#' Multi-contrast pseudobulk differential expression
#'
#' Runs a collection of two-group pseudobulk edgeR comparisons on one object, aggregating the
#' raw counts once (through [SCUBA::pseudobulk()]) and testing every contrast on the resulting
#' sample x group profiles. Contrasts that are not replicated enough for edgeR (see
#' [check_pseudobulk_design()]) fall back to the cell-level Wilcoxon test of [run_dge()]; the
#' `method` column of the result says which test produced each row. This is the natural way to
#' compare cell-level groups such as clusters across biological samples: those contrasts are
#' paired designs (`~ sample + group`), which edgeR handles.
#'
#' @param object A Seurat, SingleCellExperiment, AnnData or file-backed SCUBA object
#'   (`scuba_source`), or a pre-aggregated `scuba_pseudobulk` from [SCUBA::pseudobulk()] (then
#'   `group_by` and `sample_by` are taken from the object and no fallback is possible).
#' @param group_by Metadata column containing the groups to compare.
#' @param sample_by Metadata column containing biological sample identifiers.
#' @param contrast_mode Comparison strategy: `"pairwise"` compares every pair
#' of selected groups; `"one_vs_rest"` compares each selected group with all
#' other observed groups pooled within each sample; `"reference"` compares
#' every selected non-reference group with `reference_group`.
#' @param groups Optional character vector of groups to use. The default uses
#' every non-missing value observed in `group_by`, in observed order. In
#' `one_vs_rest` mode, this selects numerator groups; the rest side still pools
#' all other observed groups.
#' @param reference_group Reference denominator required when
#' `contrast_mode = "reference"`.
#' @param layer Raw-count layer. Defaults to `"counts"`. For AnnData objects whose counts live
#'   in `X` or `raw`, pass `"X"` or `"raw"` explicitly.
#' @param seurat_assay Assay / modality holding the raw counts (default: the object's default
#'   assay, as SCUBA defines it).
#' @param min_cells Minimum cells required per sample/comparison-side
#' pseudobulk profile.
#' @param positive_only Whether to retain only positive log2 fold changes.
#' @param remove_raw_pval Whether to remove the raw p-value column.
#' @param robust Whether to use robust empirical Bayes estimation in
#' `edgeR::glmQLFit()`.
#' @param test_use Only `"edgeR"` is supported here; cell-level tests are the per-contrast
#'   fallback, not a choice.
#' @param p_adjust_scope Multiple-testing scope. `"contrast"` retains the within-contrast
#' adjustment; `"global"` replaces `pval_adj` with BH adjustment over
#' every returned gene/contrast test; `"both"` retains within-contrast
#' `pval_adj` and adds `pval_adj_global`. Global adjustment is applied before
#' `positive_only` filtering.
#'
#' @return A tibble with the columns `group_1`, `group_2`, `contrast`, `method`, `group`,
#' `feature`, `avgExpr`, `log2FC`, `pval`, `pval_adj` (and, for rows from the Wilcoxon fallback,
#' `auc`, `pct_in`, `pct_out`; `NA` on edgeR rows). The `edger_multicontrast_details` attribute
#' holds the contrast table (with `method`, `design_type` and `reason`) and the per-contrast
#' details.
#'
#' @export
run_dge_multicontrast <- function(
  object,
  group_by,
  sample_by,
  contrast_mode = c("pairwise", "one_vs_rest", "reference"),
  groups = NULL,
  reference_group = NULL,
  layer = NULL,
  seurat_assay = NULL,
  min_cells = 10L,
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  robust = FALSE,
  p_adjust_scope = c("contrast", "global", "both"),
  test_use = "edgeR"
) {
  contrast_mode <- match.arg(contrast_mode)
  if (!.scde_is_edger(test_use)) {
    stop("`test_use` must be \"edgeR\" (cell-level tests are the per-contrast fallback).", call. = FALSE)
  }
  p_adjust_scope <- match.arg(p_adjust_scope)
  .scde_validate_min_cells(min_cells)
  .scde_validate_flag(robust, "robust")
  .scde_validate_flag(positive_only, "positive_only")
  .scde_validate_flag(remove_raw_pval, "remove_raw_pval")

  pre_aggregated <- inherits(object, "scuba_pseudobulk")
  if (pre_aggregated) {
    pb <- object
    profiles <- pb$profiles
    if (all(is.na(profiles$group))) {
      stop("The scuba_pseudobulk object has no groups; aggregate with `group_by`.", call. = FALSE)
    }
    group_by <- pb$group_by
    sample_by <- pb$sample_by
  } else {
    .scde_validate_scalar_name(group_by, "group_by")
    .scde_validate_scalar_name(sample_by, "sample_by")
    layer <- layer %||% "counts"
    .scde_check_counts_layer(layer)
    design <- SCUBA::pseudobulk_design(object, sample_by, group_by)
    profiles <- design$profiles
  }
  .scde_reject_cell_level(profiles, sample_by %||% "sample")
  observed_groups <- unique(as.character(profiles$group))

  if (is.null(groups)) groups <- observed_groups
  if (!is.character(groups) || length(groups) == 0L || anyNA(groups) ||
      any(!nzchar(groups)) || anyDuplicated(groups)) {
    stop("`groups` must contain unique, non-missing group names.", call. = FALSE)
  }
  missing_groups <- setdiff(groups, observed_groups)
  if (length(missing_groups)) {
    stop("Requested groups were not found in '", group_by, "': ",
         paste(missing_groups, collapse = ", "), ".", call. = FALSE)
  }

  contrast_table <- switch(
    contrast_mode,
    pairwise = {
      if (length(groups) < 2L) stop("Pairwise mode requires at least two groups.", call. = FALSE)
      pairs <- utils::combn(groups, 2L)
      data.frame(group_1 = pairs[1, ], group_2 = pairs[2, ], stringsAsFactors = FALSE)
    },
    reference = {
      if (!is.character(reference_group) || length(reference_group) != 1L ||
          is.na(reference_group) || !nzchar(reference_group)) {
        stop("`reference_group` is required in reference mode.", call. = FALSE)
      }
      if (!reference_group %in% observed_groups) {
        stop("Reference group '", reference_group, "' was not found in '",
             group_by, "'.", call. = FALSE)
      }
      numerators <- setdiff(groups, reference_group)
      if (!length(numerators)) {
        stop("Reference mode requires at least one non-reference group.", call. = FALSE)
      }
      data.frame(group_1 = numerators, group_2 = reference_group,
                 stringsAsFactors = FALSE)
    },
    one_vs_rest = {
      if (length(observed_groups) < 2L) stop("One-vs-rest mode requires at least two observed groups.", call. = FALSE)
      data.frame(group_1 = groups, group_2 = "rest", stringsAsFactors = FALSE)
    }
  )

  # Replication rules per contrast, from metadata alone
  verdict <- check_pseudobulk_design(profiles, contrast_table, min_cells)
  contrast_table <- verdict$contrasts
  if (verdict$any_fallback && pre_aggregated) {
    stop(
      "The following contrast(s) are not replicated enough for edgeR and a pre-aggregated ",
      "object has no cells for the Wilcoxon fallback: ",
      paste(contrast_table$contrast[contrast_table$method == "wilcoxon"], collapse = "; "),
      call. = FALSE
    )
  }

  # Aggregate once, only when some contrast will use edgeR
  if (verdict$any_edger) {
    .scde_require_edger()
    if (!pre_aggregated) {
      pb <- SCUBA::pseudobulk(object, sample_by = sample_by, group_by = group_by,
                              assay = seurat_assay, layer = layer)
      attr(pb, "scde_source_class") <- class(object)[1]
    }
    .scde_validate_pseudobulk_counts(pb)
  } else {
    pb <- NULL
  }

  details <- vector("list", nrow(contrast_table))
  names(details) <- contrast_table$contrast
  warning_messages <- vector("list", nrow(contrast_table))
  names(warning_messages) <- contrast_table$contrast
  cell_level_all <- NULL      # one-vs-rest fallbacks share one cell-level run

  result_list <- lapply(seq_len(nrow(contrast_table)), function(i) {
    group_1 <- contrast_table$group_1[[i]]
    group_2 <- contrast_table$group_2[[i]]
    selection <- verdict$selections[[i]]
    result <- tryCatch(
      withCallingHandlers(
        {
          if (selection$ok) {
            counts <- .scde_contrast_counts(pb, selection)
            fit <- .scde_test_edger(counts, selection$retained, group_1, group_2,
                                    selection$design_type, robust)
            .scde_pseudobulk_result(fit, selection, group_1, group_2, min_cells,
                                    positive_only = FALSE, remove_raw_pval = FALSE, pb, robust)
          } else {
            if (identical(group_2, "rest") && is.null(cell_level_all)) {
              cell_level_all <<- run_dge(object, group_by = group_by, seurat_assay = seurat_assay)
            }
            .scde_cell_level_contrast(
              object, group_by, group_1, group_2, assay = seurat_assay, lfc_format = "log2",
              positive_only = FALSE, remove_raw_pval = FALSE, selection = selection,
              min_cells = min_cells, cell_level = cell_level_all
            )
          }
        },
        warning = function(cnd) {
          warning_messages[[i]] <<- c(warning_messages[[i]], conditionMessage(cnd))
          invokeRestart("muffleWarning")
        },
        message = function(cnd) invokeRestart("muffleMessage")
      ),
      error = function(cnd) {
        stop(
          contrast_table$method[[i]], " contrast '", contrast_table$contrast[[i]], "' failed: ",
          conditionMessage(cnd),
          call. = FALSE
        )
      }
    )
    details[[i]] <<- attr(result, "edger_details")
    result$group_1 <- group_1
    result$group_2 <- group_2
    result$contrast <- contrast_table$contrast[[i]]
    dplyr::relocate(result, dplyr::all_of(c("group_1", "group_2", "contrast", "method")))
  })

  results <- dplyr::bind_rows(result_list)
  if (p_adjust_scope %in% c("global", "both")) {
    global_fdr <- stats::p.adjust(results$pval, method = "BH")
    if (identical(p_adjust_scope, "global")) {
      results$pval_adj <- global_fdr
    } else {
      results$pval_adj_global <- global_fdr
    }
  }
  # filters apply after the global adjustment so the set of adjusted tests is complete
  if (positive_only) results <- dplyr::filter(results, .data$log2FC > 0)
  if (remove_raw_pval) results <- dplyr::select(results, -dplyr::all_of("pval"))
  results <- dplyr::arrange(results, .data$contrast, .data$pval_adj, dplyr::desc(abs(.data$log2FC)))

  warnings_by_contrast <- warning_messages[lengths(warning_messages) > 0L]
  if (length(warnings_by_contrast)) {
    warning(
      paste(
        vapply(names(warnings_by_contrast), function(name) {
          paste0(name, ": ", paste(unique(warnings_by_contrast[[name]]), collapse = " "))
        }, character(1)),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  if (verdict$any_fallback) {
    fell_back <- contrast_table[contrast_table$method == "wilcoxon", , drop = FALSE]
    message(
      "Cell-level Wilcoxon was used for ", nrow(fell_back), " contrast(s) that edgeR could not test: ",
      paste0(fell_back$contrast, " (", fell_back$reason, ")", collapse = "; ")
    )
  }

  attr(results, "edger_multicontrast_details") <- list(
    method = "edgeR",
    contrast_mode = contrast_mode,
    contrast_table = contrast_table,
    p_adjust_scope = p_adjust_scope,
    min_cells = as.integer(min_cells),
    contrasts = details,
    warnings = warning_messages,
    aggregation = if (is.null(pb)) NULL else list(assay = pb$assay, layer = pb$layer, fun = pb$fun, engine = pb$engine,
                                                  source_class = attr(pb, "scde_source_class") %||% NA_character_)
  )
  class(results) <- c("scDE_edger_multicontrast_results", class(results))
  results
}
