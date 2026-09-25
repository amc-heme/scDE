#' Multi-contrast pseudobulk differential expression with edgeR
#'
#' Runs a collection of two-group edgeR quasi-likelihood comparisons while
#' reusing scDE's raw-count extraction, pseudobulk aggregation, sample pairing,
#' and design validation. This is particularly useful for comparing cell-level
#' groups such as clusters across biological samples.
#'
#' @param object A Seurat, SingleCellExperiment, or AnnData object.
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
#' @param layer Raw-count layer or assay. Defaults to `"counts"`. For AnnData,
#' use `"X"` or `"raw"` explicitly when counts are stored there.
#' @param seurat_assay Seurat assay containing raw counts. Defaults to `"RNA"`.
#' @param min_cells Minimum cells required per sample/comparison-side
#' pseudobulk profile.
#' @param positive_only Whether to retain only positive log2 fold changes.
#' @param remove_raw_pval Whether to remove the raw p-value column.
#' @param robust Whether to use robust empirical Bayes estimation in
#' `edgeR::glmQLFit()`.
#' @param p_adjust_scope Multiple-testing scope. `"contrast"` retains edgeR's
#' within-contrast FDR; `"global"` replaces `pval_adj` with BH adjustment over
#' every returned gene/contrast test; `"both"` retains within-contrast
#' `pval_adj` and adds `pval_adj_global`.
#'
#' @return A tibble containing the common scDE result columns plus `group_1`,
#' `group_2`, and `contrast`. Per-contrast edgeR details are stored in the
#' `edger_multicontrast_details` attribute.
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
  p_adjust_scope = c("contrast", "global", "both")
) {
  contrast_mode <- match.arg(contrast_mode)
  p_adjust_scope <- match.arg(p_adjust_scope)
  .scde_validate_scalar_name(group_by, "group_by")
  .scde_validate_scalar_name(sample_by, "sample_by")
  .scde_validate_edgeR_args(
    group_by, sample_by, ".group_1", ".group_2", min_cells, "log2", robust
  )
  for (value in list(
    positive_only = positive_only,
    remove_raw_pval = remove_raw_pval
  )) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("`positive_only` and `remove_raw_pval` must be TRUE or FALSE.",
           call. = FALSE)
    }
  }

  extracted <- .scde_extract_edgeR_data(object, layer, seurat_assay)
  metadata <- extracted$metadata
  .scde_validate_count_matrix(extracted$counts, rownames(metadata))
  if (!group_by %in% colnames(metadata)) {
    stop("Metadata column '", group_by, "' was not found.", call. = FALSE)
  }
  observed <- as.character(metadata[[group_by]])
  observed_groups <- unique(observed[!is.na(observed) & nzchar(trimws(observed))])

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
    one_vs_rest = data.frame(
      group_1 = groups,
      group_2 = "rest",
      stringsAsFactors = FALSE
    )
  )
  contrast_table$contrast <- paste(
    contrast_table$group_1, "vs", contrast_table$group_2
  )

  details <- vector("list", nrow(contrast_table))
  names(details) <- contrast_table$contrast
  warning_messages <- vector("list", nrow(contrast_table))
  names(warning_messages) <- contrast_table$contrast

  result_list <- lapply(seq_len(nrow(contrast_table)), function(i) {
    display_group_1 <- contrast_table$group_1[[i]]
    display_group_2 <- contrast_table$group_2[[i]]
    analysis_metadata <- metadata
    analysis_group_by <- group_by
    analysis_group_1 <- display_group_1
    analysis_group_2 <- display_group_2

    if (identical(contrast_mode, "one_vs_rest")) {
      analysis_group_by <- ".scDE_multicontrast_group"
      while (analysis_group_by %in% colnames(analysis_metadata)) {
        analysis_group_by <- paste0(analysis_group_by, "_")
      }
      rest_label <- ".scDE_rest"
      while (rest_label %in% observed_groups) rest_label <- paste0(rest_label, "_")
      analysis_metadata[[analysis_group_by]] <- ifelse(
        observed == display_group_1,
        display_group_1,
        ifelse(is.na(observed) | !nzchar(trimws(observed)), NA_character_, rest_label)
      )
      analysis_group_2 <- rest_label
    }

    result <- tryCatch(
      withCallingHandlers(
        .scde_run_edger(
          extracted$counts,
          analysis_metadata,
          sample_by,
          analysis_group_by,
          analysis_group_1,
          analysis_group_2,
          min_cells,
          positive_only,
          FALSE,
          robust,
          validate_counts = FALSE
        ),
        warning = function(cnd) {
          warning_messages[[i]] <<- c(warning_messages[[i]], conditionMessage(cnd))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(cnd) {
        stop(
          "edgeR contrast '", contrast_table$contrast[[i]], "' failed: ",
          conditionMessage(cnd),
          call. = FALSE
        )
      }
    )

    contrast_details <- attr(result, "edger_details")
    if (identical(contrast_mode, "one_vs_rest")) {
      contrast_details$group_2 <- "rest"
      contrast_details$sample_summary$group[
        contrast_details$sample_summary$group == analysis_group_2
      ] <- "rest"
      contrast_details$retained_profiles$group[
        contrast_details$retained_profiles$group == analysis_group_2
      ] <- "rest"
      contrast_details$excluded_profiles$group[
        contrast_details$excluded_profiles$group == analysis_group_2
      ] <- "rest"
    }
    details[[i]] <<- contrast_details

    result$group_1 <- display_group_1
    result$group_2 <- display_group_2
    result$contrast <- contrast_table$contrast[[i]]
    dplyr::relocate(
      result,
      dplyr::all_of(c("group_1", "group_2", "contrast"))
    )
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
  if (remove_raw_pval) {
    results <- dplyr::select(results, -dplyr::all_of("pval"))
  }

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

  attr(results, "edger_multicontrast_details") <- list(
    contrast_mode = contrast_mode,
    contrast_table = contrast_table,
    p_adjust_scope = p_adjust_scope,
    contrasts = details,
    warnings = warning_messages
  )
  class(results) <- c("scDE_edger_multicontrast_results", class(results))
  results
}
