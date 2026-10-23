# Pseudobulk edgeR backend ---------------------------------------------------------------------
#
# Aggregation lives in SCUBA (`SCUBA::pseudobulk()`), so nothing in this file touches a cell-level
# matrix. For a contrast, the backend
#   1. previews the sample x group profiles from metadata alone (`SCUBA::pseudobulk_design()`),
#   2. decides with `check_pseudobulk_design()` whether the design is replicated enough for edgeR
#      (the same rules cellDIVER shows before a run, so preflight and run always agree),
#   3. aggregates once, runs the edgeR quasi-likelihood test on the small features x profiles
#      matrix (independent design: ~ 0 + group with a contrast; paired design, where samples
#      contribute cells to both groups: ~ sample + group), and
#   4. falls back to the existing cell-level Wilcoxon test of `run_dge()` when the contrast is
#      not replicated, labelling every result row with its `method`.

.scde_is_edger <- function(test_use) {
  is.character(test_use) && length(test_use) == 1L && !is.na(test_use) && tolower(test_use) == "edger"
}

.scde_require_edger <- function() {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(
      "Package 'edgeR' (>= 4.0.0) is required for test_use = 'edgeR'. Install with: ",
      "BiocManager::install('edgeR')",
      call. = FALSE
    )
  }
  if (utils::packageVersion("edgeR") < "4.0.0") {
    stop(
      "edgeR >= 4.0.0 is required for test_use = 'edgeR' (found ",
      utils::packageVersion("edgeR"), "). Update with: BiocManager::install('edgeR')",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Prevent a supplied biological-sample column from looking effective when a
# cell-level backend was selected. The pseudobulk path returns before this
# helper is called.
.scde_warn_ignored_sample_by <- function(sample_by) {
  if (!is.null(sample_by)) {
    warning(
      "`sample_by` is ignored by cell-level differential-expression tests. ",
      "Use `test_use = \"edgeR\"` for sample-level (pseudobulk) testing, or ",
      "omit `sample_by` for cell-level Wilcoxon.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.scde_validate_scalar_name <- function(x, argument) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", argument, "` must be one non-missing column name.", call. = FALSE)
  }
}

.scde_validate_min_cells <- function(min_cells) {
  if (!is.numeric(min_cells) || length(min_cells) != 1L || is.na(min_cells) ||
      !is.finite(min_cells) || min_cells < 1 || min_cells != as.integer(min_cells)) {
    stop("`min_cells` must be one positive whole number.", call. = FALSE)
  }
}

.scde_validate_flag <- function(x, argument) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

.scde_validate_pseudobulk_args <- function(group_1, group_2, min_cells, lfc_format, robust) {
  for (x in list(group_1 = group_1, group_2 = group_2)) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop("`group_1` and `group_2` must each be one non-missing group name.",
           call. = FALSE)
    }
  }
  if (identical(group_1, group_2)) {
    stop("`group_1` and `group_2` must be different.", call. = FALSE)
  }
  .scde_validate_min_cells(min_cells)
  if (!identical(lfc_format, "log2")) {
    stop("Pseudobulk tests report log2 fold changes; use `lfc_format = \"log2\"`.",
         call. = FALSE)
  }
  .scde_validate_flag(robust, "robust")
}

# Layers whose name says they hold normalised values must not be summed.
.scde_check_counts_layer <- function(layer) {
  if (layer %in% c("X", "raw", "counts")) return(invisible(layer))
  if (layer %in% c("data", "scale.data", "logcounts", "normcounts", "cpm", "logcpm") ||
      grepl("log|norm|scale", layer, ignore.case = TRUE)) {
    stop(
      "Pseudobulk tests need raw counts; layer '", layer, "' holds normalised values. ",
      "Use layer = \"counts\" (or \"X\" / \"raw\" for AnnData objects whose counts live there).",
      call. = FALSE
    )
  }
  invisible(layer)
}

# Per-cell sample identifiers (median profile size below two cells) are a cell-level
# comparison in disguise; refuse them before any expression data is read.
.scde_reject_cell_level <- function(profiles, sample_by) {
  if (!nrow(profiles)) {
    stop("No cells with a valid sample identifier were found in `sample_by`.", call. = FALSE)
  }
  if (stats::median(profiles$n_cells) < 2) {
    stop(
      "`sample_by` ('", sample_by, "') yields about one cell per sample/group profile, ",
      "i.e. a cell-level comparison. Pseudobulk tests require biological sample ",
      "identifiers. For cell-level comparisons, omit `test_use` (Wilcoxon) and `sample_by`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# The aggregated matrix is small, so the integrality check that used to scan every cell block
# is now a few vector operations: a normalised layer summed over two or more cells is
# (essentially) never integer-valued.
.scde_validate_pseudobulk_counts <- function(pb) {
  counts <- pb$counts
  if (!identical(pb$fun, "sum")) {
    stop("Pseudobulk tests need summed counts; the object was aggregated with fun = \"", pb$fun, "\".", call. = FALSE)
  }
  if (!is.matrix(counts) || !is.numeric(counts) || !length(counts)) {
    stop("The pseudobulk object holds no count matrix.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || anyNA(rownames(counts)) || anyDuplicated(rownames(counts))) {
    stop("Pseudobulk counts must have unique, non-missing feature names.", call. = FALSE)
  }
  if (!identical(colnames(counts), pb$profiles$profile_id)) {
    stop("Pseudobulk count columns do not match the profile table.", call. = FALSE)
  }
  if (any(!is.finite(counts))) stop("Pseudobulk counts contain non-finite values.", call. = FALSE)
  if (any(counts < 0)) stop("Pseudobulk counts contain negative values.", call. = FALSE)
  if (any(abs(counts - round(counts)) > 1e-6)) {
    stop(
      "Pseudobulk profiles contain non-integer values, so layer '", pb$layer,
      "' does not hold raw counts. Pass layer = \"counts\" (or \"X\" / \"raw\" for AnnData).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# --- entry point ------------------------------------------------------------------------------

# `object` is any SCUBA-supported object, or a pre-aggregated `scuba_pseudobulk` (no cells, so
# no fallback). `assay` is the SCUBA assay / modality (`seurat_assay` in the public API).
.scde_run_pseudobulk <- function(object, group_by, sample_by, group_1, group_2, layer = NULL,
                                 assay = NULL, min_cells = 10L, robust = FALSE,
                                 lfc_format = "log2", positive_only = FALSE,
                                 remove_raw_pval = FALSE) {
  .scde_validate_pseudobulk_args(group_1, group_2, min_cells, lfc_format, robust)
  .scde_validate_flag(positive_only, "positive_only")
  .scde_validate_flag(remove_raw_pval, "remove_raw_pval")
  pre_aggregated <- inherits(object, "scuba_pseudobulk")
  if (pre_aggregated) {
    pb <- object
    profiles <- pb$profiles
    if (all(is.na(profiles$group))) {
      stop("The scuba_pseudobulk object has no groups; aggregate with `group_by` to test a contrast.",
           call. = FALSE)
    }
    group_by <- pb$group_by
    sample_by <- pb$sample_by
  } else {
    .scde_validate_scalar_name(group_by, "group_by")
    .scde_validate_scalar_name(sample_by, "sample_by")
    layer <- layer %||% "counts"
    .scde_check_counts_layer(layer)
    # metadata only: nothing is read from the expression matrix yet
    design <- SCUBA::pseudobulk_design(object, sample_by, group_by)
    profiles <- design$profiles
  }
  .scde_reject_cell_level(profiles, sample_by %||% "sample")
  .scde_check_groups(profiles, c(group_1, group_2), group_by)

  selection <- .scde_contrast_profiles(profiles, group_1, group_2, min_cells)
  if (!selection$ok) {
    if (pre_aggregated) {
      stop(
        "edgeR cannot run for ", group_1, " vs ", group_2, ": ", selection$reason,
        " A pre-aggregated object has no cells, so the cell-level Wilcoxon fallback is unavailable.",
        call. = FALSE
      )
    }
    message(
      "edgeR cannot run for ", group_1, " vs ", group_2, ": ", selection$reason,
      " Falling back to the cell-level Wilcoxon test."
    )
    return(.scde_cell_level_contrast(
      object, group_by, group_1, group_2, assay = assay, lfc_format = lfc_format,
      positive_only = positive_only, remove_raw_pval = remove_raw_pval,
      selection = selection, min_cells = min_cells
    ))
  }

  .scde_require_edger()
  if (!pre_aggregated) {
    pb <- SCUBA::pseudobulk(object, sample_by = sample_by, group_by = group_by,
                            assay = assay, layer = layer)
    attr(pb, "scde_source_class") <- class(object)[1]
  }
  .scde_validate_pseudobulk_counts(pb)
  counts <- .scde_contrast_counts(pb, selection)
  fit <- .scde_test_edger(counts, selection$retained, group_1, group_2, selection$design_type, robust)
  .scde_pseudobulk_result(fit, selection, group_1, group_2, min_cells, positive_only,
                          remove_raw_pval, pb, robust)
}

.scde_check_groups <- function(profiles, groups, group_by) {
  observed <- unique(profiles$group)
  groups <- setdiff(groups, "rest")
  missing_groups <- setdiff(groups, observed)
  if (length(missing_groups)) {
    stop(
      "Group(s) not found in '", group_by, "': ", paste(missing_groups, collapse = ", "),
      ". Observed groups: ", paste(observed, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# --- profile selection and replication rules --------------------------------------------------

# Decide which sample x group profiles enter a contrast and whether edgeR can run on them.
# `profiles`: data.frame(profile_id, sample, group, n_cells) from SCUBA. `group_2 = "rest"`
# pools, per sample, every profile of a group other than `group_1` (the pooled count is what
# `min_cells` is applied to). The design is "paired" when at least one sample has both sides
# after the `min_cells` filter (samples lacking a side are then dropped), "independent"
# otherwise. edgeR needs `min_replicates` samples per side (independent) or paired samples.
.scde_contrast_profiles <- function(profiles, group_1, group_2, min_cells, min_replicates = 2L) {
  p <- profiles[, c("profile_id", "sample", "group", "n_cells"), drop = FALSE]
  p$sample <- as.character(p$sample); p$group <- as.character(p$group)
  members <- as.list(p$profile_id)
  if (identical(group_2, "rest")) {
    own <- p$group == group_1
    others <- p[!own, , drop = FALSE]
    own_members <- members[own]
    p <- p[own, , drop = FALSE]
    if (nrow(others)) {
      by_sample <- split(others, others$sample)
      rest <- data.frame(
        profile_id = paste0(names(by_sample), "__rest"),
        sample = names(by_sample),
        group = "rest",
        n_cells = vapply(by_sample, function(d) sum(d$n_cells), numeric(1)),
        stringsAsFactors = FALSE
      )
      p <- rbind(p, rest)
      members <- c(own_members, lapply(by_sample, function(d) d$profile_id))
    } else {
      members <- own_members
    }
  } else {
    keep <- p$group %in% c(group_1, group_2)
    p <- p[keep, , drop = FALSE]
    members <- members[keep]
  }
  rownames(p) <- NULL
  p$n_cells <- as.integer(p$n_cells)
  p$pass <- p$n_cells >= min_cells
  samples_1 <- unique(p$sample[p$group == group_1 & p$pass])
  samples_2 <- unique(p$sample[p$group == group_2 & p$pass])
  both <- intersect(samples_1, samples_2)
  design_type <- if (length(both)) "paired" else "independent"
  retained <- if (design_type == "paired") p$pass & p$sample %in% both else p$pass
  p$retained <- retained
  p$exclusion_reason <- ifelse(
    retained, NA_character_,
    ifelse(!p$pass, paste0("fewer than ", min_cells, " cells"), "paired group profile unavailable")
  )
  n_1 <- length(unique(p$sample[retained & p$group == group_1]))
  n_2 <- length(unique(p$sample[retained & p$group == group_2]))
  pairs <- if (design_type == "paired") length(both) else 0L
  ok <- if (design_type == "paired") pairs >= min_replicates else (n_1 >= min_replicates && n_2 >= min_replicates)
  reason <- if (ok) {
    NA_character_
  } else if (design_type == "paired") {
    sprintf("Only %d sample%s at least %d cells in both '%s' and '%s' (edgeR needs %d paired samples).",
            pairs, if (pairs == 1L) " has" else "s have", min_cells, group_1, group_2, min_replicates)
  } else {
    short <- c(n_1, n_2) < min_replicates
    sprintf("%s (edgeR needs %d samples per group with at least %d cells).",
            paste(sprintf("'%s' has %d sample%s", c(group_1, group_2)[short], c(n_1, n_2)[short],
                          ifelse(c(n_1, n_2)[short] == 1L, "", "s")), collapse = " and "),
            min_replicates, min_cells)
  }
  warnings <- character()
  if (ok && (if (design_type == "paired") pairs < 3L else min(n_1, n_2) < 3L)) {
    warnings <- "Fewer than three retained biological replicates are available; three or more are recommended."
  }
  excluded <- p[!p$retained, , drop = FALSE]
  if (ok && nrow(excluded)) {
    warnings <- c(warnings, paste0(nrow(excluded), " sample/group pseudobulk profile(s) were excluded."))
  }
  list(
    group_1 = group_1, group_2 = group_2, min_cells = as.integer(min_cells),
    profiles = p, members = members,
    retained = p[p$retained, , drop = FALSE], retained_members = members[p$retained],
    excluded = excluded,
    design_type = design_type, samples_group_1 = n_1, samples_group_2 = n_2,
    paired_samples = pairs, ok = ok, reason = reason, warnings = warnings
  )
}

# Columns of the aggregated matrix for the retained profiles; a pooled "rest" profile is the
# row sum of its member profiles (exact: integer sums of integer sums), so nothing is re-read.
.scde_contrast_counts <- function(pb, selection) {
  retained <- selection$retained
  if (!nrow(retained)) stop("No pseudobulk profiles remain after filtering.", call. = FALSE)
  columns <- lapply(selection$retained_members, function(m) {
    if (length(m) == 1L) pb$counts[, m] else rowSums(pb$counts[, m, drop = FALSE])
  })
  counts <- matrix(unlist(columns, use.names = FALSE), nrow = nrow(pb$counts),
                   dimnames = list(rownames(pb$counts), retained$profile_id))
  counts
}

# --- edgeR ------------------------------------------------------------------------------------

# edgeR 4 quasi-likelihood pipeline: filter, then TMM (normLibSizes), then glmQLFit, which
# estimates the NB dispersion itself (no estimateDisp). `group_1` is always the numerator.
.scde_test_edger <- function(counts, profiles, group_1, group_2, design_type, robust = FALSE) {
  library_sizes <- colSums(counts)
  if (any(!is.finite(library_sizes) | library_sizes <= 0)) {
    bad <- colnames(counts)[!is.finite(library_sizes) | library_sizes <= 0]
    stop("Pseudobulk profiles have zero or invalid library sizes: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  group_factor <- factor(profiles$group, levels = c(group_2, group_1))
  y <- edgeR::DGEList(counts = counts, group = group_factor)
  keep <- edgeR::filterByExpr(y, group = group_factor)
  if (!any(keep)) stop("No genes passed edgeR::filterByExpr().", call. = FALSE)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::normLibSizes(y, method = "TMM")
  if (identical(design_type, "independent")) {
    design <- stats::model.matrix(~ 0 + group_factor)
    colnames(design) <- c("group_2", "group_1")
    contrast <- c(group_2 = -1, group_1 = 1)
    coef <- NULL
  } else {
    sample_factor <- factor(profiles$sample)
    design <- stats::model.matrix(~ sample_factor + group_factor)
    colnames(design) <- c("intercept", paste0("sample_", levels(sample_factor)[-1]), "group_1_vs_group_2")
    contrast <- NULL
    coef <- ncol(design)
  }
  rownames(design) <- profiles$profile_id
  if (qr(design)$rank < ncol(design)) stop("The edgeR design matrix is not full rank.", call. = FALSE)
  fit <- edgeR::glmQLFit(y, design, robust = robust)
  test <- if (is.null(coef)) edgeR::glmQLFTest(fit, contrast = contrast) else edgeR::glmQLFTest(fit, coef = coef)
  table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  list(
    table = data.frame(
      feature = rownames(table),
      avgExpr = as.numeric(table$logCPM),
      log2FC = as.numeric(table$logFC),
      pval = as.numeric(table$PValue),
      pval_adj = as.numeric(table$FDR),
      stringsAsFactors = FALSE
    ),
    design_matrix = design, contrast = contrast, coef = coef, genes_tested = sum(keep),
    normalization_factors = stats::setNames(y$samples$norm.factors, rownames(y$samples))
  )
}

.scde_plain_profiles <- function(p) {
  cols <- intersect(c("profile_id", "sample", "group", "n_cells", "exclusion_reason"), names(p))
  out <- p[, cols, drop = FALSE]
  rownames(out) <- NULL
  out
}

.scde_sample_summary <- function(selection) {
  g <- c(selection$group_1, selection$group_2)
  data.frame(
    group = g,
    retained_samples = c(selection$samples_group_1, selection$samples_group_2),
    excluded_samples = vapply(g, function(x) {
      length(unique(selection$excluded$sample[selection$excluded$group == x]))
    }, integer(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

.scde_pseudobulk_result <- function(fit, selection, group_1, group_2, min_cells, positive_only,
                                    remove_raw_pval, pb, robust) {
  results <- tibble::tibble(method = "edgeR", group = group_1, fit$table)
  if (positive_only) results <- dplyr::filter(results, .data$log2FC > 0)
  if (remove_raw_pval) results <- dplyr::select(results, -dplyr::all_of("pval"))
  results <- dplyr::arrange(results, .data$group, .data$pval_adj, dplyr::desc(abs(.data$log2FC)))
  for (w in selection$warnings) warning(w, call. = FALSE)
  attr(results, "edger_details") <- list(
    method = "edgeR",
    design_type = selection$design_type,
    group_1 = group_1,
    group_2 = group_2,
    min_cells = as.integer(min_cells),
    fallback_reason = NA_character_,
    sample_summary = .scde_sample_summary(selection),
    retained_profiles = .scde_plain_profiles(selection$retained),
    excluded_profiles = .scde_plain_profiles(selection$excluded),
    design_matrix = fit$design_matrix,
    contrast = fit$contrast,
    coef = fit$coef,
    genes_tested = fit$genes_tested,
    normalization_method = "TMM",
    normalization_factors = fit$normalization_factors,
    robust = robust,
    aggregation = list(assay = pb$assay, layer = pb$layer, fun = pb$fun, engine = pb$engine,
                       source_class = attr(pb, "scde_source_class") %||% NA_character_)
  )
  class(results) <- c("scDE_edger_results", class(results))
  results
}

# --- cell-level fallback ----------------------------------------------------------------------

# The existing cell-level `run_dge()` of the object's class, restricted to the contrast: for
# `group_2 = "rest"` the one-vs-rest rows of a full run (`cell_level`, computed once by the
# caller when several contrasts fall back); otherwise a run on the cells of the two groups,
# whose one-vs-rest result is exactly group_1 vs group_2.
.scde_cell_level_contrast <- function(object, group_by, group_1, group_2, assay = NULL,
                                      lfc_format = "log2", positive_only = FALSE,
                                      remove_raw_pval = FALSE, selection, min_cells,
                                      cell_level = NULL) {
  if (identical(group_2, "rest")) {
    res <- cell_level %||% run_dge(object, group_by = group_by, seurat_assay = assay, lfc_format = lfc_format)
  } else {
    cells <- SCUBA::fetch_cells(object, group_by, c(group_1, group_2))
    subset <- .scde_subset_cells(object, cells)
    res <- run_dge(subset, group_by = group_by, seurat_assay = assay, lfc_format = lfc_format)
  }
  res <- tibble::as_tibble(res)
  res <- res[as.character(res$group) == group_1, , drop = FALSE]
  lfc_col <- if (identical(lfc_format, "log2")) "log2FC" else "logFC"
  if (positive_only) res <- res[res[[lfc_col]] > 0, , drop = FALSE]
  if (remove_raw_pval && "pval" %in% names(res)) res$pval <- NULL
  res <- dplyr::arrange(res, .data$pval_adj, dplyr::desc(abs(.data[[lfc_col]])))
  res <- dplyr::relocate(tibble::add_column(res, method = "wilcoxon", .before = 1), "method", "group")
  attr(res, "edger_details") <- list(
    method = "wilcoxon",
    design_type = NA_character_,
    group_1 = group_1,
    group_2 = group_2,
    min_cells = as.integer(min_cells),
    fallback_reason = selection$reason,
    sample_summary = .scde_sample_summary(selection),
    retained_profiles = .scde_plain_profiles(selection$retained[0, , drop = FALSE]),
    excluded_profiles = .scde_plain_profiles(selection$excluded),
    design_matrix = NULL, contrast = NULL, coef = NULL,
    genes_tested = nrow(res),
    normalization_method = NA_character_, normalization_factors = NULL, robust = NA,
    aggregation = NULL
  )
  res
}

# Cells of an object by name, for every class run_dge() supports.
.scde_subset_cells <- function(object, cells) UseMethod(".scde_subset_cells")
.scde_subset_cells.default <- function(object, cells) {
  stop("Cannot subset an object of class ", paste(class(object), collapse = ", "),
       " for the cell-level fallback.", call. = FALSE)
}
.scde_subset_cells.Seurat <- function(object, cells) object[, cells]
.scde_subset_cells.SingleCellExperiment <- function(object, cells) object[, cells]
.scde_subset_cells.scuba_source <- function(object, cells) object[cells]
.scde_subset_cells.AnnDataR6 <- function(object, cells) object[cells, ]
