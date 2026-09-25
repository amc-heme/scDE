# Internal helpers for the edgeR pseudobulk backend.

.scde_require_edger <- function() {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(
      "Package 'edgeR' is required for test_use = 'edgeR'. Install with: ",
      "BiocManager::install('edgeR')",
      call. = FALSE
    )
  }
}

.scde_validate_scalar_name <- function(x, argument) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", argument, "` must be one non-missing column name.", call. = FALSE)
  }
}

.scde_validate_edgeR_args <- function(group_by, sample_by, group_1, group_2,
                                      min_cells, lfc_format, robust = FALSE) {
  .scde_validate_scalar_name(group_by, "group_by")
  .scde_validate_scalar_name(sample_by, "sample_by")
  for (x in list(group_1 = group_1, group_2 = group_2)) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop("`group_1` and `group_2` must each be one non-missing group name.",
           call. = FALSE)
    }
  }
  if (identical(group_1, group_2)) {
    stop("`group_1` and `group_2` must be different.", call. = FALSE)
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1L || is.na(min_cells) ||
      !is.finite(min_cells) || min_cells < 1 || min_cells != as.integer(min_cells)) {
    stop("`min_cells` must be one positive whole number.", call. = FALSE)
  }
  if (!identical(lfc_format, "log2")) {
    stop("edgeR reports log2 fold changes; use `lfc_format = \"log2\"`.",
         call. = FALSE)
  }
  if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
    stop("`robust` must be TRUE or FALSE.", call. = FALSE)
  }
}

.scde_validate_count_matrix <- function(counts, expected_cells = NULL,
                                        chunk_size = 1000L) {
  dimensions <- dim(counts)
  if (is.null(dimensions) || length(dimensions) != 2L || any(dimensions == 0L)) {
    stop("Raw counts must be a non-empty two-dimensional matrix-like object.",
         call. = FALSE)
  }
  if (is.null(rownames(counts)) || anyNA(rownames(counts)) ||
      any(!nzchar(rownames(counts))) || anyDuplicated(rownames(counts))) {
    stop("Raw counts must have unique, non-missing feature names.", call. = FALSE)
  }
  if (is.null(colnames(counts)) || anyNA(colnames(counts)) ||
      any(!nzchar(colnames(counts))) || anyDuplicated(colnames(counts))) {
    stop("Raw counts must have unique, non-missing cell names.", call. = FALSE)
  }
  if (!is.null(expected_cells) && !identical(colnames(counts), expected_cells)) {
    stop("Raw-count cell names or order do not match the object metadata.",
         call. = FALSE)
  }

  starts <- seq.int(1L, ncol(counts), by = chunk_size)
  tolerance <- sqrt(.Machine$double.eps)
  for (start in starts) {
    columns <- start:min(ncol(counts), start + chunk_size - 1L)
    block <- counts[, columns, drop = FALSE]
    values <- if (inherits(block, "sparseMatrix")) block@x else as.vector(as.matrix(block))
    if (!is.numeric(values)) stop("Raw counts must be numeric.", call. = FALSE)
    if (any(!is.finite(values))) stop("Raw counts contain non-finite values.", call. = FALSE)
    if (any(values < 0)) stop("Raw counts contain negative values.", call. = FALSE)
    if (any(abs(values - round(values)) > tolerance)) {
      stop("Raw counts contain non-integer values. edgeR requires untransformed counts.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

.scde_seurat_counts <- function(object, assay, layer) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required for Seurat count extraction.", call. = FALSE)
  }
  assay <- assay %||% "RNA"
  if (!assay %in% names(object@assays)) {
    stop("Assay '", assay, "' was not found. Available assays: ",
         paste(names(object@assays), collapse = ", "), ".", call. = FALSE)
  }
  layer <- layer %||% "counts"
  if (layer %in% c("data", "scale.data")) {
    stop("edgeR requires raw counts; layer '", layer,
         "' contains processed expression values.", call. = FALSE)
  }
  assay_object <- object[[assay]]
  if (inherits(assay_object, "Assay5")) {
    available <- SeuratObject::Layers(assay_object)
    if (!layer %in% available) {
      stop("Raw-count layer '", layer, "' was not found in assay '", assay, "'.",
           call. = FALSE)
    }
    counts <- SeuratObject::LayerData(assay_object, layer = layer)
  } else if (inherits(assay_object, "Assay")) {
    if (!layer %in% c("counts", "data", "scale.data")) {
      stop("Raw-count layer '", layer, "' is unavailable in this Seurat Assay.",
           call. = FALSE)
    }
    counts <- SeuratObject::GetAssayData(object, assay = assay, layer = layer)
  } else {
    stop("Unsupported Seurat assay class: ", class(assay_object)[1], ".", call. = FALSE)
  }
  metadata <- object[[]]
  list(counts = counts, metadata = metadata)
}

.scde_sce_counts <- function(object, layer) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required for SingleCellExperiment objects.",
         call. = FALSE)
  }
  layer <- layer %||% "counts"
  available <- SummarizedExperiment::assayNames(object)
  if (!layer %in% available) {
    stop("Raw-count assay '", layer, "' was not found. Available assays: ",
         paste(available, collapse = ", "), ".", call. = FALSE)
  }
  if (grepl("log|norm|scale", layer, ignore.case = TRUE)) {
    stop("edgeR requires a raw-count assay; requested assay '", layer,
         "' appears to contain processed values.", call. = FALSE)
  }
  list(
    counts = SummarizedExperiment::assay(object, layer),
    metadata = as.data.frame(SummarizedExperiment::colData(object))
  )
}

.scde_anndata_counts <- function(object, layer) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for AnnData objects.", call. = FALSE)
  }
  layer <- layer %||% "counts"
  as_r <- function(x) {
    if (inherits(x, "python.builtin.object")) reticulate::py_to_r(x) else x
  }
  matrix_py <- NULL
  if (identical(layer, "X")) {
    matrix_py <- object$X
  } else if (identical(layer, "raw")) {
    if (is.null(object$raw)) stop("AnnData `raw` counts are unavailable.", call. = FALSE)
    matrix_py <- object$raw$X
  } else {
    available <- unlist(as_r(object$layers$keys()), use.names = FALSE)
    if (!layer %in% available) {
      stop("Raw-count AnnData layer '", layer, "' was not found. Available layers: ",
           paste(available, collapse = ", "),
           ". Use layer = 'X' or 'raw' explicitly when appropriate.", call. = FALSE)
    }
    matrix_py <- object$layers$get(layer)
  }
  cell_by_feature <- as_r(matrix_py)
  counts <- Matrix::t(cell_by_feature)
  rownames(counts) <- if (identical(layer, "raw")) {
    as.character(as_r(object$raw$var_names))
  } else {
    as.character(as_r(object$var_names))
  }
  colnames(counts) <- as.character(as_r(object$obs_names))
  metadata <- as.data.frame(as_r(object$obs))
  rownames(metadata) <- colnames(counts)
  list(counts = counts, metadata = metadata)
}

.scde_extract_edgeR_data <- function(object, layer = NULL,
                                     seurat_assay = NULL) {
  if (inherits(object, "Seurat")) {
    .scde_seurat_counts(object, seurat_assay, layer)
  } else if (inherits(object, "SingleCellExperiment")) {
    .scde_sce_counts(object, layer)
  } else if (inherits(object, "AnnDataR6")) {
    .scde_anndata_counts(object, layer)
  } else {
    stop(
      "edgeR does not support objects of class ",
      paste(class(object), collapse = ", "),
      ". Supported classes: Seurat, SingleCellExperiment, and AnnData.",
      call. = FALSE
    )
  }
}

.scde_make_pseudobulk <- function(counts, metadata, sample_by, group_by,
                                  group_1, group_2, min_cells) {
  missing_columns <- setdiff(c(sample_by, group_by), colnames(metadata))
  if (length(missing_columns)) {
    stop("Metadata column(s) not found: ", paste(missing_columns, collapse = ", "),
         ".", call. = FALSE)
  }
  if (nrow(metadata) != ncol(counts) ||
      is.null(rownames(metadata)) || !identical(rownames(metadata), colnames(counts))) {
    stop("Count-matrix cells and metadata rows must have identical names and order.",
         call. = FALSE)
  }
  samples <- as.character(metadata[[sample_by]])
  groups <- as.character(metadata[[group_by]])
  valid_sample <- !is.na(samples) & nzchar(trimws(samples))
  valid_group <- !is.na(groups) & nzchar(trimws(groups))
  selected <- valid_group & groups %in% c(group_1, group_2)
  eligible <- valid_sample & selected
  if (!any(eligible)) {
    stop("No cells with valid sample IDs belong to `group_1` or `group_2`.",
         call. = FALSE)
  }

  sample_group_counts <- tapply(groups[eligible], samples[eligible],
                                function(x) length(unique(x)))
  design_type <- if (all(sample_group_counts == 1L)) "ordinary" else "blocked"
  key <- paste(nchar(samples[eligible]), samples[eligible], groups[eligible], sep = ":")
  unique_key <- unique(key)
  membership <- match(key, unique_key)
  first <- match(unique_key, key)
  info <- data.frame(
    sample_id = samples[eligible][first],
    group = groups[eligible][first],
    n_cells = as.integer(tabulate(membership, nbins = length(unique_key))),
    stringsAsFactors = FALSE
  )
  info$min_cells_pass <- info$n_cells >= min_cells
  paired <- if (design_type == "blocked") {
    table(info$sample_id[info$min_cells_pass],
          factor(info$group[info$min_cells_pass], levels = c(group_2, group_1)))
  } else NULL
  paired_samples <- if (is.null(paired)) character() else rownames(paired)[rowSums(paired > 0) == 2L]
  info$retained <- info$min_cells_pass &
    (design_type == "ordinary" | info$sample_id %in% paired_samples)
  info$exclusion_reason <- ifelse(
    info$retained, NA_character_,
    ifelse(!info$min_cells_pass, paste0("fewer than ", min_cells, " cells"),
           "paired group profile unavailable")
  )
  info$pseudobulk_id <- make.unique(paste(info$sample_id, info$group, sep = "__"))
  retained_rows <- which(info$retained)
  if (!length(retained_rows)) stop("No pseudobulk profiles remain after filtering.", call. = FALSE)
  retained_membership <- match(membership, retained_rows)
  keep_eligible <- !is.na(retained_membership)
  cell_positions <- which(eligible)[keep_eligible]
  aggregator <- Matrix::sparseMatrix(
    i = seq_along(cell_positions), j = retained_membership[keep_eligible], x = 1,
    dims = c(length(cell_positions), length(retained_rows))
  )
  aggregated <- counts[, cell_positions, drop = FALSE] %*% aggregator
  rownames(aggregated) <- rownames(counts)
  colnames(aggregated) <- info$pseudobulk_id[retained_rows]
  retained <- info[retained_rows, , drop = FALSE]
  rownames(retained) <- retained$pseudobulk_id
  list(
    counts = aggregated,
    metadata = retained,
    all_profiles = info,
    excluded_profiles = info[!info$retained, , drop = FALSE],
    design_type = design_type,
    cell_exclusions = c(
      missing_sample_id = sum(!valid_sample),
      missing_group = sum(valid_sample & !valid_group),
      group_not_selected = sum(valid_sample & valid_group & !selected)
    )
  )
}

.scde_edgeR_design <- function(pseudobulk, group_1, group_2) {
  meta <- pseudobulk$metadata
  if (pseudobulk$design_type == "ordinary") {
    n_by_group <- vapply(c(group_1, group_2), function(g) {
      length(unique(meta$sample_id[meta$group == g]))
    }, integer(1))
    if (any(n_by_group < 2L)) {
      bad <- paste0(names(n_by_group)[n_by_group < 2L], " (", n_by_group[n_by_group < 2L], ")")
      stop("edgeR requires at least two retained biological replicates per group. ",
           "Insufficient groups: ", paste(bad, collapse = ", "), ".", call. = FALSE)
    }
    group_factor <- factor(meta$group, levels = c(group_2, group_1))
    design <- stats::model.matrix(~ 0 + group_factor)
    colnames(design) <- c("group_2", "group_1")
    contrast <- c(group_2 = -1, group_1 = 1)
    retained_samples <- n_by_group
  } else {
    retained_pairs <- unique(meta$sample_id)
    if (length(retained_pairs) < 2L) {
      stop("edgeR requires at least two retained paired biological samples; found ",
           length(retained_pairs), ".", call. = FALSE)
    }
    sample_factor <- factor(meta$sample_id)
    group_indicator <- as.integer(meta$group == group_1)
    design <- stats::model.matrix(~ sample_factor + group_indicator)
    contrast <- stats::setNames(rep(0, ncol(design)), colnames(design))
    contrast["group_indicator"] <- 1
    retained_samples <- stats::setNames(rep(length(retained_pairs), 2L),
                                        c(group_1, group_2))
  }
  rownames(design) <- rownames(meta)
  if (qr(design)$rank < ncol(design)) stop("The edgeR design matrix is not full rank.", call. = FALSE)
  list(design = design, contrast = contrast, retained_samples = retained_samples)
}

.scde_run_edger <- function(counts, metadata, sample_by, group_by, group_1,
                            group_2, min_cells, positive_only,
                            remove_raw_pval, robust,
                            validate_counts = TRUE) {
  .scde_require_edger()
  if (isTRUE(validate_counts)) {
    .scde_validate_count_matrix(counts, rownames(metadata))
  }
  pseudobulk <- .scde_make_pseudobulk(
    counts, metadata, sample_by, group_by, group_1, group_2, min_cells
  )
  design_info <- .scde_edgeR_design(pseudobulk, group_1, group_2)
  library_sizes <- as.numeric(Matrix::colSums(pseudobulk$counts))
  if (any(!is.finite(library_sizes) | library_sizes <= 0)) {
    bad <- colnames(pseudobulk$counts)[!is.finite(library_sizes) | library_sizes <= 0]
    stop("Pseudobulk profiles have zero or invalid library sizes: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }

  y <- edgeR::DGEList(counts = as.matrix(pseudobulk$counts))
  y <- edgeR::calcNormFactors(y, method = "TMM")
  keep <- edgeR::filterByExpr(y, design = design_info$design)
  if (!any(keep)) stop("No genes passed edgeR::filterByExpr().", call. = FALSE)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::estimateDisp(y, design_info$design)
  fit <- edgeR::glmQLFit(y, design_info$design, robust = robust)
  test <- edgeR::glmQLFTest(fit, contrast = design_info$contrast)
  table <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table

  results <- tibble::tibble(
    group = group_1,
    feature = rownames(table),
    avgExpr = as.numeric(table$logCPM),
    log2FC = as.numeric(table$logFC),
    pval = as.numeric(table$PValue),
    pval_adj = as.numeric(table$FDR)
  )
  if (positive_only) results <- dplyr::filter(results, .data$log2FC > 0)
  if (remove_raw_pval) results <- dplyr::select(results, -dplyr::all_of("pval"))
  results <- dplyr::arrange(
    results, .data$group, .data$pval_adj, dplyr::desc(abs(.data$log2FC))
  )

  summary <- data.frame(
    group = c(group_1, group_2),
    retained_samples = as.integer(design_info$retained_samples[c(group_1, group_2)]),
    excluded_samples = vapply(c(group_1, group_2), function(g) {
      length(unique(pseudobulk$excluded_profiles$sample_id[
        pseudobulk$excluded_profiles$group == g
      ]))
    }, integer(1)),
    stringsAsFactors = FALSE
  )
  if (any(summary$retained_samples < 3L)) {
    warning("Fewer than three retained biological replicates are available; ",
            "three or more are recommended.", call. = FALSE)
  }
  if (nrow(pseudobulk$excluded_profiles)) {
    warning(nrow(pseudobulk$excluded_profiles),
            " sample/group pseudobulk profile(s) were excluded.", call. = FALSE)
  }
  attr(results, "edger_details") <- list(
    design_type = pseudobulk$design_type,
    group_1 = group_1,
    group_2 = group_2,
    min_cells = as.integer(min_cells),
    sample_summary = summary,
    retained_profiles = pseudobulk$metadata,
    excluded_profiles = pseudobulk$excluded_profiles,
    cell_exclusions = pseudobulk$cell_exclusions,
    design_matrix = design_info$design,
    contrast = design_info$contrast,
    genes_tested = sum(keep),
    normalization_method = "TMM",
    normalization_factors = stats::setNames(y$samples$norm.factors, rownames(y$samples))
  )
  class(results) <- c("scDE_edger_results", class(results))
  results
}
