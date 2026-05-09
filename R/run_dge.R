#' Differential Gene Expression
#'
#' Performs differential expression testing on single cell data depending on the
#' class of the object passed. See the "Methods (by class)" section for details
#' on the test used.
#'
#' @param object a single-cell object. Currently, Seurat,
#' SingleCellExperiment, and AnnData objects are supported.
#' @param group_by metadata variable to use for forming differential gene
#' expression groups. A group will be created for all values in this variable
#' present in the object passed to this function.
#' @param layer The feature matrix layer to pull data from. If NULL,
#' test-specific conventions will be used. The "data" layer is used for Seurat
#' objects, and the "logcounts" layer is used for SingleCellExperiment objects.
#' This parameter is currently ignored for anndata objects: only the "X" matrix
#' may be used. The layer parameter is also ignored in BPCells objects.
#' @param slot No longer used. Please use `layer` instead.
#' @param seurat_assay for Seurat objects, the assay to use for evaluating
#' differences in expression. This is used only for tests on Seurat objects.
#' For SingleCellExperiment objects, the assay defaults to the main experiment,
#' and for anndata objects the "X" matrix is used.
#' @param lfc_format sets the format used for reporting log-fold changes. If
#' "log2", log-fold changes will be reported in a log two fold change format,
#' and if "ln", log-fold changes will be reported as a natural log fold change.
#' "log2" is the default value. Presto uses a natural log fold change by
#' default, so `lfc_format` should be set to "ln" if consistency with Presto is
#' desired.
#' @param positive_only If TRUE, only display genes that are up-regulated in a
#' group compared to the reference (LFC > 0). Defaults to FALSE.
#' @param remove_raw_pval If TRUE, show only the adjusted p-value.
#' Defaults to FALSE.
#' @param test_use For SingleCellExperiment objects only. Controls which
#' backend is used. `NULL` (default) auto-detects: `"Presto"` for in-memory
#' matrices (`dgCMatrix` or `matrix`), `"scran"` for all other matrix classes
#' (e.g. HDF5-backed `DelayedMatrix`, `SVT_SparseMatrix`).
#' `"Presto"` forces presto via `presto::wilcoxauc` — errors if the matrix is
#' not `dgCMatrix` or `matrix`. `"scran"` forces a true one-vs-rest Wilcoxon
#' test per group via `scran::pairwiseWilcox`, indexing the result directly
#' (returns an `auc` column; LFC from group means). `"scran_pairwise"` uses
#' `scran::findMarkers` with pairwise comparisons combined across all groups
#' (`pval.type = "any"`); no `auc` column is returned and `log2FC` is computed
#' from group vs. all-other-cells means (same method as `"scran"`).
#'
#' @rdname run_dge
#'
#' @importFrom rlang %||%
#'
#' @importFrom magrittr %>%
#'
#' @export
run_dge <-
  function(
    object,
    group_by,
    layer = NULL,
    seurat_assay = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE,
    test_use = NULL,
    slot = lifecycle::deprecated(),
    ...
  ){
    # Slot parameter has been deprecated and changed to layer for consistency
    if (lifecycle::is_present(slot)){
      lifecycle::deprecate_warn(
        when = "0.2.0",
        what = "run_dge(slot)",
        details =
          paste0(
            "Please use the `layer` parameter instead. The `slot` ",
            "parameter will be removed in 1.0.0."
            )
        )

      # Direct slot to layer
      layer <- slot
      }

    UseMethod("run_dge")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
run_dge.default <-
  function(
    object,
    group_by,
    layer = NULL,
    seurat_assay = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE,
    test_use = NULL,
    slot = lifecycle::deprecated()
  ){
    warning(
      paste0(
        "run_dge does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat, SingleCellExperiment, and AnnData."
      )
    )
  }

#' @describeIn run_dge Seurat objects: uses Presto
#'
#' @export
run_dge.Seurat <-
  function(
    object,
    group_by,
    layer = NULL,
    seurat_assay = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE,
    test_use = NULL,
    slot = lifecycle::deprecated()
  ){
    # Define layer to use. If not specified by the user, use the
    # default layer, "data"
    layer <- layer %||% SCUBA::default_layer(object)
    # Seurat assay to use: "RNA" if not specified
    seurat_assay <- seurat_assay %||% "RNA"

    # If matrix at the defined layer is a BPCells matrix, use the BPCells
    # marker_feature function.
    # If a dGC or dense matrix, use presto.
    # Code to test class of expression matrix depends on the structure of the
    # assay (v3 vs. v5)
    if (inherits(object[[seurat_assay]], "Assay5")){
      if (inherits(object[[seurat_assay]]@layers[[layer]], "IterableMatrix")){
        test_use <- "BPCells"
      } else {
        test_use <- "Presto"
      }
    } else if (inherits(object[[seurat_assay]], "Assay")){
      # Pull the matrix at the specified layer (slot in Seurat 4.x.x.) and
      # test if it is a BPCells matrix (has the class IterableMatrix)
      # object[[seurat_assay]][slot] was used originally but caused #6.
      if (layer == "counts"){
        is_BPCells <- inherits(object[[seurat_assay]]@counts, "IterableMatrix")
      } else if (layer == "data"){
        is_BPCells <- inherits(object[[seurat_assay]]@data, "IterableMatrix")
      } else if (layer == "scale.data"){
        is_BPCells <-
          inherits(object[[seurat_assay]]@scale.data, "IterableMatrix")
      } else {
        stop(
          paste0(
            'The assay passed to the object (', seurat_assay, ') is a Seurat ',
            'v3 assay. For v3 assays, `layer` must be set to "counts", ',
            '"data", or "scale.data".'
            )
          )
      }

      # Determine test to use based on results above for checking
      # for a BPCells matrix
      if (is_BPCells){
        test_use <- "BPCells"
        } else {
          test_use <- "Presto"
        }
    } else {
      stop('Unsupported assay class for assay "', seurat_assay, '".')
    }

    if (test_use == "BPCells"){
      # marker_features requires vector of group labels by cell
      groups <-
        SCUBA::fetch_metadata(
          object,
          vars = group_by,
          return_class = "vector"
          )

      # Run BPCells marker_features
      dge_table <-
        BPCells::marker_features(
          # Always uses the "data" layer
          mat = object[[seurat_assay]]$data,
          groups = groups,
          method = "wilcoxon"
          )

      # Post-DGE modifications to table
      # 1. Add LFC, adjusted p-value columns
      dge_table <-
        dge_table %>%
        dplyr::mutate(
          # 1.1. Log-fold change
          # Compute LFC based on log scaling format (ln or log2)
          # Regardless of the format, the column is named logFC until the final
          # step so downstream operations can identify the column more easily
          logFC =
            # Calculation depends on if the matrix used has been natural
            # log-transformed, and on the LFC format requested by the user
            # Matrices at all layers except counts are ln-transformed
            if (layer == "data"){
              if (lfc_format == "log2"){
                # L2FC for log-transformed matrices:
                # subtract ln-transformed means and divide by ln(2) to
                # convert to log2 format
                (foreground_mean - background_mean)/log(2)
              } else if (lfc_format == "ln"){
                # ln-FC for log-transformed matrices:
                # subtract ln-transformed means
                # foreground/background = ln(foreground) - ln(background)
                foreground_mean - background_mean
              }
            } else {
              if (lfc_format == "log2"){
                # L2FC for non-transformed matrices:
                # compute a fold-change in log-2 format
                log2(foreground_mean/background_mean)
              } else if (lfc_format == "ln"){
                # ln-FC for non-transformed matrices:
                # Take the natural log in the fold change of the means
                log(foreground_mean/background_mean)
              }
            },
          # 1.2. Add Benjamini-Hochberg p-value correction
          pval_adj = stats::p.adjust(p_val_raw, method = "BH")
        )

      # Average expression in group: convert to log-2 form?
      # dge_table <-
      #   dge_table |>
      #   dplyr::mutate(
      #     # (division is by ln(2))
      #     avgExpr = foreground_mean/log(2),
      #   )

      # 2. Filter table, rename and remove columns
      # Identify columns to be renamed in this step (old = new pairs)
      rename_cols <-
        c("group" = "foreground",
          "pval" = "p_val_raw",
          "avgExpr" = "foreground_mean"
        )

      dge_table <-
        dge_table %>%
        # 2.1. Filter DGE table for non-zero LFC values, if positive_only is TRUE
        {if (positive_only == TRUE) dplyr::filter(., logFC > 0) else .} %>%
        # 2.2. Rename columns for output consistency
        dplyr::rename(any_of(rename_cols)) %>%
        # 2.3. Sort table by group, then by adjusted p-value,
        # then by descending LFC
        dplyr::arrange(group, pval_adj, desc(abs(logFC))) %>%
        # 2.4. Move columns for consistency with outputs from other DGE methods
        dplyr::relocate(
          c(avgExpr, logFC, pval, pval_adj),
          .after = feature
          ) %>%
        # 2.5. Remove columns that are not common to outputs of
        # other DGE methods
        dplyr::select(
          -any_of(c("background", "foreground_mean", "background_mean"))
        ) %>%
        # 2.6. Remove raw p-value column if remove_raw_pval == TRUE
        {if (remove_raw_pval == TRUE) dplyr::select(., -pval) else .}
    } else if (test_use == "Presto"){
      # Run presto
      dge_table <-
        presto::wilcoxauc(
          object,
          group_by = group_by,
          assay = layer,
          seurat_assay = seurat_assay
        )

      # Convert to tibble, remove wilcoxon rank sum U statistic
      dge_table <-
        dge_table %>%
        tibble::as_tibble() %>%
        dplyr::select(-statistic) %>%
        # Remove p-value column if specified (only shows adjusted p-value)
        {if (remove_raw_pval == TRUE) dplyr::select(., -pval) else .} %>%
        # Return only genes that are upregulated in each group
        # if positive_only is TRUE
        {if (positive_only == TRUE) dplyr::filter(., logFC > 0) else .} %>%
        # Sort by group, then by adjusted p-value, then by descending LFC
        dplyr::arrange(group, padj, desc(abs(logFC)))

      # Convert to log2FC format if desired by user
      if (lfc_format == "log2"){
        dge_table$logFC <-
          scDE:::ln_to_log2(dge_table$logFC)

        # Rename column to specify log2FC values
        dge_table <-
          dge_table %>%
          dplyr::rename(log2FC = logFC)
      }

      # Rename adjusted p-value column to pval_adj
      dge_table <-
        dge_table %>%
        dplyr::rename(pval_adj = padj)

      # Put group before feature column
      dge_table <-
        dplyr::relocate(dge_table, group, .before = feature)
    }

    dge_table
  }

#' @describeIn run_dge SingleCellExperiment objects: uses Presto by default for
#' in-memory matrices, with automatic fallback to scran for DelayedArray-backed
#' matrices (e.g. HDF5-backed via `HDF5Array`). The backend can be controlled
#' explicitly via `test_use`.
#'
#' **Presto** (`test_use = "Presto"` or auto for in-memory): calls
#' `presto::wilcoxauc` directly on the extracted assay matrix. Output schema
#' matches the Seurat/Presto path and includes `auc`, `pct_in`, and `pct_out`.
#'
#' **scran one-vs-rest** (`test_use = "scran"` or auto for non-presto-compatible
#' matrices): runs a true one-vs-rest Wilcoxon test per group via
#' `scran::pairwiseWilcox` on a binary factor (group vs. all others), then
#' indexes the result directly by pair position. Returns an `auc` column. LFC
#' is computed from per-group means (log-space subtraction). Requires
#' `BiocManager::install(c("scran", "SummarizedExperiment"))`.
#'
#' **scran pairwise** (`test_use = "scran_pairwise"`): uses
#' `scran::findMarkers(pval.type = "any")`. No `auc` column. `log2FC` is
#' computed from group vs. all-other-cells means (same method as `"scran"`).
#' Requires `BiocManager::install(c("scran", "SummarizedExperiment"))`.
#'
#' In all cases, `avgExpr` for scran paths is computed via
#' `MatrixGenerics::rowMeans2()`, which supports both in-memory (`dgCMatrix`,
#' `matrix`) and `DelayedArray`-backed matrices (e.g. `HDF5Matrix`,
#' `SVT_SparseMatrix`). `other_mean` is derived algebraically from per-group
#' means — no additional matrix reads per group.
#'
#' @export
run_dge.SingleCellExperiment <-
  function(
    object,
    group_by,
    layer = NULL,
    seurat_assay = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE,
    test_use = NULL,
    slot = lifecycle::deprecated()
  ){
    # SummarizedExperiment is always needed for matrix/colData extraction
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)){
      stop(
        "Package 'SummarizedExperiment' is required for SingleCellExperiment ",
        "objects. Install with: BiocManager::install('SummarizedExperiment')"
      )
    }

    # Default layer for SCE is "logcounts"
    layer <- layer %||% "logcounts"

    # Extract expression matrix and group labels
    mat    <- SummarizedExperiment::assay(object, layer)
    groups <- SummarizedExperiment::colData(object)[[group_by]]

    if (is.null(groups)){
      stop(
        "Column '", group_by, "' not found in colData. ",
        "Available columns: ",
        paste(names(SummarizedExperiment::colData(object)), collapse = ", ")
      )
    }

    groups <- as.character(groups)

    # Presto requires dgCMatrix or base matrix; anything else routes to scran
    .is_presto_compatible <- function(mat) inherits(mat, c("dgCMatrix", "matrix"))

    # Resolve test_use: auto-detect based on matrix class if not specified
    if (is.null(test_use)){
      test_use <- if (.is_presto_compatible(mat)) "Presto" else "scran"
    }

    # Normalize test_use casing (accept any case from user)
    test_use <- switch(
      tolower(test_use),
      "presto"          = "Presto",
      "scran"           = "scran",
      "scran_pairwise"  = "scran_pairwise",
      stop(
        "'test_use' must be one of: \"Presto\", \"scran\", \"scran_pairwise\". ",
        "Got: \"", test_use, "\""
      )
    )

    # Presto requires dgCMatrix or base matrix
    if (test_use == "Presto" && !.is_presto_compatible(mat)){
      stop(
        "Presto requires an in-memory sparse matrix (dgCMatrix or matrix). ",
        "The provided assay is class '", class(mat)[1], "'. ",
        "Use test_use = 'scran' or test_use = 'scran_pairwise'."
      )
    }

    # ------------------------------------------------------------------ #
    # Presto path                                                         #
    # ------------------------------------------------------------------ #
    if (test_use == "Presto"){
      dge_table <-
        presto::wilcoxauc(mat, y = groups) %>%
        tibble::as_tibble() %>%
        dplyr::select(-statistic) %>%
        {if (remove_raw_pval == TRUE) dplyr::select(., -pval) else .} %>%
        {if (positive_only == TRUE) dplyr::filter(., logFC > 0) else .} %>%
        dplyr::arrange(group, padj, desc(abs(logFC)))

      # Convert to log2FC format if desired by user
      if (lfc_format == "log2"){
        dge_table$logFC <- scDE:::ln_to_log2(dge_table$logFC)
        dge_table <- dge_table %>% dplyr::rename(log2FC = logFC)
      }

      dge_table <-
        dge_table %>%
        dplyr::rename(pval_adj = padj) %>%
        dplyr::relocate(group, .before = feature)

      return(dge_table)
    }

    # ------------------------------------------------------------------ #
    # scran paths: guard + shared setup                                   #
    # ------------------------------------------------------------------ #
    if (!requireNamespace("scran", quietly = TRUE)){
      stop(
        "Package 'scran' is required when test_use is '", test_use, "'. ",
        "Install with: BiocManager::install('scran')"
      )
    }

    group_levels <- unique(groups)
    total_n      <- ncol(mat)
    group_counts <- table(groups)

    # Compute per-group mean expression once using MatrixGenerics::rowMeans2,
    # which supports DelayedArray/HDF5-backed matrices. This avoids re-reading
    # the matrix once per group for other_mean: instead, other_mean is derived
    # algebraically as (total_sum - group_sum) / (total_n - n_g), which is pure
    # arithmetic after this single pass of N group reads.
    group_means <-
      lapply(group_levels, function(g){
        MatrixGenerics::rowMeans2(mat[, groups == g, drop = FALSE])
      })
    names(group_means) <- group_levels

    # Precompute total sum (weighted sum of group means) for other_mean derivation
    total_sum <- Reduce(`+`, mapply(
      function(m, n) m * n,
      group_means,
      as.integer(group_counts[group_levels]),
      SIMPLIFY = FALSE
    ))

    # Returns other_mean for group g without any matrix reads
    .other_mean <- function(g){
      n_g <- as.integer(group_counts[[g]])
      (total_sum - group_means[[g]] * n_g) / (total_n - n_g)
    }

    if (test_use == "scran"){
      # ---------------------------------------------------------------- #
      # True one-vs-rest Wilcoxon via pairwiseWilcox.                    #
      # One binary test per group (g vs. all others).                    #
      # Directly index pw$statistics by pair position — avoids           #
      # combineMarkers, which fails in scran >=1.36 due to unnamed       #
      # statistics list. Returns auc column; LFC from group means.       #
      # ---------------------------------------------------------------- #
      dge_table <-
        lapply(group_levels, function(g){
          # Binary factor with explicit level order: g first, "other" second.
          # This guarantees pairwiseWilcox places g in pairs$first so
          # pair_idx is always row 1 (g vs. other), regardless of how scran
          # orders factor levels internally.
          binary_groups <- factor(
            ifelse(groups == g, g, "other"),
            levels = c(g, "other")
          )

          pw <- scran::pairwiseWilcox(mat, groups = binary_groups)

          # pairs$first == g is now guaranteed to be row 1
          pair_idx <- which(pw$pairs$first == g)
          stats_df <- as.data.frame(pw$statistics[[pair_idx]])
          features <- rownames(pw$statistics[[pair_idx]])

          # LFC from log-space means: subtraction = log ratio;
          # divide by log(2) to convert ln-ratio -> log2FC.
          # other_mean derived algebraically — no additional matrix reads.
          log2fc_vals <-
            (group_means[[g]][features] - .other_mean(g)[features]) / log(2)

          tibble::tibble(
            group    = g,
            feature  = features,
            avgExpr  = group_means[[g]][features],
            auc      = stats_df[["AUC"]],
            log2FC   = log2fc_vals,
            pval     = stats_df[["p.value"]],
            pval_adj = stats_df[["FDR"]]
          )
        }) %>%
        dplyr::bind_rows()

    } else if (test_use == "scran_pairwise"){
      # ---------------------------------------------------------------- #
      # scran::findMarkers with pairwise combinations.                   #
      # p-values/ranking from all pairwise comparisons (pval.type="any"). #
      # No auc column. log2FC computed from group vs. all-other means    #
      # (same as "scran" path) — findMarkers wilcox output does not      #
      # include a reliable summary LFC (scran bug, unfixed as of 1.40).  #
      # ---------------------------------------------------------------- #
      markers <-
        scran::findMarkers(
          mat,
          groups     = groups,
          test.type  = "wilcox",
          full.stats = TRUE,
          pval.type  = "any"
        )

      dge_table <-
        lapply(names(markers), function(g){
          df       <- as.data.frame(markers[[g]])
          features <- rownames(markers[[g]])

          other_mean   <- .other_mean(g)
          log2fc_vals  <-
            (group_means[[g]][features] - other_mean[features]) / log(2)

          tibble::tibble(
            group    = g,
            feature  = features,
            avgExpr  = group_means[[g]][features],
            log2FC   = log2fc_vals,
            pval     = df[["p.value"]],
            pval_adj = df[["FDR"]]
          )
        }) %>%
        dplyr::bind_rows()

    } else {
      stop(
        "Unknown test_use value '", test_use, "' for SingleCellExperiment. ",
        "Must be one of: NULL (auto), 'Presto', 'scran', 'scran_pairwise'."
      )
    }

    # ------------------------------------------------------------------ #
    # Shared post-processing for scran paths                              #
    # ------------------------------------------------------------------ #

    # Handle lfc_format: convert log2FC -> ln if requested
    if (lfc_format == "ln"){
      dge_table <-
        dge_table %>%
        dplyr::mutate(logFC = log2FC * log(2)) %>%
        dplyr::select(-log2FC)

      lfc_col <- "logFC"
    } else {
      lfc_col <- "log2FC"
    }

    dge_table %>%
      {if (positive_only) dplyr::filter(., .data[[lfc_col]] > 0) else .} %>%
      {if (remove_raw_pval) dplyr::select(., -pval) else .} %>%
      dplyr::arrange(group, pval_adj, dplyr::desc(abs(.data[[lfc_col]])))
  }

#' @describeIn run_dge Anndata objects: uses Scanpy's rank_genes_groups
#' function, with the "wilcoxon" method
#'
#' @export
run_dge.AnnDataR6 <-
  function(
    object,
    group_by,
    layer = NULL,
    seurat_assay = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE,
    test_use = NULL,
    slot = lifecycle::deprecated()
  ){
    if (!requireNamespace("reticulate", quietly = TRUE)){
      stop(
        "Package 'reticulate' is required for AnnData objects. ",
        "Install with: install.packages('reticulate')"
      )
    }

    # Source fetch_anndata python script
    python_path =
      system.file(
        "extdata",
        "python",
        "scanpy_dge.py",
        package = "scDE"
      )

    reticulate::source_python(python_path)

    # Run DGE via scanpy.tl.rank_genes_groups,
    # Return table
    dge_table <-
      py$scanpy_dge(
        adata = object,
        group_by = group_by,
        method = "wilcoxon"
        )

    rename_cols <-
      c("feature" = "names",
        "log2FC" = "logfoldchanges",
        "pval" = "pvals",
        "pval_adj" = "pvals_adj",
        "pct_in" = "pct_nz_group",
        "pct_out" = "pct_nz_reference"
        )

    dge_table %>%
      tibble::as_tibble() %>%
      dplyr::rename(all_of(rename_cols)) %>%
      dplyr::select(!scores) %>%
      # Remove p-value column if specified (only shows adjusted p-value)
      {if (remove_raw_pval == TRUE) dplyr::select(., -pval) else .} %>%
      # Return only genes that are upregulated in each group
      # if positive_only is TRUE
      {if (positive_only == TRUE) dplyr::filter(., log2FC > 0) else .} %>%
      # Sort by group, then by adjusted p-value, then by descending LFC
      dplyr::arrange(group, pval_adj, desc(abs(log2FC)))
  }
