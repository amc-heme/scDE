#' Differential Gene Expression
#'
#' Performs differential expression testing on single cell data depending on the
#' class of the object passed. See the "Methods (by class)" section for details
#' on the test used.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param group_by metadata variable to use for forming differential gene
#' expression groups. A group will be created for all values in this variable
#' present in the object passed to this function.
#' @param slot The feature matrix slot to pull data from. If NULL, the "data"
#' slot is used for Seurat objects, and the "logcounts" slot is used for
#' SingleCellExperiment objects. This parameter is currently ignored for
#' anndata objects: only the "X" matrix may be used.
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
    seurat_assay = NULL,
    slot = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE,
    ...
  ){
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
    seurat_assay = NULL,
    slot = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE
  ){
    warning(
      paste0(
        "run_dge does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and Anndata."
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
    seurat_assay = NULL,
    slot = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE
  ){
    # Define slot to use. If not specified by the user, use the
    # default slot, "data"
    slot <- slot %||% SCUBA::default_slot(object)
    # Seurat assay to use: "RNA" if not specified
    seurat_assay <- seurat_assay %||% "RNA"

    # If matrix at the defined slot is a BPCells matrix, use the BPCells
    # marker_feature function.
    # If a dGC or dense matrix, use presto.
    # Code to test class of expression matrix depends on the structure of the
    # assay (v3 vs. v5)
    if (inherits(object[[seurat_assay]], "Assay5")){
      if (inherits(object[[seurat_assay]]@layers[[slot]], "IterableMatrix")){
        test_use <- "BPCells"
      } else {
        test_use <- "Presto"
      }
    } else if (inherits(object[[seurat_assay]], "Assay")){
      # Uses single brackets to pull via @ (@counts or @data)
      if (inherits(object[[seurat_assay]][slot], "IterableMatrix")){
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
          mat = object[[seurat_assay]][[slot]],
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
            # Matrices at all slots except counts are ln-transformed
            if (slot == "data"){
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
          padj = stats::p.adjust(p_val_raw, method = "BH"),
          # 1.3. Add placeholder NA entries for pct_in, pct_out
          # These are not computed by marker_features. These are included for
          # consistency with the outputs of other DGE methods used in
          # this package.
          pct_in = NA_character_,
          pct_out = NA_character_
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
        dplyr::arrange(group, padj, desc(abs(logFC))) %>%
        # 2.4. Move columns for consistency with outputs from other DGE methods
        dplyr::relocate(
          c(avgExpr, logFC, pval, padj),
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
          assay = slot,
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
          dplyr::rename(Log2FC = logFC)
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

# @describeIn run_dge SingleCellExperiment objects
# @export
# run_dge.SingleCellExperiment <-
#   function(
#     object
#   ){
#     # SingleCellExperiment objects: cell names are column names of object
#     colnames(object)
#   }

#' @describeIn run_dge Anndata objects: uses Scanpy's rank_genes_groups
#' function, with the "wilcoxon" method
#'
#' @export
run_dge.AnnDataR6 <-
  function(
    object,
    group_by,
    seurat_assay = NULL,
    slot = NULL,
    lfc_format = "log2",
    positive_only = FALSE,
    remove_raw_pval = FALSE
  ){
    library(reticulate)

    # Source fetch_anndata python script
    python_path =
      system.file(
        "extdata",
        "Python",
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

