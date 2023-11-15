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
#' "default", the default formats for each test used will be used. If "log2",
#' log-fold changes will be reported in a log two fold change format, and if
#' "ln", log-fold changes will be reported as a natural log fold change.
#' @param positive_only If TRUE, only display genes that are up-regulated in a
#' group compared to the reference (LFC > 0). Defaults to FALSE.
#' @param remove_raw_pval If TRUE, show only the adjusted p-value.
#' Defaults to FALSE.
#'
#' @rdname run_dge
#'
#' @export
run_dge <-
  function(
    object,
    group_by,
    seurat_assay = NULL,
    slot = NULL,
    lfc_format = "default",
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
    lfc_format = "default",
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
    lfc_format = "default",
    positive_only = FALSE,
    remove_raw_pval = FALSE
  ){
    # Define slot to use. If not specified by the user, use the
    # default slot, "data"
    slot <- slot %||% SCUBA::default_slot(object)
    # Seurat assay to use: "RNA" if not specified
    seurat_assay <- seurat_assay %||% "RNA"

    # Run presto
    dge_table <-
      presto::wilcoxauc(
        AML_Seurat,
        group_by = group_by,
        assay = slot,
        seurat_assay = seurat_assay
      )

    # Convert to tibble, remove wilcoxon rank sum U statistic
    dge_table <-
      dge_table %>%
      as_tibble() %>%
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
        scDE:::to_log2(dge_table$logFC)

      # Rename column to specify log2FC values
      dge_table <-
        dge_table %>%
        dplyr::rename(Log2FC = logFC)
    }

    # Put group before feature column
    dge_table <-
      dplyr::relocate(dge_table, group, .before = feature)

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
    lfc_format = "default",
    positive_only = FALSE,
    remove_raw_pval = FALSE
  ){
    library(reticulate)

    # Source fetch_anndata python script
    python_path =
      system.file(
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
      c("names" = "feature",
        "logfoldchanges" = "log2FC",
        "pvals" = "pval",
        "pvals_adj" = "pval_adj",
        "pct_nz_group" = "pct_in",
        "pct_nz_reference" = "pct_out")

    dge_table %>%
      as_tibble() %>%
      rename(all_of(rename_cols)) %>%
      dplyr::select(!scores) %>%
      # Remove p-value column if specified (only shows adjusted p-value)
      {if (remove_raw_pval == TRUE) dplyr::select(., -pval) else .} %>%
      # Return only genes that are upregulated in each group
      # if positive_only is TRUE
      {if (positive_only == TRUE) dplyr::filter(., logFC > 0) else .} %>%
      # Sort by group, then by adjusted p-value, then by descending LFC
      dplyr::arrange(group, pval_adj, desc(abs(log2FC)))
  }

