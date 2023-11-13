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
#' @param lfc_format sets the format used for reporting log-fold changes. If
#' "default", the default formats for each test used will be used. If "log2",
#' log-fold changes will be reported in a log two fold change format, and if
#' "ln", log-fold changes will be reported as a natural log fold change.
#'
#' @rdname run_dge
#'
#' @export
run_dge <-
  function(
    object,
    group_by,
    lfc_format = "default",
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
    lfc_format = "default"
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
    lfc_format = "default"
  ){
    # Run presto
    dge_table <-
      presto::wilcoxauc(
        AML_Seurat,
        group_by = group_by
      )

    print(colnames(dge_table))

    # Convert to tibble, remove wilcoxon rank sum U statistic
    dge_table <-
      dge_table %>%
      as_tibble() %>%
      dplyr::select(-statistic)

    # Convert to log2FC format if desired by user
    if (lfc_format == "log2"){
      dge_table$logFC <-
        scDE:::to_log2(presto_cell_type$logFC)
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
    lfc_format = "default"
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

    dge_table |>
      as_tibble() |>
      rename(all_of(rename_cols)) |>
      dplyr::select(!scores)
  }

