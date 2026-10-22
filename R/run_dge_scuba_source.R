#' @describeIn run_dge File-backed SCUBA objects (`scuba_source`, from
#' `SCUBA::open_h5ad()` / `SCUBA::open_h5mu()`): uses Presto
#' (`presto::wilcoxauc`) on the expression matrix served by SCUBA's Julia
#' engine. The matrix is pulled through `SCUBA::fetch_matrix()` in blocks of
#' features (`options(scDE.chunk_features)`, default 2000) so memory stays
#' bounded for large objects; the Wilcoxon statistics are computed per feature,
#' so blocking does not change them, and the Benjamini-Hochberg adjustment is
#' applied over all features within each group afterwards, as Presto does.
#' Groups come from `SCUBA::fetch_metadata()`, so metadata columns added in R
#' (e.g. by `SCUBA::update_object_metadata()`) are honoured. `seurat_assay`
#' selects the modality (default: the main modality, `"X"` for h5ad files);
#' `layer` defaults to `SCUBA::default_layer()` (`"X"`), which is assumed to
#' hold log-normalised values like a Seurat `data` layer. Output schema matches
#' the Seurat/Presto path: `group, feature, avgExpr, log2FC, auc, pval,
#' pval_adj, pct_in, pct_out`.
#'
#' @export
run_dge.scuba_source <-
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
  ) {
    if (!is.null(test_use) && tolower(test_use) != "presto") {
      stop(
        "run_dge() for scuba_source objects supports test_use = \"Presto\" only."
      )
    }

    assay <- seurat_assay %||% SCUBA::scuba_source_modalities(object)[1]
    layer <- layer %||% SCUBA::default_layer(object)

    # Group labels per cell (in the object's cell order), including metadata
    # columns that only exist in the R-side table.
    groups <-
      SCUBA::fetch_metadata(object, vars = group_by, return_class = "vector")
    groups <- as.character(groups)
    if (anyNA(groups)) {
      stop(
        "Column '", group_by, "' contains ", sum(is.na(groups)),
        " NA value(s). Remove or impute NA cells before running run_dge."
      )
    }
    if (length(unique(groups)) < 2) {
      stop(
        "At least 2 groups are required for differential expression. ",
        "Column '", group_by, "' contains only 1 unique value: '",
        unique(groups), "'."
      )
    }

    features <- SCUBA::features_in_assay(object, assay)
    chunk <- as.integer(getOption("scDE.chunk_features", 2000L))
    blocks <- split(features, ceiling(seq_along(features) / chunk))

    # presto::wilcoxauc on a features x cells matrix; statistics are per
    # feature, so feature blocks are independent. padj is recomputed below
    # over all features (presto adjusts within each group over the features
    # it is given).
    dge_table <-
      lapply(blocks, function(fs) {
        mat <-
          SCUBA::fetch_matrix(
            object,
            features = fs,
            assay = assay,
            layer = layer,
            sparse = TRUE
          )
        presto::wilcoxauc(mat, y = groups)
      }) %>%
      dplyr::bind_rows() %>%
      tibble::as_tibble() %>%
      dplyr::select(-statistic) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(padj = stats::p.adjust(pval, method = "BH")) %>%
      dplyr::ungroup()

    .presto_table(dge_table, lfc_format, positive_only, remove_raw_pval)
  }

#' Shared post-processing for Presto results
#'
#' Takes a `presto::wilcoxauc()` table (with `statistic` already dropped) and
#' applies the filters, the log-fold-change format and the column names /
#' order used by every Presto-based method of `run_dge()`.
#'
#' @noRd
.presto_table <- function(dge_table, lfc_format, positive_only, remove_raw_pval) {
  dge_table <-
    dge_table %>%
    tibble::as_tibble() %>%
    # Remove p-value column if specified (only shows adjusted p-value)
    {
      if (remove_raw_pval == TRUE) dplyr::select(., -pval) else .
    } %>%
    # Return only genes that are upregulated in each group
    # if positive_only is TRUE
    {
      if (positive_only == TRUE) dplyr::filter(., logFC > 0) else .
    } %>%
    # Sort by group, then by adjusted p-value, then by descending LFC
    dplyr::arrange(group, padj, desc(abs(logFC)))

  # Convert to log2FC format if desired by user
  if (lfc_format == "log2") {
    dge_table$logFC <- ln_to_log2(dge_table$logFC)
    dge_table <- dplyr::rename(dge_table, log2FC = logFC)
  }

  dge_table %>%
    dplyr::rename(pval_adj = padj) %>%
    dplyr::relocate(group, .before = feature)
}
