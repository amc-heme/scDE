# Cell-level Wilcoxon stays the default of run_dge().

test_that("Wilcoxon remains the default", {
  skip_if_not_installed("SingleCellExperiment")
  fixture <- edger_fixture(FALSE)
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = log1p(Matrix::Matrix(fixture$counts, sparse = TRUE))),
    colData = fixture$metadata
  )
  observed <- run_dge(object, group_by = "condition")
  expect_false(inherits(observed, "scDE_edger_results"))
  expect_true(all(c("group", "feature", "log2FC", "pval_adj") %in% colnames(observed)))
})
