test_that("Seurat Wilcoxon bypasses the defunct Presto slot adapter", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("presto")

  counts <- matrix(
    c(
      8, 7, 1, 0,
      0, 1, 7, 8,
      3, 2, 3, 2
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(
      paste0("gene", 1:3),
      paste0("cell", 1:4)
    )
  )
  metadata <- data.frame(
    cluster = c("A", "A", "B", "B"),
    row.names = colnames(counts)
  )
  object <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts, meta.data = metadata)
  )
  result <- expect_no_error(
    run_dge(
      object,
      group_by = "cluster",
      layer = "counts",
      seurat_assay = "RNA"
    )
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("group", "feature", "log2FC", "pval_adj") %in% names(result)))
  expect_setequal(unique(result$group), c("A", "B"))
})
