# File-backed SCUBA objects (scuba_source) must give the same Presto results as the equivalent
# in-memory Seurat object: same cells, X == Seurat "data" in the AML fixtures.
skip_if_no_julia()

h5ad <- system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA", mustWork = TRUE)
withr::local_options(list(SCUBA.engine = "julia"))
obj <- SCUBA::open_h5ad(h5ad)
seurat <- SCUBA::AML_Seurat

sort_table <- function(x) {
  x <- as.data.frame(x)
  x[order(x$group, x$feature), ]
}

test_that("run_dge on a scuba_source matches the Seurat/Presto result", {
  a <- run_dge(obj, group_by = "condensed_cell_type")
  b <- run_dge(seurat, group_by = "condensed_cell_type", seurat_assay = "RNA")
  expect_equal(names(a), names(b))
  expect_equal(
    names(a),
    c("group", "feature", "avgExpr", "log2FC", "auc", "pval", "pval_adj", "pct_in", "pct_out")
  )
  expect_equal(nrow(a), nrow(b))
  a <- sort_table(a); b <- sort_table(b)
  expect_equal(a$feature, b$feature)
  for (col in c("avgExpr", "log2FC", "pval", "pval_adj", "auc", "pct_in", "pct_out")) {
    expect_equal(a[[col]], b[[col]], tolerance = 1e-6, info = col)
  }
})

test_that("feature blocking does not change the result, including the BH adjustment", {
  full <- sort_table(run_dge(obj, group_by = "Batch"))
  blocked <- withr::with_options(
    list(scDE.chunk_features = 100L),
    sort_table(run_dge(obj, group_by = "Batch"))
  )
  expect_equal(blocked, full, tolerance = 1e-12)
})

test_that("options and group columns added in R are honoured", {
  pos <- run_dge(obj, group_by = "Batch", positive_only = TRUE, remove_raw_pval = TRUE)
  expect_false("pval" %in% names(pos))
  expect_true(all(pos$log2FC > 0))
  ln <- run_dge(obj, group_by = "Batch", lfc_format = "ln")
  expect_true("logFC" %in% names(ln) && !"log2FC" %in% names(ln))

  meta <- SCUBA::fetch_metadata(obj, full_table = TRUE)
  meta$flag <- ifelse(meta$nCount_RNA > stats::median(meta$nCount_RNA), "high", "low")
  obj2 <- SCUBA::update_object_metadata(obj, meta)
  res <- run_dge(obj2, group_by = "flag")
  expect_setequal(unique(res$group), c("high", "low"))

  # the protein modality by name
  prot <- run_dge(obj, group_by = "Batch", seurat_assay = "protein")
  expect_true(all(prot$feature %in% SCUBA::features_in_assay(obj, "protein")))
})

test_that("input validation", {
  expect_error(run_dge(obj, group_by = "Batch", test_use = "scran"), "Presto")
  one <- obj[SCUBA::fetch_cells(obj, "Batch", "BM_200AB"), ]
  expect_error(run_dge(one, group_by = "Batch"), "At least 2 groups")
})

test_that("pseudobulk edgeR on a file-backed object matches the in-memory Seurat object", {
  skip_if_no_edger()
  # paired design: B Cells vs NK Cells within the two batches (the hard minimum of two pairs)
  args <- list(group_by = "condensed_cell_type", test_use = "edgeR", sample_by = "Batch",
               group_1 = "B Cells", group_2 = "NK Cells", min_cells = 5L)
  a <- suppressWarnings(do.call(run_dge, c(list(obj), args)))
  b <- suppressWarnings(do.call(run_dge, c(list(seurat), args, list(seurat_assay = "RNA"))))
  expect_s3_class(a, "scDE_edger_results")
  expect_identical(attr(a, "edger_details")$design_type, "paired")
  expect_identical(attr(a, "edger_details")$aggregation$engine, "julia")
  expect_equal(sorted_result(a), sorted_result(b), tolerance = 1e-8)
  # one sample per batch: the contrast falls back to the chunked Presto test of scuba_source
  f <- suppressMessages(run_dge(obj, group_by = "Batch", test_use = "edgeR", sample_by = "Batch",
                                group_1 = "BM_200AB", group_2 = "PBMC_200AB"))
  expect_true(all(f$method == "wilcoxon"))
  reference <- run_dge(obj, group_by = "Batch")
  reference <- reference[reference$group == "BM_200AB", ]
  expect_equal(sorted_result(f)[, c("feature", "pval")], sorted_result(reference)[, c("feature", "pval")], tolerance = 1e-10)
})
