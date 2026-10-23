# The edgeR backend on exact, hand-aggregated profiles (no SCUBA engine involved).

test_that("independent design: ~ 0 + group with a contrast, group_1 as numerator", {
  skip_if_no_edger()
  pb <- fixture_pseudobulk(edger_fixture(FALSE))
  result <- suppressWarnings(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"))
  expect_edger_result(result, "independent")
  details <- attr(result, "edger_details")
  expect_identical(dim(details$design_matrix), c(6L, 2L))
  expect_identical(colnames(details$design_matrix), c("group_2", "group_1"))
  expect_equal(details$contrast, c(group_2 = -1, group_1 = 1))
  expect_null(details$coef)
  expect_identical(details$sample_summary$group, c("group_1", "group_2"))
  expect_identical(details$sample_summary$retained_samples, c(3L, 3L))
  expect_identical(details$sample_summary$excluded_samples, c(0L, 0L))
  expect_identical(nrow(details$retained_profiles), 6L)
  expect_identical(nrow(details$excluded_profiles), 0L)
  expect_identical(details$aggregation$layer, "counts")
  reverse <- suppressWarnings(run_dge(pb, test_use = "edgeR", group_1 = "group_2", group_2 = "group_1"))
  expect_equal(result$log2FC, -reverse$log2FC[match(result$feature, reverse$feature)], tolerance = 1e-8)
})

test_that("the edgeR sequence is filter, then TMM, then glmQLFit (no estimateDisp)", {
  skip_if_no_edger()
  fixture <- edger_fixture(FALSE)
  pb <- fixture_pseudobulk(fixture)
  result <- suppressWarnings(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"))
  # the same pipeline by hand
  group <- factor(pb$profiles$group, levels = c("group_2", "group_1"))
  y <- edgeR::DGEList(pb$counts, group = group)
  keep <- edgeR::filterByExpr(y, group = group)
  y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  design <- stats::model.matrix(~ 0 + group)
  fit <- edgeR::glmQLFit(y, design)
  tab <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast = c(-1, 1)), n = Inf, sort.by = "none")$table
  expect_equal(result$pval[match(rownames(tab), result$feature)], tab$PValue, tolerance = 1e-10)
  expect_equal(result$log2FC[match(rownames(tab), result$feature)], tab$logFC, tolerance = 1e-10)
  expect_equal(attr(result, "edger_details")$normalization_factors, stats::setNames(y$samples$norm.factors, rownames(y$samples)))
  expect_identical(attr(result, "edger_details")$genes_tested, sum(keep))
})

test_that("non-integer pseudobulk counts and normalised layer names are rejected", {
  skip_if_no_edger()
  pb <- fixture_pseudobulk(edger_fixture(FALSE))
  pb$counts[1, 1] <- 0.5
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"), "non-integer")
  pb <- fixture_pseudobulk(edger_fixture(FALSE)); pb$fun <- "mean"
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"), "summed counts")
  object <- as_sce(edger_fixture(FALSE))
  expect_error(run_dge(object, "condition", layer = "logcounts", test_use = "edgeR", sample_by = "sample_id",
                       group_1 = "group_1", group_2 = "group_2"), "normalised values")
})

test_that("argument validation", {
  skip_if_no_edger()
  pb <- fixture_pseudobulk(edger_fixture(FALSE))
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_1"), "must be different")
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1"), "group_1. and .group_2")
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2", min_cells = 0), "min_cells")
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2", lfc_format = "ln"), "log2")
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2", robust = NA), "robust")
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "nope"), "not found")
  expect_error(run_dge(pb, test_use = "Presto", group_1 = "group_1", group_2 = "group_2"), "edgeR")
  object <- as_sce(edger_fixture(FALSE))
  expect_error(run_dge(object, "condition", test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"), "sample_by")
  expect_error(run_dge(object, "condition", test_use = "edgeR", sample_by = "sample_id", group_1 = "group_1",
                       group_2 = "nope"), "not found")
})

test_that("min_cells excludes profiles and the details say why", {
  skip_if_no_edger()
  fixture <- edger_fixture(FALSE)
  keep <- !(fixture$metadata$sample_id == "sample_1" & seq_len(nrow(fixture$metadata)) %% 3L != 0L)
  fixture$counts <- fixture$counts[, keep]; fixture$metadata <- fixture$metadata[keep, ]
  pb <- fixture_pseudobulk(fixture)
  expect_lt(pb$profiles$n_cells[pb$profiles$sample == "sample_1"], 10L)
  result <- suppressWarnings(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2", min_cells = 10L))
  details <- attr(result, "edger_details")
  expect_identical(details$excluded_profiles$sample, "sample_1")
  expect_identical(details$excluded_profiles$exclusion_reason, "fewer than 10 cells")
  expect_identical(details$sample_summary$retained_samples, c(3L, 2L))
  expect_identical(details$sample_summary$excluded_samples, c(0L, 1L))
  expect_warning(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2", min_cells = 10L), "excluded")
  expect_warning(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2", min_cells = 10L), "Fewer than three")
})

test_that("filtering options preserve the common schema", {
  skip_if_no_edger()
  pb <- fixture_pseudobulk(edger_fixture(FALSE))
  observed <- suppressWarnings(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2",
                                       positive_only = TRUE, remove_raw_pval = TRUE))
  expect_false("pval" %in% colnames(observed))
  expect_true(all(observed$log2FC > 0))
  expect_identical(names(observed)[1:3], c("method", "group", "feature"))
  expect_false(is.unsorted(observed$pval_adj))
})

test_that("per-cell sample identifiers are rejected before any aggregation", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  fixture <- edger_fixture(FALSE)
  fixture$metadata$cell_id <- rownames(fixture$metadata)
  object <- as_sce(fixture)
  called <- FALSE
  testthat::local_mocked_bindings(pseudobulk = function(...) { called <<- TRUE; stop("should not aggregate") }, .package = "SCUBA")
  expect_error(run_dge(object, "condition", test_use = "edgeR", sample_by = "cell_id", group_1 = "group_1",
                       group_2 = "group_2", min_cells = 1L), "cell-level")
  expect_false(called)
})
