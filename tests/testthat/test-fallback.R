# Automatic per-contrast fallback to the cell-level Wilcoxon test when edgeR cannot run.

fallback_fixture <- function() {
  # one sample per group: never replicated
  fixture <- edger_fixture(FALSE)
  keep <- fixture$metadata$sample_id %in% c("sample_1", "sample_4")
  fixture$counts <- fixture$counts[, keep]
  fixture$metadata <- fixture$metadata[keep, ]
  fixture
}

test_that("an unreplicated single contrast returns the cell-level Wilcoxon result, labelled", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  object <- as_sce(fallback_fixture())
  expect_message(
    result <- run_dge(object, "condition", layer = "counts", test_use = "edgeR", sample_by = "sample_id",
                      group_1 = "group_1", group_2 = "group_2"),
    "Falling back"
  )
  expect_false(inherits(result, "scDE_edger_results"))
  expect_true(all(result$method == "wilcoxon"))
  expect_true(all(result$group == "group_1"))
  expect_identical(names(result)[1:3], c("method", "group", "feature"))
  expect_true(all(c("auc", "pct_in", "pct_out") %in% names(result)))
  details <- attr(result, "edger_details")
  expect_identical(details$method, "wilcoxon")
  expect_true(is.na(details$design_type))
  expect_match(details$fallback_reason, "has 1 sample")
  expect_identical(details$sample_summary$retained_samples, c(1L, 1L))
  # equal to the cell-level test on the two groups
  reference <- run_dge(object, "condition")
  reference <- reference[reference$group == "group_1", ]
  expect_equal(sorted_result(result)[, c("feature", "log2FC", "pval", "pval_adj")],
               sorted_result(reference)[, c("feature", "log2FC", "pval", "pval_adj")], tolerance = 1e-10)
})

test_that("the fallback never aggregates, and honours the filters", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  object <- as_sce(fallback_fixture())
  called <- FALSE
  testthat::local_mocked_bindings(pseudobulk = function(...) { called <<- TRUE; stop("should not aggregate") }, .package = "SCUBA")
  result <- suppressMessages(run_dge(object, "condition", layer = "counts", test_use = "edgeR", sample_by = "sample_id",
                                     group_1 = "group_1", group_2 = "group_2", positive_only = TRUE, remove_raw_pval = TRUE))
  expect_false(called)
  expect_true(all(result$log2FC > 0))
  expect_false("pval" %in% names(result))
})

test_that("the fallback restricts the comparison to the two groups of the contrast", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  data <- multicontrast_fixture_data()
  # keep only one sample so no contrast is replicated
  keep <- data$metadata$sample_id == "sample1"
  data$counts <- data$counts[, keep]; data$metadata <- data$metadata[keep, ]
  object <- as_sce(data)
  result <- suppressMessages(run_dge(object, "cluster", layer = "counts", test_use = "edgeR", sample_by = "sample_id",
                                     group_1 = "A", group_2 = "B"))
  two <- object[, object$cluster %in% c("A", "B")]
  reference <- run_dge(two, "cluster")
  reference <- reference[reference$group == "A", ]
  expect_equal(sorted_result(result)[, c("feature", "pval")], sorted_result(reference)[, c("feature", "pval")], tolerance = 1e-10)
  # gene4:6 are up in B, so down in A vs B; against the rest (B and C) they are still down
  expect_true(all(result$log2FC[match(paste0("gene", 4:6), result$feature)] < 0))
})

test_that("pre-aggregated objects cannot fall back", {
  skip_if_no_edger()
  pb <- fixture_pseudobulk(fallback_fixture())
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"), "no cells")
})

test_that("cell-level backends warn when sample_by is ignored", {
  fixture <- edger_fixture(FALSE)
  object <- as_sce(fixture)
  expect_warning(
    result <- run_dge(object, group_by = "condition", test_use = "Presto", sample_by = "sample_id"),
    "sample_by.*ignored"
  )
  expect_false(inherits(result, "scDE_edger_results"))
  expect_false("method" %in% names(result))
})
