# run_dge_multicontrast(): one aggregation, many contrasts, per-contrast fallback.

test_that("reference mode uses the reference as denominator", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  observed <- suppressWarnings(run_dge_multicontrast(
    as_sce(multicontrast_fixture_data()), group_by = "cluster", sample_by = "sample_id",
    contrast_mode = "reference", reference_group = "C"
  ))
  expect_setequal(unique(observed$contrast), c("A vs C", "B vs C"))
  expect_true(all(observed$group_2 == "C"))
  expect_true(all(observed$method == "edgeR"))
  expect_identical(names(observed)[1:5], c("group_1", "group_2", "contrast", "method", "group"))
})

test_that("one-versus-rest pools the other clusters within each sample by column sums", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  data <- multicontrast_fixture_data()
  object <- as_sce(data)
  observed <- suppressWarnings(run_dge_multicontrast(
    object, group_by = "cluster", sample_by = "sample_id", contrast_mode = "one_vs_rest", groups = c("A", "B")
  ))
  expect_setequal(unique(observed$contrast), c("A vs rest", "B vs rest"))
  expect_true(all(observed$group_2 == "rest"))
  details <- attr(observed, "edger_multicontrast_details")$contrasts[["A vs rest"]]
  expect_identical(details$design_type, "paired")
  expect_true(all(details$sample_summary$retained_samples == 3L))
  expect_setequal(unique(details$retained_profiles$group), c("A", "rest"))
  expect_identical(details$retained_profiles$n_cells[details$retained_profiles$group == "rest"], rep(24L, 3))
  # the pooled rest equals relabelling the cells and running the single contrast
  data$metadata$relabel <- ifelse(data$metadata$cluster == "A", "A", "rest")
  single <- suppressWarnings(run_dge(as_sce(data), "relabel", layer = "counts", test_use = "edgeR",
                                     sample_by = "sample_id", group_1 = "A", group_2 = "rest"))
  a <- observed[observed$contrast == "A vs rest", ]
  expect_equal(sorted_result(a), sorted_result(single), tolerance = 1e-10)
})

test_that("SCUBA::pseudobulk is called once per call, and not at all when every contrast falls back", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  object <- as_sce(multicontrast_fixture_data())
  calls <- 0L
  real <- SCUBA::pseudobulk
  testthat::local_mocked_bindings(pseudobulk = function(...) { calls <<- calls + 1L; real(...) }, .package = "SCUBA")
  suppressWarnings(run_dge_multicontrast(object, "cluster", "sample_id", contrast_mode = "pairwise"))
  expect_identical(calls, 1L)
  calls <- 0L
  one_sample <- object[, object$sample_id == "sample1"]
  suppressMessages(run_dge_multicontrast(one_sample, "cluster", "sample_id", contrast_mode = "pairwise"))
  expect_identical(calls, 0L)
})

test_that("a mixed run labels edgeR and Wilcoxon rows and NA-fills the Wilcoxon-only columns", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  data <- multicontrast_fixture_data()
  # cluster C only in sample1: A vs B paired across 3 samples (edgeR); anything vs C unreplicated
  drop <- data$metadata$cluster == "C" & data$metadata$sample_id != "sample1"
  data$counts <- data$counts[, !drop]; data$metadata <- data$metadata[!drop, ]
  object <- as_sce(data)
  expect_message(
    observed <- suppressWarnings(run_dge_multicontrast(object, "cluster", "sample_id", contrast_mode = "pairwise")),
    "Cell-level Wilcoxon was used for 2 contrast"
  )
  table <- attr(observed, "edger_multicontrast_details")$contrast_table
  expect_identical(table$method[table$contrast == "A vs B"], "edgeR")
  expect_identical(table$method[table$contrast %in% c("A vs C", "B vs C")], c("wilcoxon", "wilcoxon"))
  expect_match(table$reason[table$contrast == "A vs C"], "Only 1 sample")
  expect_true(all(observed$method[observed$contrast == "A vs B"] == "edgeR"))
  expect_true(all(observed$method[observed$contrast == "A vs C"] == "wilcoxon"))
  expect_true(all(is.na(observed$auc[observed$method == "edgeR"])))
  expect_false(anyNA(observed$auc[observed$method == "wilcoxon"]))
  expect_identical(attr(observed, "edger_multicontrast_details")$contrasts[["A vs C"]]$method, "wilcoxon")
  # one_vs_rest fallbacks share one cell-level run and equal run_dge() on the object
  ovr <- suppressMessages(suppressWarnings(run_dge_multicontrast(object, "cluster", "sample_id", contrast_mode = "one_vs_rest", min_cells = 10L)))
  reference <- run_dge(object, "cluster")
  c_rows <- ovr[ovr$contrast == "C vs rest", ]
  expect_true(all(c_rows$method == "wilcoxon"))
  expect_equal(sorted_result(c_rows)[, c("feature", "pval")],
               sorted_result(reference[reference$group == "C", ])[, c("feature", "pval")], tolerance = 1e-10)
})

test_that("global adjustment runs before positive_only filtering", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  object <- as_sce(multicontrast_fixture_data())
  both <- suppressWarnings(run_dge_multicontrast(object, "cluster", "sample_id", contrast_mode = "pairwise", p_adjust_scope = "both"))
  expect_true("pval_adj_global" %in% colnames(both))
  expect_equal(both$pval_adj_global, stats::p.adjust(both$pval, "BH"))
  positive <- suppressWarnings(run_dge_multicontrast(object, "cluster", "sample_id", contrast_mode = "pairwise",
                                                     p_adjust_scope = "global", positive_only = TRUE))
  expect_true(all(positive$log2FC > 0))
  key <- function(x) paste(x$contrast, x$feature)
  expect_equal(positive$pval_adj, both$pval_adj_global[match(key(positive), key(both))])
  global <- suppressWarnings(run_dge_multicontrast(object, "cluster", "sample_id", contrast_mode = "reference",
                                                   reference_group = "C", p_adjust_scope = "global", remove_raw_pval = TRUE))
  expect_false("pval" %in% colnames(global))
  expect_false("pval_adj_global" %in% colnames(global))
})

test_that("pre-aggregated input works when every contrast is replicated, and errors otherwise", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  data <- multicontrast_fixture_data()
  pb <- SCUBA::pseudobulk(as_sce(data), "sample_id", "cluster")
  from_pb <- suppressWarnings(run_dge_multicontrast(pb, contrast_mode = "pairwise"))
  from_object <- suppressWarnings(run_dge_multicontrast(as_sce(data), "cluster", "sample_id", contrast_mode = "pairwise"))
  # identical tables; the details differ only in the recorded source class
  expect_equal(as.data.frame(from_pb), as.data.frame(from_object), tolerance = 1e-10, ignore_attr = TRUE)
  expect_true(is.na(attr(from_pb, "edger_multicontrast_details")$aggregation$source_class))
  drop <- data$metadata$cluster == "C" & data$metadata$sample_id != "sample1"
  data$counts <- data$counts[, !drop]; data$metadata <- data$metadata[!drop, ]
  pb2 <- SCUBA::pseudobulk(as_sce(data), "sample_id", "cluster")
  expect_error(run_dge_multicontrast(pb2, contrast_mode = "pairwise"), "no cells")
})

test_that("multicontrast validation identifies invalid requests", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  object <- as_sce(multicontrast_fixture_data())
  expect_error(run_dge_multicontrast(object, "cluster", "sample_id", groups = "A"), "at least two groups")
  expect_error(run_dge_multicontrast(object, "cluster", "sample_id", contrast_mode = "reference"), "reference_group")
  expect_error(run_dge_multicontrast(object, "cluster", "sample_id", groups = c("A", "D")), "not found.*D")
  expect_error(run_dge_multicontrast(object, "cluster", "sample_id", test_use = "Presto"), "must be")
  expect_error(run_dge_multicontrast(object, "cluster", "sample_id", layer = "logcounts"), "normalised")
})
