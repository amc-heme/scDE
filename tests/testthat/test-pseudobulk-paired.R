# Paired designs: samples contribute cells to both groups (cluster vs cluster), fitted as
# ~ sample + group.

test_that("a paired design is fitted with a sample block and the group coefficient", {
  skip_if_no_edger()
  pb <- fixture_pseudobulk(edger_fixture(TRUE), group_by = "cluster")
  result <- suppressWarnings(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"))
  expect_edger_result(result, "paired")
  details <- attr(result, "edger_details")
  expect_identical(dim(details$design_matrix), c(6L, 4L))
  expect_identical(colnames(details$design_matrix), c("intercept", "sample_sample_2", "sample_sample_3", "group_1_vs_group_2"))
  expect_identical(details$coef, 4L)
  expect_null(details$contrast)
  expect_identical(details$sample_summary$retained_samples, c(3L, 3L))
  # by hand
  group <- factor(pb$profiles$group, levels = c("group_2", "group_1"))
  y <- edgeR::DGEList(pb$counts, group = group)
  keep <- edgeR::filterByExpr(y, group = group)
  y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  design <- stats::model.matrix(~ factor(pb$profiles$sample) + group)
  tab <- edgeR::topTags(edgeR::glmQLFTest(edgeR::glmQLFit(y, design), coef = 4), n = Inf, sort.by = "none")$table
  expect_equal(result$pval[match(rownames(tab), result$feature)], tab$PValue, tolerance = 1e-10)
  expect_equal(result$log2FC[match(rownames(tab), result$feature)], tab$logFC, tolerance = 1e-10)
})

test_that("a sample missing one side after min_cells is dropped; two pairs still run, one does not", {
  skip_if_no_edger()
  fixture <- edger_fixture(TRUE)
  drop <- fixture$metadata$sample_id == "sample_3" & fixture$metadata$cluster == "group_1"
  fixture$counts <- fixture$counts[, !drop]; fixture$metadata <- fixture$metadata[!drop, ]
  pb <- fixture_pseudobulk(fixture, group_by = "cluster")
  result <- suppressWarnings(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"))
  details <- attr(result, "edger_details")
  expect_identical(details$design_type, "paired")
  expect_setequal(unique(details$retained_profiles$sample), c("sample_1", "sample_2"))
  expect_identical(details$excluded_profiles$sample, "sample_3")
  expect_identical(details$excluded_profiles$exclusion_reason, "paired group profile unavailable")
  expect_identical(details$sample_summary$retained_samples, c(2L, 2L))
  expect_warning(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"), "Fewer than three")
  # one pair: not replicated (pre-aggregated objects cannot fall back)
  drop2 <- fixture$metadata$sample_id == "sample_2" & fixture$metadata$cluster == "group_2"
  fixture$counts <- fixture$counts[, !drop2]; fixture$metadata <- fixture$metadata[!drop2, ]
  pb <- fixture_pseudobulk(fixture, group_by = "cluster")
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2"), "Only 1 sample has")
})

test_that("the design type is decided after min_cells filtering", {
  skip_if_no_edger()
  # sample_1 spans both groups, but its group_1 profile is tiny: an independent design remains
  fixture <- edger_fixture(FALSE)
  stray <- which(fixture$metadata$sample_id == "sample_1")[1:2]
  fixture$metadata$condition[stray] <- "group_1"
  pb <- fixture_pseudobulk(fixture)
  expect_true(any(pb$profiles$sample == "sample_1" & pb$profiles$group == "group_1"))
  result <- suppressWarnings(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2", min_cells = 5L))
  expect_identical(attr(result, "edger_details")$design_type, "independent")
  expect_identical(attr(result, "edger_details")$excluded_profiles$exclusion_reason, "fewer than 5 cells")
  # with min_cells = 1 the same data is a paired design with one pair, plus independent samples
  # on both sides: sample_1 is the only paired sample, so edgeR needs 2 pairs and cannot run
  expect_error(run_dge(pb, test_use = "edgeR", group_1 = "group_1", group_2 = "group_2", min_cells = 1L), "paired samples")
})

test_that("multicontrast fixture: every cluster contrast is a paired edgeR test", {
  skip_if_no_edger()
  skip_if_no_scuba_pseudobulk()
  object <- as_sce(multicontrast_fixture_data())
  observed <- suppressWarnings(run_dge_multicontrast(
    object, group_by = "cluster", sample_by = "sample_id", contrast_mode = "pairwise", min_cells = 10L
  ))
  expect_s3_class(observed, "scDE_edger_multicontrast_results")
  expect_setequal(unique(observed$contrast), c("A vs B", "A vs C", "B vs C"))
  expect_true(all(observed$method == "edgeR"))
  details <- attr(observed, "edger_multicontrast_details")
  expect_true(all(details$contrast_table$method == "edgeR"))
  expect_true(all(details$contrast_table$design_type == "paired"))
  expect_true(all(vapply(details$contrasts, function(x) identical(x$design_type, "paired"), logical(1))))
  a_vs_b <- observed[observed$contrast == "A vs B", ]
  expect_true(all(a_vs_b$log2FC[match(paste0("gene", 1:3), a_vs_b$feature)] > 0))
  expect_true(all(a_vs_b$log2FC[match(paste0("gene", 4:6), a_vs_b$feature)] < 0))
})
