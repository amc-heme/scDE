multicontrast_fixture <- function() {
  set.seed(20260925)
  clusters <- c("A", "B", "C")
  samples <- paste0("sample", 1:3)
  cells_per_profile <- 12L
  profile <- expand.grid(
    cluster = clusters,
    sample_id = samples,
    stringsAsFactors = FALSE
  )
  metadata <- profile[rep(seq_len(nrow(profile)), each = cells_per_profile), ]
  rownames(metadata) <- paste0("cell", seq_len(nrow(metadata)))
  counts <- matrix(
    stats::rpois(40L * nrow(metadata), lambda = 3),
    nrow = 40L,
    dimnames = list(paste0("gene", 1:40), rownames(metadata))
  )
  counts[1:3, metadata$cluster == "A"] <-
    counts[1:3, metadata$cluster == "A"] + 10L
  counts[4:6, metadata$cluster == "B"] <-
    counts[4:6, metadata$cluster == "B"] + 8L
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(counts, sparse = TRUE)),
    colData = metadata
  )
}

test_that("pseudobulk Wilcoxon runs every paired cluster contrast", {
  skip_if_not_installed("SingleCellExperiment")
  object <- multicontrast_fixture()
  observed <- suppressWarnings(run_dge_multicontrast(
    object,
    group_by = "cluster",
    sample_by = "sample_id",
    contrast_mode = "pairwise",
    min_cells = 10L,
    test_use = "pseudobulk_wilcox"
  ))

  expect_s3_class(observed, "scDE_pseudobulk_wilcox_multicontrast_results")
  expect_setequal(unique(observed$contrast), c("A vs B", "A vs C", "B vs C"))
  expect_true(all(c("group_1", "group_2", "contrast", "feature", "log2FC") %in%
                    colnames(observed)))
  details <- attr(observed, "pseudobulk_multicontrast_details")
  expect_identical(details$contrast_mode, "pairwise")
  expect_true(all(vapply(details$contrasts, function(x) {
    identical(x$design_type, "blocked")
  }, logical(1))))
  a_vs_b <- observed[observed$contrast == "A vs B", ]
  expect_true(all(a_vs_b$log2FC[match(paste0("gene", 1:3), a_vs_b$feature)] > 0))
})

test_that("reference mode uses the reference as denominator", {
  skip_if_not_installed("SingleCellExperiment")
  observed <- suppressWarnings(run_dge_multicontrast(
    multicontrast_fixture(),
    group_by = "cluster",
    sample_by = "sample_id",
    contrast_mode = "reference",
    reference_group = "C",
    test_use = "pseudobulk_wilcox"
  ))

  expect_setequal(unique(observed$contrast), c("A vs C", "B vs C"))
  expect_true(all(observed$group_2 == "C"))
})

test_that("one-versus-rest pools non-target clusters within sample", {
  skip_if_not_installed("SingleCellExperiment")
  observed <- suppressWarnings(run_dge_multicontrast(
    multicontrast_fixture(),
    group_by = "cluster",
    sample_by = "sample_id",
    contrast_mode = "one_vs_rest",
    groups = c("A", "B"),
    test_use = "pseudobulk_wilcox"
  ))

  expect_setequal(unique(observed$contrast), c("A vs rest", "B vs rest"))
  expect_true(all(observed$group_2 == "rest"))
  details <- attr(observed, "pseudobulk_multicontrast_details")$contrasts[["A vs rest"]]
  expect_identical(details$design_type, "blocked")
  expect_true(all(details$sample_summary$retained_samples == 3L))
  expect_identical(
    sort(unique(details$retained_profiles$group)),
    c("A", "rest")
  )
})

test_that("global adjustment modes are explicit", {
  skip_if_not_installed("SingleCellExperiment")
  both <- suppressWarnings(run_dge_multicontrast(
    multicontrast_fixture(), "cluster", "sample_id",
    contrast_mode = "pairwise", p_adjust_scope = "both",
    test_use = "pseudobulk_wilcox"
  ))
  expect_true("pval_adj_global" %in% colnames(both))
  expect_equal(both$pval_adj_global, stats::p.adjust(both$pval, "BH"))

  global <- suppressWarnings(run_dge_multicontrast(
    multicontrast_fixture(), "cluster", "sample_id",
    contrast_mode = "reference", reference_group = "C",
    p_adjust_scope = "global", remove_raw_pval = TRUE,
    test_use = "pseudobulk_wilcox"
  ))
  expect_false("pval" %in% colnames(global))
  expect_false("pval_adj_global" %in% colnames(global))
})

test_that("edgeR multicontrast rejects paired cluster comparisons", {
  skip_if_not_installed("SingleCellExperiment")
  expect_error(
    run_dge_multicontrast(
      multicontrast_fixture(), "cluster", "sample_id",
      contrast_mode = "pairwise"
    ),
    "cluster-to-cluster"
  )
})

test_that("multicontrast validation identifies invalid requests and contrasts", {
  skip_if_not_installed("SingleCellExperiment")
  object <- multicontrast_fixture()
  expect_error(
    run_dge_multicontrast(object, "cluster", "sample_id", groups = "A"),
    "at least two groups"
  )
  expect_error(
    run_dge_multicontrast(
      object, "cluster", "sample_id", contrast_mode = "reference"
    ),
    "reference_group"
  )
  expect_error(
    run_dge_multicontrast(object, "cluster", "sample_id", groups = c("A", "D")),
    "not found.*D"
  )
})
