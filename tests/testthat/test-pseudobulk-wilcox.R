wilcox_fixture <- function(n_samples_per_group = 6L, blocked = FALSE) {
  set.seed(20261008)
  n_genes <- 30L
  cells_per_profile <- 12L
  if (blocked) {
    samples <- paste0("sample_", seq_len(n_samples_per_group))
    sample_id <- rep(samples, each = 2L * cells_per_profile)
    group <- rep(rep(c("group_2", "group_1"), each = cells_per_profile),
                 n_samples_per_group)
  } else {
    samples <- paste0("sample_", seq_len(2L * n_samples_per_group))
    sample_id <- rep(samples, each = cells_per_profile)
    group <- rep(rep(c("group_2", "group_1"), each = n_samples_per_group),
                 each = cells_per_profile)
  }
  counts <- matrix(
    stats::rpois(n_genes * length(sample_id), lambda = 3),
    nrow = n_genes,
    dimnames = list(
      paste0("gene", seq_len(n_genes)),
      paste0("cell_", seq_along(sample_id))
    )
  )
  counts[1:3, group == "group_1"] <- counts[1:3, group == "group_1"] + 8L
  metadata <- data.frame(
    sample_id = sample_id,
    cell_id = colnames(counts),
    condition = group,
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  list(counts = counts, metadata = metadata)
}

run_wilcox <- function(fixture, sample_by = "sample_id", min_cells = 10L) {
  scDE:::.scde_run_edger(
    fixture$counts, fixture$metadata, sample_by, "condition",
    "group_1", "group_2", min_cells, FALSE, FALSE, FALSE,
    method = "wilcox"
  )
}

test_that("pseudobulk Wilcoxon detects shifted genes in ordinary designs", {
  result <- run_wilcox(wilcox_fixture(6L))
  expect_s3_class(result, "scDE_pseudobulk_wilcox_results")
  expect_true(all(c("group", "feature", "avgExpr", "log2FC", "pval", "pval_adj") %in%
                    colnames(result)))
  top <- result[match(paste0("gene", 1:3), result$feature), ]
  expect_true(all(top$log2FC > 0))
  expect_true(all(top$pval < 0.01))
  details <- attr(result, "pseudobulk_details")
  expect_identical(details$method, "wilcox")
  expect_identical(details$design_type, "ordinary")
})

test_that("pseudobulk Wilcoxon uses signed-rank tests for blocked designs", {
  result <- run_wilcox(wilcox_fixture(6L, blocked = TRUE))
  expect_identical(attr(result, "pseudobulk_details")$design_type, "blocked")
  top <- result[match(paste0("gene", 1:3), result$feature), ]
  expect_true(all(top$log2FC > 0))
  expect_true(all(top$pval < 0.05))
})

test_that("pseudobulk Wilcoxon warns when samples are too few for significance", {
  warnings <- character()
  withCallingHandlers(
    run_wilcox(wilcox_fixture(3L)),
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("smallest attainable Wilcoxon p-value", warnings)))
})

test_that("pseudobulk tests reject per-cell sample identifiers", {
  fixture <- wilcox_fixture(6L)
  expect_error(run_wilcox(fixture, "cell_id", min_cells = 1L), "cell-level")
  expect_error(
    scDE:::.scde_run_edger(
      fixture$counts, fixture$metadata, "cell_id", "condition",
      "group_1", "group_2", 1L, FALSE, FALSE, FALSE
    ),
    "cell-level"
  )
})

test_that("cell-level pseudobulk is rejected before raw-count validation", {
  fixture <- wilcox_fixture(3L)
  fixture$counts[1, 1] <- 0.5

  expect_error(
    scDE:::.scde_run_edger(
      fixture$counts, fixture$metadata, "cell_id", "condition",
      "group_1", "group_2", 1L, FALSE, FALSE, FALSE
    ),
    "cell-level"
  )
})

test_that("run_dge and run_dge_multicontrast dispatch pseudobulk_wilcox", {
  skip_if_not_installed("SingleCellExperiment")
  fixture <- wilcox_fixture(6L)
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(fixture$counts, sparse = TRUE)),
    colData = fixture$metadata
  )
  result <- run_dge(
    object, "condition", layer = "counts", test_use = "pseudobulk_wilcox",
    sample_by = "sample_id", group_1 = "group_1", group_2 = "group_2"
  )
  expect_s3_class(result, "scDE_pseudobulk_wilcox_results")

  multi <- run_dge_multicontrast(
    object, "condition", "sample_id", contrast_mode = "pairwise",
    layer = "counts", test_use = "pseudobulk_wilcox"
  )
  expect_s3_class(multi, "scDE_pseudobulk_wilcox_multicontrast_results")
  expect_identical(attr(multi, "pseudobulk_multicontrast_details")$method, "wilcox")
  expect_error(
    run_dge_multicontrast(object, "condition", "sample_id", test_use = "Presto"),
    "must be"
  )
})

test_that("cell-level backends warn when sample_by is ignored", {
  skip_if_not_installed("SingleCellExperiment")
  fixture <- wilcox_fixture(3L)
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      logcounts = log1p(Matrix::Matrix(fixture$counts, sparse = TRUE))
    ),
    colData = fixture$metadata
  )

  expect_warning(
    result <- run_dge(
      object,
      group_by = "condition",
      test_use = "Presto",
      sample_by = "sample_id"
    ),
    "sample_by.*ignored"
  )
  expect_false(inherits(result, "scDE_pseudobulk_wilcox_results"))
})
