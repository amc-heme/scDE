edger_fixture <- function(blocked = FALSE) {
  set.seed(20260924)
  n_genes <- 30L
  cells_per_profile <- 12L
  if (blocked) {
    sample_id <- rep(paste0("sample_", 1:3), each = 2L * cells_per_profile)
    group <- rep(rep(c("group_2", "group_1"), each = cells_per_profile), 3L)
  } else {
    sample_id <- rep(paste0("sample_", 1:6), each = cells_per_profile)
    group <- rep(rep(c("group_2", "group_1"), each = 3L), each = cells_per_profile)
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
    condition = group,
    cluster = group,
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  list(counts = counts, metadata = metadata)
}

run_fixture_edger <- function(object, group_by = "condition", layer = "counts", ...) {
  suppressWarnings(
    run_dge(
      object,
      group_by = group_by,
      layer = layer,
      test_use = "edgeR",
      sample_by = "sample_id",
      group_1 = "group_1",
      group_2 = "group_2",
      min_cells = 10L,
      ...
    )
  )
}

expect_edger_result <- function(result, design_type) {
  expect_s3_class(result, "scDE_edger_results")
  expect_true(all(c("group", "feature", "avgExpr", "log2FC", "pval", "pval_adj") %in%
                    colnames(result)))
  expect_true(all(result$group == "group_1"))
  expected_features <- paste0("gene", 1:3)
  expect_true(all(result$log2FC[match(expected_features, result$feature)] > 0))
  expect_identical(attr(result, "edger_details")$design_type, design_type)
  expect_identical(attr(result, "edger_details")$normalization_method, "TMM")
}
