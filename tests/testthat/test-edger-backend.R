test_that("edgeR accepts sample-level designs and rejects cluster designs", {
  ordinary <- edger_fixture(FALSE)
  ordinary_result <- scDE:::.scde_run_edger(
    ordinary$counts, ordinary$metadata, "sample_id", "condition",
    "group_1", "group_2", 10L, FALSE, FALSE, FALSE
  )
  expect_edger_result(ordinary_result, "ordinary")
  expect_identical(dim(attr(ordinary_result, "edger_details")$design_matrix), c(6L, 2L))

  blocked <- edger_fixture(TRUE)
  expect_error(
    scDE:::.scde_run_edger(
      blocked$counts, blocked$metadata, "sample_id", "cluster",
      "group_1", "group_2", 10L, FALSE, FALSE, FALSE
    ),
    "cluster-to-cluster"
  )
})

test_that("group_1 is the numerator", {
  fixture <- edger_fixture(FALSE)
  forward <- scDE:::.scde_run_edger(
    fixture$counts, fixture$metadata, "sample_id", "condition",
    "group_1", "group_2", 10L, FALSE, FALSE, FALSE
  )
  reverse <- scDE:::.scde_run_edger(
    fixture$counts, fixture$metadata, "sample_id", "condition",
    "group_2", "group_1", 10L, FALSE, FALSE, FALSE
  )
  expect_equal(
    forward$log2FC,
    -reverse$log2FC[match(forward$feature, reverse$feature)],
    tolerance = 1e-8
  )
})

test_that("raw count and replication validation is explicit", {
  fixture <- edger_fixture(FALSE)
  fixture$counts[1, 1] <- 0.5
  expect_error(
    scDE:::.scde_run_edger(
      fixture$counts, fixture$metadata, "sample_id", "condition",
      "group_1", "group_2", 10L, FALSE, FALSE, FALSE
    ),
    "non-integer"
  )

  fixture <- edger_fixture(FALSE)
  keep <- fixture$metadata$sample_id %in% c("sample_1", "sample_4")
  expect_error(
    scDE:::.scde_run_edger(
      fixture$counts[, keep], fixture$metadata[keep, ], "sample_id", "condition",
      "group_1", "group_2", 10L, FALSE, FALSE, FALSE
    ),
    "at least two"
  )
})

test_that("pseudobulk Wilcoxon retains only adequately represented pairs", {
  fixture <- edger_fixture(TRUE)
  drop <- fixture$metadata$sample_id == "sample_3" & fixture$metadata$cluster == "group_1"
  observed <- suppressWarnings(scDE:::.scde_run_edger(
    fixture$counts[, !drop], fixture$metadata[!drop, ], "sample_id", "cluster",
    "group_1", "group_2", 10L, FALSE, FALSE, FALSE, method = "wilcox"
  ))
  details <- attr(observed, "pseudobulk_details")
  expect_identical(unique(details$retained_profiles$sample_id), c("sample_1", "sample_2"))
  expect_true(all(details$sample_summary$retained_samples == 2L))
  expect_true(any(details$excluded_profiles$sample_id == "sample_3"))
})

test_that("edgeR filtering options preserve the common schema", {
  fixture <- edger_fixture(FALSE)
  observed <- scDE:::.scde_run_edger(
    fixture$counts, fixture$metadata, "sample_id", "condition",
    "group_1", "group_2", 10L, TRUE, TRUE, FALSE
  )
  expect_false("pval" %in% colnames(observed))
  expect_true(all(observed$log2FC > 0))
})
