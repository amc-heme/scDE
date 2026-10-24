# check_pseudobulk_design(): the replication rules, from profiles alone.

design_profiles <- function(...) {
  # list of sample = c(group = n_cells, ...)
  spec <- list(...)
  rows <- do.call(rbind, lapply(names(spec), function(s) {
    data.frame(sample = s, group = names(spec[[s]]), n_cells = as.integer(spec[[s]]), stringsAsFactors = FALSE)
  }))
  rows$profile_id <- paste(rows$sample, rows$group, sep = "__")
  rows
}

test_that("independent and paired designs are classified and replicated correctly", {
  p <- design_profiles(s1 = c(A = 20), s2 = c(A = 20), s3 = c(B = 20), s4 = c(B = 20))
  v <- check_pseudobulk_design(p, data.frame(group_1 = "A", group_2 = "B"))
  expect_identical(v$contrasts$method, "edgeR")
  expect_identical(v$contrasts$design_type, "independent")
  expect_identical(v$contrasts$samples_group_1, 2L)
  expect_identical(v$contrasts$paired_samples, 0L)
  expect_true(is.na(v$contrasts$reason))
  expect_true(v$any_edger); expect_false(v$any_fallback)
  expect_match(v$warnings, "Fewer than three")
  p <- design_profiles(s1 = c(A = 20, B = 20), s2 = c(A = 20, B = 20), s3 = c(A = 20, B = 20))
  v <- check_pseudobulk_design(p, data.frame(group_1 = "A", group_2 = "B"))
  expect_identical(v$contrasts$design_type, "paired")
  expect_identical(v$contrasts$paired_samples, 3L)
  expect_identical(v$contrasts$method, "edgeR")
  expect_length(v$warnings, 0L)
})

test_that("unreplicated contrasts fall back with an actionable reason", {
  p <- design_profiles(s1 = c(A = 20), s2 = c(B = 20))
  v <- check_pseudobulk_design(p, data.frame(group_1 = "A", group_2 = "B"))
  expect_identical(v$contrasts$method, "wilcoxon")
  expect_match(v$contrasts$reason, "'A' has 1 sample and 'B' has 1 sample")
  expect_match(v$contrasts$reason, "edgeR needs 2 samples per group with at least 10 cells")
  expect_true(v$any_fallback); expect_false(v$any_edger)
  p <- design_profiles(s1 = c(A = 20, B = 20), s2 = c(A = 20, B = 3))
  v <- check_pseudobulk_design(p, data.frame(group_1 = "A", group_2 = "B"), min_cells = 10L)
  expect_identical(v$contrasts$design_type, "paired")
  expect_identical(v$contrasts$paired_samples, 1L)
  expect_match(v$contrasts$reason, "Only 1 sample has at least 10 cells in both 'A' and 'B'")
  # lowering min_cells rescues it
  v <- check_pseudobulk_design(p, data.frame(group_1 = "A", group_2 = "B"), min_cells = 2L)
  expect_identical(v$contrasts$method, "edgeR")
  # raising min_replicates is honoured
  v <- check_pseudobulk_design(p, data.frame(group_1 = "A", group_2 = "B"), min_cells = 2L, min_replicates = 3L)
  expect_identical(v$contrasts$method, "wilcoxon")
})

test_that("rest pools the other groups per sample", {
  p <- design_profiles(s1 = c(A = 20, B = 6, C = 6), s2 = c(A = 20, B = 5, C = 6), s3 = c(B = 30))
  v <- check_pseudobulk_design(p, data.frame(group_1 = c("A", "C"), group_2 = c("rest", "rest")), min_cells = 10L)
  sel <- v$selections[["A vs rest"]]
  rest <- sel$profiles[sel$profiles$group == "rest", ]
  expect_identical(rest$n_cells[match(c("s1", "s2", "s3"), rest$sample)], c(12L, 11L, 30L))
  expect_identical(v$contrasts$method[1], "edgeR")       # s1 and s2 paired (both sides >= 10)
  expect_identical(v$contrasts$design_type[1], "paired")
  expect_identical(v$contrasts$paired_samples[1], 2L)
  expect_identical(v$contrasts$method[2], "wilcoxon")    # C never reaches 10 cells
  expect_identical(sel$retained_members[[which(sel$retained$profile_id == "s1__rest")]], c("s1__B", "s1__C"))
})

test_that("input validation and accepted design inputs", {
  p <- design_profiles(s1 = c(A = 20), s2 = c(A = 20), s3 = c(B = 20), s4 = c(B = 20))
  expect_error(check_pseudobulk_design(1:3, data.frame(group_1 = "A", group_2 = "B")), "design")
  expect_error(check_pseudobulk_design(p, data.frame(x = 1)), "contrasts")
  expect_error(check_pseudobulk_design(p, data.frame(group_1 = "A", group_2 = "B"), min_cells = 0), "min_cells")
  pb <- fixture_pseudobulk(edger_fixture(FALSE))
  v <- check_pseudobulk_design(pb, data.frame(group_1 = "group_1", group_2 = "group_2"))
  expect_identical(v$contrasts$method, "edgeR")
  expect_identical(v$contrasts$samples_group_1, 3L)
  skip_if_no_scuba_pseudobulk()
  d <- SCUBA::pseudobulk_design(as_sce(edger_fixture(FALSE)), "sample_id", "condition")
  expect_identical(check_pseudobulk_design(d, data.frame(group_1 = "group_1", group_2 = "group_2"))$contrasts, v$contrasts)
})
