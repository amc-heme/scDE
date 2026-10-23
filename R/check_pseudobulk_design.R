#' Decide which contrasts can be tested with pseudobulk edgeR
#'
#' Applies scDE's replication rules to the sample x group profiles of an object, without reading
#' any expression data. For every requested contrast it reports whether edgeR can run
#' (`method = "edgeR"`) or whether `run_dge()` / `run_dge_multicontrast()` will fall back to the
#' cell-level Wilcoxon test (`method = "wilcoxon"`), and why. The same function drives the
#' tests themselves, so a preflight in an application and the eventual run always agree.
#'
#' A profile (one sample in one group) counts when it has at least `min_cells` cells. The
#' design of a contrast is *paired* when at least one sample has both groups (cluster versus
#' cluster, treated versus untreated cells of the same donor); edgeR then needs
#' `min_replicates` such samples and fits `~ sample + group`. Otherwise the design is
#' *independent* and edgeR needs `min_replicates` samples on each side (`~ 0 + group`).
#' `group_2 = "rest"` pools, per sample, every group other than `group_1`.
#'
#' @param design A `scuba_pseudobulk_design` from [SCUBA::pseudobulk_design()], a
#'   `scuba_pseudobulk` from [SCUBA::pseudobulk()], or a data.frame with columns `sample`,
#'   `group` and `n_cells`.
#' @param contrasts A data.frame with character columns `group_1` and `group_2` (one row per
#'   contrast; `group_2` may be `"rest"`).
#' @param min_cells Minimum cells for a sample x group profile to count.
#' @param min_replicates Minimum samples per side (independent design) or paired samples
#'   (paired design) for edgeR; the default of 2 is the hard minimum, three or more are
#'   recommended and fewer produce a warning when the test runs.
#'
#' @return A list with
#'   - `contrasts`: data.frame with one row per contrast: `group_1`, `group_2`, `contrast`,
#'     `method` (`"edgeR"` or `"wilcoxon"`), `design_type` (`"independent"` or `"paired"`),
#'     `samples_group_1`, `samples_group_2` (samples retained per side), `paired_samples`,
#'     and `reason` (why the contrast falls back, `NA` otherwise);
#'   - `any_edger`, `any_fallback`;
#'   - `warnings`: notes about contrasts that run with fewer than three replicates or with
#'     excluded profiles;
#'   - `selections`: the per-contrast profile tables (internal detail used by the tests).
#' @seealso [run_dge()], [run_dge_multicontrast()], [SCUBA::pseudobulk_design()]
#' @export
#'
#' @examples
#' \dontrun{
#' design <- SCUBA::pseudobulk_design(SCUBA::AML_Seurat, "Batch", "condensed_cell_type")
#' check_pseudobulk_design(
#'   design,
#'   contrasts = data.frame(group_1 = "B Cells", group_2 = c("NK Cells", "rest"))
#' )$contrasts
#' }
check_pseudobulk_design <- function(design, contrasts, min_cells = 10L, min_replicates = 2L) {
  profiles <- if (inherits(design, c("scuba_pseudobulk_design", "scuba_pseudobulk"))) {
    design$profiles
  } else if (is.data.frame(design)) {
    design
  } else {
    stop("`design` must come from SCUBA::pseudobulk_design() / SCUBA::pseudobulk(), or be a ",
         "data.frame with columns sample, group and n_cells.", call. = FALSE)
  }
  needed <- c("sample", "group", "n_cells")
  if (!all(needed %in% names(profiles))) {
    stop("`design` needs the columns ", paste(needed, collapse = ", "), ".", call. = FALSE)
  }
  if (!"profile_id" %in% names(profiles)) {
    profiles$profile_id <- make.unique(paste(profiles$sample, profiles$group, sep = "__"))
  }
  if (!is.data.frame(contrasts) || !all(c("group_1", "group_2") %in% names(contrasts)) ||
      !nrow(contrasts)) {
    stop("`contrasts` must be a data.frame with columns group_1 and group_2.", call. = FALSE)
  }
  .scde_validate_min_cells(min_cells)
  .scde_validate_min_cells(min_replicates)
  selections <- lapply(seq_len(nrow(contrasts)), function(i) {
    .scde_contrast_profiles(profiles, as.character(contrasts$group_1[i]),
                            as.character(contrasts$group_2[i]), min_cells, min_replicates)
  })
  table <- data.frame(
    group_1 = as.character(contrasts$group_1),
    group_2 = as.character(contrasts$group_2),
    stringsAsFactors = FALSE
  )
  table$contrast <- paste(table$group_1, "vs", table$group_2)
  table$method <- ifelse(vapply(selections, `[[`, logical(1), "ok"), "edgeR", "wilcoxon")
  table$design_type <- vapply(selections, `[[`, character(1), "design_type")
  table$samples_group_1 <- vapply(selections, `[[`, integer(1), "samples_group_1")
  table$samples_group_2 <- vapply(selections, `[[`, integer(1), "samples_group_2")
  table$paired_samples <- vapply(selections, `[[`, integer(1), "paired_samples")
  table$reason <- vapply(selections, `[[`, character(1), "reason")
  names(selections) <- table$contrast
  warnings <- unlist(lapply(seq_along(selections), function(i) {
    if (length(selections[[i]]$warnings)) paste0(table$contrast[i], ": ", selections[[i]]$warnings) else character()
  }), use.names = FALSE)
  list(
    contrasts = table,
    any_edger = any(table$method == "edgeR"),
    any_fallback = any(table$method == "wilcoxon"),
    warnings = warnings,
    min_cells = as.integer(min_cells),
    min_replicates = as.integer(min_replicates),
    selections = selections
  )
}
