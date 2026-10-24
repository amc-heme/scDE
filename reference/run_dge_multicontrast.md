# Multi-contrast pseudobulk differential expression

Runs a collection of two-group pseudobulk edgeR comparisons on one
object, aggregating the raw counts once (through
[`SCUBA::pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.html))
and testing every contrast on the resulting sample x group profiles.
Contrasts that are not replicated enough for edgeR (see
[`check_pseudobulk_design()`](https://amc-heme.github.io/scDE/reference/check_pseudobulk_design.md))
fall back to the cell-level Wilcoxon test of
[`run_dge()`](https://amc-heme.github.io/scDE/reference/run_dge.md); the
`method` column of the result says which test produced each row. This is
the natural way to compare cell-level groups such as clusters across
biological samples: those contrasts are paired designs
(`~ sample + group`), which edgeR handles.

## Usage

``` r
run_dge_multicontrast(
  object,
  group_by,
  sample_by,
  contrast_mode = c("pairwise", "one_vs_rest", "reference"),
  groups = NULL,
  reference_group = NULL,
  layer = NULL,
  seurat_assay = NULL,
  min_cells = 10L,
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  robust = FALSE,
  p_adjust_scope = c("contrast", "global", "both"),
  test_use = "edgeR"
)
```

## Arguments

- object:

  A Seurat, SingleCellExperiment, AnnData or file-backed SCUBA object
  (`scuba_source`), or a pre-aggregated `scuba_pseudobulk` from
  [`SCUBA::pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.html)
  (then `group_by` and `sample_by` are taken from the object and no
  fallback is possible).

- group_by:

  Metadata column containing the groups to compare.

- sample_by:

  Metadata column containing biological sample identifiers.

- contrast_mode:

  Comparison strategy: `"pairwise"` compares every pair of selected
  groups; `"one_vs_rest"` compares each selected group with all other
  observed groups pooled within each sample; `"reference"` compares
  every selected non-reference group with `reference_group`.

- groups:

  Optional character vector of groups to use. The default uses every
  non-missing value observed in `group_by`, in observed order. In
  `one_vs_rest` mode, this selects numerator groups; the rest side still
  pools all other observed groups.

- reference_group:

  Reference denominator required when `contrast_mode = "reference"`.

- layer:

  Raw-count layer. Defaults to `"counts"`. For AnnData objects whose
  counts live in `X` or `raw`, pass `"X"` or `"raw"` explicitly.

- seurat_assay:

  Assay / modality holding the raw counts (default: the object's default
  assay, as SCUBA defines it).

- min_cells:

  Minimum cells required per sample/comparison-side pseudobulk profile.

- positive_only:

  Whether to retain only positive log2 fold changes.

- remove_raw_pval:

  Whether to remove the raw p-value column.

- robust:

  Whether to use robust empirical Bayes estimation in
  [`edgeR::glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html).

- p_adjust_scope:

  Multiple-testing scope. `"contrast"` retains the within-contrast
  adjustment; `"global"` replaces `pval_adj` with BH adjustment over
  every returned gene/contrast test; `"both"` retains within-contrast
  `pval_adj` and adds `pval_adj_global`. Global adjustment is applied
  before `positive_only` filtering.

- test_use:

  Only `"edgeR"` is supported here; cell-level tests are the
  per-contrast fallback, not a choice.

## Value

A tibble with the columns `group_1`, `group_2`, `contrast`, `method`,
`group`, `feature`, `avgExpr`, `log2FC`, `pval`, `pval_adj` (and, for
rows from the Wilcoxon fallback, `auc`, `pct_in`, `pct_out`; `NA` on
edgeR rows). The `edger_multicontrast_details` attribute holds the
contrast table (with `method`, `design_type` and `reason`) and the
per-contrast details.
