# Decide which contrasts can be tested with pseudobulk edgeR

Applies scDE's replication rules to the sample x group profiles of an
object, without reading any expression data. For every requested
contrast it reports whether edgeR can run (`method = "edgeR"`) or
whether
[`run_dge()`](https://amc-heme.github.io/scDE/reference/run_dge.md) /
[`run_dge_multicontrast()`](https://amc-heme.github.io/scDE/reference/run_dge_multicontrast.md)
will fall back to the cell-level Wilcoxon test (`method = "wilcoxon"`),
and why. The same function drives the tests themselves, so a preflight
in an application and the eventual run always agree.

## Usage

``` r
check_pseudobulk_design(
  design,
  contrasts,
  min_cells = 10L,
  min_replicates = 2L
)
```

## Arguments

- design:

  A `scuba_pseudobulk_design` from
  [`SCUBA::pseudobulk_design()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk_design.html),
  a `scuba_pseudobulk` from
  [`SCUBA::pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.html),
  or a data.frame with columns `sample`, `group` and `n_cells`.

- contrasts:

  A data.frame with character columns `group_1` and `group_2` (one row
  per contrast; `group_2` may be `"rest"`).

- min_cells:

  Minimum cells for a sample x group profile to count.

- min_replicates:

  Minimum samples per side (independent design) or paired samples
  (paired design) for edgeR; the default of 2 is the hard minimum, three
  or more are recommended and fewer produce a warning when the test
  runs.

## Value

A list with

- `contrasts`: data.frame with one row per contrast: `group_1`,
  `group_2`, `contrast`, `method` (`"edgeR"` or `"wilcoxon"`),
  `design_type` (`"independent"` or `"paired"`), `samples_group_1`,
  `samples_group_2` (samples retained per side), `paired_samples`, and
  `reason` (why the contrast falls back, `NA` otherwise);

- `any_edger`, `any_fallback`;

- `warnings`: notes about contrasts that run with fewer than three
  replicates or with excluded profiles;

- `selections`: the per-contrast profile tables (internal detail used by
  the tests).

## Details

A profile (one sample in one group) counts when it has at least
`min_cells` cells. The design of a contrast is *paired* when at least
one sample has both groups (cluster versus cluster, treated versus
untreated cells of the same donor); edgeR then needs `min_replicates`
such samples and fits `~ sample + group`. Otherwise the design is
*independent* and edgeR needs `min_replicates` samples on each side
(`~ 0 + group`). `group_2 = "rest"` pools, per sample, every group other
than `group_1`.

## See also

[`run_dge()`](https://amc-heme.github.io/scDE/reference/run_dge.md),
[`run_dge_multicontrast()`](https://amc-heme.github.io/scDE/reference/run_dge_multicontrast.md),
[`SCUBA::pseudobulk_design()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk_design.html)

## Examples

``` r
if (FALSE) { # \dontrun{
design <- SCUBA::pseudobulk_design(SCUBA::AML_Seurat, "Batch", "condensed_cell_type")
check_pseudobulk_design(
  design,
  contrasts = data.frame(group_1 = "B Cells", group_2 = c("NK Cells", "rest"))
)$contrasts
} # }
```
