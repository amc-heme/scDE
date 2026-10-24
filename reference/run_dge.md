# Differential Gene Expression

Performs differential expression testing on single cell data depending
on the class of the object passed. See the "Methods (by class)" section
for details on the test used.

## Usage

``` r
run_dge(
  object,
  group_by,
  layer = NULL,
  seurat_assay = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  test_use = NULL,
  sample_by = NULL,
  group_1 = NULL,
  group_2 = NULL,
  min_cells = 10L,
  robust = FALSE,
  slot = lifecycle::deprecated(),
  ...
)

# S3 method for class 'scuba_pseudobulk'
run_dge(
  object,
  group_by = NULL,
  layer = NULL,
  seurat_assay = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  test_use = "edgeR",
  sample_by = NULL,
  group_1 = NULL,
  group_2 = NULL,
  min_cells = 10L,
  robust = FALSE,
  slot = lifecycle::deprecated(),
  ...
)

# S3 method for class 'Seurat'
run_dge(
  object,
  group_by,
  layer = NULL,
  seurat_assay = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  test_use = NULL,
  sample_by = NULL,
  group_1 = NULL,
  group_2 = NULL,
  min_cells = 10L,
  robust = FALSE,
  slot = lifecycle::deprecated(),
  ...
)

# S3 method for class 'SingleCellExperiment'
run_dge(
  object,
  group_by,
  layer = NULL,
  seurat_assay = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  test_use = NULL,
  sample_by = NULL,
  group_1 = NULL,
  group_2 = NULL,
  min_cells = 10L,
  robust = FALSE,
  slot = lifecycle::deprecated(),
  ...
)

# S3 method for class 'AnnDataR6'
run_dge(
  object,
  group_by,
  layer = NULL,
  seurat_assay = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  test_use = NULL,
  sample_by = NULL,
  group_1 = NULL,
  group_2 = NULL,
  min_cells = 10L,
  robust = FALSE,
  slot = lifecycle::deprecated(),
  ...
)

# S3 method for class 'scuba_source'
run_dge(
  object,
  group_by,
  layer = NULL,
  seurat_assay = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  test_use = NULL,
  sample_by = NULL,
  group_1 = NULL,
  group_2 = NULL,
  min_cells = 10L,
  robust = FALSE,
  slot = lifecycle::deprecated(),
  ...
)
```

## Arguments

- object:

  a single-cell object. Currently, Seurat, SingleCellExperiment, AnnData
  and file-backed SCUBA (`scuba_source`, see
  [`SCUBA::open_h5ad()`](https://amc-heme.github.io/SCUBA/reference/open_h5ad.html))
  objects are supported, as well as pre-aggregated pseudobulk profiles
  (`scuba_pseudobulk`, see
  [`SCUBA::pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.html))
  for `test_use = "edgeR"`.

- group_by:

  metadata variable to use for forming differential gene expression
  groups. A group will be created for all values in this variable
  present in the object passed to this function.

- layer:

  The feature matrix layer to pull data from. If NULL, test-specific
  conventions will be used. The "data" layer is used for Seurat objects,
  and the "logcounts" layer is used for SingleCellExperiment objects.
  This parameter is currently ignored for anndata objects: only the "X"
  matrix may be used. The layer parameter is also ignored in BPCells
  objects.

- seurat_assay:

  for Seurat objects, the assay to use for evaluating differences in
  expression. This is used only for tests on Seurat objects. For
  SingleCellExperiment objects, the assay defaults to the main
  experiment, and for anndata objects the "X" matrix is used.

- lfc_format:

  sets the format used for reporting log-fold changes. If "log2",
  log-fold changes will be reported in a log two fold change format, and
  if "ln", log-fold changes will be reported as a natural log fold
  change. "log2" is the default value. Presto uses a natural log fold
  change by default, so `lfc_format` should be set to "ln" if
  consistency with Presto is desired.

- positive_only:

  If TRUE, only display genes that are up-regulated in a group compared
  to the reference (LFC \> 0). Defaults to FALSE.

- remove_raw_pval:

  If TRUE, show only the adjusted p-value. Defaults to FALSE.

- test_use:

  For SingleCellExperiment objects only. Controls which backend is used.
  `NULL` (default) auto-detects: `"Presto"` for in-memory matrices
  (`dgCMatrix` or `matrix`), `"scran"` for all other matrix classes
  (e.g. HDF5-backed `DelayedMatrix`, `SVT_SparseMatrix`). `"Presto"`
  forces presto via
  [`presto::wilcoxauc`](https://rdrr.io/pkg/presto/man/wilcoxauc.html) —
  errors if the matrix is not `dgCMatrix` or `matrix`. `"scran"` forces
  a true one-vs-rest Wilcoxon test per group via
  [`scran::pairwiseWilcox`](https://rdrr.io/pkg/scran/man/pairwiseWilcox.html),
  indexing the result directly (returns an `auc` column; LFC from group
  means). `"scran_pairwise"` uses
  [`scran::findMarkers`](https://rdrr.io/pkg/scran/man/findMarkers.html)
  with pairwise comparisons combined across all groups
  (`pval.type = "any"`); no `auc` column is returned and `log2FC` is
  computed from group vs. all-other-cells means (same method as
  `"scran"`). For every supported object class, `test_use = "edgeR"`
  runs a sample-level (pseudobulk) edgeR quasi-likelihood test of
  `group_1` versus `group_2`: raw counts are summed per sample and group
  by
  [`SCUBA::pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.html)
  (the cell-level matrix is never materialised), profiles with fewer
  than `min_cells` cells are dropped,
  [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
  and TMM normalisation are applied, and the model is `~ 0 + group` when
  every sample belongs to one group (independent design) or
  `~ sample + group` when samples contribute cells to both groups
  (paired design, e.g. cluster versus cluster within the same donors).
  edgeR needs at least two samples per group (or two paired samples;
  three or more are recommended). When a contrast is not replicated
  enough, the same call runs the cell-level Wilcoxon test of the
  object's class on the cells of the two groups instead, says so in a
  message, and labels the rows `method = "wilcoxon"` (see
  [`check_pseudobulk_design()`](https://amc-heme.github.io/scDE/reference/check_pseudobulk_design.md)
  to predict this from metadata alone). Pseudobulk tests error if
  `sample_by` yields about one cell per profile (a cell-level comparison
  in disguise). All other values preserve the existing cell-level
  Wilcoxon behaviour.

- sample_by:

  Metadata column containing biological sample identifiers. Required
  only for `test_use = "edgeR"`.

- group_1, group_2:

  The two values of `group_by` to compare with a pseudobulk test.
  `group_1` is always the numerator, so positive `log2FC` values mean
  higher expression in `group_1`.

- min_cells:

  Minimum cells required in each sample/group pseudobulk.

- robust:

  Whether to use robust empirical Bayes estimation in
  [`edgeR::glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html).

- slot:

  No longer used. Please use `layer` instead.

- ...:

  Additional arguments reserved for methods.

## Methods (by class)

- `run_dge(scuba_pseudobulk)`: Pre-aggregated pseudobulk profiles from
  [`SCUBA::pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.html):
  edgeR on the stored sample x group counts. `group_by` and `sample_by`
  come from the object; an unreplicated contrast is an error, because no
  cells are available for the Wilcoxon fallback.

- `run_dge(Seurat)`: Seurat objects: uses Presto

- `run_dge(SingleCellExperiment)`: SingleCellExperiment objects: uses
  Presto by default for in-memory matrices, with automatic fallback to
  scran for DelayedArray-backed matrices (e.g. HDF5-backed via
  `HDF5Array`). The backend can be controlled explicitly via `test_use`.

  **Presto** (`test_use = "Presto"` or auto for in-memory): calls
  [`presto::wilcoxauc`](https://rdrr.io/pkg/presto/man/wilcoxauc.html)
  directly on the extracted assay matrix. Output schema matches the
  Seurat/Presto path and includes `auc`, `pct_in`, and `pct_out`.

  **scran one-vs-rest** (`test_use = "scran"` or auto for
  non-presto-compatible matrices): runs a true one-vs-rest Wilcoxon test
  per group via
  [`scran::pairwiseWilcox`](https://rdrr.io/pkg/scran/man/pairwiseWilcox.html)
  on a binary factor (group vs. all others), then indexes the result
  directly by pair position. Returns an `auc` column. LFC is computed
  from per-group means (log-space subtraction). Requires
  `BiocManager::install(c("scran", "SummarizedExperiment"))`.

  **scran pairwise** (`test_use = "scran_pairwise"`): uses
  `scran::findMarkers(pval.type = "any")`. No `auc` column. `log2FC` is
  computed from group vs. all-other-cells means (same method as
  `"scran"`). Requires
  `BiocManager::install(c("scran", "SummarizedExperiment"))`.

  In all cases, `avgExpr` for scran paths is computed via
  `MatrixGenerics::rowMeans2()`, which supports both in-memory
  (`dgCMatrix`, `matrix`) and `DelayedArray`-backed matrices (e.g.
  `HDF5Matrix`, `SVT_SparseMatrix`). `other_mean` is derived
  algebraically from per-group means — no additional matrix reads per
  group.

- `run_dge(AnnDataR6)`: Anndata objects: uses Scanpy's rank_genes_groups
  function, with the "wilcoxon" method

- `run_dge(scuba_source)`: File-backed SCUBA objects (`scuba_source`,
  from
  [`SCUBA::open_h5ad()`](https://amc-heme.github.io/SCUBA/reference/open_h5ad.html)
  /
  [`SCUBA::open_h5mu()`](https://amc-heme.github.io/SCUBA/reference/open_h5ad.html)):
  uses Presto
  ([`presto::wilcoxauc`](https://rdrr.io/pkg/presto/man/wilcoxauc.html))
  on the expression matrix served by SCUBA's Julia engine. The matrix is
  pulled through
  [`SCUBA::fetch_matrix()`](https://amc-heme.github.io/SCUBA/reference/fetch_matrix.html)
  in blocks of features (`options(scDE.chunk_features)`, default 2000)
  so memory stays bounded for large objects; the Wilcoxon statistics are
  computed per feature, so blocking does not change them, and the
  Benjamini-Hochberg adjustment is applied over all features within each
  group afterwards, as Presto does. Groups come from
  [`SCUBA::fetch_metadata()`](https://amc-heme.github.io/SCUBA/reference/fetch_metadata.html),
  so metadata columns added in R (e.g. by
  [`SCUBA::update_object_metadata()`](https://amc-heme.github.io/SCUBA/reference/update_object_metadata.html))
  are honoured. `seurat_assay` selects the modality (default: the main
  modality, `"X"` for h5ad files); `layer` defaults to
  [`SCUBA::default_layer()`](https://amc-heme.github.io/SCUBA/reference/default_layer.html)
  (`"X"`), which is assumed to hold log-normalised values like a Seurat
  `data` layer. Output schema matches the Seurat/Presto path:
  `group, feature, avgExpr, log2FC, auc, pval, pval_adj, pct_in, pct_out`.
  `test_use = "edgeR"` (pseudobulk) is handled before dispatch, see
  above.

## Pseudobulk results

edgeR results carry a `method` column (`"edgeR"`, or `"wilcoxon"` after
a fallback) ahead of the common columns `group`, `feature`, `avgExpr`
(average log2-CPM for edgeR), `log2FC`, `pval`, `pval_adj`, and an
`edger_details` attribute: `design_type` (`"independent"` / `"paired"`),
`fallback_reason`, `sample_summary` (retained and excluded samples per
group), `retained_profiles`, `excluded_profiles`, the design matrix and
contrast or coefficient, `genes_tested`, the TMM factors and how the
counts were aggregated.
