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
  slot = lifecycle::deprecated()
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
  slot = lifecycle::deprecated()
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
  slot = lifecycle::deprecated()
)
```

## Arguments

- object:

  a single-cell object. Currently, Seurat, SingleCellExperiment, and
  AnnData objects are supported.

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
  `"scran"`).

- slot:

  No longer used. Please use `layer` instead.

## Methods (by class)

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
