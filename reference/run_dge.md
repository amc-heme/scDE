# Differential Gene Expression

Performs differential expression testing on single cell data depending
on the class of the object passed. See the "Methods (by class)" section
for details on the test used.

## Usage

``` r
run_dge(
  object,
  group_by,
  seurat_assay = NULL,
  slot = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE,
  ...
)

# S3 method for class 'Seurat'
run_dge(
  object,
  group_by,
  seurat_assay = NULL,
  slot = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE
)

# S3 method for class 'AnnDataR6'
run_dge(
  object,
  group_by,
  seurat_assay = NULL,
  slot = NULL,
  lfc_format = "log2",
  positive_only = FALSE,
  remove_raw_pval = FALSE
)
```

## Arguments

- object:

  a single-cell object. Currently, Seurat and SingleCellExperiment
  objects are supported.

- group_by:

  metadata variable to use for forming differential gene expression
  groups. A group will be created for all values in this variable
  present in the object passed to this function.

- seurat_assay:

  for Seurat objects, the assay to use for evaluating differences in
  expression. This is used only for tests on Seurat objects. For
  SingleCellExperiment objects, the assay defaults to the main
  experiment, and for anndata objects the "X" matrix is used.

- slot:

  The feature matrix slot to pull data from. If NULL, the "data" slot is
  used for Seurat objects, and the "logcounts" slot is used for
  SingleCellExperiment objects. This parameter is currently ignored for
  anndata objects: only the "X" matrix may be used.

- lfc_format:

  sets the format used for reporting log-fold changes. If "log2",
  log-fold changes will be reported in a log two fold change format, and
  if "ln", log-fold changes will be reported as a natural log fold
  change. "log2" is the default value. Presto uses a natural log fold
  change by default, so \`lfc_format\` should be set to "ln" if
  consistency with Presto is desired.

- positive_only:

  If TRUE, only display genes that are up-regulated in a group compared
  to the reference (LFC \> 0). Defaults to FALSE.

- remove_raw_pval:

  If TRUE, show only the adjusted p-value. Defaults to FALSE.

## Methods (by class)

- `run_dge(Seurat)`: Seurat objects: uses Presto

- `run_dge(AnnDataR6)`: Anndata objects: uses Scanpy's rank_genes_groups
  function, with the "wilcoxon" method
