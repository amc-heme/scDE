# Fixtures for the pseudobulk tests: a small count matrix with a known shifted-gene signal,
# builders that wrap it in every supported object class, and an aggregation helper that
# builds the `scuba_pseudobulk` object the backend consumes.

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

# A fixture with three clusters in every one of three samples (paired designs everywhere).
multicontrast_fixture_data <- function() {
  set.seed(20260925)
  clusters <- c("A", "B", "C")
  samples <- paste0("sample", 1:3)
  cells_per_profile <- 12L
  profile <- expand.grid(cluster = clusters, sample_id = samples, stringsAsFactors = FALSE)
  metadata <- profile[rep(seq_len(nrow(profile)), each = cells_per_profile), ]
  rownames(metadata) <- paste0("cell", seq_len(nrow(metadata)))
  counts <- matrix(
    stats::rpois(40L * nrow(metadata), lambda = 3),
    nrow = 40L,
    dimnames = list(paste0("gene", 1:40), rownames(metadata))
  )
  counts[1:3, metadata$cluster == "A"] <- counts[1:3, metadata$cluster == "A"] + 10L
  counts[4:6, metadata$cluster == "B"] <- counts[4:6, metadata$cluster == "B"] + 8L
  list(counts = counts, metadata = metadata)
}

# --- object builders -------------------------------------------------------------------------

as_sce <- function(fixture, counts_layer = "counts") {
  assays <- list(Matrix::Matrix(fixture$counts, sparse = TRUE))
  names(assays) <- counts_layer
  assays$logcounts <- log1p(Matrix::Matrix(fixture$counts, sparse = TRUE))
  SingleCellExperiment::SingleCellExperiment(assays = assays, colData = fixture$metadata)
}
as_seurat <- function(fixture) {
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(
    Matrix::Matrix(fixture$counts, sparse = TRUE), meta.data = fixture$metadata
  ))
  object <- SeuratObject::SetAssayData(object, layer = "data",
                                       new.data = log1p(Matrix::Matrix(fixture$counts, sparse = TRUE)))
  object
}
as_legacy_seurat <- function(fixture) {
  object <- as_seurat(fixture)
  object[["legacy"]] <- SeuratObject::CreateAssayObject(counts = fixture$counts)
  object
}
as_bpcells_seurat <- function(fixture, dir = tempfile("scde-bpcells-")) {
  sparse_counts <- methods::as(Matrix::Matrix(fixture$counts, sparse = TRUE), "dgCMatrix")
  suppressMessages(BPCells::write_matrix_dir(sparse_counts, dir = dir))
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(
    BPCells::open_matrix_dir(dir), meta.data = fixture$metadata
  ))
  object[["RNA"]]$data <- log1p(Matrix::Matrix(fixture$counts, sparse = TRUE))
  object
}
as_hdf5_sce <- function(fixture, path = tempfile(fileext = ".h5")) {
  h5_counts <- HDF5Array::writeHDF5Array(fixture$counts, filepath = path, name = "counts", with.dimnames = TRUE)
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = h5_counts, logcounts = log1p(Matrix::Matrix(fixture$counts, sparse = TRUE))),
    colData = fixture$metadata
  )
}
# Built through reticulate rather than anndata::AnnData(): the R package (0.8) still passes a
# `dtype` argument that Python anndata >= 0.11 rejects. The result is a plain reticulate
# `anndata.AnnData` object, which is also what anndata::read_h5ad() returns with anndata >= 0.8.
as_anndata <- function(fixture) {
  ad <- reticulate::import("anndata", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)
  object <- ad$AnnData(
    X = reticulate::r_to_py(t(log1p(fixture$counts))),
    obs = reticulate::r_to_py(fixture$metadata),
    var = pd$DataFrame(index = reticulate::r_to_py(rownames(fixture$counts)))
  )
  reticulate::py_set_item(object$layers, "counts", reticulate::r_to_py(t(fixture$counts)))
  object
}
as_scuba_source <- function(fixture, path = tempfile(fileext = ".h5ad")) {
  as_anndata(fixture)$write_h5ad(path)
  SCUBA::open_h5ad(path)
}

# `scuba_pseudobulk`-shaped object straight from the fixture (no SCUBA engine involved), so
# the backend can be unit-tested on exact, hand-computed sums.
fixture_pseudobulk <- function(fixture, group_by = "condition", sample_by = "sample_id") {
  key <- paste(fixture$metadata[[sample_by]], fixture$metadata[[group_by]], sep = "__")
  counts <- t(rowsum(t(fixture$counts), key))
  storage.mode(counts) <- "double"
  profiles <- data.frame(
    profile_id = colnames(counts),
    sample = sub("__.*$", "", colnames(counts)),
    group = sub("^.*__", "", colnames(counts)),
    n_cells = as.integer(table(key)[colnames(counts)]),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  structure(list(
    counts = counts, profiles = profiles, assay = "RNA", layer = "counts", fun = "sum",
    sample_by = sample_by, group_by = group_by,
    n_cells_dropped = c(missing_sample = 0L, missing_group = 0L), engine = "r"
  ), class = "scuba_pseudobulk")
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
  expect_true(all(c("method", "group", "feature", "avgExpr", "log2FC", "pval", "pval_adj") %in%
                    colnames(result)))
  expect_true(all(result$method == "edgeR"))
  expect_true(all(result$group == "group_1"))
  expected_features <- paste0("gene", 1:3)
  expect_true(all(result$log2FC[match(expected_features, result$feature)] > 0))
  expect_identical(attr(result, "edger_details")$design_type, design_type)
  expect_identical(attr(result, "edger_details")$normalization_method, "TMM")
}

# Sort a result table for a column-wise comparison across object classes.
sorted_result <- function(result) {
  result <- as.data.frame(result)
  result <- result[order(result$feature), c("feature", "avgExpr", "log2FC", "pval", "pval_adj")]
  rownames(result) <- NULL
  result
}
