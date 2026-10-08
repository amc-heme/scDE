test_that("edgeR supports Seurat Assay5 and legacy Assay", {
  skip_if_not_installed("SeuratObject")
  fixture <- edger_fixture(FALSE)

  assay5_object <- suppressWarnings(SeuratObject::CreateSeuratObject(
    Matrix::Matrix(fixture$counts, sparse = TRUE), meta.data = fixture$metadata
  ))
  expect_s4_class(assay5_object[["RNA"]], "Assay5")
  expect_edger_result(run_fixture_edger(assay5_object), "ordinary")

  legacy <- assay5_object
  legacy[["legacy"]] <- SeuratObject::CreateAssayObject(counts = fixture$counts)
  expect_s4_class(legacy[["legacy"]], "Assay")
  result <- suppressWarnings(run_dge(
    legacy, "condition", layer = "counts", seurat_assay = "legacy",
    test_use = "edgeR", sample_by = "sample_id",
    group_1 = "group_1", group_2 = "group_2", min_cells = 10L
  ))
  expect_edger_result(result, "ordinary")
})

test_that("edgeR supports BPCells-backed Seurat counts", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("SeuratObject")
  fixture <- edger_fixture(FALSE)
  path <- tempfile("scde-bpcells-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  sparse_counts <- methods::as(Matrix::Matrix(fixture$counts, sparse = TRUE), "dgCMatrix")
  suppressMessages(BPCells::write_matrix_dir(sparse_counts, dir = path))
  disk_counts <- BPCells::open_matrix_dir(path)
  object <- SeuratObject::CreateSeuratObject(disk_counts, meta.data = fixture$metadata)
  expect_edger_result(run_fixture_edger(object), "ordinary")
})

test_that("edgeR supports in-memory SingleCellExperiment counts", {
  skip_if_not_installed("SingleCellExperiment")
  fixture <- edger_fixture(FALSE)
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(fixture$counts, sparse = TRUE)),
    colData = fixture$metadata
  )
  expect_edger_result(run_fixture_edger(object), "ordinary")
})

test_that("edgeR supports DelayedArray and HDF5-backed SCE counts", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("DelayedArray")
  fixture <- edger_fixture(FALSE)
  delayed <- DelayedArray::DelayedArray(fixture$counts)
  delayed_object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = delayed), colData = fixture$metadata
  )
  expect_edger_result(run_fixture_edger(delayed_object), "ordinary")

  skip_if_not_installed("HDF5Array")
  h5_path <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_path), add = TRUE)
  h5_counts <- HDF5Array::writeHDF5Array(fixture$counts, filepath = h5_path,
                                         name = "counts", with.dimnames = TRUE)
  h5_object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = h5_counts), colData = fixture$metadata
  )
  expect_edger_result(run_fixture_edger(h5_object), "ordinary")
})

test_that("edgeR supports AnnData count layers", {
  skip_if_not_installed("anndata")
  skip_if_not_installed("reticulate")
  skip_if_no_python <- function() {
    available <- tryCatch(
      suppressWarnings(reticulate::py_module_available("anndata")),
      error = function(cnd) FALSE
    )
    if (!isTRUE(available)) skip("Python anndata unavailable")
  }
  skip_if_no_python()
  fixture <- edger_fixture(FALSE)
  # Recent Python anndata releases no longer accept a dtype argument here.
  # Let AnnData infer the dtype from the raw-count matrix.
  object <- anndata::AnnData(
    X = t(fixture$counts),
    obs = fixture$metadata,
    var = data.frame(row.names = rownames(fixture$counts))
  )
  object$layers[["counts"]] <- object$X
  expect_edger_result(run_fixture_edger(object), "ordinary")
})

test_that("Wilcoxon remains the default", {
  skip_if_not_installed("SingleCellExperiment")
  fixture <- edger_fixture(FALSE)
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = log1p(Matrix::Matrix(fixture$counts, sparse = TRUE))),
    colData = fixture$metadata
  )
  observed <- run_dge(object, group_by = "condition")
  expect_false(inherits(observed, "scDE_edger_results"))
  expect_true(all(c("group", "feature", "log2FC", "pval_adj") %in% colnames(observed)))
})
