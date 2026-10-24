# Skip helpers shared by the pseudobulk test files. Three tiers of optional dependencies:
# edgeR (Bioconductor), a Python with anndata (reticulate), and SCUBA's Julia engine.

skip_if_no_edger <- function() {
  skip_if_not_installed("edgeR", "4.0.0")
}

skip_if_no_scuba_pseudobulk <- function() {
  skip_if_not(exists("pseudobulk", asNamespace("SCUBA")), "SCUBA without pseudobulk()")
}

skip_if_no_python_anndata <- function() {
  skip_if_not_installed("anndata")
  skip_if_not_installed("reticulate")
  available <- tryCatch(
    suppressWarnings(reticulate::py_module_available("anndata")),
    error = function(cnd) FALSE
  )
  if (!isTRUE(available)) skip("Python anndata unavailable")
}

skip_if_no_julia <- function() {
  skip_if_not_installed("SCUBA")
  skip_if_not(exists("open_h5ad", asNamespace("SCUBA")), "SCUBA without open_h5ad()")
  skip_if_not(
    isTRUE(tryCatch(SCUBA::scuba_julia_available(quiet = TRUE), error = function(e) FALSE)),
    "SCUBA's Julia engine is not available"
  )
  skip_if_not(
    isTRUE(tryCatch(JuliaCall::julia_eval("isdefined(ScubaR, :pseudobulk)"), error = function(e) FALSE)),
    "SCUBA.jl without ScubaR.pseudobulk"
  )
}
