# `testthat::skip_if_not_installed()` calls `requireNamespace()`, which loads the
# package; loading the Bioconductor stack re-emits the benign
# "replacing previous import ... when loading 'SummarizedExperiment'" warning.
# Checking with `system.file()` avoids the load, and the functions under test
# load the namespace quietly themselves (see `load_namespace_quietly()`).
skip_if_not_installed_noload <- function(pkg) {
  testthat::skip_if(!nzchar(system.file(package = pkg)), paste(pkg, "not installed"))
}
