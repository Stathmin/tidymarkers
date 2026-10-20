# Fast by default, full on demand -- without losing truth.
#
# The suite has two kinds of checks:
#   * core logic and biology: design/gate/selection contracts, DAG invariants,
#     fixtures and regression tests.  These always run.
#   * external cross-checks: independent engines (BLAST+, Rsubread, rprimer),
#     annotation databases (VariantAnnotation/TxDb/BSgenome/EnsDb) and real
#     second-process backends (SOCK).  These are valuable but pay a heavy
#     one-off cost (package loads, index builds, worker processes).
#
# `Rscript dev/test.R` runs the fast suite.  `Rscript dev/test.R --full` (or an
# `R CMD check`, which sets `_R_CHECK_PACKAGE_NAME_`) runs both.

full_tests <- function() {
  identical(Sys.getenv("TIDYMARKERS_FULL_TESTS"), "1") ||
    nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
}

skip_if_not_full <- function() {
  if (!full_tests()) {
    testthat::skip("external cross-check; set TIDYMARKERS_FULL_TESTS=1 or run R CMD check")
  }
}
