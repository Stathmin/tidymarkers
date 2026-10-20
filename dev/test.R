#!/usr/bin/env Rscript
# =============================================================================
#  dev/test.R  --  the fast test loop
# =============================================================================
#  The core suite (design/gate/selection contracts, DAG invariants, fixtures)
#  always runs.  External cross-checks (BLAST+, Rsubread FM, rprimer,
#  VariantAnnotation/TxDb/BSgenome/EnsDb, real SOCK workers) are opt-in because
#  their one-off cost dominates a routine run.
#
#    Rscript dev/test.R               # fast: core suite only
#    Rscript dev/test.R --full        # + every external cross-check
#    Rscript dev/test.R crispr        # filter one file (as test_local(filter=))
#
#  Fast mode runs files in parallel (2 workers).  Full mode is serial: the SOCK
#  tests already spawn worker processes, and nesting them under testthat's
#  workers crashes the pool.  `R CMD check` is serial too (config is false).
#
#  Tests never install.  Install once with `Rscript dev/install-deps.R`
#  (`--with-suggests` for the full external checks).
# =============================================================================

if (!nzchar(system.file(package = "tidymarkers"))) {
  stop("tidymarkers is not installed; run Rscript dev/install-deps.R first")
}

args <- commandArgs(trailingOnly = TRUE)
full <- "--full" %in% args
args <- setdiff(args, "--full")
if (full) {
  Sys.setenv(TIDYMARKERS_FULL_TESTS = "1", TESTTHAT_PARALLEL = "FALSE")
} else if (!nzchar(Sys.getenv("TESTTHAT_PARALLEL"))) {
  Sys.setenv(TESTTHAT_PARALLEL = "TRUE")
}
filter <- if (length(args)) args[[1L]] else NULL
testthat::test_local(filter = filter, stop_on_failure = TRUE)
