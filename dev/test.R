#!/usr/bin/env Rscript
# =============================================================================
#  dev/test.R  --  the routine test loop
# =============================================================================
#  Three modes, from cheapest to complete:
#
#    Rscript dev/test.R                  # lean: shared contracts only (~1 min)
#    Rscript dev/test.R --scheme crispr  # one family + shared contracts
#    Rscript dev/test.R --full           # every test, external cross-checks on
#    Rscript dev/test.R <regex>          # legacy: raw testthat filter
#
#  The shared contracts (tests/testthat/test-{criteria,decision,kasp-checks,
#  plot-methods,technology}.R) encode the invariants every scheme inherits, so
#  they are the minimum gate for a routine change.  `--scheme <regex>` adds the
#  test files whose names match the regex (e.g. crispr, lamp|rpa, indel) on top
#  of them; the regex is matched against the file name, never the file body.
#  `--full` runs every test with TIDYMARKERS_FULL_TESTS=1, which turns on the
#  external cross-checks (BLAST+, Rsubread, rprimer, annotation DBs, real SOCK
#  workers); R CMD check also sets that flag, so the full path is never weakened.
#
#  Lean and scheme modes run files in parallel (2 workers).  Full mode is serial:
#  the SOCK tests already spawn worker processes, and nesting them under
#  testthat's workers crashes the pool.
#
#  Tests never install.  Install once with `Rscript dev/install-deps.R`
#  (`--with-suggests` for the full external checks).
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

if (!nzchar(system.file(package = "tidymarkers"))) {
  stop("tidymarkers is not installed; run Rscript dev/install-deps.R first")
}

# The cross-scheme contracts.  Matching on file names keeps this independent of
# test content: a renamed contract file is a deliberate edit here, not a silent
# drop from the routine gate.
shared <- c("criteria", "decision", "kasp-checks", "plot-methods", "technology")

args <- commandArgs(trailingOnly = TRUE)
full <- "--full" %in% args
args <- setdiff(args, "--full")

scheme <- NULL
if ("--scheme" %in% args) {
  i <- match("--scheme", args)
  if (i == length(args)) {
    stop("--scheme needs a regex, e.g. --scheme crispr", call. = FALSE)
  }
  scheme <- args[[i + 1L]]
  args <- args[-c(i, i + 1L)]
}

filter <- if (length(args)) args[[1L]] else NULL
if (length(args) > 1L) {
  stop("unexpected arguments: ", paste(args[-1L], collapse = " "), call. = FALSE)
}
if (full && (!is.null(scheme) || !is.null(filter))) {
  stop("--full runs every test; drop --scheme/<filter>", call. = FALSE)
}

files <- list.files("tests/testthat", pattern = "^test-.*[.]R$")

if (full) {
  filter <- NULL
  Sys.setenv(TIDYMARKERS_FULL_TESTS = "1", TESTTHAT_PARALLEL = "FALSE")
  message("test.R: full -- all tests, external cross-checks on")
} else {
  # an inherited TIDYMARKERS_FULL_TESTS must not silently turn the routine loop
  # into the release suite; --full is the switch
  Sys.unsetenv("TIDYMARKERS_FULL_TESTS")
  if (!nzchar(Sys.getenv("TESTTHAT_PARALLEL"))) {
    Sys.setenv(TESTTHAT_PARALLEL = "TRUE")
  }
  if (!is.null(scheme)) {
    hits <- files[grepl(scheme, files)]
    # an empty match would silently run the shared contracts and look green
    if (!length(hits)) {
      stop("--scheme ", scheme, " matches no test file under tests/testthat/ ",
        "(try crispr, lamp|rpa, indel)",
        call. = FALSE
      )
    }
    filter <- paste(c(shared, scheme), collapse = "|")
    message("test.R: scheme '", scheme, "' + shared contracts (", length(hits), " scheme file(s))")
  } else if (is.null(filter)) {
    filter <- paste(shared, collapse = "|")
    message("test.R: lean -- shared contracts")
  } else if (!any(grepl(filter, files))) {
    # a legacy filter that matches nothing is a silent no-op, not a pass
    stop("filter ", filter, " matches no test file under tests/testthat/", call. = FALSE)
  }
}

testthat::test_local(filter = filter, stop_on_failure = TRUE)
