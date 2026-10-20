#!/usr/bin/env Rscript
# =============================================================================
#  dev/install-deps.R  --  fast, tiered dependency install
# =============================================================================
#  The package itself needs only `Imports`.  The heavy `Suggests`
#  (VariantAnnotation, BSgenome, TxDb, Rsubread, ensembldb, rprimer, pkgdown,
#  ...) are opt-in: installing them by default is what made setup take many
#  minutes.
#
#    Rscript dev/install-deps.R                 # Imports + dev/test tools + package
#    Rscript dev/install-deps.R --with-suggests # + every Suggests (full gate)
#
#  `pak` reads DESCRIPTION, downloads in parallel and uses binaries where the
#  platform provides them.  Do not use `renv::restore()` for routine setup: the
#  lockfile is a historical 262-package dev snapshot.
#
#  The package itself is installed last (`R CMD INSTALL`), so the full test
#  suite can run SOCK workers and subprocesses against the installed copy.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
with_suggests <- "--with-suggests" %in% args

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cloud.r-project.org")
}

deps <- c("Depends", "Imports", "LinkingTo")
if (with_suggests) deps <- c(deps, "Suggests")
pak::local_install_dev_deps(".", dependencies = deps)

# The test loop and the local gate are always needed; they are small and do not
# pull the Bioconductor annotation stack.
tools <- c("testthat", "pkgload", "roxygen2", "styler", "lintr", "snow", "svglite")
missing <- tools[!vapply(tools, function(p) nzchar(system.file(package = p)), logical(1))]
if (length(missing)) {
  pak::pkg_install(missing, ask = FALSE)
}

# Install the package itself: the full tests exercise the installed copy.
lib <- .libPaths()[1L]
message("installing tidymarkers into ", lib)
status <- system2(
  file.path(R.home("bin"), "R"),
  c(
    "CMD", "INSTALL", "--no-multiarch", "--no-byte-compile",
    paste0("--library=", lib), "."
  )
)
if (!identical(as.integer(status), 0L)) {
  stop("R CMD INSTALL failed (status ", status, ")")
}
