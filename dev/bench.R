# =============================================================================
#  dev/bench.R  --  time and memory per function (a regression baseline)
# =============================================================================
#  Times a representative call of the hot functions and reports the memory each
#  allocates, so a slow-down or a new allocation shows up as a number rather
#  than a hunch.  Backed by `bench` (median time + `mem_alloc`) and `lobstr`
#  (object sizes); both are pinned in renv.lock.
#
#  Needs OligoArrayAux for the thermodynamic rows (see
#  data-raw/install-oligoarrayaux.R); without it those rows are skipped.
#
#  Run:  Rscript dev/bench.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

for (pkg in c("bench", "lobstr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("dev/bench.R needs '", pkg, "'; run renv::restore()")
}
suppressMessages(pkgload::load_all(root, quiet = TRUE))

# --- fixture ----------------------------------------------------------------
mk_fixture <- function(seed = 7L, n = 2000L) {
  withr::with_seed(seed, {
    w <- paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
    list(genome = Biostrings::DNAStringSet(c(Chr1 = w)), primer = substr(w, 500L, 524L))
  })
}
fx <- mk_fixture()

cases <- list(
  primer_tm        = function() primer_tm(fx$primer),
  gc_frac          = function() gc_frac(fx$primer),
  snp_template     = function() snp_template(fx$genome, "Chr1", 1000L, flank = 150L),
  in_silico_pcr    = function() in_silico_pcr(fx$primer, as.character(Biostrings::reverseComplement(Biostrings::DNAString(fx$primer))), fx$genome),
  specificity_seed = function() specificity_seed(c(P = fx$primer), fx$genome)
)
if (nzchar(Sys.which("hybrid-min"))) cases$oligo_dg <- function() oligo_dg(fx$primer)

# --- run --------------------------------------------------------------------
cat("tidykasp function benchmark (bench)\n")
if (!nzchar(Sys.which("hybrid-min"))) cat("note: OligoArrayAux not on PATH - thermodynamic rows skipped\n")
cat("\n")
for (nm in names(cases)) {
  res <- bench::mark(cases[[nm]](), iterations = 20L, check = FALSE, filter_gc = FALSE)
  cat(sprintf("  %-18s %9.3f ms  %10s\n", nm, as.numeric(res$median) * 1000, format(res$mem_alloc)))
}

cat("\nobject sizes (lobstr):\n")
objects <- list(
  genome = fx$genome,
  template = snp_template(fx$genome, "Chr1", 1000L, flank = 150L)
)
for (nm in names(objects)) cat(sprintf("  %-18s %10s\n", nm, format(lobstr::obj_size(objects[[nm]]))))
cat("\n")
