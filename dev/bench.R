# =============================================================================
#  dev/bench.R  --  time and memory per function (a regression baseline)
# =============================================================================
#  Times a representative call of the hot functions and reports the memory each
#  allocates, so a slow-down or a new allocation shows up as a number rather
#  than a hunch.  Uses `bench` when installed (median time + `mem_alloc`); on a
#  bare machine it falls back to `system.time()` and base `Rprofmem()`.
#
#  Needs OligoArrayAux for the thermodynamic rows (see
#  data-raw/install-oligoarrayaux.R); without it those rows are skipped.
#
#  Run:  Rscript dev/bench.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

suppressMessages(pkgload::load_all(root, quiet = TRUE))

has_bench <- requireNamespace("bench", quietly = TRUE)

# --- measurement ------------------------------------------------------------
fmt_bytes <- function(b) {
  if (b < 1024) return(sprintf("%d B", round(b)))
  if (b < 1024^2) return(sprintf("%.1f KB", b / 1024))
  sprintf("%.1f MB", b / 1024^2)
}

# total bytes allocated by one call of fn (cumulative, so freed memory counts)
allocated <- function(fn) {
  tmp <- tempfile()
  Rprofmem(tmp, threshold = 0L)
  on.exit({
    Rprofmem(NULL)
    unlink(tmp)
  }, add = TRUE)
  fn()
  Rprofmem(NULL)
  ln <- grep("^[0-9]", readLines(tmp, warn = FALSE), value = TRUE)
  if (!length(ln)) return(0)
  sum(as.numeric(sub(" .*$", "", ln)))
}

measure <- function(label, fn, iterations = 10L) {
  if (has_bench) {
    res <- bench::mark(fn(), iterations = iterations, check = FALSE)
    ms  <- as.numeric(res$median) * 1000
    mem <- as.numeric(res$mem_alloc)
  } else {
    fn() # warm up
    ms  <- system.time(for (i in seq_len(iterations)) fn())[["elapsed"]] * 1000 / iterations
    mem <- allocated(fn)
  }
  data.frame(fn = label, ms = ms, bytes = mem)
}

# --- a small, self-contained fixture ---------------------------------------
mk_fixture <- function(seed = 7L, n = 2000L) {
  set.seed(seed)
  w <- paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  list(
    genome = Biostrings::DNAStringSet(c(Chr1 = w)),
    primer = substr(w, 500L, 524L)
  )
}

fx <- mk_fixture()
ot <- Sys.which("hybrid-min")

cases <- list(
  "primer_tm"        = function() primer_tm(fx$primer),
  "gc_frac"          = function() gc_frac(fx$primer),
  "snp_template"     = function() snp_template(fx$genome, "Chr1", 1000L, flank = 150L),
  "in_silico_pcr"    = function() in_silico_pcr(fx$primer, as.character(Biostrings::reverseComplement(Biostrings::DNAString(fx$primer))), fx$genome),
  "specificity_seed" = function() specificity_seed(c(P = fx$primer), fx$genome)
)
if (nzchar(ot)) {
  cases$oligo_dg <- function() oligo_dg(fx$primer)
}

# --- run --------------------------------------------------------------------
cat("tidykasp function benchmark (", if (has_bench) "bench" else "system.time + Rprofmem", ")\n", sep = "")
if (!nzchar(ot)) cat("note: OligoArrayAux not on PATH - thermodynamic rows skipped\n")
cat("\n")
rows <- lapply(names(cases), function(nm) measure(nm, cases[[nm]]))
tab  <- do.call(rbind, rows)
for (i in seq_len(nrow(tab))) {
  cat(sprintf("  %-18s %9.3f ms  %10s\n", tab$fn[i], tab$ms[i], fmt_bytes(tab$bytes[i])))
}
cat("\n")
