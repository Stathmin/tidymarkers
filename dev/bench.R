# =============================================================================
#  dev/bench.R  --  time and memory per function (a regression baseline)
# =============================================================================
#  Times a representative call of the hot functions and reports the memory each
#  allocates, so a slow-down or a new allocation shows up as a number rather
#  than a hunch.  Backed by `bench` (median time + `mem_alloc`) and `lobstr`
#  (object sizes); both are pinned in renv.lock.
#
#  Compares against dev/bench-baseline.csv when it exists; refresh it with
#  `--update` after an intended change.  `--gate` makes the comparison a gate:
#  exit 1 when a case is more than 2x its baseline, and an explicit `not_run`
#  (exit 0) when the baseline is missing.
#
#  Needs OligoArrayAux for the thermodynamic rows (see
#  data-raw/install-oligoarrayaux.R); without it those rows are skipped.
#
#  Run:  Rscript dev/bench.R [--update|--gate]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
update <- "--update" %in% commandArgs(trailingOnly = TRUE)
gate <- "--gate" %in% commandArgs(trailingOnly = TRUE)
if (update && gate) stop("--update rewrites the baseline; --gate compares against it", call. = FALSE)
# a timing baseline belongs to one host: the same call differs by an order of
# magnitude across CPUs, so a foreign baseline is `not_run`, not a regression
this_host <- Sys.info()[["nodename"]]
this_r <- paste(R.version$major, R.version$minor, sep = ".")

for (pkg in c("bench", "lobstr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("dev/bench.R needs '", pkg, "'; run renv::restore()")
}
suppressMessages(pkgload::load_all(root, quiet = TRUE))
source(file.path(root, "dev", "proc-mem.R"))

# --- fixture ----------------------------------------------------------------
mk_fixture <- function(seed = 7L, n = 2000L) {
  withr::with_seed(seed, {
    w <- paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
    pos <- seq(300L, 1700L, by = 60L)
    ref <- substring(w, pos, pos)
    alt <- vapply(ref, function(b) setdiff(c("A", "C", "G", "T"), b)[1L], character(1))
    list(
      genome = Biostrings::DNAStringSet(c(Chr1 = w)),
      primer = substr(w, 500L, 524L),
      variants = tibble::tibble(chrom = "Chr1", pos = pos, ref = ref, alt = alt)
    )
  })
}
fx <- mk_fixture()

cases <- list(
  primer_tm = function() primer_tm(fx$primer),
  gc_frac = function() gc_frac(fx$primer),
  snp_template = function() snp_template(fx$genome, "Chr1", 1000L, flank = 150L),
  in_silico_pcr = function() in_silico_pcr(fx$primer, as.character(Biostrings::reverseComplement(Biostrings::DNAString(fx$primer))), fx$genome),
  specificity_seed = function() specificity_seed(c(P = fx$primer), fx$genome),
  screen_snp_targets = function() screen_snp_targets(fx$variants, fx$genome)
)
if (nzchar(Sys.which("hybrid-min"))) cases$oligo_dg <- function() oligo_dg(fx$primer)

# --- measure ----------------------------------------------------------------
fmt_bytes <- function(b) {
  if (!is.finite(b)) {
    return("n/a")
  }
  if (b < 1024) sprintf("%d B", round(b)) else if (b < 1024^2) sprintf("%.1f KB", b / 1024) else sprintf("%.1f MB", b / 1024^2)
}
# mem_alloc is NA on an R build without memory profiling, so that case falls
# back to the case's heap delta (retained, not peak).
rows <- lapply(names(cases), function(nm) {
  gc()
  before <- lobstr::mem_used()
  res <- bench::mark(cases[[nm]](), iterations = 20L, check = FALSE, filter_gc = FALSE)
  after <- lobstr::mem_used()
  bytes <- as.numeric(res$mem_alloc)
  fallback <- !is.finite(bytes)
  if (fallback) bytes <- max(0, after - before)
  data.frame(fn = nm, ms = as.numeric(res$median) * 1000, bytes = bytes, fallback = fallback)
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$fn), ]
used_fallback <- tab$fallback
tab <- tab[, c("fn", "ms", "bytes")]

# --- report -----------------------------------------------------------------
baseline <- file.path(root, "dev", "bench-baseline.csv")
cat("tidymarkers function benchmark (bench)\n")
if (!nzchar(Sys.which("hybrid-min"))) cat("note: OligoArrayAux not on PATH - thermodynamic rows skipped\n")
if (gate && !file.exists(baseline)) cat("bench-gate: not_run (dev/bench-baseline.csv missing; run --update once)\n")
if (any(used_fallback)) cat("note: bench mem_alloc is unavailable on this R build - bytes are the case's heap delta (retained, not peak)\n")
cat("\n")

if (update) {
  tab$host <- this_host
  tab$r_version <- this_r
  readr::write_csv(tab, baseline)
}

if (!update && file.exists(baseline)) {
  base <- readr::read_csv(baseline, show_col_types = FALSE)
  foreign <- gate && (!"host" %in% names(base) || !identical(unique(base$host), this_host))
  if (foreign) {
    cat(sprintf(
      "bench-gate: not_run (baseline host %s, this host %s; rerun --update here)\n",
      if ("host" %in% names(base)) paste(unique(base$host), collapse = ",") else "unknown",
      this_host
    ))
    base <- base[0, ]
  }
  tab$time_x <- tab$ms / base$ms[match(tab$fn, base$fn)]
  for (i in seq_len(nrow(tab))) {
    flag <- if (!is.na(tab$time_x[i]) && tab$time_x[i] > 1.5) "   <- slower" else ""
    x <- if (foreign) "   " else sprintf("x%-5.2f", tab$time_x[i])
    cat(sprintf("  %-18s %9.3f ms  %s %10s%s\n", tab$fn[i], tab$ms[i], x, fmt_bytes(tab$bytes[i]), flag))
  }
  if (gate && !foreign) {
    slow <- tab$fn[!is.na(tab$time_x) & tab$time_x > 2]
    if (length(slow)) {
      cat("\nbench-gate: FAIL\n")
      for (fn in slow) cat(sprintf("  - %s is x%.2f the baseline\n", fn, tab$time_x[tab$fn == fn]))
      quit(status = 1L)
    }
    cat("\nbench-gate: ok (every case within 2x of dev/bench-baseline.csv)\n")
  }
} else {
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("  %-18s %9.3f ms  %10s\n", tab$fn[i], tab$ms[i], fmt_bytes(tab$bytes[i])))
  }
  if (update) cat("\nbaseline written to dev/bench-baseline.csv\n")
}

cat(sprintf("peak RSS: %s\n", fmt_peak_rss(proc_peak_rss_kb())))

cat("\nobject sizes (lobstr):\n")
objects <- list(
  genome = fx$genome,
  template = snp_template(fx$genome, "Chr1", 1000L, flank = 150L)
)
for (nm in names(objects)) cat(sprintf("  %-18s %10s\n", nm, format(lobstr::obj_size(objects[[nm]]))))
cat("\n")
