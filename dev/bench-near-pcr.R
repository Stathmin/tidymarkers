# =============================================================================
#  dev/bench-near-pcr.R  --  cost of the radius-1 near pass of in-silico PCR
# =============================================================================
#  The exact-mode call `in_silico_pcr()` warns when a <= 1-mismatch product
#  exists, and `in_silico_pcr_pool(near = TRUE)` does the same for a pool.  Both
#  pay for one genome-wide radius-1 scan.  This script prints ABSOLUTE totals for
#  every call it times; it does not subtract one call from another to invent a
#  "pass" cost, because the near pool call runs the exact pass first and so its
#  total is always >= the exact total.
#
#  The exact pool call is timed twice: `cold` is the first pool call in the
#  process and pays its setup, `warm` is the one to compare with the near call
#  that is timed after it.  The per-pair lines are absolute too.
#
#  Needs the TAIR10 cache of `ath_genome()` (120 Mb, seven sequences); when it
#  is neither cached nor downloadable the script says so and exits 0.  Run:
#
#    Rscript dev/bench-near-pcr.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

suppressMessages(pkgload::load_all(root, quiet = TRUE))

genome <- tryCatch(ath_genome(), error = function(e) NULL)
if (is.null(genome)) {
  cat("near-pass benchmark skipped: ath_genome() is unavailable (no cache, no network)\n")
  quit(save = "no", status = 0L)
}
chroms <- genome_chroms(genome)
# warm the per-chromosome sequence cache outside the timed regions
invisible(lapply(chroms, function(nm) genome_seq(genome, nm)))

clear_hits <- function() {
  rm(list = ls(.marker_state$hits, all.names = TRUE), envir = .marker_state$hits)
}

timed <- function(expr) {
  gc()
  system.time(expr)[["elapsed"]]
}

f <- "ACGTACGATCGATCGTACGT"
r <- "TTGCAGCTAGCTAGGCTAAC"
pair <- list(near_pair = c(F = f, R = r))

# warm the first-call setup (blocks, masks) outside the timed regions
invisible(in_silico_pcr("GGCATGCATGCTAGCTAGCT", "AATTCCGGAATTCCGGAATT", genome))

cat("tidymarkers near-pass benchmark (TAIR10, ", sum(Biostrings::width(genome)), " bp)\n", sep = "")
cat("\nper-pair path (in_silico_pcr())\n")
clear_hits()
t0 <- timed(out0 <- in_silico_pcr(f, r, genome))
clear_hits()
t1 <- timed(out1 <- in_silico_pcr(f, r, genome, max_mismatch = 1L))
cat(sprintf("  radius 0: %7.2f s  (%d products)\n", t0, nrow(out0)))
cat(sprintf("  radius 1: %7.2f s  (%d products)\n", t1, nrow(out1)))

cat("\npool path (in_silico_pcr_pool(), one pair)\n")
clear_hits()
p0a <- timed(pool0a <- in_silico_pcr_pool(pair, genome))
clear_hits()
p0b <- timed(pool0b <- in_silico_pcr_pool(pair, genome))
clear_hits()
p1 <- timed(pool1 <- in_silico_pcr_pool(pair, genome, near = TRUE))
cat(sprintf("  exact cold (first pool call): %7.2f s  (%d products)\n", p0a, nrow(pool0a[[1L]])))
cat(sprintf("  exact warm:                   %7.2f s  (%d products)\n", p0b, nrow(pool0b[[1L]])))
cat(sprintf("  near (after the warm exact):   %7.2f s  (%d products)\n", p1, nrow(pool1[[1L]])))

cat("\nsame exact verdict across the two paths: ",
  identical(
    as.data.frame(pool0b[[1L]][, c("chrom", "start", "end", "width")]),
    as.data.frame(out0[, c("chrom", "start", "end", "width")])
  ), "\n",
  sep = ""
)
cat("same near verdict across the two paths: ",
  identical(
    as.data.frame(pool1[[1L]][, c("chrom", "start", "end", "width")]),
    as.data.frame(out1[, c("chrom", "start", "end", "width")])
  ), "\n",
  sep = ""
)
