# =============================================================================
#  dev/bench-near-pcr.R  --  cost of the radius-1 near pass of in-silico PCR
# =============================================================================
#  The exact-mode call `in_silico_pcr()` warns when a <= 1-mismatch product
#  exists, and `in_silico_pcr_pool(near = TRUE)` does the same for a pool.  Both
#  pay for one genome-wide radius-1 scan.  This script times that scan against
#  the radius-0 call it accompanies, so the residual cost is a number:
#
#    per-pair near pass   = in_silico_pcr(r = 1) - in_silico_pcr(r = 0)
#    pool near pass       = in_silico_pcr_pool(near = TRUE) - (near = FALSE)
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
cat(sprintf("  near pass: %6.2f s\n", max(0, t1 - t0)))

cat("\npool path (in_silico_pcr_pool(), one pair)\n")
clear_hits()
p0 <- timed(pool0 <- in_silico_pcr_pool(pair, genome))
clear_hits()
p1 <- timed(pool1 <- in_silico_pcr_pool(pair, genome, near = TRUE))
cat(sprintf("  exact:     %7.2f s  (%d products)\n", p0, nrow(pool0[[1L]])))
cat(sprintf("  near:      %7.2f s  (%d products)\n", p1, nrow(pool1[[1L]])))
cat(sprintf("  near pass: %6.2f s\n", max(0, p1 - p0)))

cat("\nsame exact verdict across the two paths: ",
  identical(
    as.data.frame(pool0[[1L]][, c("chrom", "start", "end", "width")]),
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
