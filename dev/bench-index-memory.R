#!/usr/bin/env Rscript
# =============================================================================
#  dev/bench-index-memory.R  --  peak RSS of the two genome-access strategies
# =============================================================================
#  Measures the process high-water mark (VmHWM) for
#    memory <fa>            read the whole FASTA into a DNAStringSet
#    shard  <fa> [chroms]   build/query per-chromosome bowtie shards
#  The point is the peak, not the wall time: a per-chromosome shard holds one
#  chromosome, so a 16 Gb fixpur is a disk problem, not a 21.5 GB RAM problem.
#
#  Run:  Rscript dev/bench-index-memory.R shard .tmp/tair10.fa Chr1
#        Rscript dev/bench-index-memory.R memory .tmp/tair10.fa
# =============================================================================
suppressPackageStartupMessages(pkgload::load_all("."))

peak_mb <- function() {
  status <- readLines("/proc/self/status")
  as.numeric(sub("\\D+(\\d+).*", "\\1", status[grepl("^VmHWM:", status)])) / 1024
}

args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1L]]
fa <- args[[2L]]
if (identical(mode, "memory")) {
  Biostrings::readDNAStringSet(fa)
} else if (identical(mode, "shard")) {
  genome <- marker_genome(fa)
  chroms <- if (length(args) > 2L) strsplit(args[[3L]], ",")[[1L]] else names(genome$lengths)
  dir <- tempfile("bench-index-memory-")
  fs::dir_create(dir)
  index <- genome_index(genome, "bowtie", chroms = chroms, dir = dir)
  specificity_of(c(AS = "ACGTACGTACGTACGT"), index, k = SEED_K, max_mismatch = 0L)
} else {
  stop("mode must be `memory` or `shard`", call. = FALSE)
}
cat(sprintf("%s\tpeak_rss_mb=%.1f\tfasta_mb=%.1f\n", mode, peak_mb(), file.size(fa) / 1024^2))
