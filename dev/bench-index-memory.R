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
json <- NULL
if ("--json" %in% args) {
  j <- match("--json", args)
  json <- args[[j + 1L]]
  args <- args[-c(j, j + 1L)]
}
mode <- args[[1L]]
fa <- args[[2L]]
if (identical(mode, "memory")) {
  g <- Biostrings::readDNAStringSet(fa)
  bases <- sum(Biostrings::width(g))
} else if (identical(mode, "shard")) {
  genome <- marker_genome(fa)
  chroms <- if (length(args) > 2L) strsplit(args[[3L]], ",")[[1L]] else names(genome$lengths)
  bases <- sum(as.numeric(genome$lengths[names(genome$lengths) %in% chroms]))
  dir <- tempfile("bench-index-memory-")
  fs::dir_create(dir)
  index <- genome_index(genome, "bowtie", chroms = chroms, dir = dir)
  specificity_of(c(AS = "ACGTACGTACGTACGT"), index, k = SEED_K, max_mismatch = 0L)
} else {
  stop("mode must be `memory` or `shard`", call. = FALSE)
}
peak <- peak_mb()
cat(sprintf("%s\tpeak_rss_mb=%.1f\tfasta_mb=%.1f\n", mode, peak, file.size(fa) / 1024^2))
if (!is.null(json)) {
  record <- list(
    mode = mode, label = fs::path_file(fa), bases = bases,
    peak_rss_mb = peak, fasta_mb = file.size(fa) / 1024^2,
    chroms = if (mode == "shard") length(chroms) else NA_integer_
  )
  old <- if (fs::file_exists(json)) jsonlite::fromJSON(json, simplifyVector = FALSE) else list()
  writeLines(jsonlite::toJSON(c(old, list(record)), auto_unbox = TRUE, pretty = TRUE), json)
}

