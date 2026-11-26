#!/usr/bin/env Rscript
# =============================================================================
#  dev/bench-index.R  --  one locus, memory vs indexed scan
# =============================================================================
#  Measures the one-locus cost of the two answers to the same question on a
#  synthetic 3 Mb genome (three 1 Mb chromosomes): the in-memory 3'-seed scan
#  and the per-chromosome bowtie index.  Writes the table used in
#  docs/index-perf.md as markdown on stdout.
#
#  Run:  Rscript dev/bench-index.R
# =============================================================================
suppressPackageStartupMessages(pkgload::load_all("."))

set.seed(20261126)
rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
site <- rand(24L)
# one planted locus, so both strategies have exactly one place to find
genome <- Biostrings::DNAStringSet(c(
  chr1 = paste0(rand(5e5), site, rand(499976L)),
  chr2 = rand(1e6),
  chr3 = rand(1e6)
))
fa <- tempfile(fileext = ".fa")
Biostrings::writeXStringSet(genome, fa)
primers <- c(FWD = site, REV = as.character(Biostrings::reverseComplement(
  Biostrings::DNAString(site)
)))

seconds <- function(fun, times = 3L) {
  stats::median(vapply(seq_len(times), function(i) {
    unname(system.time(fun())[["elapsed"]])
  }, numeric(1)))
}

mem_idx <- genome_index(genome, "memory")
mem_s <- seconds(function() {
  specificity_of(primers, mem_idx, k = SEED_K, max_mismatch = 0L)
})
mem_hits <- specificity_of(primers, mem_idx, k = SEED_K, max_mismatch = 0L)$hits

dir <- tempfile("bench-index-")
fs::dir_create(dir)
bowtie_idx <- genome_index(marker_genome(fa), "bowtie",
  chroms = c("chr1", "chr2", "chr3"), dir = dir
)
# the build is one-off and deliberately inside the measured call
build_s <- unname(system.time(
  specificity_of(primers, bowtie_idx, max_mismatch = 0L)
)[["elapsed"]])
query_s <- seconds(function() {
  specificity_of(primers, bowtie_idx, max_mismatch = 0L)
})
idx_hits <- specificity_of(primers, bowtie_idx, max_mismatch = 0L)$hits

kb <- function(x) sprintf("%.1f", as.numeric(object.size(x)) / 1024^2)
rows <- tibble::tibble(
  strategy = c("memory (3'-seed)", "bowtie, per-chromosome (whole-primer)"),
  shards = c(1L, nrow(bowtie_idx$manifest)),
  build_s = c(NA_real_, round(build_s, 3)),
  query_s = round(c(mem_s, query_s), 4),
  hits = c(mem_hits[[1L]], idx_hits[[1L]]),
  index_object_MB = c(kb(mem_idx), kb(bowtie_idx))
)
cat("| strategy | shards | build (s) | query median (s) | hits | index object (MB) |\n")
cat("| --- | --- | --- | --- | --- | --- |\n")
for (i in seq_len(nrow(rows))) {
  cat(sprintf(
    "| %s | %d | %s | %.4f | %d | %s |\n",
    rows$strategy[[i]], rows$shards[[i]],
    if (is.na(rows$build_s[[i]])) "-" else sprintf("%.3f", rows$build_s[[i]]),
    rows$query_s[[i]], rows$hits[[i]], rows$index_object_MB[[i]]
  ))
}
cat("\n")
cat("genome bases:", sum(Biostrings::width(genome)), "\n")
cat("genome fasta MB:", sprintf("%.1f", file.size(fa) / 1024^2), "\n")
cat("manifest rows:", nrow(index_manifest(bowtie_idx)), "\n")
cat("index files MB:", sprintf(
  "%.1f",
  sum(as.numeric(fs::file_size(Sys.glob(file.path(dir, "*"))))) / 1024^2
), "\n")
