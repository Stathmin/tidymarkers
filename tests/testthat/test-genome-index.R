# ADR 0005: the genome declares its strategy and the scan names its question.

test_that("genome_index declares the capability matrix per strategy", {
  gs <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGT"))
  expected <- list(
    memory = list(
      seed_scan = TRUE, whole_primer_scan = FALSE,
      mismatch_tolerant = TRUE, max_mismatch = Inf, engine = "Biostrings"
    ),
    bowtie = list(
      seed_scan = TRUE, whole_primer_scan = TRUE,
      mismatch_tolerant = TRUE, max_mismatch = 3, engine = "Rbowtie"
    ),
    blast = list(
      seed_scan = TRUE, whole_primer_scan = TRUE,
      mismatch_tolerant = FALSE, max_mismatch = 0, engine = "BLAST+"
    ),
    fm = list(
      seed_scan = FALSE, whole_primer_scan = TRUE,
      mismatch_tolerant = TRUE, max_mismatch = 2, engine = "Rsubread"
    )
  )
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)
  for (strategy in names(expected)) {
    genome <- if (strategy == "memory") gs else marker_genome(fa)
    idx <- genome_index(genome, strategy)
    expect_identical(idx$strategy, strategy)
    expect_identical(capabilities(idx), expected[[strategy]])
    expect_identical(capabilities(idx, "engine"), expected[[strategy]]$engine)
  }
  expect_error(genome_index("not a genome", "memory"), "DNAStringSet")
})

test_that("memory and bowtie index give the same verdict on a small fixture", {
  skip_if_not_installed("Rbowtie")
  withr::local_seed(7)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  unique <- rand(25L)
  repeated <- rand(20L)
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
    rand(300L), unique, rand(300L), repeated, rand(300L), repeated, rand(300L)
  )))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(genome, fa)
  primers <- c(UNIQUE = unique, TWICE = repeated)

  mem <- specificity_of(primers, genome_index(genome, "memory"),
    k = 15L, max_mismatch = 0L
  )
  idx <- specificity_of(
    primers,
    genome_index(marker_genome(fa), "bowtie", dir = withr::local_tempdir()),
    k = 15L, max_mismatch = 0L
  )

  # the two paths answer different questions, but the locus verdict is the same
  verdict <- function(res) ifelse(res$hits == 1L, "unique", "off-target")
  expect_identical(mem$primer, idx$primer)
  expect_identical(mem$hits, idx$hits)
  expect_identical(verdict(mem), verdict(idx))
  expect_identical(verdict(mem), c("unique", "off-target"))

  # every path declares attr(x, "scan"), and the strategy names it
  expect_identical(attr(mem, "scan"), "3'-seed")
  expect_identical(attr(idx, "scan"), "whole-primer")
  expect_identical(mem$engine, rep("Biostrings", 2L))
  expect_identical(idx$engine, rep("Rbowtie", 2L))
})

test_that("the manifest names every shard and hashes its artefact set", {
  skip_if_not_installed("Rbowtie")
  withr::local_seed(11)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  gs <- Biostrings::DNAStringSet(c(
    chr1 = rand(400L), chr2 = rand(350L), chr3 = rand(300L)
  ))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)
  dir <- withr::local_tempdir()

  idx <- genome_index(marker_genome(fa), "bowtie", chroms = c("chr1", "chr2"), dir = dir)
  expect_identical(names(idx$manifest), INDEX_MANIFEST_COLUMNS)
  expect_identical(idx$manifest$chrom, c("chr1", "chr2"))
  expect_identical(idx$manifest$tool, rep("Rbowtie", 2L))
  expect_false(is.na(idx$manifest$tool_version[[1L]]))
  expect_true(all(is.na(idx$manifest$sha256))) # nothing built by the constructor

  # the object holds paths, never sequence
  expect_false(any(vapply(idx, function(x) inherits(x, "DNAStringSet"), logical(1))))
  # a lazy per-chromosome build keeps no whole chromosome in the query cache
  genome <- marker_genome(fa)
  idx <- genome_index(genome, "bowtie", chroms = c("chr1", "chr2"), dir = dir)
  primer <- as.character(Biostrings::subseq(gs[["chr1"]], 1L, 20L))
  specificity_of(c(AS = primer), idx, k = 15L, max_mismatch = 0L)
  expect_identical(length(ls(genome$cache)), 0L)

  man <- index_manifest(idx)
  expect_true(all(!is.na(man$sha256)))
  # one row per requested chromosome; the measured hash covers the shard files
  files <- sort(Sys.glob(paste0(man$path[[1L]], ".*")))
  parts <- vapply(files, function(f) {
    paste0(basename(f), ":", digest::digest(file = f, algo = "sha256"))
  }, character(1))
  expect_identical(man$sha256[[1L]], digest::digest(paste(parts, collapse = "|"), algo = "sha256"))
  expect_identical(man$bytes[[1L]], sum(as.numeric(fs::file_size(files))))
  expect_false(is.na(man$created[[1L]]))
  expect_identical(man$source_sha256[[1L]], digest::digest(file = fa, algo = "sha256"))
})

test_that("a valid shard is reused without touching its mtime", {
  skip_if_not_installed("Rbowtie")
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  gs <- Biostrings::DNAStringSet(c(chr1 = rand(200L), chr2 = rand(200L)))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)
  dir <- withr::local_tempdir()
  primer <- c(AS = as.character(Biostrings::subseq(gs[["chr1"]], 1L, 20L)))
  build <- function() {
    idx <- genome_index(marker_genome(fa), "bowtie", chroms = c("chr1", "chr2"), dir = dir)
    specificity_of(primer, idx, k = 15L, max_mismatch = 0L)
    idx
  }
  first <- build()
  shard <- paste0(first$manifest$path[[1L]], ".1.ebwt")
  before <- fs::file_info(shard)$modification_time
  build()
  expect_identical(fs::file_info(shard)$modification_time, before)
})

test_that("a missing shard names the chromosome, a changed one names the shard", {
  skip_if_not_installed("Rbowtie")
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  gs <- Biostrings::DNAStringSet(c(chr1 = rand(200L), chr2 = rand(200L)))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)
  dir <- withr::local_tempdir()
  primer <- c(AS = as.character(Biostrings::subseq(gs[["chr1"]], 1L, 20L)))
  build <- function() {
    idx <- genome_index(marker_genome(fa), "bowtie", chroms = c("chr1", "chr2"), dir = dir)
    specificity_of(primer, idx, k = 15L, max_mismatch = 0L)
    # record the built shards, as a second run would
    idx$manifest <- index_manifest(idx)
    idx
  }

  tampered_idx <- build()
  cat("x", file = paste0(tampered_idx$manifest$path[[1L]], ".1.ebwt"), append = TRUE)
  expect_error(
    index_verify(tampered_idx),
    class = "tidymarkers_index_integrity_error"
  )

  unlink(Sys.glob(paste0(tampered_idx$manifest$path[2L], ".*")))
  idx <- genome_index(marker_genome(fa), "bowtie", chroms = c("chr1", "chr2"), dir = dir)
  expect_error(
    index_verify(idx),
    regexp = "chromosome chr2",
    class = "tidymarkers_index_missing_error"
  )
})
