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
