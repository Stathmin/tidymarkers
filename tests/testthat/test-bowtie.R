# bowtie1 genome scan (Rbowtie): the indexed path that scales to large genomes.

oracle_hits <- function(primer, genome) {
  # both strands: bowtie aligns the reverse complement as a separate hit
  rc <- Biostrings::reverseComplement(Biostrings::DNAString(primer))
  sum(vapply(genome, function(seq) {
    Biostrings::countPattern(primer, seq, fixed = TRUE) +
      Biostrings::countPattern(rc, seq, fixed = TRUE)
  }, integer(1)))
}

test_that("specificity_bowtie agrees with the two-strand countPattern oracle", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(1)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  # non-palindromic, so the two strands cannot double-count one locus
  motif <- rand(25L)
  repeated <- rand(25L)
  genome <- Biostrings::DNAStringSet(c(
    chr1 = paste0(rand(500L), motif, rand(500L), repeated, rand(500L)),
    chr2 = paste0(rand(600L), repeated, rand(600L))
  ))
  idx <- marker_bowtie_index(genome, basename = fs::path(withr::local_tempdir(), "bw"))
  primers <- c(UNIQUE = motif, REPEATED = repeated)
  res <- specificity_bowtie(primers, idx)
  expect_setequal(names(res), c("primer", "seed", "hits"))
  expect_equal(res$hits, c(oracle_hits(motif, genome), oracle_hits(repeated, genome)))
  expect_equal(res$hits[res$primer == "REPEATED"], 2L)
})

test_that("an absent primer reports zero hits and input order is preserved", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(2)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  genome <- Biostrings::DNAStringSet(c(chr1 = rand(400L)))
  idx <- marker_bowtie_index(genome, basename = fs::path(withr::local_tempdir(), "bw"))
  primers <- c(Z = "AAAAAAAAAAAAAAAAAAAAAAAAA", A = rand(25L))
  res <- specificity_bowtie(primers, idx)
  # unaligned reads appear in SAM with FLAG 4; counting them would report 1
  expect_equal(res$primer, names(primers))
  expect_equal(res$hits[res$primer == "Z"], 0L)
})

test_that("specificity_of routes a bowtie-index genome to specificity_bowtie", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(3)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  motif <- rand(25L)
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">chr1", paste0(rand(300L), motif, rand(300L))), fa)
  genome <- marker_genome(fa)
  idx <- marker_bowtie_index(genome, basename = fs::path(withr::local_tempdir(), "bw"))
  genome$bowtie_index <- idx
  res <- specificity_of(c(UNIQUE = motif), genome)
  expect_equal(res$hits, 1L)
  # the bowtie path scans the whole primer; a 3'-seed request must warn
  expect_warning(
    specificity_of(c(UNIQUE = motif), genome, k = 12L),
    "whole primer"
  )
})
