# bowtie1 genome scan (Rbowtie): the indexed path that scales to large genomes.

oracle_hits <- function(primer, genome) {
  # both strands: bowtie aligns the reverse complement as a separate hit
  rc <- Biostrings::reverseComplement(Biostrings::DNAString(primer))
  sum(vapply(genome, function(seq) {
    Biostrings::countPattern(primer, seq, fixed = TRUE) +
      Biostrings::countPattern(rc, seq, fixed = TRUE)
  }, integer(1)))
}

test_that("a bowtie build worker count is a positive whole number", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGT"))
  for (x in list(0.5, 0L, -1L, Inf, "4", c(4L, 8L))) {
    expect_error(
      marker_bowtie_index(genome, threads = x),
      "positive whole number of workers"
    )
  }
})

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
  expect_setequal(names(res), c("primer", "seed", "hits", "engine"))
  expect_equal(res$engine, rep("Rbowtie", nrow(res)))
  expect_equal(res$hits, c(oracle_hits(motif, genome), oracle_hits(repeated, genome)))
  expect_equal(res$hits[res$primer == "REPEATED"], 2L)
})

test_that("a bowtie radius is checked before the engine runs", {
  skip_if_not_installed("Rbowtie")
  # the index is never opened: the radius is refused first
  expect_error(
    specificity_bowtie("ACGTACGTACGTACGT", "no-index", max_mismatches = 1.9),
    "non-negative whole number"
  )
  expect_error(
    specificity_bowtie("ACGTACGTACGTACGT", "no-index", max_mismatches = 4L),
    "at most 3"
  )
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
  expect_equal(attr(res, "scan"), "whole-primer")
  expect_equal(res$engine, "Rbowtie")
  expect_warning(
    specificity_of(c(UNIQUE = motif), genome, k = 12L),
    "whole primer"
  )
})

test_that("an indexed genome keeps the AS 3'-seed tie-break (D1)", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(4)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  ref <- rand(25L)
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">chr1", paste0(rand(300L), ref, rand(300L))), fa)
  genome <- marker_genome(fa)
  genome$bowtie_index <- marker_bowtie_index(genome, basename = fs::path(withr::local_tempdir(), "bw"))
  as_cand <- tibble::tibble(
    ref_primer = ref, discrimination = 0.9, tm = 58, dg_3p = -1,
    wobble = FALSE, n_var = 0L, w_var = 0, n_3p = 0L
  )
  res <- marker_pick_as(as_cand, genome,
    min_disc = 0.5, specificity_candidates = 1L, BPPARAM = NULL,
    snp = list(chrom = "chr1", pos = 301L), tmpl = NULL
  )
  # the indexed path used to refuse the seed question, leaving seed_hits NA
  expect_true("seed_hits" %in% names(res$as_cand))
  expect_false(is.na(res$as_cand$seed_hits[[1L]]))
  expect_equal(res$as_cand$seed_hits[[1L]], 1L)
})
