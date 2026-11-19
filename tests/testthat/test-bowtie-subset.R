# per-chromosome (pangenome) bowtie index: one index per chromosome, hits summed.

# the two-strand oracle from test-bowtie.R, kept here so the file stands alone
oracle_chrom_hits <- function(primer, genome) {
  rc <- Biostrings::reverseComplement(Biostrings::DNAString(primer))
  sum(vapply(genome, function(seq) {
    Biostrings::countPattern(primer, seq, fixed = TRUE) +
      Biostrings::countPattern(rc, seq, fixed = TRUE)
  }, integer(1)))
}

# `motif` sits in chr1 only, `repeated` in both; both are non-palindromic
two_chrom_genome <- function() {
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  motif <- rand(25L)
  repeated <- rand(25L)
  list(
    genome = Biostrings::DNAStringSet(c(
      chr1 = paste0(rand(500L), motif, rand(500L), repeated, rand(500L)),
      chr2 = paste0(rand(500L), repeated, rand(500L))
    )),
    motif = motif,
    repeated = repeated
  )
}

test_that("the monolithic index is one unprefixed basename and is reused", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(11)
  g <- two_chrom_genome()
  base <- fs::path(withr::local_tempdir(), "bw")
  idx <- marker_bowtie_index(g$genome, basename = base)
  expect_length(idx, 1L)
  expect_equal(idx, as.character(base))
  expect_null(attr(idx, "chroms"))
  expect_true(fs::file_exists(paste0(idx, ".1.ebwt")))
  # a second call reuses the built index instead of rebuilding it
  expect_equal(marker_bowtie_index(g$genome, basename = base), idx)
})

test_that("chroms builds one index per chromosome and reuses each of them", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(11)
  g <- two_chrom_genome()
  base <- fs::path(withr::local_tempdir(), "bw")
  idx <- marker_bowtie_index(g$genome, basename = base, chroms = c("chr1", "chr2"))
  expect_s3_class(idx, "marker_bowtie_index_set")
  expect_length(idx, 2L)
  expect_named(idx, c("chr1", "chr2"))
  expect_equal(unname(as.character(idx)), paste0(as.character(base), "_", c("chr1", "chr2")))
  expect_equal(attr(idx, "chroms"), c("chr1", "chr2"))
  expect_true(all(fs::file_exists(paste0(idx, ".1.ebwt"))))
  expect_equal(marker_bowtie_index(g$genome, basename = base, chroms = c("chr1", "chr2")), idx)
})

test_that("hits over the set equal hits over the monolithic index", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(11)
  g <- two_chrom_genome()
  primers <- c(UNIQUE = g$motif, REPEATED = g$repeated)
  base <- fs::path(withr::local_tempdir(), "bw")
  mono <- specificity_bowtie(primers, marker_bowtie_index(g$genome, basename = base))
  set_res <- specificity_bowtie(
    primers,
    marker_bowtie_index(g$genome, basename = base, chroms = c("chr1", "chr2"))
  )
  expect_equal(set_res$hits, mono$hits)
  expect_equal(set_res$primer, mono$primer)
  expect_equal(set_res$engine, mono$engine)
  expect_equal(attr(set_res, "engine_stamp"), attr(mono, "engine_stamp"))
  expect_equal(set_res$hits, c(oracle_chrom_hits(g$motif, g$genome), 2L))
})

test_that("a one-chromosome subset counts that chromosome alone", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(11)
  g <- two_chrom_genome()
  primers <- c(UNIQUE = g$motif, REPEATED = g$repeated)
  base <- fs::path(withr::local_tempdir(), "bw")
  idx <- marker_bowtie_index(g$genome, basename = base, chroms = c("chr1", "chr2"))
  only1 <- specificity_bowtie(primers, idx["chr1"])
  expect_equal(only1$primer, names(primers))
  # the repeated motif is in both chromosomes; one chromosome alone sees one copy
  expect_equal(only1$hits[only1$primer == "REPEATED"], 1L)
  expect_equal(only1$hits[only1$primer == "UNIQUE"], 1L)
  # a plain subset of the set needs no new argument
  expect_equal(specificity_bowtie(primers, idx[c("chr1", "chr2")])$hits, c(1L, 2L))
})

test_that("a FASTA-backed genome slices by its mapped chromosome names", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(12)
  g <- two_chrom_genome()
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(g$genome, fa)
  genome <- marker_genome(fa, chrom_map = c(chr1 = "Chr1", chr2 = "Chr2"))
  idx <- marker_bowtie_index(
    genome,
    basename = fs::path(withr::local_tempdir(), "bw"), chroms = c("Chr1", "Chr2")
  )
  expect_named(idx, c("Chr1", "Chr2"))
  res <- specificity_bowtie(c(REPEATED = g$repeated), idx)
  expect_equal(res$hits, 2L)
})

test_that("a chromosome the genome does not have aborts and names it", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(11)
  g <- two_chrom_genome()
  expect_error(
    marker_bowtie_index(
      g$genome,
      basename = fs::path(withr::local_tempdir(), "bw"), chroms = "chrX"
    ),
    "chrX"
  )
})

test_that("an empty chroms vector aborts and names the argument", {
  skip_if_not_full()
  skip_if_not_installed("Rbowtie")
  withr::local_seed(11)
  g <- two_chrom_genome()
  expect_error(
    marker_bowtie_index(
      g$genome,
      basename = fs::path(withr::local_tempdir(), "bw"), chroms = character(0)
    ),
    "chroms"
  )
})
