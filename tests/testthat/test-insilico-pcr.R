# Virtual PCR with Biostrings::matchProbePair (no external tools needed).

test_that("in_silico_pcr finds the expected product and nothing else", {
  genome <- Biostrings::DNAStringSet(
    c(chr1 = paste0(strrep("A", 30L), "ACGTACGTACGTACGTACGT", strrep("T", 30L)))
  )
  pr <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$width, 20L)
})

test_that("in_silico_pcr reports no product when a primer is absent", {
  genome <- Biostrings::DNAStringSet(c(chr1 = strrep("A", 60L)))
  expect_equal(nrow(in_silico_pcr("ACGTACGTACGTACGTACGT", "TTTTGGGGCCCCAAAATTTT", genome)), 0L)
})

test_that("in_silico_pcr is invariant to the parallel backend", {
  skip_if_not_installed("BiocParallel")
  amp <- "GGGGGGGGGGACGTACGTACGTACGTACGTGGGGGGGGGG"
  genome <- Biostrings::DNAStringSet(c(
    chr1 = paste0(amp, strrep("C", 40L), amp),
    chr2 = amp
  ))
  serial <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
  par <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_equal(par, serial)
  skip_on_os("windows")
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  mc <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome,
    BPPARAM = bp
  )
  expect_equal(mc, serial)
})

test_that("in_silico_pcr supports degenerate primers", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGGAAAATGCAAAACGT"))
  pr <- in_silico_pcr("ACGS", "TGCA", genome)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$start, 1L)
  expect_equal(pr$end, 12L)
  expect_equal(pr$width, 12L)
})

test_that("in_silico_pcr finds a product with one mismatch in a primer", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGGAAAATGCAAAACGT"))
  # ACGT vs ACGG is one mismatch at the 4th base; the reverse primer stays intact
  pr <- in_silico_pcr("ACGT", "TGCA", genome, max_mismatch = 1L)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$start, 1L)
  expect_equal(pr$end, 12L)
  expect_equal(nrow(in_silico_pcr("ACGT", "TGCA", genome)), 0L)
})

test_that("exact mode warns when a one-mismatch product is hidden", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGGAAAATGCAAAACGT"))
  expect_warning(
    in_silico_pcr("ACGT", "TGCA", genome),
    "near-match"
  )
})

test_that("in_silico_pcr finds off-target products too", {
  amp <- "GGGGGGGGGGACGTACGTACGTACGTACGTGGGGGGGGGG"
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(amp, strrep("C", 40L), amp)))
  pr <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
  # two intra-copy products plus the real cross-copy product (forward in copy 1,
  # reverse in copy 2)
  expect_equal(nrow(pr), 3L)
})

test_that("the specificity engine advertises in-silico PCR without changing the matcher", {
  eng <- marker_engines()$specificity
  expect_true(capabilities(eng, "in_silico_pcr"))
  expect_true(capabilities(eng, "iupac"))
  expect_identical(capabilities(eng, "not_a_capability"), "unknown")
  # metadata only: the matcher keeps its signature (no conditions argument)
  expect_false("conditions" %in% names(formals(in_silico_pcr)))
  expect_identical(eng$pcr, in_silico_pcr)
})
