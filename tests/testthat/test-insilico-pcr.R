# Virtual PCR with Biostrings::matchProbePair (no external tools needed).

test_that("in_silico_pcr finds the expected product and nothing else", {
  genome <- Biostrings::DNAStringSet(
    c(chr1 = paste0(strrep("A", 30L), "ACGTACGTACGTACGTACGT", strrep("T", 30L))))
  pr <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$width, 20L)
})

test_that("in_silico_pcr reports no product when a primer is absent", {
  genome <- Biostrings::DNAStringSet(c(chr1 = strrep("A", 60L)))
  expect_equal(nrow(in_silico_pcr("ACGTACGTACGTACGTACGT", "TTTTGGGGCCCCAAAATTTT", genome)), 0L)
})

test_that("in_silico_pcr finds off-target products too", {
  amp <- "GGGGGGGGGGACGTACGTACGTACGTACGTGGGGGGGGGG"
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(amp, strrep("C", 40L), amp)))
  pr <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
  expect_equal(nrow(pr), 2L)            # both copies amplified
})
