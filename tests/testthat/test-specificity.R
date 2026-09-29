# Specificity seed counting: a known seed must be counted on both strands.

test_that("specificity_seed counts seed occurrences on both strands", {
  genome  <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  primers <- c(P = "TTTTTTAAACCCGGG")          # 3'-terminal 9 nt = AAACCCGGG
  # forward "AAACCCGGG" x2 + reverse-complement "CCCGGGTTT" x2 = 4
  res <- specificity_seed(primers, genome, k = 9L, max_mismatch = 0L)
  expect_equal(res$hits, 4L)
  expect_equal(res$seed, "AAACCCGGG")
})

test_that("specificity_seed returns zero hits for an absent seed", {
  genome  <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  primers <- c(P = "TTTTTTACGACGAC")
  expect_equal(specificity_seed(primers, genome, k = 9L, max_mismatch = 0L)$hits, 0L)
})
