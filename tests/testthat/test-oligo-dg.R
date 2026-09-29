# OligoArrayAux wrapper: the thermodynamics must be finite, physical, symmetric.

test_that("oligo_dg returns a finite dG for normal oligos", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  s <- "ACGTACGTACGTACGTACGT"
  expect_length(oligo_dg(s, s), 1L)
  expect_true(is.finite(oligo_dg(s, s)))     # duplex mode
  expect_true(is.finite(oligo_dg(s)))        # hairpin mode (hybrid-ss-min)
})

test_that("a strand pairs far more strongly with its reverse complement", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  s  <- "ACGTACGTACGTACGTACGT"
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
  unrelated <- "TTTTTTTTTTTTTTTTTTTT"
  # more negative = stronger; perfect duplex must beat an unrelated strand
  expect_lt(oligo_dg(s, rc), oligo_dg(s, unrelated))
})

test_that("duplex dG is symmetric in its two arguments", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  a <- "ACGTACGTACGTACGTACGT"; b <- "TTTTGGGGCCCCAAAATTTT"
  expect_equal(oligo_dg(a, b), oligo_dg(b, a), tolerance = 1e-6)
})

test_that("a self-complementary oligo folds into a hairpin; poly-A cannot", {
  skip_if_not(nzchar(Sys.which("hybrid-ss-min")), "OligoArrayAux not on PATH")
  palindrome <- "GCGCGCGCGCGCGCGCGCGC"     # perfect self-complementary stem
  poly_a     <- "AAAAAAAAAAAAAAAAAAAA"     # nothing to pair with
  expect_lt(oligo_dg(palindrome), oligo_dg(poly_a))
})
