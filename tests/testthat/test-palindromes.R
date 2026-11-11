# Inverted-repeat / palindrome screen (Biostrings::findPalindromes).

test_that("has_palindrome flags inverted repeats, not good primers", {
  expect_true(has_palindrome("ACGTACGTACGTACGT"))
  expect_false(has_palindrome("CATGGACCTAGCATTGACGTA"))
  expect_false(has_palindrome("GCTTAGCCATGGTCAATGCTA"))
})

test_that("gc_frac counts GC over every base it was given", {
  # the denominator is `nchar()`, which is what `as.prob = TRUE` divided by: a
  # `+`/`-`/`.` counts in it but is never GC
  expect_equal(gc_frac(c("ACGT", "acgt", "SSSS", "NNNN", "GC", "ATAT")), c(0.5, 0.5, 0, 0, 1, 0))
  expect_equal(gc_frac("ACGT-"), 0.4)
  expect_true(is.nan(gc_frac("")))
  expect_null(names(gc_frac(c(a = "ACGT"))))
  expect_error(gc_frac("ACXT"), "must be DNA")
  expect_error(gc_frac(NA_character_), "must be DNA")
})

test_that("has_palindrome fast path agrees with findPalindromes", {
  # the plain-ACGT fast path (`palindromic_arms()`) must answer the same
  # question as the Biostrings screen it replaces
  set.seed(42)
  s <- vapply(1:200, function(i) paste(sample(c("A", "C", "G", "T"), 24, TRUE), collapse = ""), "")
  ref <- vapply(s, function(x) {
    length(Biostrings::findPalindromes(Biostrings::DNAString(x),
      min.armlength = 4L, max.looplength = 6L, max.mismatch = 0L
    )) > 0L
  }, logical(1))
  expect_identical(has_palindrome(s), unname(ref))
  # the loop bound is the fast path's own: arms ACGT/ACGT with a loop of 8
  expect_true(has_palindrome("ACGTAAAAAAAAACGT", max_looplength = 8L))
  expect_false(has_palindrome("ACGTAAAAAAAAACGT", max_looplength = 6L))
  # a degenerate code leaves the fast path and is still screened
  expect_true(has_palindrome("ACGTNNNNACGT"))
  # `map_int(...) > 0L` dropped names; the vector shape must not change
  expect_null(names(has_palindrome(c(a = "ACGTACGTACGTACGT"))))
})

test_that("find_palindromes returns a tidy hit table per input", {
  h <- find_palindromes(c(good = "CATGGACCTAGCATTGACGTA", bad = "ACGTACGTACGTACGT"))
  expect_length(h, 2L)
  expect_true(all(c("start", "end", "width") %in% names(h[[1]])))
  expect_equal(nrow(h[["good"]]), 0L)
  expect_gt(nrow(h[["bad"]]), 0L)
  # every reported view is inside the sequence
  expect_true(all(h[["bad"]]$end <= 16L))
})
