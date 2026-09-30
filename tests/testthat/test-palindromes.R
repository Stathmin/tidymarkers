# Inverted-repeat / palindrome screen (Biostrings::findPalindromes).

test_that("has_palindrome flags inverted repeats, not good primers", {
  expect_true(has_palindrome("ACGTACGTACGTACGT"))
  expect_false(has_palindrome("CATGGACCTAGCATTGACGTA"))
  expect_false(has_palindrome("GCTTAGCCATGGTCAATGCTA"))
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
