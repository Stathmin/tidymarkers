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

test_that("has_palindrome keeps its plain-ACGT answers", {
  expect_identical(has_palindrome(c("ACGTACGT", "AACCGG")), c(TRUE, FALSE))
  expect_true(has_palindrome("ACGTAAAACGT"))
  expect_true(has_palindrome("AACCGGTT"))
  # the arms do not pair at max_mismatch = 0 (checked against the concrete
  # expansion, which finds no member with an inverted repeat)
  expect_false(has_palindrome("AACGTAAAATGCA"))
  expect_false(has_palindrome("ACGTTTTTTTACGG"))
  expect_false(has_palindrome("CATGGACCTAGCATTGACGTA"))
})

test_that("has_palindrome treats IUPAC codes as wildcards", {
  # the N arms can be completed to ACGT...ACGT; findPalindromes() says FALSE
  # because it only pairs N with N
  expect_true(has_palindrome("ACGTNNACGT"))
  # arms ACGT / ACGT with `RR` = {A,G} in the right arm: one concrete member
  # (`ACGTACGT`) is itself the palindrome, so this is TRUE
  expect_true(has_palindrome("ACGTRRACGT"))
  expect_true(has_palindrome("ACGTRGACGT"))
  # lowercase degenerate inputs: Y = {C,T} faces the complement of `g`/`c` in
  # some member, so `acgtnynacgt` completes (member `ACGTACAACGT`)
  expect_true(has_palindrome("acgtncnacgt"))
  expect_true(has_palindrome("acgtnynacgt"))
  expect_false(has_palindrome("acgtnnacgt", min_armlength = 6L))
  # arms ACGT / ACGT with a 4-base loop of N
  expect_true(has_palindrome("ACGTAANNNNACGT"))
})

test_that("has_palindrome wildcard answers match the concrete mixture", {
  # a degenerate oligo is a mixture of concrete oligos: the wildcard answer is
  # "some member of the mixture has an inverted repeat", nothing wider
  expand <- function(s) {
    # one concrete oligo per base choice, as a full grid: a multi-base code
    # must fill exactly one position, so `Reduce(paste0, ...)` (which joins
    # element-wise) would silently shorten the sequence
    parts <- lapply(strsplit(s, "")[[1L]], function(b) {
      bases <- IUPAC_BASES[[toupper(b)]]
      if (is.null(bases)) b else bases
    })
    do.call(paste0, expand.grid(parts, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE))
  }
  oligos <- c(
    "ACGTNNACGT", "ACGTRRACGT", "ACGTRGACGT", "ACGYNNACGR",
    "CANNTG", "GTNACN", "ACGSNNACGS", "ACGTRNNACGT"
  )
  # the mixture must be small enough to enumerate here
  expect_lte(max(vapply(oligos, function(o) length(expand(o)), integer(1))), 64L)
  oracle <- vapply(oligos, function(o) any(has_palindrome(expand(o))), logical(1))
  expect_identical(has_palindrome(oligos), unname(oracle))
})

test_that("has_palindrome survives NA and keeps its shape", {
  r <- has_palindrome(c("ACGTACGT", NA_character_, "NACGT"))
  expect_type(r, "logical")
  expect_length(r, 3L)
  expect_true(r[1L])
  expect_false(r[3L])
  expect_length(has_palindrome(character()), 0L)
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
