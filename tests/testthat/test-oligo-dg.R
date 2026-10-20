# OligoArrayAux wrapper: the thermodynamics must be finite, physical, symmetric.

test_that("oligo_dg returns a finite dG for normal oligos", {
  skip_if_no_oligoarrayaux()
  s <- "ACGTACGTACGTACGTACGT"
  expect_length(oligo_dg(s, s), 1L)
  expect_true(is.finite(oligo_dg(s, s))) # duplex mode
  expect_true(is.finite(oligo_dg(s))) # hairpin mode (hybrid-ss-min)
})

test_that("a strand pairs far more strongly with its reverse complement", {
  skip_if_no_oligoarrayaux()
  s <- "ACGTACGTACGTACGTACGT"
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
  unrelated <- "TTTTTTTTTTTTTTTTTTTT"
  # more negative = stronger; perfect duplex must beat an unrelated strand
  expect_lt(oligo_dg(s, rc), oligo_dg(s, unrelated))
})

test_that("duplex dG is symmetric in its two arguments", {
  skip_if_no_oligoarrayaux()
  a <- "ACGTACGTACGTACGTACGT"
  b <- "TTTTGGGGCCCCAAAATTTT"
  expect_equal(oligo_dg(a, b), oligo_dg(b, a), tolerance = 1e-6)
})

test_that("a self-complementary oligo folds into a hairpin; poly-A cannot", {
  skip_if_not(nzchar(Sys.which("hybrid-ss-min")), "OligoArrayAux not on PATH")
  palindrome <- "GCGCGCGCGCGCGCGCGCGC" # perfect self-complementary stem
  poly_a <- "AAAAAAAAAAAAAAAAAAAA" # nothing to pair with
  expect_lt(oligo_dg(palindrome), oligo_dg(poly_a))
})

# Conditions / engines seam ---------------------------------------------------

test_that("reaction_conditions() defaults are the historical conditions", {
  cond <- reaction_conditions()
  expect_identical(names(cond), c("Na", "Mg", "dNTPs", "DMSO", "temperature"))
  expect_equal(cond$Na, 50)
  expect_equal(cond$Mg, 1.5)
  expect_equal(cond$dNTPs, 0.6)
  expect_equal(cond$DMSO, 0)
  expect_equal(cond$temperature, 37)
})

test_that("reaction_conditions() validates its fields", {
  expect_error(reaction_conditions(Na = c(50, 100)), "single finite")
  expect_error(reaction_conditions(DMSO = -1), "negative")
  expect_error(reaction_conditions(temperature = NA_real_), "single finite")
})

test_that("oligo_dg() defaults reproduce the explicit 50 mM / 37 C call", {
  skip_if_no_oligoarrayaux()
  s <- "ACGTACGTACGTACGTACGT"
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
  cond <- reaction_conditions()
  expect_equal(oligo_dg(s, rc), oligo_dg(s, rc, na = cond$Na / 1000, temp = cond$temperature))
})

test_that("oligo_dg() dG depends on the reaction temperature and salt", {
  skip_if_no_oligoarrayaux()
  s <- "ACGTACGTACGTACGTACGT"
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
  expect_false(isTRUE(all.equal(oligo_dg(s, rc, temp = 60), oligo_dg(s, rc))))
  expect_false(isTRUE(all.equal(oligo_dg(s, rc, na = 1), oligo_dg(s, rc))))
})

test_that("primer_tm() reproduces the legacy conditions and honours overrides", {
  s <- "ACGTACGTACGTACGTACGT"
  expect_equal(primer_tm(s), primer_tm(s, DMSO = 0))
  expect_equal(primer_tm(s), primer_tm(s, conditions = reaction_conditions()))
  expect_false(isTRUE(all.equal(
    primer_tm(s), primer_tm(s, conditions = reaction_conditions(Na = 200))
  )))
  expect_false(isTRUE(all.equal(
    primer_tm(s), primer_tm(s, conditions = reaction_conditions(DMSO = 5))
  )))
})

test_that("primer_tm() fills a partial conditions list from the defaults", {
  s <- "ACGTACGTACGTACGTACGT"
  expect_equal(
    primer_tm(s, conditions = list(Na = 50)),
    primer_tm(s, conditions = reaction_conditions())
  )
  expect_error(primer_tm(s, conditions = list(salt = 50)), "Unknown reaction condition")
})

test_that("the thermo engine advertises its capabilities", {
  eng <- marker_engines()$thermo
  expect_true(capabilities(eng, "tm"))
  expect_true(capabilities(eng, "temperature"))
  expect_identical(capabilities(eng, "duplex_dg"), has_oligoarrayaux())
  expect_identical(capabilities(eng, "not_a_capability"), "unknown")
})
