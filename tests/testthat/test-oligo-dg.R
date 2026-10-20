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

test_that("primer_tm() with explicit default conditions matches the pre-seam call", {
  s <- c("ACGTACGTACGTACGTACGT", "GAAGGTGACCAAGTTCATGCT", "TTTTGGCGCGCCGATATATAT")
  # the pre-seam call: the conditions were a global, DMSO was the only argument
  legacy <- primer_tm(s, DMSO = 0)
  expect_equal(primer_tm(s), legacy)
  expect_equal(primer_tm(s, conditions = NULL), legacy)
  expect_equal(primer_tm(s, conditions = reaction_conditions()), legacy)
  expect_equal(
    primer_tm(s, conditions = list(Na = 50, Mg = 1.5, dNTPs = 0.6, DMSO = 0)),
    legacy
  )
})

test_that("primer_tm() takes salt, dNTP and DMSO from the conditions", {
  s <- "GAAGGTGACCAAGTTCATGCT"
  base <- primer_tm(s)
  expect_gt(primer_tm(s, conditions = reaction_conditions(Na = 100)), base)
  expect_lt(primer_tm(s, conditions = reaction_conditions(DMSO = 5)), base)
  expect_false(isTRUE(all.equal(
    primer_tm(s, conditions = reaction_conditions(dNTPs = 0.2)), base
  )))
  # an explicit DMSO still wins over the conditions field
  expect_equal(
    primer_tm(s, DMSO = 2, conditions = reaction_conditions(DMSO = 5)),
    primer_tm(s, DMSO = 2)
  )
})

test_that("primer_tm() temperature is informational (the NN model gives Tm)", {
  s <- "GAAGGTGACCAAGTTCATGCT"
  expect_equal(
    primer_tm(s, conditions = reaction_conditions(temperature = 60)),
    primer_tm(s)
  )
})

test_that("scheme_conditions() encodes the documented reaction temperatures", {
  window <- function(scheme) attr(scheme_conditions(scheme), "temperature_range")
  # PCR-style chemistries: one 60 C point, the anneal/extension temperature the
  # package already uses for pcr_efficiency()
  expect_equal(window("pcr"), c(60, 60))
  for (scheme in c("snp", "indel", "multi", "duplex", "taqman", "hrm")) {
    expect_equal(scheme_conditions(scheme)$temperature, 60)
  }
  # isothermal windows, as recorded in R/marker_lamp.R and R/marker_rpa.R
  expect_equal(window("lamp"), c(60, 65))
  expect_equal(window("rpa"), c(37, 42))
  expect_equal(scheme_conditions("lamp")$temperature, 65)
  expect_equal(scheme_conditions("rpa")$temperature, 42)
  # room-temperature methods, and a CAPS digest is not a 60 C reaction
  expect_equal(scheme_conditions("caps")$temperature, 25)
  expect_equal(window("ambient"), c(25, 25))
})

test_that("scheme_conditions() keeps the historical buffer for every scheme", {
  for (scheme in c("pcr", "snp", "lamp", "rpa", "ambient", "not_a_scheme")) {
    cond <- scheme_conditions(scheme)
    expect_equal(cond[c("Na", "Mg", "dNTPs", "DMSO")], reaction_conditions()[c("Na", "Mg", "dNTPs", "DMSO")])
  }
  # an unregistered id is a PCR-style chemistry, not an error
  expect_equal(scheme_conditions("not_a_scheme")$temperature, 60)
  expect_error(scheme_conditions(c("lamp", "rpa")), "single scheme id")
  expect_error(scheme_conditions(NA_character_), "single scheme id")
})

test_that("RPA is not gated on a PCR-style Tm", {
  # RPA's Tm window is open by construction: the recombinase needs primer length,
  # not an annealing Tm, and its crowded buffer is not this NN model
  expect_equal(tidymarkers:::RPA_NO_TM, c(0, 200))
  expect_equal(attr(scheme_conditions("rpa"), "temperature_range"), tidymarkers:::RPA_TEMP)
  expect_equal(scheme_conditions("rpa")$temperature, max(tidymarkers:::RPA_TEMP))
  # the schemes that run their primers through RPA carry the open Tm window and
  # are flagged informational, so a design can never be failed on that Tm
  expect_true(attr(scheme_conditions("rpa"), "tm_informational"))
  expect_true(attr(scheme_conditions("crispr"), "tm_informational"))
  expect_false(attr(scheme_conditions("pcr"), "tm_informational"))
  expect_false(attr(scheme_conditions("lamp"), "tm_informational"))
  tm <- c(20, 60, 120)
  expect_true(all(dplyr::between(tm, tidymarkers:::RPA_NO_TM[1L], tidymarkers:::RPA_NO_TM[2L])))
})

test_that("the encoded temperatures are the ones the method files already use", {
  # scheme_conditions() repeats the windows as literals (RPA_TEMP collates after
  # marker_core.R); this pins them to the method files so they cannot drift
  expect_equal(attr(scheme_conditions("rpa"), "temperature_range"), tidymarkers:::RPA_TEMP)
  # LAMP runs 60-65 C (R/marker_lamp.R header); it is not an exported constant
  expect_equal(attr(scheme_conditions("lamp"), "temperature_range"), c(60, 65))
  # PCR-style = the anneal/extension temperature already used for pcr_efficiency()
  expect_equal(scheme_conditions("pcr")$temperature, tidymarkers:::PCR_EFF$temp)
  # the buffer is the historical PCR buffer for every family (no invented recipe)
  pcr <- reaction_conditions()
  for (scheme in c("pcr", "lamp", "rpa", "ambient")) {
    cond <- scheme_conditions(scheme)
    expect_equal(cond$Na, pcr$Na)
    expect_equal(cond$Mg, pcr$Mg)
    expect_equal(cond$dNTPs, pcr$dNTPs)
    expect_equal(cond$DMSO, pcr$DMSO)
  }
})

test_that("the thermo engine advertises its capabilities", {
  eng <- marker_engines()$thermo
  expect_true(capabilities(eng, "tm"))
  expect_true(capabilities(eng, "temperature"))
  expect_identical(capabilities(eng, "duplex_dg"), has_oligoarrayaux())
  expect_identical(capabilities(eng, "not_a_capability"), "unknown")
})
