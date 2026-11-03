# TaqMan / 5'-nuclease science extracts (znavsmgb2014, apol1_2019). Data-only,
# from inst/extdata; the generators in data-raw/* re-derive every number from
# the fetched JATS before writing.

tq_read <- function(name) {
  utils::read.csv(tm_extdata(name),
    stringsAsFactors = FALSE
  )
}

test_that("znavsmgb2014's Tm ladder is monotone and internally consistent", {
  t <- tq_read("taqman_znavsmgb2014_tm.csv")
  v <- stats::setNames(t$value, t$item)
  tm_c <- v[["probe allele C: nearest-neighbour Tm, unmodified (C)"]]
  tm_t <- v[["probe allele T: nearest-neighbour Tm, unmodified (C)"]]
  tm_z <- v[["Tm with ZNA (spermine) conjugation (C)"]]
  expect_identical(unname(c(tm_c, tm_t)), c(47.5, 49.8))
  expect_identical(tm_z, 67)
  expect_lt(tm_c, tm_t)
  expect_lt(tm_t, tm_z)
  expect_identical(v[["Tm rise from ZNA conjugation (C)"]], tm_z - tm_c)
  # the probe is under the package's 18-25 nt default window
  expect_lt(v[["probe allele C: FAM probe length (nt)"]], 18)
})

test_that("znavsmgb2014's design rules are weak, never outcomes", {
  r <- tq_read("taqman_znavsmgb2014_rules.csv")
  expect_identical(nrow(r), 3L)
  expect_true(all(r$evidence_tier == "weak"))
  expect_true(any(grepl("Longer probes are less sensitive", r$reported)))
})

test_that("apol1_2019 records a real assay failure with its own arithmetic", {
  f <- tq_read("taqman_apol1_2019_failure.csv")
  expect_identical(nrow(f), 8L)
  v <- stats::setNames(f$freq_pct, f$genotype)
  expect_identical(sum(v), 99.9)
  expect_identical(v[["G2,G2"]], 1.4)
  expect_identical(v[["G1GM,G2"]], 5.9)
  expect_identical(round(v[["G2,G2"]] + v[["G1GM,G2"]], 1), 7.3)
  expect_true(any(grepl("No amplification", f$predesigned_call)))
  expect_true(all(f$n_samples == 1110L))
  # a failure mode, explicitly not a mismatch discrimination
  expect_match(unique(f$outcome_definition), "failure MODE")
})

test_that("kutyavin2000 measures an MGB Tm gain and a match-mismatch dTm", {
  t <- tq_read("taqman_kutyavin2000_mgb_tm.csv")
  v <- stats::setNames(t$value, t$item)
  # same 15mer with and without 3'-MGB: 42 -> 60 C, the gain is the paper's 18
  expect_identical(v[["ODN1 15mer, unmodified, Tm (match)"]], 42)
  expect_identical(v[["ODN1 15mer, 3'-MGB, Tm (match)"]], 60)
  expect_identical(
    v[["ODN1 15mer, Tm gain from MGB"]],
    v[["ODN1 15mer, 3'-MGB, Tm (match)"]] -
      v[["ODN1 15mer, unmodified, Tm (match)"]]
  )
  # a 12mer + MGB matches a 27mer without one, so MGB buys ~15 nt
  expect_lte(abs(
    v[["12mer MGB probe, Tm (match)"]] -
      v[["no-MGB 27mer probe, Tm (match)"]]
  ), 1)
  # the measured discriminations: 20 C with MGB, 4-6 C without
  expect_identical(v[["12mer MGB probe, match-minus-mismatch dTm"]], 20)
  expect_lt(
    v[["no-MGB probes, match-minus-mismatch dTm, high end"]],
    v[["12mer MGB probe, match-minus-mismatch dTm"]]
  )
  # position dependence: the same T/G mismatch gains only when under the MGB
  expect_gt(
    v[["probe #7 T/G mismatch at 7th base from 3' end, MGB dTm"]],
    v[["probe #7 T/G mismatch at 7th base from 3' end, no-MGB dTm"]]
  )
  # the short probes are below the package's no-MGB 18-25 nt default
  expect_lt(v[["MGB probe length, lower bound"]], 18)
  # recorded as a measurement, and it names what it does not measure
  expect_true(all(t$evidence_tier == "moderate"))
  expect_false(any(c("dct", "delta_ct", "endpoint_ratio") %in% names(t)))
  expect_match(unique(t$outcome_definition), "no dCt")
})
