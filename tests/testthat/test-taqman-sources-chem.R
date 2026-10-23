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
