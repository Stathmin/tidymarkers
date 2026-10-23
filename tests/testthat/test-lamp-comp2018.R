# Competitive SNP-LAMP probes (Cell Rep Methods 2022): endpoint fluorescence,
# not kinetics. Data-only, from inst/extdata.

cp_read <- function(name) {
  utils::read.csv(tm_extdata(name),
    stringsAsFactors = FALSE
  )
}

test_that("the three comp2018 files exist with required columns", {
  for (f in c(
    "lamp_comp2018_separation.csv", "lamp_comp2018_design_space.csv",
    "lamp_comp2018_durations.csv"
  )) {
    d <- cp_read(f)
    expect_gt(nrow(d), 0L)
    expect_true(all(c("outcome_definition", "evidence_tier", "source_file") %in% names(d)))
  }
})

test_that("the separation is an endpoint ratio and says so", {
  s <- cp_read("lamp_comp2018_separation.csv")
  expect_identical(s$magnitude[s$comparison == "SNP vs non-SNP (synthetic oligos, ACTB c.776A>C)"], 2.7)
  expect_true(all(grepl("endpoint fluorescence", s$signal_type)))
  expect_match(unique(s$outcome_definition), "NO time-to-positive")
  expect_false(any(c("time_to_positive_min", "ct") %in% names(s)))
})

test_that("the design-space statistics are internally consistent", {
  d <- cp_read("lamp_comp2018_design_space.csv")
  v <- stats::setNames(d$value, d$item)
  expect_identical(v[["random probe sets simulated"]], 10000)
  expect_lt(v[["fraction reaching >70% signal and >100x specificity"]], 1)
  expect_gt(v[["median SNP/WT specificity over 10 GA runs"]], 1)
  expect_gt(v[["SD of that specificity"]], 0)
  expect_lte(v[["fraction of maximum fluorescence achieved"]], 100)
  # the one-pot probe cap is a real constraint, not a preference
  expect_gt(v[["signal:background at the cap"]], 1)
})

test_that("durations are 30/60/75 min and 'Anneal only' rows carry none", {
  d <- cp_read("lamp_comp2018_durations.csv")
  expect_identical(nrow(d), 8L)
  expect_setequal(stats::na.omit(d$minutes), c(30L, 60L, 75L))
  expect_true(all(is.na(d$minutes[d$format == "Anneal only"])))
  expect_true(all(!is.na(d$minutes[d$format == "One-pot"])))
  expect_true(all(d$format %in% c("One-pot", "Two-step", "Anneal only")))
  # the two-step SARS-CoV-2 experiment is the fastest
  expect_identical(d$minutes[d$format == "Two-step"], 30L)
})
