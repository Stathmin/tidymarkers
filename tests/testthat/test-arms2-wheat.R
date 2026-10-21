# ARMS/AS-PCR family, batch 2: wheat allele-specific marker design funnel.
#
# The wheat paper (Bhoite et al. 2023, Plant Direct 7:e520) is included as a
# CROP OUTCOME RATE, with the assay labelled KASP/PACE. These tests pin:
#   * the arithmetic the paper itself closes (39 - 7 = 32);
#   * that 18 designed markers is corroborated by the supplement, not just prose;
#   * that the 8-vs-7 polymorphic disagreement is PRESERVED, not smoothed;
#   * that nothing can later be relabelled as ARMS.

wheat_dir <- function() test_path("../../inst/extdata")

funnel_csv <- function() {
  path <- file.path(wheat_dir(), "arms_wheat_as2023_funnel.csv")
  skip_if_not(file.exists(path), "wheat funnel CSV not built")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

outcomes_csv <- function() {
  path <- file.path(wheat_dir(), "arms_wheat_as2023_outcomes.csv")
  skip_if_not(file.exists(path), "wheat outcomes CSV not built")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

test_that("the reported funnel closes arithmetically", {
  f <- funnel_csv()
  n <- stats::setNames(f$n, f$stage)
  expect_equal(n[["significant SNPs"]] - n[["unmappable"]], n[["designable"]])
  expect_equal(n[["significant SNPs"]], 39)
  expect_equal(n[["unmappable"]], 7)
  expect_equal(n[["designable"]], 32)
  expect_equal(n[["designed markers"]], 18)
})

test_that("the supplement corroborates 18 designed markers", {
  o <- outcomes_csv()
  expect_equal(nrow(o), 18L)
  expect_equal(anyDuplicated(o$marker), 0L)
  expect_equal(length(unique(o$marker)), 18L)
})

test_that("the 8-vs-7 polymorphic disagreement is preserved, not averaged", {
  f <- funnel_csv()
  n <- stats::setNames(f$n, f$stage)
  expect_equal(n[["polymorphic (abstract)"]], 8)
  expect_equal(n[["polymorphic (named in 2.3)"]], 7)
  # The gap must stay visible: it is a disagreement in the source.
  expect_false(n[["polymorphic (abstract)"]] == n[["polymorphic (named in 2.3)"]])
  # and the named ones must be a subset of the designed ones, counted in the CSV
  o <- outcomes_csv()
  expect_equal(sum(o$named_polymorphic_in_text), n[["polymorphic (named in 2.3)"]])
})

test_that("every polymorphic count is a subset of the designed markers", {
  f <- funnel_csv()
  n <- stats::setNames(f$n, f$stage)
  expect_lte(n[["polymorphic (abstract)"]], n[["designed markers"]])
  expect_lte(n[["polymorphic (named in 2.3)"]], n[["designed markers"]])
  expect_lte(n[["designed markers"]], n[["designable"]])
})

test_that("the assay is labelled KASP/PACE and never ARMS", {
  o <- outcomes_csv()
  expect_true(all(grepl("KASP", o$assay, fixed = TRUE)))
  expect_false(any(grepl("ARMS", o$assay, ignore.case = TRUE)))
  expect_true(all(grepl("NOT ARMS", o$outcome_definition, fixed = TRUE)))
})
