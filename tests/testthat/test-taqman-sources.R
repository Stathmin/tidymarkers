# TaqMan / 5'-nuclease science extracts (wheat2019). Data-only, from
# inst/extdata; data-raw/wheat2019.R re-derives every number from the fetched
# JATS. The other three sources live in the sibling *-citrus / *-chem files.

tq_read <- function(name) {
  utils::read.csv(testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE
  )
}

test_that("the TaqMan extracts all share the common schema", {
  for (f in c(
    "taqman_wheat2019_platform.csv", "taqman_wheat2019_composition.csv",
    "taqman_wheat2019_anova.csv", "taqman_wheat2019_concordance.csv",
    "taqman_citrus2020_calls.csv", "taqman_citrus2020_probes.csv",
    "taqman_citrus2020_amplicons.csv", "taqman_znavsmgb2014_tm.csv",
    "taqman_znavsmgb2014_rules.csv", "taqman_apol1_2019_failure.csv"
  )) {
    d <- tq_read(f)
    expect_gt(nrow(d), 0L)
    expect_true(all(
      c("outcome_definition", "evidence_tier", "source_file") %in% names(d)
    ))
    expect_true(all(nzchar(d$outcome_definition)))
    expect_true(all(d$evidence_tier %in% c("strong", "moderate", "weak")))
  }
})

test_that("wheat2019 is a platform comparison, never a mismatch measurement", {
  p <- tq_read("taqman_wheat2019_platform.csv")
  expect_identical(nrow(p), 3L)
  v <- stats::setNames(p$pct_amplified, p$platform)
  expect_identical(v[["TaqMan"]], 93)
  expect_identical(v[["rhAmp"]], 97)
  # the plate frame: 96 wells = 94 genotypes + 2 NTCs; 2,784 = 96 * 29
  expect_true(all(p$datapoints == 2784L))
  expect_true(all(p$assays_compared == 29L))
  expect_true(all(p$assays_designed_ok <= p$assays_attempted))
  # TaqMan amplified less than each of the other two, and cannot exceed 100%
  expect_true(all(p$pct_amplified > 0 & p$pct_amplified <= 100))
  expect_lt(v[["TaqMan"]], v[["KASP"]])
  expect_lt(v[["TaqMan"]], v[["rhAmp"]])
  expect_true(all(round(100 - p$pct_unamplified, 1) == p$pct_amplified))
  # the extraction must not claim a mismatch readout it does not have
  expect_match(unique(p$outcome_definition), "NOT a mismatch discrimination")
  expect_false(any(c("dct", "delta_ct", "mismatch_position") %in% names(p)))
})

test_that("wheat2019's ANOVA blocks close internally", {
  a <- tq_read("taqman_wheat2019_anova.csv")
  expect_identical(nrow(a), 4L)
  expect_true(all(a$df_tech == 2L))
  # F = tech MS / residual MS (within 5% of the printed 2-dp values); the
  # call-quality block prints an unparsable residual SS ("32.9.4") and so has
  # no recomputed F -- that NA must be the only one
  has <- !is.na(a$F_recomputed)
  expect_identical(sum(!has), 1L)
  expect_match(a$parameter[!has], "Call quality")
  expect_true(all(
    abs(a$F_recomputed[has] - a$F_reported[has]) / a$F_reported[has] < 0.05
  ))
  expect_true(all(a$df_resid > 2L))
  # the separation-angle block is the one with 59 df
  expect_identical(a$df_resid[grepl("Separation angle", a$parameter)], 59L)
})

test_that("wheat2019's composition row set matches its ANOVA label letters", {
  c0 <- tq_read("taqman_wheat2019_composition.csv")
  expect_identical(nrow(c0), 3L)
  # TaqMan's letters: angle 'a' (the only a), compactness 'b'
  expect_identical(c0$lsd_angle[c0$platform == "TaqMan"], "aabb")
  expect_identical(sum(c0$separation_angle_rad), 0.73)
  # the angle is a scatter-plot geometry, not a thermodynamic gap
  expect_match(unique(c0$outcome_definition), "clustering readout")
})

test_that("wheat2019's concordance pairs are percentages and centred on TaqMan", {
  k <- tq_read("taqman_wheat2019_concordance.csv")
  expect_identical(nrow(k), 3L)
  expect_true(all(k$concordance_pct > 0 & k$concordance_pct <= 100))
  # TaqMan is one side of two of the three pairs (the third is KASP vs rhAmp)
  expect_identical(sum(grepl("TaqMan", k$pair)), 2L)
})

test_that("the TaqMan extracts still match the fetched JATS when .tmp/lit is present", {
  root <- testthat::test_path("..", "..")
  lit <- file.path(root, ".tmp", "lit")
  skip_if_not(
    dir.exists(file.path(lit, "citrus2020")),
    "fetched literature cache is absent"
  )
  # re-run each generator in a child R process; they stop() on any
  # disagreement with the CSV they wrote
  for (g in c("wheat2019", "citrus2020", "znavsmgb2014", "apol1_2019")) {
    f <- file.path(root, "data-raw", paste0(g, ".R"))
    skip_if_not(file.exists(f), paste("missing generator", g))
    out <- suppressWarnings(system2("Rscript", c(shQuote(f)),
      stdout = TRUE, stderr = TRUE
    ))
    expect_null(attr(out, "status"))
  }
})
