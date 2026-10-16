# The experimenter-facing risk roll-up over the acceptance checks.

mk_assay_ev <- function(failed = character()) {
  crit <- c("discrimination", "AS GC", "longest homopolymer (bp)", "off-target exact primer hits (n)")
  list(evaluation = tibble::tibble(criterion = crit, pass = !crit %in% failed))
}

test_that("marker_risk tiers grave failures above minor ones", {
  expect_equal(as.character(marker_risk(mk_assay_ev())$risk), "low")
  expect_equal(as.character(marker_risk(mk_assay_ev("AS GC"))$risk), "medium")
  expect_equal(as.character(marker_risk(mk_assay_ev("longest homopolymer (bp)"))$risk), "medium")
  expect_equal(as.character(marker_risk(mk_assay_ev("discrimination"))$risk), "high")
  expect_equal(as.character(marker_risk(mk_assay_ev("off-target exact primer hits (n)"))$risk), "high")
})

test_that("marker_risk returns the failed criteria and handles a missing evaluation", {
  r <- marker_risk(mk_assay_ev(c("AS GC", "discrimination")))
  expect_equal(r$n_failed, 2L)
  expect_setequal(r$failed[[1L]], c("AS GC", "discrimination"))
  expect_true(is.na(marker_risk(list())$risk))
})

test_that("marker_risk grades discrimination from its margin, not a bare fail", {
  ev <- function(disc) {
    list(evaluation = tibble::tibble(
      criterion = c("discrimination", "AS GC"),
      observed = c(disc, 0.5),
      pass = c(disc >= KASP_LIMITS$discrimination, TRUE)
    ))
  }
  # >= min_disc: low, even though we do not use the pass column here
  expect_equal(as.character(marker_risk(ev(0.7))$disc_risk), "low")
  expect_equal(as.character(marker_risk(ev(0.7))$risk), "low")
  # just below min_disc (within the margin): medium, not high
  expect_equal(as.character(marker_risk(ev(0.45))$disc_risk), "medium")
  expect_equal(as.character(marker_risk(ev(0.45))$risk), "medium")
  # well below: high
  expect_equal(as.character(marker_risk(ev(0.2))$disc_risk), "high")
  expect_equal(as.character(marker_risk(ev(0.2))$risk), "high")
})
