# The experimenter-facing risk roll-up over the acceptance checks.

mk_assay_ev <- function(failed = character()) {
  crit <- c("discrimination", "AS GC", "longest homopolymer (bp)", "off-target exact primer hits (n)")
  list(evaluation = tibble::tibble(criterion = crit, pass = !crit %in% failed))
}

test_that("kasp_risk tiers grave failures above minor ones", {
  expect_equal(as.character(kasp_risk(mk_assay_ev())$risk), "low")
  expect_equal(as.character(kasp_risk(mk_assay_ev("AS GC"))$risk), "medium")
  expect_equal(as.character(kasp_risk(mk_assay_ev("longest homopolymer (bp)"))$risk), "medium")
  expect_equal(as.character(kasp_risk(mk_assay_ev("discrimination"))$risk), "high")
  expect_equal(as.character(kasp_risk(mk_assay_ev("off-target exact primer hits (n)"))$risk), "high")
})

test_that("kasp_risk returns the failed criteria and handles a missing evaluation", {
  r <- kasp_risk(mk_assay_ev(c("AS GC", "discrimination")))
  expect_equal(r$n_failed, 2L)
  expect_setequal(r$failed[[1L]], c("AS GC", "discrimination"))
  expect_true(is.na(kasp_risk(list())$risk))
})
