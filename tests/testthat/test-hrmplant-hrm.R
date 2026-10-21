# HRM in crops + the Liew physics anchor. Invariants only: these assertions pin
# the numbers that data-raw/hrmplant.R derives and stop()s on. No network.

csv <- function(name) {
  p <- test_path("..", "..", "inst", "extdata", name)
  if (!file.exists(p)) skip(paste("missing", name))
  utils::read.csv(p, stringsAsFactors = FALSE)
}

test_that("jsbbs2018 Table 3 grid is 6 bimorphisms x 16 and only 32 cells are bp-neutral", {
  g <- csv("hrm_jsbbs2018_context_grid.csv")
  tb <- table(g$bimorphism)
  expect_equal(length(tb), 6L)
  expect_true(all(tb == 16L))
  expect_equal(sum(g$bimorphism %in% c("C/T", "G/A")), 32L)
  expect_true(all(g$evidence_tier == "model"))
})

test_that("jsbbs2018 Table 1 frame is 70 markers, 38 + 32, no M##_TF block", {
  f <- csv("hrm_jsbbs2018_table1_frame.csv")
  expect_equal(nrow(f), 70L)
  expect_equal(sum(f$experiment == "Exp. 1"), 38L)
  expect_equal(sum(f$experiment == "Exp. 2"), 32L)
  expect_false(any(grepl("^M", f$marker_id)))
})

test_that("jsbbs2018 dTm classes: C/T highest mean, C/G lowest", {
  s <- csv("hrm_jsbbs2018_dtm_by_bimorphism.csv")
  expect_equal(sum(s$n), 56L)
  expect_equal(s$bimorphism[which.max(s$dtm_mean)], "C/T")
  expect_equal(s$bimorphism[which.min(s$dtm_mean)], "C/G")
})

test_that("gundry2008 Table 4: 6 loci, n sums to 1147, BX647987 is the only regression", {
  a <- csv("hrm_gundry2008_homozygote_accuracy.csv")
  expect_equal(nrow(a), 6L)
  expect_equal(sum(a$n), 1147L)
  expect_equal(sum(a$accuracy_after_pct >= 99), 5L)
  worse <- a$gene[a$accuracy_after_pct < a$accuracy_before_pct]
  expect_equal(worse, "BX647987")
})

test_that("gundry2008 Table 3: 8 amplicons x 2 homozygotes, symmetric pairs predict 0", {
  t3 <- csv("hrm_gundry2008_amplicon_tm.csv")
  expect_equal(nrow(t3), 16L)
  expect_true(all(t3$dtm_after >= 0 & t3$dtm_after <= 1))
  expect_setequal(unique(t3$amplicon[t3$symmetric_nn]), c("CPS1 51 bp", "OTC 68 bp"))
  expect_true(all(t3$dtm_predicted[t3$symmetric_nn] == 0))
  # SIX non-symmetric predictions are printed; the prose quotes only four.
  pp <- sort(unique(round(abs(t3$dtm_predicted[!t3$symmetric_nn]), 2)))
  expect_equal(pp, c(0.04, 0.11, 0.18, 0.19, 0.24, 0.31))
})

test_that("liew2004 bands are secondary_cited, never measured", {
  l <- csv("hrm_liew2004_restated_bands.csv")
  expect_equal(nrow(l), 6L)
  expect_true(all(l$evidence_tier == "secondary_cited"))
  expect_false(any(l$evidence_tier == "measured"))
  v <- stats::setNames(l$value, l$quantity)
  expect_lt(v[["dtm_at_to_gc_lower_c"]], v[["dtm_at_to_gc_upper_c"]])
  expect_lt(v[["dtm_other_substitutions_upper_c"]], v[["dtm_at_to_gc_lower_c"]])
  expect_equal(unname(v[["homozygote_resolved_pct_gundry"]] +
                      v[["nn_symmetric_bp_neutral_pct_gundry"]]), 100)
  expect_true(all(grepl("NOT MEASURED HERE", l$outcome_definition, fixed = TRUE)))
})
