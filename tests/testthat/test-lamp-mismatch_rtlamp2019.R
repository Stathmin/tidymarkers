# LAMP mismatch position ranking (DENV RT-LAMP, 12 single-position mutants).
# Data-only, from inst/extdata; the generator re-parses Table 1 from .tmp/lit.

mm_read <- function(name) {
  utils::read.csv(testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE
  )
}

test_that("both mismatch files exist with the required columns", {
  for (f in c(
    "lamp_mismatch_rtlamp2019_tt.csv",
    "lamp_mismatch_rtlamp2019_position_rules.csv",
    "lamp_mismatch_rtlamp2019_contradictions.csv"
  )) {
    d <- mm_read(f)
    expect_gt(nrow(d), 0L)
    expect_true(all(c("outcome_definition", "evidence_tier", "source_file") %in% names(d)))
    expect_true(all(d$evidence_tier == "moderate"))
  }
})

test_that("13 templates x 2 conditions = 26 Tt rows, 3 replicates each", {
  d <- mm_read("lamp_mismatch_rtlamp2019_tt.csv")
  expect_identical(nrow(d), 26L)
  expect_identical(sum(d$template == "Wild"), 2L)
  expect_identical(length(unique(d$template)), 13L)
  expect_setequal(
    unique(d$condition),
    c("novel (Bst 2.0 + 0.15 U Q5)", "conventional (Bst 2.0)")
  )
  # every replicate is a real time, and the mean lies inside them
  expect_true(all(d$t1 > 0 & d$t2 > 0 & d$t3 > 0))
  lo <- pmin(d$t1, d$t2, d$t3)
  hi <- pmax(d$t1, d$t2, d$t3)
  expect_true(all(d$mean >= lo & d$mean <= hi))
  expect_true(all(d$sd >= 0))
})

test_that("a 3'-terminal mismatch always SLOWS the conventional reaction", {
  d <- mm_read("lamp_mismatch_rtlamp2019_tt.csv")
  conv <- d[d$condition == "conventional (Bst 2.0)", ]
  wt <- conv$mean[conv$template == "Wild"]
  expect_identical(length(wt), 1L)
  mut <- conv$template != "Wild"
  expect_true(all(conv$mean[mut] > wt))
  # the delay is a delay in minutes, never a call
  expect_match(unique(d$outcome_definition), "mismatch-tolerance measure")
  expect_match(unique(d$outcome_definition), "NOT two-allele discrimination")
})

test_that("the position ranking is FIP-F2 > BLP > F1c > F3", {
  r <- mm_read("lamp_mismatch_rtlamp2019_position_rules.csv")
  get1 <- function(p) r$conventional_tt_over_wt_min[startsWith(r$item, p)]
  expect_lt(get1("mismatch at the 3' end of F3"), get1("mismatch at the 5' end of FIP"))
  expect_lt(get1("mismatch at the 3' end of FIP"), get1("mismatch at the 3' end of BLP"))
  f2 <- get1("mismatch at the 3' end of FIP")
  f3 <- get1("mismatch at the 3' end of F3")
  expect_gt(f2 / f3, 3) # the F2 position delays >= 3x more than F3
  # the novel (proofreading) method shrinks every gap, but never to zero
  expect_true(all(r$novel_tt_over_wt_max_min < r$conventional_tt_over_wt_max_min))
})

test_that("the recorded arithmetic deviations are not silently corrected", {
  k <- mm_read("lamp_mismatch_rtlamp2019_contradictions.csv")
  expect_gte(nrow(k), 3L)
  expect_true(any(grepl("conv mean", k$text_version)))
  d <- mm_read("lamp_mismatch_rtlamp2019_tt.csv")
  # the printed cells survive verbatim: F1-Mu-A conventional mean is 28.6, not 28.7
  expect_identical(
    d$mean[d$template == "F1-Mu-A" & d$condition == "conventional (Bst 2.0)"],
    28.6
  )
})
