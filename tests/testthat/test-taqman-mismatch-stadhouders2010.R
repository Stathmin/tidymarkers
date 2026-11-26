# =============================================================================
#  Stadhouders et al. 2010  --  what a 3' mismatch costs a 5'-nuclease assay
# =============================================================================
#  Two extracts, one DOI.  `taqman_mismatch_stadhouders2010.csv` is the measured
#  Delta-Ct / fold data; `taqman_mismatch_guidelines_stadhouders2010.csv` is the
#  paper's Table 1 guidance, verbatim.  The package has no TaqMan
#  mismatch-position criterion, so nothing here is a `criterion_spec()` row:
#  these tests keep the tables honest measurements, not numbers to render.
#
#  The two file names appear verbatim above on purpose -- `test-science-ledger.R`
#  greps this source for every extract its ledger rows name.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

measured <- function() ext("taqman_mismatch_stadhouders2010.csv")
guidance <- function() ext("taqman_mismatch_guidelines_stadhouders2010.csv")

test_that("mismatch measurements: one DOI, closed vocabularies", {
  m <- measured()
  expect_equal(nrow(m), 20L)
  expect_identical(unique(m$doi), "10.2353/jmoldx.2010.090035")
  expect_identical(unique(m$key), "stadhouders2010")
  expect_setequal(
    unique(m$setup),
    c("standard_taq", "taq_mmlv_rt", "rtth_rt", "across_master_mixes")
  )
  expect_setequal(unique(m$primer), c("forward", "reverse", "pooled"))
  expect_setequal(
    unique(m$mismatch_class),
    c(
      "all", "multiple", "purine_purine_cc", "pyrimidine_pyrimidine",
      "purine_pyrimidine", "mixed"
    )
  )
  expect_setequal(unique(stats::na.omit(m$position)), c(1, 2, 3, 5))
  expect_setequal(
    unique(m$quantity),
    c(
      "delta_ct_mean", "delta_ct_max", "delta_ct_class_terminal",
      "symmetry_gap_terminal_delta_ct", "fold_of_standard_taqman",
      "amplification_failure_percent", "fold_master_mix_difference",
      "fold_underestimation_max"
    )
  )
})

test_that("mismatch measurements: a reported range is ordered", {
  m <- measured()
  # a point estimate fills `value_low` only; a reported range fills both
  expect_true(all(is.na(m$value_high) | m$value_low <= m$value_high))
  expect_true(all(!is.na(m$value_low)))
})

test_that("mismatch measurements: one row per (setup, primer, class, position, quantity)", {
  m <- measured()
  key <- paste(m$setup, m$primer, m$mismatch_class, m$position, m$quantity,
    sep = "\r"
  )
  expect_false(anyDuplicated(key) > 0L)
})

test_that("mismatch measurements: cost decays away from the 3' end, but position 5 is not position 3", {
  m <- measured()
  pooled <- m[m$setup == "standard_taq" & m$quantity == "delta_ct_mean", ]
  pooled <- pooled[!is.na(pooled$position), ]
  at <- function(pos) pooled$value_low[pooled$position == pos][[1L]]
  # the paper's own reading: 5.19 / 1.87 / 0.53 Ct at the terminal, penultimate
  # and third-base positions
  expect_gt(at(1), at(2))
  expect_gt(at(2), at(3))
  # position 5 is measured separately and is NOT the smallest of the four:
  # asserting strict monotonicity here would be a fiction
  expect_gt(at(1), at(5))
  expect_lt(at(3), at(5))
})

test_that("mismatch measurements: the worst case is a maximum, not a mean", {
  m <- measured()
  top <- m[
    m$setup == "standard_taq" & is.na(m$position) &
      m$quantity == "delta_ct_max",
  ]
  expect_equal(nrow(top), 1L)
  mean_1 <- m$value_low[
    which(m$setup == "standard_taq" & m$position == 1 & m$quantity == "delta_ct_mean")
  ]
  expect_gt(top$value_low, mean_1)
})

test_that("guidelines: Table 1 is a complete 3 x 3 block", {
  g <- guidance()
  expect_equal(nrow(g), 9L)
  expect_identical(unique(g$doi), "10.2353/jmoldx.2010.090035")
  expect_setequal(
    unique(g$mismatch_class),
    c("purine_purine_cc", "pyrimidine_pyrimidine", "purine_pyrimidine")
  )
  expect_setequal(
    unique(g$position),
    c("terminal", "penultimate", "positions_3_5")
  )
  # every class x position pair appears exactly once
  pairs <- table(g$mismatch_class, g$position)
  expect_true(all(pairs == 1L))
})

test_that("guidelines: the advice vocabulary is closed", {
  g <- guidance()
  allowed <- c(
    "Avoid (both primers)",
    "Acceptable (both primers)",
    "Avoid in FWD, Acceptable in REV",
    "Acceptable in FWD, Avoid in REV"
  )
  cols <- c("standard_taq", "taq_mmlv_rt", "rtth_rt")
  for (col in cols) {
    expect_true(
      all(g[[col]] %in% allowed),
      label = paste("vocabulary of", col)
    )
  }
})

test_that("guidelines: standard Taq relaxes away from the 3' end", {
  g <- guidance()
  # distance rank: terminal is the 3' end, positions_3_5 is the farthest
  rank_of <- c(terminal = 1L, penultimate = 2L, positions_3_5 = 3L)
  g$rank <- rank_of[g$position]
  avoid <- grepl("^Avoid", g$standard_taq)
  for (cls in unique(g$mismatch_class)) {
    idx <- which(g$mismatch_class == cls)
    idx <- idx[order(g$rank[idx])]
    # "Avoid" may only form a prefix: once a position is acceptable, a farther
    # one cannot be avoid
    expect_false(
      any(diff(as.integer(avoid[idx])) > 0L),
      label = paste("prefix of Avoid for", cls)
    )
  }
  # the two classes the paper puts on top of the hierarchy are avoided at the
  # terminal base; a purine-pyrimidine mismatch is the class it calls acceptable
  expect_true(all(avoid[g$position == "terminal" & g$mismatch_class != "purine_pyrimidine"]))
  expect_false(any(avoid[g$position == "terminal" & g$mismatch_class == "purine_pyrimidine"]))
})
