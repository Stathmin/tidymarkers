# HRM resolvability: what the extracted evidence will and will not support.
#
# These are the extract-level checks for the jsbbs2018 / ijms2022 HRM tables;
# the design-side tests live in test-hrm.R.

test_that("jsbbs2018 Table 4 keeps the unusable pairs instead of dropping them", {
  fx <- utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", "hrm_jsbbs2018_nnn_effect.csv"),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(fx), 13L)
  both <- !is.na(fx$dtm_nonnns) & !is.na(fx$dtm_nnns)
  expect_equal(sum(both), 10L)
  # the three dropped pairs are the paper's own stated failures
  expect_equal(sum(!both), 3L)
  expect_equal(fx$effect[both], fx$dtm_nnns[both] - fx$dtm_nonnns[both],
    tolerance = 1e-9
  )
  expect_true(all(fx$dtm_nnns[both] > fx$dtm_nonnns[both]))
})

test_that("jsbbs2018 selection rule is a free-energy rule, not a dTm cut", {
  r <- utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", "hrm_jsbbs2018_selection_rule.csv"),
    stringsAsFactors = FALSE
  )
  thr <- r[grepl("threshold", r$criterion), ]
  expect_equal(thr$value, "1.43")
  expect_equal(thr$units, "kcal/mol")
  # and it must say so in words, because 1.43 next to a dTm column is a trap
  expect_match(thr$outcome_definition, "NOT a melting-temperature threshold", fixed = TRUE)
})

test_that("jsbbs2018 SNP-class frequencies are a composition, not an outcome", {
  f <- utils::read.csv(
    testthat::test_path(
      "..", "..", "inst", "extdata",
      "hrm_jsbbs2018_snp_class_frequencies.csv"
    ),
    stringsAsFactors = FALSE
  )
  expect_equal(f$total_snps, c(8372, 6142, 4138822))
  # AF and TF percentages do not sum to 100 in the source itself; keep that
  sums <- f$pct_class_I + f$pct_class_II + f$pct_class_III + f$pct_class_IV
  expect_equal(round(sums, 1), c(104.5, 94.9, 100.0))
})

test_that("HRM resolvability rests on the paper's per-bimorphism dTm, and says when it has none", {
  # the decision's source of truth: the lookup must BE the extract, not a copy
  # of it that can drift
  d <- utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", "hrm_jsbbs2018_dtm_by_bimorphism.csv"),
    stringsAsFactors = FALSE
  )
  expect_setequal(names(HRM_BIMORPHISM_DTM), d$bimorphism)
  for (i in seq_len(nrow(d))) {
    got <- HRM_BIMORPHISM_DTM[[d$bimorphism[[i]]]]
    expect_equal(unname(got[["min"]]), d$dtm_min[[i]], tolerance = 1e-3)
    expect_equal(unname(got[["mean"]]), d$dtm_mean[[i]], tolerance = 1e-3)
  }
  # the threshold is the source's own stated rule, not a round number of ours
  expect_equal(HRM_DTM_RESOLVABLE, 0.5)

  # C/T is the best class, and the threshold is a PASS bar: a class whose worst
  # marker still reaches 0.5 is resolvable, even if the class mean is lower
  expect_equal(hrm_resolvability("C", "T")$status, "resolvable")
  expect_equal(hrm_resolvability("T", "C")$status, "resolvable") # orientation-free
  expect_equal(hrm_resolvability("A", "C")$status, "resolvable")
  # A/G has one sub-threshold marker (0.245) but a class mean of 0.65: the class
  # is resolvable as a class, and the single weak marker is visible in dtm_min
  expect_equal(hrm_resolvability("G", "A")$status, "resolvable")
  expect_lt(hrm_resolvability("G", "A")$dtm_min, HRM_DTM_RESOLVABLE)
  # C/G's ONLY measured marker is 0.198 C: below the default, no mean to save it
  expect_equal(hrm_resolvability("C", "G")$status, "below_threshold")
  expect_equal(hrm_resolvability("C", "G")$dtm_min, 0.1976, tolerance = 1e-3)
  # G/T: the class mean (0.519) just clears 0.5 while its worst member does not,
  # so it is the one class that is resolvable but not assured
  expect_lt(hrm_resolvability("G", "T")$dtm_min, HRM_DTM_RESOLVABLE)
  # no landed dTm for other pairs => unknown, never a pass and never a fail
  for (b in list(c("A", "T"), c("A", "A"), c("A", NA), c("R", "C"))) {
    expect_equal(hrm_resolvability(b[[1L]], b[[2L]])$status, "unknown")
    expect_true(is.na(hrm_resolvability(b[[1L]], b[[2L]])$dtm_min))
    expect_equal(hrm_resolvability(b[[1L]], b[[2L]])$evidence, "none landed")
  }
  # the free-energy and length inputs are recorded, not thresholded here (the
  # paper's 1.43 kcal/mol rule was one experiment's, and is a filter only)
  r <- hrm_resolvability("C", "T", dg3bp = 1.2, product = 74)
  expect_equal(r$dg3bp, 1.2)
  expect_equal(r$product, 74)
  expect_false("resolvable" %in% names(r)) # no per-SNP yes/no column

  # a multi-locus panel is not answerable per SNP: 14 of 36 controls collided
  p <- utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", "hrm_ijms2022_unresolvable_3snp.csv"),
    stringsAsFactors = FALSE
  )
  expect_equal(unique(p$n_controls_total), 36L)
  expect_equal(unique(stats::na.omit(p$n_uniquely_identified)), 22L)
  expect_equal(sum(!is.na(p$control_a)), 7L) # the indistinguishable pairs
})

test_that("ijms2022 dTm table is the substitute for the blocked Liew 2004 figures", {
  d <- utils::read.csv(
    testthat::test_path(
      "..", "..", "inst", "extdata",
      "hrm_ijms2022_dtm_single_snp_controls.csv"
    ),
    stringsAsFactors = FALSE
  )
  expect_equal(d$dtm_c[d$n_t_alleles == 0L], 0)
  expect_true(max(d$dtm_c) < 2) # all differences are sub-2 C
  # every row must name its own basis, never a bare magic number
  expect_true(all(nzchar(d$outcome_definition)))
  expect_true(all(d$evidence_tier == "measured"))
})
