# HRM resolvability: what the extracted evidence will and will not support.
#
# These are the extract-level checks for the jsbbs2018 / ijms2022 HRM tables;
# the design-side tests live in test-hrm.R.

test_that("jsbbs2018 Table 4 keeps the unusable pairs instead of dropping them", {
  fx <- utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", "hrm_jsbbs2018_nnn_effect.csv"),
    stringsAsFactors = FALSE)
  expect_equal(nrow(fx), 13L)
  both <- !is.na(fx$dtm_nonnns) & !is.na(fx$dtm_nnns)
  expect_equal(sum(both), 10L)
  # the three dropped pairs are the paper's own stated failures
  expect_equal(sum(!both), 3L)
  expect_equal(fx$effect[both], fx$dtm_nnns[both] - fx$dtm_nonnns[both],
               tolerance = 1e-9)
  expect_true(all(fx$dtm_nnns[both] > fx$dtm_nonnns[both]))
})

test_that("jsbbs2018 selection rule is a free-energy rule, not a dTm cut", {
  r <- utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", "hrm_jsbbs2018_selection_rule.csv"),
    stringsAsFactors = FALSE)
  thr <- r[grepl("threshold", r$criterion), ]
  expect_equal(thr$value, "1.43")
  expect_equal(thr$units, "kcal/mol")
  # and it must say so in words, because 1.43 next to a dTm column is a trap
  expect_match(thr$outcome_definition, "NOT a melting-temperature threshold", fixed = TRUE)
})

test_that("jsbbs2018 SNP-class frequencies are a composition, not an outcome", {
  f <- utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata",
                        "hrm_jsbbs2018_snp_class_frequencies.csv"),
    stringsAsFactors = FALSE)
  expect_equal(f$total_snps, c(8372, 6142, 4138822))
  # AF and TF percentages do not sum to 100 in the source itself; keep that
  sums <- f$pct_class_I + f$pct_class_II + f$pct_class_III + f$pct_class_IV
  expect_equal(round(sums, 1), c(104.5, 94.9, 100.0))
})

test_that("ijms2022 dTm table is the substitute for the blocked Liew 2004 figures", {
  d <- utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata",
                        "hrm_ijms2022_dtm_single_snp_controls.csv"),
    stringsAsFactors = FALSE)
  expect_equal(d$dtm_c[d$n_t_alleles == 0L], 0)
  expect_true(max(d$dtm_c) < 2)          # all differences are sub-2 C
  # every row must name its own basis, never a bare magic number
  expect_true(all(nzchar(d$outcome_definition)))
  expect_true(all(d$evidence_tier == "measured"))
})
