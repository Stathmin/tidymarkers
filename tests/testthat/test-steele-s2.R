# =============================================================================
#  Y6  --  per-marker KASP call rates from Steele et al. 2025 (File S2)
# =============================================================================
#  The registry entry `steele2025_rice_kasp` is the only source here that
#  reaches a PER-MARKER reading of call success on a real panel.  This test is
#  the consumer gate for the ingest in data-raw/steele2025_rice_s2.R: it pins
#  the five-class partition, the measured panel totals, and three markers whose
#  counts were read off the shipped matrix by hand.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

s2 <- function() ext("steele2025_rice_kasp_calls.csv")

test_that("Y6: the S2 matrix keeps one row per marker and a full partition", {
  d <- s2()
  expect_equal(nrow(d), 3366L)
  expect_false(anyDuplicated(d$assay_id) > 0L)
  expect_false(anyNA(d$missing_calls))
  counts <- c("ref_calls", "alt_calls", "het_calls", "missing_calls", "other_calls")
  expect_equal(rowSums(as.matrix(d[counts])), rep(178L, nrow(d)))
  expect_equal(d$n_called, 178L - d$missing_calls)
  expect_true(all(d$call_rate > 0 | d$n_called == 0L))
  expect_true(all(d$call_rate <= 1))
})

test_that("Y6: the panel totals are the measured ones, not the paper's prose", {
  d <- s2()
  expect_equal(sum(d$missing_calls), 9810L)
  expect_equal(sum(d$ref_calls), 417319L)
  expect_equal(sum(d$alt_calls), 168706L)
  expect_equal(sum(d$het_calls), 2004L)
  expect_equal(sum(d$other_calls), 1309L)
  # the text claims 3,366 markers above 90 % called; the shipped matrix gives
  # 3,275, and ROADMAP Y6 records that gap rather than resolving it
  expect_equal(sum(d$passes_90), 3275L)
  expect_equal(sum(!d$passes_90), 91L)
})

test_that("Y6: three hand-read markers keep their exact counts", {
  d <- s2()
  m <- function(id) d[d$assay_id == id, ]
  expect_equal(nrow(m("RM462_SNP_nn_2")), 1L)
  expect_equal(
    unname(unlist(m("RM462_SNP_nn_2")[c("ref_calls", "alt_calls", "het_calls", "missing_calls")])),
    c(145L, 30L, 2L, 1L)
  )
  # a column printed as "?: 147" carries that allele's no-call count
  expect_true(is.na(m("PitaYL183_func_1")$ref_allele))
  expect_equal(m("PitaYL183_func_1")$alt_allele, "T:T")
  expect_equal(m("PitaYL183_func_1")$alt_calls, 31L)
  expect_equal(m("PitaYL183_func_1")$missing_calls, 147L)
  # the single marker with no call printed for either allele
  expect_true(is.na(m("RM214_SNP_nn_2")$ref_allele))
  expect_true(is.na(m("RM214_SNP_nn_2")$alt_allele))
  expect_equal(m("RM214_SNP_nn_2")$n_called, 0L)
})

test_that("Y6: the conversion registry carries the unreconciled panel totals", {
  d <- ext("conversion_rates.csv")
  s <- d[d$study == "steele2025_rice_kasp", ]
  expect_equal(nrow(s), 1L)
  expect_equal(s$n_designed, 4000L)
  expect_equal(s$n_success, 3366L)
  expect_equal(s$n_fail_other, 634L)
  expect_equal(s$n_failed, 634L)
  # the paper does not decompose the failure: no PCR stage, no no-allele stage
  expect_true(is.na(s$n_fail_pcr))
  expect_true(is.na(s$n_fail_no_allele))
  expect_equal(s$n_fail_known, 0L)
  expect_false(s$decomposes)
  expect_false(s$fully_decomposed)
  expect_equal(s$outcome_scale, "usable")
  expect_false(s$calibrated)
  expect_equal(s$doi, "10.1093/g3journal/jkae251")
  # 3,366 called + 635 failed = 4,001 against 4,000 designed: the note records
  # the arithmetic rather than inventing the missing assay
  expect_match(s$failure_note, "3,366 \\+ 635 = 4,001", fixed = FALSE)
  expect_match(s$failure_note, "15.9%", fixed = TRUE)
})

test_that("Y6: the ambiguity sweep buys designs, not assay success", {
  d <- ext("steele2025_design_sweep.csv")
  expect_equal(nrow(d), 5L)
  expect_equal(d$max_ambiguous_bases, 4:0)
  # cutting the flanking-ambiguity budget shrinks the design space and spreads
  # the candidates out
  expect_true(all(diff(d$designs) < 0))
  expect_true(all(diff(d$median_spacing_bp) > 0))
  expect_equal(d$designs[[1L]], 1392359L)
  expect_equal(d$designs[[5L]], 223783L)
  expect_equal(d$median_spacing_bp[[1L]], 68L)
  expect_equal(d$median_spacing_bp[[5L]], 580L)
  # ... while the successful-assay rate stays inside a 5-point band: this is
  # the table's point, and the guard against reading design count as quality
  expect_lt(max(d$pct_successful) - min(d$pct_successful), 5)
  expect_gte(min(d$pct_successful), 83.7)
  expect_true(all(d$pct_producing_results >= d$pct_successful))
  expect_setequal(d$doi, "10.1093/g3journal/jkae251")
})
