# =============================================================================
#  SNaPshot: detection vs RNA input, and the paper's prose claims
# =============================================================================
#  Source: PMC13644189 (DOI 10.1186/s12864-026-13173-0), Table S4.
#  The CSV is a 14 locus x 5 input grid recomputed from the workbook.  These
#  tests fail if a count is hand-edited out of range, if the grid loses a row,
#  or if the prose claims stop matching the table they were compared against.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}
proposed <- function(name) {
  # data-raw/ is build-ignored: it exists in a checkout, never in the installed
  # package R CMD check runs against
  f <- test_path("../../data-raw/proposed", name)
  skip_if_not(file.exists(f), "data-raw/proposed not present (installed check)")
  readr::read_csv(f,
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("SNaPshot sensitivity: the grid is complete and physically possible", {
  d <- ext("snapshot2026_sensitivity.csv")
  expect_equal(nrow(d), 70L) # 14 loci x 5 inputs
  expect_equal(length(unique(d$locus)), 14L)
  expect_setequal(unique(d$rna_input_label), c("10ng", "5ng", "1ng", "0.5ng", "0.1ng"))
  expect_false(anyDuplicated(d[c("locus", "rna_input_ng")]) > 0L)

  # a count cannot exceed the number of replicates, nor be negative
  expect_true(all(d$n_replicates == 3L))
  expect_true(all(d$n_detected >= 0L & d$n_detected <= 3L))
  expect_true(all(d$n_failed >= 0L & d$n_failed <= 3L))
  expect_equal(d$n_detected + d$n_failed, d$n_replicates)
  # rate is derived from the counts, not typed in
  expect_equal(d$failure_rate_pct, round(100 * d$n_failed / d$n_replicates, 1L))
  # every row carries the rule that defines "detected"
  expect_true(all(grepl("100", d$detection_rule, fixed = TRUE)))
  expect_true(all(d$evidence_tier == "strong"))
})

test_that("SNaPshot sensitivity: saturating input never fails", {
  d <- ext("snapshot2026_sensitivity.csv")
  ten <- d[d$rna_input_label == "10ng", ]
  expect_equal(nrow(ten), 14L)
  expect_true(all(ten$n_failed == 0L))
  # at 1 ng the paper's weaker statement holds: no locus is lost in >1 of 3
  one <- d[d$rna_input_label == "1ng", ]
  expect_true(all(one$n_detected >= 2L))
})

test_that("SNaPshot sensitivity: the prose-vs-table comparison is honest", {
  c <- proposed("snapshot2026_failure_claims.rows.csv")
  expect_gt(nrow(c), 10L)
  expect_true(all(c("claimed_failure_rate_pct", "tabled_failure_rate_pct", "reproduces") %in% names(c)))
  # the comparison records MISMATCHES: an all-true file would mean the
  # disagreements were quietly dropped
  expect_true(any(!c$reproduces))
  # every claim is credited to the paper or contradicted by it, never "NA"
  expect_true(all(nzchar(c$verdict)))
  # the tabulated values must be one of the four possible 3-replicate rates
  expect_true(all(c$tabled_failure_rate_pct %in% c(0, 33.3, 66.7, 100)))
  # a 0% tabled rate against a 100% claim is recorded as contradicted
  contra <- c$claimed_failure_rate_pct == 100 & c$tabled_failure_rate_pct == 0
  expect_true(all(grepl("contradicted", c$verdict[contra])))
})

test_that("SNaPshot specificity: the 500-sample matrix is self-consistent", {
  s <- ext("snapshot2026_specificity.csv")
  expect_equal(nrow(s), 95L) # 19 loci x 5 fluids
  expect_true(all(s$n_samples == 100L))
  expect_true(all(s$n_detected >= 0L & s$n_detected <= 100L))
  expect_equal(s$detection_rate_pct, round(100 * s$n_detected / s$n_samples, 1L))
  # the paper's headline: all 14 cSNPs in every PB sample
  pb <- s[s$body_fluid == "Peripheral" & s$marker_class == "cSNP", ]
  expect_equal(nrow(pb), 14L)
  expect_true(all(pb$n_detected == 100L))
  # specificity: cSNPs stay below 50% outside the two blood fluids
  other <- s[s$marker_class == "cSNP" & !s$body_fluid %in% c("Peripheral", "Menstrual"), ]
  expect_true(all(other$detection_rate_pct < 50))
  # GYPA is the shared blood marker, not a one-fluid control.  It is 84/100 in
  # PB in this table, so the test asserts the LANDED value and that it is the
  # marker's best fluid, rather than asserting a textbook 100.
  gy <- s[s$locus == "GYPA", ]
  expect_equal(nrow(gy), 5L)
  expect_equal(gy$detection_rate_pct[gy$body_fluid == "Peripheral"], 84)
  expect_equal(gy$detection_rate_pct[gy$body_fluid == "Menstrual"], 67)
  expect_lt(gy$detection_rate_pct[gy$body_fluid == "Saliva"], 10)
  # and the paper's own 76% claim is recorded as a mismatch, not overwritten
  claim <- proposed("snapshot2026_specificity_claims.rows.csv")
  expect_equal(nrow(claim), 1L)
  expect_equal(claim$claimed_detection_rate_pct, 76)
  expect_false(claim$reproduces)
  expect_equal(claim$tabled_detection_rate_pct, 67)
})

test_that("SNaPshot: the CSVs still agree with the fetched workbook", {
  wb <- test_path("../../.tmp/lit/PMC13644189/supplements/12864_2026_13173_MOESM2_ESM.xlsx")
  skip_if_not(file.exists(wb), "fetched supplement not present")
  # the workbook must still carry the two sheets the CSVs were built from
  con <- unz(wb, "xl/workbook.xml")
  on.exit(close(con), add = TRUE)
  wbx <- paste(readLines(con, warn = FALSE), collapse = "")
  expect_true(grepl("Table S4 Sensitivity", wbx, fixed = TRUE))
  expect_true(grepl("Table S8 500 sample analysis", wbx, fixed = TRUE))
})
