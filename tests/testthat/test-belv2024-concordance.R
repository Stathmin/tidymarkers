# =============================================================================
#  Y6  --  eleven-assay pairwise concordance from Belv 2024 (Table 2)
# =============================================================================
#  The corpus holds one platform's readout in `mismatch_family_*.csv` and the
#  dPCR Tm tables.  `belv2024_assay_concordance.csv` is the only extract that
#  measures how far the ASSAY CHOICE itself moves the call: eleven assays on
#  the same samples, 55 pairwise Cohen's kappa and percent-agreement cells.
#  This test is the consumer gate for data-raw/belv2024_assay_concordance.R.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

belv <- function() ext("belv2024_assay_concordance.csv")

test_that("Y6: the pairwise matrix is complete and one row per pair", {
  d <- belv()
  expect_equal(nrow(d), 55L)
  expect_equal(ncol(d), 4L)
  expect_false(anyDuplicated(paste(d$assay_a, d$assay_b)) > 0L)
  expect_false(anyNA(d$kappa))
  expect_false(anyNA(d$agreement_pct))
  expect_true(all(d$agreement_pct >= 0 & d$agreement_pct <= 100))
})

test_that("Y6: the assay set is the paper's, ddPCR6 included", {
  d <- belv()
  expect_setequal(unique(c(d$assay_a, d$assay_b)), c(
    "qPCR1", "qPCR2", "qPCR3", "qPCR4", "qPCR5", "ddPCR6",
    "qPCR7", "qPCR8", "qPCR9", "qPCR10", "qPCR11"
  ))
  # the sixth slot is a digital assay; a positional read that assumed qPCR6
  # would rename it and shift every later column
  expect_true("ddPCR6" %in% d$assay_b)
})

test_that("Y6: the concordance distribution is the measured one", {
  d <- belv()
  kappa <- as.numeric(d$kappa)
  # 30 of 55 pairs fall below kappa 0.4 -- the extract's whole point: most
  # assay pairs barely agree beyond chance
  expect_equal(sum(kappa < 0.4), 30L)
  expect_equal(sum(kappa >= 0.8), 3L)
  expect_equal(sort(kappa[kappa >= 0.8]), c(0.845, 0.927, 1))
  expect_equal(min(kappa), 0.070)
  expect_equal(max(kappa), 1)
  expect_equal(min(d$agreement_pct), 66.8)
  expect_equal(max(d$agreement_pct), 100)
  # three decimals in the file, so 1.000 and 0.070 do not round away
  raw <- readLines(tm_extdata("belv2024_assay_concordance.csv"))
  kappa_text <- vapply(strsplit(raw[-1L], ",", fixed = TRUE), `[[`, character(1), 3L)
  expect_true(all(grepl("^[01][.][0-9]{3}$", kappa_text)))
})

test_that("Y6: the paper's two traps stay resolved", {
  d <- belv()
  pair <- function(a, b) d[d$assay_a == a & d$assay_b == b, ]
  # the paper prints this kappa with a decimal comma (`0,075`); a naive
  # as.numeric() would have turned it into NA and dropped the row
  expect_equal(pair("qPCR3", "qPCR5")$kappa, 0.075)
  expect_equal(pair("qPCR3", "qPCR5")$agreement_pct, 79.5)
  # this pair sits in the same row and is NOT the comma cell
  expect_equal(pair("qPCR3", "ddPCR6")$kappa, 0.539)
  expect_equal(pair("qPCR3", "ddPCR6")$agreement_pct, 84.1)
})

test_that("Y6: the two triangles stay distinct", {
  d <- belv()
  # kappa and agreement are read from opposite triangles of the same matrix;
  # if the coordinates were swapped the two columns would be equal
  expect_false(isTRUE(all.equal(
    as.numeric(d$kappa), d$agreement_pct / 100
  )))
  # one perfect pair: identical calls, kappa 1.000 and 100 % agreement
  perfect <- d[d$kappa == 1, ]
  expect_equal(nrow(perfect), 1L)
  expect_equal(perfect$agreement_pct, 100)
})
