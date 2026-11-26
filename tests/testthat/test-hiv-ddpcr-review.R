# =============================================================================
#  Y6  --  bulk vs ddPCR detection floor, 16 Sanger-negative patients
# =============================================================================
#  `hiv_ddpcr_review_detection_floor.csv` is the corpus's measured version of
#  "the reference workflow cannot see the minority allele": the same cohort
#  called by bulk PCR and by ddPCR, with each mutation's printed VAF.  This
#  test is the consumer gate for data-raw/hiv_ddpcr_review_detection_floor.R.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

hiv <- function() ext("hiv_ddpcr_review_detection_floor.csv")

test_that("Y6: one row per patient, 16 patients", {
  d <- hiv()
  expect_equal(nrow(d), 16L)
  expect_equal(ncol(d), 8L)
  expect_false(anyDuplicated(d$case_id) > 0L)
  expect_true(all(d$rna_copies_per_ml > 0))
  # digit-group spaces in the paper's copy counts must not survive
  expect_equal(d$rna_copies_per_ml[d$case_id == 1658L], 3000000)
})

test_that("Y6: ddPCR calls more DRMs than bulk on this cohort", {
  d <- hiv()
  expect_equal(sum(d$n_bulk, na.rm = TRUE), 31L)
  expect_equal(sum(d$n_ddpcr, na.rm = TRUE), 39L)
  expect_true(sum(d$n_ddpcr, na.rm = TRUE) > sum(d$n_bulk, na.rm = TRUE))
  expect_equal(sum(d$n_ddpcr > d$n_bulk, na.rm = TRUE), 6L)
})

test_that("Y6: the detection floor is 1 % VAF", {
  d <- hiv()
  expect_equal(min(d$min_vaf, na.rm = TRUE), 0.01)
  expect_equal(max(d$max_vaf, na.rm = TRUE), 0.10)
  expect_true(all(d$min_vaf <= d$max_vaf, na.rm = TRUE))
})

test_that("Y6: discordance runs both ways, and `none` is not `NA`", {
  d <- hiv()
  # six patients get a DRM ddPCR-only, two the other way round: the paper's
  # point is that neither assay is a superset of the other
  expect_equal(
    sum(!is.na(d$ddpcr_only_drms) & d$ddpcr_only_drms != "none"), 6L
  )
  expect_equal(
    sum(!is.na(d$bulk_only_drms) & d$bulk_only_drms != "none"), 2L
  )
  expect_equal(
    d$ddpcr_only_drms[d$case_id == 1362L], "Y181C"
  )
  # `none` = no discordant variant; NA = the paper reports no calls at all.
  # readr must not collapse the two (na = c("", "NA"), never "none").
  expect_true(all(is.na(d$ddpcr_only_drms[d$case_id == 1758L])))
  expect_true(all(is.na(d$n_bulk[d$case_id == 1758L])))
  expect_equal(sum(is.na(d$n_ddpcr)), 1L)
  # 1758 is the only unreported row, so 9 + 6 + 1 and 13 + 2 + 1 partition
  expect_equal(sum(d$ddpcr_only_drms == "none", na.rm = TRUE), 9L)
  expect_equal(sum(d$bulk_only_drms == "none", na.rm = TRUE), 13L)
})

test_that("Y6: the bulk-only calls are the paper's, verbatim", {
  d <- hiv()
  # `M46IV` is printed that way in Table 2 (no such substitution exists; it is
  # almost certainly M46I/M46V).  Kept as shipped so the extract stays
  # checkable against the source, and pinned here so nobody "fixes" it.
  expect_equal(d$bulk_only_drms[d$case_id == 1931L], "M46I")
  expect_equal(d$bulk_only_drms[d$case_id == 2443L], "M46IV")
})
