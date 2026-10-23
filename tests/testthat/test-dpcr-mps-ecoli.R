# Multi-assay ddPCR validation panel: 10.3390/mps7060094 (PMC11587158).
#
# These tests assert the properties that make this file a PANEL rather than one
# more single-assay LoB: several assays on ONE dilution series, with an
# explicit working range that is narrower than the linear range.

ext <- function(name) {
  utils::read.csv(
    tm_extdata(name),
    stringsAsFactors = FALSE
  )
}

test_that("mps2024: seven assays share one five-level dilution series", {
  d <- ext("dpcr_mps_ecoli_dilution.csv")
  # 7 assays x 5 levels x 3 days = 105 cells, of which 3 are blank in the PDF
  expect_equal(nrow(d), 102L)
  expect_equal(length(unique(d$assay)), 7L)
  expect_setequal(
    unique(d$assay),
    c("uidA", "lacY", "eaeA", "rfbE", "stx1", "stx2", "Z3276")
  )
  expect_setequal(unique(d$level), 1:5)
  expect_true(all(table(d$assay, d$level) %in% 1:3))
  # the two partial cells are the ones the PDF shows blank
  expect_equal(unname(table(d$assay, d$level)["stx2", "1"]), 1L)
  expect_equal(unname(table(d$assay, d$level)["rfbE", "5"]), 2L)
})

test_that("mps2024: copy numbers fall monotonically across the dilution series", {
  d <- ext("dpcr_mps_ecoli_dilution.csv")
  # every replicate is a positive concentration, and lower levels are lower
  vals <- c(d$rep1, d$rep2, d$rep3)
  expect_false(anyNA(vals)) # blank cells were dropped, not filled
  expect_true(all(vals > 0))
  mean_by_level <- tapply(
    rowMeans(d[c("rep1", "rep2", "rep3")], na.rm = TRUE),
    d$level, median
  )
  expect_true(all(diff(as.numeric(mean_by_level)) < 0))
  # a five-level 1:10 series spans at least 1e3 in copies/uL
  expect_gt(mean_by_level[[1]] / mean_by_level[[5]], 1e3)
})

test_that("mps2024: unreadable cells are dropped, never imputed", {
  d <- ext("dpcr_mps_ecoli_dilution.csv")
  # every shipped row carries three real replicates
  expect_true(all(d$n_replicates_parsed == 3L))
  expect_false(anyNA(d[c("rep1", "rep2", "rep3")]))
  # exactly three of the 105 possible cells are absent
  expect_equal(105L - nrow(d), 3L)
  expect_true(all(d$dilution_print > 0))
})

test_that("mps2024: working range is narrower than the linear range", {
  l <- ext("dpcr_mps_ecoli_limits.csv")
  expect_equal(nrow(l), 7L)
  expect_true(all(l$linear_range_low_copies_ul < 6.6))
  expect_true(all(l$lod_copies_ul < 6.6))
  expect_true(all(l$linear_r2 > 0.998))
  # the author-stated exception is rfbE, and it must stay the unique minimum
  expect_equal(l$assay[which.min(l$lod_copies_ul)], "rfbE")
  expect_equal(min(l$lod_copies_ul), 0.17)
  expect_equal(sum(l$lod_copies_ul == 0.27), 6L)
  expect_true(all(l$evidence_tier == "measured"))
  expect_true(all(nzchar(l$outcome_definition)))
})
