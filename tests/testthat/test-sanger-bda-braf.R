# Measured Sanger VAF floor: 10.1038/s41598-021-88391-5 (PMC8079675).
#
# The point of this file is that the floor is MEASURED and that it belongs to a
# named enrichment-plus-caller pipeline.  A test that merely checked "Sanger is
# about 20%" would be satisfied by the textbook rows this replaces.

ext <- function(name) {
  utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE
  )
}

test_that("cheng2021: the 23 clinical calls partition exactly as the paper states", {
  d <- ext("sanger_bda_braf_calls.csv")
  expect_equal(nrow(d), 23L)
  expect_equal(sum(d$bda_identity == "V600E"), 8L)
  expect_equal(sum(d$bda_identity %in% c("V600K", "V600R")), 6L)
  expect_equal(sum(d$bda_identity == "V600-WT"), 8L)
  expect_equal(sum(d$is_ffpe_damage), 1L)
  # prose: "Nine samples were identified as BRAF V600-WT"
  expect_equal(sum(d$bda_identity == "V600-WT" | d$is_ffpe_damage), 9L)
  expect_true(all(d$ffpe_repaired))
})

test_that("cheng2021: '<0.1' stays a bound and is never a zero", {
  d <- ext("sanger_bda_braf_calls.csv")
  expect_true(all(d$is_bound == (d$bda_identity == "V600-WT")))
  expect_true(all(is.na(d$bda_vaf_pct[d$is_bound])))
  expect_true(all(d$is_bound[grepl("<", d$raw)]))
  # no bound row may masquerade as a measured low value
  expect_equal(sum(d$bda_vaf_pct < 0.1, na.rm = TRUE), 0L)
  expect_true(all(d$bda_vaf_pct >= 0 & d$bda_vaf_pct <= 100, na.rm = TRUE))
})

test_that("cheng2021: the lowest positive call is 0.20% at sample 129533", {
  d <- ext("sanger_bda_braf_calls.csv")
  v <- d[!is.na(d$bda_vaf_pct), ]
  low <- which.min(v$bda_vaf_pct)
  expect_equal(v$bda_vaf_pct[low], 0.20)
  # sample_id is a label; read.csv infers integer, integer vs character is not
  # the property under test
  expect_equal(as.character(v$sample_id[low]), "129533")
  expect_equal(d$patient[as.character(d$sample_id) == "129533"], 21L)
})

test_that("cheng2021: the floor rows are measured and pipeline-specific", {
  f <- ext("sanger_bda_braf_floor.csv")
  expect_equal(nrow(f), 3L)
  expect_true(all(f$evidence_tier == "measured"))
  expect_true(all(f$usable_as_calibration))
  expect_true(all(grepl("enrich", f$method, ignore.case = TRUE)))
  expect_setequal(f$value_pct, c(0.2, 0.1, 0.5))
  # the extract must not be mistaken for a Sanger-wide constant
  expect_true(all(grepl("NOT the population-Sanger floor", f$outcome_definition)))
})

test_that("sanger floor tables keep the textbook rows separate from the measured ones", {
  tf <- ext("sanger_detection_floor.csv")
  expect_true(all(tf$evidence_tier == "textbook"))
  expect_false(any(tf$usable_as_calibration))
  # every textbook row is still bounded and cited
  expect_true(all(tf$value_pct >= 5 & tf$value_pct <= 30))
  expect_true(all(nzchar(tf$source_file)))
})
