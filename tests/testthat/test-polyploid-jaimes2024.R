# =============================================================================
#  test-polyploid-jaimes2024
# =============================================================================
#  The sugarcane row is the family's AUTOPOLYPLOID anchor: no subgenome exists,
#  so the outcome axis is allele DOSAGE, not homoeolog specificity.  What these
#  tests protect:
#    * "converted" and "dosage assigned" stay different outcomes (42 of 76);
#    * the denominator is the 76 SNPs the paper tested, not the 42 that worked;
#    * the rows say "no homoeolog axis" in the package's vocabulary;
#    * the band split that did NOT reproduce from the supplement stays landed
#      as "not landed" rather than quietly turning into a number.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("polyploid/jaimes2024: 42 of 76 with the denominator kept", {
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  expect_true(all(c(
    "ploidy", "homeolog_handled", "outcome_definition",
    "evidence_tier", "source_file"
  ) %in% names(j)))
  conv <- j[j$key == "jaimes2024", ]
  expect_equal(sort(conv$n_tested), rep(76, nrow(conv)))
  expect_equal(sum(conv$n_worked[conv$metric == "dosage assigned (n)"]), 42L)
  # worked + failed reconstructs the tested count on every conversion row ...
  expect_true(all(conv$n_worked + conv$n_failed == conv$n_tested))
  # ... and the two rows that state the same 42 (fraction / count) agree
  expect_equal(unique(conv$n_worked[conv$metric %in% c(
    "dosage assigned (fraction)", "dosage assigned (n)"
  )]), 42)
  # ... and the fraction row is the same ratio
  fr <- conv$value[conv$metric == "dosage assigned (fraction)"]
  expect_length(fr, 1L)
  expect_lt(abs(fr - 42 / 76), 1e-3)
})

test_that("polyploid/jaimes2024: dosage outcome is not an amplification outcome", {
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  expect_true(all(grepl("dosage", j$outcome_definition)))
  expect_true(all(grepl("NOT an amplification", j$outcome_definition)))
  # the paper's own words for the 34 are the definition, so they must be quoted
  expect_true(all(grepl("did not cluster", j$outcome_definition)))
  # autopolyploid, and explicitly no homoeolog axis - not "not stated"
  expect_true(all(j$ploidy == "auto-polyploid"))
  expect_true(all(j$homeolog_handled == "not applicable"))
})

test_that("polyploid/jaimes2024: the unreproducible band split is not a number", {
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  b <- j[j$key == "jaimes2024_kasp_bands", ]
  expect_equal(nrow(b), 1L)
  expect_true(grepl("not landed", b$band_detail, fixed = TRUE))
  # and the mean correlation it did land is inside the printed +/- sd envelope
  expect_equal(b$value, 0.3815, tolerance = 1e-4)
  base <- j[j$key == "jaimes2024_baseline", ]
  expect_equal(base$value, 0.616, tolerance = 1e-3)
})

test_that("polyploid/jaimes2024: the CSV still matches Appendix S3", {
  lit <- test_path("../../.tmp/lit/jaimes2024")
  skip_if_not(file.exists(file.path(lit, "supp", "APS3-12-e11574-s002.xlsx")))
  source(test_path("../../data-raw/lib/xlsx-read.R"))
  s3 <- xlsx_sheets(file.path(lit, "supp", "APS3-12-e11574-s002.xlsx"))[[1L]]
  rows <- s3[7:82]
  kasp <- vapply(rows, function(r) xcell(r, "F"), character(1))
  n_assigned <- sum(!is.na(kasp))
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  expect_equal(n_assigned, 42L)
  expect_equal(j$n_worked[j$metric == "dosage assigned (n)"], as.numeric(n_assigned))
})
