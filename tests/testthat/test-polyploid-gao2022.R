# =============================================================================
#  tests/testthat/test-polyploid-gao2022.R
# =============================================================================
#  A nested outcome chain carried as three rows, so that a later reader cannot
#  collapse it into one rate.  The tests pin:
#    * strict monotonicity 77 > 49 > 24;
#    * the file-backed tier (24) equals the workbook's assay rows;
#    * the two sentence-only tiers are labelled weak, the table tier measured;
#    * the workbook's allele-prime end is recorded, not assumed.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name), show_col_types = FALSE,
                  na = c("", "NA"))
}

test_that("polyploid/gao2022: 77 -> 49 -> 24 is strict and in order", {
  g <- ext("polyploid_gao2022_sugarcane_chain.csv")
  expect_equal(nrow(g), 3L)
  expect_equal(g$tier, c("selected", "genotyped in both parents",
                         "in line with the map"))
  expect_equal(g$n, c(77, 49, 24))
  expect_true(all(diff(g$n) < 0))
  expect_true(all(c("ploidy", "homeolog_handled", "outcome_definition",
                    "evidence_tier", "source_file") %in% names(g)))
  expect_true(all(g$ploidy == "auto-polyploid"))
  expect_true(all(g$homeolog_handled %in% c("not stated", "not applicable", "no")))
})

test_that("polyploid/gao2022: only the workbook tier is measured", {
  g <- ext("polyploid_gao2022_sugarcane_chain.csv")
  expect_equal(g$evidence_tier, c("weak", "weak", "measured"))
  # the measured tier must carry a file locator, the weak tiers a text one
  expect_true(grepl("sheet S10", g$source_locator[[3L]]))
  expect_true(all(grepl("fulltext", g$source_locator[1:2])))
  # every tier carries its definition; nothing is an empty string
  expect_true(all(nzchar(g$outcome_definition)))
  expect_true(all(grepl("two parents", g$outcome_definition)))
})

test_that("polyploid/gao2022: the allele-specific end is recorded, not assumed", {
  g <- ext("polyploid_gao2022_sugarcane_chain.csv")
  # the workbook prints the SNP at the 5' end of the allele primers (only 4/48
  # at the 3' end), so the assay row must say so rather than inherit the wheat
  # "SNP at the AS 3' base" rule
  expect_equal(g$as_allele_match_3prime[[3L]], 4)
  expect_equal(g$as_allele_match_5prime[[3L]], 13)
  expect_true(all(is.na(g$as_allele_match_3prime[1:2])))
})

test_that("polyploid/gao2022: the CSV still matches Supplementary Table 10", {
  lit <- test_path("../../.tmp/lit/gao2022")
  skip_if_not(file.exists(file.path(lit, "supp", "Table_1.XLSX")))
  source(test_path("../../data-raw/lib/xlsx-read.R"))
  s10 <- xlsx_sheets(file.path(lit, "supp", "Table_1.XLSX"))[["S10"]]
  ids <- vapply(s10[-(1:2)], function(r) {
    v <- xcell(r, "A")
    if (is.na(v)) "" else v
  }, character(1))
  expect_equal(sum(nzchar(ids)), 24L)
  g <- ext("polyploid_gao2022_sugarcane_chain.csv")
  expect_equal(g$n[[3L]], 24)
})
