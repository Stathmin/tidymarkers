# =============================================================================
#  test-polyploid-durum2019
# =============================================================================
#  The tetraploid-vs-hexaploid comparison.  The tests pin the arithmetic the
#  paper's footnote implies:
#    * the two sub-tables sum to the combined total on all three columns;
#    * the failed % divide by ALL assays, the polymorphic % by NON-FAILED ones;
#    * the two conventions are genuinely different here, so a swap is visible;
#    * the D-genome contrast (AABBDD vs AABB panel origin) is labelled.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("polyploid/durum2019: sub-tables sum to the combined row", {
  d <- ext("polyploid_durum2019_hexaploid_vs_tetraploid.csv")
  expect_equal(nrow(d), 3L)
  expect_equal(d$n_assays[[3L]], sum(d$n_assays[1:2]))
  expect_equal(d$n_failed[[3L]], sum(d$n_failed[1:2]))
  expect_equal(d$n_polymorphic[[3L]], sum(d$n_polymorphic[1:2]))
  expect_equal(d$n_assays, c(73647, 7940, 81587))
  # no assay is both failed and polymorphic
  expect_true(all(d$n_failed + d$n_polymorphic <= d$n_assays))
})

test_that("polyploid/durum2019: the two percentage conventions are separate", {
  d <- ext("polyploid_durum2019_hexaploid_vs_tetraploid.csv")
  failed_pct <- 100 * d$n_failed / d$n_assays
  poly_pct <- 100 * d$n_polymorphic / (d$n_assays - d$n_failed)
  expect_equal(d$failed_pct_of_assays, failed_pct, tolerance = 1e-3)
  expect_equal(d$polymorphic_pct_of_nonfailed, poly_pct, tolerance = 1e-3)
  # the printed cells are kept verbatim next to the computed ones, and the only
  # cell where they differ is the durum 8.8 (printed) vs 8.854 (computed)
  expect_equal(d$polymorphic_pct_of_nonfailed_printed, c(30.2, 51.7, 32.3))
  gap <- abs(d$failed_pct_of_assays_printed - d$failed_pct_of_assays)
  expect_lt(max(gap), 0.06)
  # and they must NOT be interchangeable: the printed 32.3 is a non-failed rate
  expect_gt(
    abs(
      d$polymorphic_pct_of_nonfailed[[3L]] -
        d$polymorphic_pct_of_assays[[3L]]
    ),
    1
  )
  expect_true(all(grepl("NON-FAILED", d$outcome_definition)))
  expect_true(all(grepl("not KASP", d$outcome_definition)))
})

test_that("polyploid/durum2019: the D-genome axis is labelled, not inferred", {
  d <- ext("polyploid_durum2019_hexaploid_vs_tetraploid.csv")
  expect_equal(d$ploidy, c("allo-hexaploid", "allo-tetraploid", "mixed"))
  expect_true(grepl("AABBDD", d$species[[1L]]))
  expect_true(grepl("AABB", d$species[[2L]]))
  # the hexaploid assay is designed without homoeolog separation, the durum
  # catalogue is built from homoeolog-resolved SNPs - both must be explicit
  expect_equal(d$homeolog_handled, c("no", "yes", "mixed"))
  expect_true(all(d$evidence_tier == "measured"))
})

test_that("polyploid/durum2019: the CSV still matches Table 2 in the XML", {
  xml_path <- test_path("../../.tmp/lit/pmc6456968/fulltext.xml")
  skip_if_not(file.exists(xml_path))
  skip_if_not_installed("xml2")
  doc <- xml2::xml_ns_strip(xml2::read_xml(xml_path))
  tabs <- xml2::xml_find_all(doc, ".//table-wrap")
  cap <- vapply(tabs, function(t) {
    paste(xml2::xml_text(xml2::xml_find_all(t, ".//caption")), collapse = " ")
  }, character(1))
  hit <- grep("Number of failed, polymorphic and total SNPs", cap, fixed = TRUE)
  expect_length(hit, 1L)
  rows <- xml2::xml_find_all(tabs[[hit]], ".//tr")
  cells <- lapply(rows, function(r) {
    xml2::xml_text(xml2::xml_find_all(r, ".//td"))
  })
  total_row <- Filter(function(r) length(r) == 4L && r[[1L]] == "Total", cells)
  expect_length(total_row, 1L)
  expect_equal(gsub("[^0-9]", "", total_row[[1L]][[2L]]), "81587")
  expect_equal(gsub("[^0-9]", "", sub(" .*$", "", total_row[[1L]][[3L]])), "5485")
})
