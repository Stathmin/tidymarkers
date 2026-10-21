# =============================================================================
#  test-polyploid-cotton2023
# =============================================================================
#  Cotton is the family's NEGATIVE control on the homeolog axis: an AADD
#  allotetraploid screened on cluster separation alone, in which two assays
#  designed from diploid relatives worked without any homoeolog step.  The
#  tests pin the per-assay classification and the "not stated" claim.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("polyploid/cotton2023: 11 tested, 8 usable, 3 bad, per assay", {
  c <- ext("polyploid_cotton2023_kasp_screen.csv")
  expect_equal(nrow(c), 11L)
  expect_equal(anyDuplicated(c$marker), 0L)
  expect_equal(sum(c$usable), 8L)
  expect_equal(sum(c$outcome_class == "BAD"), 3L)
  expect_equal(sum(c$outcome_class == "GOOD"), 7L)
  expect_equal(sum(c$outcome_class == "GOOD*"), 1L)
  # the split is on the assay rows, so the totals are constant columns
  expect_true(all(c$n_tested == 11L))
  expect_true(all(c$n_usable == 8L))
})

test_that("polyploid/cotton2023: the screen is cluster separation, not conversion", {
  c <- ext("polyploid_cotton2023_kasp_screen.csv")
  expect_true(all(grepl("cluster separation", c$outcome_definition)))
  expect_true(all(grepl("not a SNP-to-KASP conversion", c$outcome_definition)))
  expect_true(all(c$evidence_tier == "measured"))
  # 72.72 is 8/11 truncated; the computed fraction must agree with it
  expect_true(all(abs(c$computed_success_pct - c$printed_success_pct) < 0.05))
})

test_that("polyploid/cotton2023: no homoeolog claim is made", {
  c <- ext("polyploid_cotton2023_kasp_screen.csv")
  expect_true(all(c$ploidy == "allo-tetraploid"))
  expect_true(all(c$homeolog_handled == "not stated"))
  # the diploid-relative assays are the negative datapoint: both must be usable
  rel <- c[grepl("^i[0-9]+(Gt|Gm)$", c$marker), ]
  expect_equal(nrow(rel), 2L)
  expect_true(all(rel$usable))
  expect_true(all(grepl("diploid", rel$marker_origin)))
})

test_that("polyploid/cotton2023: the CSV still matches Table 3 in the XML", {
  xml_path <- test_path("../../.tmp/lit/pmc10005503/fulltext.xml")
  skip_if_not(file.exists(xml_path))
  skip_if_not_installed("xml2")
  doc <- xml2::xml_ns_strip(xml2::read_xml(xml_path))
  tabs <- xml2::xml_find_all(doc, ".//table-wrap")
  cap <- vapply(tabs, function(t) {
    paste(xml2::xml_text(
      xml2::xml_find_all(t, ".//caption")
    ), collapse = " ")
  }, character(1))
  hit <- grep("KASP assay screening results", cap, fixed = TRUE)
  expect_length(hit, 1L)
  rows <- xml2::xml_find_all(tabs[[hit]], ".//tr")
  cells <- lapply(rows, function(r) {
    gsub("[[:space:]]+", " ", xml2::xml_text(
      xml2::xml_find_all(r, ".//td|.//th")
    ))
  })
  ids <- Filter(function(r) length(r) == 6L && grepl("^i[0-9]", r[[1L]]), cells)
  expect_equal(length(ids), 11L)
  verdicts <- vapply(ids, `[[`, character(1), 6L)
  expect_equal(sum(grepl("^GOOD$", verdicts)), 7L)
  expect_equal(sum(grepl("^BAD$", verdicts)), 3L)
})
