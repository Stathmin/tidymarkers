# =============================================================================
#  The engine stamp is complete and serialisable
# =============================================================================
#  A manifest that names three of four engines cannot support the claim "the
#  results changed because a model changed": the absent engine is
#  indistinguishable from one that was never installed.  So the stamp must
#  always carry all four measured engines, each with a real version string, and
#  must survive being flattened into a single CSV or JSON column.
#
#  Scope is deliberately narrow: this is provenance, not a Tm test.  These
#  versions are what THIS machine resolves; they pin nothing about accuracy.

REQUIRED_ENGINES <- c("DECIPHER", "OligoArrayAux", "TmCalculator", "BLAST+")

test_that("stamp names every engine the measurement layer depends on", {
  stamp <- marker_engine_stamp()
  expect_type(stamp, "character")
  expect_length(stamp, length(marker_engine_versions()$engine))
  expect_false(anyDuplicated(names(stamp)) > 0L)
  expect_true(all(REQUIRED_ENGINES %in% names(stamp)))
})

test_that("every required engine reports a non-empty version string", {
  stamp <- marker_engine_stamp()
  for (engine in REQUIRED_ENGINES) {
    v <- stamp[[engine]]
    expect_false(is.na(v), info = paste(engine, "is not installed on this machine"))
    expect_true(is.character(v) && nzchar(v),
      info = paste(engine, "reported an empty version")
    )
  }
})

test_that("a missing engine is written as NA, not dropped", {
  # the regression this guards: filtering on `available` made a four-engine
  # stamp silently become a three-engine one
  fake <- c(
    DECIPHER = "3.8.1", OligoArrayAux = NA_character_,
    TmCalculator = "1.1.2", `BLAST+` = "2.16.0+"
  )
  csv <- marker_engine_stamp_field("csv", stamp = fake)
  expect_type(csv, "character")
  expect_length(csv, 1L)
  expect_true(grepl("OligoArrayAux=NA", csv, fixed = TRUE))
  expect_true(grepl("DECIPHER=3.8.1", csv, fixed = TRUE))
})

test_that("the stamp flattens to CSV and JSON for a manifest column", {
  csv <- marker_engine_stamp_field("csv")
  expect_type(csv, "character")
  expect_length(csv, 1L)
  # one cell, no newlines, no commas to split a CSV row on
  expect_false(grepl("[,\n]", csv))
  for (engine in REQUIRED_ENGINES) {
    expect_true(grepl(paste0(engine, "="), csv, fixed = TRUE),
      info = paste(engine, "missing from the CSV field")
    )
  }

  js <- marker_engine_stamp_field("json")
  expect_type(js, "character")
  expect_length(js, 1L)
  parsed <- jsonlite::fromJSON(js)
  expect_true(all(REQUIRED_ENGINES %in% names(parsed)))
  expect_equal(unname(parsed[["BLAST+"]]), unname(marker_engine_stamp()[["BLAST+"]]))
})
