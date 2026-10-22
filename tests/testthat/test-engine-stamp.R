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

assay_with_stamp <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      g <- ath_genome()
      snp <- pick_target_snp(
        ensembl_variants("4", 1188000L, 1192000L), g, ensembl_gene("AT4G00650")
      )
      cached <<- assemble_kasp_assay(g, snp, flank = 150L, run_blast = FALSE)
    }
    cached
  }
})

test_that("an evaluation carries the engines that produced it", {
  # ZW2.2: a verdict that cannot name its engines cannot be re-derived after a
  # model bumps.  The stamp rides on the checks table, not on the session.
  stamp <- marker_evaluation_stamp(assay_with_stamp()$evaluation)
  expect_type(stamp, "character")
  expect_true(all(REQUIRED_ENGINES %in% names(stamp)))
  expect_equal(stamp, marker_engine_stamp())
})

test_that("the stamp survives the rows appended after it", {
  assay <- assay_with_stamp()
  ev <- assay$evaluation
  # a real population: marker_population_gate() has three rows to append
  pop <- tibble::tibble(
    alt_freq = 0.2, n_haplotypes = 10L, n_alt_haplotypes = 2L, min_alt_freq = 0.01
  )
  with_pop <- marker_with_population_checks(ev, list(population = pop, snp = assay$snp))
  expect_gt(nrow(with_pop), nrow(ev))
  expect_equal(marker_evaluation_stamp(with_pop), marker_evaluation_stamp(ev))
  # and through a verb that drops attributes
  out <- marker_with_stamp(dplyr::bind_rows(ev), from = ev)
  expect_equal(marker_evaluation_stamp(out), marker_engine_stamp())
})

test_that("a table built by hand has no stamp rather than an invented one", {
  ev <- marker_evaluation_table(
    stage = "assay", criterion = "x", observed = "1", limit = "1",
    pass = TRUE, level = "assay"
  )
  expect_null(marker_evaluation_stamp(ev))
})
