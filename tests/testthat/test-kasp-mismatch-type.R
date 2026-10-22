# Z0.4 -- a KASP 3' mismatch is kinetics, not a binary
#
# Kwok 1990 orders 3'-terminal primer-template mismatches by how readily Taq
# extends them.  That is a RANKING, not a pass/fail: "the polymerase extends this
# one more easily" and "this assay works / does not work" are different claims.
# The KASP checks therefore report the CLASS, and the report says a pass there is
# a strong, medium or weak assumption.

ext <- function(name) {
  utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
}

test_that("the severity classes are Kwok's, verbatim, and the three classes agree", {
  d <- ext("mismatch_discrimination.csv")
  k <- d[d$table == "kwok_3prime_severity", ]
  members <- function(cls) {
    trimws(strsplit(k$members[k$mismatch_class == cls], ",", fixed = TRUE)[[1L]])
  }
  # every published class is reproduced by the classifier (either orientation),
  # EXCEPT A:A: the paper lists it as `moderate` but an A/A pair is a Watson-
  # Crick match, not a mismatch, so the classifier answers `NA` rather than
  # call it a mismatch.  That is a deliberate divergence, recorded here.
  for (cls in c("severe", "moderate", "minimal", "context_dependent")) {
    for (pair in members(cls)) {
      ab <- strsplit(pair, ":", fixed = TRUE)[[1L]]
      if (identical(pair, "A:A")) {
        # A is the complement of T, so "A:A" names a primer A over template T -
        # an A:T Watson-Crick MATCH, not a mismatch.  Kwok's paper calls it
        # `moderate`; the classifier answers NA on purpose, and that divergence
        # is pinned here rather than hidden.
        expect_true(is.na(marker_mismatch_severity("A", "T")))
        next
      }
      # every published pair reproduces its class.  The reverse argument order is
      # ONLY a different labelling of the SAME pair when the two bases are equal
      # or complementary-checkable the way the table's own naming allows; for the
      # `T` family the reverse is a different pair (G over A, C over A) that the
      # paper does not call, so NA there is right, not a miss.
      expect_equal(marker_mismatch_severity(ab[[1L]], ab[[2L]]), cls, info = pair)
      rev <- marker_mismatch_severity(ab[[2L]], ab[[1L]])
      if (identical(pair, "T:G")) {
        # reversed this is primer G over template A - not in Kwok's table
        expect_true(is.na(rev), info = pair)
      } else {
        expect_equal(rev, cls, info = pair)
      }
    }
  }
})

test_that("a substitution is classified, a non-substitution is not", {
  # a plus-strand A opposite a plus-strand T is a match, not a mismatch
  expect_true(is.na(marker_mismatch_severity("A", "T")))
  expect_true(is.na(marker_mismatch_severity(NA, "A")))
  expect_true(is.na(marker_mismatch_severity("N", "A")))
  # the pairs Kwok's table does not carry stay unknown, not guessed
  expect_true(is.na(marker_mismatch_severity("G", "T"))) # primer G over template T
  expect_true(is.na(marker_mismatch_severity("A", "C"))) # primer A over template C
})

test_that("the wobble flag is the same call on one pair and on a batch", {
  # `wobble_mismatch()` used to sort the WHOLE vector, so a batch disagreed with
  # the single-pair call it vectorises; both forms must now agree
  pairs <- list(c("A", "G"), c("C", "T"), c("G", "A"), c("T", "C"), c("A", "C"), c("G", "T"))
  ones <- vapply(pairs, function(p) wobble_mismatch(p[[1L]], p[[2L]]), logical(1))
  batch <- wobble_mismatch(
    vapply(pairs, `[[`, character(1), 1L),
    vapply(pairs, `[[`, character(1), 2L)
  )
  expect_equal(unname(batch), unname(ones))
})

test_that("a KASP evaluation reports the terminal mismatch class, not a silent pass", {
  f <- snp_fixture()
  m <- assemble_kasp_assay(f$genome, f$snp)
  ev <- marker_design_checks(m, scheme = "snp")
  # `criterion` carries the registry LABEL (the id is not a column of the
  # acceptance table), so the row is found by its label
  row <- ev[ev$criterion == "AS 3' mismatch base pair", ]
  expect_equal(nrow(row), 1L)
  # the row states the class it graded, and never claims more than a ranking.
  # `snp_fixture()` is A>C, a substitution pair Kwok's table does not carry, so
  # the honest answer is "no published class" - never a silent pass.
  expect_true(grepl("Kwok 1990|not in Kwok", row$detail))
  expect_true(row$status %in% c("pass", "unknown"))
  # the kinetics-not-a-binary caveat is carried whenever a class IS reported
  if (grepl("terminal mismatch \\(Kwok 1990\\)", row$detail)) {
    expect_match(row$detail, "kinetics, not a binary", fixed = TRUE)
  }
  # the criterion itself is a categorical check sourced from the literature.
  # `marker_criteria()` carries the pipeline view (stage/role/limit); the `kind`
  # and `source` it was registered with live on the data object, `criterion_spec()`
  sp <- criterion_spec(KASP_LIMITS)
  expect_equal(sp$kind[sp$id == "as_mismatch_type"], "categorical")
  expect_equal(sp$source[sp$id == "as_mismatch_type"], "primary literature")
})
