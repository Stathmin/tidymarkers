# The rprimer cross-check: an independent designer/reporter on an alignment.

test_that("alignment_from builds an ungapped alignment and rejects ragged input", {
  a <- alignment_from(c(REF = "ACGTACGT", ALT = "ACGAACGT"))
  expect_s4_class(a, "DNAMultipleAlignment")
  expect_error(alignment_from(c(a = "ACGT", b = "ACG")), "equal length")
})

test_that("rprimer_assays designs and reports a match table", {
  skip_if_not_installed("rprimer")
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 120L, TRUE), collapse = ""))
  aln <- alignment_from(c(
    REF = w,
    ALT = paste0(substr(w, 1L, 60L), "A", substr(w, 62L, 120L))
  ))
  des <- rprimer_assays(aln, probe = FALSE, lengthPrimer = c(18, 22), tmPrimer = c(45, 70))
  expect_true(all(c("type", "fwd", "start", "end", "sequence") %in% names(des)))
  if (nrow(des)) {
    expect_true(all(
      c("perfectMatch", "oneMismatch", "offTargetMatch") %in% names(des)
    ))
    expect_true(all(dplyr::between(des$offTargetMatch, 0, 1)))
  }
})
