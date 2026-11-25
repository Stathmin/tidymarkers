# design_common_primers(): the geometry and its abort path.

test_that("a template with no room for a common primer aborts with the filter message", {
  # 16-base window, lengths start at 18, so no candidate row survives the
  # geometry.  The empty candidate set used to reach primer_tm(character(0)),
  # which aborted inside TmCalculator before this message could be raised.
  expect_error(
    design_common_primers(Biostrings::DNAString("ACGTACGTACGTACGT"), 8L),
    "No common primer passed the filters",
    fixed = TRUE
  )
})
