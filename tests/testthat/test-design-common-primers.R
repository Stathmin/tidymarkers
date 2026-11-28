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

test_that("the lengths, gaps and anchor are whole numbers before the grid", {
  # `lengths = 20.9` recorded `len` 20.9 and `end` 123.9 for a 20 nt primer,
  # and `snp_i = 80.9` anchored the search at 100.9 while `substring()` read
  # from 100 (ADR 0055).  The 16-base window cannot yield a candidate, so a
  # value that passes these guards reaches the documented filter message.
  w <- Biostrings::DNAString("ACGTACGTACGTACGT")
  call <- function(...) design_common_primers(w, 8L, ...)
  expect_error(call(lengths = 20.9), "positive whole numbers of nucleotides")
  expect_error(call(lengths = c(20L, 20.5)), "positive whole numbers of nucleotides")
  expect_error(call(lengths = NA_real_), "positive whole numbers of nucleotides")
  expect_error(call(gaps = 20.5), "positive whole numbers of nucleotides")
  expect_error(call(gaps = c(20L, Inf)), "positive whole numbers of nucleotides")
  expect_error(design_common_primers(w, 8.5), "positive whole number")
})
