# The fixed (given) primer constructor must answer "unknown, not pass" when the
# PCR-efficiency engine is absent, instead of aborting before the evaluation:
# discrimination used to be computed unconditionally, so a missing `hybrid-min`
# killed the exported entry point before its own unknown path could run.

test_that("fixed primers without the PCR engine report discrimination unknown, not an error", {
  f <- snp_fixture()
  len <- 20L
  pos <- f$snp$pos
  as_ref <- substring(f$W, pos - len + 1L, pos)
  as_alt <- paste0(substring(f$W, pos - len + 1L, pos - 1L), f$alt)
  com <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAString(substring(f$W, pos + 30L, pos + 30L + len - 1L))
  ))
  # the seam, not PATH: engine_has() is what the constructor gates on
  local_mocked_bindings(engine_has = function(...) FALSE, .package = "tidymarkers")
  out <- suppressWarnings(assemble_primer_assay(
    f$genome, f$snp,
    c(AS_ref = as_ref, AS_alt = as_alt, COM = com)
  ))
  expect_identical(out$efficiency$discrimination, NA_real_)
  row <- out$evaluation[out$evaluation$criterion == "discrimination", ]
  expect_equal(nrow(row), 1L)
  expect_equal(as.character(row$status), "unknown")
  expect_false(row$pass)
  expect_true(is.na(row$observed))
  expect_match(row$detail, "OligoArrayAux is unavailable", fixed = TRUE)
  # an unknown row is never a pass, so the assay cannot come out usable; the
  # verdict-ordering rule itself is covered in test-kasp-checks.R
  expect_false(identical(evaluation_verdict(out), "usable"))
})
