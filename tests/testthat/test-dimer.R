# Duplex free energies must form a symmetric matrix with a homodimer diagonal.

test_that("dimer_thermo is symmetric and includes homodimers", {
  skip_if_no_oligoarrayaux()
  p <- c(
    A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT",
    C = "GGGGCCCCAAAATTTTACGT"
  )
  d <- dimer_thermo(p)
  expect_equal(nrow(d), length(p)^2)
  m <- matrix(d$dg, nrow = length(p), dimnames = list(names(p), names(p)))
  expect_equal(m, t(m), tolerance = 1e-8)
  expect_false(any(is.na(diag(m))))
})

test_that("dimer_thermo is invariant to the parallel backend", {
  skip_if_no_oligoarrayaux()
  skip_if_not_installed("BiocParallel")
  p <- c(
    A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT",
    C = "GGGGCCCCAAAATTTTACGT"
  )
  s <- dimer_thermo(p)
  expect_equal(dimer_thermo(p, BPPARAM = BiocParallel::SerialParam()), s)
  skip_on_os("windows")
  # SOCK workers (a real second process) also exercise worker_biostrings();
  # bpstop() reaps them so the session exits cleanly
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  expect_equal(dimer_thermo(p, BPPARAM = bp), s)
})
