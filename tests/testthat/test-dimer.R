# Duplex free energies must form a symmetric matrix with a homodimer diagonal.

test_that("dimer_thermo is symmetric and includes homodimers", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  p <- c(A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT",
         C = "GGGGCCCCAAAATTTTACGT")
  d <- dimer_thermo(p)
  expect_equal(nrow(d), length(p)^2)
  m <- matrix(d$dg, nrow = length(p), dimnames = list(names(p), names(p)))
  expect_equal(m, t(m), tolerance = 1e-8)
  expect_false(any(is.na(diag(m))))
})
