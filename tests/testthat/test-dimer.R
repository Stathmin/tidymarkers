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

test_that("dimer_thermo defaults reproduce the 37 C / 50 mM legacy call", {
  skip_if_no_oligoarrayaux()
  p <- c(
    A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT",
    C = "GGGGCCCCAAAATTTTACGT"
  )
  d <- dimer_thermo(p)
  expect_equal(d, dimer_thermo(p, conditions = reaction_conditions()))
  legacy <- vapply(seq_len(nrow(d)), function(i) {
    oligo_dg(p[[d$a[i]]], p[[d$b[i]]], na = 0.05, temp = 37)
  }, numeric(1))
  expect_equal(d$dg, legacy)
})

test_that("dimer_thermo honours a custom reaction temperature and salt", {
  skip_if_no_oligoarrayaux()
  p <- c(
    A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT",
    C = "GGGGCCCCAAAATTTTACGT"
  )
  d37 <- dimer_thermo(p)
  d60 <- dimer_thermo(p, conditions = reaction_conditions(temperature = 60))
  expect_equal(d37$a, d60$a)
  expect_equal(d37$b, d60$b)
  expect_false(isTRUE(all.equal(d37$dg, d60$dg)))
  d_salt <- dimer_thermo(p, conditions = reaction_conditions(Na = 1000))
  expect_false(isTRUE(all.equal(d37$dg, d_salt$dg)))
})

test_that("dimer_thermo accepts a partial conditions list", {
  skip_if_no_oligoarrayaux()
  p <- c(A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT")
  expect_equal(
    dimer_thermo(p, conditions = list(temperature = 37)),
    dimer_thermo(p)
  )
  expect_error(dimer_thermo(p, conditions = list(tmp = 37)), "Unknown reaction condition")
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
  skip_if_not_full()
  # SOCK workers (a real second process) also exercise worker_biostrings();
  # bpstop() reaps them so the session exits cleanly
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  expect_equal(dimer_thermo(p, BPPARAM = bp), s)
})
