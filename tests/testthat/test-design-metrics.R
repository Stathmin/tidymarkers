# Design diagnostics: 3'-terminal stability and wobble-mismatch typing.

test_that("wobble_mismatch flags the weak G:T and A:C pairs, symmetrically", {
  expect_true(wobble_mismatch("G", "A")) # primer G faces complement(A)=T -> G:T
  expect_true(wobble_mismatch("A", "G"))
  expect_true(wobble_mismatch("C", "T")) # C:T ... (complement(T)=A -> C:A)
  expect_false(wobble_mismatch("A", "C")) # A faces G -> A:G, a strong mismatch
  expect_false(wobble_mismatch("A", "T"))
  # the flag must not depend on which allele we call REF
  expect_equal(wobble_mismatch("G", "A"), wobble_mismatch("A", "G"))
  expect_equal(wobble_mismatch("A", "C"), wobble_mismatch("C", "A"))
})

test_that("terminal_dg is finite and rewards a 3' GC clamp", {
  skip_if_no_oligoarrayaux()
  gc3 <- "ACGTACGTACGTACGCGG" # GC-rich 3' end
  at3 <- "ACGTACGTACGTACAAAA" # A/T-rich 3' end
  expect_true(is.finite(terminal_dg(gc3)))
  expect_lt(terminal_dg(gc3), terminal_dg(at3)) # more negative = more stable
  expect_true(is.na(terminal_dg("ACGT", n = 5L))) # shorter than the window
})

test_that("primer_thermo reports the terminal stability", {
  skip_if_no_oligoarrayaux()
  th <- primer_thermo("ACGTACGTACGTACGCGG")
  expect_true(all(c("tm", "gc", "hairpin_dg", "terminal_dg") %in% names(th)))
  expect_true(is.finite(th$terminal_dg))
})

test_that("the thermo maps are invariant to the parallel backend", {
  skip_if_no_oligoarrayaux()
  skip_if_not_installed("BiocParallel")
  p <- c("ACGTACGTACGTACGCGG", "TTTTGGGGCCCCAAAATTTT", "GAAGGTGACCAAGTTCATGCT")
  # a duplicated 3' end must also map back to the same (deduped) value
  s <- terminal_dg(c(p, p[3L]))
  expect_equal(terminal_dg(c(p, p[3L]), BPPARAM = BiocParallel::SerialParam()), s)
  expect_equal(primer_thermo(p, BPPARAM = BiocParallel::SerialParam()), primer_thermo(p))
  skip_on_os("windows")
  expect_equal(terminal_dg(c(p, p[3L]), BPPARAM = BiocParallel::MulticoreParam(2L)), s)
})
