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

test_that("a KASP assay's reported dG is its scheme-temperature evaluation (Z3)", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  a <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp, flank = 150L))
  # the assembler threads its own scheme into the reported screen
  expect_equal(a$scheme, "snp")
  set <- c(AS_ref = a$tails[["AS_ref"]], AS_alt = a$tails[["AS_alt"]], COM = a$tails[["COM"]])
  expect_equal(a$complementarity, dimer_thermo(set, conditions = scheme_conditions("snp")))
  # the verdict is unchanged: a 60 C dG that is >= 37 C is still >= 37 C
  cold <- dimer_thermo(set, conditions = reaction_conditions())
  expect_true(all(a$complementarity$dg >= cold$dg - 1e-8))
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

test_that("dimer_thermo is monotone in the reaction temperature", {
  skip_if_no_oligoarrayaux()
  p <- c(
    A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT",
    C = "GGGGCCCCAAAATTTTACGT"
  )
  temps <- c(25, 37, 50, 60, 65)
  sim <- seq_len(nrow(dimer_thermo(p)))
  at <- function(t) {
    dimer_thermo(p, conditions = reaction_conditions(temperature = t))$dg[sim]
  }
  by_temp <- vapply(temps, at, numeric(length(sim)))
  # raising T destabilises a duplex, so dG must never fall further below its
  # 25 C value, and the mean must rise over the whole window
  expect_true(all(diff(colMeans(by_temp)) > 0))
  expect_gte(min(by_temp), min(by_temp[, 1L]) - 1e-8)
  # strictly rising wherever there is something to melt: the fully-paired pairs
  # (both >= -5 kcal/mol at 25 C) and, up to 50 C, every pair.  A pair that does
  # not bind approaches 0 from below and the curve flattens or turns over at the
  # top of the window (raw hybrid-min gives the same numbers, so this is the
  # model, not the wrapper).
  paired <- by_temp[, 1L] <= -5
  expect_true(any(paired))
  expect_true(all(diff(t(by_temp[paired, ])) > 0))
  expect_true(all(diff(t(by_temp[, temps %in% c(25, 37, 50)])) > 0))
})

test_that("dimer_thermo at 37 C equals the legacy call", {
  skip_if_no_oligoarrayaux()
  p <- c(A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT")
  # the pre-seam code: deltas at 37 C, 50 mM Na, from the full matrix
  full <- tidyr::crossing(a = names(p), b = names(p))
  legacy <- vapply(seq_len(nrow(full)), function(i) {
    oligo_dg(p[[full$a[i]]], p[[full$b[i]]], na = 0.05, temp = 37)
  }, numeric(1))
  expect_equal(dimer_thermo(p, conditions = list(temperature = 37))$dg, legacy)
  expect_equal(
    dimer_thermo(p, conditions = reaction_conditions(temperature = 37)),
    dimer_thermo(p)
  )
})

test_that("dimers at LAMP temperature differ from 37 C", {
  skip_if_no_oligoarrayaux()
  p <- c(
    A = "ACGTACGTACGTACGTACGT", B = "TTTTGGGGCCCCAAAATTTT",
    C = "GGGGCCCCAAAATTTTACGT"
  )
  lamp <- scheme_conditions("lamp")
  expect_equal(attr(lamp, "temperature_range"), c(60, 65))
  d37 <- dimer_thermo(p)
  d_lamp <- dimer_thermo(p, conditions = scheme_conditions("lamp"))
  expect_equal(d37$a, d_lamp$a)
  expect_equal(d37$b, d_lamp$b)
  expect_false(isTRUE(all.equal(d37$dg, d_lamp$dg)))
  # weaker at the reaction temperature: never more negative than at 37 C, and
  # strictly so for at least one pair
  expect_true(all(d_lamp$dg >= d37$dg - 1e-8))
  expect_true(any(d_lamp$dg > d37$dg + 1e-8))
})

test_that("hairpin_dg_thermo defaults to 37 C and follows the conditions", {
  skip_if_no_oligoarrayaux()
  s <- c("GCGCGCGCGCGCGCGCGCGC", "AAAAAAAAAAAAAAAAAAAA", "GAAGGTGACCAAGTTCATGCT")
  expect_equal(hairpin_dg_thermo(s), oligo_dg_vec(s, na = 0.05, temp = 37))
  expect_equal(hairpin_dg_thermo(s, conditions = NULL), hairpin_dg_thermo(s))
  expect_false(isTRUE(all.equal(
    hairpin_dg_thermo(s), hairpin_dg_thermo(s, conditions = scheme_conditions("lamp"))
  )))
})

test_that("duplex_dg_thermo matches dimer_thermo and the scheme temperature", {
  skip_if_no_oligoarrayaux()
  a <- "GAAGGTGACCAAGTTCATGCT"
  b <- "AGCATGAACTTGGTCACCTTC"
  expect_equal(duplex_dg_thermo(a, b), oligo_dg_vec(a, b, na = 0.05, temp = 37))
  expect_equal(
    duplex_dg_thermo(a, b, conditions = reaction_conditions(temperature = 37)),
    duplex_dg_thermo(a, b)
  )
  # a 20-mer duplex is held together by more than 10 kcal/mol near 37 C, so it is
  # certainly melted by the LAMP window
  expect_lt(duplex_dg_thermo(a, b), -10)
  expect_gt(duplex_dg_thermo(a, b, conditions = scheme_conditions("lamp")), -10)
  # and the vector form agrees with the matrix engine
  d <- dimer_thermo(c(a = a, b = b), conditions = scheme_conditions("pcr"))
  expect_equal(
    duplex_dg_thermo(a, b, conditions = scheme_conditions("pcr")),
    d$dg[d$a == "a" & d$b == "b"]
  )
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
