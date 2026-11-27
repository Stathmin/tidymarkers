# Threshold provenance and the cross-table Tm check.

test_that("marker_priors names every threshold and its source", {
  p <- marker_priors()
  expect_true(all(c("parameter", "value", "source") %in% names(p)))
  expect_true("as_tm (C)" %in% p$parameter)
  expect_true("Tm model" %in% p$parameter)
  expect_false(any(is.na(p$source) | !nzchar(p$source)))
})

test_that("thermo_crosscheck reports the Tm spread across NN tables", {
  seqs <- c("GAAGGTGACCAAGTTCATGCT", "TTTTGGCGCGCCGATATATAT")
  x <- thermo_crosscheck(seqs)
  expect_equal(nrow(x), 2L)
  # the first table is our default, so it matches primer_tm()
  expect_equal(x[["DNA_NN_SantaLucia_2004"]], primer_tm(seqs))
  expect_true(all(x$spread >= 0 & is.finite(x$spread)))
  # the models really do disagree (otherwise the check is vacuous)
  expect_true(any(x$spread > 0))
})

test_that("primer_tm honours an alternative NN table", {
  s <- "GAAGGTGACCAAGTTCATGCTACGTAC"
  expect_false(isTRUE(all.equal(primer_tm(s), primer_tm(s, nn_table = "DNA_NN_Breslauer_1986"))))
})

test_that("outcomes become per-scheme counts, never a probability", {
  o <- tibble::tibble(
    id = c("a1", "a1", "a2"),
    target = c("Chr1:10 A>C", "Chr1:20 A>C", "Chr1:10 A>C"),
    scheme = c("snp", "hrm", "snp"),
    result = c("works", "fails", "works"),
    notes = c(NA, "no amplification", NA),
    calibrated = FALSE
  )
  p <- marker_priors(o)
  expect_equal(p$scheme, c("hrm", "snp"))
  expect_equal(unname(p$works), c(0L, 2L))
  expect_equal(unname(p$fails), c(1L, 0L))
  expect_true(all(p$calibrated == FALSE))
  expect_false(any(grepl("prob|p_success|expected", names(p))))
  expect_true(all(grepl("a1|a2", p$source)))

  empty <- marker_priors(o[0, ])
  expect_equal(nrow(empty), 0L)
  expect_error(marker_priors(tibble::tibble(x = 1)), "marker_outcomes")
})

test_that("a prior reorders attempts without changing the verdict", {
  genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
  loci <- utils::read.csv(
    system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
    colClasses = "character"
  )[, c("chrom", "pos", "ref", "alt")]
  loci$pos <- as.integer(loci$pos)
  priors <- tibble::tibble(
    scheme = c("snp", "hrm"), works = c(1L, 0L), fails = c(0L, 1L),
    unknown = c(0L, 0L), source = "test", calibrated = FALSE
  )
  base <- marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = c("snp", "hrm"))
  nudge <- marker_plan(genome) + marker_targets(loci) +
    marker_tech(schemes = c("snp", "hrm"), priors = priors)
  expect_equal(tidymarkers:::plan_schemes(base), c("snp", "hrm"))
  expect_equal(tidymarkers:::plan_schemes(nudge), c("snp", "hrm"))

  reverse <- tibble::tibble(
    scheme = c("snp", "hrm"), works = c(0L, 1L), fails = c(1L, 0L),
    unknown = c(0L, 0L), source = "test", calibrated = FALSE
  )
  flipped <- marker_plan(genome) + marker_targets(loci) +
    marker_tech(schemes = c("snp", "hrm"), priors = reverse)
  expect_equal(tidymarkers:::plan_schemes(flipped), c("hrm", "snp"))
  expect_error(
    marker_tech(priors = tibble::tibble(x = 1)),
    "marker_priors"
  )
  bad <- function(scheme, works, fails) {
    tibble::tibble(scheme = scheme, works = works, fails = fails)
  }
  expect_error(marker_tech(priors = bad("hmr", 1L, 0L)), "Unknown scheme")
  expect_error(marker_tech(priors = bad("snp", -1L, 0L)), "non-negative")
  expect_error(marker_tech(priors = bad("snp", NA_integer_, 0L)), "finite")
  expect_error(marker_tech(priors = bad("snp", "1", 0L)), "non-negative")
})
