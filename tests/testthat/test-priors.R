# Threshold provenance and the cross-table Tm check.

test_that("kasp_priors names every threshold and its source", {
  p <- kasp_priors()
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
