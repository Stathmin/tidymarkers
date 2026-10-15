# Readout models: what a scheme simulates, and the explicit "not simulated" flag.

test_that("kasp_readout returns a modelled readout for a hydrolysis-probe scheme", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_taqman_marker(genome, snp, flank = 150L)
  r <- kasp_readout(m)
  expect_true(r$simulated)
  expect_equal(r$scheme, "taqman")
  expect_match(r$model, "hydrolysis-probe")
  expect_setequal(r$table$allele, c("REF", "ALT"))
  # each allele-specific probe reports on/off efficiencies and their difference
  expect_equal(r$table$d_eff, r$table$on - r$table$off)
  expect_true(all(r$table$on >= 0 & r$table$on <= 1))
  expect_true(all(r$table$off >= 0 & r$table$off <= 1))
  # the probe was placed over the variant, so its own allele wins
  expect_true(all(r$table$d_eff > 0))
})

test_that("kasp_readout names the assumption for a non-simulated scheme", {
  m <- structure(list(snp = tibble::tibble()), class = "kasp_hrm")
  r <- kasp_readout(m)
  expect_false(r$simulated)
  expect_equal(r$scheme, "hrm")
  expect_true(is.na(r$model))
  expect_null(r$table)
  expect_match(r$note, "melting curve")
})
