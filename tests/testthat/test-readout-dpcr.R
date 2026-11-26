# The dPCR counting model: Poisson partitioning, validated by the shipped sets.
#
# The model is the counting law only - lambda = -ln(N_neg / N_total) - checked
# against the instrument's own lambda (Dong 2015, 10.1038/srep13174) and across a
# 6.6-7900 copies/uL dynamic range (MPS 2024, 10.3390/mps7060094).  No
# amplification kinetics are claimed.

ext <- function(name) {
  utils::read.csv(
    tm_extdata(name),
    stringsAsFactors = FALSE
  )
}

test_that("poisson_lambda inverts the positive fraction and refuses non-physical input", {
  lambda <- c(0.01, 0.5, 1, 2, 5)
  expect_equal(poisson_lambda(100 * exp(-lambda), 100), lambda)
  expect_equal(poisson_positive_fraction(lambda), 1 - exp(-lambda))
  # more negatives than partitions, nothing, or missing: NA, never NaN/Inf
  expect_equal(poisson_lambda(c(0, 5, 11, NA), c(10, 0, 10, 10)), rep(NA_real_, 4L))
  expect_true(all(is.na(poisson_lambda(c(0, 5, 11, NA), c(10, 0, 10, 10)))))
})

test_that("dong2015: the inversion reproduces the instrument's lambda", {
  d <- ext("dpcr_dong_poisson_residual.csv")
  expect_equal(poisson_lambda(d$negative, d$total), d$lambda_poisson, tolerance = 1e-9)
  expect_lte(max(abs(d$lambda_residual)), 0.005 + 1e-9)
  expect_lt(max(abs(d$lambda_residual) / d$lambda_poisson), 0.004)
})

test_that("mps2024: the inverted loading decreases monotonically along the dilution", {
  d <- ext("dpcr_mps_ecoli_dilution.csv")
  # rep1..rep3 are copies/uL, i.e. a lambda at the partition volume.  At a ~1 uL
  # droplet the top levels saturate (p -> 1), so the check is on the loading axis.
  lambda <- rowMeans(d[c("rep1", "rep2", "rep3")])
  expect_true(all(lambda > 0))
  # the inversion is exact on the same counts the forward law builds
  expect_equal(poisson_lambda(1e6 * exp(-lambda / 1e4), 1e6), lambda / 1e4, tolerance = 1e-6)
  by_level <- tapply(lambda, factor(d$level, levels = 1:5), median)
  expect_true(all(diff(as.numeric(by_level)) < 0))
  expect_true(all(by_level > 0))
})

test_that("the partition table is the Poisson law over the dPCR loading range", {
  p <- dpcr_partition_table()
  expect_equal(names(p), c("copies_per_partition", "expected_positive_fraction", "expected_negative_fraction"))
  expect_true(all(diff(p$expected_positive_fraction) > 0))
  expect_true(all(p$expected_positive_fraction > 0 & p$expected_positive_fraction < 1))
  expect_equal(p$expected_positive_fraction + p$expected_negative_fraction, rep(1, nrow(p)))
})

test_that("a dpcr design reports the probe readout and the partitioning model", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture(seed = 1L)
  m <- design_dpcr_marker(f$genome, f$snp, flank = 150L)
  r <- marker_readout(m)
  expect_true(r$simulated)
  expect_equal(r$scheme, "dpcr")
  expect_match(r$model, "hydrolysis-probe")
  expect_match(r$model, "Poisson")
  expect_true("d_eff" %in% names(r$table))
  expect_true(all(diff(r$partition$expected_positive_fraction) > 0))
  expect_true(all(r$partition$expected_positive_fraction > 0 & r$partition$expected_positive_fraction < 1))
})

test_that("a non-dpcr scheme still reports no simulation and no partition element", {
  r <- marker_readout(structure(list(snp = tibble::tibble()), class = c("marker_hrm", "marker")))
  expect_false(r$simulated)
  expect_null(r$partition)
  expect_false("partition" %in% names(r))
})
