# Poisson-model anchor: 10.1038/srep13174 (PMC4548226).
#
# The package converts a positive-partition fraction to a concentration with
# lambda = -ln(1 - p).  These tests assert that the shipped extract still
# carries the measurement that validates that conversion, and that the
# per-replicate arithmetic is internally consistent.

ext <- function(name) {
  utils::read.csv(
    tm_extdata(name),
    stringsAsFactors = FALSE
  )
}

test_that("dong2015: every replicate's counts sum to its printed total", {
  d <- ext("dpcr_dong_poisson_residual.csv")
  expect_equal(nrow(d), 60L)
  expect_equal(sum(d$negative + d$positive != d$total), 0L)
  expect_equal(sum(d$total != d$reported_total), 0L)
  expect_setequal(
    unique(d$platform),
    c("QuantStudio 12k", "BioMark", "QX100", "RainDrop")
  )
  expect_true(all(table(d$platform) == 15L)) # 5 replicates x 3 vials
  expect_true(all(table(d$vial) == 20L))
})

test_that("dong2015: lambda recomputed from the counts reproduces the printed mean", {
  d <- ext("dpcr_dong_poisson_residual.csv")
  recomputed <- -log(d$negative / d$total)
  # the CSV stores lambda to ~15 significant digits
  expect_lt(max(abs(recomputed - d$lambda_poisson)), 1e-9)
  expect_lt(max(abs(d$lambda_residual)), 0.005 + 1e-9)
  # The residual is not just rounding: the printed mean is given to 2 d.p., so a
  # pure rounding artefact would be bounded by 0.005/2.  Check the residual is
  # larger than that on at least some rows, i.e. it measures something.
  expect_true(any(abs(d$lambda_residual) > 0.0025))
  # and it is small relative to lambda itself
  expect_lt(max(abs(d$lambda_residual) / d$lambda_poisson), 0.004)
})

test_that("dong2015: the per-platform summary reproduces Table 2's lambda", {
  s <- ext("dpcr_dong_poisson_summary.csv")
  expect_equal(nrow(s), 4L)
  expect_equal(sum(s$n), 60L)
  expect_equal(sum(s$lambda_reported - s$table2_lambda > 0.01), 0L)
  expect_equal(sum(s$lambda_poisson - s$table2_lambda > 0.01), 0L)
  # the summary must agree with the rows it summarises
  d <- ext("dpcr_dong_poisson_residual.csv")
  for (i in seq_len(nrow(s))) {
    v <- d[d$platform == s$platform[i], ]
    expect_equal(round(mean(v$reported_mean), 4), s$lambda_reported[i])
    expect_equal(round(mean(v$total)), s$partitions_mean[i])
  }
  # partition counts span the whole dPCR design space in the paper
  expect_gt(max(s$partitions_mean) / min(s$partitions_mean), 1e4)
})
