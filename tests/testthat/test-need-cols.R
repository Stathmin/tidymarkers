# `need_cols()` is the one required-column guard: one message shape, one pair of
# classes, and the missing columns in the caller's order.

test_that("need_cols aborts with the contract and missing-column classes", {
  err <- rlang::catch_cnd(need_cols(data.frame(chrom = "1"), c("chrom", "pos", "ref")))
  expect_s3_class(err, "tidymarkers_contract_error")
  expect_s3_class(err, "tidymarkers_missing_columns")
  expect_identical(
    conditionMessage(err),
    "`x` is missing column(s): pos, ref"
  )
})

test_that("need_cols names its argument and keeps the given column order", {
  err <- rlang::catch_cnd(need_cols(list(a = 1), c("z", "a", "m"), arg = "loci"))
  expect_identical(
    conditionMessage(err),
    "`loci` is missing column(s): z, m"
  )
})

test_that("need_cols passes a complete frame, tibble or named list", {
  expect_true(need_cols(data.frame(chrom = "1", pos = 1L), c("chrom", "pos")))
  expect_true(need_cols(tibble::tibble(chrom = "1", pos = 1L), c("chrom", "pos")))
  expect_true(need_cols(list(chrom = "1", pos = 1L), c("chrom", "pos")))
})

test_that("need_cols returns invisibly", {
  x <- tibble::tibble(chrom = "1")
  expect_invisible(need_cols(x, "chrom"))
})
