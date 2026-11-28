test_that("is_whole_number answers alone, without an NA in the way", {
  expect_true(is_whole_number(1L))
  expect_true(is_whole_number(4))
  expect_true(is_whole_number(.Machine$integer.max))
  expect_true(is_whole_number(0, min = 0L))

  expect_false(is_whole_number(0)) # a positive count starts at one
  expect_false(is_whole_number(4.9))
  expect_false(is_whole_number(c(4, 5)))
  expect_false(is_whole_number(numeric(0)))
  expect_false(is_whole_number(NA_real_))
  expect_false(is_whole_number(NaN))
  expect_false(is_whole_number(Inf))
  expect_false(is_whole_number(-Inf))
  expect_false(is_whole_number("4"))
  expect_false(is_whole_number(TRUE))
  expect_false(is_whole_number(.Machine$integer.max + 1))
})
