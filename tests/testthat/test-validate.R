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

test_that("bad_position refuses a spelling that is not a whole number", {
  expect_false(any(bad_position(c(1L, 250, "250"))))
  expect_false(any(bad_position(c("+250", " 250 "))))

  # as.integer() truncates or parses these, so only the spelling can refuse them
  expect_true(all(bad_position(c("250.9", "1e3", "2.0", ""))))
  expect_true(all(bad_position(c(250.9, 0L, -1L, NA_real_, Inf))))
  expect_true(all(bad_position(c("abc", NA_character_))))
})
