# The dev test reporter's formatting contract.  A measurement that did not happen
# must read as `n/a` (or `—` inside a memory column that exists for other rows),
# never as `NA` and never as a silent zero.  `dev/reporter.R` is pure so this file
# needs no suite run; R CMD check does not ship `dev/`, so it skips there.

reporter_src <- test_path("..", "..", "dev", "reporter.R")

skip_if_not(file.exists(reporter_src), "dev/reporter.R is not shipped in the build")
source(reporter_src)

file_row <- function(file = "test-x.R", started_at = "2026-01-01 00:00:00",
                     duration_s = 1.5, pass = 3L, fail = 0L, mem_delta_mb = NA_real_) {
  list(
    file = file, started_at = started_at, duration_s = duration_s,
    pass = pass, fail = fail, mem_delta_mb = mem_delta_mb
  )
}

test_that("a missing second count prints as n/a, never as NA", {
  expect_identical(fmt_seconds(NA_real_), "n/a")
  expect_identical(fmt_seconds(4.2), "4.2s")
  expect_identical(fmt_seconds(68), "1m 08s")
  expect_false(grepl("NA", fmt_seconds(NA_real_), fixed = TRUE))
})

test_that("a missing start stamp prints as n/a", {
  expect_identical(fmt_stamp(NA_character_), "n/a")
  expect_identical(fmt_stamp("2026-01-01 00:00:00"), "2026-01-01 00:00:00")
})

test_that("the memory column is dropped when nothing was measured", {
  rows <- list(file_row(), file_row("test-y.R", duration_s = 3))
  lines <- slow_files_lines(rows)
  expect_true(any(grepl("^## Slowest files$", lines)))
  expect_false(any(grepl("mem delta", lines)))
  expect_length(grep("^\\|", lines), 4L) # header, rule, two rows
  expect_false(any(grepl("NA", lines, fixed = TRUE)))
})

test_that("a partly measured memory column prints an em dash, not NA", {
  rows <- list(file_row(), file_row("test-y.R", mem_delta_mb = -1.34))
  lines <- slow_files_lines(rows)
  expect_true(any(grepl("mem delta (MB, serial phase)", lines, fixed = TRUE)))
  expect_true(any(grepl("| — |", lines, fixed = TRUE)))
  expect_true(any(grepl("| -1.3 |", lines, fixed = TRUE)))
  expect_false(any(grepl("NA", lines, fixed = TRUE)))
})

test_that("a row with no timer at all still renders both cells as n/a", {
  rows <- list(file_row(started_at = NA_character_, duration_s = NA_real_))
  lines <- slow_files_lines(rows)
  expect_true(any(grepl("| n/a | n/a |", lines, fixed = TRUE)))
})

test_that("the slowest-test table escapes a pipe in a test name", {
  rows <- list(list(
    file = "test-x.R", test = "a | b", started_at = "2026-01-01 00:00:00",
    duration_s = 2
  ))
  lines <- slow_tests_lines(rows)
  expect_true(any(grepl("a \\| b", lines, fixed = TRUE)))
  expect_false(any(grepl("NA", lines, fixed = TRUE)))
})

test_that("no rows print no section at all", {
  expect_identical(slow_files_lines(list()), character(0))
  expect_identical(slow_tests_lines(list()), character(0))
  expect_identical(timing_gaps_line(list()), character(0))
})

test_that("timing gaps are counted and reported by name", {
  rows <- list(
    file_row(),
    file_row("test-y.R", started_at = NA_character_, duration_s = NA_real_)
  )
  gap <- timing_gaps_line(rows)
  expect_length(gap, 1L)
  expect_match(gap, "^- timing gaps: 1 of 2 file")
})
