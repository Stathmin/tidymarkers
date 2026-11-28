# A self-repair row is only useful if its duration is measured (ZV14).  The
# helper writes the journal row and the repair metric from one start/end pair.

repair_script <- function() {
  path <- testthat::test_path("..", "..", "dev", "repair.R")
  # dev/ is repo tooling and is not part of the built package, so R CMD check
  # has no copy; only a dev checkout can exercise the script.
  testthat::skip_if_not(file.exists(path), "dev/repair.R is not shipped in the built package")
  path
}

repair_run <- function(args) {
  quoted <- vapply(args, shQuote, character(1))
  suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(repair_script()), quoted),
    stdout = TRUE, stderr = TRUE
  ))
}

repair_log_fixture <- function(dir) {
  path <- file.path(dir, "self-repair-log.md")
  writeLines(c(
    "# Self-repair log",
    "",
    "| time (UTC) | blocker | class | action | result | commit | attempt | status | next task |",
    "| ---------- | ------- | ----- | ------ | ------ | ------ | ------- | ------ | --------- |",
    "| 2026-01-01 10:00 | old row | logic | fixed | done | abc123 | 1 | fixed | next |",
    "",
    "Status values: open, fixed, bypassed, quarantined, watchdog."
  ), path)
  path
}

test_that("a repair row is appended with its measured duration", {
  dir <- withr::local_tempdir()
  log <- repair_log_fixture(dir)
  before <- readLines(log)

  res <- repair_run(c(
    "--blocker", "an example blocker | with a pipe", "--class", "logic",
    "--action", "the row is checked", "--result", "fixed: test",
    "--commit", "def456", "--attempt", "2", "--status", "fixed",
    "--next", "next audit", "--start", "2026-01-01 10:00",
    "--end", "2026-01-01 10:30", "--log", log, "--dir", dir,
    "--run", "repair-run", "--date", "2026-01-01"
  ))
  expect_null(attr(res, "status"))

  after <- readLines(log)
  expect_equal(length(after), length(before) + 1L)
  # the new row stays inside the table, above the prose
  expect_match(after[[6L]], "an example blocker \\| with a pipe", fixed = TRUE)
  expect_match(after[[6L]], "| 2 | fixed | next audit |", fixed = TRUE)
  expect_match(after[[7L]], "^$")
  expect_match(after[[8L]], "^Status values")

  rows <- jsonlite::fromJSON(readLines(file.path(dir, "2026-01-01-repair-run.jsonl")))
  expect_length(rows$stage, 1L)
  expect_identical(rows$stage[[1L]], "repair")
  expect_equal(rows$duration_ms[[1L]], 30 * 60 * 1000)
  expect_identical(rows$status[[1L]], "ok")
  expect_identical(rows$object_id[[1L]], "an_example_blocker___with_a_pipe")
})

test_that("a repair row without an end time writes nothing", {
  dir <- withr::local_tempdir()
  log <- repair_log_fixture(dir)
  before <- readLines(log)

  res <- repair_run(c(
    "--blocker", "no end", "--class", "logic", "--action", "none",
    "--result", "none", "--commit", "abc", "--status", "fixed",
    "--start", "2026-01-01 10:00", "--log", log, "--dir", dir
  ))
  expect_false(is.null(attr(res, "status")))
  expect_match(paste(res, collapse = " "), "--end is required", fixed = TRUE)
  expect_identical(readLines(log), before)
  expect_length(list.files(dir, pattern = "[.]jsonl$"), 0L)
})

test_that("a repair row with an end before its start is refused", {
  dir <- withr::local_tempdir()
  log <- repair_log_fixture(dir)

  res <- repair_run(c(
    "--blocker", "backwards", "--class", "logic", "--action", "none",
    "--result", "none", "--commit", "abc", "--status", "fixed",
    "--start", "2026-01-01 11:00", "--end", "2026-01-01 10:00",
    "--log", log, "--dir", dir
  ))
  expect_false(is.null(attr(res, "status")))
  expect_match(paste(res, collapse = " "), "must not be before", fixed = TRUE)
  expect_length(list.files(dir, pattern = "[.]jsonl$"), 0L)
})
