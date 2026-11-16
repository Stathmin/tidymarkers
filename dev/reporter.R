# =============================================================================
#  dev/reporter.R  --  the pure part of the test timing report
# =============================================================================
#  Sourced by dev/test.R and by tests/testthat/test-dev-reporter.R; it has no
#  side effects and reads no state, the clock or the filesystem.  Every function
#  here takes the rows it renders and returns text, so the formatting contract is
#  testable without running a test suite.
#
#  The contract: a measurement that did not happen prints as `n/a` (a cell) or
#  `—` (a memory column that exists for other rows); it never prints `NA`, never
#  a silent zero.  A run with no serial phase drops the memory column entirely
#  instead of filling it with 132 identical dashes.
# =============================================================================

fmt_seconds <- function(s) {
  if (length(s) != 1L || is.na(s)) {
    return("n/a")
  }
  if (s < 60) sprintf("%.1fs", s) else sprintf("%dm %02ds", floor(s / 60), round(s) %% 60)
}

fmt_stamp <- function(x) {
  if (length(x) != 1L || is.na(x)) "n/a" else x
}

fmt_mem <- function(x) {
  if (length(x) != 1L || is.na(x)) "—" else sprintf("%+.1f", x)
}

has_mem <- function(rows) {
  any(vapply(rows, function(r) length(r$mem_delta_mb) == 1L && !is.na(r$mem_delta_mb), logical(1)))
}

slow_files_lines <- function(rows) {
  if (!length(rows)) {
    return(character(0))
  }
  mem <- has_mem(rows)
  head <- if (mem) {
    c(
      "| file | started | duration | pass | fail | mem delta (MB, serial phase) |",
      "|---|---|---|---|---|---|"
    )
  } else {
    c("| file | started | duration | pass | fail |", "|---|---|---|---|---|")
  }
  body <- vapply(rows, function(r) {
    cells <- sprintf(
      "| %s | %s | %s | %d | %d", r$file, fmt_stamp(r$started_at),
      fmt_seconds(r$duration_s), r$pass, r$fail
    )
    if (mem) paste0(cells, " | ", fmt_mem(r$mem_delta_mb), " |") else paste0(cells, " |")
  }, character(1))
  c("## Slowest files", "", head, body, "")
}

slow_tests_lines <- function(rows) {
  if (!length(rows)) {
    return(character(0))
  }
  c(
    "## Slowest tests", "", "| file | test | started | duration |", "|---|---|---|---|",
    vapply(rows, function(r) {
      sprintf(
        "| %s | %s | %s | %s |", r$file, gsub("|", "\\|", r$test, fixed = TRUE),
        fmt_stamp(r$started_at), fmt_seconds(r$duration_s)
      )
    }, character(1)),
    ""
  )
}

timing_gaps_line <- function(rows) {
  gaps <- sum(vapply(rows, function(r) is.na(r$duration_s) || is.na(r$started_at), logical(1)))
  if (!gaps) {
    return(character(0))
  }
  sprintf(
    "- timing gaps: %d of %d file(s) have no start/duration (events were not streamed)",
    gaps, length(rows)
  )
}
