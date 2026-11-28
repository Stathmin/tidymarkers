# Metrics writer/validator: a repair row carries a measured duration or it
# does not exist (ZV14, ROADMAP.md).  dev/metrics.R is not shipped, so the
# script is sourced into a bare environment instead of calling library().

metrics_writer <- function() {
  path <- testthat::test_path("..", "..", "dev", "metrics.R")
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  env
}

test_that("--emit refuses a row with no duration instead of writing 0 ms", {
  w <- metrics_writer()
  dir <- tempfile("metrics-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  script <- testthat::test_path("..", "..", "dev", "metrics.R")
  res <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      shQuote(script), "--emit", "--dir", shQuote(dir), "--date", "2026-01-01",
      "--run", "repair-run", "--stage", "repair", "--object", "R-1",
      "--status", "fixed"
    ),
    stdout = TRUE, stderr = TRUE
  ))
  expect_false(is.null(attr(res, "status"))) # non-zero exit
  expect_match(paste(res, collapse = " "), "unmeasured duration is not 0 ms", fixed = TRUE)
  expect_false(dir.exists(file.path(dir, "2026-01-01-repair-run.jsonl")))
})

test_that("a measured repair row reaches the report with its own p50/p95", {
  w <- metrics_writer()
  repair_log <- testthat::test_path("..", "..", "docs", "plans", "self-repair-log.md")
  dir <- tempfile("metrics-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  rows <- lapply(1:20, function(i) {
    w$metrics_row(
      run_id = "repair-run", stage = "repair", object_id = paste0("R-", i),
      duration_ms = i * 1000, status = "ok", language = "en"
    )
  })
  w$write_rows(rows, file.path(dir, "2026-01-01-repair-run.jsonl"))
  expect_length(w$check_date("2026-01-01", dir)$problems, 0L)

  md <- w$report_md(
    "2026-01-01", dir,
    baseline = file.path(dir, "none.json"), repair_log = repair_log
  )$md
  expect_match(md, "| repair | 20 |", fixed = TRUE)
  expect_false(grepl("repair rows: 0 of", md, fixed = TRUE))
})

test_that("the report names the untimed log lines instead of showing zero repairs", {
  w <- metrics_writer()
  repair_log <- testthat::test_path("..", "..", "docs", "plans", "self-repair-log.md")
  dir <- tempfile("metrics-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  row <- w$metrics_row(
    run_id = "scenario-run", stage = "scenario", object_id = "scenario:J1",
    duration_ms = 10, status = "ok", language = "en"
  )
  w$write_rows(list(row), file.path(dir, "2026-01-01-scenario-run.jsonl"))
  md <- w$report_md(
    "2026-01-01", dir,
    baseline = file.path(dir, "none.json"), repair_log = repair_log
  )$md
  expect_match(
    md,
    sprintf("repair rows: 0 of %d log line(s)", w$repair_stats(repair_log)$rows),
    fixed = TRUE
  )
})

test_that("a malformed or unknown journal row is named, not counted or fatal", {
  w <- metrics_writer()
  log <- tempfile("repair-log-", fileext = ".md")
  writeLines(c(
    "| time (UTC) | blocker | class | action | result | commit | attempt | status | next task |",
    "| ---------- | ------- | ----- | ------ | ------ | ------ | ------- | ------ | --------- |",
    "| 2026-01-01 00:00 | b | logic | a | fixed: x | abc | 1 | fixed | next |",
    "| 2026-01-01 00:01 | a short row |",
    "| 2026-01-01 00:02 | b | logic | a | fixed: x | abc | 1 | fixedd | next |",
    # an escaped pipe is markdown, not a column break (six real rows did this)
    "| 2026-01-01 00:03 | b | logic | `a \\| b` | fixed: x | abc | 1 | fixed | next |"
  ), log)

  st <- w$repair_stats(log)
  expect_equal(st$rows, 4L)
  expect_equal(st$fixed, 2L)
  expect_equal(st$malformed, 1L)
  expect_equal(st$unknown, 1L)
})

test_that("the validator enforces the unit the contract names", {
  w <- metrics_writer()
  row <- w$metrics_row(
    run_id = "x", stage = "ingest", object_id = "i1",
    duration_ms = 10, language = "en"
  )
  expect_length(w$validate_row(row), 0L)

  # §0.1 of the plan: duration_ms, retries, queue_depth, tokens are integers
  for (field in c("duration_ms", "retries", "queue_depth", "tokens_in", "tokens_out")) {
    bad <- row
    bad[[field]] <- 1.5
    expect_match(
      paste(w$validate_row(bad), collapse = " "),
      paste0(field, " must be .*integer"),
      info = field
    )
  }

  # a cost in USD is a finite number: jsonlite reads 1e400 as Inf, and Inf
  # passed `>= 0` and then made the stage cost Inf in the report
  bad <- row
  bad$cost_usd <- Inf
  expect_match(
    paste(w$validate_row(bad), collapse = " "),
    "cost_usd must be null or a finite non-negative number"
  )

  # gpu_util is a fraction in [0, 1], not a percentage
  ok <- row
  ok$gpu_util <- 0.75
  expect_length(w$validate_row(ok), 0L)
  bad <- row
  bad$gpu_util <- 100
  expect_match(
    paste(w$validate_row(bad), collapse = " "),
    "gpu_util must be null or a fraction in \\[0, 1\\]"
  )
})

test_that("a baseline stage without a p95 does not stop the report", {
  w <- metrics_writer()
  repair_log <- testthat::test_path("..", "..", "docs", "plans", "self-repair-log.md")
  dir <- tempfile("metrics-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  rows <- lapply(1:3, function(i) {
    w$metrics_row(
      run_id = "run", stage = "ingest", object_id = paste0("i", i),
      duration_ms = i * 10, status = "ok", language = "en"
    )
  })
  w$write_rows(rows, file.path(dir, "2026-01-01-run.jsonl"))
  writeLines(
    '{"stages":[{"stage":"ingest","n":1,"p50":1,"p95":null}]}',
    file.path(dir, "baseline.json")
  )

  md <- w$report_md(
    "2026-01-01", dir,
    baseline = file.path(dir, "baseline.json"), repair_log = repair_log
  )$md
  expect_match(md, "| ingest | 3 |", fixed = TRUE)
})
