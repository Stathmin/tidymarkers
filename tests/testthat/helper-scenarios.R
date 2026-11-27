# Scenario helpers: a committed fixture, the typed step, and an optional
# per-step metrics row.  A scenario run is a normal test run; setting
# TIDYMARKERS_METRICS=<run> also records one JSONL row per step (contract
# §0.1).  The writer lives in dev/metrics.R, which is not shipped, so an
# installed package run simply does not emit.

scenario_file <- function(name) {
  testthat::test_path("fixtures", "scenarios", name)
}

scenario_metrics_writer <- function() {
  for (up in c("../..", "../../..")) {
    path <- file.path(up, "dev", "metrics.R")
    if (file.exists(path)) {
      env <- new.env(parent = globalenv())
      sys.source(path, envir = env)
      return(list(env = env, root = normalizePath(file.path(up))))
    }
  }
  NULL
}

#' Record one scenario step
#'
#' Returns the step's value, so a test reads as the pipeline it proves.  The
#' metrics row is written only when `TIDYMARKERS_METRICS` names a run.
#' @param scenario scenario id (`"J1"`, `"J2"`, ...).
#' @param stage metrics stage (see the contract).
#' @param object_id stable id; never a path or a sequence.
#' @param code the step, evaluated once.
#' @keywords internal
#' @noRd
scenario_step <- function(scenario, stage, object_id, code) {
  run <- Sys.getenv("TIDYMARKERS_METRICS")
  t0 <- proc.time()[[3L]]
  out <- force(code)
  ms <- round(1000 * (proc.time()[[3L]] - t0))
  if (nzchar(run)) {
    w <- scenario_metrics_writer()
    if (!is.null(w)) {
      row <- w$env$metrics_row(
        run_id = run, stage = stage,
        object_id = paste0("scenario:", scenario, ":", object_id),
        duration_ms = ms, status = "ok", language = "en"
      )
      w$env$append_row(
        row, file.path(w$root, "data", "metrics", paste0(Sys.Date(), "-", run, ".jsonl"))
      )
    }
  }
  out
}

#' Assert the typed class a scenario step promised
#' @keywords internal
#' @noRd
expect_scenario_class <- function(object, class) {
  ok <- inherits(object, class)
  # the ingest tables are semantic classes (`variants`) carried by a tibble,
  # not an S3 class of their own
  if (!ok && class == "variants") ok <- inherits(object, "data.frame")
  # a FASTA without an index is a plain DNAStringSet; `genome` is the semantic class
  if (!ok && class == "genome") ok <- inherits(object, "DNAStringSet")
  testthat::expect_true(ok, info = paste("expected", class, "got", paste(class(object), collapse = "/")))
}
