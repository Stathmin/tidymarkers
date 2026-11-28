#!/usr/bin/env Rscript

#  repair.R  --  append one self-repair journal row and emit its measured metric
#
#  The journal says what happened; the repair metric says how long it took.
#  Both are written from the same start/end pair, so a fixed blocker leaves a
#  timed row instead of a markdown-only note with no MTTR.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  "dev"
}
arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) {
    return(default)
  }
  if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
    stop(name, " needs a value", call. = FALSE)
  }
  args[[i + 1L]]
}
need <- function(name) {
  value <- arg(name)
  if (is.null(value) || !nzchar(value)) {
    stop(name, " is required", call. = FALSE)
  }
  value
}

classes <- c("transient", "environment", "data", "access", "logic", "spec", "external")
statuses <- c("open", "fixed", "bypassed", "quarantined", "watchdog")
metric_status <- c(
  open = "fail", fixed = "ok", bypassed = "warn",
  quarantined = "warn", watchdog = "warn"
)

blocker <- need("--blocker")
class <- need("--class")
action <- need("--action")
result <- need("--result")
commit <- need("--commit")
attempt <- arg("--attempt", "1")
status <- need("--status")
next_task <- arg("--next", "—")
start <- need("--start")
end <- need("--end")
log <- arg("--log", file.path("docs", "plans", "self-repair-log.md"))
dir <- arg("--dir", file.path("data", "metrics"))
run <- arg("--run", "repair")
date <- arg("--date", format(Sys.Date()))

if (!class %in% classes) {
  stop("--class must be one of ", paste(classes, collapse = ", "), call. = FALSE)
}
if (!status %in% statuses) {
  stop("--status must be one of ", paste(statuses, collapse = ", "), call. = FALSE)
}

as_time <- function(x, name) {
  out <- as.POSIXct(x, tz = "UTC", format = "%Y-%m-%d %H:%M")
  if (is.na(out)) {
    out <- as.POSIXct(x, tz = "UTC", format = "%Y-%m-%dT%H:%M")
  }
  if (is.na(out)) {
    stop(name, " must be a UTC time such as 2026-11-27 18:20", call. = FALSE)
  }
  out
}
started <- as_time(start, "--start")
ended <- as_time(end, "--end")
ms <- as.numeric(difftime(ended, started, units = "secs")) * 1000
if (!is.finite(ms) || ms < 0) {
  stop("--end must not be before --start; an unmeasured duration is not 0 ms", call. = FALSE)
}

cell <- function(x) gsub("|", "\\|", gsub("[\r\n]+", " ", x), fixed = TRUE)
row <- paste0(
  "| ", format(ended, "%Y-%m-%d %H:%M"), " | ", cell(blocker), " | ", class,
  " | ", cell(action), " | ", cell(result), " | ", commit, " | ", attempt,
  " | ", status, " | ", cell(next_task), " |"
)

lines <- readLines(log, warn = FALSE)
table <- grep("^\\|", lines)
if (!length(table)) {
  stop("no markdown table in ", log, call. = FALSE)
}
writeLines(append(lines, row, after = max(table)), log)

env <- new.env(parent = globalenv())
sys.source(file.path(script_dir, "metrics.R"), envir = env)
# the journal keeps the human text; the metric key may only be an id
object <- gsub("[^A-Za-z0-9_:.@/+<>-]", "_", blocker)
object <- substr(object, 1L, 120L)
metric <- env$metrics_row(
  run_id = run, stage = "repair", object_id = object, duration_ms = ms,
  status = metric_status[[status]], language = "en"
)
path <- file.path(dir, paste0(date, "-", run, ".jsonl"))
env$append_row(metric, path)
cat("ok repair ", path, " (", round(ms), " ms)\n", sep = "")
