#!/usr/bin/env Rscript
# =============================================================================
#  dev/metrics.R  --  run metrics: writer, validator and daily report
# =============================================================================
#  Contract: docs/plans/2026-11-27-integration-3day.md §0.1 and ADR 0011 §5.
#  One JSONL file per run, one object per stage, numbers and ids only.  There is
#  no free text: `prompt`, `content` and `sequence` are not fields, so a file
#  carrying them fails `--check` instead of becoming a debugging note.
#
#      Rscript dev/metrics.R --emit --run baseline-15a3c41 --stage design \
#          --object scenario:J2 --ms 9123 --status ok
#      Rscript dev/metrics.R --check 2026-11-27
#      Rscript dev/metrics.R --report 2026-11-27
#      Rscript dev/metrics.R --baseline 2026-11-27
#
#  Base R + jsonlite only: the gate runs before the package is loaded.
# =============================================================================

STAGES <- c(
  "scenario", "ingest", "design", "specificity", "thermo", "checks",
  "artefact", "mcp", "gate", "repair"
)
STATUSES <- c("ok", "warn", "fail", "skip", "error")
LANGUAGES <- c("en", "ru", "und")
FIELDS <- c(
  "run_id", "stage", "object_id", "duration_ms", "status", "retries",
  "tokens_in", "tokens_out", "cost_usd", "cache_hit", "queue_depth",
  "gpu_util", "language"
)

root_dir <- function() {
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))))
  }
  getwd()
}

metrics_dir <- function(root = root_dir()) file.path(root, "data", "metrics")

arg_value <- function(args, key, default = NULL) {
  i <- match(key, args)
  if (is.na(i)) {
    return(default)
  }
  if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
    return(TRUE)
  }
  args[[i + 1L]]
}

as_int <- function(x, field) {
  if (is.null(x) || is.na(x) || identical(x, "null") || identical(x, "NA")) {
    return(NA_integer_)
  }
  n <- suppressWarnings(as.numeric(x))
  if (is.na(n) || n < 0 || n != floor(n)) {
    stop(field, " must be a non-negative integer", call. = FALSE)
  }
  as.integer(n)
}

as_num <- function(x) {
  if (is.null(x) || is.na(x) || identical(x, "null") || identical(x, "NA")) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x))
}

# --- one row ----------------------------------------------------------------

# Build a row in the contract's field order.  Unknown arguments are ignored on
# purpose: the command line is not the schema, `validate_row()` is.
metrics_row <- function(run_id, stage, object_id, duration_ms,
                        status = "ok", retries = 0L, tokens_in = NA_integer_,
                        tokens_out = NA_integer_, cost_usd = NA_real_,
                        cache_hit = FALSE, queue_depth = 0L, gpu_util = NA_real_,
                        language = "und") {
  list(
    run_id = as.character(run_id),
    stage = as.character(stage),
    object_id = as.character(object_id),
    duration_ms = as_int(duration_ms, "duration_ms"),
    status = as.character(status),
    retries = as_int(retries, "retries"),
    tokens_in = as_int(tokens_in, "tokens_in"),
    tokens_out = as_int(tokens_out, "tokens_out"),
    cost_usd = as_num(cost_usd),
    cache_hit = isTRUE(cache_hit),
    queue_depth = as_int(queue_depth, "queue_depth"),
    gpu_util = as_num(gpu_util),
    language = as.character(language)
  )
}

# A row is valid when it has exactly the contract's keys and each value is in
# range.  `object_id` is an id, not text: no spaces, quotes or long ACGT runs.
validate_row <- function(row, where = "row") {
  problems <- character()
  keys <- names(row)
  if (!setequal(keys, FIELDS)) {
    extra <- setdiff(keys, FIELDS)
    missing <- setdiff(FIELDS, keys)
    if (length(extra)) problems <- c(problems, paste0(where, ": unknown field(s) ", paste(extra, collapse = ", ")))
    if (length(missing)) problems <- c(problems, paste0(where, ": missing field(s) ", paste(missing, collapse = ", ")))
    return(problems)
  }
  add <- function(cond, msg) if (!isTRUE(cond)) problems <<- c(problems, paste0(where, ": ", msg))

  add(length(row$run_id) == 1L && nzchar(row$run_id) && !grepl("[[:space:]]", row$run_id),
    "run_id must be one non-empty token")
  add(length(row$stage) == 1L && row$stage %in% STAGES,
    paste0("stage must be one of ", paste(STAGES, collapse = ", ")))
  add(length(row$object_id) == 1L && grepl("^[A-Za-z0-9_:.@/+<>-]{1,120}$", row$object_id),
    "object_id must be an id (letters, digits, _ : . @ / + < > -), at most 120 chars")
  add(!(nchar(row$object_id) > 60L && grepl("^[ACGTNacgtn]+$", row$object_id)),
    "object_id looks like a nucleotide sequence")
  add(length(row$duration_ms) == 1L && !is.na(row$duration_ms) && row$duration_ms >= 0,
    "duration_ms must be a non-negative integer")
  add(length(row$status) == 1L && row$status %in% STATUSES,
    paste0("status must be one of ", paste(STATUSES, collapse = ", ")))
  add(length(row$retries) == 1L && !is.na(row$retries) && row$retries >= 0,
    "retries must be a non-negative integer")
  add(is.null(row$tokens_in) || (length(row$tokens_in) == 1L && (is.na(row$tokens_in) || row$tokens_in >= 0)),
    "tokens_in must be null or a non-negative integer")
  add(is.null(row$tokens_out) || (length(row$tokens_out) == 1L && (is.na(row$tokens_out) || row$tokens_out >= 0)),
    "tokens_out must be null or a non-negative integer")
  add(is.null(row$cost_usd) || (length(row$cost_usd) == 1L && (is.na(row$cost_usd) || row$cost_usd >= 0)),
    "cost_usd must be null or a non-negative number")
  add(is.logical(row$cache_hit) && length(row$cache_hit) == 1L && !is.na(row$cache_hit),
    "cache_hit must be TRUE or FALSE")
  add(length(row$queue_depth) == 1L && !is.na(row$queue_depth) && row$queue_depth >= 0,
    "queue_depth must be a non-negative integer")
  add(is.null(row$gpu_util) || (length(row$gpu_util) == 1L && (is.na(row$gpu_util) || (row$gpu_util >= 0 && row$gpu_util <= 100))),
    "gpu_util must be null or a percentage in [0, 100]")
  add(length(row$language) == 1L && row$language %in% LANGUAGES,
    paste0("language must be one of ", paste(LANGUAGES, collapse = ", ")))
  problems
}

row_json <- function(row) {
  as.character(jsonlite::toJSON(row, auto_unbox = TRUE, null = "null", na = "null", digits = NA))
}

write_rows <- function(rows, path) {
  bad <- unlist(lapply(seq_along(rows), function(i) validate_row(rows[[i]], paste0(basename(path), ":", i))))
  if (length(bad)) {
    stop(paste(bad, collapse = "\n"), call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(vapply(rows, row_json, character(1)), path, useBytes = TRUE)
  invisible(path)
}

append_row <- function(row, path) {
  bad <- validate_row(row, basename(path))
  if (length(bad)) stop(paste(bad, collapse = "\n"), call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  cat(row_json(row), "\n", file = path, sep = "", append = TRUE)
  invisible(path)
}

read_rows <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lapply(seq_along(lines), function(i) {
    row <- tryCatch(
      jsonlite::fromJSON(lines[[i]], simplifyVector = FALSE),
      error = function(e) structure(list(), bad_line = conditionMessage(e))
    )
    if (length(row) == 0L) {
      list(row = NULL, error = paste0(basename(path), ":", i, ": not JSON (", attr(row, "bad_line"), ")"))
    } else {
      list(row = row, error = NULL)
    }
  })
}

files_for <- function(date, dir = metrics_dir()) {
  if (!dir.exists(dir)) return(character())
  list.files(dir, pattern = paste0("^", date, "-.*[.]jsonl$"), full.names = TRUE)
}

check_date <- function(date, dir = metrics_dir()) {
  files <- files_for(date, dir)
  if (!length(files)) {
    return(list(files = files, rows = 0L, problems = paste0("no metrics file for ", date, " in ", dir)))
  }
  problems <- character()
  total <- 0L
  for (f in files) {
    parsed <- read_rows(f)
    errs <- vapply(parsed, function(p) p$error %||% NA_character_, character(1))
    problems <- c(problems, errs[!is.na(errs)])
    for (p in parsed) {
      if (!is.null(p$row)) {
        total <- total + 1L
        problems <- c(problems, validate_row(p$row, basename(f)))
      }
    }
  }
  list(files = files, rows = total, problems = unique(problems))
}

# --- report -----------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

nearest_rank <- function(x, p) {
  x <- sort(x[!is.na(x)])
  if (!length(x)) return(NA_real_)
  x[[max(1L, ceiling(p * length(x)))]]
}

pct <- function(x, n) if (!n) NA_real_ else 100 * sum(as.logical(x)) / n

rows_as_df <- function(rows) {
  if (!length(rows)) {
    return(data.frame(
      run_id = character(), stage = character(), object_id = character(),
      duration_ms = numeric(), status = character(), retries = numeric(),
      tokens_in = numeric(), tokens_out = numeric(), cost_usd = numeric(),
      cache_hit = logical(), queue_depth = numeric(), gpu_util = numeric(),
      language = character(), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    run_id = vapply(rows, function(r) as.character(r$run_id), character(1)),
    stage = vapply(rows, function(r) as.character(r$stage), character(1)),
    object_id = vapply(rows, function(r) as.character(r$object_id), character(1)),
    duration_ms = vapply(rows, function(r) as.numeric(r$duration_ms), numeric(1)),
    status = vapply(rows, function(r) as.character(r$status), character(1)),
    retries = vapply(rows, function(r) as.numeric(r$retries), numeric(1)),
    tokens_in = vapply(rows, function(r) as.numeric(r$tokens_in %||% NA_real_), numeric(1)),
    tokens_out = vapply(rows, function(r) as.numeric(r$tokens_out %||% NA_real_), numeric(1)),
    cost_usd = vapply(rows, function(r) as.numeric(r$cost_usd %||% NA_real_), numeric(1)),
    cache_hit = vapply(rows, function(r) isTRUE(r$cache_hit), logical(1)),
    queue_depth = vapply(rows, function(r) as.numeric(r$queue_depth), numeric(1)),
    gpu_util = vapply(rows, function(r) as.numeric(r$gpu_util %||% NA_real_), numeric(1)),
    language = vapply(rows, function(r) as.character(r$language), character(1)),
    stringsAsFactors = FALSE
  )
}

repair_stats <- function(log_path = file.path(root_dir(), "docs", "plans", "self-repair-log.md")) {
  if (!file.exists(log_path)) return(list(rows = 0L, resolved = 0L, quarantined = 0L))
  lines <- readLines(log_path, warn = FALSE)
  rows <- grep("^\\| [0-9]", lines, value = TRUE)
  cells <- strsplit(rows, "|", fixed = TRUE)
  status <- vapply(cells, function(x) trimws(x[[9]]), character(1))
  list(
    rows = length(rows),
    fixed = sum(status == "fixed"),
    bypassed = sum(status == "bypassed"),
    quarantined = sum(status == "quarantined"),
    open = sum(status == "open")
  )
}

summary_of <- function(df) {
  stages <- lapply(STAGES, function(s) {
    d <- df[df$stage == s, , drop = FALSE]
    n <- nrow(d)
    list(
      stage = s, n = n,
      p50 = if (n >= 20L) nearest_rank(d$duration_ms, 0.50) else NA_real_,
      p95 = if (n >= 20L) nearest_rank(d$duration_ms, 0.95) else NA_real_,
      fail_share = if (n) pct(d$status %in% c("fail", "error"), n) else NA_real_,
      cost = if (n && !all(is.na(d$cost_usd))) sum(d$cost_usd, na.rm = TRUE) else NA_real_
    )
  })
  routes <- lapply(paste0("J", 1:7), function(j) {
    d <- df[df$object_id == paste0("scenario:", j), , drop = FALSE]
    runs <- length(unique(d$run_id))
    list(
      route = j, runs = runs,
      cost_per_run = if (runs && !all(is.na(d$cost_usd))) sum(d$cost_usd, na.rm = TRUE) / runs else NA_real_
    )
  })
  list(
    stages = stages, routes = routes, rows = nrow(df),
    runs = length(unique(df$run_id)),
    fail_share = if (nrow(df)) pct(df$status %in% c("fail", "error"), nrow(df)) else NA_real_,
    cache_hit_rate = if (nrow(df)) 100 * mean(df$cache_hit) else NA_real_,
    queue_p95 = nearest_rank(df$queue_depth, 0.95)
  )
}

fmt <- function(x, digits = 1) if (is.na(x)) "n/a" else formatC(x, format = "f", digits = digits)

report_md <- function(date, dir = metrics_dir(), baseline = NULL) {
  chk <- check_date(date, dir)
  if (length(chk$problems)) {
    stop(paste(chk$problems, collapse = "\n"), call. = FALSE)
  }
  rows <- unlist(lapply(chk$files, function(f) lapply(read_rows(f), `[[`, "row")), recursive = FALSE)
  df <- rows_as_df(rows)
  s <- summary_of(df)
  rep_stats <- repair_stats()
  failures <- sum(df$status %in% c("fail", "error"))

  base <- NULL
  base_path <- baseline %||% file.path(dir, "baseline.json")
  if (file.exists(base_path)) base <- jsonlite::fromJSON(base_path, simplifyVector = FALSE)

  out <- c(
    paste0("# Run metrics — ", date),
    "",
    paste0("Files: ", length(chk$files), "; rows: ", s$rows, "; runs: ", s$runs, "."),
    ""
  )
  out <- c(out, "## Stages", "", "| stage | n | p50 ms | p95 ms | fail share % | cost USD | delta p95 |", "| --- | --- | --- | --- | --- | --- | --- |")
  for (st in s$stages) {
    delta <- "n/a"
    if (!is.null(base)) {
      b <- Filter(function(x) identical(x$stage, st$stage), base$stages)
      if (length(b) && !is.na(b[[1L]]$p95) && !is.na(st$p95)) {
        delta <- fmt(100 * (st$p95 / b[[1L]]$p95 - 1))
      } else if (is.na(st$p95)) {
        delta <- "not_run"
      }
    }
    note <- if (st$n && st$n < 20L) "insufficient_n" else fmt(st$p50, 0)
    if (is.na(st$p95)) {
      out <- c(out, sprintf("| %s | %d | %s | insufficient_n | %s | %s | %s |", st$stage, st$n, note, fmt(st$fail_share), fmt(st$cost, 2), delta))
    } else {
      out <- c(out, sprintf("| %s | %d | %s | %s | %s | %s | %s |", st$stage, st$n, fmt(st$p50, 0), fmt(st$p95, 0), fmt(st$fail_share), fmt(st$cost, 2), delta))
    }
  }
  out <- c(out, "", "## Routes J1–J7", "", "| route | runs | cost per run USD |", "| --- | --- | --- |")
  for (r in s$routes) {
    out <- c(out, sprintf("| %s | %d | %s |", r$route, r$runs, if (r$runs && !is.na(r$cost_per_run)) fmt(r$cost_per_run, 4) else "not_run"))
  }
  out <- c(
    out, "",
    sprintf("- failure share: %s%% (%d fail/error of %d rows)", fmt(s$fail_share), failures, s$rows),
    sprintf("- cache-hit rate: %s%%", fmt(s$cache_hit_rate)),
    sprintf("- p95 queue depth: %s", fmt(s$queue_p95, 0)),
    sprintf(
      "- self-repair: %d log line(s), fixed %d, bypassed %d, quarantined %d, open %d",
      rep_stats$rows, rep_stats$fixed %||% 0L, rep_stats$bypassed %||% 0L,
      rep_stats$quarantined %||% 0L, rep_stats$open %||% 0L
    ),
    ""
  )
  if (is.null(base)) {
    out <- c(out, "_No baseline.json: deltas are not_run, not zero._", "")
  }
  paste(out, collapse = "\n")
}

# --- CLI --------------------------------------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  dir <- arg_value(args, "--dir", metrics_dir())
  date <- arg_value(args, "--date", format(Sys.Date()))

  if ("--emit" %in% args || "--append" %in% args) {
    run <- arg_value(args, "--run", "run")
    row <- metrics_row(
      run_id = run,
      stage = arg_value(args, "--stage", "gate"),
      object_id = arg_value(args, "--object", "gate"),
      duration_ms = arg_value(args, "--ms", "0"),
      status = arg_value(args, "--status", "ok"),
      retries = arg_value(args, "--retries", "0"),
      tokens_in = arg_value(args, "--tokens-in", NA),
      tokens_out = arg_value(args, "--tokens-out", NA),
      cost_usd = arg_value(args, "--cost", NA),
      cache_hit = isTRUE(arg_value(args, "--cache-hit", FALSE)),
      queue_depth = arg_value(args, "--queue", "0"),
      gpu_util = arg_value(args, "--gpu", NA),
      language = arg_value(args, "--language", "und")
    )
    path <- file.path(dir, paste0(date, "-", run, ".jsonl"))
    append_row(row, path)
    cat("ok metrics ", path, "\n")
    return(invisible(NULL))
  }

  if ("--check" %in% args) {
    d <- arg_value(args, "--check", format(Sys.Date()))
    if (identical(d, TRUE)) d <- format(Sys.Date())
    chk <- check_date(d, dir)
    if (length(chk$problems)) {
      cat("FAIL metrics: ", length(chk$problems), " problem(s)\n", sep = "")
      cat(paste0("  - ", chk$problems, collapse = "\n"), "\n")
      quit(status = 1L)
    }
    cat("ok metrics ", d, ": ", length(chk$files), " file(s), ", chk$rows, " row(s)\n", sep = "")
    return(invisible(NULL))
  }

  if ("--report" %in% args) {
    d <- arg_value(args, "--report", format(Sys.Date()))
    if (identical(d, TRUE)) d <- format(Sys.Date())
    md <- report_md(d, dir, arg_value(args, "--baseline"))
    path <- file.path(dir, paste0("metrics-", d, ".md"))
    writeLines(md, path)
    cat("ok report ", path, "\n")
    return(invisible(NULL))
  }

  if ("--baseline" %in% args) {
    d <- arg_value(args, "--baseline", format(Sys.Date()))
    if (identical(d, TRUE)) d <- format(Sys.Date())
    chk <- check_date(d, dir)
    if (length(chk$problems)) stop(paste(chk$problems, collapse = "\n"), call. = FALSE)
    rows <- unlist(lapply(chk$files, function(f) lapply(read_rows(f), `[[`, "row")), recursive = FALSE)
    s <- summary_of(rows_as_df(rows))
    s$date <- d
    s$commit <- arg_value(args, "--commit", NA)
    writeLines(as.character(jsonlite::toJSON(s, auto_unbox = TRUE, na = "null", digits = NA)),
      file.path(dir, "baseline.json")
    )
    cat("ok baseline ", file.path(dir, "baseline.json"), "\n")
    return(invisible(NULL))
  }

  cat("usage: metrics.R --emit|--check|--report|--baseline [date]\n")
  quit(status = 2L)
}

if (sys.nframe() == 0L) main()
