#!/usr/bin/env Rscript
# =============================================================================
#  dev/test.R  --  the routine test loop
# =============================================================================
#  Three modes, from cheapest to complete:
#
#    Rscript dev/test.R                  # lean: shared contracts only (~1 min)
#    Rscript dev/test.R --scheme crispr  # one family + shared contracts
#    Rscript dev/test.R --full           # every test, external cross-checks on
#    Rscript dev/test.R <regex>          # legacy: raw testthat filter
#
#  The shared contracts (tests/testthat/test-{criteria,decision,kasp-checks,
#  plot-methods,technology}.R) encode the invariants every scheme inherits, so
#  they are the minimum gate for a routine change.  `--scheme <regex>` adds the
#  test files whose names match the regex (e.g. crispr, lamp|rpa, indel) on top
#  of them; the regex is matched against the file name, never the file body.
#  `--full` runs every test with TIDYMARKERS_FULL_TESTS=1, which turns on the
#  external cross-checks (BLAST+, Rsubread, rprimer, annotation DBs, real SOCK
#  workers); R CMD check also sets that flag, so the full path is never weakened.
#
#  Lean and scheme modes run files in parallel (2 workers).  Full mode is serial:
#  the SOCK tests already spawn worker processes, and nesting them under
#  testthat's workers crashes the pool.
#
#  Every run is teed to dev/logs/test-<mode>-<stamp>.log (dev/logs/latest.log
#  points at the newest one) and leaves machine-readable timings beside it:
#
#    dev/logs/test-<mode>-<stamp>-files.csv   per test file: start/end wall
#                                             clock, duration, counts, lobstr
#                                             heap delta
#    dev/logs/test-<mode>-<stamp>-tests.csv   per test: start wall clock,
#                                             duration (serial only)
#    dev/logs/test-<mode>-<stamp>-summary.md  started/finished/wall, top-N
#                                             slowest files and tests
#
#  The log carries an elapsed stamp ([+MM:SS] since the run began) and a wall
#  clock on every per-file line, so a 17-minute `--full` is readable while it
#  runs and correlatable with anything else that logged a timestamp.
#
#  These are gitignored on purpose: timings are machine-specific, so they must
#  never enter docs/artifacts.yml (whose hashes gate the committed figures).
#  `--profile` additionally writes a summaryRprof() profile (costs ~20%).
#
#  Timings are real in `--full` (serial).  In parallel modes testthat replays a
#  file's events when the file finishes, so per-test/per-file wall times are
#  written as NA instead of a silent zero.
#
#  Tests never install.  Install once with `Rscript dev/install-deps.R`
#  (`--with-suggests` for the full external checks).
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

if (!nzchar(system.file(package = "tidymarkers"))) {
  stop("tidymarkers is not installed; run Rscript dev/install-deps.R first")
}

# The cross-scheme contracts.  Matching on file names keeps this independent of
# test content: a renamed contract file is a deliberate edit here, not a silent
# drop from the routine gate.
shared <- c("criteria", "decision", "kasp-checks", "plot-methods", "technology")

args <- commandArgs(trailingOnly = TRUE)
full <- "--full" %in% args
args <- setdiff(args, "--full")
profile <- "--profile" %in% args
args <- setdiff(args, "--profile")

scheme <- NULL
if ("--scheme" %in% args) {
  i <- match("--scheme", args)
  if (i == length(args)) {
    stop("--scheme needs a regex, e.g. --scheme crispr", call. = FALSE)
  }
  scheme <- args[[i + 1L]]
  args <- args[-c(i, i + 1L)]
}

filter <- if (length(args)) args[[1L]] else NULL
if (length(args) > 1L) {
  stop("unexpected arguments: ", paste(args[-1L], collapse = " "), call. = FALSE)
}
if (full && (!is.null(scheme) || !is.null(filter))) {
  stop("--full runs every test; drop --scheme/<filter>", call. = FALSE)
}

files <- list.files("tests/testthat", pattern = "^test-.*[.]R$")
# testthat's filter is matched against the context name (file stem without the
# leading "test-"), so validation must use exactly the same string
contexts <- sub("[.]R$", "", sub("^test[-_]", "", files))
slug <- function(x) gsub("^-+|-+$", "", gsub("[^A-Za-z0-9]+", "-", x))

mode <- if (full) {
  "full"
} else if (!is.null(scheme)) {
  paste0("scheme-", slug(scheme))
} else if (is.null(filter)) {
  "lean"
} else {
  paste0("filter-", slug(filter))
}

if (full) {
  filter <- NULL
  Sys.setenv(TIDYMARKERS_FULL_TESTS = "1", TESTTHAT_PARALLEL = "FALSE")
} else {
  # an inherited TIDYMARKERS_FULL_TESTS must not silently turn the routine loop
  # into the release suite; --full is the switch
  Sys.unsetenv("TIDYMARKERS_FULL_TESTS")
  if (!nzchar(Sys.getenv("TESTTHAT_PARALLEL"))) {
    Sys.setenv(TESTTHAT_PARALLEL = "TRUE")
  }
  if (!is.null(scheme)) {
    hits <- contexts[grepl(scheme, contexts)]
    # an empty match would silently run the shared contracts and look green
    if (!length(hits)) {
      stop("--scheme ", scheme, " matches no test file under tests/testthat/ ",
        "(try crispr, lamp|rpa, indel)",
        call. = FALSE
      )
    }
    filter <- paste(c(shared, scheme), collapse = "|")
  } else if (is.null(filter)) {
    filter <- paste(shared, collapse = "|")
  } else if (!any(grepl(filter, contexts))) {
    # a legacy filter that matches nothing is a silent no-op, not a pass
    stop("filter ", filter, " matches no test file under tests/testthat/", call. = FALSE)
  }
}

selected <- if (is.null(filter)) files else files[grepl(filter, contexts)]
# testthat only goes parallel for more than one file; a single-file filter runs
# serial, so its timings are real even though TESTTHAT_PARALLEL is set
parallel <- identical(Sys.getenv("TESTTHAT_PARALLEL"), "TRUE") &&
  length(selected) > 1L

# ---------------------------------------------------------------------------
#  Logs and timings
# ---------------------------------------------------------------------------
log_dir <- file.path(root, "dev", "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
base <- file.path(log_dir, paste0("test-", mode, "-", stamp))
log_path <- paste0(base, ".log")
files_csv <- paste0(base, "-files.csv")
tests_csv <- paste0(base, "-tests.csv")
summary_path <- paste0(base, "-summary.md")
profile_path <- paste0(base, "-profile.txt")
profile_raw <- paste0(base, "-profile.out")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)

say <- function(...) cat(..., "\n", sep = "")

now <- function() unname(proc.time()[["elapsed"]])
t0 <- now()
wall_stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
# elapsed since the run began, so log lines are ordered without reading the
# wall clock of the machine (which moves, and is local time)
hms_elapsed <- function() {
  s <- as.integer(now() - t0)
  h <- s %/% 3600L
  m <- (s %% 3600L) %/% 60L
  if (h > 0L) {
    sprintf("+%d:%02d:%02d", h, m, s %% 60L)
  } else {
    sprintf("+%02d:%02d", m, s %% 60L)
  }
}

# lobstr::mem_used() forces a GC, so it is sampled once per file boundary, not
# around every test.  A missing lobstr only costs the memory column.
have_lobstr <- requireNamespace("lobstr", quietly = TRUE)
mem_used_mb <- function() {
  if (!have_lobstr) {
    return(NA_real_)
  }
  as.numeric(lobstr::mem_used()) / 1024^2
}

state <- new.env(parent = emptyenv())
state$files <- list()
state$file_order <- character()
state$tests <- list()
state$test_starts <- list()
state$mem_prev <- if (have_lobstr && !parallel) mem_used_mb() else NA_real_
state$started_at <- wall_stamp()
state$done <- FALSE
state$completed <- FALSE

# Same console look as before (ProgressReporter), plus timing bookkeeping.
TimingReporter <- R6::R6Class("TimingReporter",
  inherit = testthat::ProgressReporter,
  public = list(
    start_file = function(filename) {
      key <- basename(filename)
      if (is.null(state$files[[key]])) {
        state$files[[key]] <- list(
          file = key, started_at = wall_stamp(), started = now(),
          duration_s = NA_real_, ended_at = NA_character_, pass = 0L,
          fail = 0L, warn = 0L, skip = 0L, error = 0L, mem_delta_mb = NA_real_
        )
        state$file_order <- c(state$file_order, key)
        # in parallel the events replay when the file is already finished, so a
        # start line there would be a lie; only serial runs announce starts
        if (!parallel) {
          say(sprintf("[%s] %s start", hms_elapsed(), key))
        }
      }
      super$start_file(filename)
    },
    add_result = function(context, test, result) {
      key <- basename(self$file_name)
      f <- state$files[[key]]
      if (!is.null(f)) {
        cls <- class(result)
        if ("expectation_skip" %in% cls) {
          f$skip <- f$skip + 1L
        } else if ("expectation_error" %in% cls) {
          f$error <- f$error + 1L
        } else if ("expectation_failure" %in% cls) {
          f$fail <- f$fail + 1L
        } else if ("expectation_warning" %in% cls) {
          f$warn <- f$warn + 1L
        } else if ("expectation_success" %in% cls) {
          f$pass <- f$pass + 1L
        }
        state$files[[key]] <- f
      }
      super$add_result(context, test, result)
    },
    start_test = function(context, test) {
      # parallel replays call start_test before every event; keep the first time
      key <- paste(self$file_name, test, sep = "\r")
      if (is.null(state$test_starts[[key]])) {
        state$test_starts[[key]] <- list(at = now(), stamp = wall_stamp())
      }
      super$start_test(context, test)
    },
    end_test = function(context, test) {
      key <- paste(self$file_name, test, sep = "\r")
      started <- state$test_starts[[key]]
      state$test_starts[[key]] <- NULL
      state$tests[[length(state$tests) + 1L]] <- list(
        file = basename(self$file_name),
        test = test,
        started_at = if (parallel || is.null(started)) NA_character_ else started$stamp,
        duration_s = if (parallel || is.null(started)) NA_real_ else now() - started$at
      )
    },
    end_file = function() {
      key <- basename(self$file_name)
      f <- state$files[[key]]
      if (!is.null(f)) {
        # parallel replays a whole file's events at once, so the wall time here
        # would be the replay burst, not the run: the duration stays NA.  The
        # stamp is still honest -- this fires as the file completes.
        if (!parallel) {
          f$duration_s <- now() - f$started
        } else {
          # the replay moment is this file's completion, not its start
          f$started_at <- NA_character_
        }
        f$ended_at <- wall_stamp()
        if (isTRUE(have_lobstr) && !parallel && !is.na(state$mem_prev)) {
          m <- mem_used_mb()
          f$mem_delta_mb <- m - state$mem_prev
          state$mem_prev <- m
        }
        state$files[[key]] <- f
        dur <- if (is.na(f$duration_s)) "done (parallel)" else fmt_seconds(f$duration_s)
        mem <- if (is.na(f$mem_delta_mb)) "" else sprintf("  mem %+.1f MB", f$mem_delta_mb)
        say(sprintf(
          "[%s] %s %s%s  | pass %d / fail %d / error %d / warn %d / skip %d",
          hms_elapsed(), f$file, dur, mem, f$pass, f$fail, f$error, f$warn, f$skip
        ))
      }
      write_outputs()
    }
  )
)

# Rewritten after every file so a Ctrl-C still leaves the completed rows behind.
write_outputs <- function() {
  if (!length(state$file_order)) {
    return(invisible(NULL))
  }
  ok <- tryCatch({
    files_df <- do.call(rbind, lapply(state$file_order, function(k) {
      as.data.frame(state$files[[k]], stringsAsFactors = FALSE)
    }))
    utils::write.csv(files_df, files_csv, row.names = FALSE)
    if (length(state$tests)) {
      tests_df <- do.call(rbind, lapply(state$tests, as.data.frame, stringsAsFactors = FALSE))
      tests_df <- tests_df[order(-tests_df$duration_s, na.last = TRUE), , drop = FALSE]
      utils::write.csv(tests_df, tests_csv, row.names = FALSE)
    }
    TRUE
  }, error = function(e) {
    if (!isTRUE(state$write_error)) {
      state$write_error <- TRUE
      say("test.R: WARNING: could not write timings: ", conditionMessage(e))
    }
    FALSE
  })
  invisible(ok)
}

fmt_seconds <- function(s) {
  if (is.na(s)) {
    return("NA")
  }
  if (s < 60) sprintf("%.1fs", s) else sprintf("%dm %02ds", floor(s / 60), round(s) %% 60)
}

summary_lines <- function(status) {
  rows <- lapply(state$file_order, function(k) state$files[[k]])
  total <- function(field) sum(vapply(rows, function(r) as.numeric(r[[field]]), 0))
  slow_files <- rows[order(-vapply(rows, function(r) {
    if (is.na(r$duration_s)) -Inf else r$duration_s
  }, 0))]
  slow_files <- utils::head(slow_files, 8L)

  test_rows <- state$tests
  if (length(test_rows)) {
    test_rows <- test_rows[order(-vapply(test_rows, function(r) {
      if (is.na(r$duration_s)) -Inf else r$duration_s
    }, 0))]
    test_rows <- utils::head(test_rows, 10L)
  }

  out <- c(
    sprintf("# test.R -- %s", mode),
    "",
    sprintf("- status: %s", status),
    sprintf("- filter: %s", if (is.null(filter)) "(all)" else filter),
    sprintf("- started: %s", state$started_at),
    sprintf("- finished: %s", wall_stamp()),
    sprintf("- wall: %s", fmt_seconds(now() - t0)),
    sprintf("- files: %d", length(rows)),
    sprintf(
      "- pass %d / fail %d / error %d / warn %d / skip %d",
      total("pass"), total("fail"), total("error"), total("warn"), total("skip")
    ),
    sprintf(
      "- per-test timings: %s",
      if (parallel) "NA (parallel replay; use --full or a single-file filter)" else "measured"
    ),
    ""
  )
  if (length(slow_files)) {
    out <- c(out, "## Slowest files", "", "| file | started | duration | pass | fail | mem delta (MB) |", "|---|---|---|---|---|---|")
    out <- c(out, vapply(slow_files, function(r) {
      sprintf(
        "| %s | %s | %s | %d | %d | %s |", r$file, r$started_at, fmt_seconds(r$duration_s),
        r$pass, r$fail, if (is.na(r$mem_delta_mb)) "NA" else sprintf("%+.1f", r$mem_delta_mb)
      )
    }, character(1)))
    out <- c(out, "")
  }
  if (length(test_rows)) {
    out <- c(out, "## Slowest tests", "", "| file | test | started | duration |", "|---|---|---|---|")
    out <- c(out, vapply(test_rows, function(r) {
      start_at <- if (is.null(r$started_at)) "NA" else r$started_at
      sprintf(
        "| %s | %s | %s | %s |", r$file, gsub("|", "\\|", r$test, fixed = TRUE),
        start_at, fmt_seconds(r$duration_s)
      )
    }, character(1)))
    out <- c(out, "")
  }
  out <- c(out, sprintf("- log: `%s`", file.path("dev", "logs", basename(log_path))),
    sprintf("- timings: `%s`", file.path("dev", "logs", basename(files_csv))))
  out
}

finish <- function() {
  if (isTRUE(state$done)) {
    return(invisible(NULL))
  }
  state$done <- TRUE
  if (isTRUE(state$profiling)) {
    Rprof(NULL)
    if (file.exists(profile_raw) && file.info(profile_raw)$size > 0) {
      prof <- try(summaryRprof(profile_raw), silent = TRUE)
      if (!inherits(prof, "try-error")) {
        writeLines(capture.output(prof), profile_path)
      }
    }
  }
  write_outputs()
  status <- if (any(vapply(state$files, function(r) r$fail + r$error > 0, logical(1)))) {
    "failed"
  } else if (isTRUE(state$completed)) {
    "ok"
  } else {
    "interrupted"
  }
  lines <- summary_lines(status)
  writeLines(lines, summary_path)
  cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
  for (ext in c(".log", "-summary.md", "-files.csv", "-tests.csv", "-profile.txt")) {
    src <- paste0(base, ext)
    if (file.exists(src)) {
      file.copy(src, file.path(log_dir, paste0("latest", ext)), overwrite = TRUE)
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
#  Run
# ---------------------------------------------------------------------------
say("test.R: ", mode, " -- started ", state$started_at)
say("test.R: filter ", if (is.null(filter)) "(all)" else filter,
  " | parallel ", parallel, " | full tests ", full)
say("test.R: log ", file.path("dev", "logs", basename(log_path)))
say("")

if (isTRUE(profile)) {
  if (parallel) {
    say("test.R: --profile is unavailable in parallel mode; skipping Rprof")
    state$profiling <- FALSE
  } else {
    state$profiling <- TRUE
    Rprof(profile_raw, interval = 0.01)
    say("test.R: Rprof on (interval 0.01 s, ~20% slower)")
  }
}

# The cleanup lives in a function: on.exit() at the top level of a script is
# accepted but never runs, so errors and Ctrl-C would lose the log and timings.
run_tests <- function() {
  on.exit({
    try(finish(), silent = TRUE)
    try(sink(), silent = TRUE)
    try(close(log_con), silent = TRUE)
  })
  testthat::test_local(
    filter = filter,
    reporter = TimingReporter$new(),
    stop_on_failure = TRUE
  )
  state$completed <- TRUE
}

started <- now()
run_tests()
elapsed <- now() - started
say("")
say("test.R: finished in ", fmt_seconds(elapsed), " at ", wall_stamp())
