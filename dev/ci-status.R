#!/usr/bin/env Rscript
# =============================================================================
#  dev/ci-status.R  --  what is green on which runner, and why not
# =============================================================================
#  Reads GitHub Actions through the authenticated `gh` CLI.  No new R
#  dependency, no CRAN, no token juggling: `gh` already holds the credential.
#
#    Rscript dev/ci-status.R          # last 5 runs
#    Rscript dev/ci-status.R 12       # last 12 runs
#    Rscript dev/ci-status.R 5 --all  # also show passing jobs' wall time
#
#  For each run: per-job conclusion and wall time, then the R CMD check
#  ERROR / WARNING / NOTE lines pulled from the failed-step log.  That is the
#  "spectre of systems" view: one OS failing a check the others pass is a
#  platform issue, every OS failing the same line is a package issue.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) && grepl("^[0-9]+$", args[[1L]])) as.integer(args[[1L]]) else 5L
show_passing <- "--all" %in% args

run_gh <- function(...) {
  a <- c(...)
  # system2() does not quote arguments on Unix: a jq expression containing `|`
  # would become a shell pipe.  Quote every argument.
  if (.Platform$OS.type == "unix") a <- vapply(a, shQuote, character(1))
  out <- suppressWarnings(system2("gh", a, stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status"))) {
    warning(trimws(paste(out, collapse = " ")), call. = FALSE)
    return(character())
  }
  out
}

slug <- run_gh("repo", "view", "--json", "nameWithOwner", "-q", ".nameWithOwner")
if (!length(slug)) stop("gh is not authenticated, or this is not a GitHub repo")
slug <- trimws(slug[[1L]])

tsv <- function(path, sel) {
  lines <- run_gh("api", path, "-q", sel)
  if (!length(lines)) {
    return(NULL)
  }
  utils::read.delim(text = lines, header = FALSE, sep = "\t", quote = "",
                    stringsAsFactors = FALSE, fill = TRUE)
}

secs <- function(from, to) {
  if (!nzchar(from) || !nzchar(to)) {
    return(NA_character_)
  }
  d <- as.numeric(difftime(
    as.POSIXct(to, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    as.POSIXct(from, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    units = "secs"
  ))
  sprintf("%dm%02ds", d %/% 60, d %% 60)
}

# Only the R CMD check diagnostic lines, never the whole log: a failed job can
# be tens of MB, and the question is always "which check line failed here".
job_diagnostics <- function(slug, job_id, max_lines = 8L) {
  log <- tempfile(fileext = ".log")
  on.exit(unlink(log), add = TRUE)
  a <- c("api", sprintf("repos/%s/actions/jobs/%s/logs", slug, job_id))
  if (.Platform$OS.type == "unix") a <- vapply(a, shQuote, character(1))
  st <- suppressWarnings(system2("gh", a, stdout = log, stderr = FALSE))
  if (!identical(st, 0L) || !file.exists(log)) {
    return(character())
  }
  # R CMD check marks a failed check on the `* checking ...` line itself
  # (`... ERROR` / `... WARNING` / `... NOTE`) and repeats it in the status
  # summary.  Setup tools also print WARNING, so match only those two shapes.
  hits <- suppressWarnings(system2("grep",
    c("-E", shQuote("(\\.\\.\\. (ERROR|WARNING|NOTE)$)|(Status: [0-9]+)"), shQuote(log)),
    stdout = TRUE, stderr = FALSE))
  utils::head(hits, max_lines)
}

runs <- tsv(
  sprintf("repos/%s/actions/runs?per_page=%d", slug, n),
  ".workflow_runs[] | [.id, .name, .head_branch, .event, .status, .conclusion, .created_at] | @tsv"
)
if (is.null(runs)) {
  cat("no workflow runs for", slug, "\n")
  quit(status = 0L)
}
# id, name, head_branch, event, status, conclusion, created_at
names(runs) <- c("id", "name", "branch", "event", "status", "conclusion", "created")
cat(sprintf("%s -- last %d run(s)\n\n", slug, nrow(runs)))

for (i in seq_len(nrow(runs))) {
  r <- runs[i, ]
  cat(sprintf("run %s  %s  [%s %s]  %s  %s\n",
    r$id, r$name, r$branch, r$event, r$conclusion, r$created))
  jobs <- tsv(
    sprintf("repos/%s/actions/runs/%s/jobs?per_page=100", slug, r$id),
    ".jobs[] | [.id, .name, .conclusion, .started_at, .completed_at] | @tsv"
  )
  if (!is.null(jobs)) {
    # id, name, conclusion, started_at, completed_at
    names(jobs) <- c("id", "name", "conclusion", "started", "completed")
    for (j in seq_len(nrow(jobs))) {
      ok <- identical(jobs$conclusion[[j]], "success")
      if (!ok || show_passing) {
        cat(sprintf("  %-28s %-10s %s\n", jobs$name[[j]], jobs$conclusion[[j]],
          secs(jobs$started[[j]], jobs$completed[[j]])))
      }
    }
    failed <- jobs[jobs$conclusion %in% c("failure", "cancelled"), , drop = FALSE]
    for (j in seq_len(nrow(failed))) {
      hits <- job_diagnostics(slug, failed$id[[j]])
      if (length(hits)) {
        cat(sprintf("  --- %s: R CMD check diagnostics ---\n", failed$name[[j]]))
        cat(paste0("  ", substr(trimws(hits), 1L, 200L), collapse = "\n"), "\n")
      }
    }
  }
  cat("\n")
}
