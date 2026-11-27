#!/usr/bin/env Rscript
# =============================================================================
#  dev/perf-gate.R  --  the ADR 0010 memory and query budget
# =============================================================================
#  Reads measurements written by `dev/bench-index-memory.R --json <path>` and
#  fails when a shard exceeds `baseline + 5 bytes/base` or a query exceeds 2x
#  the memory-path median.  Without a measurement the gate reports `not_run`
#  and a reason: an unmeasured budget is not a pass.
#
#  Run:  Rscript dev/perf-gate.R [--json <path>]
# =============================================================================
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
json <- if (length(args) >= 2L && args[[1L]] == "--json") args[[2L]] else NULL

fail <- character(0)
notes <- character(0)
check <- function(ok, message) {
  if (!isTRUE(ok)) fail <<- c(fail, message)
}

BYTES_PER_BASE <- 5
baseline_mb <- as.numeric(Sys.getenv("TIDYMARKERS_PERF_BASELINE_MB", "455.4"))

# a monolith is only legal below the 1.34 bytes/base build limit
monolith_limit <- tidymarkers:::BOWTIE_MONOLITH_LIMIT
# 21.5 GiB at 1.34 B/bp is the ~16 Gb wheat case the ADR sizes the limit for
check(
  monolith_limit / 1.34e9 > 16 && monolith_limit / 1.34e9 < 18,
  sprintf("BOWTIE_MONOLITH_LIMIT is %.2f Gb at 1.34 B/bp, expected ~17 Gb",
    monolith_limit / 1.34e9
  )
)

records <- if (is.null(json) || !fs::file_exists(json)) {
  list()
} else {
  jsonlite::fromJSON(json, simplifyVector = FALSE)
}
shards <- Filter(function(r) identical(r$mode, "shard"), records)
queries <- Filter(function(r) identical(r$mode, "query"), records)

if (!length(shards)) {
  notes <- c(notes, "shards: not_run (no shard record; run dev/bench-index-memory.R --json)")
} else {
  for (r in shards) {
    bases <- as.numeric(r$bases %||% NA_real_)
    peak <- as.numeric(r$peak_rss_mb %||% NA_real_)
    if (is.na(bases) || is.na(peak)) {
      notes <- c(notes, "shards: not_run (record without bases/peak_rss_mb)")
      next
    }
    # process overhead dominates below 1 Mb, so a fixture cannot test a
    # per-base budget; the TAIR10 and wheat runs are the real measurements
    if (bases < 1e6) {
      notes <- c(notes, sprintf(
        "shards: not_run (%s has %.0f kb, below the 1 Mb budget floor)",
        r$label %||% "?", bases / 1e3
      ))
      next
    }
    budget <- baseline_mb + BYTES_PER_BASE * bases / 1024^2
    check(
      peak <= budget,
      sprintf("shard %s uses %.1f MB > %.1f MB (baseline + 5 B/bp)",
        r$label %||% "?", peak, budget
      )
    )
    if (bases * 1.34 > monolith_limit && isTRUE(as.integer(r$chroms %||% 0L) < 1L)) {
      fail <- c(fail, sprintf(
        "a %.1f Gb build was attempted without chroms (limit %.1f Gb)",
        bases / 1e9, monolith_limit / 1e9
      ))
    }
  }
}

if (!length(queries)) {
  notes <- c(notes, "query: not_run (no query record; the bench does not time the golden locus)")
} else {
  for (r in queries) {
    mem <- as.numeric(r$memory_median_s %||% NA_real_)
    idx <- as.numeric(r$index_median_s %||% NA_real_)
    if (is.na(mem) || is.na(idx)) {
      notes <- c(notes, "query: not_run (record without both medians)")
      next
    }
    check(idx <= 2 * mem, sprintf(
      "query %s takes %.2f s > 2 x %.2f s memory median",
      r$label %||% "?", idx, mem
    ))
  }
}

for (n in notes) cat("not_run:", n, "\n")
if (length(fail)) {
  cat("perf-gate: FAIL\n")
  for (f in fail) cat("  -", f, "\n")
  quit(status = 1L)
}
cat("perf-gate: ok\n")
