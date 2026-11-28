#!/usr/bin/env Rscript
# =============================================================================
#  dev/check-writers.R  --  every writer declares the names it creates
# =============================================================================
#  ADR 0029 and its amendments (0030-0034) established one rule for every
#  function that creates a file: the names it creates are declared before the
#  first write, so a link at one of them is refused instead of followed.  The
#  rule was closed by hand, one writer at a time (ZV29/ZV33-ZV44); this gate
#  makes the roster a checked ledger, not a memory.
#
#  A top-level function in R/ that calls a write primitive must also call a
#  guard in the same body, or be registered in dev/writers.yml with a reason.
#  The heuristic is deliberately narrow: only the primitives that actually
#  create a file are listed, and `dev/writers.yml` records those whose name is
#  declared by their caller (a helper that cannot know the directory).
#
#  Fails when:
#    - an unguarded write has no ledger entry        (new undeclared writer)
#    - a ledger entry no longer has a write          (stale justification)
#    - an entry has no reason
#  Passes `--self-test`: the checker is run against one unguarded and one
#  guarded synthetic body, so a scanner that stopped matching is caught.
#
#  Run:  Rscript dev/check-writers.R [--self-test]
# =============================================================================

root <- local({
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
})
setwd(root)

if (!requireNamespace("yaml", quietly = TRUE)) {
  message("skip (yaml not installed)")
  quit(status = 0L)
}

# The primitives that create a file.  `writeLines()` on a string is also a hit
# when the string is a script another function emits; a ledger entry says so.
write_pattern <- paste0(
  "(", paste(c(
    "readr::write_csv", "readr::write_lines", "readr::write_file",
    "writeLines", "saveRDS", "save", "ggplot2::ggsave", "xml2::write_xml",
    "utils::download.file", "fs::file_copy", "fs::file_move", "file.copy",
    "file.rename", "writeBin", "utils::write.csv", "utils::write.table"
  ), collapse = "|"), ")\\s*\\("
)
guard_pattern <- paste0(
  "refuse_links\\(|declared_path\\(|declared_index\\(|cache_present\\("
)

# Top-level `name <- function(...)` bodies in R/.  Parsed, not grepped, so a
# helper named in a comment does not count.
r_functions <- function(files = list.files("R", "\\.R$", full.names = TRUE)) {
  out <- list()
  for (file in files) {
    exprs <- tryCatch(parse(file), error = function(e) NULL)
    if (is.null(exprs)) next
    for (expr in as.list(exprs)) {
      if (!is.call(expr) || length(expr) < 3L) next
      if (!is.name(expr[[2]]) || !is.call(expr[[3]])) next
      if (!identical(expr[[3]][[1]], as.name("function"))) next
      out[[as.character(expr[[2]])]] <- paste(deparse(expr[[3]][[3]]), collapse = "\n")
    }
  }
  out
}

# Every write call, and whether the same body declares a name first.
scan_writers <- function(bodies) {
  out <- list()
  for (symbol in names(bodies)) {
    body <- bodies[[symbol]]
    found <- regmatches(body, gregexpr(write_pattern, body, perl = TRUE))[[1L]]
    if (!length(found)) next
    out[[symbol]] <- list(
      writes = sort(unique(sub("\\s*\\($", "", found))),
      guarded = grepl(guard_pattern, body, perl = TRUE)
    )
  }
  out
}

writer_registry <- function() {
  yaml::read_yaml("dev/writers.yml")$entries
}

check_writers <- function(registry, hits) {
  problems <- character()
  symbols <- vapply(registry, `[[`, character(1), "symbol")
  reason <- vapply(registry, function(e) e$reason %||% "", character(1))
  no_reason <- symbols[!nzchar(reason)]
  if (length(no_reason)) {
    problems <- c(problems, paste0("writer entry without a reason: ", no_reason))
  }
  unguarded <- names(hits)[!vapply(hits, `[[`, logical(1), "guarded")]
  unregistered <- setdiff(unguarded, symbols)
  if (length(unregistered)) {
    problems <- c(
      problems,
      paste0("undeclared writer (guard it, or add to dev/writers.yml): ", unregistered)
    )
  }
  stale <- setdiff(symbols, unguarded)
  if (length(stale)) {
    problems <- c(
      problems,
      paste0("dev/writers.yml is stale (the symbol no longer writes unguarded): ", stale)
    )
  }
  problems
}

self_test <- function() {
  bodies <- c(
    undeclared = "function(path) writeLines(\"x\", path)",
    declared = "function(dir) { refuse_links(dir, \"a.csv\"); readr::write_csv(x, dir) }"
  )
  hits <- scan_writers(bodies)
  empty <- check_writers(list(), hits)
  if (!identical(empty, "undeclared writer (guard it, or add to dev/writers.yml): undeclared")) {
    stop("self-test: an unguarded writer was not caught")
  }
  ok <- check_writers(
    list(list(symbol = "undeclared", reason = "self-test")), hits
  )
  if (length(ok)) {
    stop("self-test: a guarded body or a registered writer was flagged")
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
self_test()
bodies <- r_functions()
hits <- scan_writers(bodies)
registry <- writer_registry()
problems <- check_writers(registry, hits)

if (length(problems)) {
  cat("writer gate failed:\n")
  cat(paste0("  - ", problems, "\n"), sep = "")
  quit(status = 1L)
}
cat(
  "  ok   ", length(hits), " writer(s), ",
  sum(vapply(hits, `[[`, logical(1), "guarded")), " guard in place, ",
  length(registry), " registered with a reason; self-test passed\n",
  sep = ""
)
