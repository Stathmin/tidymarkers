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
#  Passes a built-in self-test on every run: the checker is run against
#  synthetic guarded and unguarded bodies, so a scanner that stopped matching
#  is caught.
#
#  Run:  Rscript dev/check-writers.R
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

# The primitives that create a file.  The call form is not part of the rule: a
# package that imports `write_csv()` may call it bare, so the namespace prefix
# is optional (ADR 0035, second hole).  The rule applies to a *call*; a body
# that only builds the text `writeLines(` writes nothing (ADR 0036).
write_primitive <- c(
  "write_csv", "write_lines", "write_file", "writeLines", "saveRDS", "save",
  "ggsave", "write_xml", "download.file", "file_copy", "file_move", "file.copy",
  "file.rename", "writeBin", "write.csv", "write.table", "writeXStringSet",
  "untar", "unzip", "tar", "zip"
)
# A guard refuses a name; it does not merely report one. `cache_present()` is a
# hit test (ZV44), so a body that asks it and then writes is still undeclared
# (ADR 0035).
guard_primitive <- c("refuse_links", "declared_path", "declared_index")

# The called name without a namespace prefix, or NA for a computed head.
call_name <- function(call) {
  if (!is.call(call)) {
    return(NA_character_)
  }
  head <- tryCatch(call[[1L]], error = function(e) NULL)
  if (is.name(head)) {
    return(as.character(head))
  }
  if (is.call(head)) {
    op <- tryCatch(as.character(head[[1L]]), error = function(e) "")
    if (op %in% c("::", ":::")) {
      return(as.character(tryCatch(head[[3L]], error = function(e) "")))
    }
  }
  NA_character_
}

# Every function definition in a parsed file, named by its assignment path:
# `outer` for a top-level one, `outer/inner` for a closure inside it (ADR 0036).
# A closure passed as an argument is `<anonymous>` under its enclosing unit.
collect_units <- function(expr, prefix, out) {
  if (!is.call(expr)) {
    return(out)
  }
  head <- call_name(expr)
  if (identical(head, "function")) {
    unit <- if (nzchar(prefix)) prefix else "<anonymous>"
    while (!is.null(out[[unit]])) unit <- paste0(unit, "#", length(out) + 1L)
    out[[unit]] <- expr[[3L]]
    return(collect_units(expr[[3L]], paste0(unit, "/"), out))
  }
  if (head %in% c("<-", "=") && length(expr) == 3L && identical(call_name(expr[[3L]]), "function")) {
    name <- if (is.name(expr[[2L]])) as.character(expr[[2L]]) else "<anonymous>"
    return(collect_units(expr[[3L]], paste0(prefix, name), out))
  }
  for (i in seq_along(expr)) out <- collect_units(expr[[i]], prefix, out)
  out
}

# Write and guard calls in one unit, in body order, not descending into a
# nested definition: that definition is its own unit and needs its own guard.
calls_in <- function(expr) {
  found <- character()
  walk <- function(e) {
    if (!is.call(e)) {
      return(invisible(NULL))
    }
    name <- call_name(e)
    if (!is.na(name) && name %in% c(write_primitive, guard_primitive)) {
      found <<- c(found, name)
    }
    if (identical(name, "function")) {
      return(invisible(NULL))
    }
    for (i in seq_along(e)) walk(e[[i]])
    invisible(NULL)
  }
  walk(expr)
  found
}

# Function units in R/, from the parse tree so a call in a comment or a string
# does not count.
r_functions <- function(files = list.files("R", "\\.R$", full.names = TRUE)) {
  out <- list()
  for (file in files) {
    exprs <- tryCatch(parse(file), error = function(e) NULL)
    if (is.null(exprs)) next
    for (expr in as.list(exprs)) out <- collect_units(expr, "", out)
  }
  out
}

# Every write call in a unit, and whether a guard call comes before the first
# one.  Position is body order, not execution order: a guard in an untaken
# branch still counts, which is why this stays a heuristic (ADR 0037).
scan_writers <- function(bodies) {
  out <- list()
  for (symbol in names(bodies)) {
    calls <- calls_in(bodies[[symbol]])
    writes <- intersect(calls, write_primitive)
    if (!length(writes)) next
    guards <- which(calls %in% guard_primitive)
    out[[symbol]] <- list(
      writes = sort(unique(writes)),
      guarded = length(guards) > 0L && min(guards) < min(which(calls %in% write_primitive))
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
    imported = "function(path) write_csv(x, path)",
    cached = "function(path) { if (cache_present(path)) return(path); writeLines(\"x\", path) }",
    declared = "function(dir) { refuse_links(dir, \"a.csv\"); readr::write_csv(x, dir) }",
    laundered = "function(d) { f <- function(p) writeLines(\"x\", p); refuse_links(d, \"a.csv\"); f(d) }",
    nested_declared = "function(d) { f <- function(p) { refuse_links(p, \"a.csv\"); writeLines(\"x\", p) }; f(d) }",
    strings = "function(x) paste0(\"writeLines(\", x, \")\")",
    sequences = "function(p) Biostrings::writeXStringSet(x, p)",
    archive = "function(d, a) utils::untar(a, exdir = d)",
    zipped = "function(d, a) utils::unzip(a, exdir = d)",
    late = "function(p) { writeLines(\"x\", p); refuse_links(p, \"a.csv\") }"
  )
  exprs <- Map(
    function(name, text) parse(text = paste0(name, " <- ", text))[[1L]],
    names(bodies), bodies
  )
  units <- list()
  for (expr in exprs) units <- collect_units(expr, "", units)
  hits <- scan_writers(units)
  empty <- check_writers(list(), hits)
  expected <- paste0(
    "undeclared writer (guard it, or add to dev/writers.yml): ",
    c("archive", "cached", "imported", "late", "laundered/f", "sequences", "undeclared", "zipped")
  )
  if (!identical(sort(empty), sort(expected))) {
    stop("self-test: an unguarded writer was not caught")
  }
  if ("strings" %in% names(hits)) {
    stop("self-test: a call built as text was reported as a writer")
  }
  ok <- check_writers(
    list(
      list(symbol = "undeclared", reason = "self-test"),
      list(symbol = "imported", reason = "self-test"),
      list(symbol = "cached", reason = "self-test"),
      list(symbol = "laundered/f", reason = "self-test"),
      list(symbol = "sequences", reason = "self-test"),
      list(symbol = "archive", reason = "self-test"),
      list(symbol = "zipped", reason = "self-test"),
      list(symbol = "late", reason = "self-test")
    ),
    hits
  )
  if (length(ok)) {
    stop("self-test: a guarded body or a registered writer was flagged")
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

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
