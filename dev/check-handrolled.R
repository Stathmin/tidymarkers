#!/usr/bin/env Rscript
# =============================================================================
#  dev/check-handrolled.R  --  the reuse ledger: register or replace
# =============================================================================
#  Adapted from GWASILIY (R/handrolled.{R,yml} + scripts/handrolled.R).  The
#  rule there and here: a function whose body reimplements something an already
#  installed package does must be *registered* with a decision, or the gate
#  fails.  Read docs/capabilities.md (GWASILIY) / the headers of R/marker_*.R
#  first; this package's whole design rule is "delegate to Biostrings,
#  TmCalculator and DECIPHER, never reimplement the primitive".
#
#  status:
#    replaced  - the hand-roll is gone; `replacement` is what replaced it
#    justified - the symbol stays, for the stated reason
#    todo      - the symbol must be replaced; the gate fails while it exists
#
#  Fails when:
#    - a heuristic hit has no registry entry              (unregistered hand-roll)
#    - a `todo` symbol still exists in R/                 (unfinished replacement)
#    - a `replaced` symbol still exists in R/             (stale registry claim)
#    - an entry has an unknown status or no reason
#    - a `package::function` in `replacement` does not exist
#
#  Run:  Rscript dev/check-handrolled.R      (exit 1 on a problem)
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

# The heuristic is deliberately narrow: these are the sequence primitives the
# package already gets from Biostrings / TmCalculator.  Ordinary `for` loops and
# `sum()` calls are not listed, because the resulting noise would make the gate
# useless.  A pattern that matches nothing today is a guard for the next change.
handrolled_patterns <- function() {
  c(
    manual_complement = "chartr\\(\\s*\"[ACGTacgtRYSWKMBDHV]+\"",
    manual_revcomp = "paste0?\\(\\s*rev\\(",
    manual_gc = "sum\\([^()]*==\\s*\"[GC]\"",
    manual_tm = "81\\.5|64\\.9"
  )
}

# Registered hand-rolls (the ledger). The decisions live in dev/handrolled.yml
# (data, not code) so dev/graphify-handrolled.R can emit the `supersedes` edges.
handrolled_registry <- function() {
  yaml::read_yaml("dev/handrolled.yml")$entries
}

# Top-level `name <- function(...)` definitions in R/.  Parsed, not grepped, so a
# helper named in a comment or a string does not trigger the gate.
r_functions <- function(files = list.files("R", "\\.R$", full.names = TRUE)) {
  out <- list()
  for (file in files) {
    exprs <- tryCatch(parse(file), error = function(e) NULL)
    if (is.null(exprs)) next
    for (expr in as.list(exprs)) {
      if (!is.call(expr) || length(expr) < 3L) next
      if (!is.name(expr[[2]]) || !is.call(expr[[3]])) next
      if (!identical(expr[[3]][[1]], as.name("function"))) next
      out[[as.character(expr[[2]])]] <- paste(deparse(expr[[3]]), collapse = "\n")
    }
  }
  out
}

scan_handrolled <- function(bodies) {
  patterns <- handrolled_patterns()
  hits <- list()
  for (symbol in names(bodies)) {
    found <- names(patterns)[
      vapply(patterns, function(p) grepl(p, bodies[[symbol]], perl = TRUE), logical(1))
    ]
    if (length(found)) {
      hits[[symbol]] <- found
    }
  }
  hits
}

# `pkg::fn` exists in an installed package.
replacement_exists <- function(replacement) {
  parts <- strsplit(replacement, "::", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) {
    return(FALSE)
  }
  requireNamespace(parts[[1L]], quietly = TRUE) &&
    exists(parts[[2L]], envir = asNamespace(parts[[1L]]), inherits = FALSE)
}

check_handrolled <- function(registry, hits, bodies) {
  problems <- character()
  symbols <- vapply(registry, `[[`, character(1), "symbol")
  status <- vapply(registry, `[[`, character(1), "status")
  reason <- vapply(registry, function(e) e$reason %||% "", character(1))
  replacement <- vapply(registry, function(e) e$replacement %||% "", character(1))

  unknown <- setdiff(unique(status), c("replaced", "justified", "todo"))
  if (length(unknown)) {
    problems <- c(problems, paste0("unknown hand-roll status: ", unknown))
  }
  no_reason <- symbols[!nzchar(reason)]
  if (length(no_reason)) {
    problems <- c(problems, paste0("hand-roll entry without a reason: ", no_reason))
  }
  still_there <- intersect(symbols[status == "todo"], names(bodies))
  if (length(still_there)) {
    problems <- c(problems, paste0("hand-roll still present (todo): ", still_there))
  }
  gone <- symbols[status == "replaced"]
  stale <- intersect(gone, names(bodies))
  if (length(stale)) {
    problems <- c(problems, paste0("registry says replaced but symbol exists: ", stale))
  }
  unregistered <- setdiff(names(hits), symbols)
  if (length(unregistered)) {
    problems <- c(
      problems,
      paste0("unregistered hand-roll (add to dev/handrolled.yml): ", unregistered)
    )
  }
  named <- replacement[nzchar(replacement)]
  missing <- named[!vapply(named, replacement_exists, logical(1))]
  if (length(missing)) {
    problems <- c(problems, paste0("replacement does not exist: ", missing))
  }
  problems
}

`%||%` <- function(x, y) if (is.null(x)) y else x

bodies <- r_functions()
hits <- scan_handrolled(bodies)
registry <- handrolled_registry()
problems <- check_handrolled(registry, hits, bodies)

if (length(problems)) {
  cat("hand-roll gate failed:\n")
  cat(paste0("  - ", problems, "\n"), sep = "")
  quit(status = 1L)
}
cat(
  "  ok   ", length(registry), " registered hand-roll(s); ",
  length(hits), " heuristic hit(s) all registered\n",
  sep = ""
)
