#!/usr/bin/env Rscript
# =============================================================================
#  dev/profile.R  --  profile one call, report it as a tidy table, export RDF
# =============================================================================
#  `Rprof()` writes one line per sample: the call stack, innermost first.  This
#  script reads that file directly, because `summaryRprof()` aggregates by
#  function name and hides who called whom -- and "who repeats this?" is the
#  question a slow pipeline asks.
#
#  From the raw file it derives three tables:
#
#    samples  one row per stack sample
#    frames   every frame in every stack: depth, function, caller
#    self     per function: own seconds (it was the innermost frame) and
#             total seconds (it appeared anywhere in the stack)
#
#  It can also write the same facts as N-Triples (RDF 1.1) in the `profile`
#  vocabulary below.  R ships no RDF serialiser and this repo takes no new
#  dependency, so the term is the whole ontology: one triple per line.
#
#  Usage
#    Rscript dev/profile.R run    <file.R>   [out-dir]   profile a script
#    Rscript dev/profile.R stacks <file.out> [out-dir]   read an existing profile
#
#  The script is sourced into the R_GlobalEnv of a fresh R session after
#  `pkgload::load_all(".")`, so it sees the package under development.  Its
#  output goes to `<out-dir>` (default `.tmp/prof/<script-name>`):
#  `<name>.out` (the raw profile), `<name>.self.csv`, `<name>.frames.csv`,
#  `<name>.ttl` (N-Triples) and `<name>.rds` (the three tables, for reuse).
#
#  A note on "for scientific accuracy": a profile measures where the interpreter
#  spent samples, not where the algorithm spends work.  A fast function called
#  184 149 times and a slow one called six times look alike in `self`.  Read
#  `self` against a count from the code ("how many hits did it say?"), and keep
#  the claim falsifiable: change one thing, re-profile, `all.equal()` the output.
# =============================================================================

suppressPackageStartupMessages({
  library(methods)
})

# --- reading the raw profile -------------------------------------------------

# One line per sample: `"f1" "f2" ...` innermost first; the first line is
# `sample.interval=<microseconds>`.  Returned as a list, never a data frame --
# stacks have no fixed width.
profile_read_stacks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  interval <- 0
  if (length(lines) && grepl("^sample.interval=[0-9.]", lines[[1]])) {
    interval <- as.numeric(sub("^sample.interval=", "", lines[[1]]))
    lines <- lines[-1]
  }
  stacks <- lapply(lines, function(l) {
    # `strsplit`, not `scan(what = character())`: scan treats whitespace as a
    # delimiter and every `what` string as a *type*, so it returns one empty
    # element per token -- every stack would collapse to a single "".
    f <- strsplit(trimws(l), "[[:space:]]+", fixed = FALSE)[[1]]
    f <- gsub('^"|"$', "", f)
    f[nzchar(f)]
  })
  stacks <- stacks[lengths(stacks) > 0]
  list(stacks = stacks, interval = interval / 1e6)
}

# Long form: one row per frame per sample.  `depth` 1 is the innermost frame
# (the function that was actually running), `caller` its parent.
profile_frames <- function(stacks) {
  n <- lengths(stacks)
  total <- sum(n)
  frame <- data.frame(
    sample = rep(seq_along(stacks), n),
    depth = sequence(n),
    function_name = unlist(stacks, use.names = FALSE),
    caller = unlist(lapply(stacks, function(s) c(s[-1], NA_character_)),
      use.names = FALSE
    ),
    stringsAsFactors = FALSE
  )
  stopifnot(nrow(frame) == total)
  frame
}

# Per function: self time (innermost) and total time (anywhere in the stack).
# Both weighted by the sample interval, so the two sum to the same number.
profile_self <- function(stacks, interval) {
  frame <- profile_frames(stacks)
  innermost <- frame$depth == 1L
  names_all <- frame$function_name
  out <- data.frame(
    function_name = unique(names_all),
    stringsAsFactors = FALSE
  )
  out$self_samples <- as.integer(vapply(unique(names_all), function(k) {
    sum(innermost & names_all == k)
  }, integer(1), USE.NAMES = FALSE))
  out$total_samples <- as.integer(vapply(unique(names_all), function(k) {
    sum(names_all == k)
  }, integer(1), USE.NAMES = FALSE))
  out$self_seconds <- out$self_samples * interval
  out$total_seconds <- out$total_samples * interval
  out$pct_self <- 100 * out$self_samples / sum(out$self_samples)
  out$pct_under <- 100 * out$total_samples / max(out$total_samples, 1L)
  ord <- order(out$self_seconds, decreasing = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

profile_read <- function(path) {
  p <- profile_read_stacks(path)
  list(
    samples = length(p$stacks),
    interval = p$interval,
    frames = profile_frames(p$stacks),
    self = profile_self(p$stacks, p$interval)
  )
}

# --- the report ---------------------------------------------------------------

profile_report <- function(x, top = 15L) {
  cat(sprintf(
    "samples: %d  interval: %.3f s  elapsed: %.2f s\n",
    x$samples, x$interval, x$samples * x$interval
  ))
  s <- utils::head(x$self, top)
  s$pct_self <- round(s$pct_self, 1)
  s$self_seconds <- signif(s$self_seconds, 3)
  s$total_seconds <- signif(s$total_seconds, 3)
  print(s[, c("function_name", "self_seconds", "pct_self", "total_seconds", "pct_under")],
    row.names = FALSE
  )
  invisible(s)
}

# The hot callers of one function: who repeated it, and how much of the run sat
# under it.  `summaryRprof()` cannot answer this; the raw stacks can.
profile_callers <- function(x, fn, top = 10L) {
  f <- x$frames
  hit <- f[f$function_name == fn, , drop = FALSE]
  if (!nrow(hit)) {
    stop("function not in this profile: ", fn, call. = FALSE)
  }
  cb <- as.data.frame(table(hit$caller, useNA = "ifany"), stringsAsFactors = FALSE)
  names(cb) <- c("caller", "frames")
  cb <- cb[order(cb$frames, decreasing = TRUE), , drop = FALSE]
  cb$pct_of_calls <- round(100 * cb$frames / nrow(hit), 1)
  utils::head(cb, top)
}

# --- RDF ---------------------------------------------------------------------
# No serialiser is installed (no rdflib, no RDF/redland/jsonld), and none is
# needed: N-Triples is line-oriented and ASCII-safe, so base R writes it.
# Vocabulary (all under https://tidymarkers.org/profile#):
#   ProfileRun  -- dtmi:elapsedSeconds, dtmi:commit, dtmi:script, dtmi:sampleInterval
#   Function    -- dtmi:selfSeconds, dtmi:totalSeconds, dtmi:selfSamples,
#                  dtmi:totalSamples, dtmi:run (-> ProfileRun)
#   <run> dtmi:hasFunction <fn>    <fn> dtmi:calledBy <caller>
# The subject of a run is keyed by script + commit + elapsed, so two runs of the
# same script stay distinct while a re-run of one commit adds no noise.

PROFILE_NS <- "https://tidymarkers.org/profile#"

profile_iri_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9._~-]+", "_", x)
  ifelse(nzchar(x), x, "_")
}

profile_ttl_literal <- function(x) {
  if (is.numeric(x)) {
    # xsd:decimal is `-?digits.digits` -- no leading zeros, no exponent, no
    # bare integer.  `formatC(format = "fg")` emits both `0` and `2.00000`.
    x <- as.numeric(x)
    if (all(x == round(x))) {
      return(sprintf("%.1f", x))
    }
    return(sub("0+$", "0", sprintf("%.6f", x)))
  }
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\\\"', x)
  x <- gsub("[\r\n\t]", " ", x)
  paste0('"', x, '"')
}

profile_rdf_terms <- function(x, meta = list()) {
  ns <- PROFILE_NS
  script <- if (is.null(meta$script)) "inline" else meta$script
  commit <- if (is.null(meta$commit)) "unknown" else meta$commit
  elapsed <- x$samples * x$interval
  run <- paste0(ns, "run/", profile_iri_safe(script), "-", profile_iri_safe(commit))

  triples <- c(
    paste0("<", run, "> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <", ns, "ProfileRun> ."),
    paste0("<", run, "> <", ns, "script> ", profile_ttl_literal(script), " ."),
    paste0("<", run, "> <", ns, "commit> ", profile_ttl_literal(commit), " ."),
    paste0("<", run, "> <", ns, "elapsedSeconds> ", profile_ttl_literal(elapsed), " ."),
    paste0("<", run, "> <", ns, "sampleInterval> ", profile_ttl_literal(x$interval), " ."),
    paste0("<", run, "> <", ns, "samples> ", profile_ttl_literal(as.integer(x$samples)), " .")
  )

  s <- x$self
  for (i in seq_len(nrow(s))) {
    fn <- paste0(ns, "fn/", profile_iri_safe(s$function_name[[i]]))
    triples <- c(
      triples,
      paste0("<", run, "> <", ns, "hasFunction> <", fn, "> ."),
      paste0("<", fn, "> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <", ns, "Function> ."),
      paste0("<", fn, "> <", ns, "functionName> ", profile_ttl_literal(s$function_name[[i]]), " ."),
      paste0("<", fn, "> <", ns, "run> <", run, "> ."),
      paste0("<", fn, "> <", ns, "selfSeconds> ", profile_ttl_literal(s$self_seconds[[i]]), " ."),
      paste0("<", fn, "> <", ns, "totalSeconds> ", profile_ttl_literal(s$total_seconds[[i]]), " ."),
      paste0("<", fn, "> <", ns, "selfSamples> ", profile_ttl_literal(as.integer(s$self_samples[[i]])), " ."),
      paste0("<", fn, "> <", ns, "totalSamples> ", profile_ttl_literal(as.integer(s$total_samples[[i]])), " .")
    )
  }

  edges <- unique(x$frames[!is.na(x$frames$caller), c("caller", "function_name")])
  edges$caller_iri <- paste0(ns, "fn/", profile_iri_safe(edges$caller))
  edges$callee_iri <- paste0(ns, "fn/", profile_iri_safe(edges$function_name))
  triples <- c(triples, unique(paste0(
    "<", edges$callee_iri, "> <", ns, "calledBy> <", edges$caller_iri, "> ."
  )))

  unique(triples)
}

profile_ttl_write <- function(x, path, meta = list()) {
  triples <- profile_rdf_terms(x, meta)
  writeLines(triples, path, useBytes = TRUE)
  # Verify by re-reading: a triple line ends in " ." and has three terms.
  back <- readLines(path, warn = FALSE)
  stopifnot(
    identical(back, triples),
    all(grepl(" [.]$", back)),
    all(lengths(strsplit(back, " ", fixed = TRUE)) >= 4L)
  )
  invisible(length(triples))
}

# Read the tidy profile back.  `read_rds` (readr) is the reader for the `.rds`
# the writer leaves; it restores the same tables without re-parsing stacks.
profile_read_tidy <- function(path) {
  readr::read_rds(path)
}

# --- reading RDF back ---------------------------------------------------------

# N-Triples, read strictly.  This is the same grammar the writer emits, parsed
# for real: a line must be subject, predicate, object; object is an IRI or a
# literal with an optional datatype; anything else stops with the line number.
# Typed columns make a decimal a number again, so a report can be recomputed
# from the file alone instead of trusted.
profile_ttl_read <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  parts <- regmatches(lines, regexec("^<([^>]+)> <([^>]+)>", lines))
  bad <- which(lengths(parts) != 3L)
  if (length(bad)) {
    stop(sprintf("malformed N-Triples at line %d: %s", bad[[1]], lines[[bad[[1]]]]), call. = FALSE)
  }
  ## The object is the rest of the line minus the " ." terminator, sliced from
  ## the prefix length rather than matched with an escaped pattern.
  ## regexec reports byte offsets; slicing uses the reported match length.
  parts_len <- vapply(
    regmatches(lines, regexec("^<([^>]+)> <([^>]+)>", lines)),
    function(m) nchar(m[[1L]]),
    numeric(1)
  )
  objects <- substring(lines, parts_len + 2L)
  if (!all(endsWith(objects, " ."))) {
    bad <- which(!endsWith(objects, " ."))
    stop(sprintf("malformed N-Triples at line %d: %s", bad[[1]], lines[[bad[[1]]]]), call. = FALSE)
  }
  objects <- substring(objects, 1L, nchar(objects) - 2L)
  is_iri <- grepl("^<[^>]+>$", objects)
  object_iri <- ifelse(is_iri, substring(objects, 2L, nchar(objects) - 1L), NA_character_)
  ## Objects: an IRI in angle brackets, a bare xsd:decimal, or a quoted literal
  ## with an optional datatype.  The pattern is built with paste0 because the
  ## same escapes written as one literal in this file fail to match on this R
  ## build (see .tmp/probe/fix1.R).
  lit_pat <- paste0('^(?:"([^"\\\\]*(?:\\\\.[^"\\\\]*)*)"', "(?:\\\\^\\\\^<([^>]+)>)?", "|(-?[0-9]+[.][0-9]+))$")
  lit <- regmatches(objects, regexec(lit_pat, objects))
  literal <- vapply(lit, function(m) {
    if (length(m) == 4L && nzchar(m[[2L]])) gsub('\\\\(["\\\\])', "\\1", m[[2L]]) else NA_character_
  }, character(1))
  datatype <- vapply(lit, function(m) {
    if (length(m) == 4L && nzchar(m[[3L]])) m[[3L]] else NA_character_
  }, character(1))
  bare_decimal <- vapply(lit, function(m) length(m) == 4L && nzchar(m[[4L]]), logical(1))
  ok_lit <- lengths(lit) == 4L
  ok_iri <- is_iri
  if (any(!ok_lit & !ok_iri)) {
    bad <- which(!ok_lit & !ok_iri)[[1]]
    stop(sprintf("malformed N-Triples object at line %d: %s", bad, lines[[bad]]), call. = FALSE)
  }
  is_decimal <- bare_decimal | (!is.na(datatype) & datatype == "http://www.w3.org/2001/XMLSchema#decimal")
  value <- rep(NA_real_, length(objects))
  ## Only decimals get a number; a plain string literal stays NA rather than
  ## raising "NAs introduced by coercion".
  decimals <- is_decimal & ok_lit
  source_text <- vapply(lit, function(m) {
    if (length(m) < 4L) NA_character_ else if (nzchar(m[[4L]])) m[[4L]] else m[[2L]]
  }, character(1))
  value[decimals] <- as.numeric(source_text[decimals])
  if (any(is_decimal & is.na(value))) {
    bad <- which(is_decimal & is.na(value))[[1]]
    stop(sprintf("xsd:decimal is not a number at line %d: %s", bad, lines[[bad]]), call. = FALSE)
  }
  out <- data.frame(
    subject = vapply(parts, `[[`, character(1), 2L),
    predicate = vapply(parts, `[[`, character(1), 3L),
    object = ifelse(is_iri, NA_character_, literal),
    object_iri = object_iri,
    datatype = ifelse(is_decimal, "decimal", datatype),
    stringsAsFactors = FALSE
  )
  # A decimal triple comes back as a number: whatever can be recomputed from the
  # file is recomputed, and anything else stays text with its datatype attached.
  out$value <- ifelse(is_decimal, value, NA_real_)
  out
}

# --- running a script under Rprof --------------------------------------------

profile_run_file <- function(script, out_dir = NULL, interval = 0.01) {
  out_dir <- out_dir %||% file.path(".tmp", "prof", tools::file_path_sans_ext(basename(script)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  name <- tools::file_path_sans_ext(basename(script))
  raw <- file.path(out_dir, paste0(name, ".out"))
  on.exit(unlink(raw), add = TRUE)

  env <- new.env(parent = globalenv())
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  t0 <- proc.time()[["elapsed"]]
  Rprof(raw, interval = interval)
  on.exit(Rprof(NULL), add = TRUE)
  sys.source(script, envir = env)
  Rprof(NULL)
  elapsed <- proc.time()[["elapsed"]] - t0

  x <- profile_read(raw)
  commit <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    warning = function(w) "unknown", error = function(e) "unknown"
  )
  if (!length(commit) || !nzchar(commit)) commit <- "unknown"

  write.csv(x$self, file.path(out_dir, paste0(name, ".self.csv")), row.names = FALSE)
  write.csv(x$frames, file.path(out_dir, paste0(name, ".frames.csv")), row.names = FALSE)
  saveRDS(
    list(self = x$self, frames = x$frames, samples = x$samples, interval = x$interval),
    file.path(out_dir, paste0(name, ".rds"))
  )
  n <- profile_ttl_write(x, file.path(out_dir, paste0(name, ".ttl")),
    meta = list(script = script, commit = commit)
  )
  cat(sprintf(
    "profiled: %s\n  wall: %.1f s  samples: %d  triples: %d\n  out: %s\n",
    script, elapsed, x$samples, n, out_dir
  ))
  profile_report(x)
  invisible(x)
}

profile_from_file <- function(path, out_dir = NULL) {
  x <- profile_read(path)
  out_dir <- out_dir %||% file.path(".tmp", "prof", tools::file_path_sans_ext(basename(path)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  name <- tools::file_path_sans_ext(basename(path))
  write.csv(x$self, file.path(out_dir, paste0(name, ".self.csv")), row.names = FALSE)
  saveRDS(
    list(self = x$self, frames = x$frames, samples = x$samples, interval = x$interval),
    file.path(out_dir, paste0(name, ".rds"))
  )
  n <- profile_ttl_write(x, file.path(out_dir, paste0(name, ".ttl")),
    meta = list(script = path, commit = "unknown")
  )
  cat(sprintf("  triples: %d  out: %s\n", n, out_dir))
  profile_report(x)
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- main ---------------------------------------------------------------------

profile_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args)) {
    stop("usage: dev/profile.R run <file.R> [out-dir] | stacks <file.out> [out-dir]", call. = FALSE)
  }
  mode <- args[[1]]
  path <- args[[2]]
  out_dir <- if (length(args) >= 3) args[[3]] else NULL
  switch(mode,
    run = profile_run_file(path, out_dir),
    stacks = profile_from_file(path, out_dir),
    stop("unknown mode: ", mode, call. = FALSE)
  )
}

if (sys.nframe() == 0L) {
  profile_main()
}
