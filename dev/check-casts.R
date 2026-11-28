#!/usr/bin/env Rscript
# =============================================================================
#  dev/check-casts.R  --  every position cast is guarded or registered
# =============================================================================
#  A position read from a file, a user table or a pointer is not necessarily a
#  whole number: `as.integer()` truncates `10.9` to `10` and turns `-5` into a
#  plausible-looking start.  ZV96-ZV104 closed that hole one reader at a time;
#  the rule behind them is ADR 0068 and its siblings: judge a position by its
#  spelling before the cast.
#
#  This gate makes the roster checked, not remembered: a top-level function in
#  R/ that casts a position-like expression with `as.integer()` must call a
#  guard *for that expression* first, or be registered in `dev/casts.yml` with
#  a reason.  A guard for `pos` does not bless a later `as.integer(start)`:
#  the unit-wide test it replaced let one checked argument hide an unchecked
#  neighbour (ZV110).  The scan is
#  deliberately narrow: only expressions that name `pos`, `POS`, `start`, `end`,
#  `offset`, `positions`, `pos_raw` or `stab_at` are position-like, and only the
#  guards below count.  A cast on an internal integer object still needs an
#  entry, so the ledger says why it cannot receive a fraction.
#
#  A base-pair or copy *threshold* is the same question: `as.numeric()` let
#  `min_band_gap = 20.9` narrow the arms floor and `NA` reach a base R `if`
#  with a message naming no argument (ADR 0081).  A direct cast of an argument
#  named `min_*`/`max_*` is therefore a cast of the same kind; a threshold
#  merely mentioned inside a comparison casts the comparison's logical result,
#  not the threshold, and is not reported.
#
#  Fails when:
#    - an unguarded cast has no ledger entry        (new unguarded cast)
#    - a ledger entry no longer has a cast          (stale justification)
#    - an entry has no reason
#  Passes a built-in self-test on every run, so a scanner that stopped matching
#  is caught.
#
#  Run:  Rscript dev/check-casts.R
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

# The names a position can be spelled with.  `start`/`end` are included: an
# interval bound comes from a file or a caller's table exactly like `pos`
# (ZV106), and the earlier exclusion hid `normalise_features()`.  A cast on an
# internal interval is registered, not exempted.
position_name <- c("pos", "POS", "pos_raw", "positions", "stab_at", "start", "end")

# A count of base pairs or copies carries the same whole-in-its-spelling rule
# when it gates a design (ADR 0081).
threshold_name <- "^(min|max)_"

# A guard refuses or normalises a position; it does not merely report one.
guard_primitive <- c(
  "bad_position", "check_target_pos", "check_depth_pos", "check_coord",
  "check_scan_radius", "check_nucleotide_length", "is_whole_number",
  "check_base_pairs", "check_worker_count"
)

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

# Every function definition in a parsed file, named by its assignment path.
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

# TRUE when the expression names a position, directly or through `$`/`[[`.
mentions_position <- function(e) {
  if (is.name(e)) {
    return(as.character(e) %in% position_name)
  }
  if (!is.call(e)) {
    return(FALSE)
  }
  head <- call_name(e)
  if (head %in% c("$", "[[") && length(e) >= 3L) {
    tag <- tryCatch(e[[3L]], error = function(err) NULL)
    if (is.name(tag) && nzchar(as.character(tag)) && as.character(tag) %in% position_name) {
      return(TRUE)
    }
    if (is.character(tag) && any(tag %in% position_name)) {
      return(TRUE)
    }
  }
  for (i in seq_along(e)) {
    if (mentions_position(e[[i]])) {
      return(TRUE)
    }
  }
  FALSE
}

# TRUE when the expression *is* a threshold name (or a `$`/`[[` of one).
# A threshold nested inside a comparison is not a cast of the threshold: the
# value cast there is the comparison's logical result.
directly_names_threshold <- function(e) {
  if (is.name(e)) {
    return(grepl(threshold_name, as.character(e)))
  }
  if (!is.call(e)) {
    return(FALSE)
  }
  head <- call_name(e)
  if (head %in% c("$", "[[") && length(e) >= 3L) {
    tag <- tryCatch(e[[3L]], error = function(err) NULL)
    if (is.name(tag)) {
      return(grepl(threshold_name, as.character(tag)))
    }
    if (is.character(tag)) {
      return(any(grepl(threshold_name, tag)))
    }
  }
  FALSE
}

# Cast and guard calls in one unit, in body order, not descending into a
# nested definition: that definition is its own unit and needs its own guard.
scan_calls <- function(expr, found = list()) {
  if (!is.call(expr)) {
    return(found)
  }
  name <- call_name(expr)
  if (!is.na(name) && name %in% guard_primitive) {
    found <- c(found, list(guard_call(expr, name)))
  }
  position_cast <- identical(name, "as.integer") && length(expr) == 2L &&
    mentions_position(expr[[2L]])
  threshold_cast <- name %in% c("as.integer", "as.numeric", "as.double") &&
    length(expr) == 2L && directly_names_threshold(expr[[2L]])
  if (position_cast || threshold_cast) {
    found <- c(found, list(list(kind = "cast", name = name, arg = deparse_text(expr[[2L]]))))
  }
  if (identical(name, "function")) {
    return(found)
  }
  for (i in seq_along(expr)) found <- scan_calls(expr[[i]], found)
  found
}

deparse_text <- function(e) paste(deparse(e), collapse = "")

# One guard call: the expressions it judges and the labels it names them by.
# `check_base_pairs(min_gap, "min_gap")` judges `min_gap` twice: as a value and
# as the name in its message, so both spellings cover a later cast.
guard_call <- function(expr, name) {
  args <- as.list(expr)[-1L]
  labels <- vapply(
    args,
    function(a) if (is.character(a) && length(a) == 1L && !is.na(a)) a else NA_character_,
    character(1)
  )
  list(
    kind = "guard", name = name,
    args = vapply(args, deparse_text, character(1)),
    labels = labels[!is.na(labels)]
  )
}

# TRUE when one guard judged the expression a cast reads.  The guard's
# expression must match the cast's text (the same value), its message label
# must name it, or the cast must index a vector the guard judged whole
# (`d$POS[idx]` after `bad_position(d$POS)`).
guard_covers <- function(guard, cast_arg) {
  if (identical(guard$kind, "cast")) {
    return(FALSE)
  }
  any(guard$args == cast_arg) ||
    any(guard$labels == cast_arg) ||
    any(startsWith(cast_arg, paste0(guard$args, "[")))
}

scan_casts <- function(bodies) {
  out <- list()
  for (symbol in names(bodies)) {
    calls <- scan_calls(bodies[[symbol]])
    kinds <- vapply(calls, `[[`, character(1), "kind")
    if (!"cast" %in% kinds) next
    cast_at <- which(kinds == "cast")
    guarded <- vapply(
      cast_at,
      function(i) {
        before <- calls[seq_len(i - 1L)]
        any(vapply(before, guard_covers, logical(1), cast_arg = calls[[i]]$arg))
      },
      logical(1)
    )
    out[[symbol]] <- list(casts = length(cast_at), guarded = all(guarded))
  }
  out
}

cast_registry <- function() {
  yaml::read_yaml("dev/casts.yml")$entries %||% list()
}

check_casts <- function(registry, hits) {
  problems <- character()
  symbols <- vapply(registry, `[[`, character(1), "symbol")
  reason <- vapply(registry, function(e) e$reason %||% "", character(1))
  no_reason <- symbols[!nzchar(reason)]
  if (length(no_reason)) {
    problems <- c(problems, paste0("cast entry without a reason: ", no_reason))
  }
  unguarded <- names(hits)[!vapply(hits, `[[`, logical(1), "guarded")]
  unregistered <- setdiff(unguarded, symbols)
  if (length(unregistered)) {
    problems <- c(
      problems,
      paste0("unguarded position/threshold cast (guard it, or add to dev/casts.yml): ", unregistered)
    )
  }
  stale <- setdiff(symbols, unguarded)
  if (length(stale)) {
    problems <- c(
      problems,
      paste0("dev/casts.yml is stale (the symbol no longer casts unguarded): ", stale)
    )
  }
  problems
}

self_test <- function() {
  bodies <- c(
    unguarded = "function(x) as.integer(x$pos)",
    threshold = "function(x) as.numeric(x$min_gap)",
    threshold_guarded = "function(x) { check_base_pairs(x$min_gap, \"min_gap\"); as.numeric(x$min_gap) }",
    thermophilic = "function(x) as.numeric(x$as_tm)",
    unguarded_interval = "function(x) as.integer(x$start)",
    fractional = "function(x) as.integer(x$pos + 0.5)",
    guarded = "function(x) { bad_position(x$pos); as.integer(x$pos) }",
    guarded_interval = "function(x) { bad_position(x$end); as.integer(x$end) }",
    guarded_index = "function(x) { bad_position(x$POS); as.integer(x$POS[idx]) }",
    cross_guard = "function(x) { bad_position(x$pos); as.integer(x$start) }",
    late = "function(x) { as.integer(x$pos); check_target_pos(x$pos) }",
    integer_internal = "function(x) as.integer(x$n)",
    labelled = "function(x) as.integer(x$count)",
    strings = "function(x) paste0(\"as.integer(x$pos)\")"
  )
  exprs <- Map(
    function(name, text) parse(text = paste0(name, " <- ", text))[[1L]],
    names(bodies), bodies
  )
  units <- list()
  for (expr in exprs) units <- collect_units(expr, "", units)
  hits <- scan_casts(units)
  empty <- check_casts(list(), hits)
  expected <- paste0(
    "unguarded position/threshold cast (guard it, or add to dev/casts.yml): ",
    c("cross_guard", "fractional", "late", "threshold", "unguarded", "unguarded_interval")
  )
  if (!identical(sort(empty), sort(expected))) {
    stop(sprintf(
      "self-test: an unguarded position cast was not caught (got: %s)",
      paste(sort(empty), collapse = ", ")
    ))
  }
  if (any(c("integer_internal", "labelled", "strings", "thermophilic") %in% names(hits))) {
    stop("self-test: a non-position/non-threshold cast was reported")
  }
  if (!isTRUE(hits[["guarded_index"]]$guarded)) {
    stop("self-test: a guard on the indexed vector did not cover its cast")
  }
  if (!isTRUE(hits[["threshold_guarded"]]$guarded)) {
    stop("self-test: a guarded threshold cast was reported as unguarded")
  }
  ok <- check_casts(
    list(
      list(symbol = "cross_guard", reason = "self-test"),
      list(symbol = "unguarded", reason = "self-test"),
      list(symbol = "unguarded_interval", reason = "self-test"),
      list(symbol = "fractional", reason = "self-test"),
      list(symbol = "late", reason = "self-test"),
      list(symbol = "threshold", reason = "self-test")
    ),
    hits
  )
  if (length(ok)) {
    stop("self-test: a guarded body or a registered cast was flagged")
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

r_functions <- function(files = list.files("R", "\\.R$", full.names = TRUE)) {
  out <- list()
  for (file in files) {
    exprs <- tryCatch(parse(file), error = function(e) NULL)
    if (is.null(exprs)) next
    for (expr in as.list(exprs)) out <- collect_units(expr, "", out)
  }
  out
}

self_test()
hits <- scan_casts(r_functions())
registry <- cast_registry()
problems <- check_casts(registry, hits)

if (length(problems)) {
  cat("cast gate failed:\n")
  cat(paste0("  - ", problems, "\n"), sep = "")
  quit(status = 1L)
}
cat(
  "  ok   ", length(hits), " position/threshold cast(s), ",
  sum(vapply(hits, `[[`, logical(1), "guarded")), " guarded in place, ",
  length(registry), " registered with a reason; self-test passed\n",
  sep = ""
)
