# =============================================================================
#  zzz.R  --  package session state
# =============================================================================
#  A holder for one-shot user hints (e.g. "OligoArrayAux missing"), so the
#  package never writes a `options()` entry or a global variable to say "already
#  told them".  A package-namespace environment is the standard place for such
#  state (memoise caches use the same idea).
# =============================================================================

#' A fresh, empty session state
#'
#' Holds one-shot user hints (e.g. "OligoArrayAux missing"), the cache of
#' genome-derived, pattern-independent facts about a sequence (the runs of
#' unambiguous bases used by `exact_pattern_hits()`, keyed by a content
#' signature), and the session memo.
#'
#' The memo holds values that are *stable within one session*: an external tool's
#' location (a handful of directory stats) and the discovered scheme registry
#' (every `scheme_def_*()` called and turned into a tibble).  Both are asked for
#' per oligo / per candidate, so the cost is paid once instead of thousands of
#' times (a 6-target KASP design asked for the scheme rows ~1000 times).
#' Invalidated in exactly two places, both of which can change an answer:
#' `marker_install_tools()` (a tool appeared) and `marker_technology()` (a scheme
#' was composed).  Anything whose inputs vary per call - Tm, dG, specificity,
#' GC content - never goes here.
#'
#' `source_carrier()` re-homes this package's functions into a plain environment
#' for a socket worker, and this state must not travel by value; the re-homed
#' functions still look the name up, so the carrier gets a fresh one (the memo is
#' cold there either way).
#' @return an environment with `$blocks` and `$memo`.
#' @keywords internal
#' @noRd
marker_state_new <- function() {
  state <- new.env(parent = emptyenv())
  state$blocks <- new.env(parent = emptyenv())
  state$memo <- new.env(parent = emptyenv())
  state
}

.marker_state <- marker_state_new()

#' Memoise a session-stable value
#' @param key a string identifying the computation; include any option or
#'   working-directory value the answer depends on.
#' @param expr computed on the first call only.
#' @keywords internal
#' @noRd
session_memo <- function(key, expr) {
  if (exists(key, envir = .marker_state$memo, inherits = FALSE)) {
    return(get(key, envir = .marker_state$memo, inherits = FALSE))
  }
  value <- expr
  # `NULL` is a real answer ("tool not installed"), so the assignment is
  # unconditional; `exists()` - not `is.null()` - is what tells the two apart.
  assign(key, value, envir = .marker_state$memo)
  value
}

#' Drop every memoised session value
#' @keywords internal
#' @noRd
session_memo_clear <- function() {
  rm(list = ls(.marker_state$memo, all.names = TRUE), envir = .marker_state$memo)
  invisible(NULL)
}
