# =============================================================================
#  marker_validate.R  --  one rule for a count argument
# =============================================================================
#  A count is a single whole number inside a range.  Four gates spelled that
#  rule out by hand (`max_markers`, `spacer_len`, `min_lines`, `min_evidence`)
#  and drifted: `marker_one_pot(max_markers = Inf)` left the `if` condition at
#  `NA` and raised base R's "missing value where TRUE/FALSE needed" instead of
#  the package message, and `min_evidence = Inf` was accepted as a count.  The
#  predicate below is the one definition (ADR 0044); a caller keeps its own
#  message and does its own `as.integer()` afterwards.
# =============================================================================

#' Is a value a single whole number inside a range?
#'
#' `NA`, `NaN`, `Inf`, a character and a vector all answer `FALSE`, so an
#' `if (!is_whole_number(x))` gate reaches its own message instead of failing
#' on `NA`.  `is.numeric()` runs first and `&&` short-circuits, so the range
#' comparisons never see a non-number.
#' @param x the candidate value.
#' @param min smallest accepted value; `1L` for a positive count, `0L` for a
#'   non-negative one.
#' @param max largest accepted value that still fits an integer, by default.
#' @return a single logical.
#' @keywords internal
#' @noRd
is_whole_number <- function(x, min = 1L, max = .Machine$integer.max) {
  is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= min && x <= max && x == round(x)
}

#' Check a whole-number search radius before the cast
#'
#' `as.integer()` truncates `1.9` to `1`, so a search silently uses a narrower
#' radius than the caller named; `-1` and `Inf` reach the engine's own error.
#' The check runs on the value the caller gave and returns the integer the
#' engine can use (ADR 0047).
#' @param x the candidate radius.
#' @param arg the argument name for the message.
#' @param min smallest accepted value; `0L` for a mismatch radius, `1L` for a
#'   seed length.
#' @return the radius as a single integer.
#' @keywords internal
#' @noRd
check_scan_radius <- function(x, arg, min = 0L) {
  if (!is_whole_number(x, min = min)) {
    abort_msg(
      "`{arg}` must be one {if (min > 0L) \"positive\" else \"non-negative\"} whole number."
    )
  }
  as.integer(x)
}

#' Check a base-pair threshold before the cast
#'
#' `as.integer()` truncates `5.9` to `5`, so a design gate silently accepts a
#' winner its caller's threshold rejects (ADR 0048).  The check runs on the
#' value the caller gave and returns the integer the designer can use.
#' @param x the candidate number of base pairs.
#' @param arg the argument name for the message.
#' @return the threshold as a single integer.
#' @keywords internal
#' @noRd
check_base_pairs <- function(x, arg) {
  if (!is_whole_number(x, min = 1L)) {
    abort_msg("`{arg}` must be one positive whole number of base pairs.")
  }
  as.integer(x)
}

#' Check a worker count before the cast
#'
#' `as.integer(0.5)` is `0` and `as.integer(Inf)` is `NA`, so a truncated worker
#' count stalls an external build or reaches it as a missing value (ADR 0049).
#' The check runs on the value the caller gave and returns the integer the
#' engine can use.
#' @param x the candidate number of workers.
#' @param arg the argument name for the message.
#' @return the count as a single integer.
#' @keywords internal
#' @noRd
check_worker_count <- function(x, arg) {
  if (!is_whole_number(x, min = 1L)) {
    abort_msg("`{arg}` must be one positive whole number of workers.")
  }
  as.integer(x)
}

#' Which positions are not positive whole numbers?
#'
#' The vectorised counterpart of `is_whole_number()` for a coordinate column.
#' `as.integer()` truncates, so a normaliser that casts first can never see a
#' fraction: `250.9` becomes `250` and matches a target at another locus
#' (ADR 0045).  A character `"250"` is accepted when it casts exactly; `NA`,
#' `Inf`, a negative and `250.9` are not.
#' @param x a coordinate vector.
#' @param min smallest accepted value.
#' @return a logical vector, `TRUE` where the value is not a usable position.
#' @keywords internal
#' @noRd
bad_position <- function(x, min = 1L) {
  iv <- suppressWarnings(as.integer(x))
  is.na(x) | is.na(iv) | iv < min |
    (is.numeric(x) & (!is.finite(x) | x != round(x)))
}
