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
