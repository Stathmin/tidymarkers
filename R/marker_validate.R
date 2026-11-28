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
#
#  A position and a proportion are judged the same way - by the value the caller
#  gave, before any cast - and live here too: `bad_position()` (ADR 0045) and
#  `bad_proportion()` (ADR 0086, with `check_variant_weights()`).
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
  # `&` is not short-circuiting: a character column must not reach `round()`,
  # or the promised message is replaced by a base R error (ADR 0045 amendment).
  # A character column is judged by its spelling: as.integer("250.9") is 250
  # and as.integer("1e3") is 1000, so only the text can refuse them (ADR 0068).
  frac <- if (is.numeric(x)) {
    !is.finite(x) | x != round(x)
  } else if (is.character(x)) {
    !grepl("^ *[+-]?[0-9]+ *$", x)
  } else {
    # A factor, a logical, a complex or a list is not a spelling of a
    # position: `as.integer()` returns a factor code, a coerced 0/1 or a
    # truncated first element, so the unit is refused rather than guessed
    # (ADR 0045 amendment).
    rep(TRUE, length(x))
  }
  is.na(x) | is.na(iv) | iv < min | frac
}

#' Which values are not proportions?
#'
#' An allele frequency or a haplotype share is a proportion: one finite number
#' in `[0, 1]`.  A table in percent (`5`), a negative share and a spelling that
#' does not cast (`"5%"`) are not frequencies, yet the weight readers accepted
#' them and they changed the scores: `variant_load_in()` summed `5` for one
#' variant, a negative weight made a variant *better* than a clean site, and the
#' haplotype reader turned an unparseable share into `0`, silently removing that
#' haplotype from the population (ADR 0086).  The vectorised predicate is the one
#' definition; a caller keeps its own message.
#'
#' `NA` answers `FALSE`: a missing value is an absence, which each reader
#' resolves in its own documented way (a variant weight counts `1`, a haplotype
#' share is refused).  `NaN` is not a spelling of a frequency and answers `TRUE`,
#' like a factor, logical or complex column, whose `as.numeric()` would be a
#' code, a `0/1` or a lost imaginary part.
#' @param x a frequency vector.
#' @return a logical vector, `TRUE` where a present value is not a proportion.
#' @keywords internal
#' @noRd
bad_proportion <- function(x) {
  n <- length(x)
  if (is.factor(x) || is.logical(x) || is.complex(x) || is.list(x)) {
    return(rep(TRUE, n))
  }
  v <- if (is.character(x)) suppressWarnings(as.numeric(x)) else as.numeric(x)
  misspelled <- is.na(v) & !is.na(x)
  is.nan(x) | (!is.na(x) & (misspelled | !is.finite(v) | v < 0 | v > 1))
}

#' Check per-variant frequency weights
#'
#' `pop_freq` weights every variant's load and expected dropout, so a value
#' outside `[0, 1]` is refused before it scores anything; a missing weight is
#' the documented "count the variant once" fallback of
#' `variant_load_in()`, not a zero cost (ADR 0086).
#' @param x the weight column, or `NULL` when the table carries none.
#' @param arg the column name for the message.
#' @return the weights as a numeric vector with `NA` replaced by `1`, or `NULL`.
#' @keywords internal
#' @noRd
check_variant_weights <- function(x, arg = "pop_freq") {
  if (is.null(x)) {
    return(NULL)
  }
  bad <- bad_proportion(x)
  if (any(bad)) {
    abort_msg(
      "`{arg}` must be one frequency in [0, 1] per variant; row {which(bad)[[1L]]} is {format(x[bad][[1L]])}."
    )
  }
  w <- suppressWarnings(as.numeric(x))
  w[is.na(w)] <- 1
  w
}

#' Check one proportion argument
#'
#' A frequency floor is a proportion like a frequency is: `min_alt_freq = -1`
#' made the population check always pass, `5` made it always fail, and `NA`
#' reached base R's `if` with a message naming no argument (ADR 0086).
#' @param x the candidate proportion.
#' @param arg the argument name for the message.
#' @return the proportion as a single number.
#' @keywords internal
#' @noRd
check_proportion <- function(x, arg) {
  if (length(x) != 1L || is.na(x) || bad_proportion(x)) {
    abort_msg("`{arg}` must be one proportion in [0, 1]; got {format(x)}.")
  }
  as.numeric(x)
}

#' Check one coordinate before the read
#'
#' `as.integer()` truncates `2.9` to `2`, so a region read silently starts at a
#' base the caller never named; `NA` and `Inf` reach the engine's own error.
#' The check runs on the value the caller gave and returns the integer the reader
#' can use (ADR 0050).
#' @param x the candidate coordinate.
#' @param arg the argument name for the message.
#' @return the coordinate as a single integer.
#' @keywords internal
#' @noRd
check_coord <- function(x, arg) {
  if (length(x) != 1L || bad_position(x)) {
    abort_msg("`{arg}` must be one positive whole number; got {format(x)}.")
  }
  as.integer(x)
}

#' Check one or more nucleotide lengths before the design grid
#'
#' `as.integer()` truncates `20.9` to `20`, so a guide or spacer is designed
#' one base shorter than the caller named — a different cut position, not a
#' cosmetic difference (ADR 0051).  The check runs on the value the caller gave
#' and returns the integer vector the designer can use.
#' @param x the candidate lengths in nucleotides.
#' @param arg the argument name for the message.
#' @param multi whether a series is accepted; a single length is the default.
#' @return the lengths as an integer vector, in the caller's order.
#' @keywords internal
#' @noRd
check_nucleotide_length <- function(x, arg, multi = FALSE) {
  if (length(x) && (multi || length(x) == 1L) && !any(bad_position(x))) {
    return(as.integer(x))
  }
  if (multi) {
    abort_msg("`{arg}` must be one or more positive whole numbers of nucleotides; got {format(x)}.")
  }
  abort_msg("`{arg}` must be one positive whole number of nucleotides; got {format(x)}.")
}
