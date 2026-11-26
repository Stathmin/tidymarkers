# =============================================================================
#  marker_bncc.R  --  banded, non-compensatory, Chebyshev ranking (ZW4.6 phase C)
# =============================================================================
#  The lexicographic path in `marker_select.R` compares one tier at a time on
#  raw margins in mixed units (kcal/mol, C, bp, counts).  That is safe but it
#  cannot say how much a degree of Tm imbalance is worth against a base pair of
#  amplicon, and it hides the trade-off it makes.  This file is the declared
#  aggregation policy that replaces that ranking when a stage opts in:
#
#    L1  normalise every tier margin to a *regret* in [0, 1]  (`bncc_regret()`)
#    L2  group tiers into priority bands, compared lexicographically  (`BNCC_BANDS`)
#    L3  inside a band, score by the weighted Chebyshev (min-max) regret
#    L4  report the Pareto front over the band scores, not only the winner
#    L5  a calibration hook that stays inert until real wet-lab outcomes exist
#
#  Band membership is declared by POSITION in a stage's `tiers` vector.  The
#  stages disagree on what they call their tier columns (`tier1..` versus
#  `m_probe_tm`), so a position is the one stable key.  Every band partition is
#  validated against the `tiers` it is used with (`bncc_bands()`); a stage whose
#  tier count changes fails loudly instead of ranking on a stale policy.
#
#  This is a POLICY, not a threshold: `CRITERION_TIER_POLICY_EVIDENCE` names it
#  `rule_of_thumb` and no row's `evidence` rung moves because of it.
# =============================================================================

# Band membership per stage, ordered band by band: A = validity (would it work
# at all?), B = robustness (will it work across the sample?), C = operational
# (cost / convenience).  Positions index the stage's `tiers` vector, i.e. the
# order the designer passes to `select_with_ladder()`.  See proposal section 5.
BNCC_BANDS <- list(
  as = list(A = c(1L, 2L), C = c(3L, 4L)),
  pair = list(A = c(1L, 2L), B = 3L, C = c(4L, 5L)),
  duplex = list(A = 1L, B = 2L),
  asq = list(A = 1L),
  caps = list(A = 1L, B = c(2L, 3L)),
  arms = list(A = 1L, B = c(2L, 3L)),
  # candidate centrality is operational; the Tm balance is validity
  hrm = list(A = 2L, B = 1L),
  castpcr = list(A = 1L, B = c(2L, 3L)),
  crispr = list(A = c(1L, 2L, 3L), B = 4L),
  dpcr = list(A = c(1L, 2L), B = 3L),
  indel_size = list(A = 1L, B = 2L, C = 3L),
  # structure is validity, the inner-primer Tm set and the amplicon are robustness
  lamp = list(A = c(1L, 3L), B = c(2L, 4L)),
  lax = list(A = 1L),
  # discrimination placement and probe geometry are validity; length and amplicon robustness
  rpa = list(A = c(1L, 3L), B = c(2L, 4L)),
  # the Tm balance is validity; centrality and amplicon are operational
  sanger = list(A = 2L, B = c(1L, 3L)),
  snapshot = list(A = c(1L, 2L), B = 3L),
  ssr = list(A = 1L, B = 2L, C = 3L),
  taqman = list(A = c(1L, 2L), B = 3L)
)

# L5 -- the calibration hook.  The only fitted quantities are the weights and
# the band boundaries (proposal section 6).  Nothing here is fitted: `fitted` is FALSE
# and every weight is 1, so `bncc_rank()` is deterministic and no criterion row
# becomes `fitted_here`.  A fit lands by replacing this list with one whose
# `fitted` is TRUE and whose `weights` / `bands` carry the dataset hash and CIs.
BNCC_CALIBRATION <- list(
  fitted = FALSE,
  protocol = "monotone ordinal regression of wet-lab outcome on regret, bootstrap CIs (proposal section 6)",
  weights = list(),
  bands = list()
)

#' Band membership for a stage, as tier column names
#'
#' Falls back to one band per tier (the identity partition) for a stage that has
#' no declared bands, so an undeclared stage keeps the plain lexicographic order.
#' @param stage the selection stage id (e.g. `"pair"`).
#' @param tiers the stage's tier column names, in rank order.
#' @return a named list of character vectors (band -> tier column names).
#' @keywords internal
#' @noRd
bncc_bands <- function(stage, tiers) {
  declared <- if (is.null(stage)) NULL else BNCC_BANDS[[stage]]
  if (is.null(declared)) {
    out <- as.list(tiers)
    names(out) <- LETTERS[seq_along(out)]
    return(out)
  }
  pos <- unlist(declared, use.names = FALSE)
  if (!identical(sort(pos), seq_along(tiers))) {
    cli::cli_abort(c(
      "The BNCC bands for stage {.val {stage}} do not partition its {length(tiers)} tier(s).",
      "x" = "Declared positions: {.val {sort(pos)}}."
    ))
  }
  lapply(declared, function(i) tiers[i])
}

#' Normalise margins to a regret in `[0, 1]`
#'
#' `x` is a margin where larger is better.  `ideal` is the best value the tier
#' ranks toward (the printed target, when available); `nadir` is the worst value
#' on the candidate set.  With no ideal/nadir supplied the current candidate set
#' supplies them, which is the "normalise on the candidate set" of proposal L1.
#' A candidate at or beyond the ideal has zero regret; one at or below the nadir
#' has regret 1.  `NA` stays `NA` (an unknown margin is not evidence).
#' @keywords internal
#' @noRd
bncc_regret <- function(x, ideal = NULL, nadir = NULL) {
  if (!is.numeric(x)) {
    abort_msg("{.arg x} must be numeric.")
  }
  if (!length(x)) {
    return(numeric(0))
  }
  if (is.null(ideal)) ideal <- suppressWarnings(max(x, na.rm = TRUE))
  if (is.null(nadir)) nadir <- suppressWarnings(min(x, na.rm = TRUE))
  if (!is.finite(ideal) || !is.finite(nadir)) {
    return(rep(NA_real_, length(x)))
  }
  if (ideal < nadir) {
    abort_msg("{.arg ideal} ({ideal}) must not be below {.arg nadir} ({nadir}).")
  }
  denom <- ideal - nadir
  if (denom <= 0) {
    # no spread in this tier: every candidate is equally good, so no regret
    return(ifelse(is.na(x), NA_real_, 0))
  }
  r <- (ideal - x) / denom
  r <- pmin(pmax(r, 0), 1)
  r[is.na(x)] <- NA_real_
  r
}

# Row-wise max that skips `NA`, and stays `NA` only when every value is `NA`.
.row_max <- function(m) {
  if (!length(m)) {
    return(numeric(0))
  }
  if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
  apply(m, 1L, function(v) {
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  })
}

#' Band scores for a candidate table
#'
#' Each tier margin becomes a regret (`bncc_regret()`), then a band is scored by
#' the weighted Chebyshev (min-max) regret: the worst weighted regret in the
#' band.  The bottleneck therefore still dominates, but weights `w_j` are now
#' explicit and a missing tier is skipped rather than dropped.
#' @param data the graded candidate table.
#' @param tiers tier column names, in rank order.
#' @param bands a `bncc_bands()` result.
#' @param ideal,nadir optional named vectors (`tier -> value`) for L1.
#' @param weights optional named vector (`tier -> weight`); default 1.
#' @return a tibble with one `.bncc_<band>` column per band.
#' @keywords internal
#' @noRd
bncc_scores <- function(data, tiers, bands, ideal = NULL, nadir = NULL,
                        weights = NULL) {
  n <- nrow(data)
  if (!length(tiers)) {
    abort_msg("{.arg tiers} must name at least one tier.")
  }
  r <- lapply(tiers, function(t) {
    bncc_regret(data[[t]], ideal[[t]] %||% NULL, nadir[[t]] %||% NULL)
  })
  names(r) <- tiers
  w <- if (is.null(weights)) {
    rlang::set_names(rep(1, length(tiers)), tiers)
  } else {
    own <- rlang::set_names(rep(1, length(tiers)), tiers)
    own[names(weights)] <- weights
    own
  }
  out <- tibble::tibble(.row = seq_len(n))
  for (b in names(bands)) {
    cols <- bands[[b]]
    m <- do.call(cbind, lapply(cols, function(t) w[[t]] * r[[t]]))
    out[[paste0(".bncc_", b)]] <- .row_max(m)
  }
  out[, -1L, drop = FALSE]
}

#' Non-dominated rows over minimised columns
#'
#' A row is dominated when another is no worse in every column and strictly
#' better in one.  `NA` is treated as the worst value, so missing evidence can
#' never dominate a measured row.
#' @keywords internal
#' @noRd
bncc_pareto <- function(scores) {
  n <- nrow(scores)
  if (!n) {
    return(logical(0))
  }
  m <- as.matrix(scores)
  m[is.na(m)] <- Inf
  keep <- rep(TRUE, n)
  for (i in seq_len(n)) {
    dominated <- any(
      apply(m, 1L, function(v) all(v <= m[i, ]) && any(v < m[i, ]))
    )
    keep[i] <- !dominated
  }
  # a row cannot dominate itself; the loop above compares every row to row i,
  # and `all(v <= m[i, ]) && any(v < ...)` is FALSE for the row itself
  keep
}

#' Rank candidates by band priority, then weighted Chebyshev, then ties
#'
#' Lexicographic across bands (a candidate that wins band A always beats one that
#' loses it), ascending so smaller regret is better, `NA` last.  The winner is
#' the first row; the Pareto front over band scores is reported beside it rather
#' than collapsed into a single number.
#' @param data the graded candidate table.
#' @param stage the stage id used to look up `BNCC_BANDS` (ignored if `bands` given).
#' @param tiers tier column names, in rank order.
#' @param ties ascending tie-breaker columns.
#' @param bands optional explicit band partition (else `bncc_bands(stage, tiers)`).
#' @param ideal,nadir,weights passed to `bncc_scores()`.
#' @return a list: `data` (ranked, with `.pareto`), `bands`, `calibration`.
#' @keywords internal
#' @noRd
bncc_rank <- function(data, stage = NULL, tiers, ties = character(0), bands = NULL,
                      ideal = NULL, nadir = NULL, weights = NULL) {
  if (!nrow(data)) {
    return(list(data = data, bands = list(), calibration = BNCC_CALIBRATION))
  }
  if (is.null(bands)) bands <- bncc_bands(stage, tiers)
  scores <- bncc_scores(data, tiers, bands,
    ideal = ideal, nadir = nadir,
    weights = weights
  )
  tie_cols <- intersect(ties, names(data))
  ord <- do.call(
    order,
    unname(c(as.list(scores), as.list(data[tie_cols])))
  )
  out <- data[ord, , drop = FALSE]
  scores <- scores[ord, , drop = FALSE]
  out$.pareto <- bncc_pareto(scores)
  list(data = out, bands = bands, calibration = BNCC_CALIBRATION)
}
