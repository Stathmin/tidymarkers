# =============================================================================
#  dev/z8-harness.R  --  the Z8.4 split and the gate-rate honesty check
# =============================================================================
#  Two plain functions, no side effects, no package internals.  The harness
#  (dev/targets/_targets.R) sources them so the split lives in one place; the
#  tests source the same file instead of re-declaring it.
#
#  The split is by SOURCE, never by row: one paper with 30 working rows and one
#  failing row would otherwise leak its own success count across the split.
#  The rate is returned only when both a `worked` and a non-`worked` terminal
#  class exist; "all worked" is a discriminable-free table, not a rate of 1.

# deterministic whole-source split: sorted unique sources, first
# ceiling(share * n) go to `unique`, the rest to `held`
z8_split_by_source <- function(tbl, share = 0.5) {
  if (!"source" %in% names(tbl)) stop("tbl has no `source` column")
  if (!is.numeric(share) || length(share) != 1L || is.na(share) ||
    share < 0 || share > 1) {
    stop("share must be one number in [0, 1]")
  }
  s <- sort(unique(as.character(tbl$source)))
  k <- seq_along(s) <= ceiling(share * length(s))
  list(unique = s[k], held = s[!k])
}

# only a table with at least one `worked` AND at least one other terminal
# class can discriminate gate-pass from gate-pass-survivors
z8_gate_rate <- function(tbl) {
  if (!"outcome" %in% names(tbl)) stop("tbl has no `outcome` column")
  o <- as.character(tbl$outcome)
  n_worked <- sum(!is.na(o) & o == "worked")
  n_failed <- sum(!is.na(o) & o == "failed")
  n_other <- sum(!is.na(o) & o != "worked" & o != "failed")
  computable <- n_worked > 0L && (n_failed > 0L || n_other > 0L)
  reason <- if (!computable) {
    if (n_worked == 0L && n_failed + n_other == 0L) {
      "no terminal outcomes: there is nothing to compute a rate from"
    } else if (n_worked == 0L) {
      "0 worked outcomes: a rate needs both a numerator and a denominator"
    } else {
      paste0(
        "0 failed outcomes: gate-pass cannot be discriminated from ",
        "gate-pass-survivors"
      )
    }
  } else {
    NA_character_
  }
  list(
    rate = if (computable) n_worked / (n_worked + n_failed + n_other) else NA_real_,
    n_worked = n_worked,
    n_failed = n_failed,
    n_other = n_other,
    computable = computable,
    reason = reason
  )
}
