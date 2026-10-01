# =============================================================================
#  kasp_select.R  --  choosing among candidate assays: gate the grave, score the minor
# =============================================================================
#  A mature KASP designer is not a weighted sum.  Some flaws are fatal (a known
#  variant under a primer, an allele that cannot be told apart), others are just
#  costly.  So candidates are:
#
#    1. filtered by *hard gates* (the grave sins - any failure rejects);
#    2. ranked *lexicographically* on tier margins, i.e. the distance to the
#       threshold, tier by tier (a tier aggregates its dimensions by the
#       bottleneck, `min`, not the mean: an assay is as good as its weakest part).
#
#  When nothing clears the gates we relax them one at a time - least grave first
#  - and return the least-bad candidate with a warning per violated gate.  The
#  designer makes the trade-off; we never hide it.
# =============================================================================

# Which gates we are willing to give up, in order.  The last entries are the
# truly grave ones (a variant under a primer); they are relaxed only if even the
# milder compromises fail to yield a candidate.
KASP_RELAX_PAIR <- c(
  "gate_amplicon", "gate_structure", "gate_tm",
  "gate_discrimination", "gate_com_footprint", "gate_as_window"
)

# Human-readable name for each gate, used in the relaxation warning.
KASP_GATE_LABEL <- c(
  gate_discrimination = "allele discrimination below min_disc",
  gate_structure = "strong COM hairpin / homodimer",
  gate_tm = "COM Tm below the AS Tm",
  gate_amplicon = "amplicon outside the accepted length window",
  gate_com_footprint = "known variant under the COM footprint",
  gate_as_window = "known variant in the AS 3' window"
)

#' Grade candidate AS x COM pairs: hard gates + tiered margins
#'
#' Turns the flat inputs into the two things the selection needs: boolean gate
#' columns (the grave sins) and numeric margins (how far *inside* its threshold a
#' dimension sits - the minor sins, comparable across dimensions).  Margins are
#' already aggregated per tier by the bottleneck (`pmin`), so a pairing is scored
#' by its weakest link.
#' @param assays one row per (AS, COM) pairing.
#' @param limits thresholds (see `KASP_LIMITS`).
#' @param as_tm_mean mean AS Tm of the chosen AS placement.
#' @return a list with `data` (graded tibble), `gates`, `tiers` and `ties`.
#' @keywords internal
#' @noRd
grade_assay_pairs <- function(assays, limits, as_tm_mean) {
  tol <- limits$footprint_dropout %||% 0.05
  use_drop <- "com_dropout" %in% names(assays) && any(!is.na(assays$com_dropout))
  # footprint load: expected dropout when DECIPHER is available, else the
  # weighted count (still weights-aware, just position-blind)
  foot <- if (use_drop) {
    as.numeric(assays$com_dropout)
  } else if ("w_var" %in% names(assays)) {
    as.numeric(assays$w_var)
  } else {
    as.numeric(assays$n_var)
  }
  gate_foot <- if (use_drop) {
    is.na(assays$com_dropout) | assays$com_dropout <= tol
  } else {
    assays$n_var == 0L
  }
  use_win_drop <- "as_dropout" %in% names(assays) && any(!is.na(assays$as_dropout))
  win_margin <- if (use_win_drop) -as.numeric(assays$as_dropout) else -as.numeric(assays$as_n_3p)
  gate_win <- if (use_win_drop) {
    is.na(assays$as_dropout) | assays$as_dropout <= tol
  } else {
    assays$as_n_3p == 0L
  }
  d <- assays |>
    dplyr::mutate(
      .row = dplyr::row_number(),
      m_disc = .data$discrimination - limits$discrimination,
      m_structure = pmin(
        .data$com_homodimer - limits$homodimer_dg,
        .data$com_hairpin - limits$hairpin_dg
      ),
      m_tm_com = .data$com_tm - as_tm_mean,
      m_tm_tgt = -abs(.data$com_tm - COM_TM_TARGET),
      m_amplicon = pmin(
        .data$amplicon - limits$amplicon[1L],
        limits$amplicon[2L] - .data$amplicon
      ),
      m_footprint = -foot,
      m_window = win_margin,
      gate_discrimination = .data$discrimination >= limits$discrimination,
      gate_structure = .data$m_structure >= 0,
      gate_tm = .data$com_tm >= as_tm_mean,
      gate_amplicon = .data$amplicon >= limits$amplicon[1L] & .data$amplicon <= limits$amplicon[2L],
      gate_com_footprint = gate_foot,
      gate_as_window = gate_win,
      # the tiers themselves: lexicographic, and within a tier a bottleneck
      tier1 = .data$m_disc,
      tier2 = .data$m_structure,
      tier3 = pmin(.data$m_tm_com, .data$m_tm_tgt),
      tier4 = .data$m_amplicon,
      tier5 = .data$m_footprint
    )
  list(
    data = d,
    gates = c(
      "gate_discrimination", "gate_structure", "gate_tm", "gate_amplicon",
      "gate_com_footprint", "gate_as_window"
    ),
    tiers = c("tier1", "tier2", "tier3", "tier4", "tier5"),
    # tie-breakers, ascending (smaller is better)
    ties = c("run", "ssr", "com_start"),
    # the margin that measures each gate, so a relaxed gate can be minimised first
    gate_margin = c(
      gate_discrimination = "m_disc", gate_structure = "m_structure",
      gate_tm = "m_tm_com", gate_amplicon = "m_amplicon",
      gate_com_footprint = "m_footprint", gate_as_window = "m_window"
    )
  )
}

# AS-candidate gates, given up least-grave first (a wobble mismatch is milder
# than an allele-specific site lost to another variant)
KASP_RELAX_AS <- c("gate_wobble", "gate_discrimination", "gate_as_footprint", "gate_as_window")

#' Grade candidate AS primers: hard gates + tiered margins
#'
#' The AS-side twin of [grade_assay_pairs()]: the same vocabulary (gates for the
#' grave sins, margins ranked lexicographically) applied to one oligo on its own.
#' @param as_cand one row per AS placement.
#' @param limits thresholds (see `KASP_LIMITS`).
#' @param min_disc discrimination floor.
#' @return a list with `data` (graded tibble, carrying `as_row`), `gates`, `tiers`,
#'   `ties`.
#' @keywords internal
#' @noRd
grade_as_candidates <- function(as_cand, limits, min_disc) {
  tol <- limits$footprint_dropout %||% 0.05
  use_drop <- "as_dropout" %in% names(as_cand) && any(!is.na(as_cand$as_dropout))
  win_margin <- if (use_drop) -as.numeric(as_cand$as_dropout) else -as.numeric(as_cand$n_3p)
  gate_win <- if (use_drop) {
    is.na(as_cand$as_dropout) | as_cand$as_dropout <= tol
  } else {
    as_cand$n_3p == 0L
  }
  d <- as_cand |>
    dplyr::mutate(
      as_row = dplyr::row_number(),
      m_disc = .data$discrimination - min_disc,
      m_tm = -abs(.data$tm - AS_TM_TARGET),
      m_stability = -as.numeric(.data$dg_3p),
      m_footprint = -as.numeric(.data$w_var),
      m_window = win_margin,
      m_wobble = -as.numeric(.data$wobble),
      gate_discrimination = .data$discrimination >= min_disc,
      gate_wobble = is.na(.data$wobble) | !.data$wobble,
      gate_as_footprint = .data$n_var == 0L,
      gate_as_window = gate_win,
      tier1 = .data$m_disc,
      tier2 = .data$m_tm,
      tier3 = .data$m_stability,
      tier4 = .data$m_footprint
    )
  list(
    data = d,
    gates = c("gate_discrimination", "gate_wobble", "gate_as_footprint", "gate_as_window"),
    tiers = c("tier1", "tier2", "tier3", "tier4"),
    ties = c("seed_hits", "as_length", "as_start"),
    gate_margin = c(
      gate_wobble = "m_wobble", gate_discrimination = "m_disc",
      gate_as_footprint = "m_footprint", gate_as_window = "m_window"
    )
  )
}

#' Sort candidates lexicographically on tiers then ties (no gating)
#'
#' Used to build shortlists (e.g. which placements to seed-screen) with the same
#' ordering the final selection uses.
#' @keywords internal
#' @noRd
rank_candidates <- function(data, tiers, ties) {
  d <- data
  for (t in tiers) d[[paste0(".neg_", t)]] <- -d[[t]]
  dplyr::arrange(d, !!!rlang::syms(paste0(".neg_", tiers)), !!!rlang::syms(intersect(ties, names(d))))
}

#' Pick the best candidate, relaxing gates only as far as necessary
#'
#' Applies the gates, and if none survives drops them in `relax_order` (least
#' grave first) until at least one does.  Among the survivors, ranks
#' lexicographically on `tiers` (descending) then `ties` (ascending).  Reports
#' which gates were relaxed and which the winner still violates.
#' @param data the graded tibble (see `grade_assay_pairs()`).
#' @param gates gate column names (all must be TRUE to pass).
#' @param relax_order gate names in the order they may be given up.
#' @param tiers margin columns, ranked high-to-low, in order.
#' @param ties column names for ascending tie-breaks.
#' @return list with `best` (one row, or `NULL`), `relaxed` and `violated`.
#' @keywords internal
#' @noRd
select_with_ladder <- function(data, gates, relax_order, tiers, ties, gate_margin = NULL) {
  active <- gates
  relaxed <- character(0)
  repeat {
    ok <- if (length(active)) {
      Reduce(`&`, lapply(active, function(g) data[[g]]))
    } else {
      rep(TRUE, nrow(data))
    }
    ok[is.na(ok)] <- FALSE
    if (nrow(data) && any(ok)) {
      surv <- data[ok, , drop = FALSE]
      # a relaxed gate is the compromise we accepted: minimise it first (gravest
      # first), so the "least-bad" candidate is the least compromised one
      extra <- if (length(relaxed) && !is.null(gate_margin)) {
        unname(gate_margin[intersect(rev(relaxed), names(gate_margin))])
      } else {
        character(0)
      }
      ranked <- c(extra, tiers)
      for (t in ranked) surv[[paste0(".neg_", t)]] <- -surv[[t]]
      ord <- dplyr::arrange(
        surv,
        !!!rlang::syms(paste0(".neg_", ranked)),
        !!!rlang::syms(intersect(ties, names(surv)))
      )
      best <- ord[1L, , drop = FALSE]
      violated <- gates[!vapply(gates, function(g) isTRUE(best[[g]]), logical(1))]
      return(list(best = best, relaxed = relaxed, violated = unique(violated)))
    }
    if (!length(relax_order)) {
      return(list(best = NULL, relaxed = relaxed, violated = relaxed))
    }
    g <- relax_order[[1L]]
    relax_order <- relax_order[-1L]
    active <- setdiff(active, g)
    relaxed <- c(relaxed, g)
  }
}

#' Warn (once) about the gates the chosen assay had to relax
#' @keywords internal
#' @noRd
warn_relaxed <- function(violated) {
  if (!length(violated)) {
    return(invisible(NULL))
  }
  labels <- unname(KASP_GATE_LABEL[violated])
  labels[is.na(labels)] <- violated[is.na(labels)]
  rlang::warn(c(
    gettext("No assay cleared every hard constraint; returned the least-bad after relaxing some."),
    stats::setNames(gettext(labels), rep("*", length(labels)))
  ))
  invisible(NULL)
}
