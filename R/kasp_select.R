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

# Which gates may be *somewhat* relaxed, least grave first.  Gates absent here
# are hard limits set by the assay's chemistry / KASP literature and are never
# lowered (a lone variant in the AS 3' window, or no allele discrimination, makes
# the assay meaningless).  NB: a *linked* variant in the AS 3' window is not such
# a violation - it rides with the ALT haplotype and is baked into the ALT primer,
# so the window count excludes it (see assemble_kasp_assay(linked=)).
KASP_RELAX_PAIR <- c("gate_amplicon", "gate_structure", "gate_tm", "gate_com_footprint")
KASP_RELAX_AS <- c("gate_wobble", "gate_as_footprint")

# The gate set: the AS stage decides discrimination and the 3' window (they are
# fixed by the AS placement), so the pair stage only gates what varies across
# COM placements - structure, Tm, amplicon and the COM footprint.  Re-gating the
# AS's constants here would be vacuous.
KASP_GATES_PAIR <- c("gate_structure", "gate_tm", "gate_amplicon", "gate_com_footprint")
KASP_GATES_AS <- c("gate_discrimination", "gate_wobble", "gate_as_footprint", "gate_as_window")

# How far a gate's margin bar may be lowered (floor, in margin units: >= 0 passes
# at the original threshold).  A floor of 0 = hard (never relaxed); a negative
# floor bounds the relaxation (never -Inf).  gate_com_footprint's floor depends on
# the limits/backend, so it is filled in by grade_assay_pairs().
KASP_GATE_FLOOR_AS <- c(
  gate_discrimination = 0, gate_wobble = -1, gate_as_footprint = -1, gate_as_window = 0
)
KASP_GATE_FLOOR_PAIR <- c(
  gate_discrimination = 0, gate_structure = -3, gate_tm = -2,
  gate_amplicon = -20, gate_com_footprint = 0, gate_as_window = 0
)

# The lexicographic tiers and their tie-breakers, named (label = meaning).  The
# grade_*() functions and the decision graph both read these, so adding a tier or
# changing the order shows up in the picture.
KASP_TIERS_PAIR <- c(
  tier1 = "structure (hairpin / homodimer)", tier2 = "Tm (COM vs AS / target)",
  tier3 = "amplicon length", tier4 = "footprint load"
)
KASP_TIES_PAIR <- c("run", "ssr", "com_start")
KASP_TIERS_AS <- c(
  tier1 = "discrimination", tier2 = "Tm", tier3 = "3' stability", tier4 = "footprint load"
)
KASP_TIES_AS <- c("seed_hits", "as_length", "as_start")

# A duplex chemistry (dKASP / mxKASP / Multi-KASP) puts two loci in one tube, so
# on top of the per-locus design it adds a compatibility stage.  The rules are
# the published ones: Multi-KASP bans >= 6 bp inter-primer complementarity, keeps
# every primer Tm in 59-65 C and each amplicon <= 150 bp; dKASP pairs homologous
# genes, mxKASP non-homologous ones (zygosity then follows from the combination).
# A shared refinement branch may widen the amplicon or the self-complementarity
# cap, but never the cross-talk rules.  Readout may be two cassettes (co-dominant
# per locus) or presence/absence (dominant) - see kasp_schemes().
KASP_GATES_DUPLEX <- c(
  "gate_duplex_targets", "gate_duplex_tm", "gate_duplex_complement",
  "gate_duplex_selfcomp", "gate_duplex_amplicon"
)
KASP_RELAX_DUPLEX <- c("gate_duplex_amplicon", "gate_duplex_selfcomp")
KASP_TIERS_DUPLEX <- c(
  tier1 = "duplex pair compatibility", tier2 = "Tm window match", tier3 = "amplicon budget"
)

# ASQ / PACE: a different cassette readout (Kalendar 2022), up to 4-plex in one
# real-time reaction.  The allele-specific core is the KASP one; the extra gates
# are the cassette budget and the real-time resolution.
KASP_GATES_ASQ <- c("gate_asq_cassette", "gate_asq_readout")
KASP_RELAX_ASQ <- character(0)
KASP_TIERS_ASQ <- c(tier1 = "allele discrimination", tier2 = "Tm")

# CAPS / dCAPS: the variant is read by a restriction digest (R/kasp_caps.R), so
# selection is "is there an enzyme that separates the alleles, and does the pair
# amplify the site cleanly" - none of the KASP 3'-terminal machinery applies.
KASP_GATES_CAPS <- c("gate_caps_site", "gate_caps_resolvable", "gate_caps_amplicon", "gate_caps_tm")
KASP_RELAX_CAPS <- c("gate_caps_amplicon", "gate_caps_tm", "gate_caps_resolvable")
KASP_TIERS_CAPS <- c(
  tier1 = "band gap (gel resolvability)", tier2 = "Tm balance", tier3 = "amplicon length"
)

# Human-readable name for each gate, used in the relaxation warning and the
# decision graph.
KASP_GATE_LABEL <- c(
  gate_discrimination = "allele discrimination below min_disc",
  gate_structure = "strong COM hairpin / homodimer",
  gate_tm = "COM Tm below the AS Tm",
  gate_amplicon = "amplicon outside the accepted length window",
  gate_com_footprint = "known variant under the COM footprint",
  gate_as_window = "known variant in the AS 3' window",
  gate_wobble = "wobble (G:T / A:C) mismatch at the 3' base",
  gate_as_footprint = "AS footprint dropout above tolerance",
  gate_duplex_targets = "the two loci do not match the chosen readout",
  gate_duplex_tm = "primers outside the shared Tm window",
  gate_duplex_complement = "inter-pair primer complementarity",
  gate_duplex_selfcomp = "primer self-complementarity above the duplex cap",
  gate_duplex_amplicon = "a duplex amplicon above 150 bp",
  gate_asq_cassette = "more targets than the ASQ cassette budget",
  gate_asq_readout = "no allele resolution in real time",
  gate_caps_site = "no enzyme separates the REF and ALT amplicons",
  gate_caps_resolvable = "the diagnostic bands are too close on a gel",
  gate_caps_amplicon = "the CAPS amplicon is outside the accepted window",
  gate_caps_tm = "a flanking primer Tm is outside the window"
)

# Short names for the decision graph (the warnings keep the full labels).
KASP_GATE_SHORT <- c(
  gate_discrimination = "discrimination", gate_structure = "structure",
  gate_tm = "COM Tm >= AS Tm", gate_amplicon = "amplicon length",
  gate_com_footprint = "COM footprint", gate_as_window = "AS 3' window",
  gate_wobble = "wobble mismatch", gate_as_footprint = "AS footprint",
  gate_duplex_targets = "duplex targets", gate_duplex_tm = "Tm window (all primers)",
  gate_duplex_complement = "inter-pair complementarity",
  gate_duplex_selfcomp = "self-complementarity",
  gate_duplex_amplicon = "amplicon <= 150 bp",
  gate_asq_cassette = "ASQ cassette", gate_asq_readout = "real-time readout",
  gate_caps_site = "restriction site", gate_caps_resolvable = "band gap",
  gate_caps_amplicon = "amplicon length (CAPS)", gate_caps_tm = "primer Tm (CAPS)"
)

# Human-readable relaxation floor per soft gate (for the graph and text export).
KASP_GATE_FLOOR_LABEL <- c(
  gate_amplicon = "widen window by 20 bp",
  gate_structure = "3 kcal/mol more stable",
  gate_tm = "allow COM 2 C below AS",
  gate_com_footprint = "dropout up to 2x tolerance",
  gate_wobble = "one weak G:T / A:C 3' mismatch",
  gate_as_footprint = "dropout up to 2x tolerance",
  gate_duplex_amplicon = "allow up to 200 bp",
  gate_duplex_selfcomp = "one extra self-complementary base",
  gate_caps_amplicon = "widen the amplicon window",
  gate_caps_resolvable = "allow a 30 bp band gap",
  gate_caps_tm = "widen the Tm window"
)

# Why a gate is hard (never lowered).
KASP_GATE_HARD_LABEL <- c(
  gate_discrimination = "no allele discrimination = no assay",
  gate_as_window = "an unlinked variant in the AS 3' window destroys specificity",
  gate_duplex_targets = "the two loci must match the chosen readout",
  gate_duplex_tm = "all primers must share the extension window",
  gate_duplex_complement = "primer cross-talk defeats multiplexing",
  gate_asq_cassette = "one cassette per target",
  gate_asq_readout = "no allele resolution in real time = no assay",
  gate_caps_site = "no discriminating enzyme = no CAPS marker"
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
      # (discrimination and the AS window are fixed at the AS stage, so they are
      # not re-ranked here)
      tier1 = .data$m_structure,
      tier2 = pmin(.data$m_tm_com, .data$m_tm_tgt),
      tier3 = .data$m_amplicon,
      tier4 = .data$m_footprint
    )
  list(
    data = d,
    gates = KASP_GATES_PAIR,
    tiers = names(KASP_TIERS_PAIR),
    # tie-breakers, ascending (smaller is better)
    ties = KASP_TIES_PAIR,
    # the margin that measures each gate, so a relaxed gate can be minimised first
    gate_margin = c(
      gate_discrimination = "m_disc", gate_structure = "m_structure",
      gate_tm = "m_tm_com", gate_amplicon = "m_amplicon",
      gate_com_footprint = "m_footprint", gate_as_window = "m_window"
    ),
    # how far each gate may be lowered (never -Inf).  Footprint tolerance scales
    # with the accepted dropout.
    gate_floor = {
      fl <- KASP_GATE_FLOOR_PAIR
      fl["gate_com_footprint"] <- if (use_drop) -2 * tol else -3
      fl
    }
  )
}

#' Grade candidate AS primers: hard gates + tiered margins
#'
#' The AS-side twin of `grade_assay_pairs()`: the same vocabulary (gates for the
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
  # the whole-footprint rule: an unlinked variant under the AS primer is a sin
  # only when it disrupts annealing, so gate on the expected dropout (position-
  # and identity-aware) - exactly like the COM footprint - and fall back to the
  # raw count only when DECIPHER is unavailable.  A linked variant is already
  # excluded from the count (as_exclude), so linkage is honoured for free.
  use_foot_drop <- "as_foot_dropout" %in% names(as_cand) && any(!is.na(as_cand$as_foot_dropout))
  foot <- if (use_foot_drop) {
    as.numeric(as_cand$as_foot_dropout)
  } else {
    as.numeric(as_cand$w_var)
  }
  gate_foot <- if (use_foot_drop) {
    is.na(as_cand$as_foot_dropout) | as_cand$as_foot_dropout <= tol
  } else {
    as_cand$n_var == 0L
  }
  d <- as_cand |>
    dplyr::mutate(
      as_row = dplyr::row_number(),
      m_disc = .data$discrimination - min_disc,
      m_tm = -abs(.data$tm - AS_TM_TARGET),
      m_stability = -as.numeric(.data$dg_3p),
      m_footprint = -foot,
      m_window = win_margin,
      m_wobble = -as.numeric(.data$wobble),
      gate_discrimination = .data$discrimination >= min_disc,
      gate_wobble = is.na(.data$wobble) | !.data$wobble,
      gate_as_footprint = gate_foot,
      gate_as_window = gate_win,
      tier1 = .data$m_disc,
      tier2 = .data$m_tm,
      tier3 = .data$m_stability,
      tier4 = .data$m_footprint
    )
  list(
    data = d,
    gates = KASP_GATES_AS,
    tiers = names(KASP_TIERS_AS),
    ties = KASP_TIES_AS,
    gate_margin = c(
      gate_wobble = "m_wobble", gate_discrimination = "m_disc",
      gate_as_footprint = "m_footprint", gate_as_window = "m_window"
    ),
    # footprint tolerance scales with the accepted dropout, like the pair stage
    gate_floor = {
      fl <- KASP_GATE_FLOOR_AS
      fl["gate_as_footprint"] <- if (use_foot_drop) -2 * tol else -1
      fl
    }
  )
}

#' Order rows lexicographically: `desc` columns high-to-low, then `asc` low-to-high
#' @keywords internal
#' @noRd
order_by_margins <- function(data, desc, asc = character(0)) {
  dplyr::arrange(
    data,
    dplyr::across(dplyr::all_of(intersect(desc, names(data))), dplyr::desc),
    dplyr::across(dplyr::all_of(intersect(asc, names(data))))
  )
}

#' Sort candidates lexicographically on tiers then ties (no gating)
#'
#' Used to build shortlists (e.g. which placements to seed-screen) with the same
#' ordering the final selection uses.
#' @keywords internal
#' @noRd
rank_candidates <- function(data, tiers, ties) order_by_margins(data, tiers, ties)

#' Pick the best candidate, lowering gate thresholds only as far as necessary
#'
#' Hard gates (floor 0) are never lowered: if nothing clears them, the design is
#' **doomed** (`best = NULL`).  Soft gates are relaxed in `relax_order` (least
#' grave first) and only **as far as the best achievable margin** - never to
#' `-Inf`, and never below the gate's floor.  A lowered gate stays active at the
#' new bound; among the survivors the choice is still the plain lexicographic tier
#' ranking (the relaxation does not override it).  Reports which thresholds were
#' lowered and which the winner still sits below.
#' @param data the graded tibble (see `grade_assay_pairs()`).
#' @param gates gate column names (all must be TRUE to pass).
#' @param relax_order soft gate names in the order they may be lowered.
#' @param tiers margin columns, ranked high-to-low, in order.
#' @param ties column names for ascending tie-breaks.
#' @param gate_margin named map `gate -> margin column`.
#' @param gate_floor named map `gate -> floor` (margin units; 0 = hard).
#' @return list with `best` (one row, or `NULL`), `relaxed` and `violated`.
#' @keywords internal
#' @noRd
select_with_ladder <- function(data, gates, relax_order, tiers, ties,
                               gate_margin = NULL, gate_floor = NULL) {
  if (!nrow(data)) {
    return(list(best = NULL, relaxed = character(0), violated = character(0)))
  }
  floors <- stats::setNames(rep(0, length(gates)), gates)
  floors[names(gate_floor)] <- gate_floor
  bars <- stats::setNames(rep(0, length(gates)), gates)
  passes <- function() {
    ok <- rep(TRUE, nrow(data))
    for (g in gates) {
      m <- data[[gate_margin[[g]]]]
      ok <- ok & (is.na(m) | m >= bars[[g]]) # NA (unknown) does not block
    }
    ok
  }
  # relax soft gates, least grave first, each only to the best achievable margin
  for (g in relax_order) {
    if (any(passes())) break
    if (floors[[g]] >= 0) next # hard
    m <- data[[gate_margin[[g]]]]
    best_achievable <- if (any(!is.na(m))) max(m, na.rm = TRUE) else -Inf
    if (is.finite(best_achievable) && best_achievable < 0) {
      bars[[g]] <- max(floors[[g]], best_achievable)
    }
  }
  ok <- passes()
  lowered <- names(bars)[bars < 0]
  if (!any(ok)) {
    return(list(best = NULL, relaxed = lowered, violated = lowered))
  }
  surv <- data[ok, , drop = FALSE]
  best <- order_by_margins(surv, tiers, intersect(ties, names(surv)))[1L, , drop = FALSE]
  violated <- gates[vapply(gates, function(g) {
    m <- best[[gate_margin[[g]]]]
    !is.na(m) && m < 0
  }, logical(1))]
  list(best = best, relaxed = lowered, violated = unique(violated))
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
