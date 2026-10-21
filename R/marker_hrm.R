# =============================================================================
#  marker_hrm.R  --  HRM (high-resolution melting) as a marker scheme
# =============================================================================
#  HRM reads the genotype from the *melting curve of a short amplicon*: the
#  alleles differ in how the duplex melts, and heterozygotes show a heteroduplex
#  shape.  It is the scheme for the loci KASP cannot take - a variant in the AS
#  3' window, a poly-A junction - because it anchors the variant in the
#  *amplicon*, not in a primer, and *tolerates neighbours* (the curve separates
#  them).  So there is no allele-specific primer: one flanking pair straddles the
#  variant, the amplicon is kept SHORT and the pair is Tm-balanced.
#  The melt curve itself is not simulated (that needs uMELT/PolyMelt-class
#  models), so this returns a well-formed amplicon, not a predicted profile.
#
#  Design rules (researched, not invented):
#   - The SHORT amplicon is HRM's core requirement, not a preference.  SNP
#     genotyping by melting works on *small* amplicons - "Genotyping of
#     single-nucleotide polymorphisms by high-resolution melting of small
#     amplicons", Liew, Pryor, Palais, et al., Clin Chem 2004, 50(7):1156-64 -
#     because a short product melts as (nearly) one domain, so the single base
#     change shifts a sharp, reproducible transition; a long product melts in
#     several domains and the variant's contribution is diluted away.  Reed &
#     Wittwer, Clin Chem 2004, 50(10):1748-54, show the same amplicon-length
#     dependence for scanning sensitivity, and Wittwer et al., Clin Chem 2003,
#     49(6):853-60, established the LCGreen closed-tube method.
#     => `gate_hrm_amplicon` is therefore HARD, and its window is exactly the
#        marker-agnostic amplicon check's window (`KASP_LIMITS$amplicon`), so a
#        marker the gate admits is never rejected by the acceptance check (C5).
#   - The flanking pair must be Tm-balanced so both alleles are amplified
#     equally: an unbalanced pair biases one allele's template and moves the
#     curve for reasons that have nothing to do with the variant.  This is a
#     real requirement but a soft one - a couple of degrees can be absorbed by
#     the annealing temperature - so it is the ladder's rung, not a hard block.
#   * The variant must sit in the melt domain and under NEITHER primer (a primer
#     whose 3' end abuts the variant starts to read the allele itself - that is
#     allele-specific PCR, not HRM).  That is enforced *structurally* by the
#     designer's flanking geometry (a clearance on each side), not by a gate.
#   * Resolvability is a DIFFERENCE of melting peaks, not a property of a SNP:
#     `hrm_resolvability()` grades the bimorphism against the source's own
#     measured dTm, and `hrm_resolvability_row()` puts that into the evaluation
#     so an HRM marker never reports "designable" without saying whether its
#     class has a curve to read (decision:
#     data-raw/proposed/hrm_resolvability.decision.md).
# =============================================================================

# ---- the HRM selection stage --------------------------------------------------
# Two gates, in application order (grave / cheap first): the amplicon (HRM's
# core, hard) then the pair's Tm balance (soft).  The ladder has one rung, the
# Tm window, widened least-grave-first; the amplicon window is never relaxed.
# The two declared orders are therefore not the same thing (C6).
KASP_GATES_HRM <- c("gate_hrm_amplicon", "gate_hrm_tm")
KASP_RELAX_HRM <- c("gate_hrm_tm")
# Tiers rank what is left once the hard amplicon gate is passed: the variant as
# central in the melt domain as possible (so neither primer is close and the
# curve difference is maximal), then the pair's Tm balance.  Tier 1 deliberately
# does *not* restate a gate (C14).
KASP_TIERS_HRM <- c(tier1 = "variant centrality (melt domain)", tier2 = "Tm balance")

# ---- resolvability (decided, Q2 / hrm_resolvability) --------------------------
# HRM's outcome is a DIFFERENCE of melting peaks, not a yes/no property of a SNP.
# `HRM_DTM_RESOLVABLE` is the smallest homozygote-vs-homozygote dTm called
# resolvable, in C.  Source: JSBBS 2018 (10.1270/jsbbs.18048) states "a dTm value
# > 0.5 C guaranteed accurate genotyping" ON A LIGHT CYCLER 96 and adds that some
# class III/IV markers still genotyped at 0.2-0.3 C - so this is a conservative
# DEFAULT FOR THAT PLATFORM, not a physical constant.  The lookup is the paper's
# own per-bimorphism measured minimum (Supplemental Table 2/5); any pair absent
# from it has no measurement and returns `NA` ("unknown"), never a guess.
HRM_DTM_RESOLVABLE <- 0.5

# Measured dTm by bimorphism: `min` = smallest dTm still called polymorphic on
# that substitution class (so a class can fail the 0.5 C default and still be
# resolvable), `mean` = the typical value.  Asserted against the extract by
# `tests/testthat/test-hrm-jsbbs2018.R`.
HRM_BIMORPHISM_DTM <- list(
  "A/C" = c(min = 0.5139, mean = 0.5886),
  "A/G" = c(min = 0.2447, mean = 0.6486),
  "C/G" = c(min = 0.1976, mean = 0.1976),
  "C/T" = c(min = 0.3915, mean = 0.6922),
  "G/T" = c(min = 0.3671, mean = 0.5187)
)

#' Will HRM resolve the two homozygotes of this SNP?
#'
#' Answers "is there a dTm to read", from the bimorphism of the pair - the one
#' input every fetched source measured it for.  `resolvable` when the class mean
#' clears the threshold (the class is what a designer can expect); `marginal`
#' when the mean falls short but the class still produced markers above the
#' threshold (it can work, it is not dependable); `below_threshold` when no
#' member of the class reached it.  A bimorphism with no landed measurement is
#' `unknown`, because no fetched source reports a dTm for it.  A multi-locus
#' panel is NOT answered per SNP: the outcome there is a partition of control
#' genotypes (3 loci left 14/36 controls in 7 indistinguishable pairs), so the
#' caller must run the panel, not ask.
#' @param ref,alt the two alleles (a base or an IUPAC code; a homozygote pair
#'   such as `"A"`/`"G"`).
#' @param dg3bp `|dG3bp|` of the 3-bp window, or `NA`; recorded, not thresholded
#'   (the free-energy filter is a filter, never a guarantee).
#' @param product amplicon length in bp, or `NA`.
#' @param threshold the pass threshold (default [HRM_DTM_RESOLVABLE]).
#' @return a classed list: `bimorphism`, `status`
#'   (`resolvable`/`marginal`/`below_threshold`/`unknown`), `dtm_min`, `dtm_mean`,
#'   `threshold`, `product`, `dg3bp`, `evidence`.
#' @keywords internal
#' @noRd
hrm_resolvability <- function(ref, alt, dg3bp = NA_real_, product = NA_real_,
                              threshold = HRM_DTM_RESOLVABLE) {
  pair <- sort(unique(toupper(c(as.character(ref)[1L], as.character(alt)[1L]))))
  # a genotype, not a set of alleles: two distinct bases make a homozygote
  # comparison; anything else (identical bases, an IUPAC code, an indel) is a
  # question this lookup has no measurement for
  pair <- if (length(pair) == 2L) sort(pair) else character(0)
  bimo <- if (length(pair) == 2L && all(pair %in% c("A", "C", "G", "T"))) {
    paste(pair, collapse = "/")
  } else {
    NA_character_
  }
  known <- !is.na(bimo) && bimo %in% names(HRM_BIMORPHISM_DTM)
  dtm <- if (known) HRM_BIMORPHISM_DTM[[bimo]] else c(min = NA_real_, mean = NA_real_)
  status <- if (!known) {
    "unknown"
  } else if (dtm[["mean"]] >= threshold) {
    "resolvable"
  } else if (dtm[["min"]] >= threshold) {
    "marginal"
  } else {
    "below_threshold"
  }
  structure(
    list(
      bimorphism = bimo,
      status = status,
      dtm_min = unname(dtm[["min"]]),
      dtm_mean = unname(dtm[["mean"]]),
      threshold = threshold,
      product = as.numeric(product),
      dg3bp = as.numeric(dg3bp),
      evidence = if (known) "measured (LightCycler 96, JSBBS 2018)" else "none landed"
    ),
    class = "marker_hrm_resolvability"
  )
}

#' Print an HRM resolvability answer
#' @param x a [hrm_resolvability()] result.
#' @param ... unused.
#' @return `x`, invisibly.
#' @keywords internal
#' @noRd
print.marker_hrm_resolvability <- function(x, ...) {
  cli::cli_inform(c(
    "i" = "HRM {x$bimorphism}: {x$status}",
    "*" = "dTm min {round(x$dtm_min, 2)} C, mean {round(x$dtm_mean, 2)} C, threshold {x$threshold} C"
  ))
  invisible(x)
}

#' The HRM resolvability row of an evaluation table (decided criterion)
#'
#' Every HRM marker is graded on whether its bimorphism has a melting difference
#' to read.  `resolvable` passes, `marginal` does not (`unknown` in the table,
#' since the class mean - not this SNP - is what falls short), `below_threshold`
#' and `unknown` report `unknown` too: no fetched source measures a dTm for the
#' pair, so nothing is known and nothing is claimed.  No design is failed by it
#' - the scheme has no per-SNP dTm to fail it with - and the reason is carried in
#' `detail` so a report can say which way it leans.
#' @param marker an HRM marker.
#' @return a one-row evaluation tibble (empty when the target alleles are
#'   unknown).
#' @keywords internal
#' @noRd
hrm_resolvability_row <- function(marker) {
  empty <- marker_evaluation_table(
    stage = character(), criterion = character(), observed = numeric(),
    limit = character(), pass = logical(), level = character()
  )
  snp <- marker$snp
  if (is.null(snp) || !nrow(snp) || anyNA(c(snp$ref[1L], snp$alt[1L]))) {
    return(empty)
  }
  r <- hrm_resolvability(snp$ref[1L], snp$alt[1L])
  if (identical(r$status, "unknown")) {
    return(empty)
  }
  amp <- marker$amplicon
  product <- if (!is.null(amp) && "product" %in% names(amp)) as.numeric(amp$product[1L]) else NA_real_
  marker_evaluation_table(
    stage = "check",
    criterion = "family_hrm_resolvable",
    observed = r$dtm_mean,
    limit = glue::glue(">= {r$threshold} C class mean dTm"),
    pass = identical(r$status, "resolvable"),
    detail = glue::glue(
      "{r$bimorphism}: {r$status}; dTm min {round(r$dtm_min, 2)} C, ",
      "mean {round(r$dtm_mean, 2)} C; {r$evidence}",
      if (is.na(product)) "" else ", amplicon {as.integer(product)} bp"
    ),
    level = "family"
  )
}

# How far the ladder may widen the flanking primers' Tm window (C, each side)
# when nothing clears the strict one: a bounded, reachable floor (C8).
HRM_TM_FLOOR <- 2
# Bases of melt domain kept between the variant and each primer's 3' end.  A
# primer 3' base on (or adjacent to) the variant reads the allele, so the pair
# must clear it on both sides - structural, not a gate.
HRM_VARIANT_CLEARANCE <- 1L

# This scheme's entries in the shared gate-label lookups (assembled in
# R/marker_select.R), kept beside the gates they name.
KASP_GATE_LABEL_HRM <- c(
  gate_hrm_amplicon = "the HRM amplicon is outside the short window",
  gate_hrm_tm = "a flanking primer Tm is outside the window"
)
KASP_GATE_SHORT_HRM <- c(
  gate_hrm_amplicon = "amplicon length (HRM)", gate_hrm_tm = "primer Tm (HRM)"
)
KASP_GATE_LIMIT_HRM <- c(
  gate_hrm_amplicon = glue::glue(
    "amplicon within [{KASP_LIMITS$amplicon[1L]}, {KASP_LIMITS$amplicon[2L]}] bp"
  ),
  gate_hrm_tm = glue::glue(
    "both flanking primers Tm within [{KASP_LIMITS$com_tm[1L]}, {KASP_LIMITS$com_tm[2L]}] C"
  )
)
KASP_GATE_FLOOR_LABEL_HRM <- c(
  gate_hrm_tm = glue::glue("widen the Tm window by {HRM_TM_FLOOR} C")
)
KASP_GATE_HARD_LABEL_HRM <- c(
  gate_hrm_amplicon = "a long amplicon melts in several domains - the curve no longer resolves the variant"
)

#' The best HRM flanking pair: shortest melt domain, then the most balanced Tm
#'
#' The HRM-specific sibling of `design_flanking_pair()`: the pair must sit
#' strictly outside `[lo, hi]`, the product must fall in `[min_product,
#' max_product]` (the hard amplicon gate), and the survivors are ranked by the
#' two declared tiers - the variant as central in the amplicon as possible, then
#' the smallest Tm difference - rather than by Tm balance alone.
#' @param seq the template (plus strand).
#' @param lo,hi plus-strand bounds the primers must stay outside (the variant
#'   plus its clearance).
#' @param lengths primer lengths to try.
#' @param min_product,max_product accepted amplicon window.
#' @param tm_range,gc_range accepted primer windows.
#' @return a pair list (see `pair_to_list()`), or `NULL`.
#' @keywords internal
#' @noRd
hrm_flanking_pair <- function(seq, lo, hi, lengths, min_product, max_product,
                              tm_range, gc_range) {
  n <- nchar(as.character(seq))
  if (lo <= 1L || hi >= n) {
    return(NULL)
  }
  fwd <- primer_candidates(seq, seq_len(lo - 1L), lengths, "+", tm_range, gc_range) |>
    dplyr::filter(.data$end < lo)
  rev <- primer_candidates(seq, seq.int(hi + 1L, n), lengths, "-", tm_range, gc_range)
  if (!nrow(fwd) || !nrow(rev)) {
    return(NULL)
  }
  centre <- (lo + hi) / 2 # the variant's position for the symmetric clearance
  pairs <- dplyr::cross_join(fwd, rev, suffix = c("_f", "_r")) |>
    dplyr::mutate(
      product = .data$end_r - .data$start_f + 1L,
      centre_dist = abs((.data$start_f + .data$end_r) / 2 - centre),
      dtm = abs(.data$tm_f - .data$tm_r)
    ) |>
    dplyr::filter(dplyr::between(.data$product, min_product, max_product))
  if (!nrow(pairs)) {
    return(NULL)
  }
  # the declared tiers, in order: variant centrality, then Tm balance
  pair_to_list(pairs[order(pairs$centre_dist, pairs$dtm), ][1L, ])
}

#' Design an HRM marker at a variant
#'
#' The melting-curve analogue of [assemble_kasp_assay()]: a single flanking pair
#' straddles the variant, which sits in the melt domain and under neither primer.
#' The design optimises the two things that govern HRM - a **short** amplicon and
#' a **Tm-balanced** pair - and, because HRM tolerates neighbouring variants,
#' imposes no clean 3' window.  This is why HRM is the fallback for the loci KASP
#' dooms (a variant in the AS 3' window, an AT-rich / poly-A junction).
#'
#' The short amplicon is the *hard* requirement: it is the marker-agnostic
#' amplicon window (`KASP_LIMITS$amplicon`), so a marker this gate admits always
#' passes the acceptance check on the same criterion.  The pair's Tm balance is
#' the soft one: when no pair clears the strict Tm window the ladder widens it by
#' `HRM_TM_FLOOR` degrees and the marker records `selection$relaxed`.
#'
#' It does **not** simulate the melting curve (that needs uMELT/PolyMelt-class
#' models), so it returns a well-formed amplicon, not a predicted melt profile.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param min_product,max_product the accepted amplicon window (HRM wants short).
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows.
#' @return a `marker_hrm` list (`scheme`, `snp`, `primers`, `amplicon`,
#'   `selection`), or `NULL` (with a message) when no flanking pair passes the
#'   filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_hrm_marker(genome, marker_example("SNP")[1, ])
#' }
#' @family KASP HRM
#' @export
design_hrm_marker <- function(genome, snp, flank = 150L,
                              min_product = KASP_LIMITS$amplicon[1L],
                              max_product = KASP_LIMITS$amplicon[2L], lengths = 18:25,
                              tm_range = KASP_LIMITS$com_tm, gc_range = KASP_LIMITS$gc) {
  ref <- toupper(as.character(snp$ref)[1L])
  alt <- toupper(as.character(snp$alt)[1L])
  if (nchar(ref) == 1L && nchar(alt) == 1L && ref %in% c("A", "T") && alt %in% c("A", "T")) {
    cli::cli_warn(gettext(
      "A/T SNP: the HRM Tm shift is small (Liew 2004); use high-resolution melt and a positive control."
    ))
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  si <- tmpl$snp_i
  # the variant must sit in the melt *domain*, under neither primer: the pair
  # flanks it with `HRM_VARIANT_CLEARANCE` bases of room on each side, so neither
  # 3' end can read the allele (that would be allele-specific PCR, not HRM)
  lo <- si - HRM_VARIANT_CLEARANCE
  hi <- si + HRM_VARIANT_CLEARANCE
  # gate_hrm_amplicon is hard: the short, single-domain amplicon is the whole
  # point, and `hrm_flanking_pair()` applies that window as a filter - never a
  # rung.  gate_hrm_tm is the soft ladder: the strict window first, then the
  # widened one (least grave first).
  primers <- hrm_flanking_pair(
    tmpl$seq, lo, hi, lengths, min_product, max_product, tm_range, gc_range
  )
  relaxed <- character(0)
  if (is.null(primers)) {
    widened <- tm_range + c(-HRM_TM_FLOOR, HRM_TM_FLOOR)
    primers <- hrm_flanking_pair(
      tmpl$seq, lo, hi, lengths, min_product, max_product, widened, gc_range
    )
    if (!is.null(primers)) relaxed <- "gate_hrm_tm"
  }
  if (is.null(primers)) {
    cli::cli_inform(gettext("No short flanking primer pair in the HRM amplicon window passed the Tm/GC filters."))
    return(NULL)
  }
  out <- list(
    scheme = "hrm",
    snp = snp, primers = primers,
    amplicon = tibble::tibble(
      start = snp$pos + primers$start - si,
      end = snp$pos + primers$end - si
    ),
    # the gate statuses are the record the graph and a report read
    selection = list(relaxed = relaxed, violated = character(0))
  )
  structure(out, class = "marker_hrm")
}

#' @export
print.marker_hrm <- function(x, ...) {
  cat(
    "<marker_hrm>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  amplicon:", x$primers$product, "bp (melt domain)\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  if (length(x$selection$relaxed)) {
    cat("  relaxed:", paste(x$selection$relaxed, collapse = ", "), "\n")
  }
  invisible(x)
}

#' Map an HRM amplicon
#'
#' The HRM analogue of a gel: the short amplicon as a bar, the two flanking
#' primers on it and the variant at the centre of the melt domain.  There is no
#' melt curve to draw here (we do not simulate it), so this shows *where* the
#' amplicon sits rather than how it melts.
#' @param marker a `marker_hrm` marker (see [design_hrm_marker()]).
#' @return a `ggplot`.
#' @examples
#' \dontrun{
#' plot_hrm_amplicon(design_hrm_marker(genome, snp))
#' }
#' @family KASP HRM
#' @export
#' @method marker_map_spec marker_hrm
#' @export
#' @noRd
marker_map_spec.marker_hrm <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  map_spec(
    marker,
    features = map_pair_rows(marker, sp[["start"]], sp[["end"]]),
    sequences = list(forward = marker$primers$forward, reverse = marker$primers$reverse)
  )
}

plot_hrm_amplicon <- function(marker) {
  if (!inherits(marker, "marker_hrm")) {
    cli::cli_abort(gettext("`marker` must be a marker_hrm object (see design_hrm_marker())."))
  }
  amp <- marker$amplicon
  pr <- marker$primers
  primers <- tibble::tibble(
    x = c(amp$start, amp$end - nchar(pr$reverse) + 1L),
    xend = c(amp$start + nchar(pr$forward) - 1L, amp$end)
  )
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = tibble::tibble(x = amp$start, xend = amp$end),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = 1, yend = 1),
      linewidth = 3, colour = "grey80", lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = primers,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = 1.3, yend = 1.3),
      linewidth = 2.2, colour = "#57534E", lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = tibble::tibble(x = marker$snp$pos),
      ggplot2::aes(x = .data$x, xend = .data$x, y = 0.65, yend = 1.35),
      colour = "firebrick", linewidth = 0.5
    ) +
    ggplot2::labs(
      x = glue::glue("{marker$snp$chrom} position (bp)"), y = NULL,
      title = "HRM amplicon (melt domain)"
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
}
