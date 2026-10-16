# =============================================================================
#  kasp_arms.R  --  tetra-primer ARMS-PCR (a gel-readout allele-specific assay)
# =============================================================================
#  Tetra-primer ARMS is the gel analogue of KASP: two outer primers give a
#  control amplicon, and two inner allele-specific primers - placed at the
#  variant, each on the opposite strand - give a REF band and an ALT band of
#  *different sizes*, so a plain agarose gel reads the genotype (AA / AB / BB).
#  Like CAPS it tolerates what KASP cannot (the readout is band presence and
#  size, not an efficiency ratio), but it shares KASP's 3'-terminal placement of
#  the variant.  The inner primers reuse [design_as_primers()]; the outer pair is
#  a plain flanking pair.
#
#  Design rules (researched, not invented):
#   * ARMS reads a 3'-terminal mismatch: an inner primer whose last base matches
#     one allele amplifies it and the other poorly (Newton et al. 1989, Nucleic
#     Acids Res 17:2503-2516, <https://doi.org/10.1093/nar/17.7.2503>).  No
#     discriminating inner primer = no ARMS assay, so that gate is hard.
#   * Tetra-primer ARMS adds the outer control pair, so AA / AB / BB are read as
#     band patterns, not presence/absence (Ye et al. 2001, Nucleic Acids Res
#     29:e88, <https://doi.org/10.1093/nar/29.17.e88>).
#   * The standard design guide (Medrano & de Oliveira 2014, Mol Biotechnol
#     56:599-610, <https://doi.org/10.1007/s12033-014-9734-4>) asks for the two
#     allele-specific products - and each of them from the control - to differ
#     enough to score on a plain gel (>= 20 bp; relaxed to a 10 bp PAGE job), a
#     PCR-able control amplicon (200-800 bp), and Tm-matched primers, since the
#     four-primer reaction runs at ONE annealing temperature.
#   * The two inner primers must sit on *opposite* strands so each competes with
#     the far outer primer; the control is always the largest product (the inner
#     primers are nested), which is what makes the control band recognisable.
# =============================================================================

# ---- the ARMS selection stage -------------------------------------------------
# Tetra-primer ARMS-PCR: an outer control pair plus two inner allele-specific
# primers on opposite strands, read as band sizes.  The two declared orders (C6):
# gates applied in the order below (grave / cheap first), then the ladder (least
# grave first: the control window is a preference, the band separation the point
# of the assay).
KASP_GATES_ARMS <- c("gate_arms_disc", "gate_arms_bands", "gate_arms_outer")
KASP_RELAX_ARMS <- c("gate_arms_outer", "gate_arms_bands")
KASP_TIERS_ARMS <- c(
  tier1 = "inner discrimination", tier2 = "band separation", tier3 = "tube Tm balance"
)

# The sizes the ARMS designer selects against.  A load-time constant so the gate
# limits can quote it (`KASP_LIMITS` lives in R/kasp_core.R, sourced *after* this
# file).  See the header for the sources; `band_gap` is the accepted closest-band
# difference, `control` the PCR-able control window, `tm_spread` the one-Tm bound
# over all four primers.
KASP_ARMS_WINDOW <- list(band_gap = 20L, control = c(200L, 800L), tm_spread = 5L)
KASP_ARMS_GAP_FLOOR <- 10L # bp: the band gap may be relaxed to a PAGE-resolvable 10 bp

#' The window the ARMS designer selects against (and the gate limits display)
#' @keywords internal
#' @noRd
arms_limits <- function() {
  list(
    band_gap = as.numeric(KASP_ARMS_WINDOW$band_gap),
    min_disc = KASP_LIMITS$discrimination,
    control = as.numeric(KASP_ARMS_WINDOW$control),
    tm_spread = as.numeric(KASP_ARMS_WINDOW$tm_spread)
  )
}

# The ladder's floors (margin units; 0 = hard) - C8.  `gate_arms_disc` is hard
# (no allele-specific inner primer = no ARMS assay).  The band gap is soft but
# *bounded*: it may be relaxed only to a 10 bp PAGE-resolvable difference - the
# floor no longer equals the 20 bp threshold (the audit's "drawn soft, can never
# relax" bug).  The control window widens by 50 bp at each end.
KASP_GATE_FLOOR_ARMS <- c(
  gate_arms_disc = 0,
  gate_arms_bands = KASP_ARMS_GAP_FLOOR - KASP_ARMS_WINDOW$band_gap, # 20 -> 10 bp
  gate_arms_outer = -50 # [200, 800] -> [150, 850]
)

# This scheme's entries in the shared gate-label lookups (assembled in
# R/kasp_select.R), kept beside the gates they name.
KASP_GATE_LABEL_ARMS <- c(
  gate_arms_disc = "an inner primer does not discriminate the allele",
  gate_arms_bands = "the diagnostic bands do not resolve on a gel",
  gate_arms_outer = "the outer control amplicon is outside the window"
)
KASP_GATE_SHORT_ARMS <- c(
  gate_arms_disc = "inner specificity", gate_arms_bands = "band gap",
  gate_arms_outer = "outer control"
)
KASP_GATE_LIMIT_ARMS <- c(
  gate_arms_disc = "inner 3' base matches the allele (discrimination >= min_disc)",
  gate_arms_bands = glue::glue(
    "REF/ALT and each vs the control >= {KASP_ARMS_WINDOW$band_gap} bp apart"
  ),
  gate_arms_outer = glue::glue(
    "outer control within [{KASP_ARMS_WINDOW$control[1L]}, {KASP_ARMS_WINDOW$control[2L]}] bp"
  )
)
KASP_GATE_FLOOR_LABEL_ARMS <- c(
  gate_arms_bands = "allow a 10 bp band gap (PAGE)",
  gate_arms_outer = "widen the control window"
)
KASP_GATE_HARD_LABEL_ARMS <- c(
  gate_arms_disc = "no allele-specific inner primer = no ARMS assay"
)

#' Candidate outer (control) primer pairs across the control window
#'
#' `design_flanking_pair()` returns the single most Tm-balanced pair, which would
#' make the control-window gate vacuous.  This keeps a *spread* of products (best
#' balanced first, ~`spread` bp apart) over a window widened by the gate's floor,
#' so the gate has something to choose between and to reject.
#' @keywords internal
#' @noRd
arms_outer_pairs <- function(seq, snp_i, limits, n_outer = 6L, spread = 50L) {
  n <- nchar(as.character(seq))
  win <- c(limits$control[1L] - spread, limits$control[2L] + spread)
  fwd <- primer_candidates(seq, seq_len(snp_i - 1L), 18:25, "+", KASP_LIMITS$com_tm, KASP_LIMITS$gc) |>
    dplyr::filter(.data$end < snp_i)
  rev <- primer_candidates(seq, seq.int(snp_i + 1L, n), 18:25, "-", KASP_LIMITS$com_tm, KASP_LIMITS$gc)
  pairs <- dplyr::cross_join(fwd, rev, suffix = c("_f", "_r")) |>
    dplyr::mutate(
      product = .data$end_r - .data$start_f + 1L,
      dtm = abs(.data$tm_f - .data$tm_r)
    ) |>
    dplyr::filter(dplyr::between(.data$product, win[1L], win[2L])) |>
    dplyr::arrange(.data$dtm, .data$product)
  # greedy spread: the best-balanced pair, then the next one >= `spread` bp away
  keep <- integer(0)
  for (k in seq_len(nrow(pairs))) {
    if (all(abs(pairs$product[k] - pairs$product[keep]) >= spread)) keep <- c(keep, k)
    if (length(keep) >= n_outer) break
  }
  pairs <- pairs[keep, , drop = FALSE]
  tibble::tibble(
    forward = pairs$primer_f, reverse = pairs$primer_r,
    start = pairs$start_f, end = pairs$end_r, product = pairs$product,
    tm_f = pairs$tm_f, tm_r = pairs$tm_r, dtm = pairs$dtm
  )
}

#' Design a tetra-primer ARMS-PCR marker at a variant
#'
#' The gel-readout allele-specific assay: an outer control pair plus two inner
#' allele-specific primers (one per allele, on opposite strands) whose products
#' differ in size, so AA / AB / BB resolve on an agarose gel.  The inner primers
#' reuse [design_as_primers()] (the same 3'-terminal placement, plus its optional
#' stabilising mismatch), the outer pair is a plain flanking pair.  Single-base
#' substitutions only.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank window half-width for the outer pair.
#' @param min_band_gap minimal size difference between the two diagnostic bands.
#' @param ... passed to [design_as_primers()] (e.g. `stab_ats`, `BPPARAM`).
#' @return a `kasp_arms` list (`snp`, `outer`, `inner`, `bands`, `amplicon`,
#'   `compatibility` - the joint dimer/hairpin screen over all four primers), or
#'   `NULL` (with a message) when a pair cannot be built.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_arms_marker(genome, kasp_example("SNP")[1, ], flank = 300)
#' }
#' @family KASP ARMS
#' @export
design_arms_marker <- function(genome, snp, flank = 300L, min_band_gap = 20L, ...) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    cli::cli_inform(gettext("tetra-primer ARMS is implemented for single-base substitutions only."))
    return(NULL)
  }
  limits <- arms_limits()
  limits$band_gap <- as.numeric(min_band_gap)
  floors <- KASP_GATE_FLOOR_ARMS
  floors[["gate_arms_bands"]] <- KASP_ARMS_GAP_FLOOR - limits$band_gap
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  outer <- arms_outer_pairs(tmpl$seq, tmpl$snp_i, limits)
  if (!nrow(outer)) {
    cli::cli_inform(gettext("No outer primer pair passed the Tm/GC filters."))
    return(NULL)
  }
  inner <- design_as_primers(tmpl$seq, tmpl$snp_i, ref, alt,
    strands = c("+", "-"), gc_range = c(0.35, 0.70),
    tm_range = KASP_LIMITS$as_tm, ...
  )
  if (!nrow(inner)) {
    cli::cli_inform(gettext("No inner allele-specific primer passed the GC/Tm filter."))
    return(NULL)
  }
  # one inner primer per strand: REF on "+" (forward), ALT on "-" (reverse), so
  # each competes with the far outer primer as the method requires
  top <- function(st) {
    utils::head(
      dplyr::arrange(dplyr::filter(inner, .data$strand == st), dplyr::desc(.data$discrimination)), 8L
    )
  }
  fwd <- top("+")
  rev <- top("-")
  if (!nrow(fwd) || !nrow(rev)) {
    cli::cli_inform(gettext("An inner primer is missing on one strand."))
    return(NULL)
  }
  # REF band: inner-forward .. outer-reverse; ALT band: outer-forward .. inner-reverse
  cand <- tidyr::crossing(o = seq_len(nrow(outer)), i = seq_len(nrow(fwd)), j = seq_len(nrow(rev))) |>
    dplyr::mutate(
      ref_primer = fwd$ref_primer[.data$i], alt_primer = rev$alt_primer[.data$j],
      ref_start = fwd$tmpl_start[.data$i], ref_end = fwd$tmpl_end[.data$i],
      alt_start = rev$tmpl_start[.data$j], alt_end = rev$tmpl_end[.data$j],
      outer_forward = outer$forward[.data$o], outer_reverse = outer$reverse[.data$o],
      outer_start = outer$start[.data$o], outer_end = outer$end[.data$o],
      control = outer$product[.data$o],
      ref_size = .data$outer_end - .data$ref_start + 1L,
      alt_size = .data$alt_end - .data$outer_start + 1L,
      disc = pmin(fwd$discrimination[.data$i], rev$discrimination[.data$j]),
      tube_spread = pmax(fwd$tm[.data$i], rev$tm[.data$j], outer$tm_f[.data$o], outer$tm_r[.data$o]) -
        pmin(fwd$tm[.data$i], rev$tm[.data$j], outer$tm_f[.data$o], outer$tm_r[.data$o])
    ) |>
    dplyr::mutate(
      m_disc = .data$disc - limits$min_disc,
      # the two diagnostic bands must resolve from *each other* and from the
      # control, so the gate is the closer of the two differences
      m_bands = pmin(
        abs(.data$ref_size - .data$alt_size),
        .data$control - pmax(.data$ref_size, .data$alt_size)
      ) - limits$band_gap,
      m_outer = pmin(.data$control - limits$control[1L], limits$control[2L] - .data$control),
      tier1 = .data$disc - limits$min_disc, tier2 = .data$m_bands,
      tier3 = -.data$tube_spread
    )
  pick <- select_with_ladder(
    cand,
    gates = KASP_GATES_ARMS, relax_order = KASP_RELAX_ARMS,
    tiers = names(KASP_TIERS_ARMS), ties = character(0),
    gate_margin = c(
      gate_arms_disc = "m_disc", gate_arms_bands = "m_bands",
      gate_arms_outer = "m_outer"
    ),
    gate_floor = floors
  )
  best <- pick$best
  if (is.null(best)) {
    cli::cli_inform(gettext("No tetra-primer set cleared even the relaxed constraints."))
    return(NULL)
  }
  off <- snp$pos - tmpl$snp_i
  out <- list(
    snp = snp,
    outer = list(
      forward = best$outer_forward, reverse = best$outer_reverse,
      start = best$outer_start + off, end = best$outer_end + off, product = best$control
    ),
    inner = list(
      ref = list(
        primer = best$ref_primer, start = best$ref_start + off,
        end = best$ref_end + off, strand = "+"
      ),
      alt = list(
        primer = best$alt_primer, start = best$alt_start + off,
        end = best$alt_end + off, strand = "-"
      )
    ),
    bands = tibble::tibble(
      band = c("control", "REF", "ALT"),
      size = c(best$control, best$ref_size, best$alt_size)
    ),
    amplicon = tibble::tibble(start = best$outer_start + off, end = best$outer_end + off),
    selection = list(relaxed = pick$relaxed, violated = pick$violated)
  )
  # all four primers share the tube: screen the joint set, not just the pair
  out$compatibility <- tube_compatibility(c(
    outer_forward = out$outer$forward, outer_reverse = out$outer$reverse,
    AS_ref = out$inner$ref$primer, AS_alt = out$inner$alt$primer
  ))
  warn_relaxed(pick$violated)
  structure(out, class = "kasp_arms")
}

#' @export
print.kasp_arms <- function(x, ...) {
  cat("<kasp_arms>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt, "\n", sep = "")
  cat("  outer:", x$outer$product, "bp\n")
  cat("  REF band:", x$bands$size[x$bands$band == "REF"], "bp (inner F:", x$inner$ref$primer, ")\n")
  cat("  ALT band:", x$bands$size[x$bands$band == "ALT"], "bp (inner R:", x$inner$alt$primer, ")\n")
  invisible(x)
}

#' Draw the expected tetra-primer ARMS gel
#'
#' Two lanes (REF and ALT homozygotes) with the control band plus the
#' allele-specific diagnostic band - the AA / AB / BB readout.
#' @param marker a `kasp_arms` object (see [design_arms_marker()]).
#' @return a `ggplot`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' plot_arms_gel(design_arms_marker(genome, kasp_example("SNP")[1, ], flank = 300))
#' }
#' @family KASP ARMS
#' @export
plot_arms_gel <- function(marker) {
  rlang::check_installed("ggplot2", reason = gettext("to draw the ARMS gel"))
  if (!inherits(marker, "kasp_arms")) {
    cli::cli_abort(gettext("`marker` must be a kasp_arms object (see design_arms_marker())."))
  }
  ctrl <- marker$bands$size[marker$bands$band == "control"]
  ref <- marker$bands$size[marker$bands$band == "REF"]
  alt <- marker$bands$size[marker$bands$band == "ALT"]
  bands <- tibble::tibble(
    allele = c("REF", "REF", "ALT", "ALT"),
    size = c(ctrl, ref, ctrl, alt),
    diagnostic = c(FALSE, TRUE, FALSE, TRUE)
  ) |>
    dplyr::mutate(
      x = as.numeric(factor(.data$allele, levels = c("REF", "ALT"))),
      col = dplyr::case_when(
        !.data$diagnostic ~ "grey65",
        .data$allele == "REF" ~ KASP_STATUS[["neutral"]],
        TRUE ~ KASP_STATUS[["broken"]]
      )
    )
  gel_plot(
    bands, c("REF/REF", "ALT/ALT"), "band size (bp)",
    glue::glue("ARMS gel: {marker$snp$chrom}:{marker$snp$pos}"),
    glue::glue("control {ctrl} bp | REF {ref} bp | ALT {alt} bp")
  )
}
