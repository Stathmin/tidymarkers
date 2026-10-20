# =============================================================================
#  marker_indel_size.R  --  InDel size-polymorphism genotyping (gel / capillary)
# =============================================================================
#  An InDel size marker reads the genotype from the *length* of one PCR product:
#  a single flanking pair spans the insertion/deletion, so the REF and ALT
#  amplicons differ by exactly the indel length and AA / AB / BB resolve by
#  fragment size - on a plain agarose gel, or (for small differences) on a
#  capillary sequencer.  The two physical requirements are therefore (a) both
#  primers sit in the conserved flanks *outside* the whole indel, so each binds
#  both alleles and the product-size difference is exactly the indel length, and
#  (b) the REF/ALT size difference is *resolvable* by the readout.
#
#  Design rules, with the sources they come from:
#   * One flanking pair, both primers in the conserved flanks *outside* the whole
#     indel, so each primer binds both alleles and the product-size difference is
#     exactly the indel length.
#   * Amplicon 80-300 bp and primer Tm ~57-63 C - the windows used for
#     agarose-resolvable cucumber InDels (Adedze et al. 2021, Sci Rep 11:3872,
#     <https://doi.org/10.1038/s41598-021-83313-x>); a companion screen used
#     100-350 bp at Tm 55-60 C (Tian et al. 2026, Sci Rep 16:23139,
#     <https://doi.org/10.1038/s41598-026-52824-w>).
#   * Resolvability: conventional agarose resolves a size difference of only
#     ~5-10 bp, so small InDels (1-2 bp) need PAGE or capillary electrophoresis
#     (Wang et al. 2026, Int J Mol Sci 27:1406,
#     <https://doi.org/10.3390/ijms27031406>).  We default `min_gap = 5` bp for
#     the agarose readout and *relax the gap gate down to a 3 bp capillary
#     difference* rather than hard-rejecting: the gap is a soft gate with a
#     bounded, reachable floor (C8), so a 3-4 bp indel still designs, with a
#     warning naming the relaxed gate.
#   * Long InDels (>30 bp) are the agarose-friendly ones, and differences
#     >= 100 bp score unambiguously (Chowdhury et al. 2026, Plant Biotechnol
#     43:227, <https://doi.org/10.5511/plantbiotechnology.26.0420a>).
#  The primer scan also drops homopolymer runs / SSRs, which is the practical
#  "unique flank" screen.  The primer structure and the genome-wide uniqueness
#  gates are hard, evaluated by [design_indel_size_marker()] (C7): the structure
#  one on the shortlist (OligoArrayAux), the uniqueness one via the batched
#  3'-seed pre-screen and a full in-silico PCR on the winner.
# =============================================================================

# The InDel-size acceptance windows.  `min_gap` is the agarose-resolution floor
# (the default starting gate); `capillary_gap` is the lower bound the ladder may
# relax it to, because a capillary/PAGE readout resolves a smaller difference.
KASP_INDEL_SIZE_LIMITS <- list(
  amplicon = c(80L, 300L),
  tm = KASP_LIMITS$com_tm,
  min_gap = 5L,
  capillary_gap = 3L
)

# How far each *soft* InDel-size gate may be relaxed (bounded, never -Inf, C8):
# the amplicon window widens 50 bp each side and the Tm window 5 C each side (the
# gap floor is `min_gap - capillary_gap`, computed per call).
KASP_INDEL_SIZE_FLOOR <- list(amplicon = 50L, tm = 5)

# The InDel-size gates, in application order (grave / cheap first, C6): the
# marker's defining criterion first (is the size difference resolvable at all),
# then the pair's cheap windows, then the expensive structure and genome-wide
# uniqueness.  `structure` and `unique` are hard: a self-structured primer does
# not prime, and a non-unique pair amplifies the wrong locus.
KASP_GATES_INDEL_SIZE <- c(
  "gate_indel_size_gap", "gate_indel_size_amplicon", "gate_indel_size_tm",
  "gate_indel_size_structure", "gate_indel_size_unique"
)
# Relaxation order - least grave first (C6): a wider amplicon or Tm window is a
# mild trade-off; the resolvable-gap floor is the last resort (a 3 bp capillary
# difference).  The structure and uniqueness gates are hard.
KASP_RELAX_INDEL_SIZE <- c("gate_indel_size_amplicon", "gate_indel_size_tm", "gate_indel_size_gap")
KASP_TIERS_INDEL_SIZE <- c(tier1 = "Tm balance", tier2 = "amplicon length", tier3 = "primer 3' end")

KASP_GATE_SHORT_INDEL_SIZE <- c(
  gate_indel_size_gap = "resolvable gap",
  gate_indel_size_amplicon = "amplicon length (InDel)",
  gate_indel_size_tm = "primer Tm (InDel)",
  gate_indel_size_structure = "primer structure",
  gate_indel_size_unique = "primer uniqueness"
)
KASP_GATE_LIMIT_INDEL_SIZE <- c(
  gate_indel_size_gap = "REF/ALT product difference >= min_gap bp",
  gate_indel_size_amplicon = glue::glue(
    "both products within [{KASP_INDEL_SIZE_LIMITS$amplicon[1L]}, {KASP_INDEL_SIZE_LIMITS$amplicon[2L]}] bp"
  ),
  gate_indel_size_tm = glue::glue(
    "both primers Tm within [{KASP_INDEL_SIZE_LIMITS$tm[1L]}, {KASP_INDEL_SIZE_LIMITS$tm[2L]}] C"
  ),
  gate_indel_size_structure = "hairpin/homodimer dG37 >= limit",
  gate_indel_size_unique = "no off-target exact primer hits; exactly one in-silico product"
)
# Human-readable limit, relaxation floor and hard-stop reason per gate (the
# `_LABEL_` spelling the other schemes use would exceed lintr's 30-char name cap
# with `INDEL_SIZE`).
KASP_GATE_FLOOR_INDEL_SIZE <- c(
  gate_indel_size_amplicon = "widen the amplicon window by 50 bp",
  gate_indel_size_tm = "widen the Tm window by 5 C",
  gate_indel_size_gap = "accept a 3 bp difference (capillary / PAGE)"
)
KASP_GATE_HARD_INDEL_SIZE <- c(
  gate_indel_size_structure = "a self-structured primer does not prime",
  gate_indel_size_unique = "a non-unique primer amplifies the wrong locus"
)

#' Design an InDel size-polymorphism marker
#'
#' The fragment-length readout: one flanking pair spans the insertion/deletion
#' and the two alleles give products that differ by the indel length, so the
#' genotype is read from band size (agarose or capillary).  Both primers sit in
#' the conserved flanks *outside* the indel, so each binds both alleles and the
#' product-size difference is exactly the indel length.  The pair is Tm-balanced
#' and kept in a short product window.  Every gate the scheme declares is
#' evaluated here (C7): the resolvability of the size difference (soft, with a
#' 3 bp capillary floor), the amplicon and Tm windows (soft), and the primer
#' structure and genome-wide uniqueness (hard).  The result records which soft
#' gates were relaxed (`$selection`).  Single-locus (1-plex) only.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`) with unequal-length
#'   `ref`/`alt` (a real insertion/deletion).
#' @param flank template half-width.
#' @param min_gap minimal REF/ALT product-size difference (bp) to call resolvable
#'   on the agarose readout; the ladder may relax it to `capillary_gap`.
#' @param min_product,max_product the accepted amplicon window (bp), applied to
#'   both alleles.
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows.
#' @param ... passed to [in_silico_pcr()] (e.g. `BPPARAM`).
#' @return a `marker_indel_size` list (`snp`, `primers`, `amplicon` with the REF
#'   and ALT product sizes, `gap`, `selection`), or `NULL` (with a message) when
#'   the indel is too short for even a capillary readout or no flanking pair
#'   passes the gates.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "ATTTTTT")
#' design_indel_size_marker(genome, snp)
#' }
#' @family KASP indel size
#' @export
design_indel_size_marker <- function(genome, snp, flank = 250L,
                                     min_gap = KASP_INDEL_SIZE_LIMITS$min_gap,
                                     min_product = KASP_INDEL_SIZE_LIMITS$amplicon[1L],
                                     max_product = KASP_INDEL_SIZE_LIMITS$amplicon[2L],
                                     lengths = 18:25,
                                     tm_range = KASP_INDEL_SIZE_LIMITS$tm,
                                     gc_range = KASP_LIMITS$gc, ...) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) == nchar(alt)) {
    cli::cli_inform(gettext("Not an InDel size polymorphism: `ref` and `alt` are the same length."))
    return(NULL)
  }
  delta <- nchar(alt) - nchar(ref)
  # the resolvable-gap gate is soft, but its floor is bounded: below a capillary
  # difference there is no marker, so fail here rather than search
  gap_floor <- -(as.integer(min_gap) - KASP_INDEL_SIZE_LIMITS$capillary_gap)
  if (abs(delta) - as.integer(min_gap) < gap_floor) {
    cli::cli_inform(gettext("The REF/ALT size difference is below `min_gap`; use capillary/PAGE or a longer indel."))
    return(NULL)
  }
  amp_pad <- KASP_INDEL_SIZE_FLOOR$amplicon
  tm_pad <- KASP_INDEL_SIZE_FLOOR$tm
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  # both primers must lie outside the whole REF allele block so each binds the
  # two alleles identically (the right flank merely shifts by `delta`)
  lo <- tmpl$snp_i
  hi <- tmpl$snp_i + nchar(ref) - 1L
  # the product window must hold BOTH alleles, widened by the gate's floor
  bar <- c(min_product - min(0L, delta), max_product - max(0L, delta))
  grid <- flanking_pair_grid(
    tmpl$seq, lo, hi, lengths,
    tm_range = c(tm_range[1L] - tm_pad, tm_range[2L] + tm_pad), gc_range = gc_range
  )
  if (is.null(grid)) {
    cli::cli_inform(gettext("No flanking primer pair passed the Tm/GC/amplicon filters for this InDel."))
    return(NULL)
  }
  grid <- grid |>
    dplyr::filter(dplyr::between(.data$product, bar[1L] - amp_pad, bar[2L] + amp_pad)) |>
    dplyr::mutate(
      m_gap = as.numeric(abs(delta) - as.integer(min_gap)),
      m_amplicon = pmin(.data$product - bar[1L], bar[2L] - .data$product),
      m_tm = pmin(
        .data$tm_f - tm_range[1L], tm_range[2L] - .data$tm_f,
        .data$tm_r - tm_range[1L], tm_range[2L] - .data$tm_r
      ),
      # the structure gate is confirmed (OligoArrayAux) and the full specificity
      # check (in-silico PCR) is confirmed on the shortlist; the seed margin is
      # the cheap pre-screen the ladder ranks by
      m_structure = 0,
      m_unique = seed_margin(.data$primer_f, .data$primer_r, genome),
      tier1 = -abs(.data$tm_f - .data$tm_r),
      tier2 = .data$m_amplicon,
      tier3 = -pmax(.data$run_f, .data$run_r)
    )
  confirm <- function(b) {
    if (nrow(in_silico_pcr(b$primer_f, b$primer_r, genome, ...)) != 1L) {
      return(FALSE)
    }
    pair_structure_ok(b$primer_f, b$primer_r)
  }
  sel <- select_flanking_pair(
    grid, KASP_GATES_INDEL_SIZE, KASP_RELAX_INDEL_SIZE, names(KASP_TIERS_INDEL_SIZE),
    margins = c(
      gate_indel_size_gap = "m_gap", gate_indel_size_amplicon = "m_amplicon",
      gate_indel_size_tm = "m_tm", gate_indel_size_structure = "m_structure",
      gate_indel_size_unique = "m_unique"
    ),
    floors = c(
      gate_indel_size_amplicon = -amp_pad, gate_indel_size_tm = -tm_pad,
      gate_indel_size_gap = gap_floor
    ),
    ties = c("dtm", "product"), confirm = confirm
  )
  if (is.null(sel$best)) {
    # name the gate that actually doomed the design, not a blanket reason
    if (any(!is.na(grid$m_unique)) && !any(grid$m_unique >= 0, na.rm = TRUE)) {
      cli::cli_inform(gettext("No flanking pair has unique 3' ends (every candidate's 3' seed repeats in the genome)."))
    } else {
      cli::cli_inform(gettext("No flanking pair passed the InDel-size gates (gap, windows, structure, uniqueness)."))
    }
    return(NULL)
  }
  primers <- pair_to_list(sel$best)
  off <- snp$pos - tmpl$snp_i
  out <- list(
    snp = snp, primers = primers,
    amplicon = tibble::tibble(
      start = primers$start + off, end = primers$end + off,
      ref_size = primers$product, alt_size = primers$product + delta
    ),
    gap = abs(delta),
    selection = list(relaxed = sel$relaxed, violated = sel$violated)
  )
  structure(out, class = "marker_indel_size")
}

#' @export
print.marker_indel_size <- function(x, ...) {
  cat(
    "<marker_indel_size>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat(
    "  amplicon: REF", x$amplicon$ref_size, "bp / ALT", x$amplicon$alt_size,
    "bp (gap", x$gap, "bp)\n"
  )
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  if (length(x$selection$relaxed)) {
    cat("  relaxed:", paste(x$selection$relaxed, collapse = ", "), "\n")
  }
  invisible(x)
}

#' Declare the InDel size-polymorphism scheme
#'
#' The metadata the shared registry (see [marker_schemes()], [marker_criteria()] and
#' `marker_decision_stages()`) needs to wire this chemistry into the decision
#' schema: its family/label, how it reads the genotype, and its selection gates,
#' refinement order and tiers.  Kept here, next to the designer, so the two
#' cannot drift.  Gate ids are prefixed `gate_indel_size_`.  Every gate declared
#' here is evaluated by [design_indel_size_marker()] (C7).
#' @return a list describing the `indel_size` scheme.
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_indel_size
#' @export
#' @noRd
marker_map_spec.marker_indel_size <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  start <- sp[["start"]]
  end <- sp[["end"]]
  map_spec(
    marker,
    features = c(
      map_pair_rows(marker, start, end),
      list(
        spec_row(
          "REF product", start, start + marker$amplicon$ref_size[1L] - 1L, "+", "amplicon",
          label = "REF"
        ),
        spec_row(
          "ALT product", start, start + marker$amplicon$alt_size[1L] - 1L, "+", "amplicon",
          label = "ALT"
        )
      )
    ),
    sequences = list(forward = marker$primers$forward, reverse = marker$primers$reverse)
  )
}

scheme_def_indel_size <- function() {
  list(
    scheme = "indel_size",
    family = "indel_size",
    label = "InDel size polymorphism (gel/capillary fragment length)",
    discrimination = "amplicon size difference equals the indel length",
    readout = "agarose / capillary fragment length",
    zygosity = "co-dominant",
    amplicon = "bounded (80-300 bp)",
    loci = 1L,
    verdict = "usable InDel size marker",
    gates = KASP_GATES_INDEL_SIZE,
    relax = KASP_RELAX_INDEL_SIZE,
    tiers = KASP_TIERS_INDEL_SIZE,
    gate_short = KASP_GATE_SHORT_INDEL_SIZE,
    gate_limit = KASP_GATE_LIMIT_INDEL_SIZE,
    gate_floor_label = KASP_GATE_FLOOR_INDEL_SIZE,
    gate_hard_label = KASP_GATE_HARD_INDEL_SIZE,
    stage_label = "InDel size marker"
  )
}
