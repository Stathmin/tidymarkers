# =============================================================================
#  kasp_ssr.R  --  SSR / microsatellite genotyping (a size marker)
# =============================================================================
#  An SSR (simple sequence repeat, microsatellite / STR) marker is not an
#  allele-specific assay: *one* flanking pair amplifies the repeat locus and the
#  genotype is read from the amplicon *size* - the alleles differ by the number
#  of repeat copies.  So the design is the generic "any primer pair" core plus
#  one extra rule: the primers must sit in the UNIQUE flanks, never over the
#  repeat.  The repeat span is located with `repeat_spans()`, the pair is placed
#  with `design_flanking_pair()` (which already rejects any primer carrying a
#  homopolymer >= 4 or an SSR via `primer_candidates()`), and the amplicon is
#  kept in the ~100-300 bp size-marker window.
#
#  Design rules, grounded on the standard guidance:
#   * a legal primer pair must flank the repeat, in UNIQUE sequence - the repeat
#     itself is excluded from the primers (Primer3 manual: SEQUENCE_TARGET for a
#     CA repeat makes a legal pair flank the target and SEQUENCE_INTERNAL_EXCLUDED
#     keeps the repeat out of the oligo; Wikipedia "Microsatellite": "using the
#     unique sequences of flanking regions as primers").
#   * product ~100-300 bp (Primer3 default PRIMER_PRODUCT_SIZE_RANGE=100-300).
#   * composition / length / Tm / GC and the no-repeat, no-run rules are the
#     generic PCR ones: Primer3 defaults are 18-27 nt, Tm 57-63 C, GC 20-80%,
#     PRIMER_MAX_POLY_X=5; Premier Biosoft recommends 18-22 nt, Tm 52-58 C,
#     GC 40-60%, <= 4 dinucleotide repeats and <= 4 bp runs per primer.
#   * long repeats stutter and short ones amplify poorly, so tetra-/penta-
#     nucleotide repeats are preferred for clean sizing (Wikipedia).  We do not
#     model stutter: it is a readout property, not a design gate here.
#   * multiplexing needs non-overlapping size ranges / different dyes; we design
#     a single locus and leave the pooling to the caller.
#  We reuse the package's composition window (KASP_LIMITS$gc, max_run,
#  has_repeat), so a primer is clean by construction; only the amplicon window
#  and the Tm window are SSR-specific inputs below.
#
#  Sources (fetched 2026-10-06):
#   * https://primer3.org/manual.html
#     (SEQUENCE_TARGET, PRIMER_PRODUCT_SIZE_RANGE, PRIMER_MAX_POLY_X,
#      PRIMER_MIN_TM/PRIMER_MAX_TM, PRIMER_MIN_GC/PRIMER_MAX_GC)
#   * https://en.wikipedia.org/wiki/Microsatellite
#     (flanking unique primers; stutter / repeat-length trade-off)
#   * https://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
#     (18-22 nt, Tm 52-58 C, GC 40-60%, <= 4 dinucleotide repeats, <= 4 bp runs)
# =============================================================================

#' The SSR span nearest (or containing) an anchor in a template
#'
#' Reuses [repeat_spans()] for the hit list and returns the widest/nearest span:
#' a span containing `anchor` wins (the smallest one), otherwise the span whose
#' midpoint is closest.
#' @param seq the template (plus strand).
#' @param anchor template coordinate of the locus anchor.
#' @param units repeat-unit lengths to search for.
#' @param min_copies tandem copies required for a repeat call.
#' @return a one-row tibble (`start`, `end`) in template coordinates, or `NULL`.
#' @keywords internal
#' @noRd
ssr_span_at <- function(seq, anchor, units = 2:4, min_copies = 3L) {
  spans <- repeat_spans(as.character(seq), units, min_copies)
  if (!nrow(spans)) {
    return(NULL)
  }
  sp <- tibble::tibble(
    start = as.integer(spans[, 1L]),
    end = as.integer(spans[, 2L])
  )
  inside <- dplyr::filter(sp, .data$start <= anchor, .data$end >= anchor)
  if (nrow(inside)) {
    return(dplyr::slice_min(inside, .data$end - .data$start, n = 1L))
  }
  sp |>
    dplyr::mutate(dist = abs((.data$start + .data$end) / 2 - anchor)) |>
    dplyr::slice_min(.data$dist, n = 1L) |>
    dplyr::select(-"dist")
}

#' Design an SSR / microsatellite marker at a repeat locus
#'
#' The size-marker analogue of [design_hrm_marker()]: a single flanking pair
#' straddles the repeat, which sits under neither primer, and the amplicon is
#' sized to the ~100-300 bp window a capillary/gel readout wants.  The repeat is
#' located around `snp$pos`, the pair is placed in the unique flanks by
#' [design_flanking_pair()] (so both primers lie outside the repeat and carry no
#' homopolymer >= 4 or SSR), and the alleles are read as the amplicon size.
#'
#' `snp$pos` is only the *anchor* of the window searched for the repeat - an SSR
#' locus is defined by its repeat, not by the REF/ALT bases.  Strong stutter is
#' expected for long repeats and multiplexing needs non-overlapping size ranges,
#' but neither is modelled here: stutter is a readout property, and pooling is
#' left to the caller.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`); `pos` anchors the
#'   window searched for the repeat.
#' @param flank template half-width.
#' @param min_product,max_product the accepted amplicon window (SSR size marker).
#' @param lengths primer lengths to try.
#' @param units repeat-unit lengths to search for.
#' @param min_copies tandem copies required for a repeat call.
#' @param tm_range,gc_range accepted primer windows.
#' @param ... reserved for future options.
#' @return a `kasp_ssr` list (`snp`, `primers`, `amplicon`, `repeat`), or `NULL`
#'   (with a message) when no repeat or flanking pair passes the filters.  The
#'   `primers`, `amplicon` and `repeat` coordinates are genomic.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_ssr_marker(genome, kasp_example("SNP")[1, ])
#' }
#' @family KASP SSR
#' @export
design_ssr_marker <- function(genome, snp, flank = 250L,
                              min_product = 100L, max_product = 300L,
                              lengths = 18:25, units = 2:4, min_copies = 3L,
                              tm_range = KASP_LIMITS$com_tm, gc_range = KASP_LIMITS$gc, ...) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  repeat_span <- ssr_span_at(tmpl$seq, tmpl$snp_i, units, min_copies)
  if (is.null(repeat_span)) {
    cli::cli_inform(gettext("No SSR / microsatellite repeat found in the flanking window."))
    return(NULL)
  }
  # the pair must flank the repeat; design_flanking_pair() also drops any primer
  # with a homopolymer >= 4 or an SSR (primer_candidates()), so the primers are
  # clean by construction
  primers <- design_flanking_pair(
    tmpl$seq,
    lo = repeat_span$start, hi = repeat_span$end, lengths = lengths,
    min_product = min_product, max_product = max_product,
    tm_range = tm_range, gc_range = gc_range
  )
  if (is.null(primers)) {
    cli::cli_inform(gettext("No flanking primer pair outside the repeat passed the Tm/GC filters."))
    return(NULL)
  }
  off <- snp$pos - tmpl$snp_i
  primers$start <- primers$start + off
  primers$end <- primers$end + off
  out <- list(
    snp = snp,
    primers = primers,
    amplicon = tibble::tibble(start = primers$start, end = primers$end)
  )
  # `repeat` is an R reserved word, so the field can only be set/read with [[
  out[["repeat"]] <- tibble::tibble(
    start = repeat_span$start + off,
    end = repeat_span$end + off
  )
  structure(out, class = "kasp_ssr")
}

#' @export
print.kasp_ssr <- function(x, ...) {
  cat(
    "<kasp_ssr>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  rep_span <- x[["repeat"]]
  cat("  repeat:", rep_span$start, "-", rep_span$end,
    "(", rep_span$end - rep_span$start + 1L, "bp )\n",
    sep = ""
  )
  cat("  amplicon:", x$primers$product, "bp\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  invisible(x)
}

#' The SSR / microsatellite scheme definition
#'
#' The registry row and decision-schema fragments for the size marker: the
#' flanking pair must sit outside the repeat, the primers must be clean, the
#' amplicon stays in the size window and the pair is Tm-balanced.  Only the
#' in-window gates are evaluated by [design_ssr_marker()]; off-target uniqueness
#' is a genome-wide check the registry runs on the chosen pair.
#' @return a list (`scheme`, `family`, `label`, `discrimination`, `readout`,
#'   `zygosity`, `amplicon`, `loci`, `verdict`, `gates`, `relax`, `tiers`,
#'   `gate_short`, `gate_limit`, `gate_floor_label`, `gate_hard_label`,
#'   `stage_label`).
#' @keywords internal
#' @noRd
scheme_def_ssr <- function() {
  gates <- c(
    "gate_ssr_amplicon", "gate_ssr_flanks", "gate_ssr_clean",
    "gate_ssr_tm", "gate_ssr_offtarget"
  )
  list(
    scheme = "ssr",
    family = "ssr",
    label = "SSR / microsatellite (simple-sequence-repeat size marker)",
    discrimination = "amplicon size tracks the repeat copy number",
    readout = "capillary / gel fragment size",
    zygosity = "co-dominant",
    amplicon = "~100-300 bp",
    loci = 1L,
    verdict = "usable SSR marker",
    gates = gates,
    relax = c("gate_ssr_amplicon", "gate_ssr_tm", "gate_ssr_offtarget"),
    tiers = c(
      tier1 = "repeat located + clean flanks",
      tier2 = "Tm balance",
      tier3 = "amplicon length"
    ),
    gate_short = c(
      gate_ssr_amplicon = "amplicon length (SSR)",
      gate_ssr_flanks = "flanking primers",
      gate_ssr_clean = "primer repeat/run",
      gate_ssr_tm = "primer Tm (SSR)",
      gate_ssr_offtarget = "off-target hits"
    ),
    gate_limit = c(
      gate_ssr_amplicon = "amplicon within [100, 300]",
      gate_ssr_flanks = "both primers outside the repeat",
      gate_ssr_clean = "no SSR / homopolymer >= 4 in the primers",
      gate_ssr_tm = "both primers Tm within [58, 63] C",
      gate_ssr_offtarget = "off-target exact primer hits <= 1"
    ),
    gate_floor_label = c(
      gate_ssr_amplicon = "widen the amplicon window",
      gate_ssr_tm = "widen the Tm window",
      gate_ssr_offtarget = "allow one extra off-target hit"
    ),
    gate_hard_label = c(
      gate_ssr_flanks = "a primer over the repeat destroys the locus",
      gate_ssr_clean = "a repeat inside a primer misprimes"
    ),
    stage_label = "SSR / microsatellite"
  )
}
