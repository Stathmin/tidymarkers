# =============================================================================
#  kasp_indel_size.R  --  InDel size-polymorphism genotyping (gel / capillary)
# =============================================================================
#  An InDel size marker reads the genotype from the *length* of one PCR product:
#  a single flanking pair spans the insertion/deletion, so the REF and ALT
#  amplicons differ by exactly the indel length and AA / AB / BB resolve by
#  fragment size - on a plain agarose gel, or (for small differences) on a
#  capillary sequencer.
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
#     ~5-10 bp, so small InDels (1-5 bp) need PAGE or capillary electrophoresis
#     (Wang et al. 2026, Int J Mol Sci 27:1406,
#     <https://doi.org/10.3390/ijms27031406>).  We default `min_gap = 5` bp and
#     treat anything shorter as a capillary job.
#   * Long InDels (>30 bp) are the agarose-friendly ones, and differences
#     >= 100 bp score unambiguously (Chowdhury et al. 2026, Plant Biotechnol
#     43:227, <https://doi.org/10.5511/plantbiotechnology.26.0420a>).
#  The primer scan also drops homopolymer runs / SSRs, which is the practical
#  "unique flank" screen.  The gates that need the model genome (primer
#  structure, off-target uniqueness) are declared by scheme_def_indel_size() and
#  evaluated by the selector, exactly as the HRM/CAPS stages do.
# =============================================================================

#' Design an InDel size-polymorphism marker
#'
#' The fragment-length readout: one flanking pair spans the insertion/deletion
#' and the two alleles give products that differ by the indel length, so the
#' genotype is read from band size (agarose or capillary).  Both primers sit in
#' the conserved flanks *outside* the indel, so each binds both alleles and the
#' product-size difference is exactly the indel length.  The pair is Tm-balanced
#' and kept in a short product window; a difference below `min_gap` is treated as
#' a capillary/PAGE job and rejected.  Single-locus (1-plex) only.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`) with unequal-length
#'   `ref`/`alt` (a real insertion/deletion).
#' @param flank template half-width.
#' @param min_gap minimal REF/ALT product-size difference (bp) to call resolvable.
#' @param min_product,max_product the accepted amplicon window (bp), applied to
#'   both alleles.
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows.
#' @param ... reserved for future arguments.
#' @return a `kasp_indel_size` list (`snp`, `primers`, `amplicon` with the REF
#'   and ALT product sizes, `gap`), or `NULL` (with a message) when the indel is
#'   too short or no flanking pair passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "ATTTTTT")
#' design_indel_size_marker(genome, snp)
#' }
#' @family KASP indel size
#' @export
design_indel_size_marker <- function(genome, snp, flank = 250L, min_gap = 5L,
                                     min_product = 80L, max_product = 300L,
                                     lengths = 18:25,
                                     tm_range = KASP_LIMITS$com_tm,
                                     gc_range = KASP_LIMITS$gc, ...) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) == nchar(alt)) {
    cli::cli_inform(gettext("Not an InDel size polymorphism: `ref` and `alt` are the same length."))
    return(NULL)
  }
  delta <- nchar(alt) - nchar(ref)
  if (abs(delta) < min_gap) {
    cli::cli_inform(gettext("The REF/ALT size difference is below `min_gap`; use capillary/PAGE or a longer indel."))
    return(NULL)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  # both primers must lie outside the whole REF allele block so each binds the
  # two alleles identically (the right flank merely shifts by `delta`)
  lo <- tmpl$snp_i
  hi <- tmpl$snp_i + nchar(ref) - 1L
  # the product window must hold BOTH alleles: alt = ref + delta, so the REF
  # product is bounded by the window shifted by delta (negative if impossible)
  win <- c(min_product - min(0L, delta), max_product - max(0L, delta))
  primers <- design_flanking_pair(
    tmpl$seq, lo, hi, lengths, win[1L], win[2L], tm_range, gc_range
  )
  if (is.null(primers)) {
    cli::cli_inform(gettext("No flanking primer pair passed the Tm/GC/amplicon filters for this InDel."))
    return(NULL)
  }
  off <- snp$pos - tmpl$snp_i
  out <- list(
    snp = snp, primers = primers,
    amplicon = tibble::tibble(
      start = primers$start + off, end = primers$end + off,
      ref_size = primers$product, alt_size = primers$product + delta
    ),
    gap = abs(delta)
  )
  structure(out, class = "kasp_indel_size")
}

#' @export
print.kasp_indel_size <- function(x, ...) {
  cat(
    "<kasp_indel_size>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat(
    "  amplicon: REF", x$amplicon$ref_size, "bp / ALT", x$amplicon$alt_size,
    "bp (gap", x$gap, "bp)\n"
  )
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  invisible(x)
}

#' Declare the InDel size-polymorphism scheme
#'
#' The metadata the shared registry (see [kasp_schemes()], [kasp_criteria()] and
#' [kasp_decision_stages()]) needs to wire this chemistry into the decision
#' schema: its family/label, how it reads the genotype, and its selection gates,
#' refinement order and tiers.  Kept here, next to the designer, so the two
#' cannot drift.  Gate ids are prefixed `gate_indel_size_`.
#' @return a list describing the `indel_size` scheme.
#' @keywords internal
#' @noRd
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
    gates = c(
      "gate_indel_size_gap", "gate_indel_size_amplicon", "gate_indel_size_tm",
      "gate_indel_size_structure", "gate_indel_size_unique"
    ),
    relax = c("gate_indel_size_amplicon", "gate_indel_size_tm", "gate_indel_size_gap"),
    tiers = c(
      tier1 = "resolvable size gap", tier2 = "Tm balance", tier3 = "amplicon length"
    ),
    gate_short = c(
      gate_indel_size_gap = "resolvable gap",
      gate_indel_size_amplicon = "amplicon length (InDel)",
      gate_indel_size_tm = "primer Tm (InDel)",
      gate_indel_size_structure = "primer structure",
      gate_indel_size_unique = "primer uniqueness"
    ),
    gate_limit = c(
      gate_indel_size_gap = "REF/ALT product difference >= min_gap bp",
      gate_indel_size_amplicon = "both products within [80, 300] bp",
      gate_indel_size_tm = "both primers Tm within [lo, hi]",
      gate_indel_size_structure = "hairpin/homodimer dG37 >= limit",
      gate_indel_size_unique = "no off-target exact primer hits"
    ),
    gate_floor_label = c(
      gate_indel_size_amplicon = "widen the amplicon window",
      gate_indel_size_tm = "widen the Tm window",
      gate_indel_size_gap = "allow a 3 bp gap (capillary)"
    ),
    gate_hard_label = c(
      gate_indel_size_structure = "a self-structured primer does not prime",
      gate_indel_size_unique = "a non-unique primer amplifies the wrong locus"
    ),
    stage_label = "InDel size marker"
  )
}
