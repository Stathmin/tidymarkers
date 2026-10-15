# =============================================================================
#  kasp_castpcr.R  --  castPCR / allele-specific qPCR (competitive AS TaqMan)
# =============================================================================
#  castPCR is KASP's real-time cousin: ONE allele-specific primer whose 3' base
#  sits on the variant, a locus-specific hydrolysis (TaqMan) probe and a common
#  reverse primer.  The 3' mismatch is still the discriminator, but the signal is
#  read as an allele-specific Ct from the probe rather than a FRET ratio; a
#  competitive MGB oligonucleotide blocker can additionally suppress the
#  non-target allele.  The amplicon is kept SHORT, so it fits a fast real-time
#  protocol exactly as a KASP product does.
#
#  Grounding (design rules verified against the primary method papers):
#   - Roma et al., BioMed Res Int 2013, 2013:385087 (PMC3864078): "An allele-
#     specific primer and a FAM dye-labelled MGB probe detect the mutant allele,
#     while an MGB oligonucleotide blocker suppresses the wild type allele ...
#     enabling the detection of as little as 0.1% mutant allele."
#   - Barbano et al., Sci Rep 2015, 5:18592 (PMC4686972): castPCR combines a
#     mutant allele-specific primer (ASP) with a wild-type allele-specific
#     blocker (ASB); "a molecular blocker completely suppresses the amplification
#     of the wild type allele"; LOD ~1:100 (0.5-1% mutant/not-mutated).
#   - Husain et al., Asian Pac J Cancer Prev 2023, 24(3):961 (PMC10334107):
#     CAST-PCR pairs real-time PCR with an MGB-blocking oligonucleotide and reads
#     Delta-Ct (mutant assay - gene reference assay).
#  Design thresholds that the method papers do not state numerically are taken
#  from the package's own KASP limits where they coincide (AS Tm, GC) and from
#  the standard hydrolysis-probe rule of thumb - probe Tm above the primers, short
#  amplicon - which is stated as an ASSUMPTION below and is a parameter the caller
#  can override.  That vendor rule (~8-10 C) is not reachable on the package's
#  compressed NN Tm scale, so `probe_margin` is a smaller practical floor.
# =============================================================================

# Hydrolysis-probe rule of thumb: a TaqMan/MGB probe must melt above the primers,
# classically ~8-10 C higher on the vendor Tm scale, so the probe carries the
# real-time readout.  The package's NN Tm at PCR conditions is compressed (see the
# TaqMan module), so a literal 8-10 C margin is not reachable here; the practical
# floor is `CASTPCR_PROBE_MARGIN` (a few degrees), capped at `CASTPCR_PROBE_TM_MAX`.
# An MGB modification raises the effective Tm further, so this lower margin is
# conservative for the unmodified oligo.
CASTPCR_PROBE_MARGIN <- 4
CASTPCR_PROBE_TM_MAX <- 80

# castPCR keeps the product short so the real-time protocol stays fast; the
# papers' assays are TaqMan products in the ~60-150 bp range.
CASTPCR_AMPLICON <- c(60L, 150L)

# The castPCR decision stage, mirroring the CAPS/ARMS/HRM stage constants in
# R/kasp_select.R.  `scheme_def_castpcr()` returns them so the shared registry
# can wire the chemistry in afterwards without a second vocabulary.
KASP_GATES_CASTPCR <- c(
  "gate_castpcr_as3", "gate_castpcr_mismatch", "gate_castpcr_as_tm",
  "gate_castpcr_probe_tm", "gate_castpcr_amplicon"
)
KASP_RELAX_CASTPCR <- c(
  "gate_castpcr_amplicon", "gate_castpcr_as_tm", "gate_castpcr_probe_tm"
)
KASP_TIERS_CASTPCR <- c(
  tier1 = "3'-terminal specificity", tier2 = "Tm ordering", tier3 = "amplicon length"
)

#' Hydrolysis-probe candidates covering the variant
#'
#' Enumerates plus-strand probes that span the variant base and pass the GC/Tm
#' and composition filters, ordered so the most centred probe comes first.  Reuses
#' [primer_candidates()] for every composition check.
#' @param seq the template (plus strand).
#' @param snp_i 1-based position of the variant on `seq`.
#' @param lengths probe lengths to try.
#' @param tm_range,gc_range accepted windows.
#' @return a tibble (`start`, `len`, `end`, `primer`, `tm`, `gc`), possibly zero
#'   rows, ordered by distance of the probe centre to the variant.
#' @keywords internal
#' @noRd
castpcr_probe_candidates <- function(seq, snp_i, lengths, tm_range, gc_range) {
  cand <- spanning_probe_candidates(seq, snp_i, lengths, tm_range, gc_range)
  if (!nrow(cand)) {
    return(cand)
  }
  cand |>
    dplyr::mutate(centre_dist = abs((.data$start + .data$end) / 2 - snp_i)) |>
    dplyr::arrange(.data$centre_dist, abs(.data$tm - mean(tm_range)))
}

#' Design a castPCR / allele-specific qPCR marker at a variant
#'
#' The real-time analogue of [assemble_kasp_assay()]: one allele-specific primer
#' whose 3' base sits on the variant ([design_as_primers()]), a locus-specific
#' hydrolysis probe that covers the variant, and a common reverse primer, all
#' packed into a short amplicon.  The 3' mismatch is the discriminator and the
#' probe is the readout, so the same variant this places in an AS 3' window is a
#' KASP placement seen through a real-time probe; the hydrolysis probe is required
#' to melt above the primers (`probe_margin`, see the file header for the
#' assumption).  Single-base substitutions only.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param as_lengths allele-specific primer lengths to try.
#' @param mate_lengths common-primer lengths to try.
#' @param probe_lengths hydrolysis-probe lengths to try.
#' @param probe_margin minimum probe Tm excess (C) over the AS primer (probe must
#'   melt above the primer); defaults to `CASTPCR_PROBE_MARGIN`.
#' @param probe_tm_max upper bound (C) on the probe Tm; defaults to
#'   `CASTPCR_PROBE_TM_MAX`.
#' @param min_product,max_product the accepted amplicon window.
#' @param tm_range accepted AS primer Tm window.
#' @param gc_range accepted GC window for every oligo.
#' @return a `kasp_castpcr` list (`snp`, `primers`, `probe`, `amplicon`,
#'   `compatibility` - the joint dimer/hairpin screen over the AS primer, its mate
#'   and the probe), or `NULL` (with a message) when no AS primer, mate or probe
#'   passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_castpcr_marker(genome, kasp_example("SNP")[1, ])
#' }
#' @family KASP castPCR
#' @export
design_castpcr_marker <- function(genome, snp, flank = 150L,
                                  as_lengths = 18:27, mate_lengths = 18:25,
                                  probe_lengths = 18:32,
                                  probe_margin = CASTPCR_PROBE_MARGIN,
                                  probe_tm_max = CASTPCR_PROBE_TM_MAX,
                                  min_product = CASTPCR_AMPLICON[1L],
                                  max_product = CASTPCR_AMPLICON[2L],
                                  tm_range = KASP_LIMITS$as_tm, gc_range = KASP_LIMITS$gc) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    cli::cli_inform(gettext("castPCR is implemented for single-base substitutions only."))
    return(NULL)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  si <- tmpl$snp_i
  n <- nchar(as.character(tmpl$seq))
  # the AS primer carries the variant on its 3' base; one placement is chosen,
  # not a pair (the probe supplies the real-time signal, the mate is common)
  as_cand <- design_as_primers(tmpl$seq, si, ref, alt,
    strands = c("+", "-"), lengths = as_lengths,
    gc_range = gc_range, tm_range = tm_range
  )
  if (!nrow(as_cand)) {
    cli::cli_inform(gettext("No allele-specific primer passed the GC/Tm filters for castPCR."))
    return(NULL)
  }
  best <- as_cand[which.max(as_cand$discrimination), ]
  # a one-row candidate set in the schema [primer_candidates()] produces, so the
  # generic [pick_primer_pair()] can pair it with the common mate
  as_one <- tibble::tibble(
    start = best$tmpl_start, end = best$tmpl_end,
    primer = best$ref_primer, tm = best$tm
  )
  if (identical(best$strand, "+")) {
    mate <- primer_candidates(
      tmpl$seq, seq.int(best$tmpl_end + 1L, n), mate_lengths, "-",
      KASP_LIMITS$com_tm, gc_range
    )
    b <- pick_primer_pair(as_one, mate, min_product, max_product)
  } else {
    if (best$tmpl_start <= 1L) {
      cli::cli_inform(gettext("No room upstream of the allele-specific primer for a common primer."))
      return(NULL)
    }
    mate <- primer_candidates(
      tmpl$seq, seq_len(best$tmpl_start - 1L), mate_lengths, "+",
      KASP_LIMITS$com_tm, gc_range
    ) |>
      dplyr::filter(.data$end < best$tmpl_start)
    b <- pick_primer_pair(mate, as_one, min_product, max_product)
  }
  if (is.null(b)) {
    cli::cli_inform(gettext("No common primer passed the Tm/GC filters within the castPCR amplicon."))
    return(NULL)
  }
  # the probe must sit inside the amplicon, cover the variant and melt above the
  # AS primer (the classical TaqMan ordering rule)
  probe_tm <- c(best$tm + probe_margin, probe_tm_max)
  probe <- castpcr_probe_candidates(tmpl$seq, si, probe_lengths, probe_tm, gc_range) |>
    dplyr::filter(.data$start >= b$start_f, .data$end <= b$end_r)
  if (!nrow(probe)) {
    cli::cli_inform(gettext("No hydrolysis probe covering the variant passed the Tm/GC filters."))
    return(NULL)
  }
  probe <- probe[1L, ]
  off <- snp$pos - si
  pair <- pair_to_list(b)
  out <- list(
    snp = snp,
    primers = c(pair, list(
      as_role = if (identical(best$strand, "+")) "forward" else "reverse",
      as_primer = best$ref_primer,
      as_tm = best$tm,
      as_discrimination = best$discrimination
    )),
    probe = list(
      sequence = probe$primer, start = probe$start + off, end = probe$end + off,
      tm = probe$tm, gc = probe$gc
    ),
    amplicon = tibble::tibble(start = pair$start + off, end = pair$end + off)
  )
  # AS primer + common mate + hydrolysis probe share the tube: joint screen
  out$compatibility <- tube_compatibility(c(
    as_primer = out$primers$as_primer, common = out$primers$forward,
    mate = out$primers$reverse, probe = out$probe$sequence
  ))
  structure(out, class = "kasp_castpcr")
}

#' @export
print.kasp_castpcr <- function(x, ...) {
  cat(
    "<kasp_castpcr>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  AS primer (", x$primers$as_role, "): ", x$primers$as_primer,
    " (Tm ", round(x$primers$as_tm, 1L), " C)\n",
    sep = ""
  )
  cat("  probe:", x$probe$sequence, " (Tm ", round(x$probe$tm, 1L), " C)\n", sep = "")
  cat("  amplicon:", x$primers$product, "bp\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  invisible(x)
}

#' The castPCR scheme definition for the shared registry
#'
#' A self-describing definition of the castPCR chemistry in the shape the marker
#' registry expects (scheme/family/label/readout, the gate/tier vocabulary and the
#' per-gate human-readable labels).  Kept beside the designer so the chemistry and
#' its thresholds travel together; it is discovered automatically (one file per
#' scheme, see `R/kasp_scheme_registry.R`).
#' @return a named list describing the scheme.
#' @keywords internal
#' @noRd
scheme_def_castpcr <- function() {
  as_tm <- KASP_LIMITS$as_tm
  margin <- CASTPCR_PROBE_MARGIN
  amplicon <- CASTPCR_AMPLICON
  list(
    scheme = "castpcr",
    family = "castpcr",
    label = "castPCR / allele-specific qPCR (competitive AS TaqMan)",
    discrimination = "3'-terminal mismatch of the AS primer (qPCR-read)",
    readout = "hydrolysis (TaqMan) probe, allele-specific Ct",
    zygosity = "co-dominant",
    amplicon = "short (60-150 bp)",
    loci = 1L,
    verdict = "usable castPCR assay",
    gates = KASP_GATES_CASTPCR,
    relax = KASP_RELAX_CASTPCR,
    tiers = KASP_TIERS_CASTPCR,
    gate_short = c(
      gate_castpcr_as3 = "AS 3' base on variant",
      gate_castpcr_mismatch = "3' mismatch (AS)",
      gate_castpcr_as_tm = "AS primer Tm",
      gate_castpcr_probe_tm = "probe Tm",
      gate_castpcr_amplicon = "amplicon length (castPCR)"
    ),
    gate_limit = c(
      gate_castpcr_as3 = "AS 3' base = the variant base",
      gate_castpcr_mismatch = "AS 3' base mismatches the non-target allele",
      gate_castpcr_as_tm = glue::glue("AS Tm within [{as_tm[1L]}, {as_tm[2L]}] C"),
      gate_castpcr_probe_tm = glue::glue("probe Tm >= AS Tm + {margin} C"),
      gate_castpcr_amplicon = glue::glue("amplicon within [{amplicon[1L]}, {amplicon[2L]}]")
    ),
    gate_floor_label = c(
      gate_castpcr_as_tm = "widen the AS Tm window",
      gate_castpcr_probe_tm = "widen the probe Tm window",
      gate_castpcr_amplicon = "widen the amplicon window"
    ),
    gate_hard_label = c(
      gate_castpcr_as3 = "the AS 3' base must sit on the variant",
      gate_castpcr_mismatch = "no 3' mismatch = no allele discrimination"
    ),
    stage_label = "castPCR / AS-qPCR"
  )
}
