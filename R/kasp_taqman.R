# =============================================================================
#  kasp_taqman.R  --  TaqMan 5'-nuclease allelic discrimination (probes)
# =============================================================================
#  TaqMan reads the genotype from two allele-specific *hydrolysis probes* in one
#  tube: one flanking primer pair amplifies a short amplicon, and the two probes
#  (one per allele, each carrying its own fluorophore - conventionally FAM/HEX)
#  are cleaved by the 5'->3' exonuclease activity of Taq only when they match
#  their template, so the fluorescence ratio is the call.  Unlike KASP there is
#  no allele-specific primer; the discrimination lives in the probe, which is why
#  a mismatch near the *middle* of the probe (maximum duplex destabilisation) is
#  the design target.
#
#  Design rules encoded here, and where they come from:
#   * one flanking pair + two allele-specific probes, each probe centred over the
#     variant; the probe anneals *within* the region the primers amplify
#     (Wikipedia, "TaqMan", fetched 2026-10-06).
#   * probe 15-30 nt, primer 15-30 nt, both GC 30-80%; canonical qPCR probes have
#     a Tm 8-10 C above the primer Tm (68-70 C probes for 50-60 C primers).
#     MacLeod, Rowley & Essex 2021, Commun. Biol. 4:227, Table 1 ("Oligonucleotide
#     Tm: 50-60 C primers / 68-70 C probes") and its "Probe design" section
#     ("Canonical qPCR probes should have a Tm 8-10 C higher than the primer Tm");
#     full text via Europe PMC fullTextXML, PMC7892852.
#   * the discriminating base sits at the probe centre (odd length) or is biased
#     toward the probe 3' terminus for even lengths.  Same source, "Probe design".
#   * a 5'-terminal G quenches the reporter even after cleavage, so probes avoid
#     it; this is the Applied Biosystems / Thermo Fisher TaqMan design guideline
#     (thermofisher.com taqman-probe-design.html; the fetch was blocked by the
#     vendor, so this is stated as the widely used vendor rule / assumption).
#
#  Assumption on the Tm scale: TmCalculator's nearest-neighbour model at the
#  package PCR conditions returns values below the vendor numbers, so the default
#  probe window (`probe_tm = c(56, 72)`) is set on that model rather than the
#  literal 68-70 C band.  The vendor ordering rule ("probe Tm 8-10 C above the
#  primers") is NOT enforced as a hard gate: on this compressed scale even a
#  30-mer probe reaches only ~+4 C over the primers, so a fixed 8 C margin would
#  doom every design.  The ordering is kept as a *ranking* preference (tier2: the
#  probe should be the hottest oligo in the tube), and a true margin gate is
#  queued in the roadmap pending a vendor-scale Tm.  The amplicon is kept short
#  (~60-150 bp) as in the usual TaqMan SNP-genotyping assays.
# =============================================================================

#' Probe candidates centred over a variant
#'
#' Reuses [primer_candidates()] for the composition screen (GC, Tm, runs, SSR),
#' keeps the plus-strand windows that span the variant and avoid a 5'-G, and adds
#' the variant's position in the probe and a `centrality` in `[0, 1]` (1 = the
#' variant is exactly central, 0 = at an end).
#' @param seq the template (plus strand).
#' @param snp_i template coordinate of the variant.
#' @param lengths probe lengths to try.
#' @param tm_range,gc_range accepted windows.
#' @return tibble (`start`, `len`, `end`, `primer`, `tm`, `gc`, `variant_pos`,
#'   `centrality`).
#' @keywords internal
#' @noRd
taqman_probe_candidates <- function(seq, snp_i, lengths, tm_range, gc_range) {
  spanning_probe_candidates(seq, snp_i, lengths, tm_range, gc_range) |>
    dplyr::filter(stringr::str_sub(.data$primer, 1L, 1L) != "G") |>
    dplyr::mutate(
      variant_pos = snp_i - .data$start + 1L,
      centrality = 1 - abs(.data$variant_pos - (.data$len + 1L) / 2) / ((.data$len - 1L) / 2)
    )
}

#' Design a TaqMan 5'-nuclease allelic-discrimination marker at a variant
#'
#' A single flanking primer pair amplifies a short amplicon and **two
#' allele-specific hydrolysis probes** (one per allele, different fluorophores,
#' FAM/HEX here) are placed over the variant.  Each probe is centred over the
#' variant so the mismatch is maximally destabilising, avoids a 5'-G, sits inside
#' the amplicon and passes the GC/Tm composition filters; the pair itself reuses
#' the marker-agnostic flanking-pair machinery.  Single-base substitutions only.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param min_product,max_product the accepted amplicon window (TaqMan wants it short).
#' @param primer_lengths,probe_lengths primer / probe lengths to try.
#' @param primer_tm accepted flanking-primer Tm window.
#' @param probe_tm accepted probe Tm window (see the Tm-scale assumption above).
#' @param gc_range accepted GC window (primers and probes).
#' @param min_centrality least variant centrality within each probe (`[0, 1]`).
#' @param ... reserved for future design options.
#' @return a `kasp_taqman` list (`snp`, `primers`, `probes` - one row per allele -,
#'   `amplicon`, `compatibility` - the joint dimer/hairpin screen over the pair and
#'   both probes, and `readout` - the per-allele probe discrimination when
#'   OligoArrayAux is available, see [kasp_readout()]), or `NULL` (with a message)
#'   when no pair or probe passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_taqman_marker(genome, kasp_example("SNP")[1, ])
#' }
#' @family KASP TaqMan
#' @export
design_taqman_marker <- function(genome, snp, flank = 150L,
                                 min_product = 60L, max_product = 150L,
                                 primer_lengths = 18:25, probe_lengths = 18:25,
                                 primer_tm = KASP_LIMITS$com_tm,
                                 probe_tm = c(56, 72),
                                 gc_range = KASP_LIMITS$gc, min_centrality = 0.5, ...) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    cli::cli_inform(gettext("TaqMan allelic discrimination is implemented for single-base substitutions only."))
    return(NULL)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  si <- tmpl$snp_i
  primers <- design_flanking_pair(
    tmpl$seq, si, si, primer_lengths, min_product, max_product, primer_tm, gc_range
  )
  if (is.null(primers)) {
    cli::cli_inform(gettext("No short flanking primer pair passed the Tm/GC filters for TaqMan."))
    return(NULL)
  }
  # the probe must lie inside the amplicon (Wikipedia: it anneals within the
  # amplified region) and keep the variant away from its ends
  cand <- taqman_probe_candidates(tmpl$seq, si, probe_lengths, probe_tm, gc_range) |>
    dplyr::filter(
      .data$start >= primers$start, .data$end <= primers$end,
      .data$centrality >= min_centrality
    )
  if (!nrow(cand)) {
    cli::cli_inform(gettext("No composition- and position-clean hydrolysis probe passed the filters."))
    return(NULL)
  }
  # tier1: the most central mismatch (maximal destabilisation); tier2: the probe
  # being the hottest oligo in the tube (the ordering rule, see the header caveat)
  primer_hot <- max(primers$tm_f, primers$tm_r)
  cand$ordered <- cand$tm >= primer_hot
  best <- cand[order(-cand$centrality, -cand$ordered, -cand$tm), ][1L, ]
  ref_probe <- best$primer
  alt_probe <- as.character(Biostrings::replaceAt(
    Biostrings::DNAString(ref_probe),
    IRanges::IRanges(best$variant_pos, best$variant_pos), alt
  ))
  off <- snp$pos - si
  primers$start <- primers$start + off
  primers$end <- primers$end + off
  probes <- tibble::tibble(
    allele = c("REF", "ALT"),
    fluorophore = c("FAM", "HEX"),
    probe = c(ref_probe, alt_probe),
    start = snp$pos + best$start - si,
    end = snp$pos + best$end - si,
    tm = best$tm, gc = best$gc,
    variant_pos = best$variant_pos, centrality = best$centrality
  )
  out <- list(
    snp = snp, primers = primers, probes = probes,
    amplicon = tibble::tibble(start = primers$start, end = primers$end, product = primers$product)
  )
  # pair plus both allele-specific probes share the tube: joint screen
  out$compatibility <- tube_compatibility(c(
    primer_forward = out$primers$forward, primer_reverse = out$primers$reverse,
    probe_ref = out$probes$probe[1L], probe_alt = out$probes$probe[2L]
  ))
  # the readout this scheme *does* model: per-allele probe discrimination
  if (has_oligoarrayaux()) {
    out$readout <- probe_discrimination(
      tmpl$seq, si,
      tibble::tibble(
        allele = c("REF", "ALT"), probe = c(ref_probe, alt_probe),
        start = best$start, end = best$end, strand = "+"
      ),
      alt
    )
  }
  structure(out, class = "kasp_taqman")
}

#' @export
print.kasp_taqman <- function(x, ...) {
  cat("<kasp_taqman>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt, "\n", sep = "")
  cat("  amplicon:", x$amplicon$product, "bp\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  cat(glue::glue(
    "  {x$probes$allele} ({x$probes$fluorophore}): {x$probes$probe}  Tm {round(x$probes$tm, 1)} C, ",
    "centrality {round(x$probes$centrality, 2)}"
  ), sep = "\n")
  cat("\n")
  invisible(x)
}

#' TaqMan scheme definition (registry entry)
#'
#' The declarative row the shared scheme registry will consume: the assay's
#' discrimination/readout, its selection gates and ranking tiers, and the gate
#' labels/limits that the decision graph and text export render.  Kept here so the
#' TaqMan chemistry stays self-contained; it is discovered automatically (one file
#' per scheme, see `R/kasp_scheme_registry.R`).
#' @return a named list with the scheme row fields plus the gate/tier vocabulary.
#' @keywords internal
#' @noRd
scheme_def_taqman <- function() {
  list(
    scheme = "taqman",
    family = "taqman",
    label = "TaqMan 5'-nuclease allelic discrimination",
    discrimination = "allele-specific hydrolysis probes (FRET ratio)",
    readout = "co-dominant two-colour real-time",
    zygosity = "co-dominant",
    amplicon = "short (60-150 bp)",
    loci = 1L,
    verdict = "usable TaqMan assay",
    gates = c(
      "gate_taqman_amplicon", "gate_taqman_primer_tm", "gate_taqman_probe_tm",
      "gate_taqman_centrality", "gate_taqman_probe_comp"
    ),
    relax = c("gate_taqman_amplicon", "gate_taqman_centrality"),
    tiers = c(
      tier1 = "probe discrimination (variant centrality)",
      tier2 = "Tm ordering (probe above primers)",
      tier3 = "amplicon length"
    ),
    gate_short = c(
      gate_taqman_amplicon = "amplicon length (TaqMan)",
      gate_taqman_primer_tm = "primer Tm (TaqMan)",
      gate_taqman_probe_tm = "probe Tm",
      gate_taqman_centrality = "probe centrality",
      gate_taqman_probe_comp = "probe composition"
    ),
    gate_limit = c(
      gate_taqman_amplicon = "amplicon within [60, 150]",
      gate_taqman_primer_tm = "both primers Tm within [58, 63] C",
      gate_taqman_probe_tm = "probe Tm the hottest in the tube (ranking)",
      gate_taqman_centrality = "variant within the probe central window",
      gate_taqman_probe_comp = "GC 40-60%, no 5'-G, no run / SSR"
    ),
    gate_floor_label = c(
      gate_taqman_amplicon = "widen the amplicon window by 20 bp",
      gate_taqman_centrality = "allow the variant nearer a probe end"
    ),
    gate_hard_label = c(
      gate_taqman_primer_tm = "the flanking pair must sit in the extension window",
      gate_taqman_probe_tm = "a probe must out-compete the primers at annealing",
      gate_taqman_probe_comp = "a 5'-G or off-GC probe reads poorly"
    ),
    stage_label = "TaqMan 5'-nuclease"
  )
}
