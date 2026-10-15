# =============================================================================
#  kasp_dpcr.R  --  dPCR (droplet digital PCR) as a marker scheme
# =============================================================================
#  dPCR quantifies a target by *partitioning* the reaction into thousands of
#  independent microreactors (droplets) and counting the fluorescence-positive
#  ones; the ratio of positive partitions gives the concentration directly by
#  Poisson statistics, so no standard curve is needed (Quan, Sauzade & Brouzes,
#  "dPCR: A Technology Review", Sensors 2018, 18(4):1271,
#  doi:10.3390/s18041271).  The dPCR data-reporting conventions are the dMIQE
#  guidelines (Huggett et al., Clin Chem 2013, 59:892-902), referenced there.
#
#  This scheme models the hydrolysis (TaqMan-style) probe assay: ONE flanking
#  primer pair plus TWO allele-specific fluorescent probes (one per allele,
#  different dyes - FAM/HEX) that anneal *between* the primers and are cleaved by
#  the polymerase 5'->3' exonuclease during extension, separating a 5' reporter
#  dye from a 3' quencher (Thermo Fisher, "How TaqMan Assays Work",
#  https://www.thermofisher.com/us/en/home/life-science/pcr/real-time-pcr/
#  real-time-pcr-learning-center/real-time-pcr-basics/how-taqman-assays-work.html).
#  A single generic probe would only report presence/absence of the locus, not the
#  allele; two allele-specific probes make the call co-dominant, as dMIQE/standard
#  ddPCR SNP assays require.
#
#    * amplicon kept SHORT (the probe must fit between the primers) - ~60-120 bp,
#      typically <= 150 bp;
#    * probe Tm ~8-10 C ABOVE the primers so the probe out-competes them at the
#      annealing step;
#    * probe 18-25 nt, GC 40-60%, no G at the 5' end (a 5' guanine quenches the
#      reporter even after cleavage) and no homopolymer run longer than 4.  We
#      reuse the primer run screen, which bans runs of 4 or more - at least as
#      strict as the heuristic.
#
#  The probe/amplicon numbers above are the conventional TaqMan design
#  heuristics: the vendor design pages could not be fetched here (HTTP 403 /
#  redirect loops), so the thresholds are recorded as ASSUMPTIONS and exposed as
#  arguments rather than hard-wired.  The designer reuses the marker-agnostic
#  primitives [snp_template()], [design_flanking_pair()], [primer_candidates()]
#  and the composition tests behind them, plus [has_palindrome()].
# =============================================================================

#' Design a dPCR hydrolysis-probe marker at a variant
#'
#' The digital-PCR analogue of [assemble_kasp_assay()]: one flanking pair plus
#' **two allele-specific hydrolysis probes** (one per allele, different dyes -
#' FAM/HEX) that span the variant between the primers, so a droplet's dye colour
#' reads the allele (co-dominant).  A single generic probe would only report the
#' locus's presence, not the genotype.  The design optimises the probe's
#' constraints - a short amplicon that still leaves room for the probe, and a
#' probe Tm about 8-10 C above the primers - then swaps the variant base to build
#' the ALT probe from the REF one.
#'
#' The probe is a hydrolysis (TaqMan-style) probe: the polymerase 5'->3'
#' exonuclease cleaves it during extension, releasing a 5' reporter from a 3'
#' quencher (Thermo Fisher).  The conventional probe composition rules (18-25 nt,
#' GC 40-60%, no 5' G, run < 4) are exposed as arguments; the Tm-above-primers
#' rule is `probe_tm_margin`, and the amplicon window is
#' `min_product`/`max_product`.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param min_product,max_product the accepted amplicon window (dPCR wants short).
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows (the dPCR annealing window).
#' @param probe_lengths probe lengths to try.
#' @param probe_gc accepted probe GC window.
#' @param probe_tm_margin minimal probe Tm above the primers (the ~8-10 C rule;
#'   8 is used as the floor).
#' @return a `kasp_dpcr` list (`snp`, `primers`, `probes` - one allele-specific
#'   probe per allele - `amplicon`, `compatibility` - the joint dimer/hairpin
#'   screen over the pair and both probes), or `NULL` (with a message) when no
#'   primer pair / probe passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_dpcr_marker(genome, kasp_example("SNP")[1, ])
#' }
#' @family KASP dPCR
#' @export
design_dpcr_marker <- function(genome, snp, flank = 150L,
                               min_product = 60L, max_product = 150L, lengths = 18:25,
                               tm_range = c(44, 50), gc_range = c(0.40, 0.60),
                               probe_lengths = 18:25, probe_gc = c(0.40, 0.60),
                               probe_tm_margin = 8) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    cli::cli_inform(gettext("dPCR allelic discrimination is implemented for single-base substitutions only."))
    return(NULL)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  si <- tmpl$snp_i
  # the probe spans the variant between the primers, so keep the pair at least
  # the longest probe away on EACH side - design_flanking_pair()'s lo/hi gap
  # reserves that room; the product window is the amplicon-length gate
  pad <- max(probe_lengths)
  primers <- design_flanking_pair(
    tmpl$seq, si - pad, si + pad, lengths, min_product, max_product, tm_range, gc_range
  )
  if (is.null(primers)) {
    cli::cli_inform(gettext("No flanking primer pair passed the Tm/GC filters for dPCR."))
    return(NULL)
  }
  probe <- pick_dpcr_probe(tmpl$seq, si, primers, probe_lengths, probe_gc, probe_tm_margin)
  if (is.null(probe)) {
    cli::cli_inform(gettext("No hydrolysis probe passed the composition/Tm filters for dPCR."))
    return(NULL)
  }
  off <- snp$pos - si
  primers$start <- primers$start + off
  primers$end <- primers$end + off
  # the selected probe is built from the reference; the ALT probe is the same
  # window with the variant base swapped, so the two dyes report the allele
  vpos <- if (identical(probe$strand, "+")) si - probe$start + 1L else probe$end - si + 1L
  alt_base <- if (identical(probe$strand, "+")) {
    alt
  } else {
    as.character(Biostrings::complement(Biostrings::DNAString(alt)))
  }
  alt_probe <- as.character(Biostrings::replaceAt(
    Biostrings::DNAString(probe$primer),
    IRanges::IRanges(vpos, vpos), alt_base
  ))
  # hoist the scalar probe fields: a tibble() column evaluates sequentially, so a
  # column named `probe` would shadow the `probe` row for the later columns
  p_start <- probe$start + off
  p_end <- probe$end + off
  p_tm <- probe$tm
  p_gc <- probe$gc
  p_strand <- probe$strand
  probes <- tibble::tibble(
    allele = c("REF", "ALT"),
    fluorophore = c("FAM", "HEX"),
    probe = c(probe$primer, alt_probe),
    start = p_start, end = p_end,
    tm = p_tm, gc = p_gc, strand = p_strand, variant_pos = vpos
  )
  out <- list(
    snp = snp,
    primers = primers,
    probes = probes,
    amplicon = tibble::tibble(start = primers$start, end = primers$end)
  )
  # pair plus both allele-specific probes share the droplet: joint screen
  out$compatibility <- tube_compatibility(c(
    primer_forward = out$primers$forward, primer_reverse = out$primers$reverse,
    probe_ref = out$probes$probe[1L], probe_alt = out$probes$probe[2L]
  ))
  structure(out, class = "kasp_dpcr")
}

#' The best hydrolysis probe inside a primer pair
#'
#' Screens both strands for probes of `probe_lengths` in the inner region between
#' the primer 3' ends with [primer_candidates()] (GC/Tm window, homopolymer runs,
#' SSRs), keeps those that span the variant, do not start with G and reach
#' `max(primer Tm) + probe_tm_margin` (the ~8-10 C rule), then prefers a
#' non-palindromic, variant-centred, short probe.
#' @keywords internal
#' @noRd
pick_dpcr_probe <- function(seq, si, primers, probe_lengths, probe_gc, probe_tm_margin) {
  inner_lo <- primers$start + nchar(primers$forward)
  inner_hi <- primers$end - nchar(primers$reverse)
  if (inner_hi - inner_lo + 1L < min(probe_lengths)) {
    return(NULL)
  }
  tm_lo <- max(primers$tm_f, primers$tm_r) + probe_tm_margin
  starts <- seq.int(inner_lo, inner_hi)
  cand <- purrr::map(c("+", "-"), function(st) {
    primer_candidates(seq, starts, probe_lengths, st, c(tm_lo, tm_lo + 20), probe_gc) |>
      dplyr::filter(.data$start >= inner_lo, .data$end <= inner_hi) |>
      dplyr::mutate(strand = st)
  }) |>
    purrr::list_rbind() |>
    dplyr::filter(.data$start <= si, .data$end >= si) |>
    dplyr::filter(stringr::str_sub(.data$primer, 1L, 1L) != "G") |>
    dplyr::mutate(
      palindrome = has_palindrome(rlang::set_names(.data$primer, seq_along(.data$primer))),
      offset = abs((.data$start + .data$end) / 2 - si)
    )
  if (!nrow(cand)) {
    return(NULL)
  }
  dplyr::arrange(cand, .data$palindrome, .data$offset, .data$len)[1L, ]
}

#' @export
print.kasp_dpcr <- function(x, ...) {
  cat("<kasp_dpcr>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt, "\n", sep = "")
  cat("  amplicon:", x$primers$product, "bp\n")
  cat("  F:", x$primers$forward, "(Tm", round(x$primers$tm_f, 1), "C)\n")
  cat("  R:", x$primers$reverse, "(Tm", round(x$primers$tm_r, 1), "C)\n")
  cat(glue::glue(
    "  {x$probes$allele} ({x$probes$fluorophore}): {x$probes$probe} (Tm {round(x$probes$tm, 1)} C)"
  ), sep = "\n")
  cat("\n")
  invisible(x)
}

#' The dPCR scheme definition
#'
#' The registry record for the dPCR hydrolysis-probe chemistry, in the shape the
#' shared scheme registry consumes (see [kasp_schemes()] and [kasp_criteria()]).
#' Four gates: the amplicon window, the primer Tm window, the probe-above-primer
#' Tm ordering, and the probe composition (GC / 5'-G / homopolymer).  Only the
#' amplicon window is relaxable; the Tm ordering and the probe composition are
#' chemistry, not preference.  Tiers rank probe/amplicon quality, then the Tm
#' ordering, then the amplicon length.
#' @return a list describing the `dpcr` scheme.
#' @keywords internal
#' @noRd
scheme_def_dpcr <- function() {
  gates <- c("gate_dpcr_amplicon", "gate_dpcr_primer_tm", "gate_dpcr_probe_tm", "gate_dpcr_probe_comp")
  list(
    scheme = "dpcr",
    family = "dpcr",
    label = "dPCR hydrolysis-probe (droplet digital)",
    discrimination = "allele-specific hydrolysis probes (FAM/HEX)",
    readout = "end-point droplet fluorescence (Poisson counting)",
    zygosity = "co-dominant",
    amplicon = "short (<= 150 bp)",
    loci = 1L,
    verdict = "usable dPCR assay",
    gates = gates,
    relax = "gate_dpcr_amplicon",
    tiers = c(
      tier1 = "probe / amplicon quality",
      tier2 = "Tm ordering (probe above primers)",
      tier3 = "amplicon length"
    ),
    gate_short = c(
      gate_dpcr_amplicon = "amplicon length (dPCR)",
      gate_dpcr_primer_tm = "primer Tm (dPCR)",
      gate_dpcr_probe_tm = "probe Tm > primer Tm",
      gate_dpcr_probe_comp = "probe composition"
    ),
    gate_limit = c(
      gate_dpcr_amplicon = "amplicon within [60, 150] bp",
      gate_dpcr_primer_tm = "both primers Tm within [lo, hi]",
      gate_dpcr_probe_tm = "probe Tm >= primer Tm + 8 C",
      gate_dpcr_probe_comp = "GC 40-60%, no 5' G, run < 4"
    ),
    gate_floor_label = c(gate_dpcr_amplicon = "widen the amplicon window"),
    gate_hard_label = c(
      gate_dpcr_primer_tm = "both primers must anneal in the extension window",
      gate_dpcr_probe_tm = "the probe must out-compete the primers",
      gate_dpcr_probe_comp = "a 5' G quenches the reporter and a GC drift kills hybridisation"
    ),
    stage_label = "dPCR (droplet digital PCR)"
  )
}
