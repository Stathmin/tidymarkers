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
#  primer pair plus ONE fluorescent probe that anneals *between* the primers and
#  is cleaved by the polymerase 5'->3' exonuclease during extension, separating a
#  5' reporter dye from a 3' quencher (Thermo Fisher, "How TaqMan Assays Work",
#  https://www.thermofisher.com/us/en/home/life-science/pcr/real-time-pcr/
#  real-time-pcr-learning-center/real-time-pcr-basics/how-taqman-assays-work.html).
#  Because the readout is end-point fluorescence, the assay is co-dominant, but
#  the design constraints are the probe's:
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
#' The digital-PCR analogue of [assemble_kasp_assay()]: one flanking pair plus a
#' single fluorescent hydrolysis probe that spans the variant between the
#' primers, so a droplet counts as positive when the target molecule is present.
#' The design optimises the probe's constraints - a short amplicon that still
#' leaves room for the probe, and a probe Tm about 8-10 C above the primers.
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
#' @return a `kasp_dpcr` list (`snp`, `primers`, `probe`, `amplicon`), or `NULL`
#'   (with a message) when no primer pair / probe passes the filters.
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
  out <- list(
    snp = snp,
    primers = primers,
    probe = list(
      sequence = probe$primer,
      start = probe$start + off, end = probe$end + off,
      tm = probe$tm, gc = probe$gc, strand = probe$strand
    ),
    amplicon = tibble::tibble(start = primers$start, end = primers$end)
  )
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
  cat(
    "  probe:", x$probe$sequence, "(Tm", round(x$probe$tm, 1), "C,",
    "probe-primer dTm", round(x$probe$tm - max(x$primers$tm_f, x$primers$tm_r), 1), "C)\n"
  )
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
    discrimination = "allele-specific hydrolysis probe",
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
