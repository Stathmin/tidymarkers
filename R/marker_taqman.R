# =============================================================================
#  marker_taqman.R  --  TaqMan 5'-nuclease allelic discrimination (probes)
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
#  primers") is therefore NOT a gate: on this compressed scale even a 30-mer probe
#  reaches only ~+4 C over the primers, so a fixed 8 C margin would doom every
#  design.  The ordering is a *ranking* preference only (tier2 here: the probe
#  should be the hottest oligo in the tube), while the physically binding rules
#  are gates: the probe must sit in its own absolute Tm window (or it never
#  anneals), it must be the only perfect probe site in the genome, its own
#  secondary structure must be acceptable, and the allele mismatch must actually
#  destabilise the probe duplex (the mismatch dG the method runs on).
#
#  Every declared gate is evaluated by the designer (C7): the candidate scans
#  enforce the hard gates and `grade_one_marker()` records them through the shared
#  relaxation ladder, so the scheme's decision graph and its engine read the same
#  thresholds.
# =============================================================================

# The two thresholds the formals do not carry.  A single mismatch that costs less
# than this much duplex stability (dG at the reaction temperature, kcal/mol) is not discriminated by the
# exonuclease step, so it is a hard gate; a probe this much more stable than the
# dG limit is a *soft* gate (a somewhat sticky probe still reports).
TAQMAN_MISMATCH_DG <- 1.5
TAQMAN_PROBE_STRUCT_RELAX <- 3
TAQMAN_CENTRALITY_RELAX <- -0.25
TAQMAN_AMPLICON_RELAX <- 20L

#' Probe candidates centred over a variant
#'
#' Reuses `primer_candidates()` for the composition screen (GC, Tm, runs, SSR),
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

#' Duplex dG cost of the probe's allele mismatch
#'
#' The physical discrimination of a 5'-nuclease probe: it is cleaved only on a
#' perfectly matched template, so the allele mismatch must cost the probe duplex
#' enough stability (`dG_off - dG_on`, kcal/mol; positive = destabilised).  The
#' REF probe is scored against its perfect complement and against the same target
#' carrying the ALT base at the variant.  Needs OligoArrayAux; `NA` otherwise.
#' @param probe the REF probe (plus-strand window).
#' @param variant_pos the variant's index within the probe.
#' @param alt the ALT base.
#' @return the mismatch dG (numeric scalar, `NA` when OligoArrayAux is absent).
#' @keywords internal
#' @noRd
taqman_mismatch_dg <- function(probe, variant_pos, alt) {
  if (!has_oligoarrayaux()) {
    return(NA_real_)
  }
  on <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(probe)))
  off <- on
  j <- nchar(probe) - variant_pos + 1L
  stringr::str_sub(off, j, j) <- as.character(Biostrings::complement(Biostrings::DNAString(alt)))
  oligo_dg(probe, off) - oligo_dg(probe, on)
}

#' Design a TaqMan 5'-nuclease allelic-discrimination marker at a variant
#'
#' A single flanking primer pair amplifies a short amplicon and **two
#' allele-specific hydrolysis probes** (one per allele, different fluorophores,
#' FAM/HEX here) are placed over the variant.  Each probe is centred over the
#' variant so the mismatch is maximally destabilising, avoids a 5'-G, sits inside
#' the amplicon and passes the GC/Tm composition filters; the pair itself reuses
#' the marker-agnostic flanking-pair machinery.  Single-base substitutions only.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param min_product,max_product the accepted amplicon window (TaqMan wants it short).
#' @param primer_lengths,probe_lengths primer / probe lengths to try.
#' @param primer_tm accepted flanking-primer Tm window.
#' @param probe_tm accepted probe Tm window (see the Tm-scale assumption above).
#' @param gc_range accepted GC window (primers and probes).
#' @param min_centrality least variant centrality within each probe (`[0, 1]`).
#' @param ... reserved for future design options.
#' @return a `marker_taqman` list (`snp`, `primers`, `probes` - one row per allele -,
#'   `amplicon`, `compatibility` - the joint dimer/hairpin screen over the pair and
#'   both probes, and `readout` - the per-allele probe discrimination when
#'   OligoArrayAux is available, see [marker_readout()]), or `NULL` (with a message)
#'   when no pair or probe passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_taqman_marker(genome, marker_example("SNP")[1, ])
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
  d <- scheme_def_taqman()
  amp <- d$amplicon_window
  # the amplicon rung of the ladder: the declared window first, widened by the
  # floor only when nothing clears it (so the soft gate is genuinely reachable)
  pair_in <- function(lo, hi) {
    design_flanking_pair(tmpl$seq, si, si, primer_lengths, lo, hi, primer_tm, gc_range)
  }
  primers <- pair_in(min_product, max_product) %||%
    pair_in(min_product - TAQMAN_AMPLICON_RELAX, max_product + TAQMAN_AMPLICON_RELAX)
  if (is.null(primers)) {
    cli::cli_inform(gettext("No short flanking primer pair passed the Tm/GC filters for TaqMan."))
    return(NULL)
  }
  # the probe must lie strictly between the primer 3' ends (not merely inside the
  # amplicon): overlapping either primer lets the primer out-compete the probe.
  # The inner bounds are the ones pick_dpcr_probe() uses for the dPCR sibling,
  # and the search still dips below the centrality gate by its floor so the ladder
  # can report the relaxation; the ranking prefers the most central probe.
  inner_lo <- primers$start + nchar(primers$forward)
  inner_hi <- primers$end - nchar(primers$reverse)
  if (inner_hi - inner_lo + 1L < min(probe_lengths)) {
    cli::cli_inform(gettext("No composition- and position-clean hydrolysis probe passed the filters."))
    return(NULL)
  }
  cand <- taqman_probe_candidates(tmpl$seq, si, probe_lengths, probe_tm, gc_range) |>
    dplyr::filter(
      .data$start >= inner_lo, .data$end <= inner_hi,
      .data$centrality >= min_centrality + TAQMAN_CENTRALITY_RELAX
    )
  if (!nrow(cand)) {
    cli::cli_inform(gettext("No composition- and position-clean hydrolysis probe passed the filters."))
    return(NULL)
  }
  # tier1: the most central mismatch (maximal destabilisation); tier2: the probe
  # being the hottest oligo in the tube (the ordering rule, see the header caveat)
  primer_hot <- max(primers$tm_f, primers$tm_r)
  cand$ordered <- cand$tm >= primer_hot
  cand$mismatch_dg <- vapply(
    seq_len(nrow(cand)),
    function(i) taqman_mismatch_dg(cand$primer[i], cand$variant_pos[i], alt),
    numeric(1)
  )
  best <- cand[order(-cand$centrality, -cand$ordered, -cand$mismatch_dg, -cand$tm), ][1L, ]
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
  ), conditions = scheme_conditions("taqman"))

  # ---- C7: the declared gates are all evaluated on the winner ----------------
  probe_gc_obs <- gc_frac(ref_probe)
  comp_ok <- substr(ref_probe, 1L, 1L) != "G" &&
    probe_gc_obs >= gc_range[1L] && probe_gc_obs <= gc_range[2L] &&
    !has_run(ref_probe) && !has_repeat(ref_probe)
  if (has_oligoarrayaux()) {
    hairpin <- oligo_dg_vec(ref_probe)
    homodimer <- oligo_dg_vec(ref_probe, ref_probe)
  } else {
    hairpin <- NA_real_
    homodimer <- NA_real_
  }
  primer_tm_all <- c(out$primers$tm_f, out$primers$tm_r)
  # Both allele-specific probes must have a single perfect genomic site.  The REF
  # probe's on-target site is in the reference genome; the ALT probe's on-target
  # site is the ALT allele, which the reference genome does not carry, so any
  # reference-genome hit is an off-target (a second site).  Count each probe in
  # its own allele context (the degenerate alt == genome base case makes the two
  # probes identical, so the on-target site is already counted).
  ref_hits <- exact_oligo_hits(ref_probe, genome)
  alt_hits <- exact_oligo_hits(alt_probe, genome)
  alt_sites <- if (alt_probe == ref_probe) alt_hits else alt_hits + 1L
  probes_unique <- ref_hits == 1L && alt_sites == 1L
  graded <- tibble::tibble(
    m_probe_comp = if (comp_ok) 1 else -1,
    m_probe_specific = if (probes_unique) 1 else -1,
    m_probe_structure = pmin(hairpin - KASP_LIMITS$hairpin_dg, homodimer - KASP_LIMITS$homodimer_dg),
    m_probe_tm = min(out$probes$tm[1L] - d$probe_tm_window[1L], d$probe_tm_window[2L] - out$probes$tm[1L]),
    m_primer_tm = min(primer_tm_all - d$primer_tm_window[1L], d$primer_tm_window[2L] - primer_tm_all),
    m_discrimination = best$mismatch_dg - TAQMAN_MISMATCH_DG,
    m_centrality = best$centrality - min_centrality,
    m_amplicon = pmin(out$amplicon$product - amp[1L], amp[2L] - out$amplicon$product)
  )
  gate_margin <- c(
    gate_taqman_probe_comp = "m_probe_comp", gate_taqman_probe_specific = "m_probe_specific",
    gate_taqman_probe_structure = "m_probe_structure", gate_taqman_probe_tm = "m_probe_tm",
    gate_taqman_primer_tm = "m_primer_tm", gate_taqman_discrimination = "m_discrimination",
    gate_taqman_centrality = "m_centrality", gate_taqman_amplicon = "m_amplicon"
  )
  sel <- grade_one_marker(
    graded,
    gates = d$gates, relax_order = d$relax,
    tiers = c("m_discrimination", "m_probe_tm", "m_amplicon"), # tier1..3 in scheme_def_taqman()
    gate_margin = gate_margin, gate_floor = d$gate_floor
  )
  if (is.null(sel)) {
    cli::cli_inform(gettext(
      "The TaqMan candidate failed a hard gate (probe composition, specificity, Tm window, primer Tm or mismatch dG)."
    ))
    return(NULL)
  }
  out$selection <- sel

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
  structure(out, class = "marker_taqman")
}

#' @export
print.marker_taqman <- function(x, ...) {
  cat("<marker_taqman>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt, "\n", sep = "")
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
#' per scheme, see `R/marker_scheme_registry.R`).
#' @return a named list with the scheme row fields plus the gate/tier vocabulary.
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_taqman
#' @export
#' @noRd
marker_map_spec.marker_taqman <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  pr <- marker$probes
  sequences <- list(forward = marker$primers$forward, reverse = marker$primers$reverse)
  if (!is.null(pr) && all(c("allele", "probe") %in% names(pr))) {
    for (i in seq_len(nrow(pr))) {
      sequences[[paste0("probe ", pr$allele[i])]] <- pr$probe[i]
    }
  } else if (!is.null(pr) && "probe" %in% names(pr)) {
    sequences$probe <- pr$probe[1L]
  }
  map_spec(
    marker,
    features = c(map_pair_rows(marker, sp[["start"]], sp[["end"]]), map_probe_rows(pr)),
    sequences = sequences
  )
}

scheme_def_taqman <- function() {
  # the windows the designer applies, read from its formals so the two cannot
  # drift (formals() gives the *default expression*, hence eval())
  fm <- formals(design_taqman_marker)
  amp <- c(eval(fm$min_product), eval(fm$max_product))
  primer_tm <- eval(fm$primer_tm)
  probe_tm <- eval(fm$probe_tm)
  gc <- eval(fm$gc_range)
  cent <- eval(fm$min_centrality)
  gates <- c(
    "gate_taqman_probe_comp", "gate_taqman_probe_specific", "gate_taqman_probe_structure",
    "gate_taqman_probe_tm", "gate_taqman_primer_tm", "gate_taqman_discrimination",
    "gate_taqman_centrality", "gate_taqman_amplicon"
  )
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
    gates = gates,
    # least grave first: allow the mismatch nearer a probe end, then a longer
    # amplicon, and only then a slightly stickier probe
    relax = c("gate_taqman_centrality", "gate_taqman_amplicon", "gate_taqman_probe_structure"),
    tiers = c(
      tier1 = "probe discrimination (mismatch dG margin)",
      tier2 = "Tm ordering (probe above primers)",
      tier3 = "amplicon length (near the ideal)"
    ),
    gate_short = c(
      gate_taqman_probe_comp = "probe composition",
      gate_taqman_probe_specific = "probe specificity",
      gate_taqman_probe_structure = "probe self-structure",
      gate_taqman_probe_tm = "probe Tm window",
      gate_taqman_primer_tm = "primer Tm (TaqMan)",
      gate_taqman_discrimination = "probe mismatch dG",
      gate_taqman_centrality = "probe centrality",
      gate_taqman_amplicon = "amplicon length (TaqMan)"
    ),
    gate_limit = c(
      gate_taqman_probe_comp = glue::glue("GC {gc[1L] * 100}-{gc[2L] * 100} %, no 5' G, no run / SSR"),
      gate_taqman_probe_specific = "exactly one perfect probe site in the genome",
      gate_taqman_probe_structure = glue::glue("probe hairpin/homodimer dG (reaction T) >= {KASP_LIMITS$hairpin_dg}"),
      gate_taqman_probe_tm = glue::glue("probe Tm within [{probe_tm[1L]}, {probe_tm[2L]}] C"),
      gate_taqman_primer_tm = glue::glue("both primers Tm within [{primer_tm[1L]}, {primer_tm[2L]}] C"),
      gate_taqman_discrimination = glue::glue("mismatch dG >= {TAQMAN_MISMATCH_DG} kcal/mol"),
      gate_taqman_centrality = glue::glue("variant centrality >= {cent}"),
      gate_taqman_amplicon = glue::glue("amplicon within [{amp[1L]}, {amp[2L]}] bp")
    ),
    # the numeric windows the engine scores against (so gate margins and the
    # declared limits cannot drift), and the ladder's floors in margin units
    amplicon_window = amp,
    primer_tm_window = primer_tm,
    probe_tm_window = probe_tm,
    gc_window = gc,
    min_centrality = cent,
    gate_floor = c(
      gate_taqman_centrality = TAQMAN_CENTRALITY_RELAX,
      gate_taqman_amplicon = -TAQMAN_AMPLICON_RELAX,
      gate_taqman_probe_structure = -TAQMAN_PROBE_STRUCT_RELAX
    ),
    gate_floor_label = c(
      gate_taqman_centrality = "allow the variant nearer a probe end",
      gate_taqman_amplicon = glue::glue("widen the amplicon window by {TAQMAN_AMPLICON_RELAX} bp"),
      gate_taqman_probe_structure = glue::glue("allow {TAQMAN_PROBE_STRUCT_RELAX} kcal/mol more stable")
    ),
    gate_hard_label = c(
      gate_taqman_probe_comp = "a 5'-G or off-GC probe reads poorly",
      gate_taqman_probe_specific = "a probe that binds elsewhere is cleaved there and lights the wrong dye",
      gate_taqman_probe_tm = "a probe outside its Tm window does not anneal, so there is no signal",
      gate_taqman_primer_tm = "the flanking pair must sit in the extension window",
      gate_taqman_discrimination = "without a destabilising mismatch the probe cannot tell REF from ALT"
    ),
    stage_label = "TaqMan 5'-nuclease"
  )
}
