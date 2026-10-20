# =============================================================================
#  marker_dpcr.R  --  dPCR (droplet digital PCR) as a marker scheme
# =============================================================================
#  dPCR quantifies a target by *partitioning* the reaction into thousands of
#  independent microreactors (droplets) and counting the fluorescence-positive
#  ones; the ratio of positive partitions gives the concentration directly by
#  Poisson statistics, so no standard curve is needed (Quan, Sauzade & Brouzes,
#  "dPCR: A Technology Review", Sensors 2018, 18(4):1271,
#  doi:10.3390/s18041271).  The dPCR data-reporting conventions are the dMIQE
#  guidelines (Huggett et al., Clin Chem 2013, 59:892-902), referenced there.
#  NB: this package designs *constructs*; it does not simulate the
#  partitioning/Poisson readout (see [marker_readout()], whose assumption note for
#  this scheme says so), so the scheme must not *claim* a counting model it
#  never computes.
#
#  This scheme models the hydrolysis (TaqMan-style) probe assay: ONE flanking
#  primer pair plus TWO allele-specific fluorescent probes (one per allele,
#  different dyes - FAM/HEX) that anneal *between* the primers and are cleaved by
#  the polymerase 5'->3' exonuclease during extension, separating a 5' reporter
#  dye from a 3' quencher (Thermo Fisher, "How TaqMan Assays Work",
#  https://www.thermofisher.com/us/en/home/life-science/pcr/real-time-pcr/
#  real-time-pcr-learning-center/real-time-pcr-basics/how-taqman-assays-work.html).
#  REF and ALT are told apart by *which dye* lights the droplet: only the probe
#  whose sequence matches the template is cleaved, so the *probe* - never a
#  primer 3' base - carries the discrimination.  A single generic probe would
#  only report presence/absence of the locus, not the allele; two allele-specific
#  probes make the call co-dominant, as dMIQE/standard ddPCR SNP assays require.
#
#    * amplicon kept SHORT (the probe must fit between the primers) - ~60-120 bp,
#      typically <= 150 bp;
#    * probe Tm ~8-10 C ABOVE the primers so the probe out-competes them at the
#      annealing step;
#    * probe 18-25 nt, GC 40-60%, no G at the 5' end (a 5' guanine quenches the
#      reporter even after cleavage) and no homopolymer run longer than 4.  We
#      reuse the primer run screen, which bans runs of 4 or more - at least as
#      strict as the heuristic;
#    * the probe must be *unique*: a second perfect binding site anywhere in the
#      genome would be cleaved too and light the wrong droplet, so the probe is
#      screened for a single perfect genomic site;
#    * the probe's own secondary structure (hairpin / homodimer) slows
#      hybridisation exactly as it does for a primer, so it is gated on the same
#      OligoArrayAux dG37 limits ([KASP_LIMITS]); it is a *soft* gate (a
#      somewhat sticky probe still reports, just less efficiently).
#
#  The probe/amplicon numbers above are the conventional TaqMan design
#  heuristics: the vendor design pages could not be fetched here (HTTP 403 /
#  redirect loops), so the thresholds are recorded as ASSUMPTIONS and exposed as
#  arguments rather than hard-wired.  The designer reuses the marker-agnostic
#  primitives [snp_template()], [design_flanking_pair()], [primer_candidates()]
#  and the composition tests behind them, plus [has_palindrome()].
#
#  Every declared gate is evaluated by the designer (C7): the candidate scans
#  below *enforce* the hard gates (composition, Tm windows, probe Tm ordering,
#  amplicon window) and [select_with_ladder()] then re-scores the winner, so the
#  scheme's decision graph and its engine read the same thresholds.
# =============================================================================

# How far the ladder may lower the two *soft* dPCR gates (margin units): the
# amplicon window may widen by this many bp, and a probe this much more stable
# than the dG37 limit is still tolerated.  Declared here and consumed by both the
# designer and scheme_def_dpcr(), so the two cannot drift.
DPCR_AMPLICON_RELAX <- 20L
DPCR_PROBE_STRUCT_RELAX <- 3

#' Exact genome-wide hit count of a short oligo (both orientations)
#'
#' A probe or extension primer can bind either strand, so its perfect genomic
#' sites are the exact occurrences of the oligo *and* of its reverse complement;
#' counting both gives the on-target site (exactly one) plus any off-target.  The
#' hydrolysis-probe and single-base-extension schemes handled together
#' (dpcr / taqman / snapshot) use this for their specificity gates, where a
#' 3'-seed count is the wrong question: the whole probe must anneal uniquely, not
#' just its 3' end.
#' @param oligo the probe / extension primer (ACGT).
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @return the number of perfect genomic binding sites (integer).
#' @keywords internal
#' @noRd
exact_oligo_hits <- function(oligo, genome) {
  # keep the two orientations as plain strings: c()/unique() on XString objects
  # recurse in this Biostrings build
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(oligo)))
  pats <- unique(c(oligo, rc))
  sum(vapply(genome_chroms(genome), function(ch) {
    subj <- genome_seq(genome, ch)
    sum(vapply(pats, function(p) as.integer(Biostrings::countPattern(p, subj, max.mismatch = 0L)), integer(1)))
  }, integer(1)))
}

#' Score a scheme's declared gates on one graded candidate (the C7 recorder)
#'
#' The scheme designers pick their marker with the usual candidate scans, whose
#' filters *are* the gates; this runs the same gates through the shared
#' `select_with_ladder()` so the winner's gate outcomes are recorded in the
#' marker (`$selection`), and a candidate that cannot clear a *hard* gate is
#' rejected (`NULL`) rather than returned silently.  The graded table holds one
#' numeric margin column per gate (a value `>= 0` passes at the original
#' threshold); a binary gate is encoded as `1` / `-1`.
#' @param graded a one-row tibble of margin columns (`m_*`).
#' @param gates the scheme's gate ids, in application order.
#' @param relax_order the soft gate ids, least grave first.
#' @param tiers the tier margin columns, best-first (lexicographic).
#' @param gate_margin named map `gate id -> margin column`.
#' @param gate_floor named map `gate id -> floor` (margin units, negative).
#' @return `list(gates, relaxed, violated)` or `NULL` when a hard gate blocks.
#' @keywords internal
#' @noRd
grade_one_marker <- function(graded, gates, relax_order, tiers, gate_margin, gate_floor) {
  sel <- select_with_ladder(
    graded,
    gates = gates, relax_order = relax_order, tiers = tiers, ties = character(0),
    gate_margin = gate_margin, gate_floor = gate_floor
  )
  if (is.null(sel$best)) {
    return(NULL)
  }
  pass <- vapply(gates, function(g) {
    m <- graded[[gate_margin[[g]]]]
    is.na(m) || m >= 0
  }, logical(1))
  list(gates = pass, relaxed = sel$relaxed, violated = sel$violated)
}

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
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param min_product,max_product the accepted amplicon window (dPCR wants short).
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows (the dPCR annealing window).
#' @param probe_lengths probe lengths to try.
#' @param probe_gc accepted probe GC window.
#' @param probe_tm_margin minimal probe Tm above the primers (the ~8-10 C rule;
#'   8 is used as the floor).
#' @return a `marker_dpcr` list (`snp`, `primers`, `probes` - one allele-specific
#'   probe per allele - `amplicon`, `compatibility` - the joint dimer/hairpin
#'   screen over the pair and both probes, and `readout` - the per-allele probe
#'   discrimination when OligoArrayAux is available, see [marker_readout()]), or
#'   `NULL` (with a message) when no primer pair / probe passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_dpcr_marker(genome, marker_example("SNP")[1, ])
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
  d <- scheme_def_dpcr()
  # the probe spans the variant between the primers, so keep the pair at least
  # the longest probe away on EACH side - design_flanking_pair()'s lo/hi gap
  # reserves that room; the product window is the amplicon-length gate
  pad <- max(probe_lengths)
  pair_in <- function(lo, hi) {
    design_flanking_pair(tmpl$seq, si - pad, si + pad, lengths, lo, hi, tm_range, gc_range)
  }
  # the amplicon rung of the ladder: the declared window first, widened by the
  # floor only when nothing clears it (so the soft gate is genuinely reachable)
  primers <- pair_in(min_product, max_product) %||%
    pair_in(min_product - DPCR_AMPLICON_RELAX, max_product + DPCR_AMPLICON_RELAX)
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

  # ---- C7: the declared gates are all evaluated on the winner ----------------
  # The scans above already *enforce* the hard gates (Tm/GC windows, the 5'-G and
  # run/repeat screens, the probe-Tm margin, the product window); record them and
  # let the ladder report any soft bound it had to bend.
  amp <- d$amplicon_window
  prim_tm_all <- c(out$primers$tm_f, out$primers$tm_r)
  probe_gc_obs <- gc_frac(probe$primer)
  comp_ok <- substr(probe$primer, 1L, 1L) != "G" &&
    probe_gc_obs >= probe_gc[1L] && probe_gc_obs <= probe_gc[2L] &&
    !has_run(probe$primer) && !has_repeat(probe$primer)
  if (has_oligoarrayaux()) {
    hairpin <- oligo_dg_vec(probe$primer)
    homodimer <- oligo_dg_vec(probe$primer, probe$primer)
  } else {
    hairpin <- NA_real_
    homodimer <- NA_real_
  }
  # Both allele-specific probes must have a single perfect genomic site.  The REF
  # probe's on-target site is in the reference genome; the ALT probe's on-target
  # site is the ALT allele, which the reference genome does not carry, so any
  # reference-genome hit is an off-target (a second site).  Count each probe in
  # its own allele context (the degenerate alt == genome base case makes the two
  # probes identical, so the on-target site is already counted).
  ref_hits <- exact_oligo_hits(probe$primer, genome)
  alt_hits <- exact_oligo_hits(alt_probe, genome)
  alt_sites <- if (alt_probe == probe$primer) alt_hits else alt_hits + 1L
  probes_unique <- ref_hits == 1L && alt_sites == 1L
  graded <- tibble::tibble(
    m_probe_comp = if (comp_ok) 1 else -1,
    m_probe_specific = if (probes_unique) 1 else -1,
    m_probe_structure = pmin(hairpin - KASP_LIMITS$hairpin_dg, homodimer - KASP_LIMITS$homodimer_dg),
    m_probe_tm = out$probes$tm[1L] - (max(prim_tm_all) + probe_tm_margin),
    m_primer_tm = min(prim_tm_all - d$primer_tm_window[1L], d$primer_tm_window[2L] - prim_tm_all),
    m_amplicon = pmin(out$primers$product - amp[1L], amp[2L] - out$primers$product)
  )
  gate_margin <- c(
    gate_dpcr_probe_comp = "m_probe_comp", gate_dpcr_probe_specific = "m_probe_specific",
    gate_dpcr_probe_structure = "m_probe_structure", gate_dpcr_probe_tm = "m_probe_tm",
    gate_dpcr_primer_tm = "m_primer_tm", gate_dpcr_amplicon = "m_amplicon"
  )
  sel <- grade_one_marker(
    graded,
    gates = d$gates, relax_order = d$relax,
    tiers = c("m_probe_tm", "m_probe_structure", "m_amplicon"), # tier1..3 in scheme_def_dpcr()
    gate_margin = gate_margin, gate_floor = d$gate_floor
  )
  if (is.null(sel)) {
    cli::cli_inform(gettext(
      "The dPCR candidate failed a hard gate (probe composition, specificity, Tm ordering or primer Tm)."
    ))
    return(NULL)
  }
  out$selection <- sel

  # the readout this scheme *does* model: per-allele probe discrimination
  if (has_oligoarrayaux()) {
    p_start <- probe$start
    p_end <- probe$end
    p_strand <- probe$strand
    p_ref <- probe$primer
    out$readout <- probe_discrimination(
      tmpl$seq, si,
      tibble::tibble(
        allele = c("REF", "ALT"), probe = c(p_ref, alt_probe),
        start = p_start, end = p_end, strand = p_strand
      ),
      alt
    )
  }
  structure(out, class = "marker_dpcr")
}

#' The best hydrolysis probe inside a primer pair
#'
#' Screens both strands for probes of `probe_lengths` in the inner region between
#' the primer 3' ends with `primer_candidates()` (GC/Tm window, homopolymer runs,
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
print.marker_dpcr <- function(x, ...) {
  cat("<marker_dpcr>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt, "\n", sep = "")
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
#' shared scheme registry consumes (see [marker_schemes()] and [marker_criteria()]).
#' Six gates: the probe composition (5'-G / GC / run), the probe's genomic
#' specificity, the probe's self-structure, the probe-above-primer Tm ordering,
#' the primer Tm window, and the amplicon window.  Only the amplicon window and
#' the probe's self-structure are relaxable - the composition, the Tm ordering
#' and the single perfect probe site are chemistry, not preference.  Tiers rank
#' the probe's Tm margin, then its self-structure, then the amplicon length.
#' @return a list describing the `dpcr` scheme.
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_dpcr
#' @export
#' @noRd
marker_map_spec.marker_dpcr <- function(marker, scheme = NULL, ...) {
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

scheme_def_dpcr <- function() {
  # the windows the designer applies, read from its formals so the two cannot
  # drift (formals() gives the *default expression*, hence eval())
  fm <- formals(design_dpcr_marker)
  amp <- c(eval(fm$min_product), eval(fm$max_product))
  tm <- eval(fm$tm_range)
  gc <- eval(fm$probe_gc)
  pmar <- eval(fm$probe_tm_margin)
  gates <- c(
    "gate_dpcr_probe_comp", "gate_dpcr_probe_specific", "gate_dpcr_probe_structure",
    "gate_dpcr_probe_tm", "gate_dpcr_primer_tm", "gate_dpcr_amplicon"
  )
  list(
    scheme = "dpcr",
    family = "dpcr",
    label = "dPCR hydrolysis-probe (droplet digital)",
    discrimination = "allele-specific hydrolysis probes (FAM/HEX)",
    # the construct is designed here; the droplet partition/Poisson count is NOT
    # simulated (see marker_readout()), so do not claim it in the readout
    readout = "end-point droplet fluorescence (FAM/HEX); partitioning not simulated",
    zygosity = "co-dominant",
    amplicon = "short (<= 150 bp)",
    loci = 1L,
    verdict = "usable dPCR assay",
    gates = gates,
    # least grave first: widen the amplicon window, then tolerate a stickier probe
    relax = c("gate_dpcr_amplicon", "gate_dpcr_probe_structure"),
    tiers = c(
      tier1 = "probe Tm margin (out-competes the primers)",
      tier2 = "probe self-structure (hairpin / homodimer)",
      tier3 = "amplicon length (near the ideal)"
    ),
    gate_short = c(
      gate_dpcr_probe_comp = "probe composition",
      gate_dpcr_probe_specific = "probe specificity",
      gate_dpcr_probe_structure = "probe self-structure",
      gate_dpcr_probe_tm = "probe Tm ordering",
      gate_dpcr_primer_tm = "primer Tm (dPCR)",
      gate_dpcr_amplicon = "amplicon length (dPCR)"
    ),
    gate_limit = c(
      gate_dpcr_probe_comp = glue::glue("GC {gc[1L] * 100}-{gc[2L] * 100} %, no 5' G, no run >= 4"),
      gate_dpcr_probe_specific = "exactly one perfect probe site in the genome",
      gate_dpcr_probe_structure = glue::glue("probe hairpin/homodimer dG37 >= {KASP_LIMITS$hairpin_dg}"),
      gate_dpcr_probe_tm = glue::glue("probe Tm >= primer Tm + {pmar} C"),
      # the annealing window the designer applies (design_dpcr_marker's tm_range
      # default): a real parameterised value, never a `[lo, hi]` (C9)
      gate_dpcr_primer_tm = glue::glue("both primers Tm within [{tm[1L]}, {tm[2L]}] C"),
      gate_dpcr_amplicon = glue::glue("amplicon within [{amp[1L]}, {amp[2L]}] bp")
    ),
    # the floors the ladder may relax the two soft gates to (margin units)
    gate_floor = c(
      gate_dpcr_amplicon = -DPCR_AMPLICON_RELAX,
      gate_dpcr_probe_structure = -DPCR_PROBE_STRUCT_RELAX
    ),
    # the numeric windows the engine scores against (so gate margins and the
    # declared limits cannot drift)
    amplicon_window = amp,
    primer_tm_window = tm,
    probe_gc_window = gc,
    gate_floor_label = c(
      gate_dpcr_amplicon = glue::glue("widen the amplicon window by {DPCR_AMPLICON_RELAX} bp"),
      gate_dpcr_probe_structure = glue::glue("allow {DPCR_PROBE_STRUCT_RELAX} kcal/mol more stable")
    ),
    gate_hard_label = c(
      gate_dpcr_probe_comp = "a 5' G quenches the reporter and a GC drift kills hybridisation",
      gate_dpcr_probe_specific = "a probe that binds elsewhere is cleaved there and lights the wrong droplet",
      gate_dpcr_probe_tm = "the probe must out-compete the primers",
      gate_dpcr_primer_tm = "both primers must anneal in the extension window"
    ),
    stage_label = "dPCR (droplet digital PCR)"
  )
}
