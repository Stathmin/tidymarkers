# =============================================================================
#  marker_castpcr.R  --  castPCR / competitive allele-specific TaqMan (mutant LOD)
# =============================================================================
#  castPCR detects a LOW-ABUNDANCE MUTANT (ALT) allele against a wild-type
#  background - its whole purpose is the mutant, never the wild type.  One assay
#  is four oligos in one tube:
#    * a MUTANT allele-specific primer (ASP) whose 3' base is the mutant base, so
#      it extends the mutant template and its 3'-terminal mismatch stalls on the
#      wild type - the same 3'-anchor mechanism KASP uses;
#    * a LOCUS-specific hydrolysis (TaqMan/MGB) probe, NOT an allele-specific one:
#      it reports whether the ASP's amplicon was made, it does not read the base
#      itself, so it must sit inside the amplicon without overlapping the ASP;
#    * a WILD-TYPE-specific MGB oligonucleotide blocker (ASB) whose 3' base is the
#      WT base, so it blocks wild-type amplification while mismatching (and so
#      sparing) the mutant.  The blocker is what pushes the LOD down to the
#      sub-percent range;
#    * a common reverse primer.
#  The readout is the mutant assay's Ct (or Delta-Ct against a reference assay),
#  NOT a two-colour ratio: one ASP cannot call REF/ALT zygosity, so castPCR is a
#  mutant-presence (dominant) assay - not co-dominant, as KASP is.
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
#  HONEST LIMIT: the package cannot simulate the competitive blocker kinetics that
#  produce that 0.1-1% LOD, and the MGB modification (which lifts a short oligo's
#  effective Tm) is not modelled either.  So the blocker is *placed and screened*
#  (WT-anchored, composition-clean) and no numeric LOD is claimed; the marker is a
#  mutant-anchored AS-qPCR design, not a predicted limit of detection.
#
#  Design thresholds the method papers do not state numerically are taken from the
#  package's own KASP limits where they coincide (AS Tm, GC) and from the standard
#  hydrolysis-probe rule of thumb - probe Tm above the primers, short amplicon -
#  which is stated as an ASSUMPTION and is a parameter the caller can override.
# =============================================================================

# Hydrolysis-probe rule of thumb: a TaqMan/MGB probe must melt above the primers
# so the probe carries the real-time readout.  This is the SAME physical rule the
# dPCR scheme hard-gates (`gate_dpcr_probe_tm`, "probe Tm >= primer Tm + 8 C") and
# is hard-gated here too - the chemistry does not change with the readout.  The
# numeric margin is smaller than dPCR's because the castPCR probe is an MGB oligo:
# the MGB modification lifts the effective Tm well above the unmodified value this
# package computes, so the unmodified-oligo floor is conservative.
CASTPCR_PROBE_MARGIN <- 8
CASTPCR_PROBE_TM_MAX <- 80

# castPCR keeps the product short so the real-time protocol stays fast.  A short
# product also matters at the method's LOD: the mutant template is a tiny
# fraction, so the amplicon must amplify at maximum efficiency.  The window is the
# marker-agnostic amplicon limit (`KASP_LIMITS$amplicon`, read lazily below so
# this file does not depend on the load order), so a marker the gate admits still
# passes the acceptance check on the same criterion (C5).
# How far the ladder may widen that window (bp each side) - a bounded floor (C8).
CASTPCR_AMPLICON_FLOOR <- 10L
# The ladder may widen the ASP Tm window by this many degrees each side (C).
CASTPCR_TM_FLOOR <- 2
# Real ASB oligos are short MGB oligos; the MGB is not modelled, so lengths are
# the unmodified-oligo window and the blocker is screened on composition.
CASTPCR_BLOCKER_LENGTHS <- 14:24
# How many ASP candidates one attempt may try (best discrimination first) before
# the ladder widens a window: the top ASP's geometry may not admit a mate inside
# the amplicon window even though a slightly weaker one does.
CASTPCR_ASP_TRIES <- 8L

# The castPCR decision stage.  Gates in application order (grave / cheap first):
# the mutant anchor (which base the ASP reads), its 3' mismatch discrimination,
# the WT blocker, then the two windows (ASP Tm, short amplicon), and the probe Tm
# ordering the readout depends on.  Hard = the assay is not castPCR without it;
# soft = the ladder may widen it, least grave first.
KASP_GATES_CASTPCR <- c(
  "gate_castpcr_as3", "gate_castpcr_mismatch", "gate_castpcr_blocker",
  "gate_castpcr_as_tm", "gate_castpcr_probe_tm", "gate_castpcr_amplicon"
)
KASP_RELAX_CASTPCR <- c(
  "gate_castpcr_amplicon", "gate_castpcr_as_tm"
)
KASP_TIERS_CASTPCR <- c(
  tier1 = "3'-terminal specificity (ASP mutant/WT margin)",
  tier2 = "Tm balance (ASP x common mate)", tier3 = "amplicon length"
)

#' Hydrolysis-probe candidates inside a primer pair (no overlap)
#'
#' The castPCR probe is *locus-specific*, not allele-specific: it only has to
#' report the amplicon, so it is placed in the inner region between the two
#' primers (never overlapping the ASP, which would make the two compete for the
#' same bases), on the plus strand the polymerase reads, avoiding a 5'-G (which
#' quenches the reporter even after cleavage) and self-complementarity.  Reuses
#' `primer_candidates()` for every composition check and orders so the probe Tm
#' closest to the middle of the accepted window comes first (robust, not marginal).
#' @param seq the template (plus strand).
#' @param lo,hi 1-based inclusive inner region (between the primers).
#' @param lengths probe lengths to try.
#' @param tm_range,gc_range accepted windows.
#' @param conditions reaction conditions for [primer_tm()] (`NULL` = the
#'   historical 50 mM Na / 37 C buffer).
#' @return a tibble (`start`, `len`, `end`, `primer`, `tm`, `gc`), possibly zero
#'   rows.
#' @keywords internal
#' @noRd
castpcr_probe_candidates <- function(seq, lo, hi, lengths, tm_range, gc_range,
                                     conditions = NULL) {
  empty <- tibble::tibble(
    start = integer(), len = integer(), end = integer(),
    primer = character(), tm = numeric(), gc = numeric()
  )
  if (hi < lo) {
    return(empty)
  }
  cand <- primer_candidates(
    seq, seq.int(lo, hi), lengths, "+", tm_range, gc_range,
    conditions = conditions
  ) |>
    dplyr::filter(.data$start >= lo, .data$end <= hi) |>
    dplyr::filter(substr(.data$primer, 1L, 1L) != "G") |>
    dplyr::filter(!has_palindrome(.data$primer))
  if (!nrow(cand)) {
    return(empty)
  }
  cand |>
    dplyr::mutate(tm_gap = abs(.data$tm - mean(tm_range))) |>
    dplyr::arrange(.data$tm_gap, .data$len)
}

#' The wild-type-specific MGB blocker for a castPCR assay
#'
#' The ASB is the wild-type counterpart of the mutant ASP: same strand and
#' geometry, its 3' base on the variant, reading the WT base.  That is exactly
#' what makes it wild-type specific - its own 3'-terminal mismatch spares the
#' mutant (Barbano 2015).  The MGB modification is not modelled, so the oligo is
#' screened on composition alone (no homopolymer run, no SSR, not
#' self-complementary, GC in the accepted window) and the longest clean option
#' (the most stable WT binding) is returned; `NULL` when none passes.
#' @param seq the template (plus strand).
#' @param si 1-based position of the variant on `seq`.
#' @param strand the ASP's strand (`"+"`/`"-"`) - the blocker competes with it.
#' @param ref,alt the WT and mutant bases.
#' @param lengths blocker lengths to try.
#' @param gc_range accepted GC window.
#' @return a one-row tibble (`len`, `primer`, `gc`), or `NULL`.
#' @keywords internal
#' @noRd
castpcr_blocker <- function(seq, si, strand, ref, alt, lengths, gc_range) {
  # the blocker's 3' base must read WT, so it mismatches the mutant at its own 3'
  # end (what makes it block the wild type and spare the mutant)
  mutant_3p <- substr(make_as_primer(seq, si, 1L, alt, strand = strand), 1L, 1L)
  cand <- purrr::map(lengths, function(L) {
    tibble::tibble(len = L, primer = make_as_primer(seq, si, L, ref, strand = strand))
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(gc = gc_frac(.data$primer)) |>
    dplyr::filter(substr(.data$primer, nchar(.data$primer), nchar(.data$primer)) != mutant_3p) |>
    dplyr::filter(
      !has_run(.data$primer), !has_repeat(.data$primer), !has_palindrome(.data$primer),
      dplyr::between(.data$gc, gc_range[1L], gc_range[2L])
    )
  if (!nrow(cand)) {
    return(NULL)
  }
  cand[which.max(cand$len), ]
}

#' The most Tm-balanced AS x common pair that still leaves room for the probe
#'
#' The castPCR sibling of `pick_primer_pair()`: the same Tm-balance ranking, but
#' the inner gap between the two primers must hold at least `min_gap` bases (the
#' shortest accepted probe), because the locus probe has to fit between them
#' without overlapping either primer.
#' @keywords internal
#' @noRd
pick_castpcr_pair <- function(fwd, rev, min_product, max_product, min_gap) {
  if (!nrow(fwd) || !nrow(rev)) {
    return(NULL)
  }
  pairs <- dplyr::cross_join(fwd, rev, suffix = c("_f", "_r")) |>
    dplyr::mutate(
      product = .data$end_r - .data$start_f + 1L,
      gap = .data$start_r - .data$end_f - 1L,
      dtm = abs(.data$tm_f - .data$tm_r)
    ) |>
    dplyr::filter(dplyr::between(.data$product, min_product, max_product), .data$gap >= min_gap)
  if (!nrow(pairs)) {
    return(NULL)
  }
  pairs[order(pairs$dtm, pairs$product), ][1L, ]
}

#' Design a castPCR / allele-specific qPCR marker at a variant
#'
#' The real-time analogue of [assemble_kasp_assay()] that targets the mutant
#' allele - a low-abundance ALT detected against a wild-type background.  One
#' mutant allele-specific primer has its 3' base on the mutant base, a
#' locus-specific hydrolysis (TaqMan/MGB) probe reports the amplicon, a
#' wild-type-specific MGB blocker (its 3' base on the WT base) suppresses
#' wild-type amplification without touching the mutant, and a common mate closes
#' a short amplicon.  The probe is *not* allele-specific: the allele specificity
#' is carried entirely by the ASP's 3' mismatch plus the blocker's, so the probe
#' only has to sit inside the amplicon and melt above the ASP.
#'
#' The returned marker records its gate outcomes in `selection`: the ASP Tm
#' window and the amplicon window are the two soft gates, any of which the ladder
#' may widen (least grave first) and report as relaxed.  The mutant anchor, the
#' ASP's mutual discrimination, the WT blocker and the probe Tm ordering are hard:
#' without them the assay is not the method (see the file header).
#'
#' HONEST LIMIT: the blocker's suppression efficiency - and so the 0.1-1% LOD the
#' method papers report - is not simulated (see the file header).
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param as_lengths mutant allele-specific primer lengths to try.
#' @param mate_lengths common-primer lengths to try.
#' @param probe_lengths hydrolysis-probe lengths to try.
#' @param blocker_lengths wild-type blocker lengths to try.
#' @param probe_margin minimum probe Tm excess (C) over the ASP (probe must
#'   melt above the primer); defaults to `CASTPCR_PROBE_MARGIN`.
#' @param probe_tm_max upper bound (C) on the probe Tm; defaults to
#'   `CASTPCR_PROBE_TM_MAX`.
#' @param min_product,max_product the accepted amplicon window.
#' @param tm_range accepted AS primer Tm window.
#' @param gc_range accepted GC window for every oligo.
#' @return a `marker_castpcr` list (`scheme`, `snp`, `primers`, `blocker`, `probe`,
#'   `amplicon`, `selection` - the gate outcomes, `compatibility` - the joint
#'   dimer/hairpin screen over the ASP, its mate, the probe and the blocker), or
#'   `NULL` (with a message) when no ASP, mate, probe or blocker passes the
#'   filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_castpcr_marker(genome, marker_example("SNP")[1, ])
#' }
#' @family KASP castPCR
#' @export
design_castpcr_marker <- function(genome, snp, flank = 150L,
                                  as_lengths = 18:27, mate_lengths = 18:25,
                                  probe_lengths = 18:32,
                                  blocker_lengths = CASTPCR_BLOCKER_LENGTHS,
                                  probe_margin = CASTPCR_PROBE_MARGIN,
                                  probe_tm_max = CASTPCR_PROBE_TM_MAX,
                                  min_product = KASP_LIMITS$amplicon[1L],
                                  max_product = KASP_LIMITS$amplicon[2L],
                                  tm_range = KASP_LIMITS$as_tm, gc_range = KASP_LIMITS$gc) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    cli::cli_inform(gettext("castPCR is implemented for single-base substitutions only."))
    return(NULL)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  conditions <- scheme_conditions("castpcr")
  si <- tmpl$snp_i
  n <- nchar(as.character(tmpl$seq))
  # fetch the ASP candidates over the *widened* Tm window in one DECIPHER batch,
  # so the ladder's Tm rung below costs nothing extra
  as_cand <- design_as_primers(tmpl$seq, si, ref, alt,
    strands = c("+", "-"), lengths = as_lengths,
    gc_range = gc_range, tm_range = tm_range + c(-CASTPCR_TM_FLOOR, CASTPCR_TM_FLOOR),
    conditions = conditions
  )
  if (!nrow(as_cand)) {
    cli::cli_inform(gettext("No allele-specific primer passed the GC/Tm filters for castPCR."))
    return(NULL)
  }
  # gate_castpcr_mismatch reads the ASP's OWN margin - the mutant may be a
  # vanishing fraction of the template, so its on-target efficiency is what a
  # positive call rests on.  Tier 1 ranks on the same margin.
  as_cand$asp_tm <- primer_tm(as_cand$alt_primer, conditions = conditions)
  as_cand$as_discrimination <- as_cand$alt_on - as_cand$alt_off
  as_cand <- as_cand[order(-as_cand$as_discrimination), ]

  # one full assay for one ASP candidate under a product window: the ASP anchor,
  # the mate pair (tier 2/3) and the locus probe (hard Tm ordering)
  build <- function(best, prod_win) {
    # gate_castpcr_as3: the ASP's 3' terminus must sit on the variant base
    on_variant <- if (identical(best$strand, "+")) {
      best$tmpl_end == si
    } else {
      best$tmpl_start == si
    }
    if (!on_variant) {
      return(NULL)
    }
    # a one-row candidate set in the schema [primer_candidates()] produces, so
    # the pair picker can treat the ASP like any other primer
    as_one <- tibble::tibble(
      start = best$tmpl_start, end = best$tmpl_end,
      primer = best$alt_primer, tm = best$asp_tm
    )
    if (identical(best$strand, "+")) {
      mate <- primer_candidates(
        tmpl$seq, seq.int(best$tmpl_end + 1L, n), mate_lengths, "-",
        KASP_LIMITS$com_tm, gc_range,
        conditions = conditions
      )
      b <- pick_castpcr_pair(as_one, mate, prod_win[1L], prod_win[2L], min_gap)
    } else {
      if (best$tmpl_start <= 1L) {
        return(NULL)
      }
      mate <- primer_candidates(
        tmpl$seq, seq_len(best$tmpl_start - 1L), mate_lengths, "+",
        KASP_LIMITS$com_tm, gc_range,
        conditions = conditions
      ) |>
        dplyr::filter(.data$end < best$tmpl_start)
      b <- pick_castpcr_pair(mate, as_one, prod_win[1L], prod_win[2L], min_gap)
    }
    if (is.null(b)) {
      return(NULL)
    }
    # the locus-specific probe: inside the amplicon, between the primers, never
    # overlapping the ASP, and melting above it (gate_castpcr_probe_tm, hard -
    # the same rule dPCR hard-gates)
    inner_lo <- b$start_f + nchar(b$primer_f)
    inner_hi <- b$end_r - nchar(b$primer_r)
    probe <- castpcr_probe_candidates(
      tmpl$seq, inner_lo, inner_hi, probe_lengths,
      c(best$asp_tm + probe_margin, probe_tm_max), gc_range,
      conditions = conditions
    )
    if (!nrow(probe)) {
      return(NULL)
    }
    # gate_castpcr_blocker: the WT-specific blocker must be placeable
    blocker <- castpcr_blocker(tmpl$seq, si, best$strand, ref, alt, blocker_lengths, gc_range)
    if (is.null(blocker)) {
      return(NULL)
    }
    list(best = best, pair = b, probe = probe[1L, ], blocker = blocker)
  }

  # one attempt under a given (Tm window, product window): the best ASP first
  # (tier 1), then the next few, so the ladder widens a window only when no
  # candidate works at all
  attempt <- function(tm_win, prod_win) {
    asp <- as_cand[
      as_cand$as_discrimination > 0 &
        dplyr::between(as_cand$asp_tm, tm_win[1L], tm_win[2L]), ,
      drop = FALSE
    ]
    if (!nrow(asp)) {
      return(NULL)
    }
    for (i in seq_len(min(nrow(asp), CASTPCR_ASP_TRIES))) {
      one <- build(asp[i, ], prod_win)
      if (!is.null(one)) {
        return(one)
      }
    }
    NULL
  }

  # the ladder (least grave first): the strict design, then the amplicon window,
  # then the ASP Tm window, then both
  min_gap <- min(probe_lengths)
  prod_open <- c(min_product - CASTPCR_AMPLICON_FLOOR, max_product + CASTPCR_AMPLICON_FLOOR)
  tm_open <- tm_range + c(-CASTPCR_TM_FLOOR, CASTPCR_TM_FLOOR)
  design <- attempt(tm_range, c(min_product, max_product))
  relaxed <- character(0)
  if (is.null(design)) {
    design <- attempt(tm_range, prod_open)
    if (!is.null(design)) relaxed <- "gate_castpcr_amplicon"
  }
  if (is.null(design)) {
    design <- attempt(tm_open, c(min_product, max_product))
    if (!is.null(design)) relaxed <- "gate_castpcr_as_tm"
  }
  if (is.null(design)) {
    design <- attempt(tm_open, prod_open)
    if (!is.null(design)) relaxed <- c("gate_castpcr_amplicon", "gate_castpcr_as_tm")
  }
  if (is.null(design)) {
    cli::cli_inform(gettext(
      "No mutant allele-specific primer, common mate, probe and WT blocker passed the castPCR filters."
    ))
    return(NULL)
  }

  best <- design$best
  pair <- pair_to_list(design$pair)
  off <- snp$pos - si
  out <- list(
    scheme = "castpcr",
    snp = snp,
    primers = c(pair, list(
      as_role = if (identical(best$strand, "+")) "forward" else "reverse",
      as_primer = best$alt_primer,
      as_tm = best$asp_tm,
      as_discrimination = best$as_discrimination
    )),
    blocker = list(
      sequence = design$blocker$primer, allele = "WT",
      len = design$blocker$len, gc = design$blocker$gc,
      modification = "MGB, non-extendable (not modelled)"
    ),
    probe = list(
      sequence = design$probe$primer, start = design$probe$start + off,
      end = design$probe$end + off, tm = design$probe$tm, gc = design$probe$gc,
      kind = "locus-specific (reports the ASP amplicon)"
    ),
    amplicon = tibble::tibble(start = pair$start + off, end = pair$end + off),
    # the gate outcomes the graph and a report read
    selection = list(relaxed = relaxed, violated = character(0))
  )
  # the ASP, its mate, the locus probe and the WT blocker share the tube: a joint
  # screen of all four, not each oligo alone
  out$compatibility <- tube_compatibility(c(
    as_primer = out$primers$as_primer, common = out$primers$forward,
    mate = out$primers$reverse, probe = out$probe$sequence, blocker = out$blocker$sequence
  ), conditions = conditions)
  structure(out, class = "marker_castpcr")
}

#' @export
print.marker_castpcr <- function(x, ...) {
  cat(
    "<marker_castpcr>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, " (detects ALT ", x$snp$alt, ")\n",
    sep = ""
  )
  cat("  ASP (", x$primers$as_role, "): ", x$primers$as_primer,
    " (Tm ", round(x$primers$as_tm, 1L), " C)\n",
    sep = ""
  )
  cat("  WT blocker:", x$blocker$sequence, "\n", sep = "")
  cat("  probe:", x$probe$sequence, " (Tm ", round(x$probe$tm, 1L), " C)\n", sep = "")
  cat("  amplicon:", x$primers$product, "bp\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  if (length(x$selection$relaxed)) {
    cat("  relaxed:", paste(x$selection$relaxed, collapse = ", "), "\n")
  }
  invisible(x)
}

#' The castPCR scheme definition for the shared registry
#'
#' A self-describing definition of the castPCR chemistry in the shape the marker
#' registry expects (scheme/family/label/readout, the gate/tier vocabulary and the
#' per-gate human-readable labels).  Kept beside the designer so the chemistry and
#' its thresholds travel together; it is discovered automatically (one file per
#' scheme, see `R/marker_scheme_registry.R`).
#' @return a named list describing the scheme.
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_castpcr
#' @export
#' @noRd
marker_map_spec.marker_castpcr <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  start <- sp[["start"]]
  end <- sp[["end"]]
  probe <- if (!is.null(marker$probe)) {
    list(spec_row(
      "probe", marker$probe$start[1L], marker$probe$end[1L], "+", "probe",
      label = "probe"
    ))
  } else {
    list()
  }
  bl <- marker$blocker
  blocker <- if (!is.null(bl) && nzchar(bl$sequence)) {
    len <- nchar(bl$sequence)
    role <- marker$primers$as_role %||% "forward"
    bs <- if (identical(role, "forward")) marker$snp$pos[1L] - len + 1L else marker$snp$pos[1L]
    list(spec_row(
      "blocker", bs, bs + len - 1L,
      if (identical(role, "forward")) "+" else "-", "probe",
      label = "blocker"
    ))
  } else {
    list()
  }
  features <- c(
    map_pair_spec(marker$primers$forward, marker$primers$reverse, start, end),
    probe, blocker
  )
  sequences <- list(
    forward = marker$primers$forward, reverse = marker$primers$reverse,
    probe = if (!is.null(marker$probe)) marker$probe$sequence else NULL,
    blocker = if (!is.null(bl)) bl$sequence else NULL
  )
  # the mutant-specific ASP is a primer in its own right; draw it unless the
  # designer reused the common reverse primer
  asp <- marker$primers$as_primer
  if (!is.null(asp) && nzchar(asp) && !asp %in% c(marker$primers$forward, marker$primers$reverse)) {
    role <- marker$primers$as_role %||% "forward"
    astart <- if (identical(role, "forward")) start else end - nchar(asp) + 1L
    aend <- if (identical(role, "forward")) start + nchar(asp) - 1L else end
    features <- c(features, list(spec_row(
      "AS", astart, aend, if (identical(role, "forward")) "+" else "-", "primer"
    )))
    sequences$AS <- asp
  }
  map_spec(marker, features = features, sequences = sequences)
}

scheme_def_castpcr <- function() {
  as_tm <- KASP_LIMITS$as_tm
  margin <- CASTPCR_PROBE_MARGIN
  amplicon <- KASP_LIMITS$amplicon
  list(
    scheme = "castpcr",
    family = "castpcr",
    label = "castPCR / competitive allele-specific TaqMan (mutant detection)",
    discrimination = "3'-terminal mismatch of the mutant-anchored AS primer, backstopped by a WT-specific MGB blocker",
    readout = "locus-specific hydrolysis (TaqMan/MGB) probe, mutant assay Ct (Delta-Ct)",
    zygosity = "dominant (mutant detection)",
    amplicon = glue::glue("short ({amplicon[1L]}-{amplicon[2L]} bp)"),
    loci = 1L,
    verdict = "usable castPCR assay",
    gates = KASP_GATES_CASTPCR,
    relax = KASP_RELAX_CASTPCR,
    tiers = KASP_TIERS_CASTPCR,
    gate_short = c(
      gate_castpcr_as3 = "ASP 3' base on variant",
      gate_castpcr_mismatch = "3' mismatch (ASP)",
      gate_castpcr_blocker = "WT blocker",
      gate_castpcr_as_tm = "ASP Tm",
      gate_castpcr_probe_tm = "probe Tm > ASP Tm",
      gate_castpcr_amplicon = "amplicon length (castPCR)"
    ),
    gate_limit = c(
      gate_castpcr_as3 = "the ASP's 3' terminus sits on the variant base",
      gate_castpcr_mismatch = "ASP 3' base = the mutant base, predicted to prefer the mutant",
      gate_castpcr_blocker = "a WT-specific blocker (3' base = the WT base) is placed",
      gate_castpcr_as_tm = glue::glue("ASP Tm within [{as_tm[1L]}, {as_tm[2L]}] C"),
      gate_castpcr_probe_tm = glue::glue("probe Tm >= ASP Tm + {margin} C"),
      gate_castpcr_amplicon = glue::glue("amplicon within [{amplicon[1L]}, {amplicon[2L]}] bp")
    ),
    gate_floor_label = c(
      gate_castpcr_as_tm = glue::glue("widen the ASP Tm window by {CASTPCR_TM_FLOOR} C"),
      gate_castpcr_amplicon = glue::glue("widen the amplicon window by {CASTPCR_AMPLICON_FLOOR} bp")
    ),
    gate_hard_label = c(
      gate_castpcr_as3 = "without the mutant anchor the assay is not allele-specific",
      gate_castpcr_mismatch = "an ASP that does not prefer the mutant has no discrimination",
      gate_castpcr_blocker = "the WT blocker is castPCR's specificity - without it the ASP leaks wild type",
      gate_castpcr_probe_tm = "the probe must out-compete the primers for the real-time signal"
    ),
    stage_label = "castPCR / AS-qPCR (mutant detection)"
  )
}
