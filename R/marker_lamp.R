# =============================================================================
#  marker_lamp.R  --  LAMP (loop-mediated isothermal amplification) as a marker scheme
# =============================================================================
#  LAMP is the isothermal, primer-heavy sibling of PCR: a strand-displacing
#  polymerase (Bst) runs at a constant 60-65 C, and 4-6 primers recognise 6-8
#  regions of a short target.  Four are mandatory - the outer F3/B3 pair and the
#  inner FIP/BIP pair (each inner primer is a fusion of a loop arm and a
#  displacement arm: FIP = F1c + F2, BIP = B1c + B2); two loop primers (LF/LB)
#  are optional accelerators.  Forward primers (F3, FIP, LF) are plus-strand
#  oligos that extend rightward; reverse primers (B3, BIP, LB) are minus-strand
#  oligos that extend leftward - the same convention as design_flanking_pair().
#
#  Because there is no thermal cycling, the assay lives or dies on primer
#  placement and thermodynamics.  For ALLELE discrimination the single governing
#  fact is that a mismatched 3'-terminal base blocks the strand-displacing
#  polymerase: the variant must sit exactly at the 3' terminus of an inner
#  primer's displacement arm (FIP's F2 here), or of a loop primer.  A variant
#  merely *inside* a primer's footprint, or left in the single-stranded loop
#  between the inner primers' arms, is read by no primer 3' end and cannot
#  discriminate the alleles - it is an unprimed loop, not a marker.  So FIP's 3'
#  arm is built to *end* at the variant, putting the variant base at FIP's
#  3'-terminal position (its 3' base pairs the variant), and the loop primers
#  LF/LB are placed in the F1c/F2 and B1c/B2 loops.
#
#  Design rules used here (web-researched; see provenance below):
#    * 4 primers (F3/B3/FIP/BIP) plus the optional loop pair; amplicon 120-300 bp
#      - "optimal 130-200 bp, sequences up to 300 bp are also acceptable (Notomi
#      2000)" (Shirshikov et al. 2019, PeerJ 7:e6801, PMC6487805; MorphoCatcher);
#      NEB quotes a 80-250 nt target and 4-6 primers
#      (https://www.neb.com/.../loop-mediated-isothermal-amplification-lamp).
#    * isothermal 60-65 C; inner FIP/BIP bind two distant sites and are longer
#      (45-49 bp; our default arm window 20:24 yields 40-48 nt, and the designer
#      lands on the long end) than the outer F3/B3 (21-24 bp); the primers must
#      not form a stable double strand at that temperature
#      (Soroka, Wasowicz & Rymaszewska 2021, Cells 10(8):1931, PMC8393631;
#      Notomi et al. 2000; Nagamine et al. 2002 for loop primers).
#    * the variant is placed at the 5'/3' ends of the inner primers and the 3'
#      ends of the outer / loop primers - "the most suitable regions ... are the
#      5'- and 3'-ends of inner (FIP and BIP) primers, as well as 3'-ends of
#      outer (F3 and B3) and loop (LF and LB) primers"
#      (Shirshikov et al. 2019, PeerJ 7:e6801).  We realise the 3'-end placement
#      (the one a mismatch actually blocks).
#
#  Allele specificity and controls: one returned set is ALT-specific.  FIP's 3'
#  base reads the ALT allele, so a REF-only template yields no signal; a REF
#  sample and a failed reaction are therefore indistinguishable.  The `paired`
#  block (Z9.1) therefore returns the REF counterpart of the same set (identical
#  positions, FIP 3' base flipped to REF), so the two tubes are a deliberate
#  control pair; `paired$status` says whether it could be built.  A pair is not
#  a discrimination model - no allele gap, Ct or call probability is predicted
#  (data-raw/proposed/lamp_status.decision.md).
#
#  Threshold provenance / assumptions:
#    * outer Tm window c(55, 62) C and inner (binding-arm) Tm window c(55, 65) C
#      are the rule-of-thumb bands above expressed in the package's TmCalculator
#      conditions (PCR_COND).  The published inner-arm target is 60-65 C, but that
#      is a different (higher-salt / different NN) Tm model: under PCR_COND a 20-24
#      nt arm with ~50% GC computes to ~50-57 C, so the literature band maps to
#      [55, 65] here - the same window the package already uses for its other
#      3'-terminal, allele-discriminating primer (KASP's AS primer,
#      `KASP_LIMITS$as_tm`).  Treat it as an explicit assumption, not a fact.
#    * GC 40-60 % and "no run >= 4 / no SSR" are the primer-agnostic KASP limits
#      (KASP_LIMITS), reused unchanged (see marker_primer_pairs.R).
#    * the structural screen is sequence-level (no inverted repeat), reusing
#      has_palindrome(); the whole primer set is then screened jointly for
#      cross-dimers/hairpins by primer_set_compatibility().
# =============================================================================

# The LAMP decision stage, declared here so scheme_def_lamp() is self-contained.
# Gate order = the order the method applies them (grave / cheap first): the
# variant placement decides discrimination and is tested first, then the inner
# (isothermal) Tm, then the outer Tm, then the sequence-level structure, the
# joint cross-dimers and finally the amplicon length.
LAMP_GATES <- c(
  "gate_lamp_variant", "gate_lamp_inner_tm", "gate_lamp_outer_tm",
  "gate_lamp_structure", "gate_lamp_cross_dimer", "gate_lamp_amplicon"
)
# relaxation order, LEAST grave first (C6): the amplicon window, then the outer
# Tm window, then a slightly longer tolerated cross-dimer.
LAMP_RELAX <- c("gate_lamp_amplicon", "gate_lamp_outer_tm", "gate_lamp_cross_dimer")
LAMP_TIERS <- c(
  tier1 = "variant discrimination", tier2 = "inner Tm set",
  tier3 = "primer structure", tier4 = "amplicon length"
)

# The variant must sit at the 3'-terminal base of FIP's displacement arm (F2),
# and the single-stranded loop between the inner primer's arms must leave room for
# a loop primer, so a loop of at least this many nt is required.
LAMP_LOOP_MIN <- 16L
# loop-primer (LF/LB) lengths to try.
LAMP_LOOP_LEN <- 15:24L
# arm sites kept per inner-primer arm before pairing (bounds the cross join; the
# arms are ranked by binding Tm first, so the kept sites are the thermodynamically
# sensible ones, spread over the region).
LAMP_ARM_TOP <- 80L
# relaxation floors for the soft gates (bounded; C8).
LAMP_AMPLICON_RELAX <- 20L # widen the amplicon window by this on each side
LAMP_OUTER_TM_RELAX <- 3 # widen the outer Tm window by this (C)
LAMP_DIMER_SLACK <- 2L # allow this many more bp in a cross-dimer run
# a 5-6 primer LAMP pool has ~15-21 pairs, so some short complementarity is
# routine; a run of 8 bp (not 6) is what actually risks a cross-dimer, matching
# the LAMP cross-dimer gate.
LAMP_CROSS_RUN <- 8L
LAMP_SELFCOMP <- 6L

# Geometry limits for the inner (FIP/BIP) pair.  These are diagnostics layered on
# top of the decision gates: the primer *choice* is unchanged, but a selected set
# whose geometry is physically implausible is reported (and warned about) instead
# of being accepted silently.
#
# F2-B2 distance = plus-strand bases from the 3' terminus of F2 (FIP's
# displacement arm) to the 3' terminus of B2 (BIP's) inclusive.  The B1 site and a
# single-stranded loop sit between the two arms, so the floor is two loop minima;
# the inner span can never exceed the amplicon ceiling (the caller passes it).
LAMP_F2B2_MIN <- 2L * LAMP_LOOP_MIN
# 3'-end stability of an inner primer: the terminal base should be a G/C clamp
# and the 3' end must not fold back on an internal region (duplex_selfcomp()).
LAMP_3P_BASES <- 5L
LAMP_3P_CLAMP <- c("G", "C")
# Deliberately weak dG floor for the terminal 5-mer (terminal_dg(), kcal/mol):
# FIP's 3'-terminal base is the ALT allele, so it cannot always be a strong clamp.
LAMP_3P_DG_MAX <- -0.5

#' The LAMP scheme definition (registry row + decision stage)
#'
#' Self-contained description of the LAMP chemistry in the shape the shared
#' registry expects ([marker_schemes()] fields plus a decision stage: `gates`,
#' `relax`, `tiers` and their human-readable labels).  Declared next to the
#' designer so a scheme can be added without touching the registry first; the
#' coordinator wires it in afterwards.
#' @return a list (see the `scheme_def_*` naming convention used for the
#'   registry's schemes).
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_lamp
#' @export
#' @noRd
marker_map_spec.marker_lamp <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  start <- sp[["start"]]
  end <- sp[["end"]]
  f3 <- nchar(marker$primers$F3)
  b3 <- nchar(marker$primers$B3)
  loop <- marker$loop
  arms <- marker$lamp_arms
  if (!is.null(arms)) {
    fip <- marker$primers$FIP
    bip <- marker$primers$BIP
    fa <- arms$FIP
    ba <- arms$BIP
    n5f <- fa$arm5[2L] - fa$arm5[1L] + 1L
    n5b <- ba$arm5[2L] - ba$arm5[1L] + 1L
    return(map_spec(
      marker,
      features = list(
        spec_row("F3", start, start + f3 - 1L, "+", "primer"),
        spec_row("B3", end - b3 + 1L, end, "-", "primer"),
        spec_row("FIP 5'", fa$arm5[1L], fa$arm5[2L], "-", "primer"),
        spec_row("FIP 3'", fa$arm3[1L], fa$arm3[2L], "+", "primer"),
        spec_row("BIP 5'", ba$arm5[1L], ba$arm5[2L], "-", "primer"),
        spec_row("BIP 3'", ba$arm3[1L], ba$arm3[2L], "-", "primer")
      ),
      sequences = list(
        F3 = marker$primers$F3, B3 = marker$primers$B3,
        "FIP 5'" = substr(fip, 1L, n5f),
        "FIP 3'" = substr(fip, n5f + 1L, nchar(fip)),
        "BIP 5'" = substr(bip, 1L, n5b),
        "BIP 3'" = substr(bip, n5b + 1L, nchar(bip))
      )
    ))
  }
  map_spec(
    marker,
    features = list(
      spec_row("F3", start, start + f3 - 1L, "+", "primer"),
      spec_row("B3", end - b3 + 1L, end, "-", "primer"),
      spec_row("FIP", start, loop$start[1L] - 1L, "+", "primer"),
      spec_row("BIP", loop$end[1L] + 1L, end, "-", "primer")
    ),
    sequences = list(
      F3 = marker$primers$F3, B3 = marker$primers$B3,
      FIP = marker$primers$FIP, BIP = marker$primers$BIP
    )
  )
}

scheme_def_lamp <- function() {
  list(
    scheme = "lamp",
    family = "lamp",
    label = "LAMP (loop-mediated isothermal amplification)",
    discrimination = "variant at the 3' terminus of FIP (F2) or a loop primer",
    readout = "turbidity / fluorescence / colorimetric / gel (isothermal)",
    zygosity = "dominant (present/absent), co-dominant with probes",
    amplicon = "120-300 bp",
    loci = 1L,
    verdict = "usable LAMP marker",
    gates = LAMP_GATES,
    relax = LAMP_RELAX,
    tiers = LAMP_TIERS,
    gate_short = c(
      gate_lamp_variant = "variant at a primer 3' end",
      gate_lamp_inner_tm = "FIP/BIP Tm",
      gate_lamp_outer_tm = "F3/B3 Tm",
      gate_lamp_structure = "primer structure",
      gate_lamp_cross_dimer = "cross-dimers",
      gate_lamp_amplicon = "amplicon length (LAMP)"
    ),
    gate_limit = c(
      gate_lamp_variant = "variant at the 3' terminus of FIP (or a loop primer)",
      gate_lamp_inner_tm = "FIP/BIP binding arm Tm within [55, 65] C",
      gate_lamp_outer_tm = "F3/B3 Tm within [55, 62] C",
      gate_lamp_structure = "no inverted repeat in any primer",
      gate_lamp_cross_dimer = "no cross-dimer run >= 8 bp in the pooled set",
      gate_lamp_amplicon = "amplicon within [120, 300] bp"
    ),
    gate_floor_label = c(
      gate_lamp_amplicon = "widen the amplicon window by 20 bp",
      gate_lamp_outer_tm = "widen the outer Tm window by 3 C",
      gate_lamp_cross_dimer = "allow a 2 bp longer cross-dimer run"
    ),
    gate_hard_label = c(
      gate_lamp_variant = "a variant away from a primer 3' end cannot be discriminated",
      gate_lamp_inner_tm = "the inner primers must anneal isothermally",
      gate_lamp_structure = "primer hairpins / inverted repeats defeat LAMP"
    ),
    stage_label = "LAMP primer set"
  )
}

#' Reverse-complement a character vector of oligos
#' @keywords internal
#' @noRd
lamp_revcomp <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

#' A guarded `primer_candidates()` for the narrow windows LAMP probes
#'
#' The shared screener fails in `primer_tm()` when every candidate is removed by
#' its homopolymer / SSR filter; LAMP arm windows can be that narrow, so an empty
#' grid is returned instead.  `strand` selects whether the returned `primer` is
#' the plus strand (`"+"`, a forward primer) or its reverse complement (`"-"`, a
#' reverse primer that extends leftward) - the outer B3 must be the latter.
#' @param tm_range,gc_range accepted windows.
#' @param strand `"+"` or `"-"` (see `primer_candidates()`).
#' @param conditions reaction conditions for [primer_tm()] (see
#'   `primer_candidates()`; `NULL` keeps the historical buffer).
#' @keywords internal
#' @noRd
lamp_candidates <- function(seq, starts, lengths, tm_range, gc_range, strand = "+",
                            conditions = NULL) {
  empty <- tibble::tibble(
    start = integer(0), len = integer(0), end = integer(0),
    primer = character(0), tm = numeric(0), gc = numeric(0)
  )
  out <- tryCatch(
    primer_candidates(seq, starts, lengths, strand, tm_range, gc_range,
      conditions = conditions
    ),
    error = function(e) NULL
  )
  if (is.null(out)) empty else out
}

#' The plus-strand binding sites of one inner-primer arm
#'
#' A set of GC/Tm-screened plus-strand sites that is also free of inverted
#' repeats *in the orientation it will be used* (the site itself for a forward
#' arm, its reverse complement for a reverse arm).  Ranked by closeness to
#' `target_tm` and capped at `top` so the arm pairing is bounded.
#' @keywords internal
#' @noRd
lamp_arm_sites <- function(seq, starts, lengths, tm_range, gc_range, reverse, top) {
  d <- lamp_candidates(seq, starts, lengths, tm_range, gc_range)
  if (!nrow(d)) {
    return(d)
  }
  arm <- if (reverse) lamp_revcomp(d$primer) else d$primer
  d <- d[!has_palindrome(arm), , drop = FALSE]
  d |>
    dplyr::arrange(abs(.data$tm - mean(tm_range))) |>
    dplyr::slice_head(n = top)
}

#' Fuse two arm sites into an inner primer (FIP or BIP)
#'
#' `five` is the 5' (loop) arm's plus-strand site, `three` the 3' (displacement)
#' arm's.  The 5' loop arm is always the reverse complement of its plus-strand
#' site, so FIP is `revcomp(F1 site) + F2 site` and BIP is
#' `revcomp(B1 site) + revcomp(B2 site)` (BIP is a minus-strand oligo, so both
#' arms are complemented).  The arms must be non-overlapping with at least
#' `min_gap` single-stranded loop between them, the fusion must be free of
#' inverted repeats, and the best is the most compact pair (shortest loop, so the
#' outer F3/B3 pair still has room to flank it), then the one whose binding Tm is
#' closest to `target_tm`.
#' @keywords internal
#' @noRd
lamp_fuse <- function(five, three, target_tm, min_gap, reverse = FALSE) {
  if (!nrow(five) || !nrow(three)) {
    return(NULL)
  }
  f <- dplyr::rename(five, s1 = "start", e1 = "end", p1 = "primer", t1 = "tm")
  t <- dplyr::rename(three, s2 = "start", e2 = "end", p2 = "primer", t2 = "tm")
  pairs <- dplyr::cross_join(f, t) |>
    dplyr::filter(
      .data$e1 < .data$s2, .data$s2 - .data$e1 - 1L >= min_gap,
      # inner-primer 5' (F1c/B1c) arm is the hot arm; the 3' displacement
      # arm (F2/B2) must be cooler so the loop opens after the 5' arm binds
      .data$t1 > .data$t2
    ) |>
    dplyr::mutate(
      a1 = lamp_revcomp(.data$p1),
      a2 = if (reverse) lamp_revcomp(.data$p2) else .data$p2,
      tm_bind = pmax(.data$t1, .data$t2)
    ) |>
    dplyr::mutate(primer = paste0(.data$a1, .data$a2)) |>
    dplyr::filter(!has_palindrome(.data$primer))
  if (!nrow(pairs)) {
    return(NULL)
  }
  best <- pairs |>
    dplyr::arrange(.data$s2 - .data$e1, abs(.data$tm_bind - target_tm)) |>
    dplyr::slice_head(n = 1L)
  list(
    primer = best$primer, tm = best$tm_bind,
    arm5 = c(best$s1, best$e1), arm3 = c(best$s2, best$e2)
  )
}

#' Build the FIP inner primer with the variant at its 3' terminus
#'
#' FIP's 3' arm (F2) is a plus-strand site *ending* at the variant, so FIP's
#' 3'-terminal base is the variant base itself (it pairs the REF template base) -
#' the allele-specific extension block.  The 5' arm (F1c) is a plus-strand site
#' upstream, leaving the F loop between them.
#' @param seq the template (plus strand).
#' @param q 1-based variant position on `seq`.
#' @param arm_lengths inner arm lengths to try.
#' @param tm_range,gc_range accepted arm windows.
#' @param target_tm binding-arm Tm to aim for.
#' @param min_loop the minimum single-stranded loop between the arms.
#' @return a list (`primer`, `tm`, `arm5`, `arm3`), or `NULL`.
#' @keywords internal
#' @noRd
lamp_fip <- function(seq, q, arm_lengths, tm_range, gc_range, target_tm, min_loop,
                     alt_base = NULL) {
  if (q <= 1L) {
    return(NULL)
  }
  # The FIP 3' arm must read the ALT allele.  The template is REF, so replace
  # the anchor base before enumerating the 3' arm; the 5' arm is unaffected.
  if (!is.null(alt_base) && nchar(alt_base) == 1L) {
    seq <- as.character(seq)
    substr(seq, q, q) <- alt_base
  }
  # 3' arm: end exactly at the variant (start = q - len + 1)
  three <- lamp_arm_sites(
    seq, q - arm_lengths + 1L, arm_lengths, tm_range, gc_range,
    reverse = FALSE, top = LAMP_ARM_TOP
  )
  three <- dplyr::filter(three, .data$end == q)
  five <- lamp_arm_sites(
    seq, seq_len(max(1L, q - 1L)), arm_lengths, tm_range, gc_range,
    reverse = FALSE, top = LAMP_ARM_TOP
  )
  lamp_fuse(five, three, target_tm, min_loop, reverse = FALSE)
}

#' Build the BIP inner primer downstream of FIP's 3' arm
#'
#' BIP is the reverse inner primer: its 5' arm (B1c) binds just downstream of FIP
#' and its 3' arm (B2) binds farther downstream, leaving the B loop; it is a
#' minus-strand oligo (extends leftward).
#' @param after plus-strand coordinate the BIP arms must start after.
#' @param ... see `lamp_fip()`.
#' @return a list (see `lamp_fip()`), or `NULL`.
#' @keywords internal
#' @noRd
lamp_bip <- function(seq, after, arm_lengths, tm_range, gc_range, target_tm, min_loop) {
  n <- nchar(as.character(seq))
  if (after + 1L >= n) {
    return(NULL)
  }
  sites <- lamp_arm_sites(
    seq, seq.int(after + 1L, n), arm_lengths, tm_range, gc_range,
    reverse = TRUE, top = LAMP_ARM_TOP
  )
  lamp_fuse(sites, sites, target_tm, min_loop, reverse = TRUE)
}

#' Geometry / physics checks for one FIP/BIP pair
#'
#' A LAMP set can pass every thermodynamic gate and still be geometrically
#' degenerate.  This reports the F2-B2 inner distance, the single-stranded loops
#' (the F1-F2 and B1-B2 intra-primer loops plus the F2-B1 loop between the inner
#' primers), and the 3' end of each inner primer.  These are diagnostics, not new
#' hard selection gates: FIP's 3' base is the ALT allele and is therefore not free
#' to be chosen as an ideal G/C clamp.  A failed check is meant to be surfaced by
#' the caller (see `design_lamp_marker()`), never silently accepted.
#' @param fip,bip inner-primer records (`primer`, `arm5`, `arm3`; plus-strand
#'   template coordinates).
#' @param max_product amplicon ceiling; the F2-B2 span cannot exceed it.
#' @return a named list with the observed values and their `*_ok` flags (see the
#'   body), plus `ok`, `failed` and `messages` summarising the failed checks.
#' @keywords internal
#' @noRd
# =============================================================================
#  Geometry / physics checks for the inner (FIP/BIP) pair
# =============================================================================
#  The F2-B2 distance, the single-stranded loops and the inner primers' 3' ends
#  are reported here; `design_lamp_marker()` warns when any of them fails.
lamp_geometry <- function(fip, bip, max_product = 300L) {
  # F2 and B2 are the inner primers' displacement (3') arms; their 3' termini
  # bound the amplifiable core.
  f2_b2 <- bip$arm3[2L] - fip$arm3[2L] + 1L
  loops <- c(
    fip_f1f2 = fip$arm3[1L] - fip$arm5[2L] - 1L,
    f2_b1 = bip$arm5[1L] - fip$arm3[2L] - 1L,
    bip_b1b2 = bip$arm3[1L] - bip$arm5[2L] - 1L
  )
  last_base <- function(p) substr(p, nchar(p), nchar(p))
  clamps <- c(FIP = last_base(fip$primer), BIP = last_base(bip$primer))
  folds <- c(
    FIP = duplex_selfcomp(fip$primer)$fold,
    BIP = duplex_selfcomp(bip$primer)$fold
  )
  # terminal_dg() is the package's 3'-stability (GC-clamp) diagnostic.  Report it
  # when the external engine is available; absence is unknown, not a failure.
  dg <- if (has_oligoarrayaux()) {
    rlang::set_names(
      terminal_dg(c(fip$primer, bip$primer), n = LAMP_3P_BASES),
      c("FIP", "BIP")
    )
  } else {
    c(FIP = NA_real_, BIP = NA_real_)
  }
  stable_3p <- function(nm) {
    !folds[[nm]] && clamps[[nm]] %in% LAMP_3P_CLAMP &&
      (is.na(dg[[nm]]) || dg[[nm]] <= LAMP_3P_DG_MAX)
  }
  short_loops <- names(loops)[loops < LAMP_LOOP_MIN]
  checks <- c(
    f2_b2_distance = f2_b2 >= LAMP_F2B2_MIN && f2_b2 <= max_product,
    loop_length = length(short_loops) == 0L,
    fip_3p_stability = stable_3p("FIP"),
    bip_3p_stability = stable_3p("BIP")
  )
  unstable_msg <- function(nm) {
    bits <- c(
      if (!clamps[[nm]] %in% LAMP_3P_CLAMP) sprintf("terminal base %s is not G/C", clamps[[nm]]),
      if (folds[[nm]]) "3' end folds back on an internal region",
      if (!is.na(dg[[nm]]) && dg[[nm]] > LAMP_3P_DG_MAX) {
        sprintf("terminal dG %.2f > %.2f kcal/mol", dg[[nm]], LAMP_3P_DG_MAX)
      }
    )
    sprintf("%s 3' end is unstable: %s", nm, paste(bits, collapse = "; "))
  }
  messages <- c(
    f2_b2_distance = sprintf(
      "F2-B2 distance %d bp is outside [%d, %d]", f2_b2, LAMP_F2B2_MIN, max_product
    ),
    loop_length = sprintf(
      "single-stranded loop(s) %s shorter than %d nt",
      paste(short_loops, collapse = "/"), LAMP_LOOP_MIN
    ),
    fip_3p_stability = unstable_msg("FIP"),
    bip_3p_stability = unstable_msg("BIP")
  )
  list(
    f2_b2_distance = f2_b2,
    f2_b2_range = c(LAMP_F2B2_MIN, max_product),
    f2_b2_ok = unname(checks[["f2_b2_distance"]]),
    loop_length = unname(loops[["f2_b1"]]),
    loops = loops,
    loop_min = LAMP_LOOP_MIN,
    loop_ok = unname(checks[["loop_length"]]),
    fip_3p_clamp = unname(clamps[["FIP"]]),
    bip_3p_clamp = unname(clamps[["BIP"]]),
    fip_3p_fold = unname(folds[["FIP"]]),
    bip_3p_fold = unname(folds[["BIP"]]),
    fip_3p_dg = unname(dg[["FIP"]]),
    bip_3p_dg = unname(dg[["BIP"]]),
    fip_3p_stable = unname(checks[["fip_3p_stability"]]),
    bip_3p_stable = unname(checks[["bip_3p_stability"]]),
    checks = checks,
    ok = all(checks),
    failed = names(checks)[!checks],
    messages = unname(messages[!checks])
  )
}

#' The first loop primer (LF or LB) fitting a loop window
#'
#' Loop primers are optional accelerators: they anneal inside a single-stranded
#' loop and prime a new strand, speeding the reaction.  One is returned only when
#' the loop is long enough for an oligo in `LAMP_LOOP_LEN` that passes the inner
#' Tm / GC screen; otherwise `NULL` (the set still works without it).  A loop
#' primer extends in the same direction as the inner primer it assists, so LF is a
#' plus-strand oligo and LB a minus-strand one.
#' @param seq the template (plus strand).
#' @param lo,hi the loop window (plus-strand inclusive coordinates).
#' @param tm_range,gc_range accepted windows.
#' @param target_tm binding Tm to aim for.
#' @param reverse `TRUE` for LB (minus-strand), `FALSE` for LF.
#' @return a loop-primer sequence, or `NULL`.
#' @keywords internal
#' @noRd
lamp_loop <- function(seq, lo, hi, tm_range, gc_range, target_tm, reverse) {
  if (hi - lo + 1L < min(LAMP_LOOP_LEN)) {
    return(NULL)
  }
  cand <- lamp_candidates(seq, seq.int(lo, hi), LAMP_LOOP_LEN, tm_range, gc_range) |>
    dplyr::filter(.data$start >= lo, .data$end <= hi)
  if (!nrow(cand)) {
    return(NULL)
  }
  p <- cand |>
    dplyr::arrange(abs(.data$tm - target_tm)) |>
    dplyr::slice_head(n = 1L)
  primer <- if (reverse) lamp_revcomp(p$primer) else p$primer
  if (has_palindrome(primer)) {
    return(NULL)
  }
  primer
}

#' Joint compatibility of a primer pool, with a configurable dimer slack
#'
#' The pooled screen behind the LAMP cross-dimer gate:
#' [primer_set_compatibility()] over the named primers, with the cross-dimer run
#' threshold optionally raised by `slack` (the ladder's soft relaxation).
#' Non-ACGT marks (a probe abasic site) are treated as `N`.
#' @param slack extra base pairs allowed on each rule (the ladder's relaxation).
#' @param conditions reaction conditions for the reported dimer dG; defaults to
#'   `scheme_conditions("lamp")`.
#' @keywords internal
#' @noRd
lamp_pool_compatibility <- function(primers, slack = 0L, conditions = scheme_conditions("lamp")) {
  oligos <- unlist(primers, use.names = TRUE)
  oligos <- gsub("[^ACGTNacgtn]", "N", toupper(oligos))
  oligos <- oligos[nzchar(oligos)]
  oligos <- oligos[!duplicated(oligos)]
  if (length(oligos) < 2L) {
    return(NULL)
  }
  primer_set_compatibility(oligos,
    min_cross_run = LAMP_CROSS_RUN + slack, max_selfcomp = LAMP_SELFCOMP + slack,
    conditions = conditions
  )
}

#' Does FIP actually carry the variant at its 3' terminus?
#'
#' The designer's contract is that FIP's 3' arm (F2) ends at the variant and its
#' 3'-terminal base is the ALT allele.  That placement is enforced while the arm
#' is built, so the gate is a constructor invariant: it re-checks the returned
#' FIP rather than assuming the invariant held.  A drift in the placement (or a
#' caller-supplied FIP) therefore fails the gate instead of being accepted.
#' @param fip an inner-primer record (`arm3`, `primer`), or `NULL`.
#' @param alt the ALT base.
#' @param pos the variant position in the sequence the FIP was built on.
#' @return a single logical.
#' @keywords internal
#' @noRd
lamp_variant_gate <- function(fip, alt, pos) {
  ok <- !is.null(fip) &&
    length(fip$arm3) >= 2L &&
    !is.null(fip$primer) &&
    nzchar(fip$primer[1L])
  if (!ok) {
    return(FALSE)
  }
  end <- as.integer(fip$arm3[[2L]])
  terminal <- substr(fip$primer[1L], nchar(fip$primer[1L]), nchar(fip$primer[1L]))
  identical(end, as.integer(pos)) && identical(terminal, as.character(alt))
}

#' Design a LAMP marker at a variant
#'
#' The isothermal analogue of [assemble_kasp_assay()]: the outer F3/B3 pair and
#' the inner FIP/BIP pair are placed so the variant sits at the 3' terminus of
#' FIP's displacement arm (F2) - where a 3' mismatch blocks the strand-displacing
#' polymerase and so discriminates the alleles.  FIP and BIP are two-arm fusions
#' (`F1c+F2`, `B1c+B2`) built by `lamp_fip()` / `lamp_bip()`; F3/B3 reuse the
#' generic pair search over the inner block; the optional loop primers LF/LB are
#' placed in the F1c/F2 and B1c/B2 loops.  The designer does **not** simulate the
#' isothermal reaction (that needs a polymerase/kinetics model), so it returns a
#' well-formed primer set, not a predicted amplification curve.
#'
#' The returned set is **ALT-specific**: FIP's 3' base reads ALT.  Because a
#' REF-only template and a failed reaction are indistinguishable in an isothermal
#' readout (no amplification in either case), the design carries a `paired` block
#' (Z9.1): a REF set that is the same primer positions with FIP's 3' base flipped
#' to REF, so the two tubes are a deliberate pair, and `paired$status` says
#' whether it could be built.  The pair is a control construct, not a
#' discrimination model: no allele gap, Ct or call probability is predicted (see
#' `data-raw/proposed/lamp_status.decision.md`).
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param arm_lengths inner (FIP/BIP) arm lengths to try; each inner primer is a
#'   two-arm fusion, so this yields 40-48 nt inner primers by default (the
#'   published inner window is 45-49 nt), and widens only with `flank`.
#' @param outer_lengths F3/B3 outer primer lengths to try (21-24 nt by default).
#' @param outer_tm accepted F3/B3 Tm window (C).
#' @param inner_tm accepted inner binding-arm Tm window (C).
#' @param gc_range accepted GC fraction window (both arm classes).
#' @param min_product,max_product accepted amplicon window (bp).
#' @param loop_primers also place the optional LF/LB loop primers (default
#'   `TRUE`).
#' @return a `marker_lamp` list (`snp`, `primers`, `amplicon`, `loop`, `tm`,
#'   `variant_primer`, `gates`, `relaxed`, `compatibility` - the joint
#'   dimer/hairpin screen over the whole primer set), or `NULL` (with a message)
#'   when no primer set passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_lamp_marker(genome, marker_example("SNP")[1, ])
#' }
#' @family KASP LAMP
#' @export
design_lamp_marker <- function(genome, snp, flank = 150L, arm_lengths = 20:24,
                               outer_lengths = 21:24,
                               outer_tm = c(55, 62), inner_tm = c(55, 65),
                               gc_range = KASP_LIMITS$gc,
                               min_product = 120L, max_product = 300L,
                               loop_primers = TRUE) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  q <- tmpl$snp_i
  target_tm <- mean(inner_tm)
  n <- nchar(as.character(tmpl$seq))
  # the ladder: relax the soft gates LEAST grave first (C6), one rung at a time,
  # never below the floors above.  The variant placement, the inner Tm and the
  # sequence structure are hard - no set is returned past them.
  rungs <- list(
    character(0),
    "gate_lamp_amplicon",
    c("gate_lamp_amplicon", "gate_lamp_outer_tm"),
    c("gate_lamp_amplicon", "gate_lamp_outer_tm", "gate_lamp_cross_dimer")
  )
  attempt <- function(relaxed) {
    # FIP carries the mutation (its 3' arm ends at the variant); BIP sits
    # downstream so the two inner primers never overlap.
    fip <- lamp_fip(
      tmpl$seq, q, arm_lengths, inner_tm, gc_range, target_tm, LAMP_LOOP_MIN,
      alt_base = if (nchar(snp$alt[1L]) == 1L) snp$alt[1L] else NULL
    )
    if (is.null(fip)) {
      return(NULL)
    }
    bip <- lamp_bip(tmpl$seq, fip$arm3[2], arm_lengths, inner_tm, gc_range, target_tm, LAMP_LOOP_MIN)
    if (is.null(bip)) {
      return(NULL)
    }
    # the outer pair flanks the whole inner block, so F3/B3 never overlap it
    lo <- fip$arm5[1]
    hi <- bip$arm3[2]
    tm_win <- if ("gate_lamp_outer_tm" %in% relaxed) {
      outer_tm + c(-LAMP_OUTER_TM_RELAX, LAMP_OUTER_TM_RELAX)
    } else {
      outer_tm
    }
    amp_lo <- if ("gate_lamp_amplicon" %in% relaxed) min_product - LAMP_AMPLICON_RELAX else min_product
    amp_hi <- if ("gate_lamp_amplicon" %in% relaxed) max_product + LAMP_AMPLICON_RELAX else max_product
    fwd <- lamp_candidates(tmpl$seq, seq_len(max(1L, lo - 1L)), outer_lengths, tm_win, gc_range) |>
      dplyr::filter(.data$end < lo, !has_palindrome(.data$primer))
    rev <- lamp_candidates(tmpl$seq, seq.int(hi + 1L, n), outer_lengths, tm_win, gc_range, strand = "-") |>
      dplyr::filter(!has_palindrome(.data$primer))
    pair <- pick_primer_pair(fwd, rev, amp_lo, amp_hi)
    if (is.null(pair)) {
      return(NULL)
    }
    outer <- pair_to_list(pair)
    primers <- list(F3 = outer$forward, B3 = outer$reverse, FIP = fip$primer, BIP = bip$primer)
    # the optional loop primers, in the F1c/F2 and B1c/B2 loops
    if (loop_primers) {
      lf <- lamp_loop(tmpl$seq, fip$arm5[2] + 1L, fip$arm3[1] - 1L, inner_tm, gc_range, target_tm, reverse = FALSE)
      lb <- lamp_loop(tmpl$seq, bip$arm5[2] + 1L, bip$arm3[1] - 1L, inner_tm, gc_range, target_tm, reverse = TRUE)
      if (!is.null(lf)) primers$LF <- lf
      if (!is.null(lb)) primers$LB <- lb
    }
    # the whole set shares one tube, so screen the *joint* pool (cross-dimers and
    # hairpins between any two), not each primer alone.  The strict screen is the
    # one reported; the rung's `slack` widens it only to test acceptance.
    slack <- if ("gate_lamp_cross_dimer" %in% relaxed) LAMP_DIMER_SLACK else 0L
    compat <- lamp_pool_compatibility(primers, 0L)
    compatible <- is.null(compat) || all(lamp_pool_compatibility(primers, slack)$pass)
    # accept at this rung only if EVERY gate clears its (possibly widened) bar, so
    # the ladder is real: a design below the strict bar is returned only by the
    # rung that lowered that bar, never silently.
    gates_rung <- c(
      gate_lamp_variant = lamp_variant_gate(fip, snp$alt[1L], q),
      gate_lamp_inner_tm = all(dplyr::between(c(fip$tm, bip$tm), inner_tm[1L], inner_tm[2L])),
      gate_lamp_outer_tm = all(dplyr::between(c(outer$tm_f, outer$tm_r), tm_win[1L], tm_win[2L])),
      gate_lamp_structure = !any(has_palindrome(unlist(primers))),
      gate_lamp_cross_dimer = compatible,
      gate_lamp_amplicon = dplyr::between(outer$product, amp_lo, amp_hi)
    )
    if (!all(gates_rung)) {
      return(NULL)
    }
    list(fip = fip, bip = bip, outer = outer, primers = primers, compat = compat)
  }
  pick <- NULL
  for (r in rungs) {
    pick <- attempt(r)
    if (!is.null(pick)) {
      break
    }
  }
  if (is.null(pick)) {
    cli::cli_inform(gettext(
      "No LAMP primer set passed the variant-placement / inner-Tm / structure / cross-dimer / amplicon filters."
    ))
    return(NULL)
  }
  fip <- pick$fip
  bip <- pick$bip
  outer <- pick$outer
  primers <- pick$primers
  off <- snp$pos - q
  # the gate table the DAG displays, evaluated on the set actually returned (C7)
  gates <- c(
    # constructor invariant: FIP's 3' arm ends at the variant and reads ALT
    gate_lamp_variant = lamp_variant_gate(fip, snp$alt[1L], q),
    gate_lamp_inner_tm = all(
      dplyr::between(c(fip$tm, bip$tm), inner_tm[1L], inner_tm[2L])
    ),
    gate_lamp_outer_tm = all(
      dplyr::between(c(outer$tm_f, outer$tm_r), outer_tm[1L], outer_tm[2L])
    ),
    gate_lamp_structure = !any(has_palindrome(unlist(primers))),
    gate_lamp_cross_dimer = if (is.null(pick$compat)) NA else all(pick$compat$pass),
    gate_lamp_amplicon = dplyr::between(outer$product, min_product, max_product)
  )
  soft <- c("gate_lamp_amplicon", "gate_lamp_outer_tm", "gate_lamp_cross_dimer")
  relaxed <- soft[vapply(soft, function(g) isFALSE(gates[[g]]), logical(1))]
  # Geometry/physics diagnostics on the set actually returned.  A failure is
  # surfaced (warning + `geometry$ok = FALSE`) rather than silently accepted; it
  # is not folded into the ladder because FIP's 3'-terminal base is the ALT allele
  # and cannot always be an ideal G/C clamp.
  geo <- lamp_geometry(fip, bip, max_product = max_product)
  if (!isTRUE(geo$ok)) {
    cli::cli_warn(c(
      gettext("LAMP primer-set geometry check failed."),
      rlang::set_names(geo$messages, rep("*", length(geo$messages))),
      i = gettext("The set is returned with `geometry$ok = FALSE`; inspect the named check(s).")
    ))
  }
  out <- list(
    snp = snp,
    primers = primers,
    amplicon = tibble::tibble(start = outer$start + off, end = outer$end + off),
    # the F2-B2 loop the inner primers read, for reference
    loop = tibble::tibble(start = fip$arm3[2] + 1L + off, end = bip$arm5[1] - 1L + off),
    # genomic arm coordinates: FIP/BIP are fusions, only these arms align to the
    # reference and can carry sequence letters on the map
    lamp_arms = list(
      FIP = list(arm5 = fip$arm5 + off, arm3 = fip$arm3 + off),
      BIP = list(arm5 = bip$arm5 + off, arm3 = bip$arm3 + off)
    ),
    tm = c(F3 = outer$tm_f, B3 = outer$tm_r, FIP = fip$tm, BIP = bip$tm),
    variant_primer = "FIP",
    # Z9.1: the REF counterpart of the same set, so a REF-only sample has a
    # positive control instead of an ambiguous no-signal
    paired = lamp_paired_set(tmpl$seq, q, fip, snp, arm_lengths, inner_tm, gc_range, target_tm),
    gates = gates,
    relaxed = relaxed,
    geometry = geo,
    compatibility = pick$compat
  )
  structure(out, class = "marker_lamp")
}

#' The REF counterpart of an ALT-specific LAMP set (Z9.1)
#'
#' An isothermal readout gives no amplification for both a REF template and a
#' failed reaction, so an ALT-only set cannot tell them apart.  This rebuilds
#' FIP against the REF allele at the same 3' arm end, keeping every other primer
#' and the amplicon geometry, so a sample can be run against both tubes.
#'
#' It deliberately returns a *construct*, not a discrimination claim: the two
#' sets differ by the single 3' base of F2, and no allele gap is predicted.
#' @param seq the template sequence the ALT set was built on.
#' @param q the variant index in `seq`.
#' @param fip the ALT FIP from `lamp_fip()`.
#' @param snp the target row (`ref`, `alt`).
#' @return a list: `status` (`complete`, `ref_not_placed`, `not_applicable`),
#'   `ref_fip` (the REF primer, or `NA`), `discriminating_base` (the 3' base pair
#'   the two sets differ at), and `note`.
#' @keywords internal
#' @noRd
lamp_paired_set <- function(seq, q, fip, snp, arm_lengths, tm_range, gc_range, target_tm) {
  ref <- as.character(snp$ref[1L])
  alt <- as.character(snp$alt[1L])
  base <- list(
    ref_fip = NA_character_,
    discriminating_base = if (nchar(ref) == 1L && nchar(alt) == 1L) paste0(ref, "/", alt) else NA_character_
  )
  # only a SNV has a single discriminating 3' base to flip
  if (!nzchar(ref) || nchar(ref) != 1L || nchar(alt) != 1L) {
    return(c(base, list(
      status = "not_applicable",
      note = "the target is not a single-base substitution, so no single 3' base is flipped"
    )))
  }
  ref_fip <- lamp_fip(seq, q, arm_lengths, tm_range, gc_range, target_tm, LAMP_LOOP_MIN,
    alt_base = ref
  )
  if (is.null(ref_fip) || is.null(ref_fip$primer) || !nzchar(ref_fip$primer[1L])) {
    return(c(base, list(
      status = "ref_not_placed",
      note = "FIP could not be rebuilt on the REF allele at this 3' arm end; run a separate REF design"
    )))
  }
  # explicit list, not c(base, list(...)): c() on two lists can silently drop an
  # earlier NA element, which made a built primer read as missing
  list(
    status = "complete",
    ref_fip = as.character(ref_fip$primer[1L]),
    discriminating_base = base$discriminating_base,
    note = "same primer positions, FIP 3' base reads REF; the two tubes are a control pair, no allele gap is predicted"
  )
}

#' @export
print.marker_lamp <- function(x, ...) {
  cat(
    "<marker_lamp>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  amplicon:", x$amplicon$end - x$amplicon$start + 1L, "bp\n")
  cat("  variant at the FIP 3' end (F2 ends at the variant)\n")
  if (!is.null(x$paired)) {
    cat("  paired REF set:", x$paired$status, "\n")
    if (identical(x$paired$status, "complete")) {
      cat("    REF FIP:", x$paired$ref_fip, "\n")
    }
  }
  cat("  F3:", x$primers$F3, "\n  B3:", x$primers$B3, "\n")
  cat("  FIP:", x$primers$FIP, "\n  BIP:", x$primers$BIP, "\n")
  if ("LF" %in% names(x$primers)) cat("  LF:", x$primers$LF, "\n")
  if ("LB" %in% names(x$primers)) cat("  LB:", x$primers$LB, "\n")
  if (length(x$relaxed)) {
    cat("  (relaxed:", paste(x$relaxed, collapse = ", "), ")\n")
  }
  if (!is.null(x$geometry) && !isTRUE(x$geometry$ok)) {
    cat("  geometry FAILED:", paste(x$geometry$failed, collapse = ", "), "\n")
  }
  invisible(x)
}
