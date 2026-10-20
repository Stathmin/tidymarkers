# =============================================================================
#  marker_rpa.R  --  RPA (recombinase polymerase amplification), isothermal
# =============================================================================
#  RPA is the isothermal cousin of PCR: a recombinase pairs LONG primers with
#  duplex DNA and a strand-displacing polymerase extends them, so the reaction
#  runs at a constant 37-42 C with no thermal cycling.  Two consequences drive
#  the design rules here and separate RPA from every other scheme in the package:
#    1. Primer LENGTH, not Tm, is the governing parameter.  The recombinase needs
#       a long filament, so primers are ~30-35 nt (a short PCR primer recombines
#       poorly and the reaction is slow/insensitive); Tm is not a gate.
#    2. The amplicon is kept SHORT (~100-200 bp), the range RPA amplifies best.
#  Allele discrimination is not as sharp as KASP's, so the variant is placed at
#  a primer's 3'-terminal base (no probe) or under the optional exo probe.  A
#  variant *inside* a primer (not at its 3' end) does not block extension and
#  would not discriminate the alleles, so it is rejected rather than accepted.  The exo probe is
#  an oligo blocked at its 3' end with a tetrahydrofuran (THF / dSpacer) abasic
#  site flanked by a dT-fluorophore and a dT-quencher; Exonuclease III cleaves
#  the THF only in a probe-target duplex, separating the pair and releasing
#  fluorescence.
#
#  Grounding (design rules verified against the primary method sources; the
#  attempts to reach duckduckgo.com's HTML search endpoint returned a bot
#  challenge, so the primary literature was read directly through Europe PMC):
#   - Lobato IM, O'Sullivan CK, Trends Analyt Chem 2018, 98:19-35
#     (PMC7112910, https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7112910/
#     fullTextXML): "initially believed that ... primers of 30-35 bases ... are
#     necessary"; "longer primers (up to 45 nucleotides) can be used, but they
#     could lead to secondary structures"; "GC content below 30% or above 70% is
#     not recommended"; "better suited to amplicons between 100 and 200 bp"; the
#     Exo probe "has a dT-fluorophore and a dT-quencher flanking a
#     tetrahydrofuran residue (dSpacer), which are separated by a maximum of 2-4
#     bases"; "operating at 37-42 C"; and (nfo probe) "46-52 oligonucleotide ...
#     internal abasic ... placed at least 30 nucleotides from the 5' end and 15
#     nucleotides from the 3' end".
#   - Piepenburg O, et al., PLoS Biol 2006, 4(7):e204 (doi:10.1371/journal.pbio.
#     0040204): the original RPA report (recombinase/SSB/strand-displacing
#     polymerase, no denaturation step).
#   - Wikipedia, "Recombinase polymerase amplification" (fetched 2026,
#     https://en.wikipedia.org/wiki/Recombinase_polymerase_amplification):
#     "Typically 30-38 base primers are needed for efficient recombinase filament
#     formation"; RPA "Operating best at temperatures of 37-42 C".
#  Thresholds the sources do not state numerically are marked as ASSUMPTIONS
#  below and are parameters the caller can override.
# =============================================================================

# The primer LENGTH window (RPA needs a long recombinase filament; Tm is not the
# gate here, so the candidate Tm window is deliberately open).  Lobato 2018 and
# Wikipedia both give ~30-35 nt (Wikipedia allows to 38); we use the review's
# 30-35 and note longer oligos are tolerated but risk secondary structure.
RPA_PRIMER_LEN <- 30:35

# RPA is "better suited to amplicons between 100 and 200 bp" (Lobato 2018); the
# reaction can reach ~1.5 kb but short products amplify faster and more evenly.
RPA_AMPLICON <- c(100L, 200L)

# "A GC content below 30% or above 70% is not recommended" (Lobato 2018).
RPA_GC <- c(0.30, 0.70)

# Tm does not gate RPA primers (length does), so the candidate window is open;
# it only exists because primer_candidates() is the shared composition screener.
RPA_NO_TM <- c(0, 200)

# Exo probes follow the abasic-probe length and geometry of the TwistAmp probe
# (Lobato 2018 gives 46-52 nt and the >=30 nt 5' / >=15 nt 3' abasic-arm rule for
# the nfo probe; applying the same geometry to the exo probe is an ASSUMPTION).
RPA_PROBE_LEN <- 46:52
RPA_PROBE_THF_5P <- 30L
RPA_PROBE_THF_3P <- 15L
# The dT-fluorophore and dT-quencher flank the THF within a maximum of 2-4 bases
# (Lobato 2018); we require a dT within this gap on each side of the THF.
RPA_PROBE_FQ_GAP <- 4L

# Without a probe the variant is read from a primer's 3'-terminal base (the
# allele-specific-extension logic: a 3' mismatch blocks the strand-displacing
# polymerase).  It must sit exactly at the 3' terminus; further in and the
# mismatch no longer stops extension, so the window is 0.
RPA_3P_WINDOW <- 0L

# The relaxation floors for the soft gates (bounded; C8).  A primer a few nt
# outside the governing length window still recombines (Wikipedia allows to 38
# nt), the amplicon may stretch, and the THF / fluorophore-quencher geometry may
# drift a little - but never without bound.
RPA_PRIMER_LEN_RELAX <- c(min(RPA_PRIMER_LEN) - 3L, max(RPA_PRIMER_LEN) + 3L)
RPA_AMPLICON_RELAX <- 20L
RPA_PROBE_THF_5P_RELAX <- 25L
RPA_PROBE_THF_3P_RELAX <- 12L
RPA_PROBE_FQ_GAP_RELAX <- 6L

# The isothermal operating window (Lobato 2018; Wikipedia): recorded for the
# readout, not gated - RPA needs no thermal cycling.
RPA_TEMP <- c(37, 42)

# The RPA decision stage, mirroring the CAPS/ARMS/HRM/castPCR stage constants in
# R/marker_select.R.  `scheme_def_rpa()` returns them so the shared registry can
# wire the chemistry in afterwards without a second vocabulary.
KASP_GATES_RPA <- c(
  "gate_rpa_variant", "gate_rpa_primer_len", "gate_rpa_structure",
  "gate_rpa_probe", "gate_rpa_amplicon"
)
# relaxation order, LEAST grave first (C6): the amplicon window is the most
# forgiving, then the probe geometry, then the governing primer length.
KASP_RELAX_RPA <- c("gate_rpa_amplicon", "gate_rpa_probe", "gate_rpa_primer_len")
KASP_TIERS_RPA <- c(
  tier1 = "discrimination placement", tier2 = "primer length",
  tier3 = "probe geometry", tier4 = "amplicon length"
)

# An empty primer_candidates() result (every candidate carries a homopolymer run or
# an SSR) makes the shared screener fail in primer_tm(); RPA probes very narrow
# windows, so guard it here and treat "no candidates" as an empty grid.
rpa_candidates <- function(seq, starts, lengths, strand, tm_range, gc_range) {
  empty <- tibble::tibble(
    start = integer(0), len = integer(0), end = integer(0),
    primer = character(0), tm = numeric(0), gc = numeric(0)
  )
  out <- tryCatch(
    primer_candidates(seq, starts, lengths, strand, tm_range, gc_range),
    error = function(e) NULL
  )
  if (is.null(out)) empty else out
}

#' Do an RPA primer pair pass the structure / dimer filters?
#'
#' The recombinase cannot work through strong primer secondary structure or
#' dimers, so each oligo must be palindrome-free and the hairpin / homodimer /
#' cross-dimer dG37 must clear the shared KASP thresholds.  When OligoArrayAux is
#' unavailable the check is skipped (the metric cannot be computed) rather than
#' failing the design.  Reuses [has_palindrome()] and [dimer_thermo()].
#' @param pair a pair list from `pair_to_list()` (`forward`, `reverse`).
#' @param BPPARAM optional `BiocParallel` backend for [dimer_thermo()].
#' @return a single logical.
#' @keywords internal
#' @noRd
rpa_structure_ok <- function(pair, BPPARAM = NULL) {
  if (!has_oligoarrayaux()) {
    return(TRUE)
  }
  p <- c(forward = pair$forward, reverse = pair$reverse)
  if (any(has_palindrome(p))) {
    return(FALSE)
  }
  dg <- dimer_thermo(p, BPPARAM = BPPARAM)
  self <- dg$dg[dg$a == dg$b]
  cross <- dg$dg[dg$a != dg$b]
  all(self >= KASP_LIMITS$homodimer_dg) && all(cross >= KASP_LIMITS$duplex_dg)
}

#' Pick an RPA primer pair (long primers; variant under a primer or flanked)
#'
#' Mirrors `design_site_pair()` / `design_flanking_pair()` but pre-filters the
#' candidate primers for palindromes (an inverted repeat is a structure defect the
#' recombinase cannot work through), which those shared helpers do not screen.
#' Reuses `primer_candidates()`, `pick_primer_pair()` and `pair_to_list()`.
#' Without a probe the variant sits under one primer (a site pair); with a probe
#' the primers flank the variant.
#' @param seq the template (plus strand).
#' @param si 1-based variant position on `seq`.
#' @param design_probe flank the variant (probe mode) rather than cover it?
#' @param lengths primer lengths to try.
#' @param min_product,max_product the accepted amplicon window.
#' @param tm_range,gc_range accepted windows.
#' @return a pair list (`pair_to_list()`), or `NULL`.
#' @keywords internal
#' @noRd
rpa_pair <- function(seq, si, design_probe, lengths, min_product, max_product,
                     tm_range, gc_range) {
  n <- nchar(as.character(seq))
  if (si <= 1L || si >= n) {
    return(NULL)
  }
  clean <- function(d) d[!has_palindrome(d$primer), , drop = FALSE]
  if (design_probe) {
    # with a probe the primers flank the variant and the probe carries the
    # discrimination, so the primers only need to bracket it
    fwd <- rpa_candidates(seq, seq_len(si - 1L), lengths, "+", tm_range, gc_range) |>
      dplyr::filter(.data$end < si) |>
      clean()
    rev <- rpa_candidates(seq, seq.int(si + 1L, n), lengths, "-", tm_range, gc_range) |>
      clean()
    b <- pick_primer_pair(fwd, rev, min_product, max_product)
    return(if (is.null(b)) NULL else pair_to_list(b))
  }
  # no probe: the variant must sit at a primer's 3' terminus (the AS-like
  # extension block), within RPA_3P_WINDOW nt of it - an internal mismatch does
  # not stop the strand-displacing polymerase and would not discriminate.
  w <- RPA_3P_WINDOW
  near_f <- seq_len(si)
  near_r <- seq.int(max(1L, si - max(lengths) + 1L), si)
  # the variant at the 3' end of the forward primer (its `end`) or the reverse
  # primer (whose 3' terminus is its plus-strand `start`)
  cover_f <- rpa_candidates(seq, near_f, lengths, "+", tm_range, gc_range) |>
    dplyr::filter(.data$end >= si, .data$end <= si + w) |>
    clean()
  cover_r <- rpa_candidates(seq, near_r, lengths, "-", tm_range, gc_range) |>
    dplyr::filter(.data$start >= si - w) |>
    clean()
  # the mates stay clear of the covering primer so they cannot overlap it
  mate_r <- rpa_candidates(seq, seq.int(si + w + 1L, n), lengths, "-", tm_range, gc_range) |>
    clean()
  mate_f <- rpa_candidates(seq, seq_len(si - w - 1L), lengths, "+", tm_range, gc_range) |>
    dplyr::filter(.data$end < si - w) |>
    clean()
  cand <- dplyr::bind_rows(
    pick_primer_pair(cover_f, mate_r, min_product, max_product),
    pick_primer_pair(mate_f, cover_r, min_product, max_product)
  )
  if (!nrow(cand)) {
    return(NULL)
  }
  pair_to_list(cand[which.min(abs(cand$tm_f - cand$tm_r)), ])
}

#' One exo probe candidate: place the THF and the dT-fluorophore/quencher
#'
#' Keeps the abasic (THF) site inside the published 5'/3' arm windows, replaces
#' that base with `"X"`, and verifies a dT within `fq_gap` bases on each side (the
#' dT-fluorophore and dT-quencher of the TwistAmp exo probe).  Returns `NULL`
#' when the geometry cannot be satisfied.  The arm/gap windows are parameters so
#' the ladder can relax them (they are soft).
#' @param row one row of `primer_candidates()` output.
#' @param si 1-based variant position on the template.
#' @param thf_5p,thf_3p the minimum 5'/3' abasic-arm lengths.
#' @param fq_gap the maximum dT-fluorophore/quencher gap on each side.
#' @return a probe list, or `NULL`.
#' @keywords internal
#' @noRd
rpa_probe_one <- function(row, si, thf_5p = RPA_PROBE_THF_5P,
                          thf_3p = RPA_PROBE_THF_3P, fq_gap = RPA_PROBE_FQ_GAP) {
  len <- row$len
  lo <- thf_5p + 1L
  hi <- len - thf_3p
  if (lo > hi) {
    return(NULL)
  }
  voff <- si - row$start + 1L
  thf <- min(max(voff, lo), hi)
  # THF on the variant would replace the allele with X and destroy the readout.
  # Move it one base aside when the arm window allows it; otherwise reject.
  if (identical(thf, voff)) {
    thf <- if (voff < hi) voff + 1L else if (voff > lo) voff - 1L else NA_integer_
    if (is.na(thf)) {
      return(NULL)
    }
  }
  left <- stringr::str_locate_all(stringr::str_sub(row$primer, 1L, thf - 1L), "T")[[1L]]
  right <- stringr::str_locate_all(stringr::str_sub(row$primer, thf + 1L, len), "T")[[1L]]
  if (!nrow(left) || !nrow(right)) {
    return(NULL)
  }
  gap5 <- thf - max(left[, "start"])
  gap3 <- min(right[, "start"])
  if (gap5 > fq_gap || gap3 > fq_gap) {
    return(NULL)
  }
  list(
    sequence = paste0(
      stringr::str_sub(row$primer, 1L, thf - 1L), "X",
      stringr::str_sub(row$primer, thf + 1L, len)
    ),
    start = row$start, end = row$end, len = len,
    thf_index = thf, thf_arm5 = thf - 1L, thf_arm3 = len - thf,
    fq_gap5 = gap5, fq_gap3 = gap3,
    tm = row$tm, gc = row$gc, blocked_3p = TRUE
  )
}

#' The first exo probe covering the variant inside the amplicon
#'
#' Enumerates plus-strand probes that span the variant and fit *between* the
#' flanking primers' 3' ends, screened for GC/composition by
#' `primer_candidates()`, then keeps the first that satisfies the THF /
#' fluorophore-quencher geometry (`rpa_probe_one()`), preferring a longer probe
#' and the THF nearest the variant.
#' @param seq the template (plus strand).
#' @param si 1-based variant position on `seq`.
#' @param pair the amplicon-bounding pair (`pair_to_list()`).
#' @param probe_len probe lengths to try.
#' @param gc_range accepted GC window.
#' @param thf_5p,thf_3p,fq_gap the (soft) geometry windows (see `rpa_probe_one()`).
#' @return a probe list, or `NULL`.
#' @keywords internal
#' @noRd
rpa_probe <- function(seq, si, pair, probe_len, gc_range,
                      thf_5p = RPA_PROBE_THF_5P, thf_3p = RPA_PROBE_THF_3P,
                      fq_gap = RPA_PROBE_FQ_GAP) {
  starts <- seq.int(max(1L, si - max(probe_len) + 1L), si)
  # the probe must not overlap either primer: stay inside the amplicon bases
  # between the flanking primers' 3' ends.
  inner_lo <- pair$start + nchar(pair$forward)
  inner_hi <- pair$end - nchar(pair$reverse)
  cand <- rpa_candidates(seq, starts, probe_len, "+", RPA_NO_TM, gc_range) |>
    dplyr::filter(
      .data$start >= inner_lo, .data$end <= inner_hi,
      .data$start <= si, .data$end >= si
    ) |>
    dplyr::arrange(dplyr::desc(.data$len), .data$start)
  probe <- purrr::map(seq_len(nrow(cand)), function(i) {
    rpa_probe_one(cand[i, ], si, thf_5p, thf_3p, fq_gap)
  })
  purrr::keep(probe, ~ !is.null(.x)) |>
    purrr::pluck(1L)
}

#' Did the exo probe clear the declared THF/FQ geometry?
#'
#' `TRUE` for a probe that passed at the strict geometry rung, `FALSE` when it
#' only appeared at the relaxed rung (so `relaxed` records `gate_rpa_probe`),
#' and `NA` when no probe was designed.  Kept separate from
#' [design_rpa_marker()] so the ladder's bookkeeping is directly testable.
#' @param probe the chosen probe list, or `NULL`.
#' @param probe_relaxed did the successful rung relax the probe geometry?
#' @return a single logical (or `NA`).
#' @keywords internal
#' @noRd
rpa_probe_gate <- function(probe, probe_relaxed) {
  if (is.null(probe)) {
    return(NA)
  }
  !isTRUE(probe_relaxed)
}

#' Is the variant actually at a discriminating position?
#'
#' Without a probe RPA discriminates through a primer's 3'-terminal base, so the
#' variant must sit within [RPA_3P_WINDOW] nt of that 3' end.  With an exo probe
#' the probe carries the allele: it must span the variant, and its THF must not
#' replace the variant base.  `design_rpa_marker()` enforces this while it picks
#' the pair/probe, so the gate re-checks the returned coordinates rather than
#' returning an unconditional `TRUE`.
#' @param under `"forward"`, `"reverse"`, `"probe"` (or anything else).
#' @param pos variant position (same coordinate system as the others).
#' @param fwd_end forward primer 3' end.
#' @param rev_start reverse primer 3' end (plus-strand start).
#' @param probe the exo probe, or `NULL`.
#' @param window allowed distance from a primer 3' end.
#' @return a single logical.
#' @keywords internal
#' @noRd
rpa_variant_gate <- function(under, pos, fwd_end = NA_real_, rev_start = NA_real_,
                             probe = NULL, window = RPA_3P_WINDOW) {
  if (identical(under, "probe")) {
    return(
      !is.null(probe) &&
        probe$start[1L] <= pos && probe$end[1L] >= pos &&
        !identical(as.integer(probe$thf_index[1L]), as.integer(pos - probe$start[1L] + 1L))
    )
  }
  if (identical(under, "forward")) {
    return(!is.na(fwd_end) && abs(fwd_end - pos) <= window)
  }
  if (identical(under, "reverse")) {
    return(!is.na(rev_start) && abs(rev_start - pos) <= window)
  }
  FALSE
}

#' Design an RPA (recombinase polymerase amplification) marker at a variant
#'
#' The isothermal analogue of a plain PCR pair: two LONG primers (default 30-35
#' nt - the recombinase needs a long filament, so length and not Tm is the gate)
#' bound a SHORT amplicon (default 100-200 bp), run at a constant 37-42 C.  The
#' variant is placed for allele discrimination either under one primer
#' (`design_probe = FALSE`) or under the optional exo probe
#' (`design_probe = TRUE`, in which case the primers flank the variant).  The exo
#' probe carries a THF abasic site flanked by a dT-fluorophore and a dT-quencher
#' and is 3'-blocked for real-time readout.  Single-base substitutions only.
#'
#' RPA discriminates alleles less sharply than KASP, so this is a placement plus a
#' well-formed probe geometry, not a predicted efficiency; the structure/dimer
#' filter is skipped when OligoArrayAux is unavailable.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param primer_len primer lengths to try (RPA wants long primers, not a Tm).
#' @param min_product,max_product the accepted amplicon window.
#' @param design_probe also design an exo probe covering the variant?  When
#'   `TRUE` the primers flank the variant and the probe carries the allele.
#' @param probe_len exo-probe lengths to try.
#' @param gc_range accepted GC window for every oligo.
#' @param tm_range accepted primer Tm window; open by default (Tm is not the RPA
#'   gate) and only used by the shared candidate screener.
#' @param BPPARAM optional `BiocParallel` backend for the dimer filter.
#' @return a `marker_rpa` list (`snp`, `primers`, `probe`, `amplicon`,
#'   `compatibility` - the joint dimer/hairpin screen over the pair and the probe -
#'   `gates` - the evaluated gate table - and `relaxed`, the soft gates it had to
#'   lower), or `NULL` (with a message) when no pair or probe passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_rpa_marker(genome, marker_example("SNP")[1, ])
#' design_rpa_marker(genome, marker_example("SNP")[1, ], design_probe = TRUE)
#' }
#' @family KASP RPA
#' @export
design_rpa_marker <- function(genome, snp, flank = 200L,
                              primer_len = RPA_PRIMER_LEN,
                              min_product = RPA_AMPLICON[1L],
                              max_product = RPA_AMPLICON[2L],
                              design_probe = FALSE,
                              probe_len = RPA_PROBE_LEN,
                              gc_range = RPA_GC,
                              tm_range = RPA_NO_TM,
                              BPPARAM = NULL) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    cli::cli_inform(gettext("RPA discrimination here is implemented for single-base substitutions only."))
    return(NULL)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  si <- tmpl$snp_i
  strict_len <- as.integer(primer_len)
  # the ladder: relax the soft gates LEAST grave first (C6), one rung at a time,
  # never below the floors above.  The variant placement (`gate_rpa_variant`) and
  # the structure/dimers (`gate_rpa_structure`) are hard - no pair is returned
  # past them, at any rung.
  rungs <- list(
    character(0),
    "gate_rpa_amplicon",
    c("gate_rpa_amplicon", "gate_rpa_probe"),
    c("gate_rpa_amplicon", "gate_rpa_probe", "gate_rpa_primer_len")
  )
  attempt <- function(relaxed) {
    lengths <- if ("gate_rpa_primer_len" %in% relaxed) {
      seq.int(min(strict_len) - 3L, max(strict_len) + 3L)
    } else {
      strict_len
    }
    amp_lo <- if ("gate_rpa_amplicon" %in% relaxed) min_product - RPA_AMPLICON_RELAX else min_product
    amp_hi <- if ("gate_rpa_amplicon" %in% relaxed) max_product + RPA_AMPLICON_RELAX else max_product
    pair <- rpa_pair(tmpl$seq, si, design_probe, lengths, amp_lo, amp_hi, tm_range, gc_range)
    if (is.null(pair) || !rpa_structure_ok(pair, BPPARAM)) {
      return(NULL)
    }
    probe <- NULL
    if (design_probe) {
      thf5 <- if ("gate_rpa_probe" %in% relaxed) RPA_PROBE_THF_5P_RELAX else RPA_PROBE_THF_5P
      thf3 <- if ("gate_rpa_probe" %in% relaxed) RPA_PROBE_THF_3P_RELAX else RPA_PROBE_THF_3P
      fq <- if ("gate_rpa_probe" %in% relaxed) RPA_PROBE_FQ_GAP_RELAX else RPA_PROBE_FQ_GAP
      probe <- rpa_probe(tmpl$seq, si, pair, probe_len, gc_range, thf5, thf3, fq)
      if (is.null(probe)) {
        return(NULL)
      }
    }
    list(pair = pair, probe = probe, probe_relaxed = "gate_rpa_probe" %in% relaxed)
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
      "No RPA primer pair passed the variant-placement / length / structure / amplicon filters."
    ))
    return(NULL)
  }
  pair <- pick$pair
  probe <- pick$probe
  off <- snp$pos - si
  if (!is.null(probe)) {
    probe$start <- probe$start + off
    probe$end <- probe$end + off
  }
  fwd_end <- pair$start + nchar(pair$forward) - 1L
  rev_start <- pair$end - nchar(pair$reverse) + 1L
  under <- if (design_probe) {
    "probe"
  } else if (abs(fwd_end - si) <= RPA_3P_WINDOW) {
    "forward"
  } else if (abs(rev_start - si) <= RPA_3P_WINDOW) {
    "reverse"
  } else {
    "none" # no discriminating placement: the gate below reports the failure
  }
  # the gate table the DAG displays, evaluated on the assay actually returned (C7)
  gates <- c(
    # constructor invariant: the variant is at a primer 3' end or under the probe
    gate_rpa_variant = rpa_variant_gate(under, si, fwd_end, rev_start, probe),
    gate_rpa_primer_len = all(
      dplyr::between(c(nchar(pair$forward), nchar(pair$reverse)), min(strict_len), max(strict_len))
    ),
    gate_rpa_structure = if (has_oligoarrayaux()) TRUE else NA,
    gate_rpa_probe = rpa_probe_gate(probe, pick$probe_relaxed),
    gate_rpa_amplicon = dplyr::between(pair$product, min_product, max_product)
  )
  soft <- c("gate_rpa_amplicon", "gate_rpa_probe", "gate_rpa_primer_len")
  relaxed <- soft[vapply(soft, function(g) isFALSE(gates[[g]]), logical(1))]
  # every coordinate exposed by the object is plus-strand/genomic (the probe and
  # amplicon below are already shifted); shift the pair's local coords to match
  primers <- pair
  primers$start <- pair$start + off
  primers$end <- pair$end + off
  out <- list(
    snp = snp,
    primers = c(primers, list(
      length_f = nchar(pair$forward), length_r = nchar(pair$reverse),
      variant_under = under, temperature = RPA_TEMP
    )),
    probe = probe,
    gates = gates,
    relaxed = relaxed,
    amplicon = tibble::tibble(start = pair$start + off, end = pair$end + off)
  )
  # the pair and the exo probe share the tube: joint screen (the probe's THF base
  # is treated as N by tube_compatibility)
  out$compatibility <- tube_compatibility(c(
    forward = out$primers$forward, reverse = out$primers$reverse,
    probe = out$probe$sequence
  ))
  structure(out, class = "marker_rpa")
}

#' @export
print.marker_rpa <- function(x, ...) {
  cat("<marker_rpa>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt, "\n", sep = "")
  cat("  F:", x$primers$forward, " (", x$primers$length_f, " nt)\n", sep = "")
  cat("  R:", x$primers$reverse, " (", x$primers$length_r, " nt)\n", sep = "")
  if (!is.null(x$probe)) {
    cat("  exo probe:", x$probe$sequence, "\n")
    cat(
      "  (THF at ", x$probe$thf_index, ": ", x$probe$thf_arm5, " nt 5' / ",
      x$probe$thf_arm3, " nt 3', 3'-blocked)\n",
      sep = ""
    )
  }
  locus <- switch(x$primers$variant_under,
    probe = "under the exo probe",
    forward = "at the forward primer 3' end",
    reverse = "at the reverse primer 3' end",
    "not at a discriminating position"
  )
  cat(
    "  amplicon:", x$primers$product, "bp; variant ", locus,
    "; isothermal ", x$primers$temperature[1L], "-",
    x$primers$temperature[2L], " C\n",
    sep = ""
  )
  if (length(x$relaxed)) {
    cat("  (relaxed:", paste(x$relaxed, collapse = ", "), ")\n")
  }
  invisible(x)
}

#' The RPA scheme definition for the shared registry
#'
#' A self-describing definition of the RPA chemistry in the shape the marker
#' registry expects (scheme/family/label/readout, the gate/tier vocabulary and the
#' per-gate human-readable labels).  Kept beside the designer so the chemistry and
#' its thresholds travel together; it is discovered automatically (one file per
#' scheme, see `R/marker_scheme_registry.R`).
#' @return a named list describing the scheme.
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_rpa
#' @export
#' @noRd
marker_map_spec.marker_rpa <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  probe <- if (!is.null(marker$probe)) {
    list(spec_row(
      "exo probe", marker$probe$start[1L], marker$probe$end[1L], "+", "probe",
      label = "exo probe"
    ))
  } else {
    list()
  }
  map_spec(
    marker,
    features = c(map_pair_rows(marker, sp[["start"]], sp[["end"]]), probe),
    sequences = list(
      forward = marker$primers$forward, reverse = marker$primers$reverse,
      "exo probe" = if (!is.null(marker$probe)) marker$probe$sequence else NULL
    )
  )
}

scheme_def_rpa <- function() {
  primer_len <- RPA_PRIMER_LEN
  amplicon <- RPA_AMPLICON
  probe_len <- RPA_PROBE_LEN
  primer_lo <- primer_len[1L]
  primer_hi <- primer_len[length(primer_len)]
  probe_lo <- probe_len[1L]
  probe_hi <- probe_len[length(probe_len)]
  amp_lo <- amplicon[1L]
  amp_hi <- amplicon[2L]
  list(
    scheme = "rpa",
    family = "rpa",
    label = "RPA (recombinase polymerase amplification, isothermal)",
    discrimination = "variant at a long primer 3' end, or under the exo probe",
    readout = "isothermal 37-42 C; real-time exo probe (or gel)",
    zygosity = "co-dominant",
    amplicon = "short (100-200 bp)",
    loci = 1L,
    verdict = "usable RPA assay",
    gates = KASP_GATES_RPA,
    relax = KASP_RELAX_RPA,
    tiers = KASP_TIERS_RPA,
    gate_short = c(
      gate_rpa_variant = "variant placement",
      gate_rpa_primer_len = "primer length",
      gate_rpa_structure = "primer structure/dimers",
      gate_rpa_probe = "probe geometry",
      gate_rpa_amplicon = "amplicon length (RPA)"
    ),
    gate_limit = c(
      gate_rpa_variant = "variant at a primer 3' end (no probe) or under the exo probe",
      gate_rpa_primer_len = glue::glue("each primer within [{primer_lo}, {primer_hi}] nt (length, not Tm)"),
      gate_rpa_structure = glue::glue("no palindrome; dimer dG37 >= {KASP_LIMITS$homodimer_dg}"),
      gate_rpa_probe = paste0(
        "probe ", probe_lo, "-", probe_hi, " nt; THF >= ", RPA_PROBE_THF_5P,
        " nt 5' / >= ", RPA_PROBE_THF_3P, " nt 3'; variant under probe (when used)"
      ),
      gate_rpa_amplicon = glue::glue("amplicon within [{amp_lo}, {amp_hi}] bp")
    ),
    gate_floor_label = c(
      gate_rpa_amplicon = "widen the amplicon window by 20 bp",
      gate_rpa_probe = "relax the THF / fluorophore-quencher geometry",
      gate_rpa_primer_len = "allow up to 3 nt outside the length window"
    ),
    gate_hard_label = c(
      gate_rpa_variant = "a variant away from a primer 3' end / the probe is not discriminated",
      gate_rpa_structure = "strong primer dimers/hairpins stall the recombinase"
    ),
    stage_label = "RPA primers/probe"
  )
}
