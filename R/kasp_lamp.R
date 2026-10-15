# =============================================================================
#  kasp_lamp.R  --  LAMP (loop-mediated isothermal amplification) as a marker
# =============================================================================
#  LAMP is the isothermal, primer-heavy sibling of PCR: a strand-displacing
#  polymerase (Bst) runs at a constant 60-65 C, and 4-6 primers recognise 6-8
#  regions of a short target.  Four are mandatory - the outer F3/B3 pair and the
#  inner FIP/BIP pair (each inner primer is a fusion of a loop arm and a
#  displacement arm: FIP = F1c + F2, BIP = B1c + B2); two loop primers (LF/LB)
#  are optional accelerators.  Because there is no thermal cycling, the assay
#  lives or dies on primer placement and thermodynamics: F3/B3 anneal first
#  (~55-60 C), the inner primers drive the cycling (~60-65 C), and the variant is
#  read from the loop that forms between F2 and B2 - so it must lie *there*, at
#  the 5'/3' ends of FIP/BIP (or under a loop primer), never under the outer
#  primers.
#
#  Design rules used here (web-researched; see provenance below):
#    * 4 primers (F3/B3/FIP/BIP); amplicon 120-300 bp - "optimal 130-200 bp,
#      sequences up to 300 bp are also acceptable (Notomi 2000)"
#      (Shirshikov et al. 2019, PeerJ 7:e6801, PMC6487805; MorphoCatcher);
#      NEB quotes a 80-250 nt target and 4-6 primers
#      (https://www.neb.com/.../loop-mediated-isothermal-amplification-lamp).
#    * isothermal 60-65 C; inner FIP/BIP bind two distant sites and are longer
#      (45-49 bp) than the outer F3/B3 (21-24 bp); the primers must not form a
#      stable double strand at that temperature
#      (Soroka, Wasowicz & Rymaszewska 2021, Cells 10(8):1931, PMC8393631;
#      Notomi et al. 2000; Nagamine et al. 2002 for loop primers).
#    * the variant is placed at the 5'/3' ends of the inner primers (and the 3'
#      ends of the outer / loop primers) - "the most suitable regions ... are the
#      5'- and 3'-ends of inner (FIP and BIP) primers, as well as 3'-ends of
#      outer (F3 and B3) and loop (LF and LB) primers"
#      (Shirshikov et al. 2019, PeerJ 7:e6801).
#
#  Threshold provenance / assumptions:
#    * outer Tm window c(55, 62) C and inner (binding-arm) Tm window c(60, 65) C
#      are the rule-of-thumb bands above rounded to the package's TmCalculator
#      conditions (PCR_COND, DMSO-corrected).  If they were authored from a
#      different salt model, treat them as an explicit assumption, not a fact.
#    * GC 40-60 % and "no run >= 4 / no SSR" are the primer-agnostic KASP limits
#      (KASP_LIMITS), reused unchanged (see kasp_primer_pairs.R).
#    * the structural screen is sequence-level (no inverted repeat), reusing
#      has_palindrome(); the whole four-primer set is then screened jointly for
#      cross-dimers/hairpins by tube_compatibility() (sequence-level run screen,
#      no OligoArrayAux needed).
# =============================================================================

# The LAMP decision stage, declared here so scheme_def_lamp() is self-contained.
# (The coordinator wires this into kasp_schemes()/kasp_criteria() afterwards.)
LAMP_GATES <- c(
  "gate_lamp_variant", "gate_lamp_inner_tm", "gate_lamp_outer_tm",
  "gate_lamp_structure", "gate_lamp_amplicon"
)
LAMP_RELAX <- c("gate_lamp_amplicon", "gate_lamp_outer_tm")
LAMP_TIERS <- c(
  tier1 = "loop-region specificity", tier2 = "Tm set", tier3 = "amplicon length"
)

#' The LAMP scheme definition (registry row + decision stage)
#'
#' Self-contained description of the LAMP chemistry in the shape the shared
#' registry expects ([kasp_schemes()] fields plus a decision stage: `gates`,
#' `relax`, `tiers` and their human-readable labels).  Declared next to the
#' designer so a scheme can be added without touching the registry first; the
#' coordinator wires it in afterwards.
#' @return a list (see the `scheme_def_*` naming convention used for the
#'   registry's schemes).
#' @keywords internal
#' @noRd
scheme_def_lamp <- function() {
  list(
    scheme = "lamp",
    family = "lamp",
    label = "LAMP (loop-mediated isothermal amplification)",
    discrimination = "the variant sits in the F2/B2 loop, at the FIP/BIP ends",
    readout = "turbidity / fluorescence / colorimetric / gel (isothermal)",
    zygosity = "dominant (present/absent), co-dominant with probes",
    amplicon = "120-300 bp",
    loci = 1L,
    verdict = "usable LAMP marker",
    gates = LAMP_GATES,
    relax = LAMP_RELAX,
    tiers = LAMP_TIERS,
    gate_short = c(
      gate_lamp_variant = "variant in loop", gate_lamp_inner_tm = "FIP/BIP Tm",
      gate_lamp_outer_tm = "F3/B3 Tm", gate_lamp_structure = "primer structure",
      gate_lamp_amplicon = "amplicon length (LAMP)"
    ),
    gate_limit = c(
      gate_lamp_variant = "variant between F2 and B2 (in the loop)",
      gate_lamp_inner_tm = "FIP/BIP binding arm Tm within [60, 65] C",
      gate_lamp_outer_tm = "F3/B3 Tm within [55, 62] C",
      gate_lamp_structure = "no inverted repeat in any primer",
      gate_lamp_amplicon = "amplicon within [120, 300] bp"
    ),
    gate_floor_label = c(
      gate_lamp_amplicon = "widen the amplicon window",
      gate_lamp_outer_tm = "widen the outer Tm window"
    ),
    gate_hard_label = c(
      gate_lamp_variant = "a variant outside the loop cannot be discriminated",
      gate_lamp_inner_tm = "the inner primers must anneal isothermally",
      gate_lamp_structure = "primer hairpins / dimers defeat LAMP"
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

#' Build one inner LAMP primer (FIP or BIP) as a two-arm fusion
#'
#' The inner primers are fusions: `FIP = F1c + F2` and `BIP = B1c + B2`, i.e. the
#' 3' "displacement" arm (`F2`/`B2`) is the reverse complement of a plus-strand
#' site, and the 5' "loop" arm (`F1c`/`B1c`) is the reverse complement of a plus
#' strand site lying between the displacement arm and the variant, so the variant
#' ends up inside the loop between the proximal arms.  Both arms must sit in the
#' inner Tm/GC window; the fusion is screened for inverted repeats.
#' @param seq the template (plus strand).
#' @param q 1-based variant position on `seq`.
#' @param side `"up"` (variant-upstream, build FIP) or `"down"` (downstream, BIP).
#' @param lengths arm lengths to try.
#' @param tm_range,gc_range accepted arm windows.
#' @param target_tm binding-arm Tm to aim for when several fusions pass.
#' @return a list (`primer`, `tm`, `gap`, `prox`, `dist`), or `NULL`.
#' @keywords internal
#' @noRd
lamp_fusion <- function(seq, q, side, lengths, tm_range, gc_range, target_tm) {
  n <- nchar(as.character(seq))
  starts <- if (identical(side, "up")) seq_len(q - 1L) else seq.int(q + 1L, n)
  sites <- primer_candidates(seq, starts, lengths, "+", tm_range, gc_range) |>
    dplyr::filter(!has_palindrome(.data$primer))
  if (!nrow(sites)) {
    return(NULL)
  }
  prox <- dplyr::rename(sites, s1 = "start", e1 = "end", p1 = "primer", t1 = "tm")
  dist <- dplyr::rename(sites, s2 = "start", e2 = "end", p2 = "primer", t2 = "tm")
  pairs <- dplyr::cross_join(prox, dist)
  pairs <- if (identical(side, "up")) {
    # distal arm upstream of the proximal arm; the variant is past the proximal end
    dplyr::filter(pairs, .data$e2 < .data$s1)
  } else {
    dplyr::filter(pairs, .data$e1 < .data$s2)
  }
  pairs <- dplyr::mutate(
    pairs,
    primer = paste0(lamp_revcomp(.data$p1), lamp_revcomp(.data$p2)),
    gap = if (identical(side, "up")) q - .data$e1 else .data$s1 - q,
    tm_bind = pmax(.data$t1, .data$t2)
  ) |>
    dplyr::filter(!has_palindrome(.data$primer))
  if (!nrow(pairs)) {
    return(NULL)
  }
  best <- pairs |>
    dplyr::arrange(.data$gap, abs(.data$tm_bind - target_tm)) |>
    dplyr::slice_head(n = 1L)
  list(
    primer = best$primer, tm = best$tm_bind, gap = best$gap,
    prox = c(best$s1, best$e1), dist = c(best$s2, best$e2)
  )
}

#' Design a LAMP marker at a variant
#'
#' The isothermal analogue of [assemble_kasp_assay()]: the outer F3/B3 pair and
#' the inner FIP/BIP pair are placed so the variant sits in the loop between F2
#' and B2 - at the 5' ends of the inner primers, where the allelic difference
#' biases loop formation.  FIP and BIP are two-arm fusions (`F1c+F2`, `B1c+B2`)
#' built by [lamp_fusion()]; F3/B3 reuse the generic [design_flanking_pair()] over
#' the inner block.  The designer does **not** simulate the isothermal reaction
#' (that needs a polymerase/kinetics model), so it returns a well-formed primer
#' set, not a predicted amplification curve.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param arm_lengths outer / inner arm lengths to try.
#' @param outer_tm accepted F3/B3 Tm window (C).
#' @param inner_tm accepted inner binding-arm Tm window (C).
#' @param gc_range accepted GC fraction window (both arm classes).
#' @param min_product,max_product accepted amplicon window (bp).
#' @return a `kasp_lamp` list (`snp`, `primers`, `amplicon`, `loop`, `tm`,
#'   `compatibility` - the joint dimer/hairpin screen over all four primers), or
#'   `NULL` (with a message) when no primer set passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_lamp_marker(genome, kasp_example("SNP")[1, ])
#' }
#' @family KASP LAMP
#' @export
design_lamp_marker <- function(genome, snp, flank = 150L, arm_lengths = 18:22,
                               outer_tm = c(55, 62), inner_tm = c(60, 65),
                               gc_range = KASP_LIMITS$gc,
                               min_product = 120L, max_product = 300L) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  q <- tmpl$snp_i
  target_tm <- mean(inner_tm)

  fip <- lamp_fusion(tmpl$seq, q, "up", arm_lengths, inner_tm, gc_range, target_tm)
  bip <- lamp_fusion(tmpl$seq, q, "down", arm_lengths, inner_tm, gc_range, target_tm)
  if (is.null(fip) || is.null(bip)) {
    cli::cli_inform(gettext("No FIP/BIP inner primer pair passed the loop/Tm/structure filters for LAMP."))
    return(NULL)
  }

  # the outer pair flanks the whole inner block, so F3/B3 never overlap FIP/BIP;
  # unlike design_flanking_pair() the candidates are also screened for hairpins
  lo <- fip$dist[1L]
  hi <- bip$dist[2L]
  n <- nchar(as.character(tmpl$seq))
  fwd <- primer_candidates(tmpl$seq, seq_len(lo - 1L), arm_lengths, "+", outer_tm, gc_range) |>
    dplyr::filter(.data$end < lo, !has_palindrome(.data$primer))
  rev <- primer_candidates(tmpl$seq, seq.int(hi + 1L, n), arm_lengths, "-", outer_tm, gc_range) |>
    dplyr::filter(!has_palindrome(.data$primer))
  pair <- pick_primer_pair(fwd, rev, min_product, max_product)
  if (is.null(pair)) {
    cli::cli_inform(gettext("No F3/B3 outer pair passed the Tm/GC, structure and amplicon filters for LAMP."))
    return(NULL)
  }
  outer <- pair_to_list(pair)

  off <- snp$pos - q
  loop_start <- fip$dist[2L] + 1L
  loop_end <- bip$dist[1L] - 1L
  if (!(loop_start <= q && q <= loop_end)) {
    cli::cli_inform(gettext("The variant does not fall between F2 and B2 (LAMP loop region)."))
    return(NULL)
  }
  out <- list(
    snp = snp,
    primers = list(F3 = outer$forward, B3 = outer$reverse, FIP = fip$primer, BIP = bip$primer),
    amplicon = tibble::tibble(start = outer$start + off, end = outer$end + off),
    loop = tibble::tibble(start = loop_start + off, end = loop_end + off),
    tm = c(F3 = outer$tm_f, B3 = outer$tm_r, FIP = fip$tm, BIP = bip$tm)
  )
  # the four primers share one tube, so screen the *joint* set (cross-dimers and
  # hairpins between them), not each primer alone
  out$compatibility <- tube_compatibility(out$primers)
  structure(out, class = "kasp_lamp")
}

#' @export
print.kasp_lamp <- function(x, ...) {
  cat(
    "<kasp_lamp>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  amplicon:", x$amplicon$end - x$amplicon$start + 1L, "bp\n")
  cat("  variant loop:", x$loop$start, "-", x$loop$end, "\n")
  cat("  F3:", x$primers$F3, "\n  B3:", x$primers$B3, "\n")
  cat("  FIP:", x$primers$FIP, "\n  BIP:", x$primers$BIP, "\n")
  invisible(x)
}
