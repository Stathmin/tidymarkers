# =============================================================================
#  marker_crispr.R  --  CRISPR-Dx (SHERLOCK / DETECTR) as a marker scheme
# =============================================================================
#  A CRISPR-Dx marker does not read the allele from a primer 3' base (KASP) or a
#  band size (CAPS/ARMS): an isothermal pre-amplification (RPA) of the
#  region around the variant feeds a Cas effector whose guide RNA (crRNA) places
#  the variant inside its spacer, and target binding licenses *collateral*
#  cleavage of a reporter - the readout is fluorescence (or a lateral-flow band),
#  not a fragment length.  So the design levers are the guide, not the primers:
#  the spacer must cover the variant, the variant must sit in the seed / critical
#  region, and - for Cas12a - a PAM must license the protospacer.
#
#  Design rules (researched, not invented):
#   * Cas12a (DETECTR) requires a 5'-TTTV PAM (V = A/C/G) immediately 5' of the
#     protospacer, and the PAM-proximal seed is the mismatch-critical region; the
#     spacer is ~20-23 nt.  Cas12a cuts its DNA target in cis and then acts as an
#     indiscriminate ssDNase in trans ("collateral" cleavage); the DETECTR readout
#     is a quenched ssDNA reporter, and the collateral activity does not degrade
#     a double-stranded reporter.
#     Chen JS et al. (2018) Science 360:436-439, doi:10.1126/science.aar6245.
#     Addgene CRISPR Guide, PAM table (AsCpf1/LbCpf1: 5' TTTV; seed at the
#     PAM-proximal end): https://www.addgene.org/guides/crispr/
#   * Cas13 (SHERLOCK) targets ssRNA, so it has no DNA PAM - but it is NOT simply
#     "PAM-free": a protospacer-flanking site (PFS) constrains the target.  Cas13a
#     has a weak 3' PFS preference for H (A/C/U), so a G immediately 3' of the
#     protospacer is avoided; that conservative rule is enforced here even though
#     LwaCas13a tolerates some 3' G in vitro.  Cas13b requires a 5' D-PFS
#     (D = A/U/G, i.e. not C); its additional 3' NAN/NNA PFS is documented but not
#     enforced (the designer warns for `cas13b`).  The LwaCas13a crRNA spacer is
#     ~28 nt, so the Cas13a window is 26-28 nt; the Cas13b window is 28-30 nt
#     (canonical ~30 nt).  Cas13 collateral cleavage is an ssRNase activity read by
#     an ssRNA reporter; an isothermal pre-amplification precedes T7 transcription.
#     Abudayyeh OO et al. (2016) Science 353:aaf5573, doi:10.1126/science.aaf5573;
#     Gootenberg JS et al. (2017) Science 356:438-442, doi:10.1126/science.aam9321;
#     Gootenberg JS et al. (2018) Science 360:439-444, doi:10.1126/science.aaq0179
#     (SHERLOCKv2); Smargon AA et al. (2017) Mol Cell 65:618-630,
#     doi:10.1016/j.molcel.2016.12.023 (Cas13b 5' D-PFS and 3' PFS).
#   * Cas13 cleaves ssRNA, so its guide must match the *transcribed* strand (the
#     one the T7 promoter turns into the target RNA), not the genomic plus strand
#     automatically; guides are searched on both strands unless `strand=` restricts
#     them, and the caller validates the choice against the T7 orientation.
#   * The pre-amplification literature works best on short amplicons; we accept a
#     bounded window (80-250 bp) rather than a single value.
#
#  Assumptions stated where the literature does not fix a number: the guide GC
#  window (0.35-0.65, around the common 40-60% ideal), the Cas13 "critical"
#  region (the central half of the spacer, the least mismatch-tolerant region),
#  and the pre-amplification pair's own windows (RPA, never a KASP PCR).
#
#  Discrimination, spelled out: the guide is matched to ONE allele and mismatched
#  to the other, so the variant must sit where a single mismatch abolishes
#  activation - the PAM-proximal seed for Cas12a, the central critical region for
#  Cas13 - and the assay must not be licensed by an off-target protospacer
#  elsewhere in the genome (a second full-match spacer would fire on the wrong
#  template).  The pre-amplification is allele-blind (it copies both alleles), so
#  it is graded only on its own chemistry's composition - never on KASP Tm/GC.
# =============================================================================

# Cas12a's PAM-proximal seed is where the variant must sit for the guide to read
# it as allele-specific; for Cas13 the critical region is central (no PAM, no
# PAM-proximal seed).  `crispr_seed_evidence()` carries what was measured about
# the magnitude, which this constant does not encode.
CRISPR_SEED_LEN <- 8L

#' What the literature measured about a seed mismatch (evidence, not a gate)
#'
#' `crispr_seed_gate()` is a POSITION rule: it says where a variant must sit, not
#' that a seed mismatch abolishes activity.  The one source that measures the
#' magnitude (Commun Biol 2024, 40 activators x 8 crRNAs, endpoint trans-cleavage
#' fluorescence, mutant/WT ratio; 1.0 = no discrimination) says:
#' * an UNMODIFIED crRNA does not discriminate at all: every single mismatch over
#'   positions 1-20 stays in 0.767-1.394;
#' * on the optimised crRNA `M39` a single seed mismatch (positions 1-8) still
#'   transmits 0.373-0.968 of WT signal - reduced, not abolished;
#' * two adjacent mismatches are usually stronger (17 of 19 below 0.5 vs 8 of 20
#'   singles) but the ordering FAILS at the spacer ends, so it is a trend and not
#'   a filter.
#' Since the package cannot predict the ratio for a given target, the honest
#' report is `unknown` on the magnitude; the gate never passes because a position
#' is merely inside the seed.  See
#' `data-raw/proposed/crispr_seed_rule.decision.md`.
#' @return a classed list of the measured evidence.
#' @keywords internal
#' @noRd
crispr_seed_evidence <- function() {
  structure(
    list(
      gate = "position only: the variant must sit in the PAM-proximal seed",
      ratio_range = c(min = 0.3728637, max = 0.9679837),
      unmodified_crrna_range = c(0.766701, 1.393679),
      n_seed_singles = 20L,
      adjacent_better = "17 of 19 adjacent doubles below 0.5 vs 8 of 20 singles; fails at the last two positions",
      magnitude = "unknown - no model predicts the ratio for a target",
      source = "10.1038/s42003-024-07173-7",
      calibrated = FALSE
    ),
    class = "marker_crispr_seed_evidence"
  )
}

#' The seed MISSES what it was described as measuring
#'
#' @param x a [crispr_seed_evidence()] result.
#' @param ... unused.
#' @return `x`, invisibly.
#' @keywords internal
#' @noRd
print.marker_crispr_seed_evidence <- function(x, ...) {
  cli::cli_inform(c(
    "i" = x$gate,
    "*" = "single seed mismatch keeps {x$ratio_range[[1L]]}-{x$ratio_range[[2L]]} of WT signal (measured)",
    "*" = "magnitude: {x$magnitude}"
  ))
  invisible(x)
}

#' Does the variant sit in the guide's discriminating region?
#'
#' A POSITION rule, not a magnitude one.  Cas12a needs the variant in the
#' PAM-proximal seed; Cas13 needs it in the central critical region.  The
#' candidate grid enforces this while it is built, so the gate re-checks the
#' chosen guide (a constructor invariant) rather than returning an unconditional
#' `TRUE`.  Being inside the seed does not mean a single mismatch abolishes
#' activity - see [crispr_seed_evidence()] for what was measured.
#' @param cas the effector.
#' @param len spacer length.
#' @param v_off the variant's offset in the spacer (0-based, PAM-proximal end).
#' @return a single logical.
#' @keywords internal
#' @noRd
crispr_seed_gate <- function(cas, len, v_off) {
  if (length(len) != 1L || length(v_off) != 1L || is.na(len) || is.na(v_off)) {
    return(FALSE)
  }
  if (identical(cas, "cas12a")) {
    as.integer(v_off) < CRISPR_SEED_LEN
  } else {
    abs(v_off - (len - 1) / 2) <= len %/% 4L
  }
}

#' Candidate guide spacers covering a variant
#'
#' Enumerates guide spacers of the accepted lengths that place the variant at
#' each offset of the *critical* window, on both strands.  For `cas = "cas12a"` a
#' 5'-TTTV PAM immediately adjacent to the protospacer is required and the
#' critical window is the PAM-proximal seed.  Cas13 has no DNA PAM (it targets
#' ssRNA) but it is not PFS-free: the target must carry a protospacer-flanking
#' site (PFS), which is recorded in `pfs` and filtered here.  For Cas13a (and its
#' `"cas13"` alias) the PFS is the base immediately 3' of the target protospacer
#' and G is rejected (H = A/C/U preference).  For Cas13b the PFS is the 5' D-PFS
#' (D = A/U/G) and C is rejected; the additional 3' NAN/NNA PFS is not enforced.
#' The Cas13 critical window is the central region of the spacer.  Candidates are
#' ordered best-first (deepest in the seed / most central, then GC closest to
#' 0.5).  PFS bases are returned in the DNA alphabet (`T` stands for `U` in the
#' target RNA) and are read from the strand the guide spacer matches.
#' @param seq the plus-strand template.
#' @param snp_i template coordinate of the variant.
#' @param cas `"cas12a"`, `"cas13"`/`"cas13a"`, or `"cas13b"`.
#' @param spacer_lengths spacer lengths to try.
#' @param gc_range accepted guide GC window.
#' @param strands strands the guide may be placed on (`"+"`/`"-"`).
#' @return a tibble (`start`, `end`, `len`, `strand`, `spacer`, `pam`, `pfs`,
#'   `pfs_side`, `gc`, `v_off`), plus-strand coordinates, or zero rows.  `pam` is
#'   `NA` for Cas13; `pfs`/`pfs_side` are `NA` for Cas12a.
#' @keywords internal
#' @noRd
crispr_guide_candidates <- function(seq, snp_i, cas, spacer_lengths, gc_range, strands = c("+", "-")) {
  strands <- match.arg(strands, c("+", "-"), several.ok = TRUE)
  cas <- match.arg(cas, c("cas12a", "cas13", "cas13a", "cas13b"))
  s <- as.character(seq)
  n <- nchar(s)
  seed_len <- CRISPR_SEED_LEN
  grid <- tidyr::crossing(len = spacer_lengths, v_off = 0:(max(spacer_lengths) - 1L)) |>
    dplyr::filter(.data$v_off < .data$len)
  if (identical(cas, "cas12a")) {
    grid <- dplyr::filter(grid, .data$v_off < seed_len)
  } else {
    grid <- dplyr::filter(grid, abs(.data$v_off - (.data$len - 1) / 2) <= .data$len %/% 4L)
  }
  plus <- dplyr::mutate(grid, strand = "+", start = snp_i - .data$v_off, end = .data$start + .data$len - 1L)
  minus <- dplyr::mutate(grid, strand = "-", end = snp_i + .data$v_off, start = .data$end - .data$len + 1L)
  cand <- dplyr::bind_rows(plus, minus) |>
    dplyr::filter(.data$strand %in% strands, .data$start >= 1L, .data$end <= n) |>
    dplyr::mutate(
      pam_lo = dplyr::if_else(.data$strand == "+", .data$start - 4L, .data$end + 1L),
      pam_hi = dplyr::if_else(.data$strand == "+", .data$start - 1L, .data$end + 4L)
    ) |>
    dplyr::filter(.data$pam_lo >= 1L, .data$pam_hi <= n)
  if (!nrow(cand)) {
    return(cand)
  }
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  pfs_side <- if (cas %in% c("cas13", "cas13a")) {
    "3'"
  } else if (identical(cas, "cas13b")) {
    "5'"
  } else {
    NA_character_
  }
  pfs_side_val <- pfs_side
  cand <- cand |>
    dplyr::mutate(
      plus_spacer = stringr::str_sub(s, .data$start, .data$end),
      pam_raw = stringr::str_sub(s, .data$pam_lo, .data$pam_hi),
      plus_pfs = if (identical(pfs_side, "3'")) {
        stringr::str_sub(s, .data$end + 1L, .data$end + 1L)
      } else if (identical(pfs_side, "5'")) {
        stringr::str_sub(s, .data$start - 1L, .data$start - 1L)
      } else {
        rep(NA_character_, dplyr::n())
      },
      minus_pfs = if (identical(pfs_side, "3'")) {
        rc(stringr::str_sub(s, .data$start - 1L, .data$start - 1L))
      } else if (identical(pfs_side, "5'")) {
        rc(stringr::str_sub(s, .data$end + 1L, .data$end + 1L))
      } else {
        rep(NA_character_, dplyr::n())
      }
    )
  cand <- cand |>
    dplyr::mutate(
      spacer = dplyr::if_else(.data$strand == "+", .data$plus_spacer, rc(.data$plus_spacer)),
      pam = dplyr::if_else(.data$strand == "+", .data$pam_raw, rc(.data$pam_raw)),
      pfs = dplyr::if_else(.data$strand == "+", .data$plus_pfs, .data$minus_pfs),
      pfs_side = pfs_side_val,
      gc = gc_frac(.data$spacer)
    )
  if (identical(cas, "cas12a")) {
    cand <- dplyr::filter(cand, stringr::str_detect(.data$pam, "^TTT[ACG]$"))
  } else if (identical(pfs_side, "3'")) {
    cand <- dplyr::filter(cand, .data$pfs != "G")
  } else {
    cand <- dplyr::filter(cand, .data$pfs != "C")
  }
  cand <- cand |>
    dplyr::filter(dplyr::between(.data$gc, gc_range[1L], gc_range[2L])) |>
    dplyr::mutate(gc_dev = abs(.data$gc - 0.5))
  if (!identical(cas, "cas12a")) {
    cand <- dplyr::mutate(cand, pam = NA_character_)
  }
  if (identical(cas, "cas12a")) {
    cand <- dplyr::arrange(cand, .data$v_off, .data$gc_dev)
  } else {
    cand <- dplyr::mutate(cand, pos_dev = abs(.data$v_off - (.data$len - 1) / 2)) |>
      dplyr::arrange(.data$pos_dev, .data$gc_dev)
  }
  dplyr::select(
    cand, "start", "end", "len", "strand", "spacer", "pam", "pfs",
    "pfs_side", "gc", "v_off"
  )
}

#' Genome-wide exact hits of a guide spacer (both strands)
#'
#' A guide is specific only when its protospacer occurs once - the on-target
#' locus.  Any other exact hit is a full-match off-target a Cas effector could be
#' licensed by, so it is counted.  Exact matching (no seed mismatch) is used on
#' purpose: a single-mismatch count in a 100+ Mb genome is dominated by chance and
#' would not mean a real off-target protospacer.
#' @param spacer the guide spacer (5'->3').
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param BPPARAM optional `BiocParallel` backend (one task per chromosome).
#' @return the both-strand hit count (integer), or `NA` when uncomputable.
#' @keywords internal
#' @noRd
crispr_spacer_hits <- function(spacer, genome, BPPARAM = NULL) {
  if (is.null(genome) || length(spacer) != 1L || is.na(spacer) || !nzchar(spacer)) {
    return(NA_integer_)
  }
  fwd <- Biostrings::DNAString(spacer)
  rev <- Biostrings::reverseComplement(fwd)
  count_chr <- function(nm) {
    chr <- genome_seq(genome, nm)
    as.integer(Biostrings::countPattern(fwd, chr) + Biostrings::countPattern(rev, chr))
  }
  chroms <- genome_chroms(genome)
  hits <- if (is.null(BPPARAM)) {
    lapply(chroms, count_chr)
  } else {
    bp_lapply(chroms, count_chr, BPPARAM)
  }
  as.integer(sum(unlist(hits)))
}

#' Default guide spacer windows per Cas effector
#' @param cas one of `"cas12a"`, `"cas13"`, `"cas13a"`, `"cas13b"`.
#' @return an integer vector of spacer lengths: 20-23 nt (Cas12a), 26-28 nt
#'   (Cas13a/`"cas13"`, around the ~28 nt LwaCas13a spacer), 28-30 nt (Cas13b).
#' @keywords internal
#' @noRd
crispr_default_spacer <- function(cas) {
  switch(cas,
    cas12a = 20:23,
    cas13 = 26:28,
    cas13a = 26:28,
    cas13b = 28:30,
    stop("unknown Cas effector: ", cas, call. = FALSE)
  )
}

#' Design a CRISPR-Dx (SHERLOCK / DETECTR) marker at a variant
#'
#' The Cas-mediated detection analogue of [assemble_kasp_assay()]: an isothermal
#' RPA pre-amplification pair straddles the variant and a guide RNA is
#' placed so its spacer covers the variant, with the variant in the seed or
#' critical region.  For Cas12a a 5'-TTTV PAM must license the protospacer.
#' Cas13 has no DNA PAM (it targets ssRNA), but it is not PFS-free: Cas13a needs
#' a non-G base immediately 3' of the target protospacer, and Cas13b needs a
#' 5' D-PFS (A/U/G, i.e. not C).  The guide must also be unique in the genome
#' (no off-target protospacer).  The collateral-cleavage readout itself is not
#' simulated, so this returns a well-formed guide plus its pre-amplification
#' amplicon, not a predicted signal.
#'
#' The guide is chosen best-first among all candidates whose pre-amplification
#' pair also passes *that chemistry's* composition filters (RPA, never a KASP
#' PCR), so a guide is only returned when its amplicon is designable too.
#' The soft gates (amplicon window, GC, spacer length) are relaxed one at a
#' time, least grave first, only when no guide clears them all; the PAM/PFS,
#' the seed/critical region and guide specificity (no off-target protospacer)
#' are hard - no guide is returned past them.  For `cas13b` the 3' NAN/NNA PFS
#' is not enforced, and the designer warns about that gap.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param cas `"cas12a"` (DETECTR; needs a 5'-TTTV PAM), `"cas13"` or
#'   `"cas13a"` (SHERLOCK; no DNA PAM, 3' non-G PFS), or `"cas13b"`
#'   (5' D-PFS; the 3' PFS is not enforced).
#' @param strand guide placement: `"both"` (default), or `"+"`/`"-"`.  Cas13
#'   cleaves ssRNA, so its guide must match the *transcribed* strand (the one the
#'   pre-amplification's T7 promoter turns into the target RNA); pass `strand`
#'   when that orientation is known, otherwise both are tried and the choice is
#'   the caller's to validate.
#' @param spacer_lengths guide spacer lengths to try (default: 20-23 nt for
#'   Cas12a, 26-28 nt for Cas13a/`"cas13"` around the ~28 nt LwaCas13a spacer,
#'   and 28-30 nt for Cas13b).
#' @param gc_range accepted guide GC window (default 0.35-0.65).
#' @param specificity_max accepted number of off-target protospacers (exact
#'   genome hits of the spacer beyond the on-target locus; default 0).  This is
#'   a hard gate: it is never relaxed to `specificity_max + 1`.
#' @param preamp the isothermal pre-amplification the guide rides on, `"rpa"`
#'   (the only implemented chemistry); it sets that pair's length / Tm / GC
#'   windows.
#' @param min_product,max_product accepted pre-amplification amplicon window.
#' @param haplotypes optional population haplotype table.  When supplied, it is
#'   attached for guide/PAM conservation warnings and the generic haplotype
#'   track; when absent, CRISPR warns that guide/PAM conservation is unchecked.
#' @param BPPARAM optional `BiocParallel` backend for the guide-specificity scan.
#' @return a `marker_crispr` list (`snp`, `cas`, `strand`, `guide`, `pam`, `pfs`,
#'   `pfs_side`, `guide_start`, `guide_end`, `seed_offset`, `gc`, `amplicon`,
#'   `off_target`, `gates`, `relaxed`), or `NULL` (with a message) when no
#'   guide / amplicon passes the hard gates.  `pam` is `NA` for Cas13;
#'   `pfs`/`pfs_side` are `NA` for Cas12a.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_crispr_marker(genome, marker_example("SNP")[1, ], cas = "cas12a")
#' design_crispr_marker(genome, marker_example("SNP")[1, ], cas = "cas13")
#' }
#' @family KASP CRISPR
#' @export
design_crispr_marker <- function(genome, snp, flank = 200L,
                                 cas = c("cas12a", "cas13", "cas13a", "cas13b"),
                                 strand = c("both", "+", "-"),
                                 spacer_lengths = NULL, gc_range = NULL,
                                 specificity_max = 0L,
                                 preamp = "rpa",
                                 min_product = 80L, max_product = 250L,
                                 haplotypes = NULL, BPPARAM = NULL) {
  cas <- match.arg(cas)
  strand <- match.arg(strand)
  preamp <- match.arg(preamp, "rpa")
  if (identical(cas, "cas13b")) {
    cli::cli_warn(gettext(
      "Cas13b: only the 5' D-PFS is enforced; the 3' PFS (NAN/NNA) is not modelled."
    ))
  }
  strands <- if (identical(strand, "both")) c("+", "-") else strand
  strict_spacer <- if (is.null(spacer_lengths)) {
    crispr_default_spacer(cas)
  } else {
    as.integer(spacer_lengths)
  }
  strict_gc <- gc_range %||% c(0.35, 0.65)
  # the relaxation floors (bounded; C8): a guide a little outside the ideal
  # window still works, so these gates are soft, but never widened without bound.
  relax_spacer <- seq.int(min(strict_spacer) - 2L, max(strict_spacer) + 2L)
  relax_gc <- c(max(0, strict_gc[1L] - 0.10), min(1, strict_gc[2L] + 0.10))
  # the pre-amplification is RPA, so its pair is screened with THAT chemistry's
  # length / Tm / GC windows - never KASP's PCR conditions.
  pre <- list(lengths = RPA_PRIMER_LEN, tm = RPA_NO_TM, gc = RPA_GC)
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  off <- snp$pos - tmpl$snp_i
  # the ladder to the first usable guide: relax the soft gates LEAST grave first
  # (C6), one rung at a time.  PAM/PFS, seed and specificity are hard: a guide
  # without a PAM (Cas12a) or a compatible PFS (Cas13), with the variant outside
  # the seed/critical region, or with an off-target protospacer cannot be
  # licensed to discriminate, so they are enforced by the enumeration / hit
  # count and never relaxed.
  rungs <- list(
    character(0),
    "gate_crispr_amplicon",
    c("gate_crispr_amplicon", "gate_crispr_gc"),
    c("gate_crispr_amplicon", "gate_crispr_gc", "gate_crispr_spacer")
  )
  attempt <- function(relaxed) {
    sl <- if ("gate_crispr_spacer" %in% relaxed) relax_spacer else strict_spacer
    gl <- if ("gate_crispr_gc" %in% relaxed) relax_gc else strict_gc
    spec_max <- specificity_max # hard: never widened, even at the last rung
    amp_lo <- if ("gate_crispr_amplicon" %in% relaxed) min_product - 20L else min_product
    amp_hi <- if ("gate_crispr_amplicon" %in% relaxed) max_product + 20L else max_product
    cand <- crispr_guide_candidates(tmpl$seq, tmpl$snp_i, cas, sl, gl, strands)
    if (!nrow(cand)) {
      return(NULL)
    }
    for (i in seq_len(nrow(cand))) {
      hits <- crispr_spacer_hits(cand$spacer[i], genome, BPPARAM)
      if (!is.na(hits) && (hits - 1L) > spec_max) next # an off-target protospacer
      amp <- design_flanking_pair(tmpl$seq, cand$start[i], cand$end[i],
        lengths = pre$lengths, min_product = amp_lo, max_product = amp_hi,
        tm_range = pre$tm, gc_range = pre$gc
      )
      if (!is.null(amp)) {
        return(list(row = cand[i, ], amp = amp, hits = hits))
      }
    }
    NULL
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
      "No guide spacer passed the PAM/PFS / seed / specificity / GC filters with a pre-amplification pair."
    ))
    return(NULL)
  }
  row <- pick$row
  amp <- pick$amp
  off_target <- if (is.na(pick$hits)) NA_integer_ else pick$hits - 1L
  # the gate table the DAG displays, evaluated on the guide actually returned
  # (C7).  `gate_crispr_pam` is the legacy key for the licensing gate: the
  # Cas12a PAM or, for Cas13, the PFS (no DNA PAM exists for an RNA target).
  gates <- c(
    gate_crispr_pam = if (identical(cas, "cas12a")) {
      grepl("^TTT[ACG]$", row$pam)
    } else if (identical(cas, "cas13b")) {
      !is.na(row$pfs) && row$pfs != "C"
    } else {
      !is.na(row$pfs) && row$pfs != "G"
    },
    # constructor invariant: the candidate grid put the variant in the seed /
    # central window; re-check the chosen guide rather than assuming it
    gate_crispr_seed = crispr_seed_gate(cas, row$len, row$v_off),
    gate_crispr_spacer = dplyr::between(row$len, min(strict_spacer), max(strict_spacer)),
    gate_crispr_specificity = if (is.na(off_target)) NA else off_target <= specificity_max,
    gate_crispr_gc = dplyr::between(row$gc, strict_gc[1L], strict_gc[2L]),
    gate_crispr_amplicon = dplyr::between(amp$product, min_product, max_product)
  )
  soft <- c(
    "gate_crispr_spacer", "gate_crispr_gc", "gate_crispr_amplicon"
  )
  relaxed <- soft[vapply(soft, function(g) isFALSE(gates[[g]]), logical(1))]
  out <- list(
    snp = snp,
    cas = cas,
    strand = row$strand,
    guide = row$spacer,
    pam = if (identical(cas, "cas12a")) row$pam else NA_character_,
    pfs = row$pfs,
    pfs_side = row$pfs_side,
    guide_start = row$start + off,
    guide_end = row$end + off,
    seed_offset = row$v_off,
    gc = row$gc,
    preamp = preamp,
    off_target = off_target,
    gates = gates,
    relaxed = relaxed,
    amplicon = list(
      forward = amp$forward, reverse = amp$reverse,
      start = amp$start + off, end = amp$end + off, product = amp$product
    )
  )
  out <- structure(out, class = "marker_crispr")
  if (is.null(haplotypes)) {
    cli::cli_warn(gettext(
      "CRISPR: no population haplotypes supplied; guide/PAM conservation is not checked."
    ))
    return(out)
  }
  marker_attach_haplotypes(out, haplotypes, snp)
}

#' @export
print.marker_crispr <- function(x, ...) {
  cat(
    "<marker_crispr>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, " [", x$cas, "]\n",
    sep = ""
  )
  pfs <- x$pfs %||% NA_character_
  lic <- if (!is.na(x$pam)) {
    paste0("PAM ", x$pam)
  } else if (identical(x$cas, "cas13b")) {
    paste0("no PAM; 5' PFS ", pfs, " (D)")
  } else {
    paste0("no PAM; 3' PFS ", pfs, " (non-G)")
  }
  cat("  guide:", x$guide, "(", x$strand, ",", lic, ", seed offset", x$seed_offset, ")\n")
  ot <- if (is.na(x$off_target)) "unscanned" else paste0(x$off_target, " off-target protospacer(s)")
  cat("  specificity:", ot, "\n")
  cat("  pre-amplification (", x$preamp, ") amplicon:", x$amplicon$product, "bp\n", sep = "")
  cat("  F:", x$amplicon$forward, "\n  R:", x$amplicon$reverse, "\n")
  if (length(x$relaxed)) {
    cat("  (relaxed:", paste(x$relaxed, collapse = ", "), ")\n")
  }
  invisible(x)
}

#' The CRISPR-Dx scheme definition (internal)
#'
#' The registry entry the shared decision schema consumes, one marker chemistry
#' per call.  It is discovered automatically (one file per scheme, see
#' `R/marker_scheme_registry.R`), so it is kept internal.
#' @return a list with the scheme identity, readout, gates, relax set, tiers and
#'   the label / limit lookups the decision graph and text export need.
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_crispr
#' @export
#' @noRd
marker_map_spec.marker_crispr <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  guide <- list(spec_row(
    "guide", marker$guide_start[1L], marker$guide_end[1L], marker$strand[1L], "guide",
    label = "guide"
  ))
  pam <- if (!is.na(marker$pam[1L])) {
    list(spec_row(
      "PAM", marker$guide_start[1L] - 4L, marker$guide_start[1L] - 1L, "+", "site",
      label = "PAM"
    ))
  } else {
    list()
  }
  map_spec(
    marker,
    features = c(map_pair_rows(marker, sp[["start"]], sp[["end"]]), guide, pam),
    sequences = list(
      forward = marker$primers$forward, reverse = marker$primers$reverse,
      guide = marker$guide, PAM = marker$pam
    )
  )
}

scheme_def_crispr <- function() {
  list(
    scheme = "crispr",
    family = "crispr",
    label = "CRISPR-Dx (SHERLOCK / DETECTR)",
    discrimination = paste0(
      "Cas collateral cleavage of a reporter, guided to the variant ",
      "(ssDNase/ssDNA reporter for Cas12a; ssRNase/ssRNA for Cas13)"
    ),
    readout = "fluorescence / lateral flow (ssDNA reporter for Cas12a; ssRNA for Cas13)",
    zygosity = "co-dominant (dominant for presence/absence)",
    amplicon = "pre-amplified (RPA)",
    loci = 1L,
    verdict = "usable CRISPR-Dx marker",
    # the gate order the method applies (grave / cheap first): a guide is
    # licensed by the PAM (Cas12a) or PFS (Cas13) and the seed first, then its
    # length, specificity, GC and finally the pre-amplification amplicon it rides
    # on.  `gate_crispr_pam` is the legacy key for that licensing gate; for Cas13
    # it reports the PFS, because an RNA target has no DNA PAM.
    gates = c(
      "gate_crispr_pam", "gate_crispr_seed", "gate_crispr_spacer",
      "gate_crispr_specificity", "gate_crispr_gc", "gate_crispr_amplicon"
    ),
    relax = c(
      "gate_crispr_amplicon", "gate_crispr_gc", "gate_crispr_spacer"
    ),
    tiers = c(
      guide = "guide / PAM/PFS quality",
      seed = "variant position in spacer",
      specificity = "guide off-target count",
      amplicon = "amplicon length"
    ),
    gate_short = c(
      gate_crispr_pam = "PAM (Cas12a) / PFS (Cas13)",
      gate_crispr_seed = "seed coverage",
      gate_crispr_spacer = "spacer length",
      gate_crispr_specificity = "guide specificity",
      gate_crispr_gc = "guide GC",
      gate_crispr_amplicon = "amplicon length"
    ),
    gate_limit = c(
      gate_crispr_pam = "5'-TTTV (Cas12a); 3' non-G (Cas13a) / 5' D (Cas13b)",
      gate_crispr_seed = "variant in the seed (Cas12a) / central region (Cas13)",
      gate_crispr_spacer = "20-23 nt (Cas12a) / 26-28 nt (Cas13a) / 28-30 nt (Cas13b)",
      gate_crispr_specificity = "no off-target protospacer (spacer hit count = 1)",
      gate_crispr_gc = "0.35-0.65",
      gate_crispr_amplicon = "80-250 bp (RPA)"
    ),
    gate_floor_label = c(
      gate_crispr_amplicon = "widen the amplicon window",
      gate_crispr_gc = "widen the GC window",
      gate_crispr_spacer = "widen the spacer length window"
    ),
    gate_hard_label = c(
      gate_crispr_pam = "no Cas12a 5'-TTTV PAM or compatible Cas13 PFS = the guide is not selected",
      gate_crispr_seed = "a variant outside the seed / critical region is not discriminated",
      gate_crispr_specificity = "an off-target protospacer would fire on the wrong template"
    ),
    stage_label = "CRISPR guide"
  )
}
