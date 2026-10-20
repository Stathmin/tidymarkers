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
#     spacer is ~20-23 nt and target binding unleashes indiscriminate ssDNase
#     ("collateral") cleavage read by a quenched reporter.
#     Chen JS et al. (2018) Science 360:436-439, doi:10.1126/science.aar6245.
#     Addgene CRISPR Guide, PAM table (AsCpf1/LbCpf1: 5' TTTV; seed at the
#     PAM-proximal end): https://www.addgene.org/guides/crispr/
#   * Cas13 (SHERLOCK) is PAM-free (it targets ssRNA, so the DNA PAM is not
#     required), tolerates a mismatch, and its crRNA spacer is ~20-28 nt; an
#     isothermal pre-amplification precedes T7 transcription and detection.
#     Gootenberg JS et al. (2017) Science 356:438-442, doi:10.1126/science.aam9321;
#     Gootenberg JS et al. (2018) Science 360:439-444, doi:10.1126/science.aaq0179
#     (SHERLOCKv2, 20-28 nt spacers).  Addgene CRISPR Guide (Cas13 targets RNA,
#     no PAM): https://www.addgene.org/guides/crispr/
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

#' Candidate guide spacers covering a variant
#'
#' Enumerates guide spacers of the accepted lengths that place the variant at
#' each offset of the *critical* window, on both strands.  For `cas = "cas12a"` a
#' 5'-TTTV PAM immediately adjacent to the protospacer is required and the
#' critical window is the PAM-proximal seed; for `cas = "cas13"` there is no PAM
#' and the critical window is the central region of the spacer.  Candidates are
#' ordered best-first (deepest in the seed / most central, then GC closest to
#' 0.5).
#' @param seq the plus-strand template.
#' @param snp_i template coordinate of the variant.
#' @param cas `"cas12a"` or `"cas13"`.
#' @param spacer_lengths spacer lengths to try.
#' @param gc_range accepted guide GC window.
#' @param strands strands the guide may be placed on (`"+"`/`"-"`).
#' @return a tibble (`start`, `end`, `len`, `strand`, `spacer`, `pam`, `gc`,
#'   `v_off`), plus-strand coordinates, or zero rows.
#' @keywords internal
#' @noRd
crispr_guide_candidates <- function(seq, snp_i, cas, spacer_lengths, gc_range, strands = c("+", "-")) {
  strands <- match.arg(strands, c("+", "-"), several.ok = TRUE)
  s <- as.character(seq)
  n <- nchar(s)
  seed_len <- 8L # Cas12a PAM-proximal seed; also the length of the offsets tried
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
  cand <- cand |>
    dplyr::mutate(
      plus_spacer = stringr::str_sub(s, .data$start, .data$end),
      pam_raw = stringr::str_sub(s, .data$pam_lo, .data$pam_hi)
    )
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  cand <- cand |>
    dplyr::mutate(
      spacer = dplyr::if_else(.data$strand == "+", .data$plus_spacer, rc(.data$plus_spacer)),
      pam = dplyr::if_else(.data$strand == "+", .data$pam_raw, rc(.data$pam_raw)),
      gc = gc_frac(.data$spacer)
    )
  if (identical(cas, "cas12a")) {
    cand <- dplyr::filter(cand, stringr::str_detect(.data$pam, "^TTT[ACG]$"))
  }
  cand <- cand |>
    dplyr::filter(dplyr::between(.data$gc, gc_range[1L], gc_range[2L])) |>
    dplyr::mutate(gc_dev = abs(.data$gc - 0.5))
  if (identical(cas, "cas12a")) {
    cand <- dplyr::arrange(cand, .data$v_off, .data$gc_dev)
  } else {
    cand <- dplyr::mutate(cand, pos_dev = abs(.data$v_off - (.data$len - 1) / 2)) |>
      dplyr::arrange(.data$pos_dev, .data$gc_dev)
  }
  dplyr::select(
    cand, "start", "end", "len", "strand", "spacer", "pam", "gc", "v_off"
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

#' Design a CRISPR-Dx (SHERLOCK / DETECTR) marker at a variant
#'
#' The Cas-mediated detection analogue of [assemble_kasp_assay()]: an isothermal
#' RPA pre-amplification pair straddles the variant and a guide RNA is
#' placed so its spacer covers the variant, with the variant in the seed or
#' critical region.  For Cas12a a 5'-TTTV PAM must license the protospacer; Cas13
#' is PAM-free.  The guide must also be unique in the genome (no off-target
#' protospacer).  The collateral-cleavage readout itself is not simulated, so this
#' returns a well-formed guide plus its pre-amplification amplicon, not a
#' predicted signal.
#'
#' The guide is chosen best-first among all candidates whose pre-amplification
#' pair also passes *that chemistry's* composition filters (RPA, never a KASP
#' PCR), so a guide is only returned when its amplicon is designable too.
#' The soft gates (amplicon window, GC, spacer length) are relaxed one at a
#' time, least grave first, only when no guide clears them all; the PAM
#' (Cas12a), the seed/critical region and guide specificity (no off-target
#' protospacer) are hard - no guide is returned past them.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param cas `"cas12a"` (DETECTR; needs a PAM) or `"cas13"` (SHERLOCK;
#'   PAM-free).
#' @param strand guide placement: `"both"` (default), or `"+"`/`"-"`.  Cas13
#'   cleaves ssRNA, so its guide must match the *transcribed* strand (the one the
#'   pre-amplification's T7 promoter turns into the target RNA); pass `strand`
#'   when that orientation is known, otherwise both are tried and the choice is
#'   the caller's to validate.
#' @param spacer_lengths guide spacer lengths to try (default: 20-23 nt for
#'   Cas12a, 22-28 nt for Cas13).
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
#' @return a `marker_crispr` list (`snp`, `cas`, `strand`, `guide`, `pam`,
#'   `guide_start`, `guide_end`, `seed_offset`, `gc`, `amplicon`, `off_target`,
#'   `gates`, `relaxed`), or `NULL` (with a message) when no guide / amplicon
#'   passes the hard gates.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_crispr_marker(genome, marker_example("SNP")[1, ], cas = "cas12a")
#' design_crispr_marker(genome, marker_example("SNP")[1, ], cas = "cas13")
#' }
#' @family KASP CRISPR
#' @export
design_crispr_marker <- function(genome, snp, flank = 200L, cas = c("cas12a", "cas13"),
                                 strand = c("both", "+", "-"),
                                 spacer_lengths = NULL, gc_range = NULL,
                                 specificity_max = 0L,
                                 preamp = "rpa",
                                 min_product = 80L, max_product = 250L,
                                 haplotypes = NULL, BPPARAM = NULL) {
  cas <- match.arg(cas)
  strand <- match.arg(strand)
  preamp <- match.arg(preamp, "rpa")
  strands <- if (identical(strand, "both")) c("+", "-") else strand
  strict_spacer <- if (is.null(spacer_lengths)) {
    if (identical(cas, "cas12a")) 20:23 else 22:28
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
  # (C6), one rung at a time.  PAM, seed and specificity are hard: a guide
  # without a PAM (Cas12a), with the variant outside the seed/critical region,
  # or with an off-target protospacer cannot be licensed to discriminate, so
  # they are enforced by the enumeration / hit count and never relaxed.
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
      "No guide spacer passed the PAM / seed / specificity / GC filters with a pre-amplification pair."
    ))
    return(NULL)
  }
  row <- pick$row
  amp <- pick$amp
  off_target <- if (is.na(pick$hits)) NA_integer_ else pick$hits - 1L
  # the gate table the DAG displays, evaluated on the guide actually returned (C7)
  gates <- c(
    gate_crispr_pam = if (identical(cas, "cas12a")) grepl("^TTT[ACG]$", row$pam) else NA,
    gate_crispr_seed = TRUE, # enforced by the seed / central-window enumeration
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
  pam <- if (is.na(x$pam)) "PAM-free" else paste0("PAM ", x$pam)
  cat("  guide:", x$guide, "(", x$strand, ",", pam, ", seed offset", x$seed_offset, ")\n")
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
    discrimination = "Cas collateral cleavage of a reporter, guided to the variant",
    readout = "fluorescence / lateral flow (collateral reporter)",
    zygosity = "co-dominant (dominant for presence/absence)",
    amplicon = "pre-amplified (RPA)",
    loci = 1L,
    verdict = "usable CRISPR-Dx marker",
    # the gate order the method applies (grave / cheap first): a guide is
    # licensed by the PAM and the seed first, then its length, specificity, GC
    # and finally the pre-amplification amplicon it rides on.
    # The single DAG describes the stricter Cas12a path for `gate_crispr_pam`:
    # the gate is hard, but it evaluates to `NA` for PAM-free Cas13.
    gates = c(
      "gate_crispr_pam", "gate_crispr_seed", "gate_crispr_spacer",
      "gate_crispr_specificity", "gate_crispr_gc", "gate_crispr_amplicon"
    ),
    relax = c(
      "gate_crispr_amplicon", "gate_crispr_gc", "gate_crispr_spacer"
    ),
    tiers = c(
      guide = "guide / PAM quality",
      seed = "variant position in spacer",
      specificity = "guide off-target count",
      amplicon = "amplicon length"
    ),
    gate_short = c(
      gate_crispr_pam = "PAM (Cas12a; Cas13 PAM-free)",
      gate_crispr_seed = "seed coverage",
      gate_crispr_spacer = "spacer length",
      gate_crispr_specificity = "guide specificity",
      gate_crispr_gc = "guide GC",
      gate_crispr_amplicon = "amplicon length"
    ),
    gate_limit = c(
      gate_crispr_pam = "5'-TTTV adjacent (Cas12a); PAM-free (Cas13)",
      gate_crispr_seed = "variant in the seed (Cas12a) / central region (Cas13)",
      gate_crispr_spacer = "20-23 nt (Cas12a) / 22-28 nt (Cas13)",
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
      gate_crispr_pam = "no 5'-TTTV PAM = Cas12a cannot be licensed (Cas13 is PAM-free)",
      gate_crispr_seed = "a variant outside the seed / critical region is not discriminated",
      gate_crispr_specificity = "an off-target protospacer would fire on the wrong template"
    ),
    stage_label = "CRISPR guide"
  )
}
