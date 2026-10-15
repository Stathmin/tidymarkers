# =============================================================================
#  kasp_crispr.R  --  CRISPR-Dx (SHERLOCK / DETECTR) as a marker scheme
# =============================================================================
#  A CRISPR-Dx marker does not read the allele from a primer 3' base (KASP) or a
#  band size (CAPS/ARMS): an isothermal pre-amplification (RPA / LAMP) of the
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
#  window (0.35-0.65, around the common 40-60% ideal) and the Cas13 "critical"
#  region (the central half of the spacer, the least mismatch-tolerant region).
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

#' Design a CRISPR-Dx (SHERLOCK / DETECTR) marker at a variant
#'
#' The Cas-mediated detection analogue of [assemble_kasp_assay()]: an isothermal
#' pre-amplification pair (RPA / LAMP) straddles the variant and a guide RNA is
#' placed so its spacer covers the variant, with the variant in the seed or
#' critical region.  For Cas12a a 5'-TTTV PAM must license the protospacer; Cas13
#' is PAM-free.  The collateral-cleavage readout itself is not simulated, so this
#' returns a well-formed guide plus its pre-amplification amplicon, not a
#' predicted signal.
#'
#' The guide is chosen best-first among all candidates whose pre-amplification
#' pair also passes the Tm/GC filters, so a guide is only returned when its
#' amplicon is designable too.
#' @param genome a `kasp_genome` or `DNAStringSet`.
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
#' @param gc_range accepted guide GC window.
#' @param min_product,max_product accepted pre-amplification amplicon window.
#' @param ... reserved for future options.
#' @return a `kasp_crispr` list (`snp`, `cas`, `strand`, `guide`, `pam`,
#'   `guide_start`, `guide_end`, `seed_offset`, `gc`, `amplicon`), or `NULL`
#'   (with a message) when no guide / amplicon passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_crispr_marker(genome, kasp_example("SNP")[1, ], cas = "cas12a")
#' design_crispr_marker(genome, kasp_example("SNP")[1, ], cas = "cas13")
#' }
#' @family KASP CRISPR
#' @export
design_crispr_marker <- function(genome, snp, flank = 200L, cas = c("cas12a", "cas13"),
                                 strand = c("both", "+", "-"),
                                 spacer_lengths = NULL, gc_range = NULL,
                                 min_product = 80L, max_product = 250L, ...) {
  cas <- match.arg(cas)
  strand <- match.arg(strand)
  strands <- if (identical(strand, "both")) c("+", "-") else strand
  if (is.null(spacer_lengths)) {
    spacer_lengths <- if (identical(cas, "cas12a")) 20:23 else 22:28
  }
  if (is.null(gc_range)) {
    gc_range <- c(0.35, 0.65)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  cand <- crispr_guide_candidates(tmpl$seq, tmpl$snp_i, cas, spacer_lengths, gc_range, strands)
  if (!nrow(cand)) {
    cli::cli_inform(gettext("No guide spacer passed the PAM / GC filters for this CRISPR-Dx marker."))
    return(NULL)
  }
  off <- snp$pos - tmpl$snp_i
  # design pre-amplification pairs best-guide-first, stopping at the first that works
  pick <- NULL
  for (i in seq_len(nrow(cand))) {
    amp <- design_flanking_pair(tmpl$seq, cand$start[i], cand$end[i],
      min_product = min_product, max_product = max_product
    )
    if (!is.null(amp)) {
      pick <- list(row = cand[i, ], amp = amp)
      break
    }
  }
  if (is.null(pick)) {
    cli::cli_inform(gettext("A guide spacer was found, but no pre-amplification pair passed the Tm/GC filters."))
    return(NULL)
  }
  row <- pick$row
  amp <- pick$amp
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
    amplicon = list(
      forward = amp$forward, reverse = amp$reverse,
      start = amp$start + off, end = amp$end + off, product = amp$product
    )
  )
  structure(out, class = "kasp_crispr")
}

#' @export
print.kasp_crispr <- function(x, ...) {
  cat(
    "<kasp_crispr>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, " [", x$cas, "]\n",
    sep = ""
  )
  pam <- if (is.na(x$pam)) "PAM-free" else paste0("PAM ", x$pam)
  cat("  guide:", x$guide, "(", x$strand, ",", pam, ", seed offset", x$seed_offset, ")\n")
  cat("  pre-amplification amplicon:", x$amplicon$product, "bp\n")
  cat("  F:", x$amplicon$forward, "\n  R:", x$amplicon$reverse, "\n")
  invisible(x)
}

#' The CRISPR-Dx scheme definition (internal)
#'
#' The registry entry the shared decision schema consumes, one marker chemistry
#' per call.  It is discovered automatically (one file per scheme, see
#' `R/kasp_scheme_registry.R`), so it is kept internal.
#' @return a list with the scheme identity, readout, gates, relax set, tiers and
#'   the label / limit lookups the decision graph and text export need.
#' @keywords internal
#' @noRd
scheme_def_crispr <- function() {
  list(
    scheme = "crispr",
    family = "crispr",
    label = "CRISPR-Dx (SHERLOCK / DETECTR)",
    discrimination = "Cas collateral cleavage of a reporter, guided to the variant",
    readout = "fluorescence / lateral flow (collateral reporter)",
    zygosity = "co-dominant (dominant for presence/absence)",
    amplicon = "pre-amplified (RPA / LAMP)",
    loci = 1L,
    verdict = "usable CRISPR-Dx marker",
    gates = c(
      "gate_crispr_spacer", "gate_crispr_pam", "gate_crispr_seed",
      "gate_crispr_gc", "gate_crispr_amplicon"
    ),
    relax = c("gate_crispr_spacer", "gate_crispr_gc", "gate_crispr_amplicon"),
    tiers = c(
      guide = "guide / PAM quality",
      seed = "variant position in spacer",
      amplicon = "amplicon length"
    ),
    gate_short = c(
      gate_crispr_spacer = "spacer length",
      gate_crispr_pam = "PAM (TTTV)",
      gate_crispr_seed = "seed coverage",
      gate_crispr_gc = "guide GC",
      gate_crispr_amplicon = "amplicon length"
    ),
    gate_limit = c(
      gate_crispr_spacer = "20-23 nt (Cas12a) / 22-28 nt (Cas13)",
      gate_crispr_pam = "5'-TTTV adjacent to the protospacer (Cas12a)",
      gate_crispr_seed = "variant within the seed / central critical region",
      gate_crispr_gc = "0.35-0.65",
      gate_crispr_amplicon = "80-250 bp (RPA / LAMP)"
    ),
    gate_floor_label = c(
      gate_crispr_spacer = "widen the spacer length window",
      gate_crispr_gc = "widen the GC window",
      gate_crispr_amplicon = "widen the amplicon window"
    ),
    gate_hard_label = c(
      gate_crispr_pam = "no 5'-TTTV PAM = Cas12a cannot be licensed",
      gate_crispr_seed = "a variant outside the seed is not discriminated"
    ),
    stage_label = "CRISPR guide"
  )
}
