# =============================================================================
#  marker_ssr.R  --  SSR / microsatellite genotyping (a size marker)
# =============================================================================
#  An SSR (simple sequence repeat, microsatellite / STR) marker is not an
#  allele-specific assay: *one* flanking pair amplifies the repeat locus and the
#  genotype is read from the amplicon *size* - the alleles differ by the number
#  of repeat copies.  So the design is the generic "any primer pair" core plus
#  the rules that make the *size* polymorphic and readable: the primers must sit
#  in the UNIQUE flanks, never over the repeat, and the repeat must be long
#  enough to vary between alleles.  The repeat span is located with
#  `ssr_span_at()` (a `repeat_spans()` lookup that also keeps the unit and the
#  copy number), the pair is placed with `flanking_pair_grid()`, and the amplicon
#  is kept in the ~100-300 bp size-marker window.
#
#  Design rules, grounded on the standard guidance:
#   * a legal primer pair must flank the repeat, in UNIQUE sequence - the repeat
#     itself is excluded from the primers (Primer3 manual: SEQUENCE_TARGET for a
#     CA repeat makes a legal pair flank the target and SEQUENCE_INTERNAL_EXCLUDED
#     keeps the repeat out of the oligo; Wikipedia "Microsatellite": "using the
#     unique sequences of flanking regions as primers").
#   * product ~100-300 bp (Primer3 default PRIMER_PRODUCT_SIZE_RANGE=100-300).
#   * composition / length / Tm / GC and the no-repeat, no-run rules are the
#     generic PCR ones: Primer3 defaults are 18-27 nt, Tm 57-63 C, GC 20-80%,
#     PRIMER_MAX_POLY_X=5; Premier Biosoft recommends 18-22 nt, Tm 52-58 C,
#     GC 40-60%, <= 4 dinucleotide repeats and <= 4 bp runs per primer.
#   * a *monomorphic* repeat is a useless marker: a dinucleotide needs many
#     copies before it varies, short tandem repeats are the least polymorphic,
#     so the locus is accepted only when the repeat reaches the polymorphic
#     floor (the detection floor is lower, so the gate can be relaxed, C8).
#   * long repeats stutter and short ones amplify poorly, so tetra-/penta-
#     nucleotide repeats are preferred for clean sizing (Wikipedia).  We do not
#     model stutter: it is a readout property, not a design gate here.
#   * multiplexing needs non-overlapping size ranges / different dyes; we design
#     a single locus and leave the pooling to the caller.
#   * the pair must amplify exactly ONE locus, and the README calls an off-target
#     with the whole 3' end fatal - so both specificity levels are hard gates:
#     the fast 3'-seed pre-screen on every candidate, and the full in-silico PCR
#     on the shortlist.
#  We reuse the package's composition window (KASP_LIMITS$gc, max_run,
#  has_repeat), so a primer is clean by construction; only the amplicon window,
#  the Tm window, the polymorphism floor and the specificity gates are
#  SSR-specific inputs below.
#
#  Sources (fetched 2026-10-06):
#   * https://primer3.org/manual.html
#     (SEQUENCE_TARGET, PRIMER_PRODUCT_SIZE_RANGE, PRIMER_MAX_POLY_X,
#      PRIMER_MIN_TM/PRIMER_MAX_TM, PRIMER_MIN_GC/PRIMER_MAX_GC)
#   * https://en.wikipedia.org/wiki/Microsatellite
#     (flanking unique primers; stutter / repeat-length trade-off)
#   * https://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
#     (18-22 nt, Tm 52-58 C, GC 40-60%, <= 4 dinucleotide repeats, <= 4 bp runs)
# =============================================================================

# The SSR acceptance windows.  `min_copies` is the detection floor (a repeat is
# called at 3 copies), `polymorphic_copies` the floor for a *usable size marker*
# (a 3-4 copy repeat is effectively monomorphic), so the two differ and the
# polymorphism gate below has a real, reachable range.
KASP_SSR_LIMITS <- list(
  amplicon = c(100L, 300L),
  tm = KASP_LIMITS$com_tm,
  gc = KASP_LIMITS$gc,
  min_copies = 3L,
  polymorphic_copies = 5L
)

# How far each *soft* SSR gate may be relaxed (bounded, never -Inf, C8): the
# amplicon window widens 50 bp each side, the Tm window 5 C each side, and the
# polymorphism floor drops to the 3-copy detection floor.
KASP_SSR_FLOOR <- list(
  amplicon = 50L, tm = 5,
  polymorphic = KASP_SSR_LIMITS$polymorphic_copies - KASP_SSR_LIMITS$min_copies
)

# The SSR gates, in application order (grave / cheap first, C6): the flanks'
# 3'-seed uniqueness (cheap, fatal), then the marker-level criteria (is the
# repeat long enough to vary; the pair's windows), then the expensive full
# in-silico PCR.  `flanks` and `offtarget` are hard: a primer with a duplicated
# 3' end, or a pair that amplifies several loci, gives an unreadable size peak.
KASP_GATES_SSR <- c(
  "gate_ssr_flanks", "gate_ssr_polymorphic", "gate_ssr_amplicon",
  "gate_ssr_tm", "gate_ssr_offtarget"
)
# Relaxation order - least grave first (C6): a wider amplicon or Tm window is a
# mild trade-off, then a shorter (less variable) repeat; the uniqueness gates are
# hard and never relaxed.
KASP_RELAX_SSR <- c("gate_ssr_amplicon", "gate_ssr_tm", "gate_ssr_polymorphic")
KASP_TIERS_SSR <- c(tier1 = "Tm balance", tier2 = "amplicon length", tier3 = "primer 3' end")

KASP_GATE_SHORT_SSR <- c(
  gate_ssr_flanks = "unique flanks",
  gate_ssr_polymorphic = "repeat variability",
  gate_ssr_amplicon = "amplicon length (SSR)",
  gate_ssr_tm = "primer Tm (SSR)",
  gate_ssr_offtarget = "off-target products"
)
KASP_GATE_LIMIT_SSR <- c(
  gate_ssr_flanks = "primers flank the repeat in unique sequence (3' seed unique)",
  gate_ssr_polymorphic = glue::glue(
    "repeat >= {KASP_SSR_LIMITS$polymorphic_copies} copies (length-polymorphic)"
  ),
  gate_ssr_amplicon = glue::glue(
    "amplicon within [{KASP_SSR_LIMITS$amplicon[1L]}, {KASP_SSR_LIMITS$amplicon[2L]}] bp"
  ),
  gate_ssr_tm = glue::glue(
    "both primers Tm within [{KASP_SSR_LIMITS$tm[1L]}, {KASP_SSR_LIMITS$tm[2L]}] C"
  ),
  gate_ssr_offtarget = "exactly one in-silico product"
)
KASP_GATE_FLOOR_LABEL_SSR <- c(
  gate_ssr_amplicon = "widen the amplicon window by 50 bp",
  gate_ssr_tm = "widen the Tm window by 5 C",
  gate_ssr_polymorphic = "accept a 3-copy repeat (monomorphic risk)"
)
KASP_GATE_HARD_LABEL_SSR <- c(
  gate_ssr_flanks = "a duplicated 3' end misprimes and splits the size peak",
  gate_ssr_offtarget = "a pair with a second product reads as a false allele"
)

#' The SSR span nearest (or containing) an anchor in a template
#'
#' Reuses `repeat_spans()` for the hit list and returns the widest/nearest span,
#' keeping the repeat `unit` and its `copies` so the polymorphism gate can be
#' scored: a span containing `anchor` wins (the smallest one), otherwise the span
#' whose midpoint is closest.
#' @param seq the template (plus strand).
#' @param anchor template coordinate of the locus anchor.
#' @param units repeat-unit lengths to search for.
#' @param min_copies tandem copies required for a repeat call.
#' @return a one-row tibble (`start`, `end`, `unit`, `copies`) in template
#'   coordinates, or `NULL`.
#' @keywords internal
#' @noRd
ssr_span_at <- function(seq, anchor, units = 2:4, min_copies = KASP_SSR_LIMITS$min_copies) {
  s <- as.character(seq)
  hits <- purrr::map(units, function(u) {
    m <- repeat_spans(s, u, min_copies)
    if (!nrow(m)) {
      return(NULL)
    }
    tibble::tibble(start = as.integer(m[, 1L]), end = as.integer(m[, 2L]), unit = as.integer(u))
  }) |>
    purrr::list_rbind()
  if (!nrow(hits)) {
    return(NULL)
  }
  hits$copies <- (hits$end - hits$start + 1L) %/% hits$unit
  inside <- dplyr::filter(hits, .data$start <= anchor, .data$end >= anchor)
  if (nrow(inside)) {
    return(dplyr::slice_min(dplyr::arrange(inside, .data$unit), .data$end - .data$start, n = 1L))
  }
  hits |>
    dplyr::mutate(dist = abs((.data$start + .data$end) / 2 - anchor)) |>
    dplyr::slice_min(.data$dist, n = 1L) |>
    dplyr::select(-"dist")
}

#' Design an SSR / microsatellite marker at a repeat locus
#'
#' The size-marker analogue of [design_hrm_marker()]: a single flanking pair
#' straddles the repeat, which sits under neither primer, and the amplicon is
#' sized to the ~100-300 bp window a capillary/gel readout wants.  The repeat is
#' located around `snp$pos`, the pair is placed in the unique flanks (so both
#' primers lie outside the repeat and carry no homopolymer >= 4 or SSR), and the
#' alleles are read as the amplicon size.  Every gate the scheme declares is
#' evaluated here (C7): the flanks' 3'-seed uniqueness, the repeat's
#' length-polymorphism, the amplicon and Tm windows and the full in-silico PCR
#' uniqueness.  A soft gate may be relaxed to its floor, least grave first, and
#' the result records which were (`$selection`).
#'
#' `snp$pos` is only the *anchor* of the window searched for the repeat - an SSR
#' locus is defined by its repeat, not by the REF/ALT bases.  Strong stutter is
#' expected for long repeats and multiplexing needs non-overlapping size ranges,
#' but neither is modelled here: stutter is a readout property, and pooling is
#' left to the caller.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`); `pos` anchors the
#'   window searched for the repeat.
#' @param flank template half-width.
#' @param min_product,max_product the accepted amplicon window (SSR size marker).
#' @param lengths primer lengths to try.
#' @param units repeat-unit lengths to search for.
#' @param min_copies tandem copies required for a repeat call.
#' @param tm_range,gc_range accepted primer windows.
#' @param ... passed to [in_silico_pcr()] (e.g. `BPPARAM`).
#' @return a `marker_ssr` list (`snp`, `primers`, `amplicon`, `repeat`,
#'   `selection`), or `NULL` (with a message) when no repeat or no flanking pair
#'   passes the gates.  The `primers`, `amplicon` and `repeat` coordinates are
#'   genomic.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_ssr_marker(genome, marker_example("SNP")[1, ])
#' }
#' @family KASP SSR
#' @export
design_ssr_marker <- function(genome, snp, flank = 250L,
                              min_product = KASP_SSR_LIMITS$amplicon[1L],
                              max_product = KASP_SSR_LIMITS$amplicon[2L],
                              lengths = 18:25, units = 2:4,
                              min_copies = KASP_SSR_LIMITS$min_copies,
                              tm_range = KASP_SSR_LIMITS$tm, gc_range = KASP_SSR_LIMITS$gc, ...) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  span <- ssr_span_at(tmpl$seq, tmpl$snp_i, units, min_copies)
  if (is.null(span)) {
    cli::cli_inform(gettext("No SSR / microsatellite repeat found in the flanking window."))
    return(NULL)
  }
  amp_pad <- KASP_SSR_FLOOR$amplicon
  tm_pad <- KASP_SSR_FLOOR$tm
  # the pair must flank the repeat; primer_candidates() also drops any primer
  # with a homopolymer >= 4 or an SSR, so a returned primer is clean by
  # construction (the "clean" rule is not restated as a gate - that would be
  # decorative, C7)
  grid <- flanking_pair_grid(
    tmpl$seq,
    lo = span$start, hi = span$end, lengths = lengths,
    tm_range = c(tm_range[1L] - tm_pad, tm_range[2L] + tm_pad), gc_range = gc_range
  )
  if (is.null(grid)) {
    cli::cli_inform(gettext("No flanking primer pair outside the repeat passed the Tm/GC filters."))
    return(NULL)
  }
  grid <- grid |>
    dplyr::filter(dplyr::between(.data$product, min_product - amp_pad, max_product + amp_pad)) |>
    dplyr::mutate(
      m_flanks = seed_margin(.data$primer_f, .data$primer_r, genome),
      m_polymorphic = as.numeric(span$copies - KASP_SSR_LIMITS$polymorphic_copies),
      m_amplicon = pmin(.data$product - min_product, max_product - .data$product),
      m_tm = pmin(
        .data$tm_f - tm_range[1L], tm_range[2L] - .data$tm_f,
        .data$tm_r - tm_range[1L], tm_range[2L] - .data$tm_r
      ),
      # the full specificity check needs the genome, so it is confirmed on the
      # shortlist; the seed margin above is the cheap pre-screen
      m_offtarget = 0,
      tier1 = -abs(.data$tm_f - .data$tm_r),
      tier2 = .data$m_amplicon,
      tier3 = -pmax(.data$run_f, .data$run_r)
    )
  confirm <- function(b) nrow(in_silico_pcr(b$primer_f, b$primer_r, genome, ...)) == 1L
  sel <- select_flanking_pair(
    grid, KASP_GATES_SSR, KASP_RELAX_SSR, names(KASP_TIERS_SSR),
    margins = c(
      gate_ssr_flanks = "m_flanks", gate_ssr_polymorphic = "m_polymorphic",
      gate_ssr_amplicon = "m_amplicon", gate_ssr_tm = "m_tm",
      gate_ssr_offtarget = "m_offtarget"
    ),
    floors = c(
      gate_ssr_amplicon = -amp_pad, gate_ssr_tm = -tm_pad,
      gate_ssr_polymorphic = -as.numeric(KASP_SSR_FLOOR$polymorphic)
    ),
    ties = c("dtm", "product"), confirm = confirm
  )
  if (is.null(sel$best)) {
    # name the gate that actually doomed the design, not a blanket reason
    if (any(!is.na(grid$m_flanks)) && !any(grid$m_flanks >= 0, na.rm = TRUE)) {
      cli::cli_inform(gettext("No flanking pair sits in unique enough sequence (every candidate's 3' seed repeats)."))
    } else {
      cli::cli_inform(gettext("No flanking pair passed the SSR gates (windows, variability, one product)."))
    }
    return(NULL)
  }
  primers <- pair_to_list(sel$best)
  off <- snp$pos - tmpl$snp_i
  primers$start <- primers$start + off
  primers$end <- primers$end + off
  out <- list(
    snp = snp,
    primers = primers,
    amplicon = tibble::tibble(start = primers$start, end = primers$end),
    selection = list(relaxed = sel$relaxed, violated = sel$violated)
  )
  # the motif itself is part of the marker: a plot must show `(AC)5`, not an
  # opaque "2nt x5", and the sequence is what the user can verify in the genome
  motif <- substr(as.character(tmpl$seq), span$start, span$end)
  # `repeat` is an R reserved word, so the field can only be set/read with [[
  out[["repeat"]] <- tibble::tibble(
    start = span$start + off,
    end = span$end + off,
    unit = span$unit,
    copies = span$copies,
    motif = motif,
    unit_seq = substr(motif, 1L, span$unit)
  )
  structure(out, class = "marker_ssr")
}

#' @export
print.marker_ssr <- function(x, ...) {
  cat(
    "<marker_ssr>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  rep_span <- x[["repeat"]]
  cat(
    "  repeat:", rep_span$start, "-", rep_span$end,
    "(", rep_span$end - rep_span$start + 1L, "bp )\n",
    sep = ""
  )
  if (!is.null(rep_span$copies)) {
    unit_seq <- if ("unit_seq" %in% names(rep_span)) rep_span$unit_seq[1L] else NULL
    motif <- if ("motif" %in% names(rep_span)) rep_span$motif[1L] else NULL
    formula <- if (!is.null(unit_seq) && !is.na(unit_seq)) {
      sprintf("(%s)%d", unit_seq, rep_span$copies)
    } else {
      sprintf("%d nt x %d", rep_span$unit, rep_span$copies)
    }
    seq <- if (!is.null(motif) && !is.na(motif)) paste0("  ", motif) else ""
    cat("  motif:", formula, seq, "(size marker)\n")
  }
  cat("  amplicon:", x$primers$product, "bp\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  if (length(x$selection$relaxed)) {
    cat("  relaxed:", paste(x$selection$relaxed, collapse = ", "), "\n")
  }
  invisible(x)
}

#' The SSR / microsatellite scheme definition
#'
#' The registry row and decision-schema fragments for the size marker: the
#' flanking primers must sit in *unique* sequence outside the repeat, the repeat
#' must be long enough to vary between alleles, the amplicon stays in the size
#' window, the pair is Tm-balanced and it must amplify exactly one locus.  Every
#' gate declared here is evaluated by [design_ssr_marker()] (C7).
#' @return a list (`scheme`, `family`, `label`, `discrimination`, `readout`,
#'   `zygosity`, `amplicon`, `loci`, `verdict`, `gates`, `relax`, `tiers`,
#'   `gate_short`, `gate_limit`, `gate_floor_label`, `gate_hard_label`,
#'   `stage_label`).
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_ssr
#' @export
#' @noRd
marker_map_spec.marker_ssr <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  tract <- if (!is.null(marker[["repeat"]])) {
    rp <- marker[["repeat"]]
    motif <- if ("motif" %in% names(rp)) rp$motif[1L] else NULL
    unit_seq <- if ("unit_seq" %in% names(rp)) rp$unit_seq[1L] else NULL
    formula <- ssr_formula(unit_seq, rp$copies[1L])
    label <- if (is.null(motif)) formula else paste0(formula, "  ", motif)
    list(spec_row(
      "repeat tract", rp$start[1L], rp$end[1L], "+", "tract",
      label = label
    ))
  } else {
    list()
  }
  map_spec(
    marker,
    features = c(map_pair_rows(marker, sp[["start"]], sp[["end"]]), tract),
    sequences = list(forward = marker$primers$forward, reverse = marker$primers$reverse)
  )
}

scheme_def_ssr <- function() {
  list(
    scheme = "ssr",
    family = "ssr",
    label = "SSR / microsatellite (simple-sequence-repeat size marker)",
    discrimination = "amplicon size tracks the repeat copy number",
    readout = "capillary / gel fragment size",
    zygosity = "co-dominant",
    amplicon = "~100-300 bp",
    loci = 1L,
    verdict = "usable SSR marker",
    gates = KASP_GATES_SSR,
    relax = KASP_RELAX_SSR,
    tiers = KASP_TIERS_SSR,
    gate_short = KASP_GATE_SHORT_SSR,
    gate_limit = KASP_GATE_LIMIT_SSR,
    gate_floor_label = KASP_GATE_FLOOR_LABEL_SSR,
    gate_hard_label = KASP_GATE_HARD_LABEL_SSR,
    stage_label = "SSR / microsatellite"
  )
}
