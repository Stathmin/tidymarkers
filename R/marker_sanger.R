# =============================================================================
#  marker_sanger.R  --  Sanger sequencing (dideoxy read-out) as a marker scheme
# =============================================================================
#  Sanger reads the variant as a *base in a chromatogram*, not as a ratio or a
#  band, so there is no allele-specific primer and no restriction site: one
#  flanking primer pair delimits an amplicon that is read with a sequencing
#  primer.  The design rules below are the consensus of the sequencing-core
#  guidance surveyed for this file:
#
#   * amplicon ~400-800 bp around the variant, with the variant >= ~50 bp from
#     the primer binding sites so it sits inside the high-quality read and never
#     under a primer;
#     https://duckduckgo.com/html/?q=Sanger+sequencing+amplicon+design+primer+400-800+bp+variant+distance+from+primer
#     (Sanger validation guide: amplicon 400-800 bp, variant >= 50 bp from the
#     primer binding sites, specificity checks against the reference genome).
#   * primer length 17-25 nt, Tm ~55-65 C (50-70 acceptable), GC ~50% (we keep
#     the package GC window 40-60%), no hairpins / dimers, avoid long runs and
#     >2 identical 3' bases (Cornell BRC Sanger primer-design guide);
#   * the first ~50 bases of a read are low quality (a 650 bp read with 50 bp
#     low-quality edges leaves a 550 bp usable window), which is the ">= 50 bp
#     from the read start" rule above;
#     https://molbiologytools.com/sanger_primer_designer_help
#   * the pair must amplify exactly one product genome-wide, or the chromatogram
#     is a mixture of two templates.
#
#  The designer itself never simulates a chromatogram: it returns a well-formed,
#  unique amplicon with the variant correctly placed, plus the sequencing primers
#  that read into it.  The read *physics* that are modelled there are geometric:
#  the variant must sit between `min_dist` and `read_length` bp of the sequencing
#  primer's 5' end.
#
#  The read-out half is separate and works on an *observed* 4-channel trace
#  supplied by the caller (this file never invents intensities):
#   * sanger_deconvolve() baseline-normalises the channels, subtracts an optional
#     bleed-through matrix, calls a primary and a secondary base per position and
#     reports the minor allele's fraction and a call quality;
#   * sanger_allele_fraction() reads one position;
#   * sanger_indel_estimate() is a deliberately simple heuristic for a heterozygous
#     indel: the point where the trace turns mixed fixes the indel site, the
#     length is the shift that best realigns the secondary channel (no ICE/TIDE).
#  The detection floor is one constant, SANGER_DETECTION: a minor allele below
#  floor (~15%) is never called a heterozygote (it is `unknown`), a call above
#  mixed (~20%) is `mixed`, and the band in between is reported as `unknown`.
#  The floor is a literature default, not a calibration of the caller's
#  instrument, hence `calibrated = FALSE` - every read-out record carries it.
#
#  Selection is real (not decorative): every gate declared by scheme_def_sanger()
#  is computed here and applied by select_with_ladder() (C7).  The cheap margins
#  (variant distance, amplicon window, Tm, 3'-end composition, 3'-seed
#  uniqueness) are computed for every candidate pair; the expensive ones (primer
#  hairpin/homodimer, the full in-silico PCR) are confirmed on the shortlist, so
#  a design never pays for an OligoArrayAux call per candidate.
# =============================================================================

# The Sanger acceptance windows (the constraints the gates below encode).  Tm is
# the sequencing-core 55-65 C rather than the tighter KASP COM window.  A read
# is usable from `min_dist` (the low-quality lead) to `read_length` bp past the
# sequencing primer's 5' end; past that the base quality / peak spacing is no
# longer reliable enough for a heterozygous call.
KASP_SANGER_LIMITS <- list(
  amplicon = c(400L, 800L),
  tm = c(55, 65),
  gc = c(0.40, 0.60),
  min_dist = 50L,
  read_length = 700L
)

# The chromatogram read-out limits - the single source of truth for every
# sanger_*() read-out helper, docs and tests.  `floor` is the minor-allele
# fraction below which a secondary peak is not distinguishable from noise: a
# mixed call there would be a false heterozygote, so it is `unknown`.  `mixed` is
# the fraction from which a heterozygote is called (between the default ~15-20%
# visibility limit of a minor peak in a dideoxy trace).  `calibrated = FALSE`:
# both are a literature default for a generic capillary instrument, not a
# calibration of the caller's own trace (a real calibration needs panels of
# known mixes, e.g. TIDE/ICE's dilution series).
SANGER_DETECTION <- list(
  floor = 0.15,
  mixed = 0.20,
  calibrated = FALSE
)

# How far each *soft* Sanger gate may be relaxed - the bounded floor the ladder
# stops at (never -Inf, C8).  The amplicon window widens 200 bp each side, the Tm
# window 5 C each side, and the 3'-end rule tolerates a 3-base 3' homopolymer;
# the primer structure confirmation tolerates 3 kcal/mol more stable a duplex.
KASP_SANGER_FLOOR <- list(amplicon = 200L, tm = 5, primer = 1L, structure = 3)

# The Sanger selection gates, in *application* order (grave / cheap first, C6):
# the variant must be readable (a hard, cheap placement test), then the pair's
# cheap windows, then the 3'-end / structure quality, then the expensive
# genome-wide uniqueness.  `dist` and `unique` are hard: a variant under a primer
# and a non-unique amplicon cannot be read at all.
KASP_GATES_SANGER <- c(
  "gate_sanger_dist", "gate_sanger_amplicon", "gate_sanger_tm",
  "gate_sanger_primer", "gate_sanger_unique"
)
# The *relaxation* order - least grave first (C6, and different from the gate
# order on purpose): widening the amplicon or Tm window is a mild trade-off, then
# a tolerated 3' homopolymer (mispriming risk); the placement and uniqueness
# gates are hard and never relaxed.
KASP_RELAX_SANGER <- c("gate_sanger_amplicon", "gate_sanger_tm", "gate_sanger_primer")
KASP_TIERS_SANGER <- c(tier1 = "variant centrality", tier2 = "Tm balance", tier3 = "amplicon length")

# Short names for the decision graph.
KASP_GATE_SHORT_SANGER <- c(
  gate_sanger_dist = "variant from read start",
  gate_sanger_amplicon = "amplicon length (Sanger)",
  gate_sanger_tm = "primer Tm (Sanger)",
  gate_sanger_primer = "primer 3' end / structure",
  gate_sanger_unique = "unique amplicon"
)

# Human-readable limit, relaxation floor and hard-stop reason per gate.
KASP_GATE_LIMIT_SANGER <- c(
  gate_sanger_dist = glue::glue("variant >= {KASP_SANGER_LIMITS$min_dist} bp from each read start"),
  gate_sanger_amplicon = glue::glue(
    "amplicon within [{KASP_SANGER_LIMITS$amplicon[1L]}, {KASP_SANGER_LIMITS$amplicon[2L]}] bp"
  ),
  gate_sanger_tm = glue::glue(
    "both primers Tm within [{KASP_SANGER_LIMITS$tm[1L]}, {KASP_SANGER_LIMITS$tm[2L]}] C"
  ),
  gate_sanger_primer = "no 3' homopolymer >= 3 / no hairpin or homodimer",
  gate_sanger_unique = "exactly one in-silico amplicon; no off-target exact primer hits"
)
KASP_GATE_FLOOR_LABEL_SANGER <- c(
  gate_sanger_amplicon = "widen the amplicon window by 200 bp",
  gate_sanger_tm = "widen the Tm window by 5 C",
  gate_sanger_primer = "tolerate a 3-base 3' homopolymer"
)
KASP_GATE_HARD_LABEL_SANGER <- c(
  gate_sanger_dist = "a variant under a primer or in the low-quality read lead is unreadable",
  gate_sanger_unique = "a non-unique amplicon cannot be called"
)

# ---- shared selection for the flanking-pair schemes (sanger/ssr/indel_size) ---
# All three pick ONE flanking pair around a forbidden core and grade it the same
# way: enumerate forward x reverse pairs, give each declared gate a margin
# ("how far inside its threshold", >= 0 passes) and let select_with_ladder()
# apply the two declared orders; the expensive gates (structure, in-silico PCR)
# are confirmed on the shortlist.  The helpers live here, beside the first of the
# three, because they are shared by all three and the shared package files must
# not grow.

# The cap on the candidate pool per side: a >= 1 kb flank can offer ~500 starts x
# 8 lengths, so the full forward x reverse grid would be millions of rows.  Every
# position is kept (the placement tiers must see the whole flank); only the
# near-duplicate lengths are dropped, and the pool is capped by Tm closeness to
# the scheme's target - a bounded, documented pool instead of an O(starts^2) one.
KASP_MARKER_MAX_PER_SIDE <- 600L

#' Keep one length per start position, capped by Tm closeness to a target
#'
#' The bound that keeps `flanking_pair_grid()` affordable on a long flank:
#' positions are preserved (so the placement tiers still see the whole flank) and
#' for each only the length whose Tm is closest to `target` is kept.
#' @param cand a candidate primer tibble (see `primer_candidates()`).
#' @param target the preferred Tm (C).
#' @param max_per_side pool cap.
#' @return a thinned tibble.
#' @keywords internal
#' @noRd
trim_primers <- function(cand, target, max_per_side = KASP_MARKER_MAX_PER_SIDE) {
  if (!nrow(cand)) {
    return(cand)
  }
  cand <- cand[order(abs(cand$tm - target), cand$start, cand$len), ]
  utils::head(cand[!duplicated(cand$start), ], max_per_side)
}

#' Every forward x reverse pair flanking a core, with its raw metrics
#'
#' The candidate grid the flanking-pair schemes select from.  `fwd` must end
#' before `lo`, `rev` must start after `hi`, so both bind outside the core and the
#' variant (or repeat / indel) sits between them under neither primer.
#' @param seq the template (plus strand).
#' @param lo,hi inclusive plus-strand bounds of the forbidden core.
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted windows (already widened to the gate floor).
#' @param target the preferred Tm (C) the pool is thinned toward.
#' @return a tibble (one row per pair) with the suffix-joined primer metrics and
#'   `product`, `dtm`, or `NULL` when no pair can flank the core.
#' @keywords internal
#' @noRd
flanking_pair_grid <- function(seq, lo, hi, lengths, tm_range, gc_range,
                               target = mean(tm_range)) {
  n <- nchar(as.character(seq))
  if (lo < 1L || hi > n) {
    return(NULL)
  }
  fwd <- primer_candidates(seq, seq_len(max(lo - 1L, 0L)), lengths, "+", tm_range, gc_range) |>
    dplyr::filter(.data$end < lo) |>
    trim_primers(target) |>
    dplyr::mutate(run = terminal_run(.data$primer))
  rev <- primer_candidates(seq, seq.int(hi + 1L, n), lengths, "-", tm_range, gc_range) |>
    dplyr::filter(.data$start > hi) |>
    trim_primers(target) |>
    dplyr::mutate(run = terminal_run(.data$primer))
  if (!nrow(fwd) || !nrow(rev)) {
    return(NULL)
  }
  dplyr::cross_join(fwd, rev, suffix = c("_f", "_r")) |>
    dplyr::mutate(
      product = .data$end_r - .data$start_f + 1L,
      dtm = abs(.data$tm_f - .data$tm_r)
    )
}

#' Length of a primer's 3'-terminal run of one base
#'
#' The Sanger / generic primer rule forbids a 3' homopolymer longer than two
#' bases (a stable 3' run misprimes).  Reported as a count, so a gate can rank it.
#' @param primer oligo(s) 5'->3'.
#' @return an integer vector, one per primer.
#' @keywords internal
#' @noRd
terminal_run <- function(primer) {
  m <- stringr::str_match(primer, "([ACGT])\\1*$")[, 1L]
  ifelse(is.na(m), 0L, nchar(m))
}

#' The 3'-seed uniqueness margin of a primer pair
#'
#' The fast specificity pre-screen of the README: each primer's 3'-terminal seed
#' is counted genome-wide, and a pair whose both seeds occur exactly once scores
#' 0 (passes); each extra genome hit costs one margin unit.  `NA` when the genome
#' is an indexed `marker_genome` (`specificity_seed()` needs an in-memory
#' `DNAStringSet`), which the ladder reads as "unknown" and lets through to the
#' full in-silico PCR.
#' @param primer_f,primer_r the pair's oligos.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param k 3'-seed length.
#' @return a numeric vector, one per pair.
#' @keywords internal
#' @noRd
seed_margin <- function(primer_f, primer_r, genome, k = SEED_K) {
  if (inherits(genome, "marker_genome")) {
    return(rep(NA_real_, length(primer_f)))
  }
  primers <- unique(c(primer_f, primer_r))
  s <- specificity_seed(stats::setNames(primers, primers), genome, k = k, max_mismatch = 0L)
  hits <- rlang::set_names(s$hits, s$primer)
  1 - pmax(hits[primer_f], hits[primer_r])
}

#' Is a flanking pair free of a hairpin / homodimer (and palindromes)?
#'
#' The structural half of the primer-quality rule.  `TRUE` when OligoArrayAux is
#' absent (the pass is then "unknown", exactly as a missing `NA` margin would be,
#' so a bare machine is never blocked by a tool it cannot run).
#' @param primer_f,primer_r the pair's oligos.
#' @param limits thresholds (see `KASP_LIMITS`).
#' @return a single logical.
#' @keywords internal
#' @noRd
pair_structure_ok <- function(primer_f, primer_r, limits = KASP_LIMITS) {
  if (!has_oligoarrayaux()) {
    return(TRUE)
  }
  p <- c(primer_f, primer_r)
  if (any(has_palindrome(p))) {
    return(FALSE)
  }
  hairpin <- oligo_dg_vec(p)
  homodimer <- c(oligo_dg(p[[1L]], p[[1L]]), oligo_dg(p[[2L]], p[[2L]]))
  all(hairpin >= limits$hairpin_dg, na.rm = TRUE) &&
    all(homodimer >= limits$homodimer_dg, na.rm = TRUE)
}

#' Select one flanking pair through the ladder, then confirm its expensive gates
#'
#' Phase 1 is the shared `select_with_ladder()` over the *cheap* margins, so the
#' declared gate order, relaxation order and floors are the ones the DAG shows.
#' Phase 2 walks the tier-ordered survivors and applies the *expensive* gate(s)
#' via `confirm()` (OligoArrayAux structure, in-silico PCR): the first candidate
#' the confirmation accepts wins, so a design never runs the expensive test on
#' more than `max_alternatives` pairs, and a candidate the confirmation rejects is
#' dropped and the ladder re-run (which can also relax a soft gate further).
#' @param grid the candidate grid (see `flanking_pair_grid()`).
#' @param gates,relax,tiers the scheme's declared stage (C6).
#' @param margins,floors named maps gate -> margin column / gate -> floor.
#' @param ties ascending tie-break columns.
#' @param confirm `function(row)` -> `TRUE` (accept), `FALSE` (drop) or `NA`
#'   (unknown, accept); `NULL` skips the phase.
#' @param max_alternatives how many rejected candidates to walk past.
#' @return the `select_with_ladder()` result (`best`, `relaxed`, `violated`),
#'   with `best = NULL` when nothing passes (doomed).
#' @keywords internal
#' @noRd
select_flanking_pair <- function(grid, gates, relax, tiers, margins, floors,
                                 ties = character(0), confirm = NULL, max_alternatives = 8L) {
  grid <- dplyr::mutate(grid, .row = dplyr::row_number())
  pick <- function(d) {
    select_with_ladder(d, gates, relax, tiers, ties, gate_margin = margins, gate_floor = floors)
  }
  res <- pick(grid)
  if (is.null(confirm)) {
    return(res)
  }
  # walk the tier-ordered winners: the first the confirmation accepts wins; a
  # rejected one is dropped and the ladder re-run (which may relax a soft gate
  # further for the survivors), so the expensive test runs on few pairs
  for (i in seq_len(max_alternatives + 1L)) {
    if (is.null(res$best)) {
      return(res)
    }
    if (!isFALSE(confirm(res$best))) {
      return(res) # TRUE or NA (unknown) both accept
    }
    if (i > max_alternatives) {
      res$best <- NULL
      return(res)
    }
    grid <- dplyr::filter(grid, .data$.row != res$best$.row[[1L]])
    res <- pick(grid)
  }
  res
}

#' A sequencing primer that reads into the variant
#'
#' A Sanger read only sees the variant when it starts >= `min_dist` bp away from
#' it (the low-quality lead), so a candidate is kept only when `keep()` holds.
#' @param seq the template (plus strand).
#' @param side `"+"` (reads the plus strand, 3' end at the right) or `"-"`.
#' @param starts plus-strand start positions the primer may take.
#' @param lengths,tm_range,gc_range accepted windows.
#' @param keep the predicate the candidate must satisfy (the min_dist rule).
#' @return a one-row tibble or `NULL`.
#' @keywords internal
#' @noRd
seq_primer_at <- function(seq, side, starts, lengths, tm_range, gc_range, keep) {
  cand <- primer_candidates(seq, starts, lengths, side, tm_range, gc_range)
  cand <- cand[keep(cand), , drop = FALSE]
  if (!nrow(cand)) {
    return(NULL)
  }
  # the best-reading nested primer: Tm closest to the window centre, then the
  # lowest 3'-terminal run (cleanest 3' anchor)
  mid <- mean(tm_range)
  cand[order(abs(cand$tm - mid), terminal_run(cand$primer)), ][1L, ]
}

#' Design a sequencing primer pair that reads into the variant
#'
#' Returns the forward and reverse sequencing primers (each inside the amplicon,
#' at least `min_dist` from the variant, Tm-screened and free of a 3' homopolymer),
#' preferring a nested primer and falling back to the flanking pair (the routine
#' workflow).  A nested primer is kept only when its 3' seed is unique, so the
#' read is not a mixture.  The nested search window is deliberately narrow (<= 40
#' bp inside the PCR primer), so the design stays a handful of Tm calls.
#' @param seq the template (plus strand); `si` the variant's template position.
#' @param start_f,end_r plus-strand coordinates of the flanking pair's ends.
#' @param flank_f,flank_r the flanking PCR primers (the fallback sequencing ones).
#' @param min_dist minimal read-start distance.
#' @param read_length usable read span from the primer's 5' end; a variant past
#'   it cannot be called reliably from that read.
#' @param genome a `marker_genome` or `DNAStringSet` (for the seed screen).
#' @return a list (`forward`, `reverse`, `nested`, `usable`).
#' @keywords internal
#' @noRd
sanger_seq_primers <- function(seq, si, start_f, end_r, flank_f, flank_r, min_dist, genome,
                               lengths, tm_range, gc_range,
                               read_length = KASP_SANGER_LIMITS$read_length) {
  nested_lens <- intersect(c(20L, 22L), lengths)
  if (!length(nested_lens)) {
    nested_lens <- min(lengths)
  }
  window <- 40L
  fwd <- seq_primer_at(
    seq, "+", start_f + seq.int(0L, window), nested_lens, tm_range, gc_range,
    function(cand) cand$end <= si - min_dist & si - cand$start <= read_length
  )
  rev <- seq_primer_at(
    seq, "-", seq.int(end_r - window - max(nested_lens) + 1L, end_r - min(nested_lens) + 1L),
    nested_lens, tm_range, gc_range,
    function(cand) {
      cand$start >= si + min_dist & cand$end <= end_r & cand$end - si <= read_length
    }
  )
  if (!is.null(fwd) && !is.null(rev) && !inherits(genome, "marker_genome")) {
    p <- c(fwd$primer, rev$primer)
    if (any(seed_margin(p, p, genome) < 0)) {
      fwd <- rev <- NULL
    }
  }
  forward <- fwd$primer %||% flank_f
  reverse <- rev$primer %||% flank_r
  forward_start <- fwd$start %||% start_f
  forward_end <- fwd$end %||% (start_f + nchar(flank_f) - 1L)
  reverse_start <- rev$start %||% (end_r - nchar(flank_r) + 1L)
  reverse_end <- rev$end %||% end_r
  # A read is usable from `min_dist` (low-quality lead) to `read_length` bp of
  # the primer's 5' end.  The nested search already enforces the upper bound; a
  # *flanking* fallback can still be too far away, so mark the side instead of
  # presenting a pair that silently cannot call the variant.
  dist <- c(forward = si - forward_start, reverse = reverse_end - si)
  usable <- dist >= min_dist & dist <= read_length
  if (!all(usable)) {
    sides <- paste(names(usable)[!usable], collapse = " and ")
    msg <- gettext("No usable Sanger read reaches the variant (%s): more than %s bp from the read primer.")
    cli::cli_warn(c(
      sprintf(msg, sides, read_length),
      "i" = gettext("The amplicon is returned, but that side cannot call the variant.")
    ))
  }
  list(
    forward = forward,
    reverse = reverse,
    nested = !is.null(fwd) && !is.null(rev),
    # template coordinates of the actual read primers (the PCR pair when a
    # nested read primer was not found), converted to genomic by the designer
    forward_start = forward_start,
    forward_end = forward_end,
    reverse_start = reverse_start,
    reverse_end = reverse_end,
    usable = usable
  )
}

# =============================================================================
#  Read-out: deconvolving an observed 4-channel trace
# =============================================================================
#  A chromatogram is logged by the caller as one row per position and four
#  integer intensities; nothing here invents data.  The model is the standard
#  additive one used by basecallers:
#
#    I_c(pos) ~ baseline_c + bleed(x) + sum over alleles of f_a * peak_a(pos)
#
#  so after a local baseline subtraction and an optional fixed bleed-through
#  matrix the channel intensities are proportional to allele amount: the
#  primary (largest) channel is the major allele, the largest of the rest is the
#  minor allele, and the minor-allele fraction is second / (primary + second).
#  The fraction is a *linear* mixture, which is exactly what sanger_allele_fraction()
#  reports; the detection floor then decides what may be claimed from it.
#
#  Assumptions (all documented, none estimated):
#   * peaks are already tracked, so `pos` is one base of the read, not a time
#     sample, and a mixed position has both bases in the same `pos` row;
#   * the local baseline is the rolling minimum of each channel (the mean of the
#     peaks' minima when the trace is shorter than the window);
#   * bleed-through is a fixed, caller-supplied matrix (identity by default): it
#     is a property of the dye set / instrument, not something a single trace can
#     identify;
#   * `noise` (if not supplied) is the median of the local maxima of the
#     baseline-corrected channels - a per-trace proxy for the minor peak a
#     basecaller would ignore.

# The bases in chromatogram channel order.
SANGER_BASES <- c("A", "C", "G", "T")

# Rolling minimum of a numeric vector, centred, `half` positions each side.
sanger_rolling_min <- function(x, half) {
  n <- length(x)
  i <- seq_len(n)
  vapply(i, function(k) min(x[pmax(1L, k - half):pmin(n, k + half)]), numeric(1))
}

# The caller's trace as an n x 4 numeric matrix, with the position labels.
sanger_channels <- function(trace) {
  if (is.matrix(trace)) {
    trace <- as.data.frame(trace)
  }
  if (!is.data.frame(trace) && !is.list(trace)) {
    cli::cli_abort(gettext("`trace` must be a data.frame or a matrix (one row per position)."))
  }
  missing <- setdiff(SANGER_BASES, names(trace))
  if (length(missing)) {
    cli::cli_abort(gettext("`trace` is missing channel column(s): {paste(missing, collapse = ', ')}."))
  }
  x <- vapply(SANGER_BASES, function(b) as.numeric(trace[[b]]), numeric(nrow(trace)))
  x <- matrix(x, ncol = 4L)
  if (!nrow(x)) {
    cli::cli_abort(gettext("`trace` has no positions."))
  }
  if (any(!is.finite(x))) {
    cli::cli_abort(gettext("`trace` channel intensities must all be finite."))
  }
  list(
    x = x,
    pos = if ("pos" %in% names(trace)) as.integer(trace$pos) else seq_len(nrow(x)),
    base = if ("base" %in% names(trace)) toupper(as.character(trace$base)) else NULL
  )
}

# Baseline-corrected channels folded for the local baseline of each channel.
# `baseline` is one of "local" (rolling minimum at `baseline_window`), "none", or
# a numeric scalar subtracted from every channel (a dark-current offset).
sanger_correct_channels <- function(ch, baseline = "local", baseline_window = 5L,
                                    bleed = NULL) {
  x <- ch$x
  # a constant per-channel baseline (the channel's minimum over the read) is the
  # dark-current floor; a *rolling* baseline would subtract a local minor peak
  # together with its own baseline and erase exactly what we are looking for.
  # With a bleed matrix the correction comes first: subtracting the floor before
  # it would mix the floor of one channel into the next.
  # with a bleed matrix the floor is subtracted first (it is a per-channel
  # additive offset, not a dye signal) and subtracted again after the correction
  # in case the leakage left a residue
  observed <- correct_minus_baseline(if (is.null(bleed)) x else pmax(x, 0), baseline)
  corrected <- observed
  if (!is.null(bleed)) {
    if (!is.matrix(bleed) || !identical(dim(bleed), c(4L, 4L))) {
      cli::cli_abort(gettext("`bleed` must be a 4 x 4 matrix (the correction that undoes the dye overlap)."))
    }
    # corrected = observed x bleed, row by row with plain vectors (matrix
    # subassignment is avoided: `corr[row, ] <- vals` silently keeps only the
    # last element on this build, and `rbind` of row vectors is unambiguous)
    rows <- vector("list", nrow(observed))
    for (row in seq_len(nrow(observed))) {
      obs <- observed[row, ]
      vals <- numeric(4L)
      for (col in seq_len(4L)) {
        vals[col] <- sum(obs * bleed[, col])
      }
      rows[[row]] <- vals
    }
    corrected <- do.call(rbind, rows)
    corrected <- pmax(corrected, 0)
  }
  # baseline subtraction, row by row with plain vectors: matrix arithmetic on
  # this build's matrices is not reliable, and `pb` is a vector by construction
  corrected <- correct_minus_baseline(corrected, baseline)
  corrected
}

# `max(x - baseline_row, 0)` applied row by row (see sanger_correct_channels).
correct_minus_baseline <- function(x, baseline) {
  n <- nrow(x)
  known <- identical(baseline[[1L]], "local") || identical(baseline[[1L]], "none")
  numeric_ok <- is.numeric(baseline) && length(baseline) == 1L && is.finite(baseline)
  if (!known && !numeric_ok) {
    cli::cli_abort(gettext("`baseline` must be \"local\", \"none\" or one finite number."))
  }
  mins <- if (identical(baseline[[1L]], "local")) {
    # the instrument floor is the smallest value any channel reaches; a
    # per-channel minimum can *be* the minor peak (a channel that is minor at
    # most positions has its smallest value at a real peak)
    rep(min(x), ncol(x))
  } else if (identical(baseline[[1L]], "none")) {
    rep(0, ncol(x))
  } else {
    rep(as.numeric(baseline), ncol(x))
  }
  rows <- vapply(seq_len(n), function(row) {
    pmax(x[row, ] - mins, 0)
  }, numeric(ncol(x)))
  t(rows)
}

# The per-trace noise proxy: the median of the channels' *interquartile* level is
# too high (it would swallow a 20% peak), so use the smallest above-baseline
# signal a channel ever carries - the floor a basecaller would ignore.  The
# `noise` argument of sanger_deconvolve() overrides it.
sanger_noise <- function(x) {
  above <- apply(x, 2L, min)
  stats::median(above) + 1e-9
}

#' Deconvolve an observed Sanger trace into per-position allele calls
#'
#' The read-out half of the Sanger scheme.  Takes a 4-channel chromatogram
#' (one row per position; `A`/`C`/`G`/`T` intensities), removes a local baseline
#' and an optional bleed-through, then calls a primary (major) and a secondary
#' (minor) base per position and reports the minor allele's fraction.  The
#' fraction is a linear mixture, so no deconvolution of base composition is
#' needed - the only hard part is deciding *what the fraction licenses*.
#'
#' The detection floor is [SANGER_DETECTION] (`floor` ~15%):
#' * `status = "mixed"` when the minor fraction is at or above `mixed` (~20%)
#'   and the minor peak stands above the trace noise - a heterozygote;
#' * `status = "unknown"` for a minor peak below `mixed` (it may be a real minor
#'   allele under the visibility limit, or noise).  It is never reported as a
#'   confident homozygote, and a trace with no secondary peak at all is
#'   `status = "homozygous"`;
#' * `status = "noise"` when the secondary peak is not above the trace noise.
#'
#' `ref`/`alt` are the two called bases in conventional order (primary first),
#' not the reference genome's alleles: the trace carries no genome position, so
#' callers wanting a REF/ALT frame should compare `ref`/`alt` with the target
#' row themselves.  `calibrated = FALSE` says the floor is a literature default,
#' not a calibration of the caller's instrument.
#' @param trace a data.frame or matrix with `A`, `C`, `G`, `T` intensity columns
#'   and an optional `pos` or `base` column (the label; `pos` defaults to
#'   `seq_len(nrow(trace))`).
#' @param baseline `"local"` (a constant per-channel minimum), `"none"`, or one
#'   number subtracted from every channel (a dark-current offset).
#' @param baseline_window unused; kept so callers that pass it do not break.
#' @param bleed optional 4 x 4 *correction* matrix (dye by observed channel):
#'   `corrected = observed %*% bleed`, i.e. the inverse of the instrument's
#'   leakage matrix; `NULL` (no bleed-through) by default.
#' @param noise the minor-peak height treated as noise; defaults to the median of
#'   the per-channel maxima of the corrected trace.
#' @param floor,mixed the minor-allele fraction limits; default [SANGER_DETECTION].
#' @return a tibble (`pos`, `base`, `ref`, `alt`, `alt_fraction`, `mixed`,
#'   `call_quality`, `status`, `calibrated`) with one row per position.
#' @examples
#' # a 50/50 A/G heterozygote: every position carries both peaks at half height
#' het <- sanger_example_trace("AG", fraction = 0.5)
#' sanger_deconvolve(het)
#' @family KASP Sanger
#' @export
sanger_deconvolve <- function(trace, baseline = "local", baseline_window = 5L,
                              bleed = NULL, noise = NULL,
                              floor = SANGER_DETECTION$floor,
                              mixed = SANGER_DETECTION$mixed) {
  ch <- sanger_channels(trace)
  x <- sanger_correct_channels(ch, baseline, baseline_window, bleed)
  nz <- if (is.null(noise)) sanger_noise(x) else as.numeric(noise)[1L]
  # an explicit argmax: `max.col()` is not reliable on this build's 4-column
  # double matrices, and this loop is also trivially readable
  top <- vapply(seq_len(nrow(x)), function(i) which.max(x[i, ]), integer(1))
  second <- vapply(seq_len(nrow(x)), function(i) {
    v <- x[i, ]
    v[top[i]] <- -Inf
    which.max(v)
  }, integer(1))
  p <- vapply(seq_len(nrow(x)), function(i) x[i, top[i]], numeric(1))
  # A 50/50 mixture ties the two channels, and the argmax then has no basis to
  # pick a primary base.  The trace's own `base` (the basecaller's primary call)
  # breaks the tie when the trace carries one.
  if (!is.null(ch$base)) {
    k <- match(ch$base, SANGER_BASES)
    for (i in seq_len(nrow(x))) {
      if (is.na(k[i]) || top[i] == k[i] || x[i, k[i]] != p[i]) next
      top[i] <- k[i]
      v <- x[i, ]
      v[top[i]] <- -Inf
      second[i] <- which.max(v)
    }
    p <- vapply(seq_len(nrow(x)), function(i) x[i, top[i]], numeric(1))
  }
  s <- vapply(seq_len(nrow(x)), function(i) x[i, second[i]], numeric(1))
  tot <- p + s
  frac <- ifelse(tot > 0, s / tot, NA_real_)
  # `above`: the secondary peak stands above the trace noise (a real second
  # peak, not a bump).  `visible`: it also reaches the detection floor.  Only a
  # `visible` peak is eligible to be called `mixed`; a visible peak under the
  # mixed threshold is `unknown` - never `homozygous`.
  above <- s > nz & p > nz
  # A minor peak that clears the noise but stays under `floor` is not visible:
  # it cannot be told from a bump, so the position is `noise` - a minor allele
  # below the floor is simply not callable from this trace.  From `floor` up to
  # `mixed` the peak is real but too low for a heterozygote, so the call is
  # `unknown` (never a confident homozygote); from `mixed` up it is `mixed`.
  visible <- above & frac >= floor & frac <= 1 - floor
  mixed_ok <- visible & frac >= mixed & frac <= 1 - mixed
  status <- rep("homozygous", nrow(x))
  status <- ifelse(above & !visible, "noise", status)
  status <- ifelse(visible & !mixed_ok, "unknown", status)
  status <- ifelse(mixed_ok, "mixed", status)
  frac <- ifelse(status == "homozygous", 0, frac)
  alt <- ifelse(status == "homozygous", "", SANGER_BASES[second])
  tibble::tibble(
    pos = ch$pos,
    base = if (is.null(ch$base)) NA_character_ else ch$base,
    ref = SANGER_BASES[top],
    alt = alt,
    alt_fraction = round(frac, 4L),
    mixed = mixed_ok,
    call_quality = round(ifelse(tot > 0, (p - s) / tot, 0), 4L),
    status = status,
    calibrated = SANGER_DETECTION$calibrated
  )
}

#' The minor-allele fraction at one position of a Sanger trace
#'
#' The one-position read of [sanger_deconvolve()], with the detection floor
#' spelled out: `above_floor` is `TRUE` only when the minor peak clears the trace
#' noise *and* reaches [SANGER_DETECTION]`$mixed`; between `floor` and `mixed` the
#' status is `unknown` (a real minor allele cannot be excluded, a confident
#' homozygote cannot be claimed), below `floor` it is `homozygous` only when no
#' secondary peak stands above the noise at all.  Never a confident homozygous
#' call for a visible minor peak under the floor.
#' @param trace a 4-channel trace, as for [sanger_deconvolve()].
#' @param pos the position to read (matched against the trace's `pos` column; a
#'   plain index when the trace has no `pos`).
#' @param ... passed to [sanger_deconvolve()] (`baseline`, `bleed`, ...).
#' @return a tibble (`pos`, `ref`, `alt`, `alt_fraction`, `call_quality`,
#'   `status`, `above_floor`, `floor`, `calibrated`), one row, or `NULL` with a
#'   message when the position is not in the trace.
#' @examples
#' sanger_allele_fraction(sanger_example_trace("AG", fraction = 0.2), 1L)
#' @family KASP Sanger
#' @export
sanger_allele_fraction <- function(trace, pos, ...) {
  if (length(pos) != 1L) {
    cli::cli_abort(gettext("`pos` must be one position."))
  }
  d <- sanger_deconvolve(trace, ...)
  row <- d[d$pos == as.integer(pos), , drop = FALSE]
  if (!nrow(row)) {
    cli::cli_inform(gettext("Position {pos} is not in the Sanger trace."))
    return(NULL)
  }
  # the plain minor share, independent of the call: a 5% minor allele keeps its
  # 0.05 even though it is not called a heterozygote
  tibble::tibble(
    pos = row$pos, ref = row$ref, alt = row$alt,
    alt_fraction = minor_fraction(trace, row$pos, ...),
    call_quality = row$call_quality,
    status = row$status,
    above_floor = row$status == "mixed",
    floor = SANGER_DETECTION$floor,
    calibrated = SANGER_DETECTION$calibrated
  )
}

# The minor-channel share second / (primary + second) at one position, with no
# floor applied (the number, before the detection rule interprets it).
minor_fraction <- function(trace, pos, ...) {
  ch <- sanger_channels(trace)
  x <- sanger_correct_channels(ch, ...)
  i <- match(as.integer(pos), ch$pos)
  v <- x[i, ]
  top <- which.max(v)
  s <- max(v[-top])
  round(s / (v[top] + s), 4L)
}

#' Estimate a heterozygous indel from a Sanger trace
#'
#' A deliberate heuristic, not a deconvolution (no TIDE / ICE / Poly Peak
#' Parser).  A heterozygous indel superposes two sequences from the indel
#' onwards, so the trace stops being a single base per position at the indel and
#' the secondary channel starts to follow a *shifted* copy of the primary one.
#' The function finds the first position after the clean lead where the trace is
#' mixed ([sanger_deconvolve()]`$status`), then scores every signed shift `d` in
#' `-max_size:max_size` by the correlation between the minor-channel fraction of
#' the interior positions and the primary-channel fraction shifted by `d` (each
#' channel scaled by the calling strength, so the two halves compare like for
#' like).  The best `d` is the indel size: negative = deletion, positive =
#' insertion.
#'
#' The assumptions are strong and stated on purpose: the two alleles are equally
#' amplified from `start` on (a fraction well below the detection floor cannot be
#' seen at all), the trace is not re-phased by the caller (some basecallers infer
#' a new spacing after a mixed call, which defeats the realignment), and the
#' lead is clean so the first mixed position is the indel site.  The result says
#' `confidence = "heuristic"`; it is a hypothesis to test, not a call.
#' @param trace a 4-channel trace, as for [sanger_deconvolve()].
#' @param start first position of the *clean lead* the indel follows; the search
#'   for the indel site starts here (default: the first position of the trace).
#' @param max_size largest indel size to score, in bases.
#' @param window number of positions after the site used to realign.
#' @param floor the minor fraction below which a position is not mixed; default
#'   [SANGER_DETECTION]`$floor`.
#' @param ... passed to [sanger_deconvolve()].
#' @return a tibble (`pos`, `indel_size`, `heterozygous_fraction`,
#'   `n_positions`, `correlation`, `confidence`, `calibrated`), one row; `pos` and
#'   `indel_size` are `NA` (with a message) when no mixed stretch is found.
#' @examples
#' trace <- sanger_example_trace(c("ACGTACGTAAAA", "ACGTACGTCCCC"), fraction = 0.5, shift = -4L)
#' sanger_indel_estimate(trace)
#' @family KASP Sanger
#' @export
sanger_indel_estimate <- function(trace, start = NULL, max_size = 10L, window = 40L,
                                  floor = SANGER_DETECTION$floor, ...) {
  d <- sanger_deconvolve(trace, floor = floor, mixed = floor, ...)
  # the indel site is the first mixed position of the *longest* clean lead: a
  # low-level noise bump before the real indel would otherwise be picked
  hit <- d$status != "homozygous" & d$status != "noise"
  if (!is.null(start)) {
    hit <- hit & d$pos >= start[1L]
  }
  if (!any(hit)) {
    cli::cli_inform(gettext("No mixed stretch above the detection floor: the trace has no visible heterozygous indel."))
    return(tibble::tibble(
      pos = NA_integer_, indel_size = NA_integer_, heterozygous_fraction = NA_real_,
      n_positions = 0L, correlation = NA_real_, confidence = "heuristic",
      calibrated = SANGER_DETECTION$calibrated
    ))
  }
  i0 <- which(hit)[1L]
  # score over a window that starts `max_size` before the site, so every candidate
  # shift has real overlap; the reported site stays the first mixed position
  lo <- max(1L, i0 - as.integer(max_size))
  idx <- seq.int(lo, min(nrow(d), i0 + as.integer(window) - 1L))
  # the minor allele's share, read only where the trace is actually mixed (a
  # deletion leaves homologous positions homozygous at 0, which would drag a
  # plain median of the window to 0)
  mixed_idx <- which(hit & seq_along(hit) >= i0 & seq_along(hit) <= i0 + as.integer(window) - 1L)
  frac <- if (length(mixed_idx)) stats::median(d$alt_fraction[mixed_idx], na.rm = TRUE) else NA_real_
  sizes <- seq.int(-as.integer(max_size), as.integer(max_size))
  # for a shift `sz` the minor allele's base at j equals the major allele's base
  # at j - sz, so count how many mixed positions agree with the major base at
  # the shifted source: the largest count is the indel size
  score <- vapply(sizes, function(sz) {
    src <- idx - sz
    ok <- idx >= i0 & idx <= nrow(d) & src >= 1L & src <= nrow(d) & nzchar(d$alt[idx])
    if (sum(ok) < 3L) {
      return(0)
    }
    sum(d$alt[idx[ok]] == d$ref[src[ok]]) / sum(ok)
  }, numeric(1))
  best <- sizes[which.max(score)]
  tibble::tibble(
    pos = as.integer(d$pos[i0]),
    indel_size = as.integer(best),
    heterozygous_fraction = round(frac, 4L),
    n_positions = sum(hit),
    correlation = round(max(score), 4L),
    confidence = "heuristic",
    calibrated = SANGER_DETECTION$calibrated
  )
}

#' A deterministic example Sanger trace
#'
#' Builds a 4-channel trace in the linear-mixture model [sanger_deconvolve()]
#' assumes, for tests and examples: a primary sequence, an optional second
#' sequence (the minor allele) superimposed from its first position, and an
#' optional shift of the second sequence from a given position (a heterozygous
#' insertion, positive, or deletion, negative - the second sequence is the one
#' with the indel).  Intensities are deterministic (a flat peak per base with a
#' small fixed noise floor), so every caller sees the same numbers.
#' @param primary the primary (major) base sequence.
#' @param secondary an optional superimposed sequence; `NULL` for a plain read.
#' @param fraction the minor allele's fraction of the mixture.
#' @param shift signed length of the second sequence's indel, in bases.
#' @param shift_at position of `secondary` (in `primary` coordinates) where the
#'   minor allele's own bases start; the minor copy is
#'   `secondary[shift_at:length(secondary)]` placed from `shift_at + shift`
#'   (positive = the minor allele carries an insertion, negative = a deletion).
#' @param peak the major-allele peak height; every channel also carries a flat
#'   1% of it as an instrument floor, so no position is exactly zero.
#' @return a tibble with `pos`, `base` and the four channel columns.
#' @examples
#' sanger_example_trace("ACGTACGT", fraction = 0.3)
#' @family KASP Sanger
#' @export
sanger_example_trace <- function(primary, secondary = NULL, fraction = 0.5, shift = 0L,
                                 shift_at = 1L, peak = 1000) {
  primary <- toupper(strsplit(paste(primary, collapse = ""), "")[[1L]])
  second <- if (is.null(secondary)) {
    character(0)
  } else {
    toupper(strsplit(paste(secondary, collapse = ""), "")[[1L]])
  }
  # the minor copy is the tail of `secondary` from `shift_at`, displaced by
  # `shift`: before `shift_at` the two alleles are identical, so the trace is
  # clean there and turns mixed exactly at the indel
  second <- second[seq_along(second) >= shift_at]
  at <- as.integer(shift_at) + as.integer(shift) + seq_along(second) - 1L
  # the read has exactly the primary's length: a minor allele longer than the
  # primary (an insertion) is clipped past the read end
  second <- second[at <= length(primary)]
  at <- at[at <= length(primary)]
  n <- length(primary)
  x <- matrix(0, n, 4L)
  add <- function(mat, letters, positions, weight) {
    for (k in seq_along(positions)) {
      pos <- positions[k]
      code <- match(letters[k], SANGER_BASES)
      if (is.na(code) || pos < 1L || pos > n) next
      old <- mat[pos, code]
      mat[pos, code] <- old + weight
    }
    mat
  }
  frac_s <- if (length(second)) fraction else 0
  frac_p <- if (length(second)) 1 - fraction else 1
  x <- add(x, primary, seq_along(primary), frac_p * peak)
  x <- add(x, second, at, frac_s * peak)
  # a per-position flat floor of 1% of the major peak (below the detection floor
  # by a wide margin): present in every channel, so it never distorts a ratio
  x <- x + 0.01 * (if (frac_p > 0) peak else 0)
  out <- tibble::as_tibble(as.data.frame(x, stringsAsFactors = FALSE))
  names(out) <- SANGER_BASES
  out <- tibble::add_column(out, pos = seq_len(n), .before = 1L)
  out <- tibble::add_column(out, base = c(primary, rep(NA_character_, n - length(primary))), .after = 1L)
  out
}

#' Design a Sanger sequencing marker at a variant
#'
#' The amplicon analogue of [assemble_kasp_assay()]: a single flanking pair
#' delimits a ~400-800 bp product, the variant sits at least `min_dist` bp
#' inside both ends (clear of the low-quality read lead and of either primer),
#' the pair is checked for genome-wide uniqueness with [in_silico_pcr()], and the
#' sequencing primers that read into the variant are returned.  There is no
#' allele-specific primer - Sanger reads the variant directly - so the design
#' optimises the amplicon window, the variant placement and the Tm balance, then
#' gates on the primer quality and off-target uniqueness.  Every gate the scheme
#' declares is evaluated here (C7); a soft one may be relaxed to its floor, in the
#' declared order, and the result records which were (`$selection`).
#'
#' It does **not** simulate the chromatogram, so it returns a well-formed,
#' unique amplicon, not a predicted trace; the trace side is
#' [sanger_deconvolve()].
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width (must leave room for the amplicon).
#' @param min_dist minimal distance from the variant to each amplicon end (the
#'   read start / primer binding site).
#' @param read_length usable read span from the sequencing primer's 5' end.  A
#'   variant farther than this is not reliably callable; the design still
#'   returns, but warns and flags that side in `$seq_primers$usable`.
#' @param min_product,max_product the accepted amplicon window.
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows.
#' @param check_unique when `TRUE`, require exactly one [in_silico_pcr()] product.
#' @param population optional population table already used by the package's
#'   haplotype-aware markers (`haplotype`, `freq`, `chrom`, `pos`, `ref`, `alt`,
#'   see [marker_attach_haplotypes()]).  A known variant under a PCR or read
#'   primer interferes with that primer's binding (allele dropout), so it is
#'   **warned and flagged** in `$allele_dropout` + `$allele_dropout_risk`, not
#'   passed silently.  When absent the design still returns, but
#'   `$allele_dropout_risk` is `NA` and a warning says the check was not run.
#' @param ... passed to [in_silico_pcr()] (e.g. `BPPARAM`).
#' @return a `marker_sanger` list (`snp`, `primers`, `amplicon`, `seq_primers`,
#'   `selection`, `allele_dropout`, `allele_dropout_risk`), or `NULL` (with a
#'   message) when no unique pair passes the gates.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_sanger_marker(genome, marker_example("SNP")[1, ], flank = 500)
#' }
#' @family KASP Sanger
#' @export
design_sanger_marker <- function(genome, snp, flank = 500L,
                                 min_dist = KASP_SANGER_LIMITS$min_dist,
                                 read_length = KASP_SANGER_LIMITS$read_length,
                                 min_product = KASP_SANGER_LIMITS$amplicon[1L],
                                 max_product = KASP_SANGER_LIMITS$amplicon[2L],
                                 lengths = 18:25, tm_range = KASP_SANGER_LIMITS$tm,
                                 gc_range = KASP_SANGER_LIMITS$gc,
                                 check_unique = TRUE, population = NULL, ...) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  si <- tmpl$snp_i
  amp_pad <- KASP_SANGER_FLOOR$amplicon
  tm_pad <- KASP_SANGER_FLOOR$tm
  # the pair flanks the variant (fwd ends before it, rev starts after it); the
  # hard `dist` gate below then keeps only the pairs that leave the variant clear
  # of both read starts, so no guard band is applied here (it would be vacuous)
  grid <- flanking_pair_grid(
    tmpl$seq,
    lo = si, hi = si, lengths = lengths,
    tm_range = c(tm_range[1L] - tm_pad, tm_range[2L] + tm_pad), gc_range = gc_range
  )
  if (is.null(grid)) {
    cli::cli_inform(gettext("No flanking pair kept the variant clear of both read starts for Sanger."))
    return(NULL)
  }
  grid <- grid |>
    dplyr::filter(dplyr::between(.data$product, min_product - amp_pad, max_product + amp_pad)) |>
    dplyr::mutate(
      m_dist = pmin(si - .data$end_f, .data$start_r - si) - as.integer(min_dist),
      m_amplicon = pmin(.data$product - min_product, max_product - .data$product),
      m_tm = pmin(
        .data$tm_f - tm_range[1L], tm_range[2L] - .data$tm_f,
        .data$tm_r - tm_range[1L], tm_range[2L] - .data$tm_r
      ),
      m_primer = 2L - pmax(.data$run_f, .data$run_r),
      m_unique = seed_margin(.data$primer_f, .data$primer_r, genome),
      tier1 = -abs(si - (.data$start_f + .data$end_r) / 2),
      tier2 = -abs(.data$tm_f - .data$tm_r),
      tier3 = .data$m_amplicon
    )
  confirm <- function(b) {
    if (check_unique && nrow(in_silico_pcr(b$primer_f, b$primer_r, genome, ...)) != 1L) {
      return(FALSE)
    }
    pair_structure_ok(b$primer_f, b$primer_r)
  }
  sel <- select_flanking_pair(
    grid, KASP_GATES_SANGER, KASP_RELAX_SANGER, names(KASP_TIERS_SANGER),
    margins = c(
      gate_sanger_dist = "m_dist", gate_sanger_amplicon = "m_amplicon",
      gate_sanger_tm = "m_tm", gate_sanger_primer = "m_primer",
      gate_sanger_unique = "m_unique"
    ),
    floors = c(
      gate_sanger_amplicon = -amp_pad, gate_sanger_tm = -tm_pad,
      gate_sanger_primer = -as.numeric(KASP_SANGER_FLOOR$primer)
    ),
    ties = c("dtm", "product"), confirm = confirm
  )
  if (is.null(sel$best)) {
    # name the gate that actually doomed the design, not a blanket reason
    if (any(!is.na(grid$m_unique)) && !any(grid$m_unique >= 0, na.rm = TRUE)) {
      cli::cli_inform(gettext("No flanking pair has unique 3' ends (every candidate's 3' seed repeats in the genome)."))
    } else {
      cli::cli_inform(gettext("No flanking pair kept the variant clear of both read starts for Sanger."))
    }
    return(NULL)
  }
  primers <- pair_to_list(sel$best)
  seq_primers <- sanger_seq_primers(
    tmpl$seq, si, primers$start, primers$end,
    primers$forward, primers$reverse, min_dist, genome, lengths, tm_range, gc_range,
    read_length = read_length
  )
  amp0 <- snp$pos - si
  for (nm in c("forward_start", "forward_end", "reverse_start", "reverse_end")) {
    seq_primers[[nm]] <- seq_primers[[nm]] + amp0
  }
  out <- list(
    snp = snp, primers = primers,
    amplicon = tibble::tibble(start = primers$start + amp0, end = primers$end + amp0),
    seq_primers = seq_primers,
    selection = list(relaxed = sel$relaxed, violated = sel$violated)
  )
  out <- structure(out, class = "marker_sanger")
  if (is.null(population)) {
    out$allele_dropout <- tibble::tibble()
    out$allele_dropout_risk <- NA
    cli::cli_warn(c(
      gettext("Sanger: no population variants supplied; primer-over-variant allele dropout is not checked."),
      "i" = gettext("`allele_dropout_risk` is unknown, not a pass.")
    ))
    return(out)
  }
  h <- normalise_haplotypes(population)
  hits <- sanger_primer_variants(out, h)
  out$allele_dropout <- hits
  out$allele_dropout_risk <- nrow(hits) > 0L
  if (nrow(hits)) {
    cli::cli_warn(c(
      gettext("Sanger: {nrow(hits)} population variant(s) fall under a primer (allele-dropout risk)."),
      "i" = gettext(
        "A variant in a primer binding site can lose that allele (or the whole product); see `$allele_dropout`."
      )
    ))
  }
  marker_attach_haplotypes(out, population, snp)
}

#' Population variants that fall under a designed Sanger primer
#'
#' Allele dropout is a *primer-binding* problem: a variant in the 17-25 nt a
#' primer anneals to can weaken that primer for the haplotype carrying it, so the
#' amplicon (or that allele) disappears in a fraction of the population.  Every
#' primer's binding site is mapped to plus-strand coordinates (the reverse
#' primers cover `start..start + length - 1` on the minus strand, the same span
#' on the plus strand) and every population variant inside a span is reported.
#' The target position itself is excluded: the variant is placed at least
#' `min_dist` from both read starts and never under a primer by construction.
#' @param marker a `marker_sanger` object.
#' @param haplotypes a normalised population table.
#' @return a tibble (`primer`, `start`, `end`, `pos`, `haplotype`, `freq`, `risk`),
#'   one row per variant per primer, or an empty tibble.
#' @keywords internal
#' @noRd
sanger_primer_variants <- function(marker, haplotypes) {
  h <- haplotypes[haplotypes$chrom == marker$snp$chrom[1L] & haplotypes$pos != marker$snp$pos[1L], ,
    drop = FALSE
  ]
  if (!nrow(h)) {
    return(tibble::tibble())
  }
  amp <- marker$amplicon
  sp <- marker$seq_primers
  spans <- tibble::tibble(
    primer = c("PCR_F", "PCR_R", "read_F", "read_R"),
    start = c(
      amp$start[1L],
      amp$end[1L] - nchar(marker$primers$reverse) + 1L,
      sp$forward_start %||% amp$start[1L],
      sp$reverse_start %||% (amp$end[1L] - nchar(sp$reverse) + 1L)
    )
  ) |>
    dplyr::mutate(end = .data$start + c(
      nchar(marker$primers$forward), nchar(marker$primers$reverse),
      nchar(sp$forward), nchar(sp$reverse)
    ) - 1L) |>
    dplyr::distinct(.data$primer, .data$start, .data$end)
  hits <- lapply(seq_len(nrow(spans)), function(i) {
    inside <- h$pos >= spans$start[i] & h$pos <= spans$end[i]
    if (!any(inside)) {
      return(NULL)
    }
    tibble::tibble(
      primer = spans$primer[i], start = spans$start[i], end = spans$end[i],
      pos = h$pos[inside], haplotype = h$haplotype[inside], freq = h$freq[inside],
      risk = "allele dropout"
    )
  })
  out <- purrr::list_rbind(hits)
  if (is.null(out)) tibble::tibble() else out
}

#' @export
print.marker_sanger <- function(x, ...) {
  cat(
    "<marker_sanger>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  amplicon:", x$primers$product, "bp (Sanger read)\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  if (!is.null(x$seq_primers)) {
    how <- if (isTRUE(x$seq_primers$nested)) "nested" else "flanking"
    cat("  sequencing primers (", how, "):\n", sep = "")
    cat("    F:", x$seq_primers$forward, "\n    R:", x$seq_primers$reverse, "\n")
    if (!is.null(x$seq_primers$usable) && !all(x$seq_primers$usable)) {
      bad <- names(x$seq_primers$usable)[!x$seq_primers$usable]
      cat(
        "  unusable read:", paste(bad, collapse = ", "),
        "(variant beyond read reach)\n"
      )
    }
  }
  if (length(x$selection$relaxed)) {
    cat("  relaxed:", paste(x$selection$relaxed, collapse = ", "), "\n")
  }
  if (isTRUE(x$allele_dropout_risk)) {
    cat("  allele-dropout risk:", nrow(x$allele_dropout), "population variant(s) under a primer\n")
  } else if (identical(x$allele_dropout_risk, NA)) {
    cat("  allele-dropout: not checked (no population variants supplied)\n")
  }
  invisible(x)
}

#' The Sanger scheme definition (internal, wired into the shared registry)
#'
#' One scheme's declarative row plus its decision-tree vocabulary, in the shape
#' the registry expects: the chemistry's attributes, its selection gates (hard /
#' soft), its ranking tiers, and the short / limit / floor / hard labels the
#' decision graph and text export read.
#' @return a list (`scheme`, `family`, `label`, `discrimination`, `readout`,
#'   `zygosity`, `amplicon`, `loci`, `verdict`, `gates`, `relax`, `tiers`,
#'   `gate_short`, `gate_limit`, `gate_floor_label`, `gate_hard_label`,
#'   `stage_label`).
#' @keywords internal
#' @noRd
#' @method marker_map_spec marker_sanger
#' @export
#' @noRd
marker_map_spec.marker_sanger <- function(marker, scheme = NULL, ...) {
  sp <- map_amplicon(marker)
  start <- sp[["start"]]
  end <- sp[["end"]]
  fwd <- marker$primers$forward
  rev <- marker$primers$reverse
  features <- list(
    spec_row("PCR F", start, start + nchar(fwd) - 1L, "+", "primer"),
    spec_row("PCR R", end - nchar(rev) + 1L, end, "-", "primer")
  )
  sequences <- list("PCR F" = fwd, "PCR R" = rev)
  sp_read <- marker$seq_primers
  if (!is.null(sp_read)) {
    if (!identical(sp_read$forward, fwd)) {
      rs <- sp_read$forward_start %||% start
      re <- sp_read$forward_end %||% (rs + nchar(sp_read$forward) - 1L)
      features <- c(features, list(spec_row("read F", rs, re, "+", "seq")))
      sequences[["read F"]] <- sp_read$forward
    }
    if (!identical(sp_read$reverse, rev)) {
      rs <- sp_read$reverse_start %||% (end - nchar(sp_read$reverse) + 1L)
      re <- sp_read$reverse_end %||% end
      features <- c(features, list(spec_row("read R", rs, re, "-", "seq")))
      sequences[["read R"]] <- sp_read$reverse
    }
  }
  map_spec(marker, features = features, sequences = sequences)
}

scheme_def_sanger <- function() {
  list(
    scheme = "sanger",
    family = "sanger",
    label = "Sanger sequencing (dideoxy read-out)",
    discrimination = "direct base call of the variant",
    readout = "dideoxy chromatogram (capillary)",
    zygosity = "co-dominant (direct read)",
    amplicon = "400-800 bp",
    loci = 1L,
    verdict = "usable Sanger amplicon",
    gates = KASP_GATES_SANGER,
    relax = KASP_RELAX_SANGER,
    tiers = KASP_TIERS_SANGER,
    gate_short = KASP_GATE_SHORT_SANGER,
    gate_limit = KASP_GATE_LIMIT_SANGER,
    gate_floor_label = KASP_GATE_FLOOR_LABEL_SANGER,
    gate_hard_label = KASP_GATE_HARD_LABEL_SANGER,
    stage_label = "Sanger amplicon"
  )
}
