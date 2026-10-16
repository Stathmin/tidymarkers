# =============================================================================
#  kasp_sanger.R  --  Sanger sequencing (dideoxy read-out) as a marker scheme
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
#  The read-out is not simulated (no chromatogram model), so this returns a
#  well-formed, unique amplicon with the variant correctly placed, plus the
#  sequencing primers that read into it, not a trace.
#
#  Selection is real (not decorative): every gate declared by scheme_def_sanger()
#  is computed here and applied by select_with_ladder() (C7).  The cheap margins
#  (variant distance, amplicon window, Tm, 3'-end composition, 3'-seed
#  uniqueness) are computed for every candidate pair; the expensive ones (primer
#  hairpin/homodimer, the full in-silico PCR) are confirmed on the shortlist, so
#  a design never pays for an OligoArrayAux call per candidate.
# =============================================================================

# The Sanger acceptance windows (the constraints the gates below encode).  Tm is
# the sequencing-core 55-65 C rather than the tighter KASP COM window; min_dist
# is the low-quality lead of a read.
KASP_SANGER_LIMITS <- list(
  amplicon = c(400L, 800L),
  tm = c(55, 65),
  gc = c(0.40, 0.60),
  min_dist = 50L
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
#' The bound that keeps [flanking_pair_grid()] affordable on a long flank:
#' positions are preserved (so the placement tiers still see the whole flank) and
#' for each only the length whose Tm is closest to `target` is kept.
#' @param cand a candidate primer tibble (see [primer_candidates()]).
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
#' is an indexed `kasp_genome` (`specificity_seed()` needs an in-memory
#' `DNAStringSet`), which the ladder reads as "unknown" and lets through to the
#' full in-silico PCR.
#' @param primer_f,primer_r the pair's oligos.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param k 3'-seed length.
#' @return a numeric vector, one per pair.
#' @keywords internal
#' @noRd
seed_margin <- function(primer_f, primer_r, genome, k = SEED_K) {
  if (inherits(genome, "kasp_genome")) {
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
#' Phase 1 is the shared [select_with_ladder()] over the *cheap* margins, so the
#' declared gate order, relaxation order and floors are the ones the DAG shows.
#' Phase 2 walks the tier-ordered survivors and applies the *expensive* gate(s)
#' via `confirm()` (OligoArrayAux structure, in-silico PCR): the first candidate
#' the confirmation accepts wins, so a design never runs the expensive test on
#' more than `max_alternatives` pairs, and a candidate the confirmation rejects is
#' dropped and the ladder re-run (which can also relax a soft gate further).
#' @param grid the candidate grid (see [flanking_pair_grid()]).
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
#' >= `min_dist` from the variant, Tm-screened and free of a 3' homopolymer),
#' preferring a nested primer and falling back to the flanking pair (the routine
#' workflow).  A nested primer is kept only when its 3' seed is unique, so the
#' read is not a mixture.  The nested search window is deliberately narrow (<= 40
#' bp inside the PCR primer), so the design stays a handful of Tm calls.
#' @param seq the template (plus strand); `si` the variant's template position.
#' @param start_f,end_r plus-strand coordinates of the flanking pair's ends.
#' @param flank_f,flank_r the flanking PCR primers (the fallback sequencing ones).
#' @param min_dist minimal read-start distance.
#' @param genome a `kasp_genome` or `DNAStringSet` (for the seed screen).
#' @return a list (`forward`, `reverse`, `nested`).
#' @keywords internal
#' @noRd
sanger_seq_primers <- function(seq, si, start_f, end_r, flank_f, flank_r, min_dist, genome,
                               lengths, tm_range, gc_range) {
  nested_lens <- intersect(c(20L, 22L), lengths)
  if (!length(nested_lens)) {
    nested_lens <- min(lengths)
  }
  window <- 40L
  fwd <- seq_primer_at(
    seq, "+", start_f + seq.int(0L, window), nested_lens, tm_range, gc_range,
    function(cand) cand$end <= si - min_dist
  )
  rev <- seq_primer_at(
    seq, "-", seq.int(end_r - window - max(nested_lens) + 1L, end_r - min(nested_lens) + 1L),
    nested_lens, tm_range, gc_range,
    function(cand) cand$start >= si + min_dist & cand$end <= end_r
  )
  if (!is.null(fwd) && !is.null(rev) && !inherits(genome, "kasp_genome")) {
    p <- c(fwd$primer, rev$primer)
    if (any(seed_margin(p, p, genome) < 0)) {
      fwd <- rev <- NULL
    }
  }
  list(
    forward = fwd$primer %||% flank_f,
    reverse = rev$primer %||% flank_r,
    nested = !is.null(fwd) && !is.null(rev)
  )
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
#' unique amplicon, not a predicted trace.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width (must leave room for the amplicon).
#' @param min_dist minimal distance from the variant to each amplicon end (the
#'   read start / primer binding site).
#' @param min_product,max_product the accepted amplicon window.
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows.
#' @param check_unique when `TRUE`, require exactly one [in_silico_pcr()] product.
#' @param ... passed to [in_silico_pcr()] (e.g. `BPPARAM`).
#' @return a `kasp_sanger` list (`snp`, `primers`, `amplicon`, `seq_primers`,
#'   `selection`), or `NULL` (with a message) when no unique pair passes the
#'   gates.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_sanger_marker(genome, kasp_example("SNP")[1, ], flank = 500)
#' }
#' @family KASP Sanger
#' @export
design_sanger_marker <- function(genome, snp, flank = 500L,
                                 min_dist = KASP_SANGER_LIMITS$min_dist,
                                 min_product = KASP_SANGER_LIMITS$amplicon[1L],
                                 max_product = KASP_SANGER_LIMITS$amplicon[2L],
                                 lengths = 18:25, tm_range = KASP_SANGER_LIMITS$tm,
                                 gc_range = KASP_SANGER_LIMITS$gc,
                                 check_unique = TRUE, ...) {
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
    primers$forward, primers$reverse, min_dist, genome, lengths, tm_range, gc_range
  )
  amp0 <- snp$pos - si
  out <- list(
    snp = snp, primers = primers,
    amplicon = tibble::tibble(start = primers$start + amp0, end = primers$end + amp0),
    seq_primers = seq_primers,
    selection = list(relaxed = sel$relaxed, violated = sel$violated)
  )
  structure(out, class = "kasp_sanger")
}

#' @export
print.kasp_sanger <- function(x, ...) {
  cat(
    "<kasp_sanger>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  amplicon:", x$primers$product, "bp (Sanger read)\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  if (!is.null(x$seq_primers)) {
    how <- if (isTRUE(x$seq_primers$nested)) "nested" else "flanking"
    cat("  sequencing primers (", how, "):\n", sep = "")
    cat("    F:", x$seq_primers$forward, "\n    R:", x$seq_primers$reverse, "\n")
  }
  if (length(x$selection$relaxed)) {
    cat("  relaxed:", paste(x$selection$relaxed, collapse = ", "), "\n")
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
