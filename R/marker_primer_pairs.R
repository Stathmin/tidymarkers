# =============================================================================
#  marker_primer_pairs.R  --  generic PCR primer-pair primitives (marker-agnostic)
# =============================================================================
#  Every marker family here (KASP, CAPS/dCAPS, tetra-primer ARMS) ends up
#  scanning for primers that pass the same composition filters and pairing them
#  across a product window.  That scan used to be copy-pasted three times; this
#  is the one place it lives.  It is deliberately marker-agnostic: no KASP
#  3'-anchor, no restriction logic - just "candidate primers of these lengths at
#  these positions, GC/Tm-screened, paired by Tm balance".
# =============================================================================

#' A GC/Tm-screened grid of candidate primers along one strand
#'
#' Builds every (start, length) primer on `seq`, takes the reverse complement for
#' the `"-"` strand, and keeps those free of homopolymer runs / SSRs and inside
#' the GC and Tm windows.  `start`/`end` are always plus-strand coordinates, so a
#' reverse primer's `start` is still left of its `end`.
#' @param seq the template (plus strand).
#' @param starts 1-based plus-strand start positions to try.
#' @param lengths primer lengths to try.
#' @param strand `"+"` (primer equals the plus strand) or `"-"` (reverse
#'   complement).
#' @param tm_range,gc_range accepted windows.
#' @param conditions reaction conditions for [primer_tm()] (`NULL` = the
#'   historical 50 mM Na / 37 C buffer).  Only the buffer fields of `conditions`
#'   enter the nearest-neighbour Tm.
#' @return tibble with `start`, `len`, `end`, `primer`, `tm`, `gc`.
#' @keywords internal
#' @noRd
primer_candidates <- function(seq, starts, lengths, strand = c("+", "-"), tm_range, gc_range,
                              conditions = NULL) {
  strand <- match.arg(strand)
  s <- as.character(seq)
  n <- nchar(s)
  d <- tidyr::crossing(start = starts, len = lengths) |>
    dplyr::filter(.data$start >= 1L, .data$start + .data$len - 1L <= n) |>
    dplyr::mutate(
      primer = stringr::str_sub(s, .data$start, .data$start + .data$len - 1L),
      end = .data$start + .data$len - 1L
    )
  if (identical(strand, "-")) {
    d$primer <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(d$primer)))
  }
  d <- dplyr::filter(d, !has_run(.data$primer), !has_repeat(.data$primer))
  d$tm <- primer_tm(d$primer, conditions = conditions)
  d$gc <- gc_frac(d$primer)
  dplyr::filter(
    d, dplyr::between(.data$tm, tm_range[1L], tm_range[2L]),
    dplyr::between(.data$gc, gc_range[1L], gc_range[2L])
  )
}

#' A GC/Tm-screened grid of plus-strand probes spanning a position
#'
#' The shared enumeration behind the hydrolysis-probe schemes (TaqMan, castPCR):
#' every (start, length) plus-strand window whose span covers `snp_i`, GC/Tm/
#' composition-screened by `primer_candidates()`.  Callers add their own extra
#' filters (e.g. no 5'-G) and ranking.
#' @param seq the template (plus strand).
#' @param snp_i 1-based variant position on `seq`.
#' @param lengths probe lengths to try.
#' @param tm_range,gc_range accepted windows.
#' @param conditions reaction conditions for [primer_tm()] (see
#'   `primer_candidates()`; `NULL` keeps the historical buffer).
#' @return a tibble (see `primer_candidates()`), possibly zero rows.
#' @keywords internal
#' @noRd
spanning_probe_candidates <- function(seq, snp_i, lengths, tm_range, gc_range,
                                      conditions = NULL) {
  n <- nchar(as.character(seq))
  starts <- seq.int(max(1L, snp_i - max(lengths) + 1L), min(n, snp_i))
  primer_candidates(seq, starts, lengths, "+", tm_range, gc_range,
    conditions = conditions
  ) |>
    dplyr::filter(.data$start <= snp_i, .data$end >= snp_i)
}

#' The best Tm-balanced forward x reverse pair in a product window
#'
#' `fwd` must be a `"+"` candidate set, `rev` a `"-"` set (plus-strand
#' coordinates).  Returns the one with the closest Tm within the product window,
#' or `NULL`.
#' @keywords internal
#' @noRd
pick_primer_pair <- function(fwd, rev, min_product, max_product) {
  if (!nrow(fwd) || !nrow(rev)) {
    return(NULL)
  }
  pairs <- dplyr::cross_join(fwd, rev, suffix = c("_f", "_r")) |>
    dplyr::mutate(
      product = .data$end_r - .data$start_f + 1L,
      dtm = abs(.data$tm_f - .data$tm_r)
    ) |>
    dplyr::filter(dplyr::between(.data$product, min_product, max_product))
  if (!nrow(pairs)) {
    return(NULL)
  }
  pairs[order(pairs$dtm, pairs$product), ][1L, ]
}

#' One pair row as the marker-facing list
#' @keywords internal
#' @noRd
pair_to_list <- function(b) {
  list(
    forward = b$primer_f, reverse = b$primer_r,
    start = b$start_f, end = b$end_r, product = b$product,
    tm_f = b$tm_f, tm_r = b$tm_r
  )
}

#' A flanking primer pair with both primers outside a region
#'
#' The generic outer pair (used by [design_arms_marker()]): a forward primer
#' ending before `lo` and a reverse primer starting after `hi`.
#' @keywords internal
#' @noRd
design_flanking_pair <- function(seq, lo, hi, lengths = 18:25,
                                 min_product = 200L, max_product = 800L,
                                 tm_range = KASP_LIMITS$com_tm, gc_range = KASP_LIMITS$gc) {
  n <- nchar(as.character(seq))
  if (lo <= 1L || hi >= n) {
    return(NULL)
  }
  fwd <- primer_candidates(seq, seq_len(lo - 1L), lengths, "+", tm_range, gc_range) |>
    dplyr::filter(.data$end < lo)
  rev <- primer_candidates(seq, seq.int(hi + 1L, n), lengths, "-", tm_range, gc_range)
  b <- pick_primer_pair(dplyr::filter(fwd, .data$end < hi), rev, min_product, max_product)
  if (is.null(b)) {
    return(NULL)
  }
  pair_to_list(b)
}

#' A primer pair around a short restricted region, both orientations
#'
#' One primer of the pair *covers* `[lo, hi]` (the site a primer introduces or the
#' cut an enzyme makes) and extends to a flanking mate on the other side; both
#' orientations (forward/covering, reverse/covering) are tried and the best
#' Tm-balanced valid pair is returned.  Used by CAPS/dCAPS.
#' @keywords internal
#' @noRd
design_site_pair <- function(seq, lo, hi, site_lengths = 18:28, mate_lengths = 18:25,
                             min_product = 120L, max_product = 1200L,
                             tm_range = KASP_LIMITS$as_tm, gc_range = KASP_LIMITS$gc) {
  n <- nchar(as.character(seq))
  if (lo <= 1L || hi >= n) {
    return(NULL)
  }
  s_starts <- seq.int(max(1L, lo - max(site_lengths) + 1L), lo)
  cover_f <- primer_candidates(seq, s_starts, site_lengths, "+", tm_range, gc_range) |>
    dplyr::filter(.data$end >= hi)
  cover_r <- primer_candidates(seq, s_starts, site_lengths, "-", tm_range, gc_range) |>
    dplyr::filter(.data$end >= hi)
  mate_r <- primer_candidates(seq, seq.int(hi + 1L, n), mate_lengths, "-", tm_range, gc_range)
  mate_f <- primer_candidates(seq, seq_len(lo - 1L), mate_lengths, "+", tm_range, gc_range) |>
    dplyr::filter(.data$end < lo)
  p1 <- pick_primer_pair(cover_f, mate_r, min_product, max_product)
  p2 <- pick_primer_pair(mate_f, cover_r, min_product, max_product)
  cand <- dplyr::bind_rows(p1, p2)
  if (!nrow(cand)) {
    return(NULL)
  }
  pair_to_list(cand[which.min(abs(cand$tm_f - cand$tm_r)), ])
}
