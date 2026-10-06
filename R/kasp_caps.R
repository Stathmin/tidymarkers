# =============================================================================
#  kasp_caps.R  --  CAPS / dCAPS markers (the PCR family beyond KASP)
# =============================================================================
#  Where KASP needs the variant at an allele-specific primer's 3' end, a CAPS
#  marker anchors the variant *inside the amplicon* and reads it with a
#  restriction digest: a biallelic site whose recognition sequence is created or
#  destroyed by the variant (or introduced by a primer, dCAPS) yields a
#  different set of fragment lengths for each allele - a binary, gel-scorable
#  call.  That is why CAPS tolerates the very things that doom a KASP: a variant
#  in the 3' window, an AT-rich junction, or marginal discrimination (see
#  ROADMAP section U).  The restriction-reaction recipe (enzyme + expected
#  fragments) is stored on the object, like the KASP mix recipe is.
# =============================================================================

#' Bundled restriction-enzyme catalogue (name, site, cut)
#'
#' A small, curated table of common, **palindromic** restriction enzymes: the
#' recognised site (IUPAC) and the top-strand cut position (bases from the 5' end
#' of the site, 0 = before the first base).  Palindromic sites are read on the
#' plus strand alone, so one match per site is enough.  The catalogue is data
#' (`inst/extdata/restriction_enzymes.csv`), not code, so it can be extended
#' without touching the engine.
#' @return tibble with `enzyme`, `site`, `cut`.
#' @examples
#' head(restriction_enzymes())
#' @family KASP CAPS
#' @export
restriction_enzymes <- function() {
  f <- system.file("extdata", "restriction_enzymes.csv", package = "tidykasp")
  if (!nzchar(f)) rlang::abort(gettext("Bundled restriction-enzyme catalogue not found."))
  readr::read_csv(f, show_col_types = FALSE)
}

#' Digest fragment lengths of one sequence, per enzyme
#'
#' For each enzyme, finds the palindromic recognition sites (IUPAC-aware) and
#' returns the top-strand cut boundaries and the resulting linear fragment
#' lengths.  An enzyme with no site returns a single band (no cut).
#' @param seq a DNAString or character sequence (the amplicon, plus strand).
#' @param enzymes enzyme table (see [restriction_enzymes()]).
#' @return tibble with `enzyme`, `n_sites`, `cuts` (list of 0-based boundaries)
#'   and `fragments` (list of integer band lengths).
#' @examples
#' caps_fragments("GAATTCTTAAGGAATTC")
#' @family KASP CAPS
#' @export
caps_fragments <- function(seq, enzymes = restriction_enzymes()) {
  s <- Biostrings::DNAString(as.character(seq))
  n <- length(s)
  purrr::pmap(enzymes, function(enzyme, site, cut) {
    hits <- Biostrings::matchPattern(site, s, fixed = FALSE)
    starts <- IRanges::start(hits) - 1L
    cut_at <- sort(unique(starts + as.integer(cut)))
    cut_at <- cut_at[cut_at > 0L & cut_at < n]
    frag <- as.integer(diff(c(0L, cut_at, n)))
    tibble::tibble(
      enzyme = enzyme, n_sites = length(starts),
      cuts = list(cut_at), fragments = list(frag)
    )
  }) |>
    purrr::list_rbind()
}

#' The fragment-length gap between two digests (how resolvable on a gel)
#'
#' The smallest absolute difference between a REF band and an ALT band (ignoring
#' equal pairs) - the closest two bands that a gel would have to separate.  `Inf`
#' when the two digest patterns are identical.
#' @keywords internal
#' @noRd
digest_gap <- function(a, b) {
  if (length(a) == length(b) && setequal(a, b)) {
    return(Inf)
  }
  d <- abs(outer(a, b, "-"))
  d <- d[d > 0]
  if (length(d)) min(d) else Inf
}

#' Enzymes whose digest distinguishes the REF and ALT amplicons (CAPS scan)
#'
#' Copies the target into a window, builds the two allele amplicons (REF vs ALT),
#' digests each with every enzyme and keeps the enzymes whose fragment-length
#' pattern differs - i.e. a natural CAPS marker.  A pattern is reported as
#' *resolvable* when the closest differing bands are at least `min_gap` apart
#' (`min_gap = 0` keeps everything that differs at all).
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank window half-width (the prospective amplicon).
#' @param enzymes enzyme table (see [restriction_enzymes()]).
#' @param min_gap minimal band separation to call the digest gel-resolvable.
#' @return tidy tibble, one row per enzyme with `fragments_ref`, `fragments_alt`,
#'   `gap` and `resolvable`, sorted by decreasing `gap`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' snp <- kasp_example("SNP")[1, ]
#' find_caps_enzymes(genome, snp, flank = 250)
#' }
#' @family KASP CAPS
#' @export
find_caps_enzymes <- function(genome, snp, flank = 250L, enzymes = restriction_enzymes(),
                              min_gap = 50L) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  ref_amp <- tmpl$seq
  ref_len <- nchar(as.character(snp$ref)[1L])
  alt_amp <- Biostrings::replaceAt(
    ref_amp, IRanges::IRanges(tmpl$snp_i, width = ref_len), as.character(snp$alt)[1L]
  )
  fr <- caps_fragments(ref_amp, enzymes)
  fa <- caps_fragments(alt_amp, enzymes)
  out <- dplyr::inner_join(
    dplyr::rename(fr, n_sites_ref = "n_sites", cuts_ref = "cuts", fragments_ref = "fragments"),
    dplyr::rename(fa, n_sites_alt = "n_sites", cuts_alt = "cuts", fragments_alt = "fragments"),
    by = "enzyme"
  )
  out$gap <- purrr::map2_dbl(out$fragments_ref, out$fragments_alt, digest_gap)
  out$resolvable <- is.finite(out$gap) & out$gap >= min_gap
  out <- dplyr::filter(out, is.finite(.data$gap))
  out[order(-out$gap), ]
}

#' Flanking primer pair around a restricted region
#'
#' A compact, GC/Tm-screened scan: forward primers from the left of `[lo, hi]`,
#' reverse primers to the right, both clear of the site, spanning an amplicon in
#' `[min_product, max_product]` with balanced Tm.  Reuses the package's Tm/GC/run
#' primitives; no restriction logic here.
#' @keywords internal
#' @noRd
design_caps_primers <- function(seq, lo, hi, target_tm = AS_TM_TARGET,
                                lengths = 18:25, min_product = 120L, max_product = 1000L,
                                tm_range = KASP_LIMITS$as_tm, gc_range = KASP_LIMITS$gc) {
  s <- as.character(seq)
  n <- nchar(s)
  if (lo <= 1L || hi >= n) {
    return(NULL)
  }
  fwd <- tidyr::crossing(start = seq_len(max(1L, lo - 1L)), length = lengths) |>
    dplyr::filter(.data$start + .data$length - 1L < lo) |>
    dplyr::mutate(primer = stringr::str_sub(s, .data$start, .data$start + .data$length - 1L)) |>
    dplyr::filter(!has_run(.data$primer), !has_repeat(.data$primer))
  rev <- tidyr::crossing(end = seq.int(hi + 1L, n), length = lengths) |>
    dplyr::filter(.data$end - .data$length + 1L > hi) |>
    dplyr::mutate(primer = as.character(Biostrings::reverseComplement(
      Biostrings::DNAStringSet(stringr::str_sub(s, .data$end - .data$length + 1L, .data$end))
    ))) |>
    dplyr::filter(!has_run(.data$primer), !has_repeat(.data$primer))
  if (!nrow(fwd) || !nrow(rev)) {
    return(NULL)
  }
  fwd$tm <- primer_tm(fwd$primer)
  fwd$gc <- gc_frac(fwd$primer)
  rev$tm <- primer_tm(rev$primer)
  rev$gc <- gc_frac(rev$primer)
  fwd <- dplyr::filter(
    fwd, dplyr::between(.data$tm, tm_range[1L], tm_range[2L]),
    dplyr::between(.data$gc, gc_range[1L], gc_range[2L])
  )
  rev <- dplyr::filter(
    rev, dplyr::between(.data$tm, tm_range[1L], tm_range[2L]),
    dplyr::between(.data$gc, gc_range[1L], gc_range[2L])
  )
  if (!nrow(fwd) || !nrow(rev)) {
    return(NULL)
  }
  pairs <- dplyr::cross_join(fwd, rev, suffix = c("_f", "_r")) |>
    dplyr::mutate(product = .data$end - .data$start + 1L, dtm = abs(.data$tm_f - .data$tm_r)) |>
    dplyr::filter(dplyr::between(.data$product, min_product, max_product)) |>
    dplyr::arrange(.data$dtm, .data$product)
  if (!nrow(pairs)) {
    return(NULL)
  }
  best <- pairs[1L, ]
  list(
    forward = best$primer_f, reverse = best$primer_r,
    start = best$start, end = best$end,
    product = best$product, tm_f = best$tm_f, tm_r = best$tm_r
  )
}

#' Design a CAPS / dCAPS marker at a variant
#'
#' The CAPS analogue of [assemble_kasp_assay()]: find the enzyme whose digest
#' separates the alleles (see [find_caps_enzymes()]), then design the flanking
#' primer pair around the discriminating site.  Returns a `kasp_caps` object with
#' the enzyme, both digest patterns and the primer pair - the restriction recipe
#' travels with the marker, as the KASP mix recipe does.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank window half-width for the scan and primer design.
#' @param enzymes enzyme table (see [restriction_enzymes()]).
#' @param min_gap minimal band separation to require (see [find_caps_enzymes()]).
#' @return a `kasp_caps` list (`snp`, `enzyme`, `fragments_ref`, `fragments_alt`,
#'   `gap`, `primers`, `amplicon`), or `NULL` (with a message) when no enzyme
#'   separates the alleles.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_caps_marker(genome, kasp_example("SNP")[1, ], flank = 250)
#' }
#' @family KASP CAPS
#' @export
design_caps_marker <- function(genome, snp, flank = 250L, enzymes = restriction_enzymes(),
                               min_gap = 50L) {
  hits <- find_caps_enzymes(genome, snp, flank, enzymes, min_gap)
  hits <- dplyr::filter(hits, .data$resolvable)
  if (!nrow(hits)) {
    rlang::inform(gettext("No enzyme separates the REF and ALT amplicons - no natural CAPS marker here."))
    return(NULL)
  }
  best <- hits[1L, ]
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  # the discriminating cut(s): where the two digest patterns first diverge
  cuts <- setdiff(best$cuts_ref[[1L]], best$cuts_alt[[1L]])
  if (!length(cuts)) cuts <- setdiff(best$cuts_alt[[1L]], best$cuts_ref[[1L]])
  lo <- min(cuts) + 1L
  hi <- max(cuts) + 1L
  primers <- design_caps_primers(tmpl$seq, lo, hi)
  if (is.null(primers)) {
    rlang::inform(gettext("A discriminating enzyme was found, but no flanking primer pair passed the Tm/GC filters."))
    return(NULL)
  }
  out <- list(
    snp = snp, enzyme = best$enzyme,
    fragments_ref = best$fragments_ref[[1L]], fragments_alt = best$fragments_alt[[1L]],
    gap = best$gap, primers = primers,
    amplicon = tibble::tibble(
      start = snp$pos + primers$start - tmpl$snp_i,
      end = snp$pos + primers$end - tmpl$snp_i
    )
  )
  structure(out, class = "kasp_caps")
}

#' @export
print.kasp_caps <- function(x, ...) {
  cat("<kasp_caps>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, " via ", x$enzyme, "\n",
    sep = ""
  )
  cat("  REF digest:", paste(x$fragments_ref, collapse = " + "), "bp\n")
  cat("  ALT digest:", paste(x$fragments_alt, collapse = " + "), "bp\n")
  cat("  amplicon:", x$primers$product, "bp (gap", x$gap, "bp)\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  invisible(x)
}
