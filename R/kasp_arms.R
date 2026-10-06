# =============================================================================
#  kasp_arms.R  --  tetra-primer ARMS-PCR (a gel-readout allele-specific assay)
# =============================================================================
#  Tetra-primer ARMS is the gel analogue of KASP: two outer primers give a
#  control amplicon, and two inner allele-specific primers - placed at the
#  variant, each on the opposite strand - give a REF band and an ALT band of
#  *different sizes*, so a plain agarose gel reads the genotype (AA / AB / BB).
#  Like CAPS it tolerates what KASP cannot (the readout is band presence and
#  size, not an efficiency ratio), but it shares KASP's 3'-terminal placement of
#  the variant.  The inner primers reuse [design_as_primers()]; the outer pair is
#  a plain flanking pair.
# =============================================================================

#' A flanking primer pair with both primers outside a region
#'
#' The generic outer pair: a forward primer ending before `lo` and a reverse
#' primer starting after `hi`, spanning an amplicon in the product window with
#' balanced Tm.  Used by [design_arms_marker()].
#' @keywords internal
#' @noRd
design_flanking_pair <- function(seq, lo, hi, lengths = 18:25,
                                 min_product = 200L, max_product = 800L,
                                 tm_range = KASP_LIMITS$com_tm, gc_range = KASP_LIMITS$gc) {
  s <- as.character(seq)
  n <- nchar(s)
  if (lo <= 1L || hi >= n) {
    return(NULL)
  }
  ok <- function(p) !has_run(p) & !has_repeat(p)
  screen <- function(d) {
    d$tm <- primer_tm(d$primer)
    d$gc <- gc_frac(d$primer)
    dplyr::filter(
      d, dplyr::between(.data$tm, tm_range[1L], tm_range[2L]),
      dplyr::between(.data$gc, gc_range[1L], gc_range[2L])
    )
  }
  fwd <- tidyr::crossing(p1 = seq_len(max(1L, lo - 1L)), len = lengths) |>
    dplyr::filter(.data$p1 + .data$len - 1L < lo) |>
    dplyr::mutate(primer = stringr::str_sub(s, .data$p1, .data$p1 + .data$len - 1L)) |>
    dplyr::filter(ok(.data$primer)) |>
    screen()
  rev <- tidyr::crossing(p1 = seq.int(hi + 1L, n), len = lengths) |>
    dplyr::filter(.data$p1 + .data$len - 1L <= n) |>
    dplyr::mutate(primer = as.character(Biostrings::reverseComplement(
      Biostrings::DNAStringSet(stringr::str_sub(s, .data$p1, .data$p1 + .data$len - 1L))
    ))) |>
    dplyr::filter(ok(.data$primer)) |>
    screen()
  if (!nrow(fwd) || !nrow(rev)) {
    return(NULL)
  }
  pairs <- dplyr::cross_join(fwd, rev, suffix = c("_f", "_r")) |>
    dplyr::mutate(
      product = .data$p1_r + .data$len_r - 1L - .data$p1_f + 1L,
      dtm = abs(.data$tm_f - .data$tm_r)
    ) |>
    dplyr::filter(dplyr::between(.data$product, min_product, max_product))
  if (!nrow(pairs)) {
    return(NULL)
  }
  b <- pairs[which.min(pairs$dtm), ]
  list(
    forward = b$primer_f, reverse = b$primer_r,
    start = b$p1_f, end = b$p1_r + b$len_r - 1L, product = b$product,
    tm_f = b$tm_f, tm_r = b$tm_r
  )
}

#' Design a tetra-primer ARMS-PCR marker at a variant
#'
#' The gel-readout allele-specific assay: an outer control pair plus two inner
#' allele-specific primers (one per allele, on opposite strands) whose products
#' differ in size, so AA / AB / BB resolve on an agarose gel.  The inner primers
#' reuse [design_as_primers()] (the same 3'-terminal placement, plus its optional
#' stabilising mismatch), the outer pair is a plain flanking pair.  Single-base
#' substitutions only.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank window half-width for the outer pair.
#' @param min_band_gap minimal size difference between the two diagnostic bands.
#' @param ... passed to [design_as_primers()] (e.g. `stab_ats`, `BPPARAM`).
#' @return a `kasp_arms` list (`snp`, `outer`, `inner`, `bands`, `amplicon`), or
#'   `NULL` (with a message) when a pair cannot be built.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_arms_marker(genome, kasp_example("SNP")[1, ], flank = 300)
#' }
#' @family KASP ARMS
#' @export
design_arms_marker <- function(genome, snp, flank = 300L, min_band_gap = 20L, ...) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    rlang::inform(gettext("tetra-primer ARMS is implemented for single-base substitutions only."))
    return(NULL)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  outer <- design_flanking_pair(tmpl$seq, tmpl$snp_i, tmpl$snp_i)
  if (is.null(outer)) {
    rlang::inform(gettext("No outer primer pair passed the Tm/GC filters."))
    return(NULL)
  }
  inner <- design_as_primers(tmpl$seq, tmpl$snp_i, ref, alt,
    strands = c("+", "-"), gc_range = c(0.35, 0.70),
    tm_range = KASP_LIMITS$as_tm, ...
  )
  if (!nrow(inner)) {
    rlang::inform(gettext("No inner allele-specific primer passed the GC/Tm filter."))
    return(NULL)
  }
  best <- function(st) {
    d <- inner[inner$strand == st, , drop = FALSE]
    if (!nrow(d)) {
      return(NULL)
    }
    d[which.max(d$discrimination), ]
  }
  fwd <- best("+") # REF inner: forward, 3' base = REF
  rev <- best("-") # ALT inner: reverse, 3' base = ALT
  if (is.null(fwd) || is.null(rev)) {
    rlang::inform(gettext("An inner primer is missing on one strand."))
    return(NULL)
  }
  # REF band: inner-forward .. outer-reverse; ALT band: outer-forward .. inner-reverse
  ref_size <- outer$end - fwd$tmpl_start + 1L
  alt_size <- rev$tmpl_end - outer$start + 1L
  if (abs(ref_size - alt_size) < min_band_gap) {
    rlang::inform(gettext("The two diagnostic bands are closer than `min_band_gap`."))
    return(NULL)
  }
  off <- snp$pos - tmpl$snp_i
  out <- list(
    snp = snp,
    outer = list(
      forward = outer$forward, reverse = outer$reverse,
      start = outer$start + off, end = outer$end + off, product = outer$product
    ),
    inner = list(
      ref = list(primer = fwd$ref_primer, start = fwd$tmpl_start + off, end = fwd$tmpl_end + off, strand = "+"),
      alt = list(primer = rev$alt_primer, start = rev$tmpl_start + off, end = rev$tmpl_end + off, strand = "-")
    ),
    bands = tibble::tibble(
      band = c("control", "REF", "ALT"),
      size = c(outer$product, ref_size, alt_size)
    ),
    amplicon = tibble::tibble(start = outer$start + off, end = outer$end + off)
  )
  structure(out, class = "kasp_arms")
}

#' @export
print.kasp_arms <- function(x, ...) {
  cat("<kasp_arms>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt, "\n", sep = "")
  cat("  outer:", x$outer$product, "bp\n")
  cat("  REF band:", x$bands$size[x$bands$band == "REF"], "bp (inner F:", x$inner$ref$primer, ")\n")
  cat("  ALT band:", x$bands$size[x$bands$band == "ALT"], "bp (inner R:", x$inner$alt$primer, ")\n")
  invisible(x)
}

#' Draw the expected tetra-primer ARMS gel
#'
#' Two lanes (REF and ALT homozygotes) with the control band plus the
#' allele-specific diagnostic band - the AA / AB / BB readout.
#' @param marker a `kasp_arms` object (see [design_arms_marker()]).
#' @return a `ggplot`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' plot_arms_gel(design_arms_marker(genome, kasp_example("SNP")[1, ], flank = 300))
#' }
#' @family KASP ARMS
#' @export
plot_arms_gel <- function(marker) {
  rlang::check_installed("ggplot2", reason = gettext("to draw the ARMS gel"))
  if (!inherits(marker, "kasp_arms")) {
    rlang::abort(gettext("`marker` must be a kasp_arms object (see design_arms_marker())."))
  }
  ctrl <- marker$bands$size[marker$bands$band == "control"]
  ref <- marker$bands$size[marker$bands$band == "REF"]
  alt <- marker$bands$size[marker$bands$band == "ALT"]
  bands <- tibble::tibble(
    allele = c("REF", "REF", "ALT", "ALT"),
    size = c(ctrl, ref, ctrl, alt),
    diagnostic = c(FALSE, TRUE, FALSE, TRUE)
  )
  bands$x <- as.numeric(factor(bands$allele, levels = c("REF", "ALT")))
  bands$col <- ifelse(bands$diagnostic, ifelse(bands$allele == "REF", "#4477aa", "#d7301f"), "grey65")
  ggplot2::ggplot(bands) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$x - 0.32, xend = .data$x + 0.32, y = .data$size, yend = .data$size),
      linewidth = 2.4, lineend = "round", colour = bands$col
    ) +
    ggplot2::scale_x_continuous(breaks = 1:2, labels = c("REF/REF", "ALT/ALT"), limits = c(0.4, 2.6)) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      x = NULL, y = "band size (bp)",
      title = glue::glue("ARMS gel: {marker$snp$chrom}:{marker$snp$pos}"),
      subtitle = glue::glue("control {ctrl} bp | REF {ref} bp | ALT {alt} bp")
    ) +
    ggplot2::theme_minimal(base_size = 10)
}
