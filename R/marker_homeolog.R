# =============================================================================
#  marker_homeolog.R  --  the two 3' axes of a KASP candidate, made explicit
# =============================================================================
#  Roadmap Y3.1: a KASP construct has TWO 3'-terminal axes that must never be
#  conflated.  The allele-specific primers' 3' base carries the ALLELE; the
#  COMMON primer's 3' base carries the SUBGENOME (it sits on a homoeolog-
#  diagnostic position).  The diagnostic position cannot be inferred from one
#  genome - it is the PolyMarker mask, produced by comparing subgenomes - so the
#  input is explicit: a caller hands in a mask, and the candidate reports
#  `unique` / `multicopy` / unknown instead of guessing.
# =============================================================================

#' Normalise a subgenome-diagnostic mask to `chrom` / `pos` (internal)
#'
#' A mask is per-position, never per-interval: a `GRanges` row contributes its
#' `start` only (a range does not expand).  A data frame (or list) is preferred
#' in the docs.
#'
#' @param x `NULL`, a data frame / list with `chrom` and `pos`, a `GRanges`, or a
#'   named list `chrom -> integer positions`.
#' @return a tibble `chrom` / `pos` (chrom normalised with `ath_chrom()`), or
#'   `NULL` for a `NULL` input.
#' @keywords internal
#' @noRd
normalise_diagnostic_positions <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "GRanges")) {
    chrom <- as.character(GenomicRanges::seqnames(x))
    pos <- as.integer(GenomicRanges::start(x))
  } else if (is.data.frame(x)) {
    miss <- setdiff(c("chrom", "pos"), names(x))
    if (length(miss)) {
      cli::cli_abort(gettext(
        "`diagnostic_positions` is missing column(s): {paste(miss, collapse = ', ')}."
      ))
    }
    chrom <- as.character(x$chrom)
    pos <- as.integer(x$pos)
  } else if (is.list(x) && all(c("chrom", "pos") %in% names(x))) {
    chrom <- as.character(x$chrom)
    pos <- as.integer(x$pos)
  } else if (is.list(x) && !is.null(names(x))) {
    chrom <- rep(names(x), lengths(x))
    pos <- as.integer(unlist(x, use.names = FALSE))
  } else {
    cli::cli_abort(gettext(paste0(
      "`diagnostic_positions` must be NULL, a data frame / list with `chrom` and `pos`, ",
      "a GRanges, or a named list `chrom -> positions`."
    )))
  }
  out <- tibble::tibble(chrom = ath_chrom(chrom), pos = pos)
  if (anyNA(out$chrom) || anyNA(out$pos) || any(out$pos < 1L)) {
    cli::cli_abort(gettext(
      "`diagnostic_positions` needs a non-NA chromosome and a positive integer position in every row."
    ))
  }
  out
}

#' The two 3' axes of a KASP candidate, made explicit
#'
#' Reads the chosen assay of a [assemble_kasp_assay()] or
#' [assemble_primer_assay()] result and returns the allele axis (the
#' allele-specific primers' plus-strand 3' base) and the subgenome axis (the
#' common primer's plus-strand 3' base and its genomic coordinate), together with
#' a verdict against the diagnostic mask the candidate was built with.
#'
#' @param assay a `marker_assay`.
#' @return a one-row tibble: `allele_base_ref`, `allele_base_alt` (the AS primers'
#'   3'-terminal base, plus strand), `subgenome_anchor` (the COM primer's
#'   3'-terminal base, plus strand), `anchor_pos` (its plus-strand genomic
#'   coordinate), `diagnostic` (`TRUE`/`FALSE`/`NA`), `specificity`
#'   (`"unique"`/`"multicopy"`/`NA`) and `reason` (a short English string for the
#'   `NA`/`FALSE` cases).
#' @examples
#' \dontrun{
#' genome <- Biostrings::DNAStringSet(c(Chr1 = paste(rep("ACGT", 150), collapse = "")))
#' snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "G")
#' assay <- assemble_kasp_assay(genome, snp, pair_screen = 0L, cores = 1L)
#' homeolog_axes(assay)
#' }
#' @family KASP specificity
#' @export
homeolog_axes <- function(assay) {
  a <- assay$assay[1L, ]
  strand <- assay$strand
  allele <- plus_3p_base(c(a$as_ref, as_alt_sequence(assay)), strand)
  com_strand <- if (identical(strand, "+")) "-" else "+"
  subgenome_anchor <- plus_3p_base(a$com, com_strand)
  anchor_pos <- if (identical(strand, "+")) a$com_start else a$com_end

  mask <- assay$diagnostic_positions
  if (is.null(mask)) {
    diagnostic <- NA
    specificity <- NA_character_
    reason <- "no subgenome mask supplied; specificity not tested"
  } else {
    chrom <- ath_chrom(assay$snp$chrom)
    hit <- any(mask$chrom == chrom & mask$pos == anchor_pos)
    diagnostic <- hit
    specificity <- if (hit) "unique" else "multicopy"
    reason <- if (hit) {
      NA_character_
    } else {
      "common primer 3' base is not on a subgenome-diagnostic position"
    }
  }

  tibble::tibble(
    allele_base_ref = allele[[1L]],
    allele_base_alt = allele[[2L]],
    subgenome_anchor = subgenome_anchor,
    anchor_pos = as.integer(anchor_pos),
    diagnostic = diagnostic,
    specificity = specificity,
    reason = reason
  )
}

#' Stamp the two 3' axes onto an assay result and warn when the anchor is off (internal)
#'
#' Reads the mask already stored on `out`, fills `diagnostic_positions`,
#' `allele_base`, `subgenome_anchor` and `homeolog_diagnostics`, and warns when
#' the common primer's 3' base is not on a subgenome-diagnostic position.
#'
#' @param out an assay result under construction.
#' @param mask a normalised diagnostic mask (`NULL` allowed).
#' @return `out`, with the homeolog axes set.
#' @keywords internal
#' @noRd
set_homeolog_axes <- function(out, mask) {
  out$diagnostic_positions <- mask
  dx <- homeolog_axes(out)
  out$allele_base <- c(AS_ref = dx$allele_base_ref, AS_alt = dx$allele_base_alt)
  out$subgenome_anchor <- if (
    !is.na(dx$subgenome_anchor) && !is.na(match(dx$subgenome_anchor, c("A", "C", "G", "T")))
  ) {
    dx$subgenome_anchor
  } else {
    NA_character_
  }
  out$homeolog_diagnostics <- dx
  if (isTRUE(dx$diagnostic == FALSE)) {
    cli::cli_warn(gettext(paste0(
      "The common primer's 3' base is not on a subgenome-diagnostic position: ",
      "the assay genotypes but does not localise."
    )))
  }
  out
}
