# =============================================================================
#  kasp_rprimer.R  --  an independent oligo-design cross-check via rprimer
# =============================================================================
#  rprimer (Bioconductor) designs degenerate oligos from a *multiple* alignment
#  and reports how each oligo matches the target sequences (on-target perfect /
#  one / two / three+ mismatches, and off-target matches outside the binding
#  region).  It is a completely separate codebase and algorithm from our
#  Biostrings/BLAST path, so it is a genuine second opinion on a locus: "does an
#  independent designer find a clean, well-covered assay here, and how broad is
#  its off-target match?".
#
#  Scope (honest): rprimer validates oligos *it* designed - its S4 classes are
#  locked, so injecting our own primers and re-running checkMatch() is not
#  supported.  This wrapper therefore designs on the alignment and reports the
#  match table; use it to cross-check a locus/alignment, not an existing assay.
# =============================================================================

#' Multiple alignment from a set of sequences (for rprimer)
#'
#' Builds an ungapped `DNAMultipleAlignment` from equal-length sequences (e.g. a
#' REF and an ALT amplicon, or a haplotype set), the input [rprimer_assays()]
#' needs.
#' @param seqs named character vector of equal-length DNA sequences.
#' @return a `Biostrings::DNAMultipleAlignment`.
#' @examples
#' alignment_from(c(REF = "ACGTACGT", ALT = "ACGAACGT"))
#' @family KASP specificity
#' @export
alignment_from <- function(seqs) {
  if (length(unique(nchar(seqs))) != 1L) {
    cli::cli_abort(gettext("`seqs` must be equal length (an ungapped alignment)."))
  }
  Biostrings::DNAMultipleAlignment(Biostrings::DNAStringSet(seqs))
}

#' Independent oligo design + match report on an alignment (rprimer)
#'
#' A thin wrapper over rprimer: builds the consensus profile, designs oligos
#' ([rprimer::designOligos()]) and reports how they match the alignment
#' ([rprimer::checkMatch()]).  Use it as a second opinion about a locus - a small
#' on-target mismatch burden and no off-target match mean the region admits a
#' clean assay, from an independent algorithm.
#' @param alignment a `Biostrings::DNAMultipleAlignment` (see [alignment_from()]).
#' @param check also run `checkMatch()` and attach the match proportions.
#' @param ... passed to [rprimer::designOligos()].
#' @return tidy tibble, one row per designed oligo: `type`, `fwd`, `start`, `end`,
#'   `sequence`, `tm`, `gc` and, when `check`, `perfectMatch`, `oneMismatch`,
#'   `twoMismatches`, `threeMismatches`, `fourOrMoreMismatches` and
#'   `offTargetMatch`.
#' @examples
#' \dontrun{
#' alignment <- alignment_from(c(REF = "...", ALT = "..."))
#' rprimer_assays(alignment)
#' }
#' @family KASP specificity
#' @export
rprimer_assays <- function(alignment, check = TRUE, ...) {
  rlang::check_installed("rprimer", reason = gettext("for an independent oligo-design cross-check"))
  prof <- rprimer::consensusProfile(alignment)
  oligos <- rprimer::designOligos(prof, ...)
  out <- tibble::tibble(
    type = as.character(oligos$type),
    fwd = as.logical(oligos$fwd),
    start = as.integer(oligos$start),
    end = as.integer(oligos$end),
    sequence = purrr::map_chr(oligos$sequence, 1L),
    tm = purrr::map_dbl(oligos$tm, 1L),
    gc = purrr::map_dbl(oligos$gcContent, 1L)
  )
  if (check && nrow(out)) {
    md <- tibble::as_tibble(rprimer::checkMatch(oligos, target = alignment))
    cols <- c(
      "perfectMatch", "oneMismatch", "twoMismatches", "threeMismatches",
      "fourOrMoreMismatches", "offTargetMatch"
    )
    out <- dplyr::bind_cols(out, md[cols])
  }
  out
}
