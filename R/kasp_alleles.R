# =============================================================================
#  kasp_alleles.R  --  allele sets and multi-allelic KASP sites
# =============================================================================
#  A KASP assay is normally biallelic: one allele-specific primer per allele.
#  A site with more than one ALT allele is therefore not automatically
#  un-designable - it falls into one of a few well-defined classes (see
#  `kasp_site_class()`), and only the single complementary ALT is truly hopeless.
#  The IUPAC consensus comes from `Biostrings::IUPAC_CODE_MAP`, so no ambiguity
#  table is reimplemented here.
# =============================================================================

#' IUPAC consensus code for a set of bases
#'
#' The single-letter IUPAC code whose members are exactly `bases`.  The lookup
#' is the inverse of `Biostrings::IUPAC_CODE_MAP` (`"AG" -> "R"`, `"CGT" -> "B"`,
#' ...), so the ambiguity alphabet is not redefined here.
#' @param bases character vector of single bases (only `A`, `C`, `G`, `T`).
#' @return a one-letter code, or `NA` when the set is empty or has no exact code
#'   (e.g. it contains `N`).
#' @keywords internal
#' @noRd
iupac_code <- function(bases) {
  bases <- sort(unique(as.character(bases)))
  bases <- bases[!is.na(bases) & nzchar(bases)]
  if (!length(bases) || !all(bases %in% c("A", "C", "G", "T"))) {
    return(NA_character_)
  }
  key <- paste(bases, collapse = "")
  hit <- names(Biostrings::IUPAC_CODE_MAP)[Biostrings::IUPAC_CODE_MAP == key]
  if (length(hit)) hit[[1L]] else NA_character_
}

#' Classify a site's allele set for KASP design
#'
#' Decides whether a variant site is designable despite its ALT roster, so a
#' multi-allelic site is not discarded just because a single ALT cannot be told
#' from the reference.  The classes are:
#' * `"biallelic"` — a single usable ALT whose base is not complementary to REF;
#' * `"degenerate"` — two or more distinct SNP ALTs: one AS primer whose
#'   3'-terminal base is the IUPAC consensus [`iupac_code()`] of the ALT set
#'   matches all of them and mismatches REF (the degenerate-primer option);
#' * `"complementary"` — a single ALT that is the Watson-Crick complement of REF
#'   (A/T or C/G), which a 3'-terminal base cannot separate;
#' * `"indel"` — only non-SNP ALTs (handled by the indel-junction path);
#' * `"none"` — no usable ALT.
#'
#' At a single site the ALTs differ only in the variant base, so a degenerate
#' code over them always excludes REF; the criterion therefore only bites when
#' the ALT set is not encodable or when an ALT is not a SNP.
#' @param ref REF allele (a single base).
#' @param alts character vector of ALT alleles at the site.
#' @return a one-row tibble: `class`, `disc` (designable at all?), `code` (the
#'   IUPAC ALT code for `"degenerate"`, else `NA`), `n_alt`.
#' @examples
#' kasp_site_class("T", c("G", "C")) # degenerate over S = C/G, T left out
#' kasp_site_class("A", "T") # complementary, doomed
#' @family KASP design
#' @export
kasp_site_class <- function(ref, alts) {
  ref <- as.character(ref)[1L]
  alts <- as.character(alts)
  alts <- alts[!is.na(alts) & nzchar(alts)]
  alts <- alts[alts != ref]
  snp <- unique(alts[nchar(alts) == 1L & grepl("^[ACGT]$", alts)])
  out <- function(class, disc, code = NA_character_) {
    tibble::tibble(class = class, disc = disc, code = code, n_alt = length(alts))
  }
  if (!length(snp)) {
    return(out(if (length(alts)) "indel" else "none", FALSE))
  }
  if (length(snp) > 1L) {
    code <- iupac_code(snp)
    return(out("degenerate", !is.na(code), code))
  }
  comp <- paste0(pmin(ref, snp), pmax(ref, snp)) %in% c("AT", "CG")
  out(if (comp) "complementary" else "biallelic", !comp)
}
