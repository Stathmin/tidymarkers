# =============================================================================
#  kasp_schemes.R  --  the marker chemistries tidykasp can describe
# =============================================================================
#  KASP is a family, not one assay.  Every member shares the same skeleton - an
#  allele-specific primer that anchors the variant and a common primer - but they
#  differ in *what* separates the alleles (a 3'-terminal base, an amplicon
#  length, an IUPAC code over a degenerate ALT set) and in *how many loci* share
#  one tube.  A duplex chemistry (dKASP, mxKASP, Multi-KASP) puts two loci in one
#  reaction, so it adds a compatibility stage on top of the per-locus design;
#  ASQ / PACE is a different cassette readout but the same allele-specific core.
#
#  This registry names those chemistries.  `kasp_criteria()` tags each criterion
#  with the schemes it applies to, so the decision schema, the text export and
#  the picture can be filtered per chemistry without a second vocabulary to keep
#  in step.  The marker-agnostic primer checks (Tm, GC, runs, repeats, dimers,
#  off-targets) belong to every scheme - including a bare `"primer"` scheme that
#  is *not* a KASP at all but the reusable "any primer pair" core.
# =============================================================================

# Which chemistries a criterion serves.  `primer` is the marker-agnostic core
# (any PCR primer pair); the rest are allele-specific assays that inherit it.
KASP_SCHEME_PRIMER <- "primer"
KASP_SCHEME_KASP <- c("snp", "indel", "multi")
KASP_SCHEME_CAPS <- "caps"
KASP_SCHEME_ARMS <- "arms"
KASP_SCHEME_CORE <- c(KASP_SCHEME_KASP, "duplex", "asq")
KASP_SCHEME_ALL <- c(KASP_SCHEME_PRIMER, KASP_SCHEME_CORE, KASP_SCHEME_CAPS, KASP_SCHEME_ARMS)

# The acceptance checks that are true of any primer pair, KASP or not.  Their
# complement is the allele-specific machinery (discrimination, the hierarchical
# AS/COM Tm rule, the variant footprints, the ALT PCR) - the part that *is* KASP.
KASP_CHECK_PRIMER <- c(
  "amplicon length (bp)", "AS Tm (C)", "COM Tm (C)", "AS GC", "COM GC",
  "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
  "inverted repeat (arms>=4, loop<=6)", "worst hairpin dG37",
  "worst cross-dimer dG37", "worst homodimer dG37",
  "off-target exact primer hits (n)"
)

#' The marker chemistries this package can describe
#'
#' One row per chemistry: how it tells the alleles apart, how the signal is read,
#' the zygosity it can call, whether amplicon length is critical, and how many
#' loci share a reaction.  Each scheme is realised in the pipeline as a piece of
#' the decision schema ([kasp_decision_schema()]): the marker-agnostic primer
#' checks are shared by all of them, the allele-specific stages by
#' `family != "primer"`, and a duplex chemistry adds its compatibility stage.
#'
#' The `"primer"` scheme is deliberately **not** a KASP: it is the reusable
#' "any primer pair" subset (composition, structure, repeats, off-targets),
#' useful for plain PCR / allele-specific PCR primer QC where no allele
#' discrimination is claimed.
#' @return tibble with one row per scheme: `scheme`, `family`, `label`,
#'   `discrimination`, `readout`, `zygosity`, `amplicon`, `loci`, `verdict`.
#' @examples
#' kasp_schemes()
#' @family KASP checks
#' @export
kasp_schemes <- function() {
  tibble::tribble(
    ~scheme, ~family, ~label, ~discrimination, ~readout, ~zygosity, ~amplicon, ~loci, ~verdict,
    "primer", "primer", "any primer pair (marker-agnostic)",
    "none - not an allele assay", "gel / melt / any", "n/a", "bounded", 1L,
    "usable primer pair",
    "snp", "kasp", "SNP KASP (biallelic)",
    "3'-terminal base + FRET ratio", "co-dominant two-colour", "co-dominant", "critical", 1L,
    "usable KASP assay",
    "indel", "kasp", "indel KASP",
    "allele-specific amplicon length", "co-dominant two-colour", "co-dominant", "critical", 1L,
    "usable KASP assay",
    "multi", "kasp", "multi-allelic / degenerate KASP",
    "3'-terminal IUPAC code over the ALT set", "co-dominant two-colour", "co-dominant", "critical", 1L,
    "usable KASP assay",
    "duplex", "multiplex", "duplex KASP (dKASP / mxKASP / Multi-KASP)",
    "two independent loci in one tube", "two cassettes (or presence/absence)",
    "co-dominant or dominant", "<= 150 bp", 2L,
    "usable duplex assay",
    "asq", "asq", "ASQ / PACE",
    "allele-specific real-time extension", "real-time, up to 4-plex", "co-dominant", "bounded", 1L,
    "usable ASQ assay",
    "caps", "caps", "CAPS / dCAPS (restriction marker)",
    "the allele creates or destroys a restriction site", "gel fragment lengths", "co-dominant",
    "bounded", 1L, "usable CAPS marker",
    "arms", "arms", "tetra-primer ARMS-PCR",
    "allele-specific inner primers, sized gel bands", "agarose band sizes (AA/AB/BB)",
    "co-dominant", "bounded", 1L, "usable ARMS marker"
  )
}

#' One scheme's row, validated
#' @param scheme a scheme id (see [kasp_schemes()]).
#' @return a one-row tibble.
#' @keywords internal
#' @noRd
kasp_scheme <- function(scheme = "snp") {
  s <- kasp_schemes()
  i <- match(scheme, s$scheme)
  if (is.na(i)) {
    rlang::abort(gettextf(
      "Unknown marker scheme %s; see kasp_schemes().",
      encodeString(as.character(scheme), quote = '"')
    ))
  }
  s[i, , drop = FALSE]
}

#' The schemes each check applies to (a list, one element per label)
#' @keywords internal
#' @noRd
check_schemes <- function(label) {
  out <- rep(list(KASP_SCHEME_CORE), length(label))
  out[label %in% KASP_CHECK_PRIMER] <- list(KASP_SCHEME_ALL)
  out
}
