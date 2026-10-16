# =============================================================================
#  marker_schemes.R  --  the marker chemistries tidymarkers can describe
# =============================================================================
#  KASP is a family, not one assay.  Every member shares the same skeleton - an
#  allele-specific primer that anchors the variant and a common primer - but they
#  differ in *what* separates the alleles (a 3'-terminal base, an amplicon
#  length, an IUPAC code over a degenerate ALT set) and in *how many loci* share
#  one tube.  A duplex chemistry (dKASP, mxKASP, Multi-KASP) puts two loci in one
#  reaction, so it adds a compatibility stage on top of the per-locus design;
#  ASQ / PACE is a different cassette readout but the same allele-specific core.
#
#  This registry names those chemistries.  `marker_criteria()` tags each criterion
#  with the schemes it applies to, so the decision schema, the text export and
#  the picture can be filtered per chemistry without a second vocabulary to keep
#  in step.  The marker-agnostic primer checks (Tm, GC, runs, repeats, dimers,
#  off-targets) belong to every scheme - including a bare `"primer"` scheme that
#  is *not* a KASP at all but the reusable "any primer pair" core.
#
#  Membership is decided by the scheme's `family`, not by hand: the acceptance
#  checks fall into three classes (marker-agnostic, KASP-selector *site*, and
#  allele-specific) and KASP_FAMILY_CHECKS says which classes a chemistry obeys.
#  So a non-KASP chemistry (PROBE/LAMP/CRISPR/Sanger/...) never inherits the
#  AS/COM/3'-window/discrimination checks that only make sense with an
#  allele-anchoring primer - the bug C2 fixes.
# =============================================================================

# Scheme groups.  `primer` is the marker-agnostic core (any PCR primer pair);
# `KASP_SCHEME_CORE` is the built-in set whose selector runs the `as` / `pair`
# stages.  Check *membership* is family-based (see KASP_FAMILY_CHECKS below), so
# these vectors only group schemes that share a stage.
KASP_SCHEME_PRIMER <- "primer"
KASP_SCHEME_KASP <- c("snp", "indel", "multi")
KASP_SCHEME_CAPS <- "caps"
KASP_SCHEME_ARMS <- "arms"
KASP_SCHEME_HRM <- "hrm"
KASP_SCHEME_CORE <- c(KASP_SCHEME_KASP, "duplex", "asq")
KASP_SCHEME_ALL <- c(
  KASP_SCHEME_PRIMER, KASP_SCHEME_CORE,
  KASP_SCHEME_CAPS, KASP_SCHEME_ARMS, KASP_SCHEME_HRM
)

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

# The 5 target-*site* checks: the KASP *selector* notions - the ALT-allele set,
# how many of them the ALT primer covers, and whether the design leans on linkage
# (multi / lax).  They apply only to a chemistry whose selector models an ALT
# primer and linkage, never to a probe/LAMP/CRISPR/flanking-pair assay.
KASP_CHECK_SITE <- c(
  "target site class", "ALT alleles (n)", "linked variants skipped (n)",
  "ALT primer covers all ALTs", "design relies on linkage"
)

# The 9 allele-specific assay checks: they read the AS primer (its Tm ordering,
# its 3' window / footprint) and the allele-specific in-silico PCR, so they apply
# only to a chemistry that actually has an allele-anchoring primer.  The three
# classes partition the 26 registered checks (test-criteria guards that).
KASP_CHECK_AS <- c(
  "COM Tm >= AS Tm",
  "AS 3' window variants (n)", "AS footprint variants (n)", "COM footprint variants (n)",
  "AS 3' window dropout", "AS footprint dropout", "COM footprint dropout",
  "discrimination", "in-silico PCR (REF 1 / ALT 0)"
)

# family -> the non-primer check classes that chemistry obeys.  Keyed by `family`
# (not by scheme), so a discovered chemistry inherits its family's set with no
# per-scheme edit here.  `site` = the KASP ALT-primer/linkage selector applies
# (the `as`/`pair` stages, plus the discovered `kasp` family such as `lax`);
# `as` = the chemistry has an allele-anchoring primer (`kasp`/`multiplex`/`asq`
# and castPCR's 3'-terminal AS primer).  Everything else is the marker-agnostic
# primer subset only.
KASP_FAMILY_CHECKS <- list(
  primer = character(0),
  kasp = c("site", "as"),
  multiplex = c("site", "as"),
  asq = c("site", "as"),
  castpcr = c("site", "as"),
  caps = character(0),
  arms = character(0),
  hrm = character(0),
  crispr = character(0),
  dpcr = character(0),
  indel_size = character(0),
  lamp = character(0),
  rpa = character(0),
  sanger = character(0),
  snapshot = character(0),
  ssr = character(0),
  taqman = character(0)
)

#' The marker chemistries this package can describe
#'
#' One row per chemistry: how it tells the alleles apart, how the signal is read,
#' the zygosity it can call, whether amplicon length is critical, and how many
#' loci share a reaction.  Each scheme is realised in the pipeline as a piece of
#' the decision schema ([marker_decision_schema()]): every scheme carries the
#' marker-agnostic primer checks, and the allele-specific / target-site checks
#' follow the scheme's `family` (see `KASP_FAMILY_CHECKS`), so a non-KASP
#' chemistry is never graded on an AS primer it does not have.  A duplex
#' chemistry adds its compatibility stage.
#'
#' The `"primer"` scheme is deliberately **not** a KASP: it is the reusable
#' "any primer pair" subset (composition, structure, repeats, off-targets),
#' useful for plain PCR / allele-specific PCR primer QC where no allele
#' discrimination is claimed.
#' @return tibble with one row per scheme: `scheme`, `family`, `label`,
#'   `discrimination`, `readout`, `zygosity`, `amplicon`, `loci`, `verdict`.
#' @examples
#' marker_schemes()
#' @family KASP checks
#' @export
marker_schemes <- function() {
  out <- tibble::tribble(
    ~scheme, ~family, ~label, ~discrimination, ~readout, ~zygosity, ~amplicon, ~loci, ~verdict,
    "primer", "primer", "any primer pair (marker-agnostic)",
    "none - not an allele assay", "gel / melt / any", "n/a", "bounded", 1L,
    "usable primer pair",
    "snp", "kasp", "SNP KASP (biallelic)",
    "3'-terminal base + FRET ratio", "co-dominant two-colour", "co-dominant", "critical", 1L,
    "usable KASP assay",
    "indel", "kasp", "indel KASP",
    "3'-terminal junction base + FRET ratio", "co-dominant two-colour", "co-dominant", "critical", 1L,
    "usable KASP assay",
    "multi", "kasp", "multi-allelic / degenerate KASP",
    "3'-terminal IUPAC code over the ALT set", "co-dominant two-colour", "co-dominant (REF vs non-REF)", "critical", 1L,
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
    "co-dominant", "bounded", 1L, "usable ARMS marker",
    "hrm", "hrm", "HRM (high-resolution melting)",
    "amplicon melting curve (neighbour-tolerant)", "melt curve", "co-dominant",
    "short (50-100 bp)", 1L, "usable HRM amplicon"
  )
  # discovered schemes (R/marker_<scheme>.R via scheme_def_*()) are appended here,
  # so a new chemistry is one file, not an edit to this registry
  dplyr::bind_rows(out, extra_scheme_rows())
}

#' One scheme's row, validated
#' @param scheme a scheme id (see [marker_schemes()]).
#' @return a one-row tibble.
#' @keywords internal
#' @noRd
marker_scheme <- function(scheme = "snp") {
  s <- marker_schemes()
  i <- match(scheme, s$scheme)
  if (is.na(i)) {
    cli::cli_abort(gettext("Unknown marker scheme {scheme}; see marker_schemes()."))
  }
  s[i, , drop = FALSE]
}

#' The schemes each check applies to (a list, one element per label)
#'
#' Membership is family-based ([KASP_FAMILY_CHECKS]): a `KASP_CHECK_PRIMER`
#' criterion belongs to every scheme, a `KASP_CHECK_SITE` criterion to the
#' families whose selector models the ALT-allele set / linkage, and a
#' `KASP_CHECK_AS` criterion to the families with an allele-anchoring primer.  A
#' discovered scheme is tagged from its `family`, never per-scheme, so a new
#' chemistry is one row in `KASP_FAMILY_CHECKS`.
#' @param label the check labels (see [marker_criteria()]).
#' @return a list of scheme-id vectors, one per label.
#' @keywords internal
#' @noRd
check_schemes <- function(label) {
  s <- marker_schemes()
  family <- function(f) KASP_FAMILY_CHECKS[[f]] %||% character(0)
  in_class <- function(cls) {
    s$scheme[vapply(s$family, function(f) cls %in% family(f), logical(1))]
  }
  lapply(label, function(l) {
    if (l %in% KASP_CHECK_PRIMER) {
      s$scheme
    } else if (l %in% KASP_CHECK_SITE) {
      in_class("site")
    } else if (l %in% KASP_CHECK_AS) {
      in_class("as")
    } else {
      character(0)
    }
  })
}
