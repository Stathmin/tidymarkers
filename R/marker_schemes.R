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
#  in step.  The marker-agnostic primer checks (amplicon length, runs, repeats,
#  dimers, off-targets) belong to every scheme - including a bare `"primer"`
#  scheme that is *not* a KASP at all but the reusable "any primer pair" core.
#  The AS/COM oligo Tm and GC windows belong to the allele-specific class instead:
#  they only make sense where both of those primers exist.
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
# complement is the allele-specific machinery (the AS/COM Tm and GC windows, the
# hierarchical Tm rule, discrimination, the variant footprints, the ALT PCR) - the
# checks that need an allele-anchoring primer.
KASP_CHECK_PRIMER <- c(
  "amplicon length (bp)", "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
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

# The 13 allele-specific assay checks: they read the AS/COM primers (their Tm
# and GC windows, the Tm ordering, the 3' window / footprint) and the
# allele-specific in-silico PCR, so they apply only to a chemistry that actually
# has an allele-anchoring primer.  The three classes partition the 26 registered
# checks (test-criteria guards that).
KASP_CHECK_AS <- c(
  "AS Tm (C)", "COM Tm (C)", "AS GC", "COM GC",
  "COM Tm >= AS Tm",
  "AS 3' window variants (n)", "AS footprint variants (n)", "COM footprint variants (n)",
  "AS 3' window dropout", "AS footprint dropout", "COM footprint dropout",
  "discrimination", "in-silico PCR (REF 1 / ALT 0)"
)

# Family-specific checks.  These are not KASP checks, so they use the MARKER_
# prefix.  The generic primer QC checks above still apply to every scheme.
#
# Membership is held as criterion **ids** (see `criterion_spec()`), never as
# labels: the id is the stable key, the label is presentation.  Each vector is
# named by its label, the way `tibble::tribble()` names its columns - a printf
# aid only, never a source of truth.  `as.vector()` is the id view; the display
# label is asked of the registry through `marker_family_check_labels()`, so a
# criterion can be reworded in one place without orphaning a family list.
MARKER_CHECK_CAPS <- c(
  "restriction site defined" = "family_caps_site_defined",
  "restriction site separates alleles" = "family_caps_site_separates",
  "band gap (bp)" = "family_caps_band_gap",
  "digest sums to amplicon" = "family_caps_digest_spans",
  "restriction site conserved (population)" = "family_caps_site_conserved"
)
MARKER_CHECK_LAMP <- c(
  "FIP 3' ALT" = "family_lamp_fip_alt",
  "inner arm Tm balance (C)" = "family_lamp_inner_tm_balance",
  "inner arms inside amplicon" = "family_lamp_inner_inside",
  "variant at FIP 3' end" = "family_lamp_fip_variant"
)
MARKER_CHECK_RPA <- c(
  "primer length (nt)" = "family_rpa_primer_length",
  "THF off variant" = "family_rpa_thf_off",
  "RPA probe spans variant" = "family_rpa_probe_spans",
  "variant under primer/probe" = "family_rpa_variant_under"
)
MARKER_CHECK_CRISPR <- c(
  "PAM valid" = "family_crispr_pam",
  "variant in seed" = "family_crispr_seed",
  "guide spans variant" = "family_crispr_guide_spans",
  "guide/PAM conserved (population)" = "family_crispr_conserved",
  "off-target protospacers (n)" = "family_crispr_offtarget"
)
MARKER_CHECK_PROBE <- c(
  "probe spans variant" = "family_probe_spans",
  "probe Tm - primer Tm (C)" = "family_probe_tm_delta",
  "ASP 3' discrimination" = "family_probe_asp_discrimination",
  "extension 3' at variant" = "family_probe_extension"
)
MARKER_CHECK_SANGER <- c(
  "variant-read distance (bp)" = "family_sanger_read_distance",
  "read primers inside amplicon" = "family_sanger_primers_inside",
  "read primer pair defined" = "family_sanger_pair_defined"
)
MARKER_CHECK_SIZE <- c(
  "repeat unit (nt)" = "family_size_repeat_unit",
  "repeat copies" = "family_size_repeat_copies",
  "REF/ALT size difference (bp)" = "family_size_delta",
  "fragment size (bp)" = "family_size_fragment"
)

MARKER_FAMILY_CHECKS <- list(
  caps = MARKER_CHECK_CAPS,
  lamp = MARKER_CHECK_LAMP,
  rpa = MARKER_CHECK_RPA,
  crispr = MARKER_CHECK_CRISPR,
  taqman = MARKER_CHECK_PROBE,
  dpcr = MARKER_CHECK_PROBE,
  snapshot = MARKER_CHECK_PROBE,
  castpcr = MARKER_CHECK_PROBE,
  sanger = MARKER_CHECK_SANGER,
  ssr = MARKER_CHECK_SIZE,
  indel_size = MARKER_CHECK_SIZE
)

#' The display labels of a family's checks
#'
#' `MARKER_FAMILY_CHECKS` / `MARKER_CHECK_*` hold criterion ids; text is printed
#' from labels.  This is the one place the two meet, so a caller that compares or
#' prints a family's checks - the plot check column, a test, a report - asks here
#' instead of carrying its own label vector.
#' @param family a scheme family (see `marker_schemes()`); `NULL`/unknown yields
#'   `character(0)`.
#' @return the registry labels of that family's checks, in registry order.
#' @keywords internal
#' @noRd
marker_family_check_labels <- function(family) {
  ids <- MARKER_FAMILY_CHECKS[[family]] %||% character(0)
  if (!length(ids)) {
    return(character(0))
  }
  # the registry's own table, keyed by id, holds the display label: read it there
  # rather than from the constant's `names`, which are a printf aid only
  fam <- criterion_spec()
  label <- fam$label[match(as.character(ids), fam$id)]
  label[!is.na(label)]
}

# family -> the non-primer check classes that chemistry obeys.  Keyed by `family`
# (not by scheme), so a discovered chemistry inherits its family's set with no
# per-scheme edit here.  `site` = the KASP ALT-primer/linkage selector applies
# (the `as`/`pair` stages, plus the discovered `kasp` family such as `lax`);
# `as` = the chemistry grades an AS/COM primer pair (`kasp`/`multiplex`/`asq`).
# Everything else is the marker-agnostic primer subset only.
KASP_FAMILY_CHECKS <- list(
  primer = character(0),
  kasp = c("site", "as"),
  multiplex = c("site", "as"),
  asq = c("site", "as"),
  # castPCR is mutant-only (one ASP, no REF primer / COM allele rule): its own
  # gates enforce mutant specificity, so it takes none of the KASP `as` checks.
  castpcr = character(0),
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

# The specificity policy of a chemistry (roadmap Y3): how strictly a genome-level
# off-target is treated.  The default is `multicopy` - tolerant: a pairing that
# fails the genome-level screen (an off-target exact hit, a second in-silico PCR
# product or an inverted repeat) is warned and the ranked leader is kept (qPCR and
# deliberately multi-copy / CNV loci are legitimate; "don't be heavy-handed").
# `unique` / `homeolog-aware` are strict: the marker is NOT returned, the target is
# doomed with the reason - for chemistries where an off-target is fatal (e.g. a
# homoeolog-specific AS primer in wheat A/B/D).  `off` is the explicit no-op.
KASP_SPECIFICITY_MODES <- c("unique", "homeolog-aware", "multicopy", "off")

# Modes that doom instead of warning.
KASP_SPECIFICITY_STRICT <- c("unique", "homeolog-aware")

#' A scheme's specificity policy (roadmap Y3)
#'
#' The mode [assemble_kasp_assay()] uses to decide whether a pairing that fails
#' the genome-level screen is doomed or merely warned.  Read from the scheme's
#' `scheme_def_<scheme>()$specificity` when set, else the tolerant `"multicopy"`
#' (an off-target is reported, not fatal - "don't be heavy-handed").
#' @param scheme a scheme id (see [marker_schemes()]).
#' @return one of `KASP_SPECIFICITY_MODES`.
#' @keywords internal
#' @noRd
marker_specificity <- function(scheme = "snp") {
  mode <- extra_scheme_defs()[[scheme]]$specificity %||% "multicopy"
  if (!mode %in% KASP_SPECIFICITY_MODES) {
    cli::cli_abort(gettext("Unknown specificity mode {mode} for scheme {scheme}."))
  }
  mode
}

#' The marker chemistries this package can describe
#'
#' One row per chemistry: how it tells the alleles apart, how the signal is read,
#' the zygosity it can call, whether amplicon length is critical, and how many
#' loci share a reaction.  Each scheme is realised in the pipeline as a piece of
#' the decision schema (`marker_decision_schema()`): every scheme carries the
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
#' Membership is family-based (`KASP_FAMILY_CHECKS`): a `KASP_CHECK_PRIMER`
#' criterion belongs to every scheme, a `KASP_CHECK_SITE` criterion to the
#' families whose selector models the ALT-allele set / linkage, and a
#' `KASP_CHECK_AS` criterion to the families with an allele-anchoring primer.  A
#' discovered scheme is tagged from its `family`, never per-scheme, so a new
#' chemistry is one row in `KASP_FAMILY_CHECKS`.
#' @param id the check one row at a time - its `label`, the key the three KASP
#'   classes are declared with.  A family check's label is mapped back to its id
#'   (`check_schemes()` is also called with `MARKER_CHECK_*` ids directly), so a
#'   family vector's own `names` are dropped first.
#' @return a list of scheme-id vectors, one per check.
#' @keywords internal
#' @noRd
check_schemes <- function(id) {
  # every key here is a *value*: a family `MARKER_CHECK_*` vector carries its
  # label as a diagnostic name, never as content
  id <- unname(as.character(id))
  # the probe and size checks serve families that share one vector, so their
  # membership is decided where it already lives, in MARKER_FAMILY_CHECKS, and
  # read back as ids
  member <- function(fam) unname(as.character(MARKER_FAMILY_CHECKS[[fam]] %||% character(0)))
  family_ids <- unique(unlist(lapply(names(MARKER_FAMILY_CHECKS), member)))
  # a label the family table carries -> its id (checked before the three KASP
  # classes, whose labels are never family labels).  Built from MARKER_FAMILY_CHECKS
  # alone, never through `criterion_spec()`, which calls this function.
  pairs <- unlist(unname(lapply(names(MARKER_FAMILY_CHECKS), function(fam) {
    x <- MARKER_FAMILY_CHECKS[[fam]] %||% character(0)
    stats::setNames(unname(as.character(x)), unname(names(x)))
  })))
  label <- names(pairs)[match(id, names(pairs))]
  is_family_label <- !is.na(label)
  id[is_family_label] <- unname(pairs[label[is_family_label]])
  s <- marker_schemes()
  family <- function(fam) KASP_FAMILY_CHECKS[[fam]] %||% character(0)
  in_class <- function(klass) {
    s$scheme[vapply(s$family, function(fam) klass %in% family(fam), logical(1))]
  }
  lapply(id, function(l) {
    if (l %in% KASP_CHECK_PRIMER) {
      s$scheme
    } else if (l %in% KASP_CHECK_SITE) {
      in_class("site")
    } else if (l %in% KASP_CHECK_AS) {
      in_class("as")
    } else if (l %in% family_ids) {
      s$scheme[vapply(s$family, function(fam) l %in% member(fam), logical(1))]
    } else {
      character(0)
    }
  })
}
