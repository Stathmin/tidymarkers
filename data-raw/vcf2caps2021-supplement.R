# =============================================================================
#  data-raw/vcf2caps2021-supplement.R  --  VCF2CAPS: what the paper actually
#                                           reported, and what it did not
# =============================================================================
#  Source (identifier resolution is PROVENANCE, not evidence):
#
#    Wesołowski M, Domnicz A, Augustynowicz J, et al.  VCF2CAPS - A
#    high-throughput CAPS marker design from VCF files and its test-use on a
#    genotyping-by-sequencing (GBS) dataset.  PLoS Computational Biology 2021;
#    17(5):e1008980.  DOI 10.1371/journal.pcbi.1008980
#
#  Fetched: .tmp/lit/vcf2caps2021/fulltext.xml (Europe PMC fullTextXML).
#
#  ===========================================================================
#  READ THIS BEFORE ADDING A NUMBER FROM THIS PAPER TO ANYTHING
#  ===========================================================================
#  The registry row says this paper "reports NO wet-lab outcomes".  THAT IS
#  WRONG and this file records the correction.  The paper has a Results
#  subsection, "Experimental verification of the identified CAPS markers",
#  reporting 10 single-cut markers assayed on 30 plants each (Fig 8), and it
#  reports the outcome for each marker IN PROSE.
#
#  What is still true, and is the reason NO numeric outcome rate is written to
#  inst/extdata/ for this paper:
#
#   * The per-marker result is a qualitative prose verdict, not a count.  The
#     paper says things like "distribution ... exactly matched genotypic
#     predictions" (1_028, 3_579), "only the discrepancy was ... a single
#     supposed heterozygote" (3_714, 3_300, 6_022), "five supposed
#     non-cleavable homozygotes were in fact heterozygotes" (3_363).  It never
#     gives a table of per-plant calls, and Fig 8 is a gel image.  Any
#     "x/30 correct" we computed would be INVENTED, so we do not.
#   * The paper's own headline numbers for this section are the DESIGN counts
#     (16,093 CAPS markers; 3,331 single-cut; Table 4), which are site-
#     detection outputs, not wet-lab conversions.  Those are the only numbers
#     VCF2CAPS is entitled to anchor, and they are already the kind of thing
#     the package computes rather than reads.
#
#  So this file writes ONE artifact: a verbatim, per-marker outcome-definition
#  ledger with NO counts, plus an explicit `extractable = FALSE` flag with the
#  reason.  It exists so that (a) nobody re-derives a fake concordance rate,
#  and (b) the registry's claim is corrected in the repo rather than in a
#  commit message.
#
#  WHAT IS ASSERTED: the ledger contains exactly the markers that the Results
#  text names, 10 in total, and they are exactly the 10 markers in Table 1
#  (which is the table of "CAPS markers selected for experimental verification").
#  If the paper is revised and the two lists diverge, this script stops.
#
#  Run:  Rscript data-raw/vcf2caps2021-supplement.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

xml_path <- fs::path(root, ".tmp", "lit", "vcf2caps2021", "fulltext.xml")
if (!file.exists(xml_path)) stop("Missing ", xml_path)
doc <- xml2::read_xml(xml_path)
out_dir <- fs::path(root, "inst", "extdata")
DOI <- "10.1371/journal.pcbi.1008980"

## --- Table 1: the 10 markers taken forward to the bench ----------------------
tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
cap <- vapply(tw, function(t) xml2::xml_text(xml2::xml_find_first(
  t, ".//*[local-name()='caption']")), character(1))
tbl1 <- tw[grepl("selected for experimental verification", cap)]
if (length(tbl1) != 1L) stop("Could not locate the verification Table 1")
trs <- xml2::xml_find_all(tbl1[[1L]], ".//*[local-name()='tr']")
cells <- lapply(trs, function(r)
  xml2::xml_text(xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']")))
data_rows <- Filter(function(v) length(v) == 6L && grepl("^[0-9]", v[1L]), cells)
if (length(data_rows) != 10L) stop("Expected 10 marker rows in verification Table 1, got ", length(data_rows))
table1 <- tibble::tibble(
  marker_original = vapply(data_rows, `[[`, character(1), 1L),
  marker_abbrev   = vapply(data_rows, `[[`, character(1), 2L),
  chromosome      = vapply(data_rows, `[[`, character(1), 3L),
  forward_primer  = vapply(data_rows, `[[`, character(1), 4L),
  reverse_primer  = vapply(data_rows, `[[`, character(1), 5L),
  enzyme          = vapply(data_rows, `[[`, character(1), 6L)
)

## --- the per-marker verdicts, quoted verbatim from the Results --------------
# Each string is the paper's own wording, shortened only by "..." where a
# sentence is elided; no number was added.
verdicts <- c(
  `1_028` = "distribution of their restriction profiles exactly matched genotypic predictions from the VCF2CAPS output file",
  `3_579` = "distribution of their restriction profiles exactly matched genotypic predictions from the VCF2CAPS output file",
  `3_351` = "similar situation ... although the band of the cleavable allele was weaker than that of the non-cleavable allele - especially in heterozygotes",
  `3_714` = "distribution very close to predictions, the only discrepancy was shown for single supposed heterozygotes which marker-wise turned out to be homozygotes - cleavable",
  `3_300` = "distribution very close to predictions, the only discrepancy was shown for single supposed heterozygotes which marker-wise turned out to be homozygotes - non-cleavable",
  `6_022` = "distribution very close to predictions, the only discrepancy was shown for single supposed heterozygotes which marker-wise turned out to be homozygotes - non-cleavable",
  `3_363` = "predicted profiles were observed for five supposed non-cleavable homozygotes as well as for all supposed heterozygotes and cleavable homozygotes; the remaining five supposed non-cleavable homozygotes were in fact heterozygotes",
  `1_017` = "confirmed the character of the supposed non-cleavable homozygotes and heterozygotes; profiles for all supposed cleavable homozygotes were similar to those of heterozygotes, although with the cleavage products stronger and the non-cleaved DNA fragment weaker",
  `2_525` = "predictions confirmed in all supposed homozygotes (both non-cleavable and cleavable), but all supposed heterozygotes looked like the non-cleavable homozygotes",
  `s_256` = "predictions confirmed in all supposed homozygotes (both non-cleavable and cleavable), but all supposed heterozygotes looked like the non-cleavable homozygotes"
)

## --- the assertion ----------------------------------------------------------
# The verdict keys must be exactly the abbreviated marker ids of Table 1
# (en dashes normalised, because the table sets them as "1–028").
abbrev <- gsub("\u2013", "_", table1$marker_abbrev)
abbrev <- gsub("-", "_", abbrev)
if (nrow(table1) != 10L || !setequal(names(verdicts), abbrev)) {
  stop(
    "Table 1 and the Results prose disagree: ", nrow(table1),
    " markers in Table 1 (", paste(abbrev, collapse = ", "),
    ") vs ", length(verdicts), " in the verdict ledger (",
    paste(names(verdicts), collapse = ", "), ")"
  )
}

OUTCOME_DEF <- paste(
  "VCF2CAPS 'success' = for a marker, the observed restriction profiles of 30",
  "plants (10 predicted non-cleavable homozygotes, 10 heterozygotes, 10",
  "cleavable homozygotes) matched the genotype classes VCF2CAPS predicted. It",
  "is a CLASS-ASSIGNMENT agreement judged on a gel, reported in prose per",
  "marker; the paper publishes no per-plant call table and no count of",
  "'correct' plants, so no numeric rate is recorded here."
)

ledger <- table1 |>
  dplyr::mutate(
    outcome_definition = OUTCOME_DEF,
    reported_verdict_verbatim = unname(verdicts[abbrev]),
    n_plants_per_marker = 30L,
    n_plants_class_breakdown = "10 non-cleavable hom / 10 het / 10 cleavable hom (design, not result)",
    extractable = FALSE,
    not_extractable_reason = paste(
      "outcome is a qualitative prose verdict; no per-plant or per-class count",
      "and no supplementary data table exist, so any success rate would be",
      "invented"
    ),
    doi = DOI
  )

readr::write_csv(ledger, fs::path(out_dir, "vcf2caps2021_verification.csv"))

message("vcf2caps2021: ", nrow(ledger),
        " markers carried a prose verdict, matching Table 1 exactly; ",
        "NO counts written (extractable = FALSE) (assertion passed)")
