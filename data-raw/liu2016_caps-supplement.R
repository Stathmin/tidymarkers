# =============================================================================
#  data-raw/liu2016_caps-supplement.R  --  watermelon CAPS conversion funnel
# =============================================================================
#  Source (identifier resolution is PROVENANCE, not evidence; every number
#  below is re-read from the fetched full text):
#
#    Liu C, Gao L, Zhu Q, et al.  Development of genome-wide cleaved amplified
#    polymorphic sequence (CAPS) markers in watermelon (Citrullus lanatus) and
#    their use in linkage map construction.  Breeding Science 2016; 66(2):
#    244-259.  DOI 10.1270/jsbbs.66.244  (PMC4785002)
#
#  Fetched: .tmp/lit/liu2016_caps/fulltext.xml (Europe PMC fullTextXML).  The
#  full text contains NO <supplementary-material> node and no "Table S" / "suppl"
#  mention at all, so the paper's tables ARE the record - there is no supplement
#  we are failing to parse.
#
#  WHY this paper: it is the only large CAPS-conversion study in the registry
#  that reports the full amplification-and-digestion funnel for a
#  genome-wide design round (532 primer pairs), which is what the package needs
#  to show that "CAPS conversion" is a chain, not a single percentage.
#
#  WHAT THIS SCRIPT DOES NOT MODEL  (read before quoting a rate):
#
#   * The 271 is "distinct bands of the expected length AND polymorphic between
#     P1/P2/F1", scored on gels.  It is NOT "assays that worked": 516/532
#     amplified (16 failed to amplify), and of the 532 only 296 gave the
#     expected amplicon size.  The paper then reports 271 as the polymorphic
#     count with a "polymorphic rate of 50.94%".  NOTE THE ARITHMETIC
#     DISCREPANCY, which we record rather than resolve: 271/532 = 50.94%,
#     but the text also says "approximately 42.64% of the primer pairs
#     amplified additional fragments or failed in enzyme digestion reaction",
#     which is not 100% - 50.94% = 49.06%, since it is presumably taken over
#     the 296 expected-size amplicons (271/296 = 91.6% would be different, and
#     the paper's own 1 - 0.4264 = 57.36% only matches 271/532 if rounding is
#     generous). The paper does not state the denominator for the 42.64%
#     figure, so we extract the counts and flag that the two percentages
#     cannot both be right against the same denominator of 532.
#   * There is NO per-primer table.  The paper's Tables 1-7 are SNP density,
#     SNP types, linkage map, phenotype means, QTLs, candidate-gene annotation
#     and trait correlations - none lists the 532 primer pairs or their
#     outcome.  So the per-marker CSV cannot be written from the paper; only
#     the stage totals exist.  The 301 mapped markers and 1836.51 cM are
#     linkage-map facts, not conversion outcomes, and stay out.
#   * No rate is converted into a probability for internal use.
#
#  Run:  Rscript data-raw/liu2016_caps-supplement.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

xml_path <- fs::path(root, ".tmp", "lit", "liu2016_caps", "fulltext.xml")
if (!file.exists(xml_path)) stop("Missing ", xml_path)
doc <- xml2::read_xml(xml_path)
out_dir <- fs::path(root, "inst", "extdata")
DOI <- "10.1270/jsbbs.66.244"

## --- Table 2: the SNP-type split, used as an arithmetic cross-check ---------
# Re-read from the table, not typed.  Rows are read by CONTENT: blank spacer
# rows are dropped, a family row carries its own total in the last cell, and
# continuation rows belong to the family above them.
read_table <- function(doc, k) {
  trs <- xml2::xml_find_all(
    xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")[[k]],
    ".//*[local-name()='tr']"
  )
  lapply(trs, function(r) {
    v <- xml2::xml_text(xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']"))
    # KEEP empty cells as position placeholders (continuation rows leave the
    # family cell blank); drop only rows that are empty throughout.
    if (!any(nzchar(trimws(v)))) character(0) else trimws(v)
  })
}
num <- function(x) suppressWarnings(as.integer(gsub("[^0-9]", "", x)))

tbl2 <- read_table(doc, 2L)
# Rows are (family?, change, amount, ratio, total?) with EMPTY STRINGS kept as
# placeholders, so we name the cells by position instead of by length.
# Drop the header row (its second cell is the literal word "Type").
tbl2 <- Filter(function(v) length(v) == 5L && !identical(v[2L], "Type"), tbl2)

fam_of <- c("Transversions", "Transitions", "Indels")
fam_total <- stats::setNames(integer(0), character(0))
cur <- NA_character_
recs <- list()
for (v in tbl2) {
  fam <- v[1L]
  if (nzchar(fam)) {
    if (!fam %in% fam_of) stop("Unexpected Table 2 family label: ", fam)
    cur <- fam
    fam_total[[cur]] <- num(v[5L])
  }
  if (!nzchar(v[2L])) next
  recs[[length(recs) + 1L]] <- tibble::tibble(
    family = cur, change = v[2L], n = num(v[3L]), pct = v[4L]
  )
}
tbl2_long <- dplyr::bind_rows(recs)
fam_total <- as.integer(fam_total[fam_of])

## --- the conversion funnel, as the Section "CAPS polymorphism" states it ----
funnel <- tibble::tibble(
  stage = c("CAPS loci excavated by SNP2CAPS with 16 enzymes",
            "candidate sequence pairs suitable for primer design",
            "primer pairs chosen for design and polymorphism analysis",
            "primer pairs that amplified a clear PCR product",
            "primer pairs with the expected amplicon size",
            "primer pairs giving expected-length AND polymorphic bands"),
  n = c(14246L, 3192L, 532L, 516L, 296L, 271L)
)

## --- the assertion that makes this a *verified* extract --------------------
sum_fam <- vapply(fam_of, function(f) sum(tbl2_long$n[tbl2_long$family == f]), integer(1))
# Table 1's Total row prints the genome-wide SNP and indel counts; Table 2's
# three family totals must add up to exactly those, which ties the two tables.
tbl1 <- read_table(doc, 1L)
tot <- Filter(function(v) length(v) == 5L && v[1L] == "Total", tbl1)
if (length(tot) != 1L) stop("Table 1 has no single Total row")
total_nums <- num(tot[[1L]])
grand <- sum(total_nums[!is.na(total_nums)])

## --- the assertion that makes this a *verified* extract --------------------
sum_fam <- vapply(fam_of, function(f) sum(tbl2_long$n[tbl2_long$family == f]), integer(1))
funnel_ok <- all(diff(funnel$n) <= 0) && all(funnel$n > 0)
pct_ok <- abs(round(100 * 271 / 532, 2) - 50.94) < 0.005

if (!identical(unname(sum_fam), unname(fam_total)) ||
    sum(fam_total) != grand ||
    !funnel_ok || !identical(as.integer(table(tbl2_long$family)), c(4L, 2L, 4L)) ||
    !pct_ok) {
  stop(
    "Table 2 / funnel transcription failed: family sums ",
    paste(sum_fam, collapse = "/"), " vs printed ",
    paste(fam_total, collapse = "/"), "; grand total ", sum(fam_total),
    " vs Table 1 Total = ", grand, "; funnel monotone: ", funnel_ok,
    "; 271/532 = ", round(100 * 271 / 532, 2), "% vs stated 50.94%"
  )
}

OUTCOME_DEF <- paste(
  "CAPS 'success' in this paper = the primer pair gave a distinct band of the",
  "EXPECTED length that was POLYMORPHIC across the three genotypes tested",
  "(P1 LSW-177, P2 COS, F1), scored after PCR and restriction digestion.",
  "Amplification success (516/532) and expected amplicon size (296) are",
  "reported as separate stages and are not the same outcome."
)

readr::write_csv(
  dplyr::mutate(funnel, outcome_definition = OUTCOME_DEF, doi = DOI),
  fs::path(out_dir, "liu2016_caps_funnel.csv")
)
readr::write_csv(
  dplyr::mutate(tbl2_long, outcome_definition = NA_character_, doi = DOI),
  fs::path(out_dir, "liu2016_caps_snp_types.csv")
)

message("liu2016_caps: Table 2 families ",
        paste(sum_fam, collapse = "/"), " sum to ", sum(fam_total),
        "; funnel ", paste(funnel$n, collapse = " -> "),
        "; 271/532 = ", round(100 * 271 / 532, 2), "% (assertion passed)")
