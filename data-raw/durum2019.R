# =============================================================================
#  data-raw/durum2019.R  --  tetraploid vs hexaploid wheat, the same array panel
# =============================================================================
#  Vendramin V, Ormanbekova D, Scalabrin S, et al.  "Genomic tools for durum
#  wheat breeding: de novo assembly of Svevo transcriptome and SNP discovery in
#  elite germplasm", BMC Genomics 20:268 (2019), 10.1186/s12864-019-5645-x,
#  PMC6456968.
#
#  Why this source is the prize of the family: it puts T. aestivum (hexaploid,
#  AABBDD) SNP assays and T. durum (tetraploid, AABB) SNP assays on the SAME
#  90K iSelect array, scored on the SAME panel of 288 durum accessions.  So the
#  D genome's presence in the assay design, and the A/B homoeologs' abundance in
#  the design, are held apart from every other difference between studies.
#
#  Outcome here is NOT "KASP conversion".  It is the classification the array
#  software returns: failed (no call) vs polymorphic (a usable call that varies)
#  vs monomorphic/informative.  The paper's own footnote fixes the denominator
#  of the polymorphic percentage to the NON-FAILED assays - i.e. it is a
#  yield-of-informative-rate among assays that returned a call, which is not the
#  array conversion rate.  Both are written, labelled.
#
#  Table 2 of the article (not a supplement) is the source.  It is also re-read
#  from the XML: the numbers below are parsed out of the row cells, and the run
#  stops if the five printed percentages stop reproducing from the counts.
#
#  Fetch:
#    curl -sS -L -m 90 -o .tmp/lit/pmc6456968/fulltext.xml \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6456968/fullTextXML'
#  (the supplementary XLSX files carry S1-S12, no Table 2)
#
#  Run:  Rscript data-raw/durum2019.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- Sys.getenv("DURUM2019_DIR", file.path(root, ".tmp", "lit", "pmc6456968"))
xml_path <- Sys.getenv("DURUM2019_XML", file.path(lit, "fulltext.xml"))
if (!file.exists(xml_path)) {
  stop("Fetch fulltext.xml into ", lit, " (see the header), or set DURUM2019_XML.")
}

DOI <- "10.1186/s12864-019-5645-x"
out_dir <- fs::path(root, "inst", "extdata")

## --- parse Table 2 out of the JATS XML --------------------------------------
# The table markup is namespaced; stripping it makes .//table Row-class selectors
# work without a d1: prefix, which is what data-raw/lib does elsewhere.
doc <- xml2::read_xml(xml_path)
doc <- xml2::xml_ns_strip(doc)
tabs <- xml2::xml_find_all(doc, ".//table-wrap")
cap <- vapply(tabs, function(t) paste(xml2::xml_text(
  xml2::xml_find_all(t, ".//caption")), collapse = " "), character(1))
hit <- grep("Number of failed, polymorphic and total SNPs", cap, fixed = TRUE)
if (length(hit) != 1L) stop("Could not find exactly one Table 2 in the XML.")
rows <- xml2::xml_find_all(tabs[[hit]], ".//tr")
cells <- lapply(rows, function(r) {
  trimws(gsub("[[:space:]]+", " ",
              xml2::xml_text(xml2::xml_find_all(r, ".//td|.//th"))))
})

## --- assertions: values asserted against the paper's own totals -------------
# counts parsed from the cells, integers only; no count is typed by hand
int_of <- function(x) {
  # "4,782 (6.5%)" -> 4782: the count is everything before the parenthesis,
  # minus the thousands separator.  Stripping every non-digit would give
  # 478,265 by swallowing the percentage.
  as.integer(gsub("[^0-9]", "", sub("[[:space:]]*\\(.*$", "", x)))
}
pct_of <- function(x) {
  as.numeric(sub("^.*\\(([0-9.]+)%\\).*$", "\\1", x))
}

body <- Filter(function(r) length(r) >= 4L && grepl("SNP|Total", r[[1L]]), cells)
if (length(body) != 3L) {
  stop("Table 2 should have three data rows (T. aestivum, T. durum, Total); got ",
       length(body))
}
# NB: `body` is a list of 4-CHARACTER vectors (one per cell), so a vapply()
# over cells is type-stable but a do.call(rbind, ...) on an unnamed list
# silently column-binds.  Rows are built by NAME, one assertion per cell.
row_of <- function(r) {
  list(group = r[[1L]], total = int_of(r[[2L]]), failed = int_of(r[[3L]]),
       failed_pct = pct_of(r[[3L]]), polymorphic = int_of(r[[4L]]),
       polymorphic_pct = pct_of(r[[4L]]))
}
rows <- lapply(body, row_of)
num <- function(nm) vapply(rows, function(x) x[[nm]], numeric(1))
tab <- data.frame(
  group = vapply(rows, function(x) x$group, character(1)),
  total = num("total"), failed = num("failed"),
  failed_pct = num("failed_pct"), polymorphic = num("polymorphic"),
  polymorphic_pct = num("polymorphic_pct"),
  stringsAsFactors = FALSE
)
# the header row must name the columns exactly as read
hdr <- paste(Filter(length, cells)[[1L]], collapse = " ")
if (!grepl("Failed SNPs", hdr) || !grepl("Polymorphic SNPs", hdr)) {
  stop("Table 2 header does not match the expected column names: ", hdr)
}

n_aes <- tab[grepl("aestivum", tab$group), ]
n_dur <- tab[grepl("durum", tab$group), ]
n_tot <- tab[grepl("^Total", tab$group), ]
if (nrow(n_aes) != 1L || nrow(n_dur) != 1L || nrow(n_tot) != 1L) {
  stop("Table 2 row labels are not one each of aestivum / durum / Total.")
}

## --- the paper's own totals ------------------------------------------------
if (n_tot$total != 81587) stop("Table 2 total is not 81,587: ", n_tot$total)
if (n_aes$total + n_dur$total != n_tot$total) {
  stop("73,647 + 7,940 does not equal the table's 81,587.")
}
if (n_aes$failed + n_dur$failed != n_tot$failed) {
  stop("failed SNPs do not add up to the table's 5,485: ",
       n_aes$failed + n_dur$failed, " vs ", n_tot$failed)
}
if (n_aes$polymorphic + n_dur$polymorphic != n_tot$polymorphic) {
  stop("polymorphic SNPs do not add up to the table's 24,545.")
}
non_failed <- n_tot$total - n_tot$failed
if (non_failed != 76102) stop("non-failed total is not 76,102.")
# the text states the same two sub-counts; they must partition the non-failed set
if (7237L + 68865L != non_failed) {
  stop("7,237 + 68,865 (the text's non-failed bases) is not 81,587 - 5,485.")
}

## --- the two percentage conventions must both reproduce ---------------------
# 1. failed % = failed / total (scored assays in the design).  Tolerance is the
#    printed precision: 703/7,940 = 8.854%, which the table prints as 8.8
#    (a 0.054 pp move, i.e. truncation at one decimal, not rounding); the other
#    four cells round correctly.  The guard is 0.055 pp, so a real drift in the
#    counts (which would move a cell by >= 0.1 pp) still fails.
calc_fail_pct <- 100 * tab$failed / tab$total
if (any(abs(calc_fail_pct - tab$failed_pct) > 0.055)) {
  stop("failed % are not failed/total within 0.055 pp: recomputed ",
       paste(round(calc_fail_pct, 2), collapse = "/"), " vs table ",
       paste(tab$failed_pct, collapse = "/"))
}
# 2. polymorphic % = polymorphic / NON-FAILED (the paper's footnote)
calc_poly_pct <- 100 * tab$polymorphic / (tab$total - tab$failed)
if (any(abs(calc_poly_pct - tab$polymorphic_pct) > 0.055)) {
  stop("polymorphic % are not polymorphic/non-failed within 0.055 pp: recomputed ",
       paste(round(calc_poly_pct, 2), collapse = "/"), " vs table ",
       paste(tab$polymorphic_pct, collapse = "/"))
}
# and the OTHER convention must NOT reproduce the print, or the footnote would
# be doing no work and the row below would be mislabelled
wrong <- 100 * tab$polymorphic / tab$total
if (any(abs(wrong - tab$polymorphic_pct) <= 0.055)) {
  stop("polymorphic % also reproduce on the full denominator; the footnote " ,
       "would be redundant - re-check the convention before landing.")
}

OUTCOME_DEF <- paste(
  "Array assay OUTCOME on a 288-accession durum (AABB) panel, as classified by",
  "the 90K iSelect pipeline: 'failed' = no usable call; 'polymorphic' = a call",
  "that segregates in the panel. Printed polymorphic % uses the NON-FAILED",
  "denominator (paper's footnote); the failed % uses all scored assays.",
  "Neither is a KASP conversion rate (this is an array, not KASP)."
)
PROVENANCE <- paste0(
  "BMC Genomics 20:268, Table 2 (article XML, not a supplement): ",
  "\"Number of failed, polymorphic and total SNPs in the assay\""
)

out <- data.frame(
  key = "durum2019",
  species = c("Triticum aestivum (AABBDD)", "Triticum durum (AABB)",
              "combined array"),
  ploidy = c("allo-hexaploid", "allo-tetraploid", "mixed"),
  assay = "Illumina 90K iSelect SNP array",
  snp_panel_origin = c("T. aestivum SNP catalogue", "T. durum SNP catalogue",
                       "both"),
  n_assays = tab$total,
  n_failed = tab$failed,
  n_polymorphic = tab$polymorphic,
  # the printed cells (kept verbatim: 703/7,940 prints as 8.8 but computes
  # 8.854) and the recomputed values, both carried so a reader sees the gap
  failed_pct_of_assays_printed = tab$failed_pct,
  polymorphic_pct_of_nonfailed_printed = tab$polymorphic_pct,
  failed_pct_of_assays = round(100 * tab$failed / tab$total, 4),
  polymorphic_pct_of_nonfailed =
    round(100 * tab$polymorphic / (tab$total - tab$failed), 4),
  polymorphic_pct_of_assays = round(100 * tab$polymorphic / tab$total, 2),
  # homeolog handling: the durum catalogue was built from homoeolog-resolved
  # SNPs (the paper's Table 1 splits its 7,940 chip SNPs into 4,461 hemi-SNPs
  # and 3,479 "simple" (single-dose, locus-specific) SNPs), so the SNP ORIGIN
  # axis was homeolog-aware even though the assay design itself was not.
  homeolog_handled = c("no", "yes", "mixed"),
  homeolog_detail = c(
    "A/D homoeologs are not separated by the assay; D-genome probes simply fail on an AABB panel",
    "durum SNPs were selected from a homoeolog-resolved transcriptome (4,461 hemi + 3,479 simple SNPs)",
    NA_character_
  ),
  outcome_definition = OUTCOME_DEF,
  evidence_tier = "measured",
  source_file = PROVENANCE,
  source_locator = "Table 2, rows 1-3; footnote on the polymorphic denominator",
  doi = DOI,
  stringsAsFactors = FALSE
)

## --- final re-derivation of every written cell ------------------------------
stopifnot(
  sum(out$n_assays[1:2]) == out$n_assays[3L],
    sum(out$n_failed[1:2]) == out$n_failed[3L],
    sum(out$n_polymorphic[1:2]) == out$n_polymorphic[3L],
    all(out$n_failed + out$n_polymorphic <= out$n_assays))
if (!identical(out$n_assays, tab$total)) stop("written n_assays drifted.")

fs::dir_create(out_dir)
readr::write_csv(out, fs::path(out_dir, "polyploid_durum2019_hexaploid_vs_tetraploid.csv"),
                 na = "")
cat("wrote polyploid_durum2019_hexaploid_vs_tetraploid.csv: failed",
    paste0(tab$failed, " (", tab$failed_pct, "%)", collapse = " / "),
    "; polymorphic ", paste0(tab$polymorphic, " (", tab$polymorphic_pct, "%)",
                             collapse = " / "), "\n")
