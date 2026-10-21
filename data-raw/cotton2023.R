# =============================================================================
#  data-raw/cotton2023.R  --  allotetraploid cotton, per-assay KASP outcomes
# =============================================================================
#  Schoonmaker AN, Hulse-Kemp AM, Youngblood RC, et al.  "Detecting Cotton Leaf
#  Curl Virus Resistance Quantitative Trait Loci in Gossypium hirsutum and
#  iCottonQTL a New R/Shiny App to Streamline Genetic Mapping", Plants 12:1153
#  (2023), 10.3390/plants12051153, PMC10005503.
#
#  Why this source: upland cotton is an ALLOTETRAPLOID (AADD) with no D-genome
#  model organism to design against, and this is one of the few cotton papers
#  that prints a PER-ASSAY KASP outcome (Table 3: 11 assays, GOOD / BAD), so it
#  gives a real denominator and a real failure count instead of "we developed
#  N markers".  It never mentions homoeologs: the assays are screened on
#  cluster separation only.
#
#  It is also the family's only *negative* datapoint on the homeolog axis:
#  among the assays that worked are two designed from DIPLOID relatives
#  (i60979Gt, G. tomentosum; i64854Gm, G. mustelinum) and used, unre-designed,
#  on an AADD upper panel.  Those primer pairs amplified across the subgenome
#  boundary with no homeolog-aware design step, which is evidence that
#  mis-priming was not the binding problem in this system.
#
#  Fetch:
#    curl -sS -L -m 90 -o .tmp/lit/pmc10005503/fulltext.xml \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10005503/fullTextXML'
#
#  What this script refuses to do
#  ------------------------------
#  1. It does not recompute the paper's 72.72%.  8/11 = 72.727%, printed 72.72,
#     so the printed value is a truncation of the same ratio; the script
#     asserts that and stores the printer's own string, not a fresh rounding.
#  2. It does not call the one GOOD* assay "polymorphic".  Table 3 marks
#     i35622Gh "Might be best used as a dominant marker" - a usable assay with
#     a degraded call class - so the outcomes are kept as GOOD / GOOD* / BAD.
#  3. It does not claim the 11 were designed and 8 converted in the same
#     experiment: the text says the eleven were "tested ... on a new set of
#     samples obtained in the US", i.e. this is a re-screen of existing assays.
#
#  Run:  Rscript data-raw/cotton2023.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- Sys.getenv("COTTON2023_DIR", file.path(root, ".tmp", "lit", "pmc10005503"))
xml_path <- Sys.getenv("COTTON2023_XML", file.path(lit, "fulltext.xml"))
if (!file.exists(xml_path)) {
  stop("Fetch fulltext.xml into ", lit, " (see the header), or set COTTON2023_XML.")
}

DOI <- "10.3390/plants12051153"
out_dir <- fs::path(root, "inst", "extdata")

doc <- xml2::xml_ns_strip(xml2::read_xml(xml_path))
tabs <- xml2::xml_find_all(doc, ".//table-wrap")
cap <- vapply(tabs, function(t) paste(xml2::xml_text(
  xml2::xml_find_all(t, ".//caption")), collapse = " "), character(1))
hit <- grep("KASP assay screening results", cap, fixed = TRUE)
if (length(hit) != 1L) stop("Could not find exactly one Table 3 in the XML.")

rows <- xml2::xml_find_all(tabs[[hit]], ".//tr")
cells <- lapply(rows, function(r) {
  trimws(gsub("[[:space:]]+", " ",
              xml2::xml_text(xml2::xml_find_all(r, ".//td|.//th"))))
})
# drop the header row and any footnote rows (a single long cell)
cells <- Filter(function(r) length(r) >= 6L, cells)
# the header row shares the six-cell shape; an assay row's first cell is a marker
# id (i04503Gh / i60979Gt / i64854Gm), so the marker pattern is the discriminator
cells <- Filter(function(r) grepl("^i[0-9]+[A-Z][a-z]?$", r[[1L]]), cells)
if (length(cells) != 11L) {
  stop("Table 3 should carry 11 assays (the text says eleven were tested); got ",
       length(cells))
}

# columns: Marker | Cross | QTL | Linkage Group | Position (cM) | Good/Bad
tab <- lapply(cells, function(r) {
  list(marker = r[[1L]], cross = r[[2L]], qtl = r[[3L]],
       lg = r[[4L]], pos = as.numeric(r[[5L]]), verdict = r[[6L]])
})
marker <- vapply(tab, `[[`, character(1), "marker")
verdict <- vapply(tab, `[[`, character(1), "verdict")
if (anyDuplicated(marker)) stop("Duplicate marker IDs in Table 3.")

## --- the three outcome classes, read from the verdict cells -----------------
good <- grepl("^GOOD$", verdict)
good_star <- grepl("^GOOD\\s*\\*$", verdict)
bad <- grepl("^BAD$", verdict)
if (!all(good | good_star | bad)) {
  stop("Unrecognised Good/Bad cell in Table 3: ",
       paste(unique(verdict[!(good | good_star | bad)]), collapse = ", "))
}
n_good <- sum(good)
n_star <- sum(good_star)
n_bad <- sum(bad)
n_tested <- length(tab)

## --- the text's own headline, re-read from the fulltext ---------------------
xml <- paste(readLines(xml_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
txt <- gsub("[[:space:]]+", " ", gsub("<[^>]+>", " ", xml))
for (sent in c(
  "we tested eleven KASP marker assays",
  "found eight to be of sufficient quality",
  "seven were found to behave in a co-dominant fashion and only one was found to behave in a dominant fashion"
)) {
  if (!grepl(sent, txt, fixed = TRUE)) {
    stop("The fulltext no longer contains the sentence: ", sent)
  }
}

## --- assertions against the paper's own totals -----------------------------
if (n_tested != 11L) stop("Table 3 does not carry the eleven assays of the text.")
if (n_good + n_star != 8L) {
  stop("Table 3 gives ", n_good + n_star, " usable assays; the text says eight.")
}
if (n_good != 7L) stop("Table 3 gives ", n_good,
                        " plain GOOD assays; the text says seven co-dominant.")
if (n_star != 1L) stop("Table 3 gives ", n_star,
                        " GOOD* assays; the text says one dominant.")
if (n_bad != 3L) stop("Table 3 gives ", n_bad,
                       " BAD assays; 11 - 8 should be 3.")
if (n_good + n_star + n_bad != n_tested) stop("verdicts do not partition the 11.")
# the printed 72.72% is a truncation of 8/11, not a different ratio
printed_rate <- 72.72
if (abs(100 * (n_good + n_star) / n_tested - printed_rate) > 0.05 ||
    printed_rate > 100 * (n_good + n_star) / n_tested) {
  stop("The printed 72.72% is not 8/11 truncated to two decimals.")
}
# the two diploid-relatives markers must be GOOD, or the negative datapoint the
# report makes rests on a failed assay
rel <- grepl("Gt$|Gm$", marker)
if (!any(rel)) stop("The G. tomentosum / G. mustelinum markers are not in Table 3.")
if (!all(good[rel])) stop("A diploid-relative marker is not GOOD in Table 3.")

OUTCOME_DEF <- paste(
  "Per-assay KASP SCREEN done on a new US sample panel: GOOD = clusters",
  "separate sufficiently (X:X, X:Y, Y:Y); GOOD* = usable but best as a dominant",
  "marker; BAD = clusters do not separate. The screen is cluster separation,",
  "not agreement with a reference genotype, and it re-screens EXISTING assays,",
  "so it is not a SNP-to-KASP conversion count."
)
PROVENANCE <- paste0(
  "Plants 12:1153, Table 3 'KASP assay screening results for QTL-associated ",
  "and flanking markers' (", n_tested, " assay rows); text sentence ",
  "\"we tested eleven KASP marker assays ... found eight to be of sufficient ",
  "quality\""
)

out <- data.frame(
  key = "cotton2023",
  species = "Gossypium hirsutum (AADD)",
  ploidy = "allo-tetraploid",
  ploidy_detail = paste(
    "A and D subgenomes; the D subgenome has diploid congeners (G. raimondii),",
    "which is what allowed cross-species markers to be designed from them"
  ),
  homeolog_handled = "not stated",
  homeolog_detail = paste(
    "no homoeolog/subgenome screen is described; the screen is cluster",
    "separation only, on an AADD panel"
  ),
  assay = "KASP",
  marker,
  marker_origin = ifelse(rel, "designed from a diploid Gossypium relative",
                         "other"),
  outcome_class = ifelse(good, "GOOD", ifelse(good_star, "GOOD*", "BAD")),
  usable = good | good_star,
  outcome_definition = OUTCOME_DEF,
  evidence_tier = "measured",
  source_file = PROVENANCE,
  source_locator = "Table 3, rows 1-11 (Good/Bad column)",
  n_tested = n_tested,
  n_usable = n_good + n_star,
  printed_success_pct = printed_rate,
  computed_success_pct = round(100 * (n_good + n_star) / n_tested, 4),
  doi = DOI,
  stringsAsFactors = FALSE
)

## --- final re-derivation ----------------------------------------------------
if (sum(out$usable) != n_good + n_star) stop("usable flag disagrees with the cells.")
if (sum(out$outcome_class == "BAD") != n_tested - (n_good + n_star)) {
  stop("BAD count does not close the tested total.")
}
if (nrow(out) != 11L) stop("the written table does not have 11 rows.")

fs::dir_create(out_dir)
readr::write_csv(out, fs::path(out_dir, "polyploid_cotton2023_kasp_screen.csv"),
                 na = "")
cat("wrote polyploid_cotton2023_kasp_screen.csv:", n_tested, "tested,",
    n_good + n_star, "usable (", n_good, "co-dominant +", n_star, "dominant),",
    n_bad, "bad; diploid-relative markers:", paste(marker[rel], collapse = ", "),
    "\n")
