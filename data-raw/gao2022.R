# =============================================================================
#  data-raw/gao2022.R  --  a three-tier conversion chain in autopolyploid sugarcane
# =============================================================================
#  Gao Y, Zhou S, Huang Y, et al.  "Quantitative Trait Loci Mapping and
#  Development of KASP Marker Smut Screening Assay ... in Sugarcane
#  (Saccharum spp.)", Frontiers in Plant Science 12:796189 (2022),
#  10.3389/fpls.2021.796189, PMC8766830.
#
#  The value here is a *nested* outcome chain the package can carry as three
#  rows without inventing a rate:
#
#    77 SNPs selected for KASP
#     -> 49 "successfully genotyped in the two parents"
#     -> 24 "in line with genotypes in our genetic map"   (Supplementary Table 10)
#
#  Twenty-four is a countable file (sheet S10 of Table_1.XLSX), so that tier is
#  MEASURED.  The intermediate 49 is a sentence in the text and the workbook
#  has no per-assay record of the 28 siblings that were genotyped but did not
#  reproduce the map genotype, so it is carried with the sentence as its only
#  locator.  Only Supplementary Table 10 is claimed as a file-backed count.
#
#  Fetch:
#    curl -sS -L -m 90 -o .tmp/lit/gao2022/supplements.zip \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8766830/supplementaryFiles'
#    unzip -o .tmp/lit/gao2022/supplements.zip -d .tmp/lit/gao2022/supp
#    curl -sS -L -m 90 -o .tmp/lit/gao2022/fulltext.xml \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8766830/fullTextXML'
#  The per-marker docx tables the site advertises (SupplementaryTables 8-10)
#  are 404 on every URL tried; Table_1.XLSX inside the EPMC zip carries S1-S10.
#
#  What this script refuses to do
#  ------------------------------
#  1. It does not call the primers "PolyMarker-style".  The workbook gives the
#     two allele primers and the common primer but NOT which end is the SNP:
#     the allele-specific 3' base equals its AlleleX/AlleleY base in only 7 of
#     24 (X) and 6 of 24 (Y) rows, so the workbook is printed with the SNP at
#     the 5' end of the allele primers.  Noted, not used as an axis.
#  2. It does not compute a rate across the three tiers: 24/77 would mix
#     "genotyped in two parents" (49) with "reproduced the map" (24).
#
#  Run:  Rscript data-raw/gao2022.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
source("data-raw/lib/xlsx-read.R")

lit <- Sys.getenv("GAO_DIR", file.path(root, ".tmp", "lit", "gao2022"))
xlsx <- Sys.getenv("GAO_XLSX", file.path(lit, "supp", "Table_1.XLSX"))
xml_path <- Sys.getenv("GAO_XML", file.path(lit, "fulltext.xml"))
if (!file.exists(xlsx) || !file.exists(xml_path)) {
  stop("Fetch Table_1.XLSX and fulltext.xml into ", lit,
       " (see the header), or set GAO_XLSX / GAO_XML.")
}

DOI <- "10.3389/fpls.2021.796189"
out_dir <- fs::path(root, "inst", "extdata")

## --- Supplementary Table 10: one row per assay in line with the map ---------
sheets <- xlsx_sheets(xlsx)
if (!"S10" %in% names(sheets)) stop("Sheet S10 is not in the workbook.")
s10 <- sheets[["S10"]]
header <- paste(unlist(s10[[1L]]), collapse = " ")
if (!grepl("KASP", header)) stop("Sheet S10 does not look like Table S10: ", header)

cols <- c(id = "A", ax = "B", ay = "C", px = "D", py = "E", pc = "F")
cells <- s10[-(1:2)]                      # row 1 = title, row 2 = column names
get <- function(row, k) {
  v <- xcell(row, k)
  if (is.na(v)) "" else v
}
tab <- as.data.frame(lapply(cols, function(k) {
  vapply(cells, function(r) get(r, k), character(1))
}), stringsAsFactors = FALSE)
names(tab) <- names(cols)
tab <- tab[nzchar(tab$id), ]

if (nrow(tab) != 24L) {
  stop("Supplementary Table 10 lists ", nrow(tab),
       " assays; the text says 24 SNPs were in line with the genetic map.")
}
if (anyDuplicated(tab$id)) stop("Duplicate assay IDs in Table S10.")
if (!all(tab$ax %in% c("A", "C", "G", "T") & tab$ay %in% c("A", "C", "G", "T"))) {
  stop("Table S10 carries an allele base outside ACGT.")
}
if (any(tab$ax == tab$ay)) stop("An assay in Table S10 has identical alleles.")
if (!all(nzchar(tab$pc))) stop("A row of Table S10 has no common primer.")

## --- the text chain, re-read from the fulltext ------------------------------
xml <- paste(readLines(xml_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
txt <- gsub("[[:space:]]+", " ", gsub("<[^>]+>", " ", xml))
chain <- "Among the 77 SNPs selected for KASP assays, 49 were successfully genotyped in the two parents, and 24 SNPs were in line with genotypes in our genetic map"
if (!grepl(chain, txt, fixed = TRUE)) {
  stop("The 77 -> 49 -> 24 sentence is no longer in the fulltext.")
}
n_selected <- 77L
n_genotyped <- 49L
n_in_map <- 24L
if (nrow(tab) != n_in_map) stop("Table S10 count disagrees with the sentence.")
if (!(n_selected > n_genotyped && n_genotyped > n_in_map)) {
  stop("The conversion chain is not strictly decreasing.")
}

## --- how the workbook prints the allele-specific end (reported, not used) ---
as_3p_match <- c(
  X = sum(substr(tab$px, nchar(tab$px), nchar(tab$px)) == tab$ax),
  Y = sum(substr(tab$py, nchar(tab$py), nchar(tab$py)) == tab$ay)
)
as_5p_match <- c(
  X = sum(substr(tab$px, 1L, 1L) == tab$ax),
  Y = sum(substr(tab$py, 1L, 1L) == tab$ay)
)

OUTCOME_DEF <- paste(
  "Tier 1 'genotyped': clear fluorescence signal in BOTH parents of GT21 x",
  "ROC25. Tier 2 'in line with the map': the KASP genotype of the two parents",
  "reproduces the SNP genotype in the genetic map (Supplementary Table 10).",
  "Both are a two-parent assay outcome, not a population call rate."
)
PROVENANCE <- paste0(
  "Supplementary Table 10 (sheet S10 of Table_1.XLSX, ",
  nrow(tab), " assay rows); the 77 and 49 tiers exist only as the sentence ",
  "\"Among the 77 SNPs selected for KASP assays, 49 were successfully genotyped ",
  "in the two parents, and 24 SNPs were in line with genotypes in our genetic ",
  "map\" (Results)"
)

out <- data.frame(
  key = "gao2022",
  species = "Saccharum spp. (hybrid)",
  ploidy = "auto-polyploid",
  ploidy_detail = paste(
    "~8-10x (modern sugarcane hybrid); the paper works on a single biparental",
    "cross and does not filter subgenomes, so copy number is not modelled"
  ),
  homeolog_handled = "not stated",
  assay = "KASP",
  tier = c("selected", "genotyped in both parents", "in line with the map"),
  n = c(n_selected, n_genotyped, n_in_map),
  outcome_definition = OUTCOME_DEF,
  evidence_tier = c("weak", "weak", "measured"),
  source_file = PROVENANCE,
  source_locator = c("fulltext Results sentence", "fulltext Results sentence",
                     "sheet S10 rows 3-26"),
  as_allele_match_3prime = as.numeric(NA),
  as_allele_match_5prime = as.numeric(NA),
  doi = DOI,
  stringsAsFactors = FALSE
)
# the primer-end note belongs on the assay row (tier 3), where the workbook is
out$as_allele_match_3prime[[3L]] <- sum(as_3p_match)
out$as_allele_match_5prime[[3L]] <- sum(as_5p_match)

## --- final re-derivation ----------------------------------------------------
if (out$n[[3L]] != nrow(tab)) stop("tier 3 count is not the workbook row count.")
if (out$n[[1L]] - out$n[[2L]] != n_selected - n_genotyped) {
  stop("tier 1 - tier 2 is not the 77 - 49 gap.")
}
if (sum(as_3p_match) + sum(as_5p_match) < 1L) stop("no allele-end reading at all.")

fs::dir_create(out_dir)
readr::write_csv(out, fs::path(out_dir, "polyploid_gao2022_sugarcane_chain.csv"),
                 na = "")
cat("wrote polyploid_gao2022_sugarcane_chain.csv:", n_selected, "->",
    n_genotyped, "->", nrow(tab), "; allele at 3' in", sum(as_3p_match),
    "/48 allele primers, at 5' in", sum(as_5p_match), "/48\n")
