# =============================================================================
#  data-raw/gao2022_kasp_pvalues.R  --  the P-values of six used KASP assays
# =============================================================================
#  Gao Y, Zhou S, Huang Y, et al.  "Quantitative Trait Loci Mapping and
#  Development of KASP Marker Smut Screening Assay ... in Sugarcane
#  (Saccharum spp.)", Frontiers in Plant Science 12:796189 (2022),
#  10.3389/fpls.2021.796189, PMC8766830.
#
#  Table 7 ("KASP markers developed in this study and their P-values" in the
#  fulltext) is the only place the paper ties a named assay to the QTL it
#  screens and to the P-value of that association.  The sibling extract
#  `polyploid_gao2022_sugarcane_chain.csv` carries the *conversion* chain
#  (77 -> 49 -> 24); this file carries the six assays that made it into the
#  paper's smut-screening panel, with the QTL and P-value the paper prints.
#
#  Values are transcribed from Table 7; the stopifnot block is the transcription
#  guard, and when the fulltext is present the script re-reads the table and
#  asserts every marker / QTL / P-value triple it expects.
#
#  Fetch:
#    curl -sS -L -m 90 -o .tmp/lit/gao2022/fulltext.xml \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8766830/fullTextXML'
#
#  Run:  Rscript data-raw/gao2022_kasp_pvalues.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

marker <- c("A014653", "A014666", "A014669", "A014680", "A014692", "A014713")
qtl <- c("qSR51-1", "qSR23-1", "qSR23-1", "qSR59-2", "qSR59-2", "qSR59-2")
# the paper's own printed values, kept as text so 0.000387 never becomes 3.87e-4
p_text <- c("0.018357", "0.012546", "0.019141", "0.044136", "0.000387", "0.017852")
p_value <- as.numeric(p_text)

pvals <- data.frame(
  marker = marker,
  qtl = qtl,
  p_value = p_text,
  species = "Saccharum spp. (hybrid)",
  doi = "10.3389/fpls.2021.796189",
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(pvals) == 6L,
  !anyDuplicated(pvals$marker),
  # all six associations are significant at the paper's own 0.05 threshold
  all(p_value > 0 & p_value < 0.05),
  # the panel is not one QTL repeated six times
  length(unique(qtl)) == 3L,
  all(tabulate(match(qtl, unique(qtl))) == c(1L, 2L, 3L)),
  # the lone sub-1e-3 association is the one the paper highlights
  which.min(p_value) == 5L,
  marker[[which.min(p_value)]] == "A014692"
)

## --- re-read Table 7 from the fulltext when it is on disk -------------------
xml_path <- Sys.getenv("GAO_XML", file.path(root, ".tmp", "lit", "gao2022", "fulltext.xml"))
if (file.exists(xml_path)) {
  xml <- paste(readLines(xml_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  txt <- gsub("[[:space:]]+", " ", gsub("<[^>]+>", " ", xml))
  for (i in seq_len(nrow(pvals))) {
    triple <- paste(pvals$marker[[i]], pvals$qtl[[i]], pvals$p_value[[i]])
    if (!grepl(triple, txt, fixed = TRUE)) {
      stop("Table 7 no longer prints: ", triple)
    }
  }
  cat("fulltext re-read: all", nrow(pvals), "marker/QTL/P triples present\n")
} else {
  cat("fulltext not found at", xml_path, "-- transcription guard only\n")
}

readr::write_csv(pvals, "inst/extdata/gao2022_kasp_marker_pvalues.csv")
cat(
  "gao2022_kasp_marker_pvalues:", nrow(pvals), "assays on",
  length(unique(pvals$qtl)), "QTLs;",
  sprintf("P %s-%s\n", pvals$p_value[[which.min(p_value)]], pvals$p_value[[which.max(p_value)]])
)
