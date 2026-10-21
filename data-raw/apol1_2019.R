# =============================================================================
#  data-raw/apol1_2019.R  --  a TaqMan assay that FAILS on one allele: the
#                            deletion allele of the G2 haplotype
# =============================================================================
#  "Taqman Assay for Genotyping CKD-Associated APOL1 SNP rs60910145: A
#  Cautionary Note". Kidney Int Rep 2019. doi:10.1016/j.ekir.2018.08.003
#  PMCID PMC6308993 (OA).
#
#  WHY THIS FILE EXISTS
#  R/marker_taqman.R assumes the 5'-nuclease probe reports the allele it is
#  designed for, with a mismatch-dG gate as the only failure mode it models.
#  Here a PREDESIGNED TaqMan assay fails outright (no amplification) on one
#  haplotype and mis-calls another, on 1,110 samples.  That is an observed
#  failure MODE of the readout -- not a mismatch-discrimination measurement.
#
#  WHAT IT DOES NOT GIVE
#  * No dCt, no signal ratio, and no probe sequence: the deformation is a 6-bp
#    deletion that removes bases under the probe, and the assay id is the only
#    probe information the paper prints.
#  * The percentages are printed to 0.1% and are NOT recomputed from the 1,110
#    samples (the paper never prints the per-genotype counts).
#
#  ASSERTIONS (any failure -> stop(), nothing written)
#   (C1) Table 2 parses to 8 genotype rows; the custom-assay column is
#        T/GT/G (one of the three printed strings), never empty.
#   (C2) The predesigned-assay column carries the word "No amplification" for
#        exactly the G2,G2 row, and "Ga" (the wrong call) for G1GM,G2.
#   (C3) The eight frequencies sum to 99.9, exactly as printed -- a 0.1%
#        rounding shortfall, asserted as printed and NOT rescaled.
#   (C4) The text states 7.3% failed-or-wrong = 1.4% no-amplification + 5.9%
#        mis-typed, and 1,110 samples.
#   (C5) Row 1 of Table 1 says the predesigned assay does not amplify the G2
#        haplotype ("No" in the amplification column), and only that row does.
#
#  Run:  Rscript data-raw/apol1_2019.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- file.path(root, ".tmp", "lit", "apol1_2019")
xml_p <- file.path(lit, "fulltext.xml")
if (!file.exists(xml_p)) stop("Missing ", xml_p)
out_dir <- file.path(root, "inst", "extdata")

DOI <- local({
  # DOI from the fetched JATS, not from memory
  d <- xml2::read_xml(xml_p)
  v <- xml2::xml_text(xml2::xml_find_all(d, "//article-id[@pub-id-type='doi']"))
  if (!length(v)) stop("apol1_2019: no DOI in the JATS")
  v[[1L]]
})
PMCID <- "PMC6308993"

doc <- xml2::read_xml(xml_p)
x <- gsub("\u00a0", " ", xml2::xml_text(doc), fixed = TRUE)
x <- gsub("[[:space:]]+", " ", x)

problems <- character(0)
need <- function(s, what) {
  if (!grepl(gsub("[[:space:]]+", " ", s), x, fixed = TRUE)) {
    problems <<- c(problems, paste0("verbatim text not found (", what, "): ",
                                    substr(s, 1, 70)))
  }
  invisible(TRUE)
}

cells <- function(tbl) {
  rows <- lapply(xml2::xml_find_all(tbl, ".//table//tr"), function(r) {
    trimws(xml2::xml_text(xml2::xml_find_all(r, "./td|./th")))
  })
  rows[vapply(rows, length, 0L) >= 2L]
}

## -- C1/C2: Table 2 ------------------------------------------------------------
# the label uses a non-breaking space ("Table\u00a02"): find it in R, not XPath
find_table <- function(d, n) {
  tw <- xml2::xml_find_all(d, "//table-wrap")
  nm <- gsub("\u00a0", " ", xml2::xml_text(xml2::xml_find_first(tw, "label")),
             fixed = TRUE)
  tw[!is.na(nm) & startsWith(trimws(nm), paste0("Table ", n))]
}
t2 <- find_table(doc, 2L)
if (length(t2) != 1L) stop("apol1_2019: expected exactly one Table 2")
r2 <- cells(t2)
# the genotype rows start with a genotype pair; the two summary rows of the
# header block do not
geno <- r2[vapply(r2, function(v) grepl("^[A-Za-z0-9+]+,[A-Za-z0-9+]+$", v[[1L]]) &&
                    length(v) >= 4L, FALSE)]
if (length(geno) != 8L) {
  stop("apol1_2019: Table 2 parsed ", length(geno), " genotype rows, expected 8")
}
tab2 <- data.frame(
  genotype = vapply(geno, `[[`, "", 1L),
  custom_call = vapply(geno, `[[`, "", 2L),
  predesigned_call = vapply(geno, `[[`, "", 3L),
  freq_pct = as.numeric(vapply(geno, `[[`, "", 4L)),
  stringsAsFactors = FALSE)
tab2$genotype <- gsub("[[:space:]]", "", tab2$genotype)
if (any(!tab2$custom_call %in% c("T", "G", "GT"))) {
  problems <- c(problems, "a custom-assay call is not one of T/GT/G")
}
no_amp <- tab2$genotype[grepl("No amplification", tab2$predesigned_call, fixed = TRUE)]
if (!identical(no_amp, "G2,G2")) {
  problems <- c(problems, paste0("no-amplification row(s): ",
                                 paste(no_amp, collapse = ", ")))
}
wrong <- tab2$genotype[tab2$predesigned_call == "Ga"]
if (!identical(wrong, "G1GM,G2")) {
  problems <- c(problems, paste0("mis-called row(s): ", paste(wrong, collapse = ", ")))
}

## -- C3: the frequencies sum to 100 -------------------------------------------
freq_sum <- round(sum(tab2$freq_pct), 1L)
# as printed the eight rows sum to 99.9, not 100.0: a rounding artifact of the
# papers own 0.1% values, recorded as printed and never rescaled here
if (!isTRUE(freq_sum == 99.9)) {
  problems <- c(problems, paste0("Table 2 frequencies sum to ", freq_sum,
                                 ", not the printed 99.9"))
}

## -- C4: the text's own decomposition -----------------------------------------
need("1110 African American samples", "1,110 samples")
need("7.3% of the rs60910145 genotypes", "7.3% failed/wrong")
need("1.4% of the samples tested had a G2/G2 genotype", "1.4% no amplification")
need("5.9% of the samples had the G1", "5.9% mis-typed")
n_samples <- 1110L
pct_fail <- tab2$freq_pct[tab2$genotype == "G2,G2"]
pct_wrong <- tab2$freq_pct[tab2$genotype == "G1GM,G2"]
if (abs((pct_fail + pct_wrong) - 7.3) > 0.05) {
  problems <- c(problems, paste0("1.4 + 5.9 = ", pct_fail + pct_wrong, ", not the 7.3% printed"))
}

## -- C5: Table 1's amplification column ---------------------------------------
t1 <- find_table(doc, 1L)
if (length(t1) != 1L) stop("apol1_2019: expected exactly one Table 1")
r1 <- cells(t1)
hap <- r1[vapply(r1, function(v) grepl("^G[0-9]", v[[1L]]) && length(v) >= 4L, FALSE)]
if (length(hap) != 4L) {
  stop("apol1_2019: Table 1 parsed ", length(hap), " haplotype rows, expected 4")
}
hap_no <- vapply(hap, function(v) v[[length(v)]], "")
hap_id <- vapply(hap, function(v) v[[1L]], "")
if (!identical(hap_id[hap_no == "No"], "G2")) {
  problems <- c(problems, paste0("Table 1 non-amplifying haplotype(s): ",
                                 paste(hap_id[hap_no == "No"], collapse = ", ")))
}
if (!all(hap_no[hap_id != "G2"] == "Yes")) {
  problems <- c(problems, "a non-G2 haplotype is marked as not amplifying")
}

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("apol1_2019: assertion failed; nothing written")
}

## -- write --------------------------------------------------------------------
out <- tab2
out$assay = "rs60910145 predesigned (C__89555688_10) vs custom (AHWR1JA)"
out$outcome_definition <- paste(
  "observed genotype-call outcome of a PREDESIGNED TaqMan SNP assay against",
  "the custom assay on 1,110 African American samples: the deletion (G2)",
  "haplotype removes bases under the probe, so G2/G2 gives no amplification",
  "and G1GM/G2 is mis-typed; a failure MODE, not a mismatch discrimination")
out$evidence_tier <- "moderate"
out$source_file <- "fulltext.xml :: Table 2 (+ text: 1,110 samples, 7.3%)"
out$doi <- DOI
out$pmcid <- PMCID
out$n_samples <- n_samples

write.csv(out[, c("genotype", "custom_call", "predesigned_call", "freq_pct",
                  "assay", "n_samples", "outcome_definition", "evidence_tier",
                  "source_file", "doi", "pmcid")],
          file.path(out_dir, "taqman_apol1_2019_failure.csv"),
          row.names = FALSE, na = "")

message(
  "apol1_2019 (", DOI, "): ", nrow(out), " genotype rows, frequencies sum to ",
  freq_sum, "%; no amplification ", pct_fail, "% + mis-typed ", pct_wrong,
  "% = ", pct_fail + pct_wrong, "% (text 7.3%); n = ", n_samples
)
