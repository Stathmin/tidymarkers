# =============================================================================
#  data-raw/sanger_bda_braf.R  --  a MEASURED Sanger floor, with its caller
# =============================================================================
#  Source (DOI is PROVENANCE; every number is re-read from the fetched files):
#
#    High sensitivity Sanger sequencing detection of BRAF mutations in
#    metastatic melanoma FFPE tissue specimens.  Sci Rep 2021.
#    DOI 10.1038/s41598-021-88391-5   EPMC PMC8079675
#    Fetched: .tmp/lit/PMC8079675/fulltext.xml   (Europe PMC fullTextXML)
#             .tmp/lit/PMC8079675/supplements/41598_2021_88391_MOESM1_ESM.pdf
#               (Europe PMC supplementaryFiles -> Supplementary Table S3)
#
#  WHY THIS MOVES THE SANGER ROW OFF THE 'textbook' TIER
#  inst/extdata/sanger_detection_floor.csv held only two restatements of "about
#  20%" that neither paper measured.  This paper reports a TITRATION (V600E
#  reference standard diluted to 10/5/3/1/0.5/0.2/0.1% VAF) AND a per-sample call
#  table on 23 clinical FFPE specimens.  Crucially the floor belongs to
#  Sanger *plus* a variant-enrichment caller (BDA + VarTrace), not to Sanger
#  chemistry alone - which is exactly why this package must not carry one
#  "Sanger detects X%" constant.
#
#  WHAT IS LANDED
#   * SANGER-S3 TABLE: 23 patients, one BDA-VAF row each, read from Table S3.
#     'BDA VAF <0.1' is kept as a BOUND (is_bound = TRUE), never coerced.
#   * The lowest POSITIVE call among repaired FFPE samples (the paper's stated
#     0.20% VAF, sample 129533).
#   * The paper's claim that the assay's reference-material LoD is 0.1% VAF and
#     that non-V600E LoDs were "no worse than 0.5% VAF" - both as tier
#     'measured', since they are the authors' own titration, not a restatement.
#
#  WHAT IS NOT LANDED
#   * Main-text Table 1 (IHC/ddPCR/NGS concordance) is published as an IMAGE in
#     this build; it could not be read, so no 4-method concordance number here.
#   * The "<20%" textbook rows stay in sanger_detection_floor.csv unchanged -
#     this file is an ADDITION, not a replacement.
#
#  ASSERTIONS
#   * Exactly 23 rows; 8 V600E, 6 non-V600E (V600K/V600R), 8 V600-WT,
#     1 FFPE-damage row - these are the counts the paper states in prose.
#   * The lowest positive VAF among the 8 V600E calls is 0.20 and belongs to
#     sample 129533, matching the abstract's "0.20% VAF from FFPE lymph node".
#   * Every V600-WT row carries the '<0.1' bound and no numeric value above it.
#
#  Run:  Rscript data-raw/sanger_bda_braf.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

xml_path <- fs::path(root, ".tmp", "lit", "PMC8079675", "fulltext.xml")
sup_pdf <- fs::path(
  root, ".tmp", "lit", "PMC8079675", "supplements",
  "41598_2021_88391_MOESM1_ESM.pdf"
)
for (p in c(xml_path, sup_pdf)) if (!file.exists(p)) stop("Missing ", p)

doc <- xml2::read_xml(xml_path)
if (!grepl("10.1038/s41598-021-88391-5", xml2::xml_text(doc), fixed = TRUE)) {
  stop("fulltext.xml does not mention the DOI")
}

# --- Table S3 from the supplementary PDF ------------------------------------
txt <- tempfile(fileext = ".txt")
if (system2("pdftotext", c("-layout", shQuote(sup_pdf), shQuote(txt))) != 0L) {
  stop("pdftotext failed on the supplementary PDF")
}
L <- readLines(txt, warn = FALSE)
# 'Sample ID' also appears in the earlier clinical table and in the contents,
# so anchor on the Section S4 heading and take the first header after it.
s4 <- grep("Section S4:", L, fixed = TRUE)
if (length(s4) != 1L) stop("Section S4 heading not found exactly once")
hdr <- grep("Sample ID", L)
hdr <- hdr[hdr > s4]
if (!length(hdr)) stop("Table S3 header not found in Section S4")
start <- hdr[1] + 2L # skip the wrapped 'dissection' header
rows <- L[start:(start + 40L)]
# one data row: Patient SampleID Type Macrodissect Repair Identity VAF
pat <- "^[[:space:]]*([0-9]+)[[:space:]]+([0-9]+)[[:space:]]+(SLNB|CLND)[[:space:]]+(Yes|No)[[:space:]]+(Yes|No)[[:space:]]+([A-Za-z0-9 ./+-]+?)[[:space:]]+(\\(?<?[0-9.]+%?\\)?)[[:space:]]*$"
m <- regmatches(rows, regexec(pat, rows))
m <- Filter(function(z) length(z) == 8L, m)
if (length(m) != 23L) stop("expected 23 rows in Table S3, parsed ", length(m))

d <- data.frame(
  patient = as.integer(vapply(m, `[`, "", 2L)),
  # keep as character: read.csv() would strip the leading "1" pattern of
  # 129533 into an integer, and the sample ids are labels, not numbers.
  sample_id = as.character(vapply(m, `[`, "", 3L)),
  tissue = vapply(m, `[`, "", 4L),
  macrodissected = vapply(m, `[`, "", 5L) == "Yes",
  ffpe_repaired = vapply(m, `[`, "", 6L) == "Yes",
  bda_identity = trimws(vapply(m, `[`, "", 7L)),
  raw = vapply(m, `[`, "", 8L),
  stringsAsFactors = FALSE
)

# The FFPE-damage row is parenthesised by the authors to mark it as NOT an
# actionable-mutation call; that distinction is preserved, not dropped.
d$is_ffpe_damage <- grepl("^\\(", d$raw)
d$is_bound <- grepl("<", d$raw)
d$bda_vaf_pct <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", d$raw)))
d$bda_vaf_pct[d$is_bound] <- NA_real_ # '<0.1' is a bound, not a value

if (!all(d$ffpe_repaired)) stop("every Table S3 sample is repaired; a non-repaired row appeared")
n_v600e <- sum(d$bda_identity == "V600E")
n_non <- sum(d$bda_identity %in% c("V600K", "V600R"))
n_wt <- sum(d$bda_identity == "V600-WT")
n_dmg <- sum(d$is_ffpe_damage)
# Table S3 splits the nine non-mutated samples into 8 explicit V600-WT rows and
# 1 row where the caller quantitated FFPE damage; the prose calls all nine
# "V600-WT".  Assert against the TABLE, and that the two agree on 9.
if (n_v600e != 8L || n_non != 6L || n_wt != 8L || n_dmg != 1L) {
  stop(
    "identity counts are ", n_v600e, "/", n_non, "/", n_wt, "/", n_dmg,
    ", expected 8 V600E / 6 non-V600E / 8 WT / 1 damage"
  )
}
if (n_wt + n_dmg != 9L) stop("WT + damage must be the nine samples the prose calls V600-WT")
if (!grepl("Nine samples were identified as BRAF V600-WT",
  paste(gsub("[[:space:]]+", " ", xml2::xml_text(
    xml2::xml_find_all(doc, "//*[local-name()='p']")
  )), collapse = " "),
  fixed = TRUE
)) {
  stop("the prose no longer says nine samples were V600-WT")
}
if (any(!is.na(d$bda_vaf_pct) & (d$bda_vaf_pct < 0 | d$bda_vaf_pct > 100))) {
  stop("a BDA VAF is outside 0-100")
}
if (any(d$is_bound != (d$bda_identity == "V600-WT"))) {
  stop("the '<0.1' bound must appear exactly on the V600-WT rows")
}

lowest <- d[which.min(replace(d$bda_vaf_pct, d$bda_identity != "V600E", NA)), ]
if (lowest$bda_vaf_pct != 0.20 || lowest$sample_id != "129533") {
  stop(
    "lowest V600E call is ", lowest$sample_id, " at ", lowest$bda_vaf_pct,
    "%, expected 129533 at 0.20%"
  )
}

# --- the authors' own titration statements ---------------------------------
paras <- gsub(
  "[[:space:]]+", " ",
  xml2::xml_text(xml2::xml_find_all(doc, "//*[local-name()='p']"))
)
abs0 <- paste(paras, collapse = " ")
if (!grepl("down to 0.20% VAF", abs0, fixed = TRUE)) {
  stop("abstract no longer states 'down to 0.20% VAF'")
}
ref_lod <- Filter(function(s) grepl("0.1% VAF limit of detection", s, fixed = TRUE), paras)
if (length(ref_lod) != 1L) stop("reference-material LoD sentence not found exactly once")
non_lod <- Filter(function(s) grepl("no worse than 0.5% VAF", s, fixed = TRUE), paras)
if (length(non_lod) != 1L) stop("non-V600E LoD sentence not found exactly once")

d$outcome_definition <- paste(
  "Outcome = variant allele fraction reported by a BDA-enriched Sanger assay",
  "(NuProbe VarTrace) for ONE patient sample; '<0.1' means no variant peak above",
  "the caller's floor and is kept as a bound, not a zero. The floor belongs to",
  "Sanger PLUS this enrichment/caller, not to Sanger chemistry alone."
)
d$evidence_tier <- "measured"
d$doi <- "10.1038/s41598-021-88391-5"
d$source_file <- "supplements/41598_2021_88391_MOESM1_ESM.pdf (Table S3)"
d$caller <- "BDA (Blocker Displacement Amplification) + VarTrace BRAF assay"
readr::write_csv(d, fs::path(root, "inst", "extdata", "sanger_bda_braf_calls.csv"))

lim <- data.frame(
  metric = c("detection_floor_pct", "detection_floor_pct", "detection_floor_pct"),
  value_pct = c(lowest$bda_vaf_pct, 0.1, 0.5),
  lower_bound = c(FALSE, FALSE, FALSE),
  scope = c(
    "clinical FFPE sample, lowest positive V600E call",
    "reference-material titration (V600E), stated assay LoD",
    "worst-case LoD across the six non-V600E spike-ins"
  ),
  method = c("BDA-enriched Sanger", "BDA-enriched Sanger", "BDA-enriched Sanger"),
  evidence_tier = "measured",
  usable_as_calibration = TRUE,
  doi = "10.1038/s41598-021-88391-5",
  source_file = c(
    "supplements/41598_2021_88391_MOESM1_ESM.pdf (Table S3)",
    "fulltext.xml", "fulltext.xml"
  ),
  outcome_definition = paste(
    "Outcome = variant allele fraction, for a Sanger assay that ENRICHES the",
    "variant before sequencing. This is NOT the population-Sanger floor; it is",
    "the floor of a specific enrichment-plus-caller pipeline, so it cannot be",
    "used as a Sanger-wide constant."
  ),
  stringsAsFactors = FALSE
)
if (!all(lim$usable_as_calibration == (lim$evidence_tier == "measured"))) {
  stop("usable_as_calibration disagrees with evidence_tier")
}
readr::write_csv(lim, fs::path(root, "inst", "extdata", "sanger_bda_braf_floor.csv"))

message(
  "sanger_bda_braf: ", nrow(d), " clinical calls (", n_v600e, " V600E, ", n_non,
  " non-V600E, ", n_wt, " WT, ", n_dmg, " damage); lowest positive ",
  lowest$bda_vaf_pct, "% at sample ", lowest$sample_id,
  "; reference LoD 0.1%, non-V600E <=0.5% -- assertions passed"
)
