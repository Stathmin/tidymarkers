# =============================================================================
#  data-raw/stutter_woerner2017.R  --  -1 stutter ratio vs repeat length
# =============================================================================
#  Source (identifier = PROVENANCE; numbers come from the file):
#
#    "Flanking Variation Influences Rates of Stutter in Simple Repeats".
#    Woerner AE, King JL, Budowle B.  Genes (Basel) 2017;8(11):329.
#    DOI 10.3390/genes8110329.  PMID 29149052.  PMCID PMC5704242.
#
#  Fetched: .tmp/lit/PMC5704242/fulltext.xml (EPMC fullTextXML, 137192 bytes,
#  Crossref confirms title/volume/pages; Unpaywall reports the version of
#  record as closed, and EPMC holds it as an author manuscript).
#
#  CORRECTION TO THE BRIEF -- READ THIS
#  ------------------------------------
#  The brief asked for 10.1016/j.fsigen.2018.04.003, "Stutter analysis of
#  complex STR MPS data" (Vilsen et al., FSI Genetics 2018;35:107-112).  A
#  prior session cached the wrong article under the PMC id: EPMC search finds
#  one hit for that title (PMID 29702330) with isOpenAccess=N and NO PMCID, so
#  that paper is CLOSED and was not fetched.  The file at PMC5704242 is a
#  DIFFERENT, open paper (Woerner 2017), and it is the one whose tables are
#  parsed below.  Nothing here may be attributed to Vilsen 2018.
#
#  Two other stutter sources in the brief are also closed and were not bypassed:
#    10.1520/jfs2003312           (J Forensic Sci 2004)  dl.astm.org -> HTTP 403
#    10.1016/j.fsigen.2011.02.001 (FSI Genetics 2012)   Elsevier -> meta-refresh
#
#  WHY THIS PAPER STILL MAKES A GOOD SOURCE
#  ----------------------------------------
#  It publishes a *model* of stutter versus repeat length rather than a bare
#  per-locus percentage: robust linear regressions of the stutter ratio (SR)
#  on the longest uninterrupted stretch (LUS) or the parental allele length
#  (PAL), one pair of flanking haplotypes per locus, 10 loci (Tables 3, 4).
#
#  WHAT THIS LANDS  (inst/extdata/stutter_woerner2017.csv)
#  ----------------------------------------------------
#  Tables 3 and 4 of the paper, one row per locus per model:
#    - intercepts for the two flanking haplotypes (and the p-value for their
#      equality);
#    - the SLOPE of SR on the repeat-length measure (PAL for Table 3, LUS for
#      Table 4).  A slope is "fraction of stutter gained per extra repeat
#      unit", which is exactly the shape the package needs, and it is what the
#      paper actually tabulates;
#    - D5S818 is the only locus with a significandy better interaction model,
#      so it alone carries two slopes (0.007 / 0.009);
#    - VIF, R2, and the coefficient on the number of correct reads (C).
#
#  MANDATORY SCOPE FLAG (also carried in the CSV, on every row)
#  ------------------------------------------------------------
#  This is FORENSIC HUMAN STR data.  It can validate the *shape* of a stutter
#  response (more repeat units -> more stutter, and inter-haplotype spread),
#  and it must NOT be read as plant SSR stutter behaviour: no plant locus is
#  measured here, the repeat motifs are human autosomal, and the measurement
#  is MPS read counts (S-1 / C), not capillary peak areas.
#
#  Run:  Rscript data-raw/stutter_mps2018.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

lit <- fs::path(root, ".tmp", "lit", "PMC5704242")
xml_path <- fs::path(lit, "fulltext.xml")

DOI <- "10.3390/genes8110329"
PMCID <- "PMC5704242"
PMID <- "29149052"

xml_path <- fs::path(lit, "fulltext.xml")
for (p in xml_path) if (!file.exists(p)) stop("Missing ", p)

doc <- xml2::read_xml(xml_path)
title <- xml2::xml_text(xml2::xml_find_first(doc, "//*[local-name()='article-title']"))
if (!grepl("Flanking Variation Influences Rates of Stutter in Simple Repeats", title, fixed = TRUE)) {
  stop("Not the Woerner stutter paper: ", title)
}
ids <- xml2::xml_text(xml2::xml_find_all(doc, "//*[local-name()='article-id']"))
if (!DOI %in% ids) stop("XML does not carry DOI ", DOI, "; found: ", paste(ids, collapse = ", "))
if (!PMID %in% ids) stop("XML does not carry PMID ", PMID)
# the journal must be Genes, not FSI Genetics: this is the file-identity check
# that the earlier session missed.
if (!grepl("Genes", xml2::xml_text(xml2::xml_find_first(doc, "//*[local-name()='journal-title']")), fixed = TRUE)) {
  stop("Fetched XML is not from Genes")
}

## --- the tables -------------------------------------------------------------
tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
if (length(tw) != 4L) stop("Expected 4 table-wraps, got ", length(tw))

read_table <- function(node) {
  rows <- xml2::xml_find_all(node, ".//*[local-name()='tr']")
  # headers are split across <italic>/<sup> children ("R 2", "p -Value",
  # "10 −06"), so ALL whitespace is collapsed out of the header only.  Cells
  # are left with single spaces, because "5.7 × 10 −06" needs them.
  cells <- lapply(rows, function(r) {
    vapply(xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']"),
           function(c) gsub("[[:space:]]+", " ",
                            trimws(paste(xml2::xml_text(xml2::xml_find_all(c, ".//text()")), collapse = " "))),
           character(1L))
  })
  hdr <- gsub("[[:space:]]+", "", cells[[1L]])
  if (any(!nzchar(hdr))) stop("A table header cell is empty: ", paste(hdr, collapse = "|"))
  body <- cells[-1L]
  if (length(unique(vapply(body, length, integer(1L)))) != 1L) {
    stop("Ragged table: rows have " , paste(unique(vapply(body, length, integer(1L))), collapse = "/"))
  }
  if (any(vapply(body, length, integer(1L)) != length(hdr))) stop("Row width differs from header width")
  # build BY COLUMN: do.call(rbind) drops names when rows are named vectors
  m <- matrix(unlist(body, use.names = FALSE), ncol = length(hdr), byrow = TRUE,
              dimnames = list(NULL, hdr))
  m
}

# labels/captions: the caption of each table names the model, so it is the
# caption ("... as a function of the PAL") that decides which slope column is
# which -- not the position of the table.
caption <- function(node) {
  paste(xml2::xml_text(xml2::xml_find_all(
    xml2::xml_find_first(node, ".//*[local-name()='caption']"), ".//text()")), collapse = " ")
}
cap3 <- caption(tw[[3L]])
cap4 <- caption(tw[[4L]])
if (!grepl("as a function of the PAL", cap3, fixed = TRUE)) stop("Table 3 caption is not the PAL model: ", cap3)
if (!grepl("as a function of the LUS", cap4, fixed = TRUE)) stop("Table 4 caption is not the LUS model: ", cap4)
t3 <- read_table(tw[[3L]])
t4 <- read_table(tw[[4L]])

# --- numbers are printed with a Unicode minus and "x10-super" ---------------
#  "−0.031", "5.7 × 10−06".  Parsing must be explicit; as.numeric() silently
#  returns NA for these.
num <- function(x) {
  x <- gsub("\u2212", "-", x, fixed = TRUE)          # U+2212 MINUS SIGN
  x <- gsub("\u00d7", "x", x, fixed = TRUE)          # U+00D7 MULTIPLICATION
  x <- gsub("[[:space:]]+", "", x)
  if (identical(x, "--") || !nzchar(x)) return(NA_real_)
  if (grepl("^-?[0-9.]+x10-?[0-9]+$", x)) {
    parts <- strsplit(x, "x10", fixed = TRUE)[[1L]]
    return(as.numeric(parts[1L]) * 10^as.numeric(parts[2L]))
  }
  as.numeric(x)
}

build <- function(m, measure) {
  cols <- colnames(m)
  if (is.null(cols) || !length(cols)) stop("Table lost its headers")
  pick <- function(prefix) {
    hit <- cols[cols == prefix]
    if (length(hit) != 1L) stop("Column '", prefix, "' is not exactly one header; headers are: ",
                               paste(cols, collapse = " | "))
    hit
  }
  want <- paste0(measure, "Slope,Haplotype1")
  if (!want %in% cols) stop("Header ", want, " absent: ", paste(cols, collapse = ", "))
  out <- data.frame(
    locus = m[, pick("Locus")],
    measure = measure,
    max_vif = vapply(m[, pick("MaxVIF")], num, numeric(1L)),
    r2 = vapply(m[, pick("R2")], num, numeric(1L)),
    intercept_hap1 = vapply(m[, pick("Intercept,Haplotype1")], num, numeric(1L)),
    intercept_hap2 = vapply(m[, pick("Intercept,Haplotype2")], num, numeric(1L)),
    p_same_intercept = vapply(m[, pick("p-Value(SameIntercept)")], num, numeric(1L)),
    slope_hap1 = vapply(m[, want], num, numeric(1L)),
    slope_hap2 = vapply(m[, pick(paste0(measure, "Slope,Haplotype2"))], num, numeric(1L)),
    p_same_slope = vapply(m[, pick("p-Value(SameSlope)")], num, numeric(1L)),
    slope_correct = vapply(m[, pick("Slope,Correct")], num, numeric(1L)),
    p_correct = vapply(m[, pick("p-Value(Correct)")], num, numeric(1L)),
    stringsAsFactors = FALSE
  )
  out
}

d3 <- build(t3, "PAL")
d4 <- build(t4, "LUS")
out <- rbind(d3, d4)

## --- assertions -------------------------------------------------------------
# 1. 10 loci, the paper's own count ("10 loci were identified that had at
#    least two ... flanking haplotypes")
if (nrow(out) != 20L) stop("Expected 10 loci x 2 measures = 20 rows, got ", nrow(out))
if (length(unique(out$locus)) != 10L) stop("Expected 10 loci")
if (!setequal(d3$locus, d4$locus)) stop("Tables 3 and 4 do not cover the same loci")

# 2. no NA survived the Unicode parse: a silent NA is the failure mode here
for (col in c("max_vif", "r2", "intercept_hap1", "intercept_hap2",
              "p_same_intercept", "slope_hap1", "slope_correct", "p_correct")) {
  if (anyNA(out[[col]])) {
    stop("Unparsed value in ", col, " for: ", paste(out$locus[is.na(out[[col]])], collapse = ", "))
  }
}

# 3. slopes are physically bounded: SR rises by < 1 per repeat unit
if (any(out$slope_hap1 < 0 | out$slope_hap1 > 1)) stop("Slope outside 0..1")
if (any(abs(out$intercept_hap1) > 1)) stop("Stutter intercept outside -1..1")

# 4. p-values are probabilities where a test was run; NA marks "--" (the
#    non-interaction model, which has no slope test).
for (col in c("p_same_intercept", "p_correct")) {
  v <- out[[col]]
  if (any(!is.na(v) & (v < 0 | v > 1))) stop(col, " is not a probability")
}
v <- out$p_same_slope
if (any(!is.na(v) & (v < 0 | v > 1))) stop("p_same_slope is not a probability")
if (any(is.na(v) != (out$locus != "D5S818"))) {
  stop("p_same_slope is only printed for D5S818, but NA pattern does not match")
}

# 5. D5S818 is the ONLY locus with a second slope (the paper: "Only one locus,
#    D5S818, showed evidence for different intercepts and different slopes").
two_slopes <- out$locus[!is.na(out$slope_hap2)]
if (!identical(unique(two_slopes), "D5S818")) {
  stop("Expected D5S818 alone to carry two slopes, got: ", paste(unique(two_slopes), collapse = ", "))
}
if (any(out$slope_hap2[out$locus == "D5S818"] != 0.009)) stop("D5S818 hap2 slope is not 0.009")

# 6. slopes agree within a locus across the two measures to within the printed
#    resolution, EXCEPT where the flank matters (the paper's point).  We assert
#    only the invariant part: both measures are positive at every locus.
if (any(d3$slope_hap1 <= 0) || any(d4$slope_hap1 <= 0)) stop("A stutter slope is not positive")

# 7. the paper's own worked example, in the Table 3 foot note:
#    "the two flanking haplotypes in the locus D9S1122 have difference in
#     intercepts that is ~0.040 ... roughly equals the effect of having an
#     additional three repeat units (0.013 x 3)"
#    Checked to the printed resolution: the gap is 0.040 and the slope x 3 is
#    0.039, and the paper calls that "roughly equal", so the tolerance is the
#    last printed digit, not zero.
i <- match("D9S1122", d3$locus)
diff_i <- abs(d3$intercept_hap1[i] - d3$intercept_hap2[i])
if (abs(diff_i - 0.04) > 0.001) stop("D9S1122 intercept gap is not ~0.040 (got ", round(diff_i, 3L), ")")
three_units <- d3$slope_hap1[i] * 3
if (abs(three_units - diff_i) > 0.002) {
  stop("D9S1122: 3 x slope (", round(three_units, 3L), ") is not ~equal to the intercept gap (",
       round(diff_i, 3L), ")")
}

out$scope_flag <- paste(
  "FORENSIC HUMAN STR (autosomal, MPS read-count stutter ratio).",
  "Validates stutter SHAPE only; not plant SSR behaviour."
)
out$outcome_definition <- paste(
  "slope = change in stutter ratio (SR = S-1 / C) per additional repeat unit of",
  "the measure named in `measure` (PAL = parental allele length, LUS = longest",
  "uninterrupted stretch), from a robust linear model SR ~ measure + flank + C"
)
out$measurement_conditions <- paste(
  "MPS read counts, 10 loci with exactly two flanking haplotypes and",
  ">=25 haplotype calls each; robust linear regression (lmrob), R 3.4.2"
)
out$evidence_tier <- "strong"
out$evidence_basis <- "Tables 3 and 4 of the paper, re-read and re-parsed here"
out$source_file <- "PMC5704242/fulltext.xml#table-wrap[3,4]"
out$doi <- DOI
out$pmcid <- PMCID

out <- out[order(out$measure, out$locus), ]
utils::write.csv(out, fs::path(root, "inst", "extdata", "stutter_woerner2017.csv"),
                 row.names = FALSE, na = "")
cat("stutter_woerner2017.csv:", nrow(out), "rows;",
    "slope range", min(out$slope_hap1), "-", max(out$slope_hap1),
    "per repeat unit; D9S1122 gap", round(diff_i, 3), "= 3 x", d3$slope_hap1[i], "\n")
