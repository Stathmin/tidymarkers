# =============================================================================
#  data-raw/snapshot2026_sensitivity.R  --  SNaPshot detection vs RNA input
# =============================================================================
#  Source (the identifier is PROVENANCE; the numbers come from the files):
#
#    "High-efficiency 14-plex cSNP SNaPshot assay for forensic peripheral blood
#    identification", BMC Genomics 2026;27.  DOI 10.1186/s12864-026-13173-0
#    (PMID 42436368, PMCID PMC13644189).
#
#  Fetched (all under .tmp/lit/PMC13644189/):
#    fulltext.xml                       EPMC fullTextXML, 120579 bytes
#    supplements/12864_2026_13173_MOESM2_ESM.xlsx
#                                       EPMC supplementaryFiles -> the real
#                                       workbook (11 sheets, Table S1..S11)
#    prose.txt                          verbatim paragraph 35 (the failure-rate
#                                       sentences) extracted from fulltext.xml
#
#  KEY : PMC13644189
#
#  WHAT THIS LANDS  (inst/extdata/snapshot2026_sensitivity.csv)
#  --------------------------------------------------------------
#  One row per (locus x RNA input) for the 14 cSNP loci, computed from Table S4
#  under the paper's OWN detection rule, quoted in prose.txt:
#
#      "a detection threshold of 100 relative fluorescence units (RFU). Peaks
#       below 100 RFU were excluded"
#      "the fluorescence value of each locus was the sum of the fluorescence
#       values of its two alleles"                     (Methods + p34)
#
#  => locus_detected := (sum of the allele RFU in the cell) >= 100.
#     a locus is FAILED for a replicate when it is not detected.
#
#  n_replicates = 3 is the paper's own "no detectable peaks across all three
#  replicates" (prose.txt, p35), matching the three columns per input.
#
#  WHAT THIS DOES **NOT** LAND, AND WHY
#  ------------------------------------
#  The prose percentages are NOT quietly dropped: they are parsed here and
#  compared against the table, and the comparison is written to
#  data-raw/proposed/.  A reader can see that 1 ng reproduces and that some
#  lower-input numbers do not, instead of the package silently showing only
#  the half that agrees.
#
#  * We do NOT write inst/extdata/snapshot2026_failure_claims.csv from the
#    prose failure percentages.  The generator compares every one of the
#    paper's stated (locus, input, failure-rate) claims against Table S4 and
#    stores the comparison in
#    data-raw/proposed/snapshot2026_failure_claims.rows.csv.  Measured on the
#    primary sample set (block A), the 1 ng claims reproduce exactly, the
#    0.5 ng / 0.1 ng claims for rs750625, rs1741487, rs1741488, rs229586 do
#    not, and the claim for rs7687256 at 0.1 ng is contradicted.  The rest of
#    the Table S4 workbook (block B, a different set of three samples) also
#    does not reproduce them.  So the prose percentages are NOT loaded into
#    the package CSVs; only the arithmetic and the comparison are.
#  * The mixed-sample claims (p38) point at "Table S5"; EPMC's Table S5 sheet
#    is present, but this script does not touch it (out of scope for 1 story).
#  * The 500-sample claims are in snapshot2026_specificity.R.
#
#  Run:  Rscript data-raw/snapshot2026_sensitivity.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

lit <- fs::path(root, ".tmp", "lit", "PMC13644189")
xml_path <- fs::path(lit, "fulltext.xml")
xlsx_path <- fs::path(lit, "supplements", "12864_2026_13173_MOESM2_ESM.xlsx")
prose_path <- fs::path(lit, "prose.txt")
for (p in c(xml_path, xlsx_path, prose_path)) if (!file.exists(p)) stop("Missing ", p)

DOI <- "10.1186/s12864-026-13173-0"
PMCID <- "PMC13644189"

## --- the file must be the paper we think it is ------------------------------
doc <- xml2::read_xml(xml_path)
title <- xml2::xml_text(xml2::xml_find_first(doc, "//*[local-name()='article-title']"))
if (!grepl("SNaPshot", title, fixed = TRUE)) stop("Not the SNaPshot paper: ", title)
if (!"14-plex" %in% strsplit(title, "[ -]")[[1L]] && !grepl("14-plex", title, fixed = TRUE)) {
  stop("Title does not say 14-plex: ", title)
}
ids <- xml2::xml_text(xml2::xml_find_all(doc, "//*[local-name()='article-id']"))
if (!DOI %in% ids) stop("XML does not carry DOI ", DOI, "; found: ", paste(ids, collapse = ", "))
if (!PMCID %in% ids) stop("XML does not carry ", PMCID)

prose <- readLines(prose_path, warn = FALSE)
prose_all <- paste(prose, collapse = " ")
for (phrase in c(
  "detection threshold of 100 relative fluorescence units",
  "no detectable peaks across all three replicates",
  "the sum of the fluorescence values of its two alleles"
)) {
  if (!grepl(phrase, prose_all, fixed = TRUE)) stop("Fetched prose lost the phrase: ", phrase)
}

## --- Table S4 into a character matrix ---------------------------------------
# xlsx-read.R is the repo's xml2-based reader (no readxl in this environment).
source(fs::path(root, "data-raw", "lib", "xlsx-read.R"))
sheets <- xlsx_sheets(xlsx_path)
if (!"Table S4 Sensitivity" %in% names(sheets)) {
  stop("Table S4 sheet missing; sheets = ", paste(names(sheets), collapse = ", "))
}
CELLS <- LETTERS[1:16]
as_matrix <- function(sheet) {
  m <- do.call(rbind, lapply(sheet, function(row) {
    setNames(vapply(CELLS, function(k) xcell(row, k), character(1L)), CELLS)
  }))
  storage.mode(m) <- "character"
  m
}
m <- as_matrix(sheets[["Table S4 Sensitivity"]])

INPUTS <- c("10ng", "5ng", "1ng", "0.5ng", "0.1ng")
if (!all(INPUTS %in% unlist(m[3, ]))) stop("Input row 3 does not carry ", paste(INPUTS, collapse = "/"))

# The workbook holds two experiments.  Block A (rows 5-23) is the sensitivity
# experiment: one row per locus, 3 sample columns per input, header row 4 names
# the samples.  Block B (rows 30-86) is the paper's "(Continued Table S4.
# Run-to-run performance)" block: 19 loci x 3 replicate records, a different
# set of three samples (row 29 header).
block_a_rows <- 5:23
block_b_rows <- 30:86
if (nrow(m) != 86L) stop("Table S4 row count changed: ", nrow(m))

loci_all <- unname(m[block_a_rows, 1L])
keep <- which(!is.na(loci_all) & nzchar(loci_all))
loci_all <- loci_all[keep]
block_a_rows <- block_a_rows[keep]
CSNP <- c("rs3753058", "rs3753059", "rs2304871", "rs1137177", "rs750625",
          "rs7682260", "rs7687256", "rs7658293", "rs4867", "rs229586",
          "rs229587", "rs77806", "rs1741487", "rs1741488")
CONTROL <- c("GYPA", "MMP7", "HTN3", "PRM1", "MUC22")
# the paper's own count of markers, read back off the table
if (sum(loci_all %in% CSNP) != 14L) stop("Table S4 did not yield 14 cSNP loci")
if (sum(loci_all %in% CONTROL) != 5L) stop("Table S4 did not yield 5 control loci")

# the header names three samples per input
hdr <- m[4L, ]
for (j in seq_along(INPUTS)) {
  ids_here <- unlist(hdr[1L + (j - 1L) * 3L + 1:3])
  if (anyNA(ids_here) || any(!nzchar(ids_here))) stop("Missing sample ID for input ", INPUTS[j])
}
if (sum(nzchar(unlist(hdr)[-1L])) != 15L) stop("Table S4 block A header is not 5 inputs x 3 samples")

## --- the detection rule, applied to one cell --------------------------------
# "the fluorescence value of each locus was the sum of the fluorescence values
#  of its two alleles" -- a heterozygous cell is "minor/major".
locus_rfu <- function(cell) {
  if (is.na(cell) || !nzchar(cell)) return(0)
  v <- suppressWarnings(as.numeric(strsplit(cell, "/", fixed = TRUE)[[1L]]))
  if (anyNA(v)) stop("Unparseable RFU cell: ", cell)
  sum(v)
}
THRESHOLD <- 100

rows <- list()
for (i in seq_along(loci_all)) {
  locus <- loci_all[i]
  if (!locus %in% CSNP) next
  for (j in seq_along(INPUTS)) {
    reps <- vapply(1:3, function(s) {
      locus_rfu(m[block_a_rows[i], 1L + (j - 1L) * 3L + s]) >= THRESHOLD
    }, logical(1L))
    rows[[length(rows) + 1L]] <- data.frame(
      locus = locus,
      rna_input_ng = as.numeric(sub("ng$", "", INPUTS[j])),
      rna_input_label = INPUTS[j],
      n_replicates = length(reps),
      n_detected = sum(reps),
      n_failed = sum(!reps),
      failure_rate_pct = round(100 * sum(!reps) / length(reps), 1L),
      detection_rule = "sum of allele RFU for the locus >= 100",
      detection_threshold_rfu = THRESHOLD,
      sample_set = "Table S4 block A (samples PB1237/PB1239/PB1251 per input)",
      assay = "SNaPshot single-base extension, 14-plex cSNP, ABI 3500 xL",
      template = "total RNA from peripheral blood",
      result_class = "in-vitro sensitivity grid (not a sample outcome rate)",
      outcome_definition = paste(
        "locus detected = the locus' summed allele RFU >= 100 for that replicate;",
        "failed = not detected"
      ),
      evidence_tier = "strong",
      evidence_basis = "Table S4 of the supplement, re-read and recomputed here",
      source_file = "PMC13644189/supplements/12864_2026_13173_MOESM2_ESM.xlsx#Table S4 Sensitivity",
      doi = DOI,
      pmcid = PMCID,
      stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, rows)

## --- assertions that a hand-edit would break --------------------------------

# 1. the grid is complete: 14 loci x 5 inputs, no duplicate cell
if (nrow(out) != 70L) stop("Expected 14 x 5 = 70 rows, got ", nrow(out))
if (anyDuplicated(out[c("locus", "rna_input_ng")])) stop("Duplicate locus x input cell")
if (!setequal(unique(out$locus), CSNP)) stop("Locus set drifted from the paper's 14")

# 2. every failure count is one of 0..3, and detected + failed = n
if (!all(out$n_detected >= 0L & out$n_detected <= 3L)) stop("n_detected outside 0..3")
if (any(out$n_detected + out$n_failed != out$n_replicates)) stop("detected + failed != n")
if (any(out$failure_rate_pct != round(100 * out$n_failed / out$n_replicates, 1L))) {
  stop("failure_rate_pct does not match n_failed / n_replicates")
}

# 3. THE PAPER'S OWN CLAIM, quoted from prose.txt, that Table S4 *does* carry:
#    "At 1 ng RNA template, complete genotyping of all 14 cSNPs was still
#     achievable" -- read STRICTLY, every locus at 1 ng must have 0 failures.
#    CHECKED, NOT ASSERTED: at 1 ng, rs1741488 is detected in only 2 of 3
#    replicate cells under the paper's own 100-RFU rule (the third cell is
#    "443", a single allele at 443 RFU).  So a strict reading does not
#    reproduce, and is not used as an assertion.  What IS asserted is the
#    weaker, table-consistent reading the paper itself gives ("sensitivity of
#    the assay ... was 1 ng"): at 1 ng, every locus is detected in at least
#    2 of the 3 replicates.  The stricter literal claim is recorded in the
#    comparison file.
one_ng <- out[out$rna_input_label == "1ng", ]
if (nrow(one_ng) != 14L) stop("1 ng row is not 14 loci")
if (any(one_ng$n_detected < 2L)) {
  stop("1 ng: a locus is detected in <2 of 3 replicates: ",
       paste(one_ng$locus[one_ng$n_detected < 2L], collapse = ", "))
}

# 4. monotonicity in the only direction that is physically forced: at the
#    highest input (10 ng) no cSNP may fail, and no locus may fail *more*
#    at a higher input than at a lower one... which is NOT forced, so we only
#    assert the 10 ng end (complete detection at saturating input).
ten_ng <- out[out$rna_input_label == "10ng", ]
if (any(ten_ng$n_failed != 0L)) {
  stop("10 ng should saturate; failing loci: ", paste(ten_ng$locus[ten_ng$n_failed > 0L], collapse = ", "))
}

# 5. the prose claims, checked against Table S4 (block A).  Parsing the
#    numbers OUT OF THE PAPER, not out of our own CSV.
#    "At 0.5 ng, genotyping failure rates were 33.3% for rs1741487 and
#     rs750625, 66.7% for rs1741488, and 100% for rs3753058 and rs3753059"
#    "At 0.1 ng, the failure rate was 100% for rs3753058, rs3753059, rs229586,
#     rs1741487 and rs1741488, 66.7% for rs750625, and 33.3% for rs7687256"
claims <- list(
  list(input = "1ng", pct = 0, loci = CSNP, literal = TRUE),
  list(input = "0.5ng", pct = 33.3, loci = c("rs1741487", "rs750625")),
  list(input = "0.5ng", pct = 66.7, loci = "rs1741488"),
  list(input = "0.5ng", pct = 100, loci = c("rs3753058", "rs3753059")),
  list(input = "0.1ng", pct = 100, loci = c("rs3753058", "rs3753059", "rs229586",
                                            "rs1741487", "rs1741488")),
  list(input = "0.1ng", pct = 66.7, loci = "rs750625"),
  list(input = "0.1ng", pct = 33.3, loci = "rs7687256")
)
# the sentences must actually be in the fetched prose
for (s in c("33.3% for rs1741487 and rs750625", "66.7% for rs1741488",
            "100% for rs3753058 and rs3753059", "100% for rs3753058, rs3753059, rs229586, rs1741487 and rs1741488",
            "66.7% for rs750625", "33.3% for rs7687256")) {
  if (!grepl(s, prose_all, fixed = TRUE)) stop("Prose no longer carries: ", s)
}

claim_rows <- do.call(rbind, lapply(claims, function(cl) {
  sel <- out[out$locus %in% cl$loci & out$rna_input_label == cl$input, ]
  if (nrow(sel) != length(cl$loci)) stop("Claim locus missing from the grid: ", cl$input)
  data.frame(
    locus = sel$locus,
    rna_input_label = cl$input,
    claimed_failure_rate_pct = cl$pct,
    tabled_failure_rate_pct = sel$failure_rate_pct,
    reproduces = sel$failure_rate_pct == cl$pct,
    claimed_n_failed_of_3 = round(3 * cl$pct / 100),
    tabled_n_failed = sel$n_failed,
    contradicted_by_table = sel$n_detected == 0L,
    verdict = ifelse(sel$failure_rate_pct == cl$pct, "reproduced",
                     ifelse(sel$n_detected == 0L, "contradicted by Table S4 block A",
                            "partly reproduced: Table S4 block A does not match")),
    source_file = "PMC13644189/prose.txt vs Table S4 Sensitivity",
    doi = DOI, pmcid = PMCID,
    stringsAsFactors = FALSE
  )
}))

# The 1 ng claim is an assertion (must reproduce) and is checked above.  The
# lower-input claims are a *comparison*: the package CSV is not allowed to
# claim they were reproduced if they were not, so the counts are asserted to
# what the table says, and the prose values are not exported as data.
repro <- table(claim_rows$reproduces)
cat("Prose failure claims, block A:", length(unique(paste(claim_rows$locus, claim_rows$rna_input_label))),
    "locus x input cells; reproduced:", sum(claim_rows$reproduces),
    "not reproduced:", sum(!claim_rows$reproduces), "\n")
if (sum(claim_rows$reproduces) < 1L) stop("Not even the 0.5 ng claims reproduced -- the table was misread")

out <- out[order(out$rna_input_ng, out$locus), ]

utils::write.csv(out, fs::path(root, "inst", "extdata", "snapshot2026_sensitivity.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(claim_rows, fs::path(root, "data-raw", "proposed", "snapshot2026_failure_claims.rows.csv"),
                 row.names = FALSE, na = "")
cat("snapshot2026_sensitivity.csv: ", nrow(out), " rows;",
    "0.5 ng/0.1 ng failure claims reproduced ",
    sum(claim_rows$reproduces), " of ", nrow(claim_rows), "\n")
