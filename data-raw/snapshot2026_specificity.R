# =============================================================================
#  data-raw/snapshot2026_specificity.R  --  SNaPshot detection by body fluid
# =============================================================================
#  Source (identifier = PROVENANCE; numbers come from the file):
#
#    "High-efficiency 14-plex cSNP SNaPshot assay for forensic peripheral blood
#    identification", BMC Genomics 2026;27.  DOI 10.1186/s12864-026-13173-0
#    (PMID 42436368, PMCID PMC13644189).
#
#  Fetched: .tmp/lit/PMC13644189/supplements/12864_2026_13173_MOESM2_ESM.xlsx
#  sheet "Table S8 500 sample analysis":
#    - 500 sample rows, 100 per body fluid (PB, MB, SE, SA, VA);
#    - 5 positive-control loci (GYPA, MMP7, HTN3, PRM1, MUC22) in columns C-G;
#    - the 14 cSNP loci in columns I,J,L,M,N,P,Q,R,S,U,V,W,X,Y
#      (the intervening columns are spacers, carrying no value in any row).
#
#  WHAT THIS LANDS  (inst/extdata/snapshot2026_specificity.csv)
#  -------------------------------------------------------------
#  One row per (locus x body fluid): how many of the 100 samples in that fluid
#  had the locus detected, using the paper's own rule (Methods, quoted in
#  prose.txt): "a detection threshold of 100 relative fluorescence units ...
#  Peaks below 100 RFU were excluded for evaluating allelic dropout".
#
#  => detected := sum of the allele RFU in the cell >= 100.
#
#  THE ASSERTIONS THAT MAKE THIS THE PAPER'S TABLE AND NOT OURS
#  ------------------------------------------------------------
#  * 100 + 100 + 100 + 100 + 100 = 500, and every row's threshold sum is <= 100.
#  * PB: 14/14 loci detected in 100/100 samples (the paper's headline: PB
#    separates from every other fluid).
#  * Each positive control is detected in 100/100 samples of exactly the fluid
#    it marks: MMP7 -> MB (but also 7/100 VA), HTN3 -> SA, PRM1 -> SE,
#    MUC22 -> VA.  GYPA is NOT that kind of control: it is the blood marker
#    shared by PB and MB, and the table has it at 84/100 in PB and 67/100 in
#    MB.  Its PB value is landed, not asserted to 100.
#  * The paper says (p46, quoted in prose.txt) "the positive control locus
#    GYPA had a detection rate of 76% in menstrual blood samples".  The table
#    gives 67%.  The own-file number is what is landed; the paper's 76% is
#    recorded in data-raw/proposed/ as a mismatch, NOT overwritten.
#
#  Run:  Rscript data-raw/snapshot2026_specificity.R
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

doc <- xml2::read_xml(xml_path)
title <- xml2::xml_text(xml2::xml_find_first(doc, "//*[local-name()='article-title']"))
if (!grepl("SNaPshot", title, fixed = TRUE)) stop("Not the SNaPshot paper: ", title)
ids <- xml2::xml_text(xml2::xml_find_all(doc, "//*[local-name()='article-id']"))
if (!DOI %in% ids) stop("XML does not carry DOI ", DOI)

prose_all <- paste(readLines(prose_path, warn = FALSE), collapse = " ")
for (phrase in c("detection threshold of 100 relative fluorescence units",
                 "detection rate of 76% in menstrual blood")) {
  if (!grepl(phrase, prose_all, fixed = TRUE)) stop("Fetched prose lost: ", phrase)
}

source(fs::path(root, "data-raw", "lib", "xlsx-read.R"))
sheets <- xlsx_sheets(xlsx_path)
SHEET <- "Table S8 500 sample analysis "
if (!SHEET %in% names(sheets)) {
  stop("Table S8 sheet missing; sheets = ", paste(names(sheets), collapse = ", "))
}
CELLS <- LETTERS[1:25]
m <- do.call(rbind, lapply(sheets[[SHEET]], function(row) {
  setNames(vapply(CELLS, function(k) xcell(row, k), character(1L)), CELLS)
}))
storage.mode(m) <- "character"
if (nrow(m) != 503L) stop("Table S8 row count changed: ", nrow(m))

FLUIDS <- c("Peripheral Blood (PB)", "Menstrual blood (MB)", "Seminal fluid (SE)",
            "Saliva (SA)", "Vaginal secretions (VA)")

# header: row 3 carries the locus names; row 2 names the genes over each block
hdr <- m[3L, ]
CSNP <- c("rs3753058", "rs3753059", "rs2304871", "rs1137177", "rs750625",
          "rs7682260", "rs7687256", "rs7658293", "rs4867", "rs229586",
          "rs229587", "rs77806", "rs1741487", "rs1741488")
CTRL_NAMES <- c("GYPA", "MMP7", "HTN3", "PRM1", "MUC22")
# locate every non-spacer column by the name printed in the header.
# Column B is the "Sample ID" label, not a locus, so start at C.
loc_cols <- list()
for (col in LETTERS[3:25]) {
  nm <- hdr[[col]]
  if (is.na(nm) || !nzchar(nm)) next
  loc_cols[[nm]] <- col
}
if (!all(CSNP %in% names(loc_cols))) {
  stop("cSNP loci missing from row 3 header: ",
       paste(setdiff(CSNP, names(loc_cols)), collapse = ", "))
}
if (!all(CTRL_NAMES %in% names(loc_cols))) stop("Control loci missing from header")
if (length(loc_cols) != 19L) stop("Header did not give 19 named loci, got ", length(loc_cols))

# body: rows 4..503; column A carries the fluid on the first row of each
# 100-sample run only (rows 1, 101, 201, 301, 401 of the body), so the label
# is carried forward.  The label spans two lines in-cell ("Peripheral Blood"
# / "(PB)"), so whitespace is collapsed.
body <- m[4:503, , drop = FALSE]
if (nrow(body) != 500L) stop("Table S8 body is not 500 rows: ", nrow(body))
labels_at <- which(!is.na(body[, "A"]) & nzchar(trimws(body[, "A"])))
if (!identical(labels_at, c(1L, 101L, 201L, 301L, 401L))) {
  stop("Body-fluid label rows are not 1/101/201/301/401: ", paste(labels_at, collapse = ","))
}
label_txt <- trimws(gsub("[[:space:]]+", " ", body[labels_at, "A"]))
if (!all(label_txt %in% FLUIDS)) {
  stop("Unrecognised body-fluid label: ", paste(setdiff(label_txt, FLUIDS), collapse = ", "))
}
fluid <- rep(label_txt, each = 100L)
if (any(!nzchar(body[, "B"]))) stop("A sample row has no sample ID")

per_fluid <- table(factor(fluid, levels = FLUIDS))
if (length(per_fluid) != 5L) stop("Expected 5 body fluids")
if (!all(per_fluid == 100L)) stop("Expected 100 samples per fluid, got ", paste(per_fluid, collapse = "/"))
if (sum(per_fluid) != 500L) stop("Expected 500 samples in total, got ", sum(per_fluid))

locus_sum <- function(cell) {
  if (is.na(cell) || !nzchar(cell)) return(0)
  v <- suppressWarnings(as.numeric(strsplit(cell, "/", fixed = TRUE)[[1L]]))
  if (anyNA(v)) stop("Unparseable RFU cell: ", cell)
  sum(v)
}
THRESHOLD <- 100

loci <- lapply(names(loc_cols), function(nm) {
  vals <- vapply(body[, loc_cols[[nm]]], locus_sum, numeric(1L))
  list(name = nm, detected = vals >= THRESHOLD)
})
names(loci) <- names(loc_cols)

rows <- list()
for (nm in names(loci)) {
  for (f in FLUIDS) {
    sel <- fluid == f
    n_det <- sum(loci[[nm]]$detected[sel])
    is_ctrl <- nm %in% CTRL_NAMES
    rows[[length(rows) + 1L]] <- data.frame(
      locus = nm,
      body_fluid = sub(" .*$", "", f),
      body_fluid_full = f,
      marker_class = if (is_ctrl) "positive control" else "cSNP",
      n_samples = sum(sel),
      n_detected = n_det,
      detection_rate_pct = round(100 * n_det / sum(sel), 1L),
      detection_rule = "sum of allele RFU for the locus >= 100",
      detection_threshold_rfu = THRESHOLD,
      sample_set = "Table S8, 100 single-source samples per body fluid",
      assay = "SNaPshot single-base extension, 14-plex cSNP, ABI 3500 xL",
      template = "RNA from a single-source body fluid",
      result_class = "specificity matrix (detection, not genotyping accuracy)",
      outcome_definition = paste(
        "detected = the locus' summed allele RFU >= 100 for that sample;",
        "rate is the share of the fluid's 100 samples where it was detected"
      ),
      evidence_tier = "strong",
      evidence_basis = "Table S8 of the supplement, re-read and recomputed here",
      source_file = "PMC13644189/supplements/12864_2026_13173_MOESM2_ESM.xlsx#Table S8 500 sample analysis",
      doi = DOI, pmcid = PMCID,
      stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, rows)

## --- assertions -------------------------------------------------------------

if (nrow(out) != 95L) stop("Expected 19 loci x 5 fluids = 95 rows, got ", nrow(out))
if (anyDuplicated(out[c("locus", "body_fluid_full")])) stop("Duplicate locus x fluid cell")
if (!all(out$n_samples == 100L)) stop("A fluid does not have exactly 100 samples")
if (any(out$n_detected < 0L | out$n_detected > out$n_samples)) stop("n_detected out of range")
if (any(out$detection_rate_pct > 100 | out$detection_rate_pct < 0)) stop("rate out of range")

# 1. the paper's headline: every cSNP is detected in every PB sample
pb_csnp <- out[out$body_fluid == "Peripheral" & out$marker_class == "cSNP", ]
if (nrow(pb_csnp) != 14L) stop("PB cSNP slice is not 14 loci")
if (any(pb_csnp$n_detected != 100L)) {
  stop("PB should detect all 14 cSNPs in 100/100; failing: ",
       paste(pb_csnp$locus[pb_csnp$n_detected != 100L], collapse = ", "))
}

# 2. specificity: outside PB and MB, no cSNP is detected in even half a fluid.
#    (PB and MB both carry the blood-specific markers - the paper says so.)
not_blood <- out[out$marker_class == "cSNP" & !out$body_fluid %in% c("Peripheral", "Menstrual"), ]
if (any(not_blood$detection_rate_pct >= 50)) {
  stop("A cSNP is detected in >=50% of a non-blood fluid: ",
       paste(not_blood$locus[not_blood$detection_rate_pct >= 50], collapse = ", "))
}

# 3. each positive control is fluid-specific: 100/100 in exactly the fluid it
#    marks (seminal/saliva/vaginal), and absent (0) from those it does not.
#    MMP7 also fires in 7% of vaginal samples; those are the table's numbers
#    and are landed as-is, so the assertion is on the marker's OWN fluid
#    (100), not on its absence elsewhere.  GYPA is the PB/MB marker, not a
#    one-fluid control: its own rate is landed, and only its PB dominance is
#    asserted.
own_fluid <- c(MMP7 = "Menstrual", HTN3 = "Saliva", PRM1 = "Seminal", MUC22 = "Vaginal")
for (mk in names(own_fluid)) {
  r <- out[out$locus == mk & out$body_fluid == own_fluid[[mk]], ]
  if (nrow(r) != 1L) stop("Control ", mk, " missing for ", own_fluid[[mk]])
  if (r$n_detected != 100L) stop("Control ", mk, " is not 100/100 in ", own_fluid[[mk]],
                                 " (got ", r$n_detected, ")")
}
gy <- out[out$locus == "GYPA", ]
if (nrow(gy) != 5L) stop("GYPA is not in all five fluids")
# GYPA is a blood marker: present in PB and MB, essentially absent elsewhere
if (gy$detection_rate_pct[gy$body_fluid == "Peripheral"] < 50) stop("GYPA below 50% in PB")
if (gy$detection_rate_pct[gy$body_fluid == "Menstrual"] < 50) stop("GYPA below 50% in MB")
for (f in c("Seminal", "Saliva", "Vaginal")) {
  if (gy$detection_rate_pct[gy$body_fluid == f] > 10) stop("GYPA leaks into ", f)
}

# 4. the paper's own GYPA sentence, compared not asserted
gy_mb <- out$detection_rate_pct[out$locus == "GYPA" & out$body_fluid == "Menstrual"]
claim <- data.frame(
  locus = "GYPA", body_fluid = "Menstrual",
  claimed_detection_rate_pct = 76,
  tabled_detection_rate_pct = gy_mb,
  reproduces = gy_mb == 76,
  verdict = if (gy_mb == 76) "reproduced" else "Table S8 gives a different rate",
  source_file = "PMC13644189/prose.txt vs Table S8 500 sample analysis",
  doi = DOI, pmcid = PMCID, stringsAsFactors = FALSE
)

out <- out[order(out$marker_class, out$locus, out$body_fluid), ]
utils::write.csv(out, fs::path(root, "inst", "extdata", "snapshot2026_specificity.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(claim, fs::path(root, "data-raw", "proposed", "snapshot2026_specificity_claims.rows.csv"),
                 row.names = FALSE, na = "")
cat("snapshot2026_specificity.csv:", nrow(out), "rows;",
    "PB cSNP 100/100 for 14 loci; GYPA in MB tabled", gy_mb, "% vs paper's 76%\n")
