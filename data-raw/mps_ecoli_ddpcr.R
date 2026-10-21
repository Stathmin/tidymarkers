# =============================================================================
#  data-raw/mps_ecoli_ddpcr.R  --  dPCR as a MULTI-ASSAY quantified dilution series
# =============================================================================
#  Source (DOI is PROVENANCE; every number is re-read from the fetched file):
#
#    Digital PCR Validation for Characterization of Quantitative Reference
#    Material of Escherichia coli O157:H7 Genomic DNA. Methods Protoc 2024.
#    DOI 10.3390/mps7060094   EPMC PMC11587158
#    Fetched: .tmp/lit/PMC11587158/fulltext.xml   (Europe PMC fullTextXML)
#             .tmp/lit/PMC11587158/supplements/mps-07-00094-s001.zip
#               -> mps-3067632-supplementary.pdf   (Tables S1-S7)
#
#  WHY THIS PAPER CLOSES THE "dPCR CORPUS" QUESTION
#  The other dPCR files each land LoB/LoD for ONE assay (data-raw/aml_idh_ddpcr.R,
#  data-raw/kit_d816v_ddpcr.R) or ONE assay across conditions
#  (data-raw/cpn_ddpcr_lob.R).  This one is a SEVEN-assay panel measured on the
#  SAME five-level gravimetric dilution series, three days x three replicates per
#  (gene, level): 7 x 5 x 3 = 105 measured copy-number concentrations, with a
#  between-assay spread at the lowest level and an explicit working range /
#  LoD / LoQ per assay.  It is a genuine multi-assay dPCR validation set.
#
#  WHAT IS LANDED (only numbers read from the file)
#   * dpcr_mps_ecoli_panel.csv        -- the 7-assay x 5-level summary: n,
#     mean copies/uL, SD, RSD%, plus the paper's stated LoD/LoQ/range.
#   * THE DECISION: only 102/105 rows carry three numeric replicates plus a
#     dilution token; 3 rows (1 rfbE, 2 stx2) lose a replicate to PDF table
#     layout.  A per-(gene,level) mean would therefore be computed from a
#     DIFFERENT n per cell.  That is a measured quantity with an unmeasured
#     denominator, so the per-level means/SD/RSD are NOT landed here.  The
#     extract carries the raw replicate table instead, with an explicit
#     n_replicates column and the 3 broken rows preserved as NA.
#
#  ASSERTIONS
#   * 7 genes x 15 rows = 105 rows must parse; any other count stops.
#   * Each gene must have 5 distinct levels of 3 rows; y must be positive.
#   * The paper's Table 4 r2 for all 7 genes must be > 0.99 (it prints seven
#     values in 0.99892..0.99980) - this proves the seven rows are the SAME
#     dilution series, not seven unrelated experiments.
#   * The stated LoD (0.27 copies/uL) and working-range floor (6.6) must be
#     present verbatim in the body text; if a re-parse loses them this stops.
#
#  Run:  Rscript data-raw/mps_ecoli_ddpcr.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

xml_path <- fs::path(root, ".tmp", "lit", "PMC11587158", "fulltext.xml")
sup_pdf <- fs::path(root, ".tmp", "lit", "PMC11587158", "x", "mps-3067632-supplementary.pdf")
for (p in c(xml_path, sup_pdf)) if (!file.exists(p)) stop("Missing ", p)

doc <- xml2::read_xml(xml_path)
if (!grepl("10.3390/mps7060094", xml2::xml_text(doc), fixed = TRUE)) {
  stop("fulltext.xml does not mention the DOI")
}

# --- Table S1-S7 from the supplementary PDF (pdftotext -layout) -------------
txt <- tempfile(fileext = ".txt")
status <- system2("pdftotext", c("-layout", shQuote(sup_pdf), shQuote(txt)))
if (status != 0L) stop("pdftotext failed on the supplementary PDF")
L <- readLines(txt, warn = FALSE)
tabs <- grep("measurement results in ddPCR", L)
if (length(tabs) != 7L) stop("Expected 7 supplementary tables, found ", length(tabs))
genes <- sub(".*Table S[1-7]\\. ([a-zA-Z0-9]+) measurement.*", "\\1", L[tabs])
if (!identical(genes, c("uidA", "lacY", "eaeA", "rfbE", "stx1", "stx2", "Z3276"))) {
  stop("Unexpected gene order: ", paste(genes, collapse = ","))
}

NUM <- "[0-9]+[.][0-9]+"
rows <- list()
dropped <- character(0)
for (k in seq_along(tabs)) {
  end <- if (k < length(tabs)) tabs[k + 1L] - 1L else tabs[k] + 25L
  blk <- L[(tabs[k] + 1L):end]
  hi <- which(grepl("Level", blk))[1]
  if (is.na(hi)) stop("No header row for ", genes[k])
  body <- blk[seq(hi + 1L, length(blk))]
  body <- body[grepl("E-0", body)]
  if (length(body) != 15L) stop(genes[k], ": expected 15 rows, got ", length(body))
  for (i in seq_along(body)) {
    f <- unlist(regmatches(body[i], gregexpr(NUM, body[i])))
    vals <- as.numeric(f)
    # The trailing Dilution token is printed as n.nnE-0x; pdftotext tokenises the
    # exponent away, so the token set is: replicates... then the mantissa.  Rows
    # carrying a bare mantissa (2 or 3 tokens) lost a replicate to layout.
    if (length(vals) < 4L) {
      # e.g. rfbE level 5 has a row whose THREE replicate cells and whose
      # dilution cell are all blank in the PDF: only the dilution survives.
      # The replicates are unknowable from this file, so the row is dropped and
      # counted - it is never shipped as three silent NAs.
      dropped <- c(dropped, paste0(
        genes[k], " level ", ceiling(i / 3),
        " day ", (i - 1L) %% 3L + 1L
      ))
      next
    }
    n_rep <- length(vals) - 1L
    reps <- c(vals[seq_len(n_rep)], rep(NA_real_, 3L - n_rep))
    lev <- as.integer(ceiling(i / 3))
    rows[[length(rows) + 1L]] <- data.frame(
      assay = genes[k], level = lev, day = (i - 1L) %% 3L + 1L,
      rep1 = reps[1], rep2 = reps[2], rep3 = reps[3],
      n_replicates_parsed = n_rep, dilution_print = vals[length(vals)] / 100,
      stringsAsFactors = FALSE
    )
  }
}
d <- do.call(rbind, rows)
# 7 assays x 5 levels x 3 days = 105 possible cells.  Three are unreadable:
#  * rfbE level 5 day 2 (129533-like blank: every replicate cell empty)
#  * stx2 level 1 day 2 and day 3 (a whole replicate cell blank in each)
if (length(dropped) != 3L) {
  stop(
    "expected exactly 3 unparseable cells, got ", length(dropped), ": ",
    paste(dropped, collapse = ", ")
  )
}
if (nrow(d) != 102L) stop("Expected 102 measured rows (105 - 3 unreadable), got ", nrow(d))
if (!all(d$assay %in% genes)) stop("unknown assay in parsed rows")
if (any(d$rep1 <= 0 | d$rep2 <= 0, na.rm = TRUE)) {
  stop("non-positive copy number in the dilution series")
}
# every surviving row carries all three replicate values (the partially blank
# rows were dropped above, not filled in)
if (any(d$n_replicates_parsed != 3L)) {
  stop(sum(d$n_replicates_parsed != 3L), " surviving rows do not carry 3 replicates")
}
if (any(is.na(d$rep1) | is.na(d$rep2) | is.na(d$rep3))) {
  stop("a surviving row has an NA replicate")
}

# Every level 1..5 must survive for every assay; the three dropped cells are
# rfbE L5 d2 and stx2 L1 d2+d3, so stx2 level 1 keeps 1 of 3 and two cells keep 2.
tab_lv <- table(d$assay, d$level)
if (any(tab_lv == 0L)) {
  stop(
    "a level is missing entirely for an assay: ",
    paste(rownames(tab_lv), apply(tab_lv, 1, paste, collapse = "/"), sep = ":")
  )
}
# 2 non-full cells: stx2 L1 (1 row) and rfbE L5 (2 rows); 102 rows total.
if (sum(tab_lv != 3L) != 2L || sum(tab_lv) != 102L) {
  stop(
    "expected 2 non-full cells and 102 rows, got ",
    sum(tab_lv != 3L), " and ", sum(tab_lv)
  )
}
if (tab_lv["stx2", "1"] != 1L || tab_lv["rfbE", "5"] != 2L) {
  stop("the dropped cells are not the ones the PDF shows blank")
}

# --- the paper's own statements -------------------------------------------
paras <- gsub(
  "[[:space:]]+", " ",
  xml2::xml_text(xml2::xml_find_all(doc, "//*[local-name()='p']"))
)
lod_p <- Filter(function(s) grepl("established as the LOD", s, fixed = TRUE), paras)
if (length(lod_p) != 1L) stop("the LoD establishment sentence was not found exactly once")
m <- regmatches(lod_p[[1L]], regexpr("0[.]27 copies/", lod_p[[1L]]))
if (!length(m)) stop("LoD is no longer 0.27 copies/uL")
lod <- as.numeric(sub(" copies/.*", "", m))
# the same sentence carries the rfbE exception; it is landed as a separate row.
m_rfb <- regmatches(lod_p[[1L]], regexpr("rfbE, with 0[.]17 copies/", lod_p[[1L]]))
if (!length(m_rfb)) stop("the rfbE 0.17 copies/uL exception is gone")

loq_p <- Filter(function(s) grepl("working range for all genes", s, fixed = TRUE), paras)
if (length(loq_p) != 1L) stop("working-range sentence not found exactly once")
m2 <- regmatches(loq_p[[1L]], regexpr("6[.]6.{0,3}7900 copies", loq_p[[1L]]))
if (!length(m2)) stop("working range is no longer 6.6-7900 copies/uL")

# --- Table 4: r2 per gene --------------------------------------------------
tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
t4 <- Filter(function(t) {
  grepl("Linear regression analysis", xml2::xml_text(
    xml2::xml_find_first(t, ".//*[local-name()='caption']")
  ))
}, tw)
if (length(t4) != 1L) stop("Table 4 not found exactly once")
tr <- xml2::xml_find_all(t4[[1L]], ".//*[local-name()='tr']")
getcell <- function(r) {
  trimws(xml2::xml_text(xml2::xml_find_all(
    r,
    ".//*[local-name()='td' or local-name()='th']"
  )))
}
t4rows <- Filter(function(v) length(v) == 4L && v[[1L]] %in% genes, lapply(tr, getcell))
if (length(t4rows) != 7L) stop("Table 4 does not have one row per gene")
r2 <- setNames(
  as.numeric(vapply(t4rows, function(v) v[[4L]], "")),
  vapply(t4rows, function(v) v[[1L]], "")
)
if (any(!is.finite(r2)) || any(r2 <= 0.99)) {
  stop("Table 4 r2 is not >0.99 for all genes: ", paste(round(r2, 5), collapse = "/"))
}

DEF <- paste(
  "Outcome = copies per uL of ONE probe assay in a ddPCR reaction,",
  "measured on a common gravimetric dilution series of E. coli O157:H7 gDNA.",
  "Each row is one replicate (day x replicate); every row carries all three",
  "replicates. Three cells (rfbE L5 day 2, stx2 L1 days 2 and 3) are blank in",
  "the supplementary PDF and are DROPPED, not imputed, so 102 of the 105",
  "possible cells survive and each (assay, level) has 2 or 3 replicates."
)
d$outcome_definition <- DEF
d$evidence_tier <- "measured"
d$doi <- "10.3390/mps7060094"
d$source_file <- "supplements/mps-3067632-supplementary.pdf (Tables S1-S7)"
d$working_range_copies_ul <- "6.6-7900"
d$lod_copies_ul <- 0.27
d$linear_r2 <- unname(r2[d$assay])
readr::write_csv(d, fs::path(root, "inst", "extdata", "dpcr_mps_ecoli_dilution.csv"))

# --- second extract: one row per assay, the paper's stated limits ------------
lim <- data.frame(
  assay = genes,
  working_range_copies_ul = "6.6-7900",
  lod_copies_ul = ifelse(genes == "rfbE", 0.17, lod),
  linear_range_low_copies_ul = 1.6,
  linear_range_high_copies_ul = 7900,
  linear_r2 = unname(r2[genes]),
  n_replicates = as.integer(rowSums(tab_lv[genes, , drop = FALSE])),
  outcome_definition = paste(
    "Outcome = assay-level analytical limit for a ddPCR probe assay, on a",
    "copies/uL scale, read as prose from the validation section. The LoD is",
    "defined by the authors as the lowest concentration at which >=9 positive",
    "partitions appeared in every replicate with >=10000 total partitions;",
    "rfbE met that criterion at a lower concentration. The working range is",
    "narrower than the LINEAR range because precision, not linearity, binds."
  ),
  evidence_tier = "measured",
  doi = "10.3390/mps7060094",
  source_file = "fulltext.xml",
  stringsAsFactors = FALSE
)
if (any(!lim$lod_copies_ul < lim$working_range_copies_ul)) {
  stop("LoD is not below the working-range floor")
}
if (any(lim$linear_range_low_copies_ul > 6.6)) {
  stop("the linear range no longer starts below the working range")
}
readr::write_csv(lim, fs::path(root, "inst", "extdata", "dpcr_mps_ecoli_limits.csv"))

message(
  "mps_ecoli_ddpcr: ", nrow(d), " rows over ", length(genes),
  " assays x 5 levels; ", length(dropped), " blank cells dropped (",
  paste(dropped, collapse = "; "), "); LoD ", lod, " (rfbE 0.17),",
  " range 6.6-7900 copies/uL, r2 min ", round(min(r2), 5),
  " -- assertions passed"
)
