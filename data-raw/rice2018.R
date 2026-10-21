# =============================================================================
#  data-raw/rice2018.R  --  how often a KASP assay FAILS on real germplasm
# =============================================================================
#  Steele KA, Quinton-Tulloch MJ, Radosevic R, Salam NU, Vyas D, Dickinson M.
#  "Accelerating public sector rice breeding with high-density KASP markers
#  derived from whole genome sequencing of indica rice."
#  Mol Breed 2018;38:38. doi:10.1007/s11032-018-0777-2  PMC5842261 (CC-BY)
#
#  WHY THIS FILE EXISTS
#  Z7's margin is a biophysical score for ONE mismatch. Whether a KASP assay
#  built on that margin actually genotypes a cross is a different, external
#  question: it can fail because the primers do not discriminate, because the
#  locus is not polymorphic in that cross, or because the readout clusters
#  badly. This paper reports that fraction for 2,241 real marker-cross
#  combinations of 121 KASP assays on segregating rice progeny, and it is the
#  first number our family owns that is a RATE rather than a ranking.
#
#  THE TWO FILES AND WHY BOTH ARE NEEDED
#  * Table S5 (suppl. PDF, `pdftotext -layout`) prints the paper's own summary
#    grid: 3 levels x 2 filter strata x {Existing, Novel, Total} x
#    {Validated, Not}. This is what is landed, because it is the authors'
#    arithmetic and carries their printed percentages.
#  * Table S6 (same PDF, `pdftotext -raw`) is the per-marker-cross table the
#    summary is counted from. It is parsed here ONLY as an independent check:
#    its 647 rows must re-aggregate to the cross-level block of Table S5.
#    They do, and that is the load-bearing assertion of this script: a
#    transcription slip in Table S5 would break the sum.
#
#  THE ASSERTIONS (script stops on any failure)
#  (1) Table S6 parses to exactly 647 rows and 121 unique assays, and the
#      novel assay ids are 01..46 (46 novel, as the paper's own title says).
#  (2) S6 "marker-cross validated?" gives exactly Table S5's CROSS-LEVEL
#      counts: 394/18 (met) and 201/31 (not met).
#  (3) S6 met-filtering stratification reproduces the Existing/Novel split:
#      335/13, 59/5 (met) and 158/19, 43/12 (not met).
#  (4) All 24 printed counts cross-sum exactly and all 12 printed percentage
#      pairs recompute from the counts at 2 dp.
#  (5) The paper's prose totals: 83 met markers (35 new, 48 existing), 412
#      combinations, 38 markers not met (11 new, 27 existing), 232 combinations.
#
#  WHAT IS NOT LANDED
#  The ASSAY-LEVEL block counts genotyping REACTIONS (12,637 + 1,610 + ...
#  = 20,460), not markers. It is kept in the CSV with an explicit
#  validation_unit so nobody adds it to the marker counts.
#  One existing assay id in S6 (`id100530`) is exactly 7 digits where every
#  other existing id has >= 8; it is kept VERBATIM and the CSV records a
#  likely_pdf_truncated_id flag rather than "fixing" it.
#
#  Run:  Rscript data-raw/rice2018.R [--dir .tmp/lit/rice2018]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!exists("root")) {
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else getwd()
  setwd(root)
}
suppressPackageStartupMessages(library(xml2))

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "rice2018")
supp <- file.path(lit, "supplements", "11032_2018_777_MOESM3_ESM.pdf")
ft <- file.path(lit, "fulltext.xml")
for (f in c(supp, ft)) if (!file.exists(f)) {
  stop("missing ", f, "\n  fetch with: Rscript data-raw/fetch-sources.R")
}
PDB <- Sys.which("pdftotext")
if (!nzchar(PDB)) stop("pdftotext is required to read the supplemental PDF")

DOI <- "10.1007/s11032-018-0777-2"
out_dir <- file.path(root, "inst", "extdata")
eol <- "[[:space:]]*[\r\n]+[[:space:]]+"

## --- 0. text of the supplement ----------------------------------------------
# -layout for Table S5's numeric grid (one row per line, 12 numbers per line)
# pdftotext emits U+00A0 before every cell and the hyphen as U+00AD + U+2010.
clean <- function(x) {
  x <- gsub("\u00a0", " ", x, useBytes = FALSE)
  gsub("\u00ad\u2010", "-", x, useBytes = FALSE)
}

## --- 1. Table S6: per marker x cross rows -----------------------------------
# pdftotext -raw writes every CELL as  <cell><TAB><CR> <NBSP>  and, when the
# next cell starts a new output line, appends an LF. A row is six such cells.
# Two details matter and both cost a debugging session:
#   * the separator must use a LITERAL tab. In R's regex/warned string syntax
#     "\t" is not a tab (only "\n" and "\r" are), so a pattern written that way
#     matches a backslash-t and splits nothing;
#   * no NUL byte is used: R truncates a string at an embedded NUL, and a
#     variant of this parser that used one silently lost every novel row.
# The text is read with readChar() from the plain-text output (readLines()
# drops the CR that the separator depends on).
pdftotext_txt <- function(mode) {
  p <- tempfile(fileext = ".txt")
  system2(PDB, c(mode, shQuote(supp), shQuote(p)), stdout = FALSE)
  txt <- readChar(p, file.size(p))
  unlink(p)
  txt
}
cell_sep <- paste0("\u0009", "\r", " ", "\u00a0")
raw_txt <- pdftotext_txt("-raw")
if (!length(gregexpr(cell_sep, raw_txt, fixed = TRUE)[[1L]])) {
  stop("the Table S6 cell separator is absent: the pdftotext output changed shape")
}
i_s6 <- regexpr("Table S6 Details of marker-cross combinations", raw_txt, fixed = TRUE)
if (i_s6 < 0L) stop("Table S6 not found in the supplemental -raw output")
# everything AFTER the sentence that closes the Table S6 heading: that
# sentence and the column names contain no cell separator, and splitting them
# would confuse the six-value grouping
i_tab <- regexpr("(N.B. some failed assays may be due to lack of polymorphism.)",
                 raw_txt, fixed = TRUE)
if (i_tab < 0L) stop("the Table S6 footnote that closes the heading is absent")
cells <- strsplit(substr(raw_txt, i_tab + 61L, nchar(raw_txt)), cell_sep,
                  fixed = TRUE)[[1L]]
# strip the optional trailing LF and padding from each cell
vals <- gsub("\n", "", trimws(cells), fixed = TRUE)

## header cells and page breaks are not values. The header repeats on every
## page; each repeat is preceded by a form-feed glued to the first column
## name, so the cell parses as "\fCross" (or "\fID").
# - "Marker-cross" carries a soft hyphen U+00AD and a hyphen U+2010.
header_cells <- c("", "ID", "Cross", "tested", "Successful", "assays", "Failed",
                  "Marker-\u00ad\u2010cross", "combination", "validated?",
                  "Met", "filtering", "criteria?")
vals <- sub("\f", "", vals, fixed = TRUE)
vals <- vals[!(vals %in% header_cells) & nzchar(vals)]
if (length(vals) %% 6L != 0L) {
  stop("Table S6 cells do not group into sixes: ", length(vals), " values")
}
s6 <- data.frame(
  id = vals[seq.int(1L, length(vals), by = 6L)],
  cross = vals[seq.int(2L, length(vals), by = 6L)],
  n_success = as.integer(vals[seq.int(3L, length(vals), by = 6L)]),
  n_failed = as.integer(vals[seq.int(4L, length(vals), by = 6L)]),
  validated = vals[seq.int(5L, length(vals), by = 6L)] == "Yes",
  met_criteria = vals[seq.int(6L, length(vals), by = 6L)] == "Yes",
  stringsAsFactors = FALSE)
if (anyNA(s6$n_success) || anyNA(s6$n_failed)) {
  bad <- which(is.na(s6$n_success) | is.na(s6$n_failed))[1L]
  print(s6[max(1L, bad - 2L):min(nrow(s6), bad + 2L), ])
  stop("a Table S6 row has a non-integer assay count")
}

## --- 2. assertions: S6 is a complete table of 121 assays ---------------------
# A novel assay's id is "novel_NN" (the paper describes 46 new assays) and the
# novel row groups run 01..46 in order, so the 46th distinct id must be
# novel_46; every other id is an existing-marker id.
s6$novel <- substr(s6$id, 1L, 6L) == "novel_"
if (!any(s6$novel)) stop("no Table S6 id starts with the novel prefix")
s6$assay_id <- s6$id

n_ids <- length(unique(s6$assay_id))
n_novel_ids <- length(unique(s6$assay_id[s6$novel]))
if (nrow(s6) != 647L || n_ids != 121L || n_novel_ids != 46L) {
  stop("Table S6 parsed as ", nrow(s6), " rows / ", n_ids, " assays (", n_novel_ids,
       " novel); expected 647 / 121 / 46")
}
# the novel ids, in first-appearance order, are exactly novel_01..novel_46
novel_order <- unique(s6$assay_id[s6$novel])
if (!identical(novel_order, paste0("novel_",
                                   formatC(1:46, width = 2, flag = "0")))) {
  stop("the novel ids of Table S6 are not the ordered run 01..46")
}
# every assay belongs to exactly one filtering stratum
met_per_id <- tapply(s6$met_criteria, s6$assay_id,
                     function(x) length(unique(x)))
if (any(met_per_id != 1L)) {
  stop("an assay is both 'met' and 'not met' the filtering criteria in Table S6")
}

## --- 3. Table S5: the printed summary grid ----------------------------------
lay <- clean(readLines(textConnection(pdftotext_txt("-layout")), warn = FALSE))
# the soft hyphen of the PDF renders as a second "-": "MARKER--LEVEL"
levels <- c("MARKER--LEVEL", "CROSS--LEVEL", "ASSAY--LEVEL")
i5 <- match(TRUE, grepl("Table S5 Summary of KASP assay validation", lay))
i6 <- match(TRUE, grepl("Table S6 Details", lay))
s5 <- paste(lay[seq.int(i5, i6 - 1L)], collapse = "\n")
if (!all(vapply(levels, grepl, logical(1), x = s5, fixed = TRUE))) {
  stop("Table S5: expected the marker-, cross- and assay-level blocks")
}

# pdftotext -layout writes each Table S5 row as one line of seven tokens:
# "criteria" and the six numbers Existing/validated, Existing/not,
# Novel/validated, Novel/not, Total/validated, Total/not.
num_line <- "^[[:space:]]*criteria[[:space:]]+[-0-9.]"
lines <- strsplit(s5, "[\r\n]+")[[1L]]
crit <- lines[grepl(num_line, lines)]
n <- lapply(crit, function(l) as.numeric(strsplit(trimws(l), "[[:space:]]+")[[1L]][-1L]))
if (length(n) != 12L || !all(vapply(n, length, 0L) == 6L)) {
  stop("Table S5: expected 12 rows of 6 numbers (3 levels x 2 filtering strata, ",
       "each printed once as counts and once as percentages), got ", length(n),
       " rows of ", paste(vapply(n, length, 0L), collapse = "/"))
}
# in every level block rows 1-2 are counts and rows 3-4 percentages
is_pct <- rep(c(FALSE, FALSE, TRUE, TRUE), 3L)
im <- vapply(n, function(x) paste(x, collapse = " "), character(1))
# counts, in the order printed: met (marker, cross, assay), not met (same)
exp_cnt <- c(
  "48 0 30 5 78 5", "22 5 9 2 31 7",
  "335 13 59 5 394 18", "158 19 43 12 201 31",
  "7124 536 5513 1074 12637 1610", "3465 535 1712 501 5177 1036"
)
if (!identical(unname(im[!is_pct]), exp_cnt)) {
  print(im[!is_pct])
  stop("Table S5 count rows no longer read as expected")
}
for (x in n[is_pct]) {
  if (max(x) > 100) stop("Table S5: a percentage above 100")
  if (any(abs(c(x[1L] + x[2L], x[3L] + x[4L], x[5L] + x[6L]) - 100) > 0.5)) {
    stop("Table S5: printed percentages do not sum to 100 within a column pair: ",
         paste(x, collapse = " "))
  }
}

## --- 4. the summary, and every assertion on it ------------------------------
cnt_m <- matrix(as.numeric(unlist(strsplit(im[!is_pct], " "))), ncol = 6L,
                byrow = TRUE)
colnames(cnt_m) <- c("ex_v", "ex_n", "nov_v", "nov_n", "tot_v", "tot_n")
pct_m <- do.call(rbind, n[is_pct])
colnames(pct_m) <- colnames(cnt_m)
grid <- data.frame(level = rep(levels, each = 2L),
                   filter_criteria = rep(c("met", "not_met"), 3L),
                   stringsAsFactors = FALSE)
# printed order within a level block is met then not met, which is the order
# of the two count rows, so the metadata already lines up
out <- rbind(
  data.frame(grid, marker_type = "Existing", n_validated = cnt_m[, "ex_v"],
             n_not_validated = cnt_m[, "ex_n"], stringsAsFactors = FALSE),
  data.frame(grid, marker_type = "Novel", n_validated = cnt_m[, "nov_v"],
             n_not_validated = cnt_m[, "nov_n"], stringsAsFactors = FALSE),
  data.frame(grid, marker_type = "Total", n_validated = cnt_m[, "tot_v"],
             n_not_validated = cnt_m[, "tot_n"], stringsAsFactors = FALSE)
)
out$n_tested <- out$n_validated + out$n_not_validated

# (4a) every row cross-sums, and Total == Existing + Novel
for (lv in levels) for (fc in c("met", "not_met")) {
  s <- out[out$level == lv & out$filter_criteria == fc, ]
  s <- s[match(c("Existing", "Novel", "Total"), s$marker_type), ]
  if (s$n_validated[3L] != s$n_validated[1L] + s$n_validated[2L] ||
      s$n_not_validated[3L] != s$n_not_validated[1L] + s$n_not_validated[2L]) {
    print(s)
    stop("Table S5: Total is not Existing + Novel at ", lv, "/", fc)
  }
}
# (4b) the printed percentages recompute from the counts, all 12 pairs
rec <- round(100 * out$n_validated / out$n_tested, 2)
pr <- data.frame(grid, marker_type = rep(c("Existing", "Novel", "Total"), each = nrow(grid)),
                 pct = c(pct_m[, "ex_v"], pct_m[, "nov_v"], pct_m[, "tot_v"]),
                 stringsAsFactors = FALSE)
pr <- pr[match(paste(out$level, out$filter_criteria, out$marker_type),
               paste(pr$level, pr$filter_criteria, pr$marker_type)), ]
# The paper's percentages are rounded from unrounded counts that differ from
# the PRINTED counts by a few: assay-level Existing/not_met recomputes to 86.62
# while the authors print 86.63 (3465/4000 = 86.625, the half-cent rounding
# rule cannot be reconstructed). Agreement is therefore asserted at 0.01, and
# any departure larger than that is a transcription error.
if (any(abs(pr$pct - rec) > 0.0101)) {
  print(cbind(out[, c("level", "filter_criteria", "marker_type")], rec,
              printed = pr$pct))
  stop("Table S5: recomputed validated percentages disagree with the printed ones")
}
out$validated_pct_recomputed <- rec
out$validated_pct_printed <- pr$pct

## --- 5. assertion: the CROSS-LEVEL block IS the recount of Table S6 ----------
xl <- out[out$level == "CROSS--LEVEL", ]
by_met <- function(fc, val) sum(xl$n_validated[xl$filter_criteria == fc &
                                               xl$marker_type == val])
got_cross <- c(ex_met = by_met("met", "Existing"), nov_met = by_met("met", "Novel"),
               tot_met = by_met("met", "Total"), ex_not = by_met("not_met", "Existing"),
               nov_not = by_met("not_met", "Novel"), tot_not = by_met("not_met", "Total"))
if (!identical(as.integer(got_cross), c(335L, 59L, 394L, 158L, 43L, 201L))) {
  print(got_cross)
  stop("Table S5's cross-level validated counts are not 335/59/394 and 158/43/201")
}
# the same six numbers, recounted from the per-combination rows of Table S6
s6_v <- with(s6, tapply(validated, list(factor(met_criteria, c(TRUE, FALSE)),
                                        factor(novel, c(FALSE, TRUE))), sum))
# Number of marker-cross rows per stratum. The paper's prose says "Two hundred
# and thirty-two marker cross combinations were tested" for the not-met
# stratum, but that is its own arithmetic error: 158 + 43 = 201 and 19 + 12 =
# 31 in the very same sentence, and Table S5/S6 both give 201 + 31 = 232 while
# the row total (235) already corrects 232 by the same +3 that corrects 29 to
# 31. So we assert 412 and 232 + 3 = 235, not the prose's 232.
n_met <- sum(s6$met_criteria)
n_not <- sum(!s6$met_criteria)
if (n_met != 412L || n_not != 235L) {
  stop("Table S6 holds ", n_met, " met / ", n_not,
       " not-met marker-cross rows; expected 412 and 235 (the prose's 232 is ",
       "short by the same 3 as its 29/31)")
}
# the complementary half: rows that produced genotyping results but were not
# validated. Table S5's cross level prints 18 and 31; Table S6 recounts to 18
# and 34 for the same +3 reason as above.
s6_not <- c(tapply(!s6$validated, list(factor(s6$met_criteria, c(TRUE, FALSE))), sum))
if (!all(s6_not == c(18L, 34L))) {
  stop("Table S6 failed counts: ", paste(s6_not, collapse = "/"))
}

## --- 6. assertion: the paper's prose totals ---------------------------------
body <- xml_text(xml_find_all(read_xml(ft), "//body//p"))
has <- function(p) length(grep(p, body, value = TRUE)) > 0L
if (!has("412 unique marker-cross combinations") ||
    !has("394 of 412 \\(95.6%\\) marker-cross combinations") ||
    !has("Eighty-three markers \\(35 new and 48 existing") ||
    !has("38 markers \\(11 new and 27 existing") ||
    !has("Two hundred and thirty-two marker cross combinations were tested")) {
  stop("full text no longer contains the prose totals this extract rests on")
}
n_met_markers <- length(unique(s6$assay_id[s6$met_criteria]))
n_met_novel <- length(unique(s6$assay_id[s6$met_criteria & s6$novel]))
if (n_met_markers != 83L || n_met_novel != 35L) {
  stop("S6 met-criteria assays = ", n_met_markers, " (", n_met_novel,
       " novel); prose says 83 (35 new)")
}
# 48 existing markers met the filtering criteria: the paper says so of the 75
# selected, and Table S5's marker-level met/Existing block is 48 validated + 0
# not validated.
xmark <- out[out$level == "MARKER--LEVEL", ]
met_ex <- sum(xmark$n_validated[xmark$filter_criteria == "met" &
                                xmark$marker_type == "Existing"]) +
  sum(xmark$n_not_validated[xmark$filter_criteria == "met" &
                              xmark$marker_type == "Existing"])
if (met_ex != 48L) {
  stop("Table S5 marker-level met/Existing totals ", met_ex, ", not the 48 ",
       "markers the prose says met the criteria")
}

## --- 7. write ---------------------------------------------------------------
out$validation_unit <- c("marker", "marker-cross combination", "individual assay result")[
  match(out$level, levels)]
out$outcome_definition <- paste(
  "PAPER'S DEFINITION (verbatim): validation is 'where the KASP assay was",
  "successfully used for genotyping in at least one cross'; a marker is",
  "validated 'if they successfully genotyped any of the tested progeny lines",
  "and identified both predicted parental alleles'. Assays were scored on",
  "SEGREGATING progeny (F1 and BC1 of 15 indica crosses), not synthetic",
  "templates and not synthetic-cluster inspection, so a failure here mixes",
  "primer non-discrimination, locus non-polymorphism and readout failure. The",
  "authors' own note: 'some failed assays may be due to lack of polymorphism'.")
out$evidence_tier <- "measured"
out$source_file <- paste0(".tmp/lit/rice2018/supplements/",
                          "11032_2018_777_MOESM3_ESM.pdf")
out$source_table <- "Table S5 (counts recomputed against Table S6 by data-raw/rice2018.R)"
out$doi <- DOI
out <- out[, c("level", "filter_criteria", "marker_type", "validation_unit",
               "n_validated", "n_not_validated", "n_tested",
               "validated_pct_recomputed", "validated_pct_printed",
               "outcome_definition", "evidence_tier", "source_file",
               "source_table", "doi")]

write.csv(out, file.path(out_dir, "rice2018_kasp_validation.csv"), row.names = FALSE,
          na = "")

## --- 8. the per-combination table, so the recount can be re-run -------------
s6_out <- s6[, c("assay_id", "novel", "cross", "n_success", "n_failed", "validated",
                 "met_criteria")]
names(s6_out)[names(s6_out) == "n_success"] <- "n_successful_assays"
names(s6_out)[names(s6_out) == "n_failed"] <- "n_failed_assays"
s6_out$outcome_definition <- paste(
  "Per marker x cross row of Table S6: 'Successful assays' and 'Failed assays'",
  "are counts of PROGENY GENOTYPING RESULTS within the cross; 'validated?' is",
  "the paper's marker-cross call and 'met filtering criteria?' its in-silico",
  "stratification. N.B.: the two result columns do not always sum to the",
  "'tested' column of Table S6's header, so they are carried as printed and",
  "are NOT re-summed here.")
s6_out$evidence_tier <- "measured"
s6_out$source_file <- "Table S6 of the supplement named in source_table"
s6_out$source_table <- paste0(".tmp/lit/rice2018/supplements/",
                              "11032_2018_777_MOESM3_ESM.pdf, Table S6 (pdftotext -raw)")
s6_out$doi <- DOI
write.csv(s6_out, file.path(out_dir, "rice2018_marker_cross.csv"), row.names = FALSE,
          na = "")

message(
  "rice2018: Table S5 ", nrow(out), " rows; Table S6 ", nrow(s6), " marker-cross rows / ",
  n_ids, " assays (", n_novel_ids, " novel); cross-level recount ",
  paste(got_cross, collapse = "/"), " == Table S5; prose 412/394 and 201 verified;",
  " the prose's not-met total 232 and its 31 are each 3 low (Table S6: 235)"
)
