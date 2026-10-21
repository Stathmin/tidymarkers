# =============================================================================
#  data-raw/asrpa_sickle.R  --  allele-specific RPA calls against IEF/HPLC
# =============================================================================
#  Allele-specific RPA for sickle-cell genotyping (HbS, rs334).
#  Anal Chem 2020. doi:10.1021/acs.analchem.0c04191  PMC7992048 (CC-BY)
#
#  WHY THIS FILE EXISTS
#  A second, independent field check of the same design idea as asrpa_mdr: a
#  primer placed so that its 3' region discriminates betaA from betaS. The
#  quantity is again a CALL RATE against a reference method, not a mismatch
#  effect, and the CSV says so.
#
#  The paper's design statement is that a mismatch in the PENULTIMATE position
#  carries the discrimination when the terminal base cannot (A/T transition:
#  the terminal position of one allele's primer is a perfect match to the
#  other allele). We record what the assay achieved; we do NOT derive a
#  penultimate weight from it -- one genotype cannot calibrate a weight, and
#  rpa_mismatch_terminal.csv is the calibration set.
#
#  THE ASSERTIONS
#  * both Table 1 denominators are 34 + 57 = 91;
#  * Table 1's percentages reproduce from its printed fractions;
#  * Table 2's genotyping matrix has 91 calls, and its two per-genotype rows
#    reproduce the same 33/34, 54/57 the sensitivity table prints;
#  * the reported accuracy 87/91 agrees with the diagonal of Table 2.
#
#  Run:  Rscript data-raw/asrpa_sickle.R [--dir .tmp/lit/asrpa_sickle]
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
  file.path(root, ".tmp", "lit", "asrpa_sickle")
ft <- file.path(lit, "fulltext.xml")
if (!file.exists(ft)) {
  stop("missing ", ft, "\n  fetch with: Rscript data-raw/fetch-sources.R")
}
DOI <- "10.1021/acs.analchem.0c04191"

doc <- read_xml(ft)
rows_of <- function(label) {
  t <- xml_find_all(doc, sprintf("//body//table-wrap[label='%s']", label))
  if (length(t) != 1L) stop("full text: expected one ", label)
  lapply(xml_find_all(t, ".//tr"), function(r) {
    trimws(gsub("[[:space:]]+", " ",
                xml_text(xml_find_all(r, "./td|./th"))))
  })
}
# "100.0% (34/34)" -> c(pct = 100.0, n = 34, d = 34)
frac <- function(cell) {
  # The publisher uses a non-ASCII space between the % and the bracket, so
  # match on the token shape rather than on an exact ASCII string.
  m <- regmatches(cell, regexec(
    "([0-9.]+)%[^0-9]+\\(([0-9]+)/([0-9]+)\\)", cell))[[1L]]
  if (length(m) != 4L) stop("unparsable fraction cell: ", cell)
  c(pct = as.numeric(m[2L]), n = as.integer(m[3L]), d = as.integer(m[4L]))
}
t1 <- rows_of("Table 1")
bA <- lapply(t1[[2L]][2:4], frac)
bS <- lapply(t1[[3L]][2:4], frac)
names(bA) <- names(bS) <- c("sens", "spec", "acc")

## --- Table 1 arithmetic -----------------------------------------------------
for (a in list(list(lab = "betaA", x = bA), list(lab = "betaS", x = bS))) {
  x <- a$x
  # betaA: sensitivity over 34 alleles, specificity over 57; betaS: 57 / 34.
  d_expect <- if (a$lab == "betaA") c(34L, 57L) else c(57L, 34L)
  if (x$sens[["d"]] != d_expect[1L] || x$spec[["d"]] != d_expect[2L]) {
    stop("Table 1 ", a$lab, ": denominators are not ", d_expect[1L], "/",
         d_expect[2L])
  }
  if (x$sens[["d"]] + x$spec[["d"]] != 91L) stop("34 + 57 != 91")
  for (k in names(x)) {
    if (abs(round(100 * x[[k]][["n"]] / x[[k]][["d"]], 1) -
            x[[k]][["pct"]]) > 0.05) {
      stop("Table 1 ", a$lab, " ", k, ": printed % does not reproduce")
    }
  }
  if (x$acc[["d"]] != 91L) stop("Table 1 ", a$lab, ": accuracy is not over 91")
}

## --- Table 2: the genotyping matrix -----------------------------------------
t2 <- rows_of("Table 2")
# Row 1 is the IEF/HPLC header plus the per-genotype summary columns; the last
# four rows are the 3x3 called-vs-reference matrix and its total row.
raw <- rows_of("Table 2")
# The matrix rows are ragged (the AA row also carries the first summary block),
# so slice them positionally from the row list rather than by label: row 3 is
# AA, row 4 AS, row 5 SS, row 6 the total row, after the two header rows.
if (length(raw) < 6L) stop("Table 2: fewer rows than expected")
# AA is called 33 times, so its row stays aligned with raw[[4]] (AS) when the
# two leading label cells are dropped.
first3 <- function(r) suppressWarnings(as.integer(r[2:4]))
aa <- suppressWarnings(as.integer(raw[[3L]][3:5]))
cm <- rbind(AA = aa, AS = first3(raw[[4L]]), SS = first3(raw[[5L]]))
colnames(cm) <- c("AA", "AS", "SS")
rows <- raw[[6L]]
if (!identical(rows[1L], "total")) stop("Table 2: row 6 is not the total row")
cm_total <- suppressWarnings(as.integer(rows[2:5]))
if (anyNA(cm) || anyNA(cm_total)) stop("Table 2: unparsable matrix cell")
if (sum(cm) != 91L) stop("Table 2: matrix sums to ", sum(cm), ", not 91")
# identical() compares names, and colSums keeps the matrix's dimnames: compare
# the numbers.
if (!all(colSums(cm) == unname(cm_total[1:3])) ||
    cm_total[4L] != 91L) {
  stop("Table 2: column totals do not reproduce from the matrix")
}
if (any(diag(cm) != c(33L, 0L, 54L)) || sum(diag(cm)) != 87L) {
  stop("Table 2: diagonal is not 33 / 0 / 54 = 87")
}
# The per-genotype summary columns of Table 2 must repeat Table 1's numbers.
# The AA block is the tail of matrix row 1; the SS block the tail of row 3.
s_AA <- lapply(raw[[3L]][8:10], frac)   # AA: sens, spec, accuracy
s_SS <- lapply(raw[[5L]][7:8], frac)    # SS: sens, spec
if (s_AA[[1L]][["n"]] != 33L || s_AA[[1L]][["d"]] != 34L ||
    s_AA[[2L]][["n"]] != 57L || s_AA[[2L]][["d"]] != 57L ||
    s_AA[[3L]][["n"]] != 87L || s_AA[[3L]][["d"]] != 91L) {
  stop("Table 2 summary row AA does not repeat Table 1")
}
if (s_SS[[1L]][["n"]] != 54L || s_SS[[1L]][["d"]] != 57L ||
    s_SS[[2L]][["n"]] != 34L || s_SS[[2L]][["d"]] != 34L) {
  stop("Table 2 summary row SS does not repeat Table 1")
}

## --- output -----------------------------------------------------------------
out <- data.frame(
  allele = c("betaA (HbA)", "betaS (HbS)"),
  sensitivity_pct = c(bA$sens[["pct"]], bS$sens[["pct"]]),
  sensitivity_n = c(bA$sens[["n"]], bS$sens[["n"]]),
  sensitivity_d = c(bA$sens[["d"]], bS$sens[["d"]]),
  specificity_pct = c(bA$spec[["pct"]], bS$spec[["pct"]]),
  specificity_n = c(bA$spec[["n"]], bS$spec[["n"]]),
  specificity_d = c(bA$spec[["d"]], bS$spec[["d"]]),
  accuracy_pct = c(bA$acc[["pct"]], bS$acc[["pct"]]),
  accuracy_n = c(bA$acc[["n"]], bS$acc[["n"]]),
  accuracy_d = c(bA$acc[["d"]], bS$acc[["d"]]),
  n_samples = 91L,
  n_samples_excluded = 109L - 91L,
  stringsAsFactors = FALSE)
out$outcome_definition <- paste(
  "CALL AGREEMENT, not a mismatch effect: allele-level sensitivity and",
  "specificity of an allele-specific RPA assay against IEF/HPLC genotype,",
  "reported by the paper as n / d and recomputed here. There is no fold change",
  "and no time-to-positive difference to land. n_samples_excluded = 109",
  "enrolled minus 91 genotyped in Table 2; the paper's own denominator for",
  "accuracy is 91.")
out$evidence_tier <- "primary table (Tables 1-2, calls vs reference genotype)"
out$source_file <- "fulltext.xml"
out$doi <- DOI
write.csv(out, file.path(root, "inst", "extdata", "rpa_asrpa_sickle_calls.csv"),
          row.names = FALSE, na = "")

message("asrpa_sickle: betaA sens ", out$sensitivity_pct[1L], "% (",
        out$sensitivity_n[1L], "/", out$sensitivity_d[1L], "), betaS spec ",
        out$specificity_pct[2L], "% (", out$specificity_n[2L], "/",
        out$specificity_d[2L], "); diagonal ", sum(diag(cm)), "/",
        sum(cm))
