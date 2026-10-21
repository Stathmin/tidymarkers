# =============================================================================
#  data-raw/asrpa_mdr.R  --  what allele-specific RPA actually called in a field
#  trial (MDR-TB), and how often it missed a resistant isolate
# =============================================================================
#  Wang Y, et al. "Allele-specific recombinase polymerase amplification for
#  rapid detection of multidrug-resistant tuberculosis" PLoS One 2021.
#  doi:10.1371/journal.pone.0253235  PMC8195408 (CC-BY)
#
#  WHY THIS FILE EXISTS
#  rpa_mismatch_class tells us what a GIVEN mismatch does to the reaction. This
#  paper tells us what a designed allele-specific RPA assay achieved on real
#  isolates: a call rate against sequencing and against phenotype. That is a
#  DIFFERENT QUANTITY and is labelled as such in the CSV; it is the field
#  counterpart of the model, not a second measurement of a mismatch effect.
#
#  THE ASSERTIONS
#  * every 2x2 sums to its printed N (141) and its margins;
#  * the printed sensitivity/specificity reproduce from the 2x2 at 2 dp;
#  * the rpoB "either of three" row is NOT the sum of the three rpoB rows on
#    the phenotype side (the assay is a combined call), which is exactly the
#    kind of transcription error this file is meant to catch;
#  * the full text's own statement about where the discriminating mismatch
#    sits (third base from the 3' terminus) is quoted, not paraphrased.
#
#  Run:  Rscript data-raw/asrpa_mdr.R [--dir .tmp/lit/asrpa_mdr]
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
  file.path(root, ".tmp", "lit", "asrpa_mdr")
ft <- file.path(lit, "fulltext.xml")
if (!file.exists(ft)) {
  stop("missing ", ft, "\n  fetch with: Rscript data-raw/fetch-sources.R")
}

DOI <- "10.1371/journal.pone.0253235"

doc <- read_xml(ft)
tbl <- function(label) {
  t <- xml_find_all(doc, sprintf("//body//table-wrap[label='%s']", label))
  if (length(t) != 1L) stop("full text: expected one ", label)
  rows <- lapply(xml_find_all(t, ".//tr"), function(r) {
    trimws(gsub("[[:space:]]+", " ",
                xml_text(xml_find_all(r, "./td|./th"))))
  })
  rows
}
t6 <- tbl("Table 6")
t7 <- tbl("Table 7")

## --- 2x2 reader -------------------------------------------------------------
# The tables print an assay NAME on its own line, then the "Mutation detected"
# line carrying TP/FP and the printed sensitivity/specificity, then the
# "Mutation not detected" line carrying FN/TN. Numbers are read from the cells;
# nothing is re-typed.
int_at <- function(row, k) {
  if (length(row) < k) NA_integer_
  else suppressWarnings(as.integer(row[[k]]))
}
parse_2x2 <- function(rows, names_expected) {
  out <- list()
  for (i in seq_along(rows)) {
    lab <- rows[[i]][1L]
    if (!lab %in% names_expected) next
    det <- rows[[i + 1L]]; ndet <- rows[[i + 2L]]
    if (det[1L] != "Mutation detected") stop("layout: ", lab)
    if (ndet[1L] != "Mutation not detected") stop("layout: ", lab)
    out[[lab]] <- c(
      TP = int_at(det, 2L), FP = int_at(det, 3L),
      FN = int_at(ndet, 2L), TN = int_at(ndet, 3L),
      sens_printed = as.numeric(sub("%.*$", "", det[4L])),
      spec_printed = as.numeric(sub("%.*$", "", det[5L])))
  }
  out
}
t6_names <- c("rpoB516", "rpoB526", "rpoB531", "katG315")
m6 <- parse_2x2(t6, t6_names)
m7 <- parse_2x2(t7, c("Either rpoB516 or rpoB526 or rpoB531", "katG315"))
if (length(m6) != 4L) stop("Table 6: not all four assays parsed")
if (length(m7) != 2L) stop("Table 7: not both assays parsed")

## --- arithmetic assertions --------------------------------------------------
# Every 2x2 is over the same 141 isolates (Table 4/5 report N = 141 and the
# margins below are checked against that).
N <- 141L
for (nm in names(m6)) {
  x <- m6[[nm]]
  if (sum(x[1:4]) != N) stop("Table 6 ", nm, ": 2x2 sums to ", sum(x[1:4]))
  if (abs(100 * x[["TP"]] / (x[["TP"]] + x[["FN"]]) - x[["sens_printed"]]) > 0.006) {
    stop("Table 6 ", nm, ": printed sensitivity does not reproduce")
  }
  if (abs(100 * x[["TN"]] / (x[["TN"]] + x[["FP"]]) - x[["spec_printed"]]) > 0.006) {
    stop("Table 6 ", nm, ": printed specificity does not reproduce")
  }
}
if (sum(m7[[1L]][1:4]) != N || sum(m7[[2L]][1:4]) != N) {
  stop("Table 7: a 2x2 does not sum to 141")
}
for (nm in names(m7)) {
  x <- m7[[nm]]
  if (abs(100 * x[["TP"]] / (x[["TP"]] + x[["FN"]]) - x[["sens_printed"]]) > 0.006) {
    stop("Table 7 ", nm, ": printed sensitivity does not reproduce")
  }
  if (abs(100 * x[["TN"]] / (x[["TN"]] + x[["FP"]]) - x[["spec_printed"]]) > 0.006) {
    stop("Table 7 ", nm, ": printed specificity does not reproduce")
  }
}
# Table 4's phenotype split and Table 7's margins must agree: the three rpoB
# 2x2s describe three mutations, but the combined assay is scored once, so the
# rpoB resistant margin is 54 = Table 4's total resistant count, NOT the sum of
# the three rpoB TP values.
t4 <- tbl("Table 4")
tot <- t4[[length(t4)]]
if (!any(grepl("Total", tot))) stop("Table 4: no total row")
r4 <- c(suppressWarnings(as.integer(tot[3L])), suppressWarnings(as.integer(tot[4L])))
if (!identical(r4, c(87L, 54L))) stop("Table 4 total row is not 87 / 54: ",
                                      paste(tot, collapse = " | "))
x <- m7[[1L]]
if (x[["TP"]] + x[["FN"]] != r4[2L]) stop("rpoB resistant margin != Table 4")
if (x[["FP"]] + x[["TN"]] != r4[1L]) stop("rpoB susceptible margin != Table 4")
rpo_sum <- sum(vapply(m6[c("rpoB516", "rpoB526", "rpoB531")],
                      function(z) z[["TP"]], numeric(1)))
if (!rpo_sum <= x[["TP"]]) stop("per-mutation rpoB TP exceeds the combined call")

## --- the discriminating-mismatch position, quoted from the full text -------
body <- gsub("[[:space:]]+", " ", xml_text(xml_find_all(doc, "//body//p")))
q_pos <- grep("third base from the 3", body, fixed = TRUE, value = TRUE)
if (length(q_pos) != 1L) {
  stop("full text: the 'third base from the 3' terminus' statement moved")
}

## --- output -----------------------------------------------------------------
mk <- function(x, assay, reference, mutation) {
  data.frame(
    assay = assay,
    reference_standard = reference,
    mutation = mutation,
    TP = x[["TP"]], FP = x[["FP"]], FN = x[["FN"]], TN = x[["TN"]],
    N = sum(x[1:4]),
    sensitivity_pct = x[["sens_printed"]],
    specificity_pct = x[["spec_printed"]],
    # Recomputed from the 2x2, not copied: if these disagree with the printed
    # columns above the assertion would already have failed.
    sensitivity_pct_recomputed = round(100 * x[["TP"]] /
                                         (x[["TP"]] + x[["FN"]]), 2),
    specificity_pct_recomputed = round(100 * x[["TN"]] /
                                         (x[["TN"]] + x[["FP"]]), 2),
    stringsAsFactors = FALSE)
}
out <- rbind(
  mk(m6[["rpoB516"]], "AS-RPA/SYBR, rpoB516", "Sanger sequencing", "rpoB 516"),
  mk(m6[["rpoB526"]], "AS-RPA/SYBR, rpoB526", "Sanger sequencing", "rpoB 526"),
  mk(m6[["rpoB531"]], "AS-RPA/SYBR, rpoB531", "Sanger sequencing", "rpoB 531"),
  mk(m6[["katG315"]], "AS-RPA/SYBR, katG315", "Sanger sequencing", "katG 315"),
  mk(m7[[1L]], "AS-RPA/SYBR, rpoB516/526/531 (combined)",
     "phenotypic drug-susceptibility test", "rifampicin resistance"),
  mk(m7[[2L]], "AS-RPA/SYBR, katG315",
     "phenotypic drug-susceptibility test", "isoniazid resistance"))
out$missed_resistant_calls <- out$FN
out$missed_pct <- round(100 * out$FN / (out$TP + out$FN), 1)
out$outcome_definition <- paste(
  "CALL AGREEMENT, not a mismatch effect: TP/FP/FN/TN of an allele-specific",
  "RPA reaction against ", out$reference_standard,
  ". Sensitivity = TP/(TP+FN), specificity = TN/(TN+FP). Every assay here is",
  "designed so that the discriminating primer carries a deliberately placed",
  "mismatch; the paper reports the resulting call rate, NOT the log-odds cost",
  "of any single mismatch position. Do not use these rows to weight a mismatch;",
  "use rpa_mismatch_terminal.csv for that.")
out$evidence_tier <- "primary table (Tables 6-7, 2x2 against a reference)"
out$source_file <- "fulltext.xml"
out$doi <- DOI
write.csv(out, file.path(root, "inst", "extdata", "rpa_asrpa_mdr_calls.csv"),
          row.names = FALSE, na = "")

message("asrpa_mdr: ", nrow(out), " 2x2 tables; combined rpoB missed ",
        out$FN[5L], " of ", out$TP[5L] + out$FN[5L],
        " phenotypically resistant isolates (", out$missed_pct[5L], "%)")
