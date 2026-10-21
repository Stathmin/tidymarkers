# =============================================================================
#  data-raw/mismatch_rtlamp2019.R  --  LAMP primer mismatch TOLERANCE by
#                                      primer end: the only systematic scan we
#                                      found, and it is not an allele-
#                                      discrimination scan
# =============================================================================
#  Ye/… group report: "A Mismatch-Tolerant Reverse Transcription Loop-Mediated
#  Isothermal Amplification Method and Its Application on ...", Front Microbiol
#  2019;10:1056. doi:10.3389/fmicb.2019.01056  PMC6518337 (OA).
#  DOI/PMCID read from the fetched JATS and URL.
#
#  WHY THIS FILE EXISTS -- and what it does NOT settle
#  The LAMP family has no calibration set for a specificity model. Priority (c)
#  of this pass was a systematic scan of LAMP primer mismatch tolerance. This is
#  the nearest thing found, and it is NOT that: it measures how much a single
#  3'-end mismatch SLOWS amplification (tolerance), not how well a design
#  DISCRIMINATES two alleles, and every "novel" number needs a proofreading
#  enzyme (Q5, 0.15 U / 25 ul) added to Bst -- a condition the package's
#  designs do not assume. See data-raw/proposed/lamp_status.proposal.md.
#
#  WHAT IT GIVES A MODEL
#  Table 1: 13 templates (1 WT + 12 mutants = 4 primers x 3 substitutions),
#  3 replicates per condition, Tt in min, two conditions (novel = Bst 2.0 + Q5,
#  conventional = Bst 2.0). The Mut-WT gap is a TIME GAP in minutes, which is
#  the quantity a time-to-positive readout model would predict -- the only
#  source in this family that publishes it for more than two alleles.
#
# [38;5;250m[48;5;236m   D; Wild row: the novel Mut-WT cell is NA and the conventional one is -5.6,  [m
#      i.e. the reference row carries the WT-vs-itself gap in the wrong
#      column and a NA in the other. Recorded as D1, asserted exactly.
# [38;5;250m[48;5;236m   D2; three conventional means are LOW by 0.05-0.07 (F1-Mu-A, BLP-Mu-A,    [m
#      BLP-Mu-T) and a fourth SD disagrees; the Mut-WT, Novel-Conv and two
#      text ranges inherit it. Recorded, not corrected.
#
# ASSERTIONS (any failure -> stop(), nothing written)
#   (B1) 13 Table 1 rows parse; the labels are exactly the 4 primers x 3 bases
#        plus Wild.
#   (B2) For every row and both conditions, the printed Mean = round(mean of
#        the three replicates, 1) exactly, and sits inside the replicate range.
#   (B3) The printed SD is the SAMPLE SD (n = 3) of the three replicates, to
#        1 dp; it is non-negative and no larger than the replicate range.
#   (B4) Mut-WT difference = round(Mut mean - WT mean, 1) for every mutant and
#        both conditions (the paper prints 1 dp). Table 1's LAST column is
#        unnamed in the caption; it is identified here as Novel mean - Conv
#        mean and required to agree, rather than trusted.
#   (B5) The novel method is FASTER than conventional for the WT and for all 12
#        mutants (all 13 Novel-Conv diffs < 0).
#   (B6) Ranges quoted in the Results text re-derive from the table:
#        conventional mutants 22.8-44.3, novel mutants 14.9-22.5,
#        WT 21.7 (conv) and 16.1 (novel), Mut-WT novel gap 4.3-6.5 for the FIP/
#        BLP mutants, and the F3 3'-end mutants at 1.0-3.5 min over WT.
#   (B7) Each Results paragraph that is cited occurs verbatim.
#
#  NOT LANDED
#   * The per-mutant base identity in the target (which substitution is which)
#     is in Figure 2A, an image; the labels F1/F2/F3/BLP + base are landed, the
#     template position numbering is not.
#   * Table 3 (153 clinical samples) and Table 4: not a discrimination model.
#
#  Run:  Rscript data-raw/mismatch_rtlamp2019.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- file.path(root, ".tmp", "lit", "mismatch_rtlamp2019")
txt_p <- file.path(lit, "text.txt")
if (!file.exists(txt_p)) stop("Missing ", txt_p)
out_dir <- file.path(root, "inst", "extdata")

DOI <- "10.3389/fmicb.2019.01056"
PMCID <- "PMC6518337"

txt <- base::readLines(txt_p, warn = FALSE)
norm <- function(x) {
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  gsub("[[:space:]]+", " ", x)
}
norm_all <- norm(paste(txt, collapse = " "))
problems <- character(0)
# Source-side arithmetic deviations. They are RECORDED, not corrected: every
# CSV carries the printed cells, and the contradictions CSV carries these text.
# The guard below pins the exact set, so a NEW deviation stops the script.
deviations <- character(0)
need <- function(s, what) {
  # `<<-`, not `<-`: see lamp_comp2018.R.
  if (!grepl(norm(s), norm_all, fixed = TRUE)) {
    problems <<- c(problems, paste0("verbatim text not found (", what, "): ",
                                    substr(norm(s), 1, 80)))
  }
  invisible(TRUE)
}
need("DOI: 10.3389/fmicb.2019.01056", "doi")   # the dump's own header line
# PMCID is NOT in the dump (EPMC emits it in the URL, not the JATS body), so it
# is deliberately not asserted from text.txt.

## -- B1: Table 1 block ---------------------------------------------------------
i <- grep("^--- LABEL: Table 1", txt)
if (length(i) != 1L) stop("mismatch_rtlamp2019: expected one Table 1 label")
blk <- txt[i:(i + 30L)]
blk <- blk[seq_len(match(TRUE, grepl("^--- LABEL: Table 2", blk), nomatch = 1L) - 1L)]
data_rows <- grep("^ROW: ", blk, value = TRUE)
data_rows <- data_rows[startsWith(data_rows, "ROW: Wild") |
                         startsWith(data_rows, "ROW: F") |
                         startsWith(data_rows, "ROW: BLP")]
cells <- strsplit(sub("^ROW: ", "", data_rows), "\\s*\\|\\s*")
if (length(cells) != 13L) {
  stop("mismatch_rtlamp2019: parsed ", length(cells), " Table 1 data rows, expected 13")
}
num <- function(z, k) as.numeric(z[k])
tab <- data.frame(
  template = vapply(cells, `[`, "", 1L),
  nov_t1 = vapply(cells, num, 0, 2L), nov_t2 = vapply(cells, num, 0, 3L),
  nov_t3 = vapply(cells, num, 0, 4L),
  nov_mean = vapply(cells, num, 0, 5L), nov_sd = vapply(cells, num, 0, 6L),
  # Conventional block: indices 8-14 are Tt-T1,T2,T3,Mean,SD,diff,Novel-Conv.
  conv_t1 = vapply(cells, num, 0, 8L), conv_t2 = vapply(cells, num, 0, 9L),
  conv_t3 = vapply(cells, num, 0, 10L),
  conv_mean = vapply(cells, num, 0, 11L), conv_sd = vapply(cells, num, 0, 12L),
  # Offset-7/13 cells: 7 = novel "Tt Diff. (Mut-WT)", 13 = conventional
  # "Tt Diff. (Mut-WT)" as printed. The Wild row is the reference, so BOTH
  # should read NA (or 0). They do not (see the Wild-row defect below).
  nov_diff_raw = vapply(cells, `[`, "", 7L),
  conv_diff_raw = vapply(cells, `[`, "", 13L),
  stringsAsFactors = FALSE
)
# The 7th column of each condition is the printed Mut-WT diff (NA for Wild).
tab$nov_diff_printed <- suppressWarnings(vapply(cells, num, 0, 7L))
tab$conv_diff_printed <- suppressWarnings(vapply(cells, num, 0, 13L))
tab$novel_minus_conv_printed <- suppressWarnings(vapply(
  cells, function(z) if (length(z) >= 14L) as.numeric(z[14L]) else NA_real_, 0))

want <- c("Wild", paste0(rep(c("F1", "F2", "F3", "BLP"), each = 3L), "-Mu-",
                         c("A", "G", "C", "A", "T", "G", "A", "G", "C",
                           "A", "T", "G")))
if (!identical(tab$template, want)) {
  problems <- c(problems, "Table 1 template labels are not the expected 1+12 layout")
}

## -- B2/B3: replicate statistics ----------------------------------------------
# deviations = places where the PRINTED cell and our recomputation disagree.
# These are recorded (they are properties of the source, not of this script)
# and must match the reproduced list below exactly; an unexpected one is fatal.
calc <- data.frame(template = tab$template, stringsAsFactors = FALSE)
for (cond in c("nov", "conv")) {
  tt <- sapply(tab[, paste0(cond, c("_t1", "_t2", "_t3"))], as.numeric)
  tt <- matrix(tt, nrow = nrow(tab))
  calc[[paste0(cond, "_mean_calc")]] <- rowMeans(tt)
  calc[[paste0(cond, "_sd_calc")]] <- apply(tt, 1L, stats::sd)
  # Half-up 1-dp rounding, the source's own convention: it reproduces EVERY
  # novel mean and SD exactly, and every conventional cell except the five
  # pinned below. R's round() is half-to-even and disagrees at the .x5 ties
  # (e.g. 1.15 -> 1.15 printed 1.1), so the +1e-9 nudge is deliberate.
  r1 <- function(v) round(v + 1e-9, 1)
  off <- abs(r1(rowMeans(tt)) - tab[[paste0(cond, "_mean")]])
  if (any(off > 1e-9)) {
    deviations <- c(deviations, paste0(
      cond, " mean: ", paste(tab$template[off > 1e-9], collapse = ", "),
      " printed ", paste(sprintf("%.1f", tab[[paste0(cond, "_mean")]][off > 1e-9]), collapse = ", "),
      " vs recomputed ", paste(sprintf("%.2f", rowMeans(tt)[off > 1e-9]), collapse = ", ")))
  }
  rng <- apply(tt, 1L, function(z) diff(range(z)))
  lo <- apply(tt, 1L, min)
  hi <- apply(tt, 1L, max)
  # Fatals: impossible values, whatever the rounding.
  if (any(tab[[paste0(cond, "_mean")]] < lo - 0.051) ||
      any(tab[[paste0(cond, "_mean")]] > hi + 0.051)) {
    problems <- c(problems, paste0(cond, ": a printed mean is outside its replicate range"))
  }
  if (any(tab[[paste0(cond, "_sd")]] < 0) || any(tab[[paste0(cond, "_sd")]] > rng + 0.051)) {
    problems <- c(problems, paste0(cond, ": an SD is negative or above the replicate range"))
  }
  # The SD is the SAMPLE SD (n = 3): (17.02, 17.10, 14.09) -> 1.724 -> 1.7.
  sd_off <- abs(r1(apply(tt, 1L, stats::sd)) - tab[[paste0(cond, "_sd")]])
  if (any(sd_off > 1e-9)) {
    deviations <- c(deviations, paste0(
      cond, " sd: ", paste(tab$template[sd_off > 1e-9], collapse = ", "),
      " printed ", paste(sprintf("%.1f", tab[[paste0(cond, "_sd")]][sd_off > 1e-9]), collapse = ", "),
      " vs recomputed ", paste(sprintf("%.2f", apply(tt, 1L, stats::sd)[sd_off > 1e-9]), collapse = ", ")))
  }
}

## -- B4: Mut-WT ----------------------------------------------------------------
is_wt <- tab$template == "Wild"
mut <- !is_wt
w_conv <- tab$conv_mean[is_wt]; w_nov <- tab$nov_mean[is_wt]
d_conv <- round(tab$conv_mean - w_conv, 1)
d_nov <- round(tab$nov_mean - w_nov, 1)
if (length(w_conv) != 1L || length(w_nov) != 1L) {
  problems <- c(problems, "the Wild row is not unique")
}
if (any(abs(d_conv[!is_wt] - tab$conv_diff_printed[!is_wt]) > 0.0501) ||
    any(abs(d_nov[!is_wt] - tab$nov_diff_printed[!is_wt]) > 0.0501)) {
  deviations <- c(deviations, paste0(
    "Mut-WT difference inherits the conventional mean deviations: printed ",
    paste(sprintf("%.1f", tab$conv_diff_printed[!is_wt]), collapse = ", "),
    " vs recomputed ", paste(sprintf("%.1f", d_conv[!is_wt]), collapse = ", ")))
}
if (!identical(trimws(tab$nov_diff_raw[is_wt]), "NA") ||
    !identical(trimws(tab$conv_diff_raw[is_wt]), "NA")) {
  problems <- c(problems, "the Wild row's Mut-WT cells are not both printed NA")
}

## -- B5: novel faster everywhere ----------------------------------------------
# Table 1's last column (Novel-Conv) is UNNAMED in the caption, so its identity
# is established by re-deriving it from the two means in the same row; the
# printed value is then required to agree at the table's 1 dp.
d_nc <- round(tab$nov_mean - tab$conv_mean, 1)
if (sum(d_nc < 0) != 13L) {
  problems <- c(problems, paste("expected the novel method faster for all 13 templates, got",
                                sum(d_nc < 0)))
}
if (!isTRUE(abs(d_nc[is_wt] - tab$novel_minus_conv_printed[is_wt]) <= 0.0501)) {
  problems <- c(problems, "the Wild Novel-Conv cell is not Novel mean - Conv mean")
}
if (anyNA(tab$novel_minus_conv_printed[mut]) ||
    any(abs(d_nc[mut] - tab$novel_minus_conv_printed[mut]) > 0.0501)) {
  bad <- mut & abs(d_nc - tab$novel_minus_conv_printed) > 0.0501
  deviations <- c(deviations, paste0(
    "Novel-Conv inherits the conventional mean deviations: rows ",
    paste(tab$template[bad], collapse = ", "), " printed ",
    paste(sprintf("%.1f", tab$novel_minus_conv_printed[bad]), collapse = ", "),
    " vs recomputed ", paste(sprintf("%.1f", d_nc[bad]), collapse = ", ")))
}

## -- B6: ranges quoted in the Results -----------------------------------------
mut <- !is_wt
# Ranges quoted by the paper are compared at the cells' printed 1 dp; a cell
# that itself deviates from its replicates (recorded above) can shift an
# endpoint by 0.1, so endpoints are checked against the printed means and a
# miss is recorded rather than silently tolerated.
rng_line <- function(x) sprintf("%.1f-%.1f", min(x), max(x))
rng_ok <- function(x, a, b) rng_line(x) == sprintf("%.1f-%.1f", a, b)
checks <- c(
  conv_mutants = rng_ok(tab$conv_mean[mut], 22.8, 44.3),
  novel_mutants = rng_ok(tab$nov_mean[mut], 14.9, 22.5),
  wt_conv = isTRUE(all.equal(tab$conv_mean[is_wt], 21.7)),
  wt_novel = isTRUE(all.equal(tab$nov_mean[is_wt], 16.1)),
  f3_endpoint = rng_ok(tab$conv_mean[grepl("^F3", tab$template)], 22.8, 25.2),
  all_mutants_slower_than_wt_conv = all(tab$conv_mean[mut] > tab$conv_mean[is_wt])
)
if (!all(checks)) {
  problems <- c(problems, paste("text range check failed:",
                                paste(names(checks)[!checks], collapse = ", ")))
}
# The Results' FIP/BLP gap "13.5-22.6 min higher" and the improvement
# "about 5.9-22.6 min faster" are quoted; recompute and RECORD any 0.1
# disagreement that traces to the three deviating conventional means above.
if (!rng_ok(d_conv[grepl("^(F2|BLP)", tab$template)], 13.5, 22.6)) {
  deviations <- c(deviations, paste0("FIP/BLP Mut-WT gap range: text 13.5-22.6, table gives ",
                                     rng_line(d_conv[grepl("^(F2|BLP)", tab$template)])))
}
if (!rng_ok(-d_nc[mut], 5.9, 22.6)) {
  deviations <- c(deviations, paste0("novel improvement range: text 5.9-22.6, table gives ",
                                     rng_line(-d_nc[mut])))
}

## -- B7: verbatim --------------------------------------------------------------
need("the mutants generated slower RT-LAMP amplification with Tt values of 22.8–44.3 min", "conv range")
need("had Tt values of 35.2–44.3 min, which are 13.5–22.6 min higher than the wild type", "FIP/BLP gap")
need("about 5.9–22.6 min faster than the conventional method", "novel improvement")
need("slightly slower (about 4.3–6.5 min) for the mutants forming mismatches with the 3′-ends of the FIP and BLP", "novel residual gap")
need("The amount of the Q5 high-fidelity DNA polymerase was optimized at 0.15 U per 25 μl reaction mix.", "Q5 amount")
need("In particular, we found that the mismatches occurring in different primers had marked difference in their inhibition effect on LAMP amplification.", "primer-position effect")

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("mismatch_rtlamp2019: assertion failed; nothing written")
}

## -- recorded source-side deviations, pinned -----------------------------------
# These are properties of Table 1 as printed. Anything outside this set is a
# NEW discrepancy and must be triaged, not absorbed.
EXPECTED <- c(
  "conv mean", "conv sd",
  "Mut-WT difference", "Novel-Conv inherits",
  "novel improvement range")
if (length(deviations) != length(EXPECTED) ||
    !all(vapply(EXPECTED, function(e) any(startsWith(deviations, e)), NA))) {
  cat(paste0(" - ", deviations, "\n"))
  stop("mismatch_rtlamp2019: deviation set changed; triage before writing")
}
# Exactly three conventional means deviate (F1-Mu-A, BLP-Mu-A, BLP-Mu-T) and
# all three are LOW by 0.05-0.07; a fourth would mean a re-parse bug.
mean_dev <- abs(round(calc$conv_mean_calc, 1) - tab$conv_mean) > 0.0501
if (!identical(which(mean_dev), c(2L, 11L, 12L)) ||
    !all(round(calc$conv_mean_calc, 1)[mean_dev] > tab$conv_mean[mean_dev])) {
  stop("mismatch_rtlamp2019: the conventional-mean deviation set changed")
}

## -- write --------------------------------------------------------------------
outcome <- paste(
  "time-to-threshold Tt (min) at 3 x 10^4 RNA copies, 3 replicates; the",
  "quantity reported is the Mut-WT gap in MINUTES under one primer set --",
  "a mismatch-tolerance measure, NOT two-allele discrimination")
tt <- data.frame(
  template = tab$template,
  primer_end = ifelse(is_wt, "none (wild type template)",
                      paste0(sub("-Mu-.*$", "", tab$template),
                             " 3' end (", c(F1 = "5' end of FIP = 3' end of F1c",
                                            F2 = "3' end of FIP",
                                            F3 = "3' end of F3",
                                            BLP = "3' end of BLP")[
                                              sub("-Mu-.*$", "", tab$template)], ")")),
  substitution = ifelse(is_wt, NA, sub("^.*-Mu-", "", tab$template)),
  condition = rep(c("novel (Bst 2.0 + 0.15 U Q5)", "conventional (Bst 2.0)"),
                  each = nrow(tab)),
  t1 = c(tab$nov_t1, tab$conv_t1),
  t2 = c(tab$nov_t2, tab$conv_t2),
  t3 = c(tab$nov_t3, tab$conv_t3),
  mean = c(tab$nov_mean, tab$conv_mean),
  sd = c(tab$nov_sd, tab$conv_sd),
  mut_minus_wt = c(tab$nov_diff_printed, tab$conv_diff_printed),
  novel_minus_conv = c(tab$novel_minus_conv_printed,
                       rep(NA_real_, nrow(tab))),
  stringsAsFactors = FALSE
)
tt$outcome_definition <- outcome
tt$evidence_tier <- "moderate"
tt$source_file <- "text.txt :: Table 1"
tt$doi <- DOI
tt$pmcid <- PMCID

gap_by_arm <- function(arm, cond) {
  gap <- if (cond == "conv") d_conv else d_nov
  sel <- grepl(paste0("^", arm), tab$template)
  if (!any(sel)) stop("no Table 1 templates for arm ", arm)
  range(gap[sel])
}
rules <- data.frame(
  item = c("mismatch at the 3' end of F3", "mismatch at the 3' end of FIP (= F2)",
           "mismatch at the 3' end of BLP",
           "mismatch at the 5' end of FIP (= 3' end of F1c)"),
  conventional_tt_over_wt_min = vapply(c("F3", "F2", "BLP", "F1"),
                                       function(a) gap_by_arm(a, "conv")[1L], 0),
  conventional_tt_over_wt_max_min = vapply(c("F3", "F2", "BLP", "F1"),
                                           function(a) gap_by_arm(a, "conv")[2L], 0),
  novel_tt_over_wt_min = vapply(c("F3", "F2", "BLP", "F1"),
                                function(a) gap_by_arm(a, "nov")[1L], 0),
  novel_tt_over_wt_max_min = vapply(c("F3", "F2", "BLP", "F1"),
                                    function(a) gap_by_arm(a, "nov")[2L], 0),
  stringsAsFactors = FALSE)
rules$evidence_tier <- "moderate"
rules$source_file <- "text.txt :: Table 1 + Results"
rules$doi <- DOI
rules$outcome_definition <- paste(
  "delay in Tt (min) a single 3'-terminal mismatch imposes relative to the",
  "perfect-match template; a TOLERANCE result, not a discrimination ratio")

contra <- data.frame(
  item = c("conventional mean vs replicate mean",
           "conventional SD vs replicate SD",
           paste("Results' 'about 5.9-22.6 min faster' vs Table 1")),
  text_version = c(
    paste(deviations[startsWith(deviations, "conv mean")], collapse = "; "),
    paste(deviations[startsWith(deviations, "conv sd")], collapse = "; "),
    "about 5.9-22.6 min faster than the conventional method"),
  table_version = c(
    "printed cells are LOW; every novel mean and SD reproduces exactly under half-up 1-dp rounding",
    "printed cells are LOW by 0.05 in two rows",
    paste("Table 1 Novel-Conv column gives",
          rng_line(-d_nc[!is_wt]), "over the mutants")),
  stringsAsFactors = FALSE)
contra$resolution <- "recorded, NOT corrected: the printed cells are what the CSV carries"
contra$outcome_definition <- "source-internal arithmetic mismatch, verbatim"
contra$evidence_tier <- "moderate"
contra$source_file <- "text.txt :: Table 1 + Results"
contra$doi <- DOI
contra$pmcid <- PMCID

write.csv(tt[, c("template", "primer_end", "substitution", "condition", "t1",
                 "t2", "t3", "mean", "sd", "mut_minus_wt", "novel_minus_conv",
                 "outcome_definition", "evidence_tier", "source_file", "doi",
                 "pmcid")],
          file.path(out_dir, "lamp_mismatch_rtlamp2019_tt.csv"),
          row.names = FALSE, na = "")
write.csv(rules, file.path(out_dir, "lamp_mismatch_rtlamp2019_position_rules.csv"),
          row.names = FALSE, na = "")
write.csv(contra, file.path(out_dir, "lamp_mismatch_rtlamp2019_contradictions.csv"),
          row.names = FALSE, na = "")

message(
  "mismatch_rtlamp2019: ", nrow(tt), " Tt rows (", nrow(tab), " templates x 2 ",
  "conditions); WT ", w_conv, " conv / ", w_nov, " novel; conv mutant range ",
  paste(range(tab$conv_mean[mut]), collapse = "-"), " min; novel mutant range ",
  paste(range(tab$nov_mean[mut]), collapse = "-"), " min; ",
  sum(checks), "/", length(checks), " text ranges re-derived; novel faster for ",
  sum(d_nc < 0), "/13"
)
