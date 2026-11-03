# =============================================================================
#  data-raw/readout_lamp.R  --  LAMP kinetics: conditions -> time-to-positive
# =============================================================================
#  Source (DOI is PROVENANCE; every number is re-read from the fetched full text):
#
#    Barone MM, Moizio M, Choudhary R, D'Errico C, Trkulja V, Torta L, Davino S,
#    Mati...  Development of a Field-Deployable Loop-Mediated Isothermal
#    Amplification Assay for the Rapid Detection of Erysiphe corylacearum.
#    J Fungi (Basel) 2026.  DOI 10.3390/jof12010079   EPMC PMC12842720
#    Fetched: .tmp/a6/PMC12842720.xml (Europe PMC fullTextXML), Tables 2 and 3
#    read from the raw <table-wrap> bodies, not from the flattened text.
#
#  WHY THIS SET EXISTS
#  ROADMAP §V marks LAMP/RPA kinetics unbuildable (Q1: no "conditions -> outcome"
#  set).  This paper prints two such grids from ONE real-time LAMP run:
#   * Table 2 - Tp (min) for 3 templates x 3 isothermal temperatures (60/63/65 C),
#     all nine cells detected, so the temperature dependence is not censored;
#   * Table 3 - Tp +/- SD for 3 templates x a stated ten-fold mass dilution
#     series (20 ng/uL down to 10^-6), with the 10^-6 cells printed "ND".
#  Reaction: 12.5 uL, 6.25 uL 2x SuperScript IV RT-LAMP master mix, 1.5 uL DNA
#  (35 ng), SYTO dye, run 60 min, 2 technical replicates per sample.
#
#  WHAT STAYS CENSORED / UNKNOWN
#   * "ND" (not determined) at 10^-6 is NOT a time; it is recorded with
#     censoring = "not_detected_within_run" and time_to_positive_min = NA, so a
#     kinetic model must treat it as >= 60 min (the stated run length), never as
#     a point.
#   * Table 3's dilutions are MASS dilutions of an unquantified extract.  The
#     copies per reaction are NOT printed, so this set constrains Tp vs relative
#     template amount, not Tp vs absolute copies.  That is the Q1 gap and it is
#     not closed here.
#   * Mg2+ is not printed (it lives in the master mix); salt is NA.
#   * SD is printed but n is not; n = 2 is stated in Methods, so the SD is a
#     2-replicate SD and is stored as printed, not recomputed.
#
#  ASSERTIONS (the source's own arithmetic)
#   * Table 2: nine Tp cells for 3 x 3; every template x temperature present.
#   * Table 2: for every template, Tp(60) < Tp(63) - the paper's own "higher
#     annealing temperature slows amplification" ordering.  Its 63 -> 65 C step
#     is NOT monotone (8.2_PdIT: 37 then 35 min), and is kept as printed.
#   * Table 3: within every template the Tp is non-decreasing along the
#     dilution axis (undiluted -> 10^-5), and 10^-6 is ND in all three rows.
#   * Table 3: the paper's own stated limit of detection (10^-5 for the two
#     isolates) means the last detected dilution of each isolate is 10^-5.
#
#  Run:  Rscript data-raw/readout_lamp.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

doi <- "10.3390/jof12010079"
pmcid <- "PMC12842720"
source_file <- "fulltext.xml :: Table 2 (Tp vs isothermal temperature), Table 3 (Tp vs ten-fold mass dilution)"

run_min <- 60 # Methods: "LAMP assay was performed for 60 min"

# --- Table 2: Tp (min), 3 templates x 3 temperatures ------------------------
t2 <- data.frame(
  template = c("positive control", "8.2_PdIT", "8.3_PdIT"),
  tp_60 = c(5, 6, 6.5),
  tp_63 = c(32, 37, 40),
  tp_65 = c(34, 35, 37),
  stringsAsFactors = FALSE
)
stopifnot(nrow(t2) == 3L, !anyNA(t2[, 2:4]))
# 60 -> 63 C slows every template; 63 -> 65 C does not (8.2_PdIT is 37 then
# 35), so only the first ordering is asserted and the paper's own non-monotone
# cell is kept as printed.
stopifnot(all(t2$tp_60 < t2$tp_63), all(t2$tp_63 != t2$tp_65))
stopifnot(all(unlist(t2[, 2:4]) <= run_min))

# --- Table 3: Tp +- SD vs dilution ------------------------------------------
# Printed column order: 20 ng/uL, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6.
t3 <- data.frame(
  template = rep(c("1.1_PdIT", "7.1_PdIT", "positive control"), each = 7L),
  dilution_label = rep(c("20 ng/uL", "1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6"), 3L),
  dilution_factor = rep(c(1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6), 3L),
  tp_min = c(
    4.05, 5.04, 36.24, 37.37, 45.43, 47.43, NA,
    4.13, 4.90, 36.57, 41.06, 46.99, 49.87, NA,
    4.70, 5.98, 32.32, 36.11, 39.74, 49.41, NA
  ),
  tp_sd_min = c(
    0.34, 0.80, 0.80, 1.20, 2.82, 1.26, NA,
    0.23, 1.09, 2.32, 0.82, 1.47, 0.18, NA,
    0.15, 0.05, 1.14, 4.22, 0.91, 3.42, NA
  ),
  stringsAsFactors = FALSE
)
stopifnot(nrow(t3) == 21L)
# 10^-6 is ND in all three rows and nowhere else is ND.
stopifnot(all(is.na(t3$tp_min) == (t3$dilution_label == "1e-6")))
# Every detected Tp lies inside the stated 60-min run.
stopifnot(all(t3$tp_min <= run_min, na.rm = TRUE))
# Tp is non-decreasing along the dilution axis within each template.
for (tpl in unique(t3$template)) {
  v <- t3$tp_min[t3$template == tpl]
  v <- v[!is.na(v)]
  stopifnot(all(diff(v) > 0))
}
# The paper's stated LOD: 10^-5 for the two isolates.
for (tpl in c("1.1_PdIT", "7.1_PdIT")) {
  det <- t3[t3$template == tpl & !is.na(t3$tp_min), ]
  stopifnot(identical(det$dilution_label[[nrow(det)]], "1e-5"))
}

# --- OUTPUT -----------------------------------------------------------------
dilution_axis <- c("20 ng/uL" = 0L, "1e-1" = 1L, "1e-2" = 2L, "1e-3" = 3L, "1e-4" = 4L, "1e-5" = 5L, "1e-6" = 6L)
t3$template_ng_ul <- 20 * t3$dilution_factor

out <- rbind(
  data.frame(
    chemistry = "lamp",
    experiment = "temperature",
    template = t2$template,
    temperature_c = rep(c(60, 63, 65), each = 3L),
    tp_min = c(t2$tp_60, t2$tp_63, t2$tp_65),
    tp_sd_min = NA_real_,
    dilution_log10 = NA_integer_,
    template_ng_ul = NA_real_,
    copies_per_reaction = NA_real_,
    censoring = "detected",
    stringsAsFactors = FALSE
  ),
  data.frame(
    chemistry = "lamp",
    experiment = "dilution",
    template = t3$template,
    temperature_c = 60, # Methods: all successive assays at 60 C
    tp_min = t3$tp_min,
    tp_sd_min = t3$tp_sd_min,
    dilution_log10 = as.integer(dilution_axis[t3$dilution_label]),
    template_ng_ul = t3$template_ng_ul,
    copies_per_reaction = NA_real_,
    censoring = ifelse(is.na(t3$tp_min), "not_detected_within_run", "detected"),
    stringsAsFactors = FALSE
  )
)
out$run_min <- run_min
out$n_replicates <- 2L
out$mg_mm <- NA_real_
out$outcome_definition <- paste0(
  "Outcome = the time to positive Tp (min) in ONE real-time LAMP run (SYTO dye, Rotor-Gene Q, readings every 15 s, ",
  "60 min total, 2 technical replicates). experiment='temperature' is Table 2 (Tp vs 60/63/65 C); ",
  "experiment='dilution' is Table 3 (Tp +/- SD vs a stated ten-fold MASS dilution of an extract, 20 ng/uL down to 1e-6). ",
  "ND at 1e-6 is censored at the 60-min run length and carries no time (tp_min = NA), never a point. ",
  "copies_per_reaction is NA because the paper prints no concentration of the extract, so the set constrains Tp vs RELATIVE ",
  "template amount, not vs absolute copies - the Q1 gap. mg_mm is NA (Mg is inside the master mix). ",
  "tp_sd_min is the paper's printed 2-replicate SD, stored as printed."
)
out$evidence_tier <- "measured"
out$source_file <- source_file
out$doi <- doi
out$pmcid <- pmcid

stopifnot(nrow(out) == 30L, sum(out$censoring == "not_detected_within_run") == 3L)

dir.create(file.path(root, "inst", "extdata"), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, file.path(root, "inst", "extdata", "readout_lamp_erysiphe2026.csv"), row.names = FALSE)

cat(sprintf(
  "readout_lamp_erysiphe2026.csv: %d rows (%d temperature cells, %d dilution cells, %d censored); copies/reaction not printed.\n",
  nrow(out), sum(out$experiment == "temperature"), sum(out$experiment == "dilution"),
  sum(out$censoring == "not_detected_within_run")
))
