# =============================================================================
#  data-raw/lamp_melting_allele_spectrum2025.R  --  the one LAMP study that
#                                                   prints a QUENCHING-PROBE
#                                                   melting peak per allele
# =============================================================================
#  Imai K, et al.  "Melting temperature-based genotyping of Treponema pallidum
#  by loop-mediated isothermal amplification with quenching probes",
#  Microbiology Spectrum 2025, 10.1128/spectrum.01565-25, PMC12671098.
#  CC-BY 4.0.
#
#  WHY THIS FILE EXISTS
#  `inst/extdata/mismatch_family_lamp.csv` already carries this paper's tp0548
#  melt peaks (SS14 92.46 +/- 0.13 C, Nichols 93.71 +/- 0.21 C).  What it does
#  not carry is the rest of the readout: Supplementary Table 1 is 133 clinical
#  samples x three quenching probes (tp0705 506 / 625 / 708), each sample with
#  a measured Tm and a genotype letter.  That is the only per-sample melting
#  table in the corpus, and it is what lets a LAMP readout model be checked
#  against real calls instead of a single aggregate peak.
#
#  TWO FILES, TWO JOBS
#   * lamp_tp0705_melting_calls.csv  -- the 133 samples.  Calls are kept
#     VERBATIM, including the literal string "NA" (the paper's own cell value
#     for an unassembled call) and "-" (not typed).  Every such row still
#     carries a numeric Tm, so Tm and call are independent columns.
#   * lamp_tp0705_reaction_conditions.csv  -- the Methods paragraph: 0.2 M
#     betaine, 6 mM MgSO4, 68 C, QProbe 80 nM, EvaGreen 1.25 uL, melt ramp
#     0.2 C/s.  Same nine-column schema as
#     `lamp_snp_lamp_lfd2025_conditions.csv`; all rows `evidence_tier = "weak"`
#     because a reaction condition is not a readout outcome.
#     This is the concrete magnitude behind `SCHEME_UNMODELLED$lamp = "betaine"`
#     in R/marker_core.R, which currently names betaine without a number.
#
#  WHAT THE PAPER DOES NOT REPORT
#   * No per-probe allele-gap table; the separation has to be read off the
#     sample Tm columns (that is what the test does).
#   * No raw amplification curves, no fluorescence units.
#   * The supplement is an .xlsx; `readxl`/`openxlsx` are not installed under
#     R 4.6.1, so the sheet is converted once and read as TSV here.
#
#  FETCH (verified routes, 2026-11)
#    curl -sS -L -m 120 -o .tmp/lit/lamp2025/supp.zip \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12671098/supplementaryFiles'
#    unzip -o .tmp/lit/lamp2025/supp.zip -d .tmp/lit/lamp2025/
#    python3 -c 'import openpyxl,sys; w=openpyxl.load_workbook(sys.argv[1]).active; \
#      open(sys.argv[2],"w",encoding="utf-8",newline="").write("\n".join( \
#      "\t".join("" if c.value is None else str(c.value) for c in row) \
#      for row in w.iter_rows()))' \
#      .tmp/lit/lamp2025/spectrum.01565-25-s0003.xlsx .tmp/lit/lamp2025/t1.tsv
#
#  Run:  Rscript data-raw/lamp_melting_allele_spectrum2025.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

tsv_path <- Sys.getenv("LAMP_TSV", ".tmp/lit/lamp2025/t1.tsv")
if (!file.exists(tsv_path)) {
  stop("missing ", tsv_path, "; see the FETCH block at the top of this script")
}

raw <- readr::read_tsv(tsv_path,
  col_types = readr::cols(.default = "c"),
  na = character(),
  locale = readr::locale(encoding = "UTF-8")
)

# The header carries "Tm (C )" with a space before the paren, so select by
# position and refuse to go on if the sheet's column order has moved.
expected_names <- c(
  "No.", "Sample type", "Whole genome sequencing", "Strain name",
  "tp0136 allele", "tp0548 allele", "tp0705 allele", "ST (Grillova, et al)",
  "tp0548 sequence", "tp0856 sequence", "Species", "Lineage",
  "tp0705 506-LAMP, Tm (C )", "tp0705 506-LAMP",
  "tp0705 625-LAMP, Tm (C )", "tp0705 625-LAMP",
  "tp0705 708-LAMP, Tm (C )", "tp0705 708-LAMP",
  "tp0548 Nichols-LAMP", "tp0548 SS14-LAMP"
)
observed_names <- names(raw)
observed_names <- gsub("\u00b0", "", observed_names, fixed = TRUE)
observed_names <- gsub("Grillov\u00e1", "Grillova", observed_names, fixed = TRUE)
stopifnot(
  length(observed_names) == length(expected_names),
  identical(observed_names, expected_names),
  nrow(raw) == 133L
)

pick <- function(i) raw[[i]]

calls <- data.frame(
  sample_no = as.integer(pick(1)),
  sample_type = pick(2),
  strain = pick(4),
  lineage = pick(12),
  tp0705_allele = pick(7),
  tm_506 = as.numeric(pick(13)),
  call_506 = pick(14),
  tm_625 = as.numeric(pick(15)),
  call_625 = pick(16),
  tm_708 = as.numeric(pick(17)),
  call_708 = pick(18),
  tp0548_nichols_lamp = pick(19),
  tp0548_ss14_lamp = pick(20),
  stringsAsFactors = FALSE
)

# ---- the measured totals this file must keep -------------------------------
stopifnot(
  nrow(calls) == 133L,
  !anyDuplicated(calls$sample_no),
  # sample composition
  sum(calls$sample_type == "Swab") == 105L,
  sum(calls$sample_type == "Saliva") == 28L,
  sum(calls$lineage == "SS14-East Asia") == 110L,
  sum(calls$lineage == "Nichols") == 9L,
  sum(calls$lineage == "SS14\u03a9-A") == 7L,
  sum(calls$lineage == "-") == 5L,
  sum(calls$lineage == "SS14\u03a9-B 1A") == 2L,
  # the tp0705 allele each call is supposed to reflect
  sum(calls$tp0705_allele == "8") == 110L,
  sum(calls$tp0705_allele == "3") == 11L,
  sum(calls$tp0705_allele == "1") == 7L,
  sum(calls$tp0705_allele == "9") == 5L,
  # 506: V/T/A partition, and each allele letter maps to one call
  sum(calls$call_506 == "V") == 106L,
  sum(calls$call_506 == "T") == 11L,
  sum(calls$call_506 == "A") == 12L,
  sum(calls$call_506 == "NA") == 4L,
  all(calls$call_506[calls$tp0705_allele == "8" & calls$call_506 != "NA"] == "V"),
  all(calls$call_506[calls$tp0705_allele == "3"] == "T"),
  all(calls$call_506[calls$tp0705_allele %in% c("1", "9") & calls$call_506 != "NA"] == "A"),
  # 625 and 708 keep their own partitions
  sum(calls$call_625 == "M") == 123L,
  sum(calls$call_625 == "V") == 7L,
  sum(calls$call_625 == "NA") == 3L,
  sum(calls$call_708 == "G") == 117L,
  sum(calls$call_708 == "S") == 12L,
  sum(calls$call_708 == "NA") == 4L,
  # an NA call still has a number in its Tm cell
  all(!is.na(calls$tm_506)),
  all(!is.na(calls$tm_625)),
  all(!is.na(calls$tm_708)),
  # the calls the paper could not resolve, by sample number
  identical(calls$sample_no[calls$call_506 == "NA"], c(112L, 113L, 123L, 128L)),
  identical(calls$sample_no[calls$call_625 == "NA"], c(112L, 123L, 128L)),
  identical(calls$sample_no[calls$call_708 == "NA"], c(112L, 113L, 123L, 128L)),
  # tp0548 dipstick: 122 untested, 9 Nichols-positive, 2 SS14-positive
  sum(calls$tp0548_nichols_lamp == "-" & calls$tp0548_ss14_lamp == "-") == 122L,
  sum(calls$tp0548_nichols_lamp == "Positive" & calls$tp0548_ss14_lamp == "Negative") == 9L,
  sum(calls$tp0548_nichols_lamp == "Negative" & calls$tp0548_ss14_lamp == "Positive") == 2L
)

# the within-probe separation the readout model would have to reproduce
rng <- function(tm, call, letter) {
  x <- tm[call == letter]
  c(min = min(x), max = max(x))
}
stopifnot(
  all(rng(calls$tm_506, calls$call_506, "V") == c(60.08, 61.52)),
  all(rng(calls$tm_506, calls$call_506, "T") == c(63.64, 64.03)),
  all(rng(calls$tm_506, calls$call_506, "A") == c(71.52, 71.91)),
  all(rng(calls$tm_625, calls$call_625, "M") == c(63.09, 63.86)),
  all(rng(calls$tm_625, calls$call_625, "V") == c(69.96, 70.94)),
  all(rng(calls$tm_708, calls$call_708, "G") == c(72.94, 74.74)),
  all(rng(calls$tm_708, calls$call_708, "S") == c(67.38, 68.00))
)

readr::write_csv(calls, "inst/extdata/lamp_tp0705_melting_calls.csv", na = "NA")

# ---- the reaction conditions, one row per Methods value --------------------
doi <- "10.1128/spectrum.01565-25"
section <- "Methods / Melting curve analysis"
weak <- "reaction condition, not a readout outcome; no allele gap is attached to these values"
condition <- data.frame(
  item = c(
    "reaction volume", "FIP/BIP", "LF/LB", "F3/B3", "isothermal buffer",
    "MgSO4", "dNTPs", "DNA input", "betaine", "Bst 3.0 polymerase",
    "QProbe", "EvaGreen", "amplification temperature", "amplification time (tp0705)",
    "amplification time (tp0548)", "melt start", "melt hold", "melt ramp",
    "melt range", "Tm decision margin", "detection limit"
  ),
  value = c(
    "25", "1.6", "0.8", "0.2", "2.5", "6", "200", "5", "0.2", "0.32",
    "80", "1.25", "68", "80", "20", "98", "40", "0.2", "40-95", "+/-1", "20"
  ),
  unit = c(
    "uL", "uM each", "uM each", "uM each", "uL of 10x Isothermal Amplification Buffer II",
    "mM", "uM each dNTP", "uL DNA sample", "M", "U/uL", "nM", "uL EvaGreen Dye",
    "C", "min", "min", "C for 3 min", "C for 1 min", "C/s", "C", "C", "copies of plasmid per reaction"
  ),
  source_section = section,
  quote = c(
    "The final reaction volume was 25 uL",
    "1.6 uM FIP and BIP primers",
    "0.8 uM LF and LB primers",
    "0.2 uM F3 and B3 primers",
    "2.5 uL of 10x Isothermal Amplification Buffer II",
    "6 mM MgSO4",
    "200 uM each dNTP",
    "5 uL DNA sample",
    "0.2 M betaine",
    "0.32 U/uL Bst 3.0 polymerase",
    "80 nM QProbe (for tp0705)",
    "1.25 uL EvaGreen Dye (for tp0548)",
    "The LAMP reactions for tp0705 and tp0548 were carried out at 68 C",
    "for 80 and 20 min, respectively",
    "for 80 and 20 min, respectively",
    "conducted at 98 C for 3 min",
    "40 C for 1 min",
    "from 40 C to 95 C at a ramp rate of 0.2 C/s",
    "from 40 C to 95 C at a ramp rate of 0.2 C/s",
    "the Tm margin for each variant was set to +/-1 C",
    "The detection limit of the LAMP assays was 20 copies of plasmid per reaction"
  ),
  outcome_definition = weak,
  evidence_tier = "weak",
  source_file = "spectrum.01565-25 Methods / Melting curve analysis",
  doi = doi,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(condition) == 21L,
  all(condition$evidence_tier == "weak"),
  all(nzchar(condition$quote))
)

# spot-check the two values R/marker_core.R names without a number
stopifnot(
  condition$value[condition$item == "betaine"] == "0.2",
  condition$value[condition$item == "MgSO4"] == "6"
)

utils::write.csv(condition, "inst/extdata/lamp_tp0705_reaction_conditions.csv", row.names = FALSE)

cat(
  "wrote 133 melting calls and", nrow(condition), "reaction-condition rows\n",
  "506: V/T/A =", sum(calls$call_506 == "V"), "/", sum(calls$call_506 == "T"),
  "/", sum(calls$call_506 == "A"), " NA =", sum(calls$call_506 == "NA"), "\n",
  "625: M/V =", sum(calls$call_625 == "M"), "/", sum(calls$call_625 == "V"),
  " NA =", sum(calls$call_625 == "NA"), "\n",
  "708: G/S =", sum(calls$call_708 == "G"), "/", sum(calls$call_708 == "S"),
  " NA =", sum(calls$call_708 == "NA"), "\n"
)
