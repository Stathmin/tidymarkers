# =============================================================================
#  data-raw/kutyavin2000.R  --  Z9.3: measured Tm / dTm of 3'-MGB probes
# =============================================================================
#  Kutyavin et al. 2000 (Nucleic Acids Res 28:655-661, 10.1093/nar/28.2.655,
#  PMC102528) measured UV melting temperatures for the SAME probe sequences with
#  and without a 3'-minor groove binder, then ran the same probes in the
#  5'-nuclease (TaqMan) assay.  This is the source Z9.3 was missing: it reports
#  an MGB Tm gain and a match-minus-mismatch dTm, not a design rule.
#
#  Every number below is printed in the paper's Results/Discussion:
#    * ODN1 15mer: 60 C with MGB vs 42 C unmodified (+18 C)
#    * 12mer MGB probe (66 C) vs no-MGB 27mer (65 C)
#    * fluorogenic MGB 12-18mers: 66-70 C; same-sequence no-MGB: 44-56 C
#    * 12mer MGB match-minus-mismatch dTm: 20 C; no-MGB 4-6 C
#    * probe #7 T/G mismatch (7th base from 3' end): dTm 6.5 -> 11 C
#    * probe #2 same mismatch 11th from 3' end: dTm 8.5 C, no MGB effect
#    * MGB probe length range 12-20 nt; 13-18mers Tm near 65 C at GC 17-69 %
#
#  What it does NOT report: dCt or endpoint ratio between matched and
#  mismatched probes, and no LNA arm.  Those gaps are named in the rows.
#
#  Run:  Rscript data-raw/kutyavin2000.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

od <- paste(
  "UV melting of the same probe sequence with and without a 3'-MGB (DPI3/CDPI3);",
  "the 20 C and 4-6 C values are match-minus-mismatch dTm, not Tm; no dCt and",
  "no endpoint ratio is reported, and no LNA probe was tested"
)

rows <- data.frame(
  item = c(
    "ODN1 15mer, unmodified, Tm (match)",
    "ODN1 15mer, 3'-MGB, Tm (match)",
    "ODN1 15mer, Tm gain from MGB",
    "12mer MGB probe, Tm (match)",
    "no-MGB 27mer probe, Tm (match)",
    "fluorogenic MGB probes 12-18 nt, lowest Tm (match)",
    "fluorogenic MGB probes 12-18 nt, highest Tm (match)",
    "same-sequence no-MGB probes 12-18 nt, lowest Tm (match)",
    "same-sequence no-MGB probes 12-18 nt, highest Tm (match)",
    "12mer MGB probe, match-minus-mismatch dTm",
    "no-MGB probes, match-minus-mismatch dTm, low end",
    "no-MGB probes, match-minus-mismatch dTm, high end",
    "probe #7 T/G mismatch at 7th base from 3' end, no-MGB dTm",
    "probe #7 T/G mismatch at 7th base from 3' end, MGB dTm",
    "probe #2 T/G mismatch at 11th base from 3' end, dTm (MGB and no-MGB)",
    "MGB probe length, lower bound",
    "MGB probe length, upper bound"
  ),
  value = c(
    42, 60, 18, 66, 65, 66, 70, 44, 56, 20, 4, 6, 6.5, 11, 8.5, 12, 20
  ),
  unit = c(
    rep("degC", 9), rep("degC (match-minus-mismatch)", 6), rep("nt", 2)
  ),
  outcome_definition = od,
  evidence_tier = "moderate",
  source_file = "PMC102528 HTML full text :: Results and Discussion",
  doi = "10.1093/nar/28.2.655",
  pmcid = "PMC102528",
  stringsAsFactors = FALSE
)

## the paper's own arithmetic, not ours
v <- stats::setNames(rows$value, rows$item)
stopifnot(
  unname(v[["ODN1 15mer, 3'-MGB, Tm (match)"]] -
    v[["ODN1 15mer, unmodified, Tm (match)"]]) ==
    unname(v[["ODN1 15mer, Tm gain from MGB"]]),
  unname(v[["ODN1 15mer, Tm gain from MGB"]]) == 18,
  abs(unname(v[["12mer MGB probe, Tm (match)"]]) -
    unname(v[["no-MGB 27mer probe, Tm (match)"]])) <= 1,
  unname(v[["12mer MGB probe, match-minus-mismatch dTm"]]) == 20,
  unname(v[["no-MGB probes, match-minus-mismatch dTm, low end"]]) == 4,
  unname(v[["no-MGB probes, match-minus-mismatch dTm, high end"]]) == 6,
  unname(v[["probe #7 T/G mismatch at 7th base from 3' end, MGB dTm"]]) ==
    unname(v[["probe #7 T/G mismatch at 7th base from 3' end, no-MGB dTm"]]) + 4.5,
  unname(v[["MGB probe length, lower bound"]]) == 12,
  unname(v[["MGB probe length, upper bound"]]) == 20,
  # the short MGB probes sit below the package's 18-25 nt no-MGB window
  unname(v[["MGB probe length, lower bound"]]) < 18
)

out <- file.path(root, "inst", "extdata", "taqman_kutyavin2000_mgb_tm.csv")
utils::write.csv(rows, out, row.names = FALSE)
cat("kutyavin2000: 15mer MGB 60 C vs 42 C (+18); 12mer MGB dTm 20 C vs 4-6 C\n")
