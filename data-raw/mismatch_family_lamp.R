# =============================================================================
#  data-raw/mismatch_family_lamp.R  --  Z5.2: LAMP mismatch MAGNITUDE
# =============================================================================
#  The package already holds the LAMP tolerance magnitude by PRIMER END
#  (lamp_mismatch_rtlamp2019_position_rules.csv, Front Microbiol 2019):
#  a 3'-terminal mismatch retards Tt by 1.1-22.6 min depending on the end.
#
#  This generator adds the other LAMP shape: what a mismatch does to the
#  ENDPOINT, after the amplification has happened.
#
#   a) Melting-peak separation between two insertion alleles amplified by the
#      SAME LAMP chemistry (Microbiol Spectr 2025, PMC12671098, T. pallidum
#      tp0548 MLST): SS14-LAMP 92.46 +/- 0.13 C, Nichols-LAMP
#      93.71 +/- 0.21 C.  The paper treats +/-1 C as the variant margin.
#   b) A QProbe over a mismatch melts lower than over the match, so the paper
#      reads the allele from the peak, not from the amplification time
#      (same source, tp0705 MLST alleles).
#
#  Run:  Rscript data-raw/mismatch_family_lamp.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

tm_ss14 <- 92.46; sd_ss14 <- 0.13
tm_nich <- 93.71; sd_nich <- 0.21
margin <- 1.0

stopifnot(
  # the printout is to two decimals, so the difference must reproduce from them
  abs((tm_nich - tm_ss14) - 1.25) < 1e-9,
  # the separation exceeds both standard deviations together
  tm_nich - tm_ss14 > sd_ss14 + sd_nich,
  # the separation clears the paper's own +/- 1 C margin, which is why the
  # two lineage peaks are called separable
  tm_nich - tm_ss14 > margin
)

rows <- data.frame(
  family = rep("lamp", 3L),
  n_mismatches = c(NA, NA, NA),
  mismatch_position = c(
    "LAMP primer-site insertion (9-bp), SS14 vs Nichols allele",
    "LAMP primer-site insertion (9-bp), SS14 vs Nichols allele",
    "QProbe over a mismatched MLST base (tp0705 alleles)"
  ),
  metric = c("melting_tm_c", "delta_melting_tm_c", "tm_shift_direction"),
  value = c(tm_ss14, tm_nich - tm_ss14, -1),
  outcome_definition = c(
    "Melting peak after a 20-min 68 C LAMP reaction with the SS14-specific primer set (92.46 +/- 0.13 C); the Nichols-specific set peaks at 93.71 +/- 0.21 C, so the pair of peaks, not the amplification time, carries the allele call",
    "Separation between the two lineage-specific melting peaks, 93.71 - 92.46 = 1.25 C; the paper sets a +/- 1 C margin for variant assignment, so this separation is just outside it",
    "Direction only: a QProbe whose sequence mismatches the amplicon melts at a LOWER Tm than the match; the paper prints no per-mismatch delta, so -1 records the sign (negative shift) and not a magnitude"
  ),
  evidence_tier = c("primary text (Results, mean +/- SD)", "derived from the two printed means",
                    "primary text (qualitative, figure 3)"),
  source_file = rep("PMC12671098", 3L),
  doi = rep("10.1128/spectrum.01565-25", 3L),
  pmcid = rep("PMC12671098", 3L),
  stringsAsFactors = FALSE
)

stopifnot(
  identical(names(rows), c("family", "n_mismatches", "mismatch_position", "metric",
                           "value", "outcome_definition", "evidence_tier",
                           "source_file", "doi", "pmcid")),
  all(rows$family == "lamp"),
  !anyNA(rows$value),
  # the derived row is the difference of the two printed peaks (93.71 - 92.46)
  abs(rows$value[2] - 1.25) < 1e-9,
  rows$value[2] == (tm_nich - tm_ss14)
)

write.csv(rows, "inst/extdata/mismatch_family_lamp.csv", row.names = FALSE)
message("wrote inst/extdata/mismatch_family_lamp.csv (", nrow(rows), " rows)")
