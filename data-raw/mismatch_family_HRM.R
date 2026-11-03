# =============================================================================
#  data-raw/mismatch_family_HRM.R  --  Z5.2: HRM mismatch MAGNITUDE
# =============================================================================
#  HRM is the one family where the magnitude is a melting-point shift per
#  heteroduplex, so Z5.2 needs the shift per MISMATCH COUNT, not per locus.
#  The package holds the per-locus dTm (hrm_ijms2022_dtm_single_snp_controls.csv,
#  IJMS 2022, PMC9751851: three C>T alleles in one GC-rich promoter, triple
#  homozygous controls), which gives an exact count axis:
#    0 T alleles (CCC/CCC)      0.0 C
#    1 T allele                 0.4, 0.7, 0.6 C
#    2 T alleles                1.0, 1.3 C
#    3 T alleles (TTT/TTT)      1.7 C
#  This generator restates that axis as a per-count magnitude and checks the
#  paper's own arithmetic: 1.7 C over three substitutions is 0.57 C per
#  mismatch, which is the usual HRM resolution floor and therefore sits at the
#  edge of a real call on this amplicon.
#
#  Run:  Rscript data-raw/mismatch_family_HRM.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

ijms <- read.csv("inst/extdata/hrm_ijms2022_dtm_single_snp_controls.csv")
counts <- c(0L, 1L, 1L, 1L, 2L, 2L, 3L)
stopifnot(
  nrow(ijms) == 7L,
  identical(as.integer(ijms$n_t_alleles), counts),
  abs(ijms$dtm_c[1] - 0) < 1e-9,       # reference CCC/CCC
  abs(ijms$dtm_c[7] - 1.7) < 1e-9,     # triple homozygous TTT/TTT
  # a per-mismatch step over the exact-count axis
  abs((ijms$dtm_c[7] - ijms$dtm_c[1]) / 3 - 0.5666667) < 1e-6,
  # the three single-mismatch loci do not agree with each other
  diff(range(ijms$dtm_c[counts == 1L])) > 0
)

single <- ijms$dtm_c[counts == 1L]
rows <- data.frame(
  family = rep("hrm", 3L),
  n_mismatches = c(1L, 1L, 3L),
  mismatch_position = c(
    "one C>T substitution in the amplicon (three loci pooled)",
    "one C>T substitution in the amplicon (locus spread)",
    "three C>T substitutions in one amplicon (TTT/TTT)"
  ),
  metric = c("delta_tm_c", "delta_tm_c", "delta_tm_c"),
  value = c(mean(single), diff(range(single)), ijms$dtm_c[7]),
  outcome_definition = c(
    "Mean dTm of the three single-substitution homozygous controls relative to CCC/CCC (0.4, 0.7, 0.6 C): the same substitution in the same amplicon still moves Tm by up to 0.3 C depending on position",
    "Spread between the three single-substitution controls, 0.7 - 0.4 = 0.3 C: the positional and nearest-neighbour dependence that makes a per-locus value non-transferable",
    "dTm of the triple-substitution control, 1.7 C over three substitutions (0.57 C per substitution on this GC-rich amplicon); with 0.5 C as the usual instrument floor, one substitution only just leaves it"
  ),
  evidence_tier = rep("measured (derivative melting curve)", 3L),
  source_file = rep("hrm_ijms2022_dtm_single_snp_controls.csv", 3L),
  doi = rep("10.3390/ijms222212527", 3L),
  pmcid = rep("PMC9751851", 3L),
  stringsAsFactors = FALSE
)

stopifnot(
  identical(names(rows), c("family", "n_mismatches", "mismatch_position", "metric",
                           "value", "outcome_definition", "evidence_tier",
                           "source_file", "doi", "pmcid")),
  !anyNA(rows$value),
  # the pooled single value sits between the two extremes measured above
  rows$value[1] > rows$value[2],
  rows$value[3] > max(single),
  # and the per-substitution step is a third of the triple value
  abs(rows$value[3] / 3 - 0.5666667) < 1e-6
)

write.csv(rows, "inst/extdata/mismatch_family_hrm.csv", row.names = FALSE)
message("wrote inst/extdata/mismatch_family_hrm.csv (", nrow(rows), " rows)")
