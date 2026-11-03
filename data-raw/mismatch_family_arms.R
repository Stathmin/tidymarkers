# =============================================================================
#  data-raw/mismatch_family_arms.R  --  Z5.2: ARMS / KASP 3'-anchor MAGNITUDE
# =============================================================================
#  `specificity_policy()` applies the 3'-terminal penalty the same way to every
#  scheme, so Z5.2 needs the size of that penalty for the ARMS-style assays.
#  Lefever et al. 2019 (Sci Rep, PMC6377641) print it directly, and the
#  package's own extract (lefever_single_mismatch.csv) already holds the
#  per-SNP ddCq.  What was missing is the DUAL boundary the paper states in
#  prose:
#    * a single 3'-terminal mismatch "is not sufficient for reliable genotype
#      calling";
#    * with a second mismatch at position 3 (fourth base from the 3' end) the
#      correct genotype was called in 376/384 (98 %) of sample x concentration
#      combinations at a 5-dCq cut-off;
#    * position 2 and 3 give the highest average ddCq, position 3 wins once the
#      penalty on the match allele is counted, so the optimum is -3.
#  This generator adds those rows to the AXIS summary rather than repeating the
#  six per-SNP ddCq values (they stay in lefever_single_mismatch.csv).
#
#  Run:  Rscript data-raw/mismatch_family_arms.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

correct <- 376L; combos <- 384L; cut <- 5
stopifnot(
  abs(correct / combos - 0.979) < 0.005,   # printed "98%"
  cut == 5L,
  # the paper's conservative calling threshold is below the cut-off it used
  2 < cut
)

# The package extract must still show that one terminal mismatch alone is weak:
# 3 of 6 SNPs separate at a 2-cycle cut-off.
lef <- read.csv("inst/extdata/lefever_single_mismatch.csv")
stopifnot(
  nrow(lef) == 6L,
  sum(lef$separates_at_2cycles) == 3L,
  abs(max(lef$min_ddCq) - 4.33) < 0.005,   # printed ddCq for rs6734469
  abs(min(lef$min_ddCq) - 0.21) < 0.005    # printed ddCq for rs187115 (< 1 cycle)
)

rows <- data.frame(
  family = rep("arms", 3L),
  n_mismatches = c(1L, 2L, 2L),
  mismatch_position = c(
    "3' terminal base only (position 0)",
    "3' terminal base plus an artificial mismatch at position -3",
    "artificial mismatch position, averaged over positions -1 to -5"
  ),
  metric = c("genotype_call_rate_pct", "genotype_call_rate_pct", "delta_cq_cycles"),
  value = c(NA, 100 * correct / combos, NA),
  outcome_definition = c(
    "Not a magnitude: the paper states a single 3'-terminal mismatch 'is not sufficient for reliable genotype calling' and only rs4245739 and rs1027154 separated the three genotypes at all, so no call rate is printed and the value stays empty",
    "Share of 384 sample x input-concentration combinations whose genotype was called correctly with a 5-dCq cut-off in the double-mismatch (DMAS) design, 376/384 = 98 %; down to 0.25 ng (~80 copies) input",
    "Positions -2 and -3 give the highest average match-minus-mismatch dCq; position -3 is chosen because its standard deviation is low and it barely costs the match allele. The paper prints the ordering, not a per-position cycle value here (the per-SNP ddCq are in lefever_single_mismatch.csv), so the value stays empty"
  ),
  evidence_tier = c("primary text (Results)", "primary text (Results + figure 3)",
                    "primary text (figure 2 A/B, written in prose)"),
  source_file = c("PMC6377641 :: Results", "PMC6377641 :: Results",
                  "PMC6377641 :: Results + lefever_single_mismatch.csv"),
  doi = rep("10.1038/s41598-019-38581-z", 3L),
  pmcid = rep("PMC6377641", 3L),
  stringsAsFactors = FALSE
)

stopifnot(
  identical(names(rows), c("family", "n_mismatches", "mismatch_position", "metric",
                           "value", "outcome_definition", "evidence_tier",
                           "source_file", "doi", "pmcid")),
  nrow(rows) == 3L,
  all(rows$family == "arms"),
  !is.na(rows$value[2]),
  # the two rows without a printed number are left empty on purpose
  is.na(rows$value[1]), is.na(rows$value[3])
)

write.csv(rows, "inst/extdata/mismatch_family_arms.csv", row.names = FALSE)
message("wrote inst/extdata/mismatch_family_arms.csv (", nrow(rows), " rows)")
