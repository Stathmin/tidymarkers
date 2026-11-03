# =============================================================================
#  data-raw/mismatch_family_crispr.R  --  Z5.2: CRISPR mismatch MAGNITUDE
# =============================================================================
#  Z5.2 is blocked because `specificity_policy()` uses ONE mismatch policy for
#  RPA, CRISPR and KASP.  The existing CRISPR extract
#  (communbiol2024_crrna_ratio.csv) gives per-position crRNA ratios that do not
#  predict the allele, so `crispr_seed_rule` reports `unknown`.  This generator
#  adds the CRISPR magnitudes that DO predict something measurable:
#
#   A. Spacer LENGTH and target concentration gate whether a single mismatch
#      matters at all (Nucleic Acids Res 2024, PMC10810210, LbuCas13a):
#      * 28-nt spacer: "no single mismatch discriminatory capacity" at the
#        concentrations tested; with a 20-nt spacer a mismatch at +7/+8/+9 gives
#        "up to a 70 % reduction in activity" versus perfect match;
#      * mismatch at position 7 costs ~40 % of signal-to-noise in the 20-nt
#        spacer; positions 4 and 11 stay within 20 % of the perfect match;
#      * a 16-nt spacer needs a perfect match for RNase activity.
#
#   B. Mismatch NUMBER in an RfxCas13d spacer (Sci Adv 2024, PMC11654686):
#      3 mismatches with the WT K-RAS transcript exceed the silencing limit
#      (11.0 % residual WT silencing), 2 with the SNV still silence it
#      (54.8 %).  This is in-cell knockdown, not a diagnostic fluorescence.
#
#   C. Cas12a allele CONTRAST (fold between SNP and WT fluorescence):
#      * one designed extra mismatch near the SNP moves the signal ratio from
#        "only twice as strong" to "nearly ten times" (Chem Sci 2021,
#        PMC8179484);
#      * one-pot Cas12a, canonical TTT PAM: 2.8-6.0-fold SNP/WT contrast;
#        suboptimal ATT PAM: 5.9-18.4-fold (Nat Commun 2026, PMC13473493).
#
#   D. Occluding the Cas13 crRNA improves the WT/mismatch discrimination ratio
#      ~50-fold and makes 24/24 single mismatches significant, against 13/24
#      without the occluder; without the occluder only 74/96 tested mismatches
#      lower activity at all (Nat Biotechnol 2026, PMC12584110).  A ratio
#      improvement, not an absolute fold-signal.
#
#  Not rows (no magnitude printed): the Dual-CRISPR "80 %-85 % fluorescence
#  reduction" (NAR 2026, PMC12956341) is approximate and paired with no table;
#  the PAM-independent Cas12a paper (JACS Au 2025, PMC11938009) reports 0.1 %
#  allele fraction but no per-mismatch number.  Both are named in
#  data-raw/proposed/z5.2-family-magnitudes.md.
#
#  Run:  Rscript data-raw/mismatch_family_crispr.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

# --- the papers' own arithmetic, checked before any row is built ------------
vm_70 <- 70; vm_mid <- 40; vm_edge <- 20
stopifnot(vm_70 > vm_mid, vm_mid > vm_edge, vm_70 <= 100)

sa_d10 <- 54.8   # % silencing of G12D (2 mismatches) with crD10C11-MM
sa_wt <- 11.0    # % residual WT silencing (3 mismatches), same crRNA
stopifnot(
  sa_d10 > sa_wt,                   # selectivity means the SNV arm wins
  sa_d10 > 0, sa_d10 <= 100, sa_wt >= 0, sa_wt <= 100,
  3 - 2 == 1                        # one mismatch past the stated limit
)

cs_matched <- 2; cs_extra <- 10
stopifnot(abs((cs_extra / cs_matched) - 5) < 1e-9, cs_extra > cs_matched)

tt <- c(2.8, 6.0)     # canonical TTT PAM, SNP/WT contrast range
att <- c(5.9, 18.4)   # suboptimal ATT PAM
stopifnot(
  all(diff(tt) > 0), all(diff(att) > 0),
  min(att) > min(tt)              # the ATT PAM floor beats the TTT floor,
  # but the two ranges overlap (5.9-6.0), which the rows keep
)

nb_tested <- 24L; nb_sig_noc <- 13L; nb_sig_occ <- 24L
stopifnot(abs(nb_sig_noc / nb_tested - 0.54) < 0.005, nb_sig_occ == nb_tested)

nb_any_reduction <- 74L; nb_loci <- 96L
stopifnot(
  abs(nb_any_reduction / nb_loci - 0.77) < 0.005,   # printed "77%"
  nb_loci == 3L * 4L * 2L * 4L,                     # 3 targets x 4 positions x 2 mutations x 4 replicates
  nb_any_reduction < nb_loci                        # not every mismatch even lowers activity
)

rows <- data.frame(
  family = rep("crispr", 11L),
  n_mismatches = c(1, 1, 1, 1, 1, 3, 2, 2, 1, 1, 1),
  mismatch_position = c(
    "single mismatch, 20-nt spacer (middle: +7 to +9)",
    "single mismatch at spacer position +7, 20-nt spacer",
    "single mismatch at positions +4 and +11, 20-nt spacer",
    "single mismatch anywhere, 28-nt spacer",
    "single mismatch anywhere, 16-nt spacer",
    "central region of a 23-nt spacer, with the WT K-RAS transcript",
    "central region of a 23-nt spacer, with the G12D transcript",
    "one designed A>U mismatch next to the SNP site (plus the 1 mismatch at the SNP)",
    "SNP versus WT, canonical TTT PAM (no engineered mismatch)",
    "SNP versus WT, suboptimal ATT PAM (no engineered mismatch)",
    "single mismatches at 4 spacer positions x 2 substitutions x 3 targets, no occluder"
  ),
  metric = c(
    "activity_reduction_pct", "snr_reduction_pct", "snr_reduction_pct",
    "activity_reduction_pct", "activity_reduction_pct",
    "residual_silencing_pct", "silencing_pct", "fold_signal",
    "fold_allele_contrast", "fold_allele_contrast", "activity_reduction_pct"
  ),
  value = c(vm_70, vm_mid, vm_edge, 0, 100, sa_wt, sa_d10, cs_extra,
            max(tt), max(att), 100 * nb_any_reduction / nb_loci),
  value_low = c(NA, NA, NA, NA, NA, NA, NA, cs_matched, min(tt), min(att), NA),
  unit = c(rep("percent of perfect-match activity", 5L),
           "percent of transcript removed", "percent of transcript removed",
           "fold fluorescence ratio", "fold SNP/WT fluorescence",
           "fold SNP/WT fluorescence", "percent of mismatches"),
  outcome_definition = c(
    "Upper bound of the loss of LbuCas13a RNase activity for a single mismatch in the middle of a 20-nt spacer ('up to a 70 % reduction in activity, compared to a perfect-match' system); edge positions 1-4 and 17-20 stay active",
    "Loss of signal-to-noise for a single mismatch at spacer position +7, 20-nt spacer; positions +4 and +11 stay within 20 % of the perfect match",
    "Signal-to-noise still within 20 % of the perfect-match system for mismatches at positions +4 and +11 of a 20-nt spacer",
    "A 28-nt LbuCas13a spacer yields 'no single mismatch discriminatory capacity at the target-RNA concentrations tested', so the reduction is recorded as 0, not as a measurement",
    "A 16-nt LbuCas13a spacer 'still required a perfect match for RNase activity', the strongest single-mismatch penalty the paper reports (100 = silencing abolished)",
    "Consequence of 3 mismatches with the WT K-RAS transcript: 'exceeding the upper limit of mismatches that allow for efficient silencing'",
    "Consequence of 2 mismatches with the G12D transcript: the same crRNA still silences the mutant, 54.8 % (second guide 39.9 %)",
    "Cas12a fluorescence with a designed A>U mismatch next to the SNP site, divided by the fully complementary guide: 'nearly ten times'",
    "Fluorescence ratio between SNP and WT alleles with unmodified canonical PAM; the paper prints the range 2.8-6.0-fold, value/value_low keep both ends",
    "Fluorescence ratio between SNP and WT alleles with a suboptimal PAM; range 5.9-18.4-fold, which is why PAM tuning is called 'tunable genotype separation'. The two ranges overlap between 5.9 and 6.0, so a tuned PAM does not beat every canonical PAM",
    "Share of 96 mismatch experiments (3 targets x 4 positions x 2 substitutions x 4 replicates) in which a single mismatch lowered Cas13 activity at all without an occluder, 74/96 (77 %): a quarter of single mismatches are simply inert. Discrimination was significant in only 54 % (13/24) without the occluder and 24/24 with it"
  ),
  evidence_tier = c(
    rep("primary text (prose + figure)", 5L),
    rep("primary text (in-cell silencing, 23-nt spacer)", 2L),
    "primary text (figure axis, written in prose)",
    rep("primary text (figure range, written in prose)", 2L),
    "primary text (counts, Results)"
  ),
  source_file = c(
    rep("PMC10810210", 5L), rep("PMC11654686", 2L), "PMC8179484",
    "PMC13473493", "PMC13473493", "PMC12584110"
  ),
  doi = c(
    rep("10.1093/nar/gkad1132", 5L), rep("10.1126/sciadv.adl0731", 2L),
    "10.1039/d0sc05717g", rep("10.1038/s41467-026-75358-1", 2L),
    "10.1038/s41587-025-02868-6"
  ),
  pmcid = c(
    rep("PMC10810210", 5L), rep("PMC11654686", 2L), "PMC8179484",
    "PMC13473493", "PMC13473493", "PMC12584110"
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  identical(names(rows), c("family", "n_mismatches", "mismatch_position", "metric",
                           "value", "value_low", "unit", "outcome_definition",
                           "evidence_tier", "source_file", "doi", "pmcid")),
  nrow(rows) == 11L,
  all(rows$family == "crispr"),
  !anyNA(rows$value),
  all(rows$n_mismatches >= 1L),
  all(is.na(rows$value_low) | rows$value_low <= rows$value)
)

write.csv(rows, "inst/extdata/mismatch_family_crispr.csv", row.names = FALSE)
message("wrote inst/extdata/mismatch_family_crispr.csv (", nrow(rows), " rows)")
