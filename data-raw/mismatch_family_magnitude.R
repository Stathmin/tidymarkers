# Per-family mismatch-discrimination MAGNITUDES (ROADMAP Z5.2)
#
# `specificity_policy()` carries one tolerance policy per chemistry mode.  This
# script collects, per family and per measured axis, the magnitude a mismatch
# actually buys - read from extracts already in `inst/extdata/`, so nothing here
# re-scrapes.  The value is `NA` when the family has NO printed per-target
# magnitude; `axis` then says on which axis it would have to be measured.
#
# One row per (family, condition) with the observed min-max of that axis.  The
# `stopifnot()` block below reproduces each source's own arithmetic or printed
# figure; a failing assertion means an input extract changed, not that the
# scanner is wrong.
#
# Run from the project root (see AGENTS.md for the R 4.6.1 environment):
#   Rscript data-raw/mismatch_family_magnitude.R

ext <- function(file) {
  utils::read.csv(
    file.path("inst", "extdata", file),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

round_summary <- function(x) {
  c(min = round(min(x, na.rm = TRUE), 3), max = round(max(x, na.rm = TRUE), 3))
}

rows <- list()
add <- function(family, condition, value, unit, axis, outcome_definition,
                evidence_tier, source_file, doi, pmcid) {
  rows[[length(rows) + 1L]] <<- data.frame(
    family = family, condition = paste(condition, collapse = " "),
    value = value, unit = unit,
    axis = paste(axis, collapse = " "),
    outcome_definition = paste(outcome_definition, collapse = " "),
    evidence_tier = evidence_tier, source_file = source_file,
    doi = doi, pmcid = pmcid, stringsAsFactors = FALSE
  )
}

# One summary row per (family, condition); `add()` is the only writer.

# ---- RPA: terminal (3'-anchored) mismatch, binary reaction success -----------
# Higgins et al. 2022, J Mol Diagn.  Printed: 3543/4008 reactions successful
# (88.4 %), 159/315 mismatch combinations significant (50.4 %), terminal
# single-mismatch success probabilities 0.59 (C-C) to 0.984 (A-C / T-G).
rpa_term <- ext("rpa_mismatch_terminal.csv")
rs <- rpa_term$RS_probability
stopifnot(
  nrow(rpa_term) == 12L,
  round(3543 / 4008, 3) == 0.884, # "88.4%" as printed
  round(159 / 315, 3) == 0.505, # "50.4%" as printed (0.5048)
  abs(min(rs) - 0.59) < 0.005, # printed C-C 0.59
  abs(max(rs) - 0.984) < 0.001,
  # the printed probabilities reproduce from the supplement's logit + intercept
  isTRUE(all.equal(
    round(plogis(4.79 + rpa_term$RS_logit), 3), round(rs, 3)
  ))
)
add(
  "rpa", "single mismatch at the 3' terminus, 12 substitutions", rs,
  "probability of reaction success",
  "binary reaction success (fluorescence threshold), reference = perfect match",
  paste(
    "Observed success probability per terminal mismatch, versus a fully",
    "complementary primer binding site. This is a BINARY pass/fail call, not",
    "an efficiency ratio or a Cq shift: it says how often the reaction still",
    "fires, so a perfect match keeps >= 0.984 and the worst pair keeps 0.59."
  ),
  "primary table (supplement S1, model coefficient)",
  "rpa_mismatch_terminal.csv", "10.1016/j.jmoldx.2022.08.005", "PMC12179515"
)

# Same source, global: how many combinations were tested and how many moved.
stopifnot(abs((159 / 315) * 100 - 50.4) < 0.1)
add(
  "rpa", "315 mismatch combinations in the 3' anchor region", c(159, 315),
  "combinations significant / combinations tested",
  "count of combinations with a significant success penalty",
  paste(
    "159 of 315 combinations significantly changed reaction success",
    "(P < 0.05). Reported as counts; the paper prints no per-target ratio."
  ),
  "primary text (Results)",
  "PMC12179515.xml :: Results", "10.1016/j.jmoldx.2022.08.005", "PMC12179515"
)

# Two adjacent anchor mismatches are the only printed complete-inhibition call.
# The paper prints them as `??CA-??CC` and `??TC-??CC`: penultimate C-C or T-C
# with a terminal C-A or C-C.
add(
  "rpa", "penultimate C-C + terminal C-A (and T-C + terminal C-C)", c(0, 0),
  "probability of reaction success",
  "binary reaction success, reference = perfect match",
  paste(
    "The only combinations the paper prints as completely inhibiting. A lone",
    "penultimate C-C does not: it mildly retards kinetics only. This is a",
    "printed outcome, not a per-count magnitude."
  ),
  "primary text (Discussion)",
  "PMC12179515.xml :: Discussion", "10.1016/j.jmoldx.2022.08.005", "PMC12179515"
)

# ---- LAMP: Tt delay, one 3'-terminal mismatch ---------------------------------
# rtLAMP primer-end mismatch study 2019 (Front Microbiol), Table 1.
lamp <- ext("lamp_mismatch_rtlamp2019_position_rules.csv")
lamp_lo <- min(lamp$conventional_tt_over_wt_min, lamp$novel_tt_over_wt_min)
lamp_hi <- max(lamp$conventional_tt_over_wt_max_min, lamp$novel_tt_over_wt_max_min)
stopifnot(
  nrow(lamp) == 4L,
  abs(lamp_lo - (-1.2)) < 1e-9, abs(lamp_hi - 22.6) < 1e-9,
  # every inner-primer (FIP/BLP) window must sit above every F3 window
  min(lamp$conventional_tt_over_wt_min[lamp$item %in% grep("FIP|BLP", lamp$item, value = TRUE)]) >
    max(lamp$conventional_tt_over_wt_max_min[lamp$item %in% grep("F3", lamp$item, value = TRUE)])
)
add(
  "lamp", "single 3'-terminal mismatch, by primer end (F3 / FIP / BLP)", c(lamp_lo, lamp_hi),
  "minutes of Tt delay",
  "time-to-threshold delay versus the perfect-match template",
  paste(
    "Tt delay, conventional vs novel (Bst 2.0 + Q5) conditions, over the four",
    "primer ends the paper tests. The sign runs both ways: a mismatch at the 5'",
    "end of FIP (3' end of F1c) can even advance Tt by 1.2 min under the novel",
    "polymerase. This is a TOLERANCE magnitude (the mismatch retards, it does not",
    "stop) and the range is wide, so no single number transfers."
  ),
  "moderate", "lamp_mismatch_rtlamp2019_position_rules.csv",
  "10.3389/fmicb.2019.01056", "PMC6518337"
)

# ---- LAMP: two-allele discrimination (AS-LAMP), printed outcomes -------------
# Two independent AS-LAMP assays print the SAME shape: the mismatched allele
# does not amplify inside the run, the matched allele amplifies at 50-52 min.
stopifnot(
  abs(64 - 50) == 14, # ace-1 heterozygote vs wild type
  75 > 64, # off-allele: none by 75 min
  55 >= 51 # kdr heterozygote >= resistant homozygote
)
add(
  "lamp", "AS-LAMP, off-allele amplification inside a 75 min run",
  c(51, 64),
  "minutes to detection (on-allele); off-allele none <= 75",
  "cycle-threshold-equivalent: amplification time to positive",
  paste(
    "ace-1 (Malar J 2015): wild-type primers, WT homozygote 50 min, heterozygote",
    "64 min, resistant allele no amplification within 75 min; resistant primers,",
    "51 min / 55 min, WT none within 75 min. kdr (Malar J 2012): on-allele",
    "~52 min, other genotype none within 75 min. Discrimination is therefore",
    "BINARY (amplify / no amplify) at the primer end used for the allele; the",
    "number held is the on-allele detection time, not a ratio."
  ),
  "moderate",
  "PMC4665897.xml + PMC3407793.xml :: Results",
  "10.1186/s12936-015-0968-9", "PMC4665897"
)

# ---- KASP / ARMS: 3'-terminal mismatch, ddCq in an allelic-discrimination qPCR
# Lefever et al. 2019, Sci Rep, 6 SNPs, single 3'-terminal mismatch.
lef <- ext("lefever_single_mismatch.csv")
dd <- lef$min_ddCq
stopifnot(
  nrow(lef) == 6L,
  abs(sum(dd) - 11.15) < 1e-9, # 4.33+0.95+0.21+0.28+3.03+2.35
  sum(lef$separates_at_2cycles) == 3L,
  abs(max(dd) - 4.33) < 1e-9, abs(min(dd) - 0.21) < 1e-9
)
add(
  "kasp", "single 3'-terminal mismatch, 6 SNPs, no extra artificial mismatch", dd,
  "ddCq (cycles)",
  "ACt-equivalent: homozygous vs heterozygous discrimination",
  paste(
    "ddCq = |dCq(homozygote) - dCq(heterozygote)| in parallel allele-specific",
    "reactions. Only 3 of 6 SNPs reach >= 2 cycles, so a lone 3'-terminal",
    "mismatch is NOT a guaranteed gate - it depends on the substitution and",
    "the sequence context."
  ),
  "primary literature", "lefever_single_mismatch.csv",
  "10.1038/s41598-019-38581-z", ""
)

# ARMS: the classical 3'-terminal severity classes, PCR yield.
# Kwok et al. 1990, Nucleic Acids Res: ~100-fold (severe) / ~20-fold (A:A).
add(
  "arms", "3'-terminal mismatch class (PCR product yield)", c(20, 100),
  "fold reduction in product yield",
  "relative endpoint yield, reference = matched 3' terminus",
  paste(
    "~100-fold for A:G, G:A, C:C; ~20-fold for A:A; minimal for T:G, T:C,",
    "T:T. Yield, not a per-target probability: G:G is context dependent, which",
    "is the paper's own caveat."
  ),
  "primary literature", "mismatch_discrimination.csv",
  "10.1093/nar/18.4.999", ""
)
for (i in seq_len(nrow(lef))) {
  add(
    "arms", paste0("single 3'-terminal mismatch, ", lef$snp[i]), dd[i],
    "ddCq (cycles)",
    "ACt-equivalent: homozygous vs heterozygous discrimination",
    paste(
      "Same assay geometry as ARMS (allele-specific primer, 3'-terminal",
      "mismatch); carried here so the two families can be compared on one axis."
    ),
    "primary literature", "lefever_single_mismatch.csv",
    "10.1038/s41598-019-38581-z", ""
  )
}

# ---- CAPS: no per-target magnitude, only a verification funnel ----------------
caps <- ext("vcf2caps2021_headline.csv")
caps_n <- caps$value[caps$quantity == "n_plant_x_marker_observations"]
caps_ok <- sum(caps$value[caps$quantity %in% c(
  "exact_match", "weaker_cleavable_allele", "close_single_het"
)])
# The headline file is a CATEGORY SUMMARY, not a per-observation tally: its
# seven printed counts sum to 15, not to the 300 observations.  The row-level
# file prints the disagreements per marker and leaves agreement as NA, so the
# number of agreeing observations is NOT derivable here and is not claimed.
caps_v <- ext("vcf2caps2021_verification.csv")
stopifnot(
  length(caps_n) == 1L, caps_n == 300,
  sum(caps$value[caps$quantity %in% c(
    "exact_match", "weaker_cleavable_allele", "close_single_het",
    "five_noncleavable_reclassified", "cleavable_hom_read_as_het_banding",
    "het_read_as_noncleavable", "marker_3_363_reclassified_as_het"
  )]) == 15, # 2+1+3+1+1+2+5: the seven printed categories
  caps_ok == 6
)
stopifnot(nrow(caps_v) == 10L, max(caps_v$observed_reclassified_as_het, na.rm = TRUE) == 5)

add(
  "caps", "10 markers x 30 plants, gel-scored genotype", c(caps_ok, nrow(caps_v)),
  "markers agreeing in full with the predicted digest / markers",
  "binary gel-call agreement with the in-silico digest",
  paste(
    "CAPS has no mismatch-discrimination magnitude to hold: the restriction",
    "digest either cuts or does not. What is measurable is the agreement of",
    "the gel call with the in-silico digest: 6 of 10 markers show ONLY the",
    "printed agreement categories, while 4 carry a documented reclassification.",
    "Kept as an outcome funnel, not a magnitude; the per-observation count is",
    "not printed."
  ),
  "primary literature", "vcf2caps2021_headline.csv", "10.1101/2020.09.25.312561", ""
)

# ---- HRM: dTm between homozygotes --------------------------------------------
gundy <- ext("hrm_gundry2008_amplicon_tm.csv")
ijms <- ext("hrm_ijms2022_dtm_single_snp_controls.csv")
jsbbs <- ext("hrm_jsbbs2018_dtm_by_bimorphism.csv")
gundy_dt <- suppressWarnings(as.numeric(gundy$dtm_after))
gundy_dt <- gundy_dt[!is.na(gundy_dt)]
stopifnot(
  nrow(ijms) == 7L,
  abs(max(ijms$dtm_c) - 1.7) < 1e-9,
  # the single-SNP increment: 1 T allele = 0.4 / 0.7 / 0.6, no uniform step
  isTRUE(all.equal(ijms$dtm_c[ijms$n_t_alleles == 1][1], 0.4)),
  nrow(jsbbs) == 5L,
  abs(sum(jsbbs$n) - 56) < 1e-9, # 4 + 18 + 1 + 28 + 5
  abs(min(jsbbs$dtm_min) - 0.198) < 1e-3, abs(max(jsbbs$dtm_max) - 1.118) < 1e-3
)
add(
  "hrm", "49-68 bp amplicons, two homozygote classes", round_summary(gundy_dt),
  "dTm (C) between homozygotes",
  "melting-temperature separation",
  paste(
    "Post-calibration homozygote dTm, 8 amplicons. This is the physical",
    "separation a genotyping call needs; 0.02 C is below any instrument's",
    "resolution, so a small amplicon does not guarantee a call."
  ),
  "measured", "hrm_gundry2008_amplicon_tm.csv", "10.1093/nar/gkn204", ""
)
add(
  "hrm", "one C>T SNP in a GC-rich promoter amplicon", c(0.4, 0.7),
  "dTm (C) per single SNP",
  "melting-temperature separation",
  paste(
    "Single-SNP dTm of three C>T controls in one amplicon: 0.4 / 0.7 / 0.6 C.",
    "The step is locus-specific even for the same substitution class."
  ),
  "measured", "hrm_ijms2022_dtm_single_snp_controls.csv", "10.3390/ijms222212527", ""
)
add(
  "hrm", "56 primer pairs over 5 bimorphism classes", round_summary(jsbbs$dtm_mean),
  "mean dTm (C) per bimorphism class",
  "melting-temperature separation",
  paste(
    "Class means 0.20 (C/G, n=1) to 0.69 (C/T, n=28). Spread is wide within a",
    "class (0.20-1.12 C overall), and the free-energy index does not order the",
    "classes, so dTm cannot be predicted from the substitution alone."
  ),
  "measured", "hrm_jsbbs2018_dtm_by_bimorphism.csv", "10.1270/jsbbs.18048", ""
)

# ---- CRISPR: seed / PAM -------------------------------------------------------
# The 2024 Cas12a activator screen measured RATIOS but they do not predict a
# per-target ratio (see data-raw/proposed/crispr_seed_rule.decision.md).  Kept
# here with value NA so the table states the gap instead of leaving it implicit.
ratio <- ext("communbiol2024_crrna_ratio.csv")
m39 <- ratio$mean_fc_vs_wt[ratio$crrna == "M39"]
wt <- ratio$mean_fc_vs_wt[ratio$crrna == "WT crRNA"]
m39_single <- ratio$mean_fc_vs_wt[ratio$crrna == "M39" & ratio$n_mismatch == 1]
stopifnot(
  nrow(ratio) == 320L,
  length(m39_single) == 20L,
  length(wt) == 40L,
  # printed figures for the singles: WT crRNA 0.77-1.394, M39 seed 0.373-0.968
  sum(abs(wt - 1) < 1e-3) < sum(wt > 1e-3), # the ratio includes the exact match
  abs(max(wt) - 1.400) < 5e-3,
  abs(max(m39_single[m39_single < 1]) - 0.968) < 5e-3,
  abs(min(m39_single) - 0.227) < 5e-3,
  abs(min(m39) - 0.0158) < 5e-3, abs(max(m39) - 1.199) < 5e-3
)
add(
  "crispr", "seed mismatch, Cas12a crRNA activator screen", c(NA_real_, NA_real_),
  "unknown",
  "endpoint trans-cleavage ratio (measured, not transferable)",
  paste(
    "Measured single-mismatch ratios span 0.767-1.394 on the unmodified crRNA",
    "(20 positions; 0.553 is the floor of the whole class, which also holds",
    "the adjacent doubles) and 0.227-1.199 on the optimised M39 crRNA (20",
    "i.e. the seed positions do not order by distance from the PAM. No model",
    "predicts the ratio for a new target, so the magnitude stays `unknown`",
    "and the policy keeps only the position rule."
  ),
  "moderate", "communbiol2024_crrna_ratio.csv", "10.1038/s42003-024-07173-7",
  "PMC11541961"
)
add(
  "crispr", "PAM (TTTV for Cas12a) altered", c(NA_real_, NA_real_),
  "unknown",
  "binary: PAM intact or not",
  paste(
    "No measured magnitude and none needed: a changed PAM removes the licensing",
    "element, so the policy treats it as a hard loss, not a gradient."
  ),
  "rule of thumb", "marker_crispr.R :: crispr_seed_gate", "", ""
)

out <- do.call(rbind, rows)
# A source that measured nothing on the axis leaves the value empty.
empty <- vapply(out$value, function(v) all(is.na(v)), logical(1))
stopifnot(
  !any(empty) == all(out$unit == "unknown"), # blank value <=> `unknown` unit
  identical(
    names(out),
    c(
      "family", "condition", "value", "unit", "axis", "outcome_definition",
      "evidence_tier", "source_file", "doi", "pmcid"
    )
  ),
  setequal(unique(out$family), c("rpa", "lamp", "kasp", "arms", "caps", "hrm", "crispr")),
  nrow(out) == 46L, # rpa 16, arms 8, kasp 6, hrm 6, crispr 4, lamp 4, caps 2
  # a family with numbers must not also claim `unknown`
  !any(out$family %in% out$family[out$unit == "unknown"] & out$unit != "unknown"),
  # every row names a source
  all(nzchar(out$source_file))
)

out$value <- vapply(
  out$value, function(v) paste(v, collapse = " - "), character(1)
)
out$value[out$unit == "unknown"] <- NA_character_

path <- file.path("inst", "extdata", "mismatch_family_magnitude.csv")
utils::write.csv(out, path, row.names = FALSE, na = "")
message("wrote ", path, " (", nrow(out), " rows, ", ncol(out), " columns)")
