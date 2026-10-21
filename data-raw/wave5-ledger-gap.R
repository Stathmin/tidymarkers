# =============================================================================
#  data-raw/wave5-ledger-gap.R
# =============================================================================
#  Wave 4 registered 21 sources and landed their extracts, but no ledger row
#  named most of them: `data-raw/merge-proposals.R` printed them as ORPHAN
#  extracts and, until this wave, only as a NOTE. This script seeds those rows
#  ONCE, and data-raw/science_predictions.R carries them from here on.
#
#  Every row is a claim about code that exists, and every claim is checked
#  before it is written: the named test file must exist, must name the extract,
#  and the extract must be on disk. Same checks as test-science-ledger.R, made
#  here so a wrong row cannot be written at all.
#
#  Run:  Rscript data-raw/wave5-ledger-gap.R      (idempotent)
# =============================================================================

a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(a)) dirname(dirname(normalizePath(sub("^--file=", "", a[[1L]])))) else getwd()
setwd(root)
ext_dir <- file.path(root, "inst", "extdata")
led_path <- file.path(ext_dir, "science_ledger.csv")

led <- utils::read.csv(led_path, stringsAsFactors = FALSE, na.strings = c("NA", ""))

# one row per (source, extract, consumer). Written one at a time so counts
# cannot drift out of alignment across columns.
row <- function(source_key, extract, use, test_file, verdict_moved, reason = NA_character_) {
  data.frame(source_key = source_key, extract = extract, use = use,
             test_file = test_file, verdict_moved = verdict_moved,
             reason = reason, stringsAsFactors = FALSE)
}

seed <- do.call(rbind, list(
  # --- LAMP family -----------------------------------------------------------
  row("badolo_kdr2012", "lamp_badolo_kdr2012_time_to_positive.csv", "test",
      "test-lamp-badolo_kdr2012.R",
      "the AS-LAMP readout is a TIME WINDOW with censoring, not a binary call: plasmid Tt ~50 min (WT primers) and ~52 min (kdr-w), non-target not amplified within 75 min",
      NA),
  row("badolo_kdr2012", "lamp_badolo_kdr2012_primers.csv", "config",
      "test-lamp-badolo_kdr2012.R",
      "the two 3'-terminal BIP B2 dinucleotides are design INPUT, so the time window in the outcome table can be read back to the base that produced it",
      "no verdict is asserted from the primer table alone"),
  row("lamp_comp2018", "lamp_comp2018_durations.csv", "test",
      "test-lamp-comp2018.R",
      "LAMP multiplex competition: a shared primer pool moves each target's time to positive, so assays cannot be compared outside their pool",
      NA),
  row("mismatch_rtlamp2019", "lamp_mismatch_rtlamp2019_tt.csv", "test",
      "test-lamp-mismatch_rtlamp2019.R",
      "real-time LAMP discrimination is a time-to-positive RATIO against wild type, not a per-cycle fold change",
      NA),
  row("mismatch_rtlamp2019", "lamp_mismatch_rtlamp2019_position_rules.csv", "test",
      "test-lamp-mismatch_rtlamp2019.R",
      "the position rules the paper states and their own Tt-over-WT bands, which is the only position evidence LAMP has",
      NA),
  row("mismatch_rtlamp2019", "lamp_mismatch_rtlamp2019_contradictions.csv", "report",
      "test-lamp-mismatch_rtlamp2019.R",
      "where the paper's prose and its own numbers disagree, recorded rather than reconciled: the claim is weakened, not averaged",
      "a reported-claim table, kept so a contradiction cannot be silently resolved in the paper's favour"),
  row("snp_lamp_lfd2025", "lamp_snp_lamp_lfd2025_detection.csv", "test",
      "test-lamp-snp_lamp_lfd2025.R",
      "a lateral-flow SNP-LAMP readout: a dipstick band at a time, so 'no band' before the window is censored, not negative",
      NA),
  row("snp_lamp_lfd2025", "lamp_snp_lamp_lfd2025_primers.csv", "config",
      "test-lamp-snp_lamp_lfd2025.R",
      "the primer set as design INPUT, so the detection table's allele column traces to an oligo",
      "no verdict is asserted from the primer table alone"),
  row("lamp_mismatch_scan_negative", "none", "blocked", NA_character_,
      "a documented NEGATIVE search: no paper in the corpus reports a position-wise mismatch scan for conventional LAMP",
      "the search itself is the result; there is no file to extract, and manufacturing one would invent evidence for a position rule LAMP does not have"),

  # --- ARMS / probes family --------------------------------------------------
  row("celegans2018", "celegans2018_arms_primers.csv", "config",
      "test-probes-celegans2018.R",
      "the primer table records the mismatch each allele-specific primer carries and the offset it sits at, so the Sanger observation about mismatch incorporation can be read back to a design",
      "no verdict is asserted from the primer table alone"),
  row("barley_arms2004", "arms_barley_arms2004_blocked.csv", "report",
      "test-arms2-negative.R",
      "the barley tetra-primer ARMS source is closed access, so the only thing that can be asserted is that nothing was landed and why",
      "closed access: no PMCID, isOpenAccess=false, Science CDN returns HTTP 403; the row records the block, not a number"),

  # --- CRISPR family ---------------------------------------------------------
  row("communbiol2024_crrna", "communbiol2024_crrna_activators.csv", "test",
      "test-crispr-communbiol.R",
      "Cas12a seed tolerance per position and per multiplicity: a double mismatch is not simply stronger than a single one, and single seed mismatches leak 37-97% of signal",
      NA),
  row("communbiol2024_crrna", "communbiol2024_crrna_ratio.csv", "test",
      "test-crispr-communbiol.R",
      "the mutant/WT endpoint fluorescence ratio the tolerance classes are computed from, so the class cannot be quoted without its ratio",
      NA),

  # --- KASP family -----------------------------------------------------------
  row("maize_kasp2020", "maize_kasp2020_assays.csv", "config",
      "test-kasp-platforms.R",
      "the maize KASP assay table as design INPUT: which SNP each assay targets, so a platform rate has a denominator",
      "no verdict is asserted from the assay table alone"),
  row("maize_kasp2020", "maize_kasp2020_outcomes.csv", "test",
      "test-kasp-platforms.R",
      "plant KASP conversion as per-assay outcome classes, so failure modes are not pooled into one rate",
      NA),
  row("rice2018", "rice2018_kasp_validation.csv", "test",
      "test-probes-rice2018.R",
      "KASP call accuracy against a reference genotype on a rice panel: the mismatched calls, not just the total",
      NA),
  row("rice2018", "rice2018_marker_cross.csv", "test",
      "test-probes-rice2018.R",
      "the marker x sample cross the validation rates are counted over, so a rate cannot be quoted without its matrix",
      NA),
  row("rice2024", "rice2024_kasp_error_rates.csv", "test",
      "test-probes-rice2024.R",
      "per-assay error rates on a second rice panel: how much of a KASP failure is design rather than genotype",
      NA),
  row("wheat_platform2019", "wheat_platform2019_platform_outcomes.csv", "test",
      "test-kasp-platforms.R",
      "platform-level outcomes on a plant panel: which failure class each assay landed in",
      NA),
  row("wheat_platform2019", "wheat_platform2019_cluster_metrics.csv", "test",
      "test-kasp-platforms.R",
      "the cluster metrics a genotype call is made from: a call is a clustering decision, so its thresholds are data",
      NA),

  # --- reported-claim tables -------------------------------------------------
  # A `report` row is what a paper ASSERTS, kept beside its measured tables so a
  # prose claim can never be promoted to a measurement. It carries a reason that
  # says what the claim is not (test-science-ledger.R requires one on `report`).
  row("rpa_mismatch_class", "rpa_mismatch_claims.csv", "report",
      "test-rpa-sources.R",
      "what the RPA literature CLAIMS about mismatch tolerance, kept separate from what it measured",
      "a reported-claim table, not an outcome table: its rows are what the papers assert, and none is a measured ratio"),
  row("sharplab2021", "sharplab2021_outcomes.csv", "report",
      "test-kasp-platforms.R",
      "the SHARP-lab outcome counts as REPORTED, with their own denominators; no test recomputes these rates yet",
      "a reported-claim table read only by the KASP platform tests, so its rates are claims, not measurements"),

  # --- ARMS/dPCR readouts ----------------------------------------------------
  row("sars2022", "sars2022_dmas_ct.csv", "test",
      "test-probes-sars2022.R",
      "Ct-based ARMS discrimination on SARS-CoV-2: the readout is a delta Ct, so 'got a band' cannot be read as a ratio",
      NA),
  row("sars2022", "sars2022_dmas_primers.csv", "config",
      "test-probes-sars2022.R",
      "the DMAS primer triples as design INPUT: which mismatch each allele-specific primer carries",
      "no verdict is asserted from the primer table; it makes the Ct table traceable to a design"),
  row("sharplab2021", "sharplab2021_assays.csv", "config",
      "test-kasp-platforms.R",
      "the SHARP-lab assay table as design INPUT, so an outcome rate has a denominator",
      "no verdict is asserted from the assay table alone"),
  row("vcf2caps2021", "vcf2caps2021_headline.csv", "test",
      "test-typing-vcf2caps2021.R",
      "the pipeline's own reported conversion headline, kept beside the recomputation so the two can be compared",
      NA),
  row("vcf2caps2021", "vcf2caps2021_verification.csv", "test",
      "test-typing-vcf2caps2021.R",
      "the per-marker verification that decides whether the headline is reproducible",
      NA)
))

# Refuse a row whose named test does not mention its extract, or whose extract
# is not on disk. The same checks test-science-ledger.R makes.
for (i in seq_len(nrow(seed))) {
  tf <- seed$test_file[[i]]
  if (is.na(tf)) next
  f <- file.path(root, "tests", "testthat", tf)
  if (!file.exists(f)) stop("seed names a test file that does not exist: ", tf)
  src <- readLines(f, warn = FALSE)
  for (one in trimws(strsplit(seed$extract[[i]], ",", fixed = TRUE)[[1L]])) {
    if (one == "none") next
    if (!file.exists(file.path(ext_dir, one))) stop("seed names a missing extract: ", one)
    if (!any(grepl(one, src, fixed = TRUE))) stop(tf, " does not read ", one)
  }
}

key <- function(d) paste(d$source_key, d$extract)
fresh <- seed[!key(seed) %in% key(led), ]
led <- rbind(led, fresh[, names(led), drop = FALSE])
utils::write.csv(led, led_path, row.names = FALSE, na = "")
message("wave5 ledger gap: ", nrow(fresh), " row(s) added, ", nrow(led), " total")
