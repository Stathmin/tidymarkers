# =============================================================================
#  data-raw/jsbbs2018.R  --  HRM marker funnel, observed dTm, selection rule
# =============================================================================
#  Source (identifier is PROVENANCE; every number is re-read):
#    "Selection criteria for SNP loci to maximize robustness of high-resolution
#    melting analysis for variant genotyping in soybean"
#    Breed Sci 2018;68(4).  DOI 10.1270/jsbbs.18048, PMID 30369824, PMC6198901.
#
#  Fetched: .tmp/lit/jsbbs2018/fulltext.xml (EPMC fullTextXML, 119,717 B) and
#  .tmp/lit/jsbbs2018/supplements/68_18048_2.xlsx (8 sheets), read with
#  data-raw/lib/xlsx-read.R (no spreadsheet package installed).
#
#  WHY: a PLANT (soybean / Oryza) HRM marker-development paper that publishes
#  per-marker observed dTm for every marker that worked.  That distribution,
#  not a textbook "1 C", is the empirical basis for a dTm threshold.
#
#  WHAT THIS SCRIPT DOES NOT MODEL
#   * No funnel stage becomes "P(HRM marker works)": denominators differ
#     (48 / 32 / 24), the stages are not nested and the population is one
#     soybean cross plus one Oryza comparison.
#   * |dG-deg3bp| is a FREE-ENERGY index in kcal/mol - never a dTm threshold.
#   * No dTm below is a threshold WE chose; the observed minimum is reported
#     as descriptive, with no stated false-negative rate.
#
#  KNOWN INTERNAL INCONSISTENCY (recorded, not forced):
#   Prose: "average dTmobs values of non-NNNs-HRM and NNNs-HRM markers were
#   0.43 C and 0.63 C ... (11 combinations ... df = 10, t value = 7.6)".
#   Supplemental Table 4 has 13 marker pairs, of which ONE has no polymorphism
#   (M50/M66), ONE has a PCR failure (M56/M72) and ONE has no polymorphism
#   without the substitution (M64/M80) -> at most 11 pairs have both dTm values
#   and only 10 have a non-zero non-NNNs dTm.  From the file the 10 usable
#   pairs give mean 0.431 / 0.689 C and t = 8.15 on df = 9.  Only the 0.43 C
#   reproduces; 0.63 C, df = 10 and t = 7.6 do not, so they are carried with
#   matched_from_file = FALSE instead of being adjusted.
#
#  ASSERTIONS: A1 Table 1 totals + percentages sum to 100; A2 funnel 38/48,
#  24/32, 18/24 against the printed percentages; A3 Table 2 dTmobs = |Tm_a-Tm_f|
#  and r(dTmobs, |dG3bp|) = 0.647; A4 Table 5 dTmobs = |Tm_F-Tm_T|, 54..100 bp,
#  |dG3bp| >= 1.43, mean 0.74 C; A5 Table 4 effect = NNNs - non-NNNs and
#  10/13 usable; A6 Table 6 672 genotypes, 10 discordant, 9 unknown, 653 same.
#
#  Run:  Rscript data-raw/jsbbs2018.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

xml_path <- fs::path(root, ".tmp", "lit", "jsbbs2018", "fulltext.xml")
xlsx_path <- fs::path(root, ".tmp", "lit", "jsbbs2018", "supplements", "68_18048_2.xlsx")
for (p in c(xml_path, xlsx_path)) if (!file.exists(p)) stop("Missing ", p)

DOI <- "10.1270/jsbbs.18048"
out_dir <- fs::path(root, "inst", "extdata")

doc <- xml2::read_xml(xml_path)
paras <- vapply(
  xml2::xml_find_all(doc, "//*[local-name()='p']"),
  function(p) gsub("[[:space:]]+", " ", xml2::xml_text(p)), character(1))
if (!any(grepl(DOI, xml2::xml_text(doc), fixed = TRUE)))
  stop("fulltext.xml does not mention ", DOI)

source("data-raw/lib/xlsx-read.R")
sh <- xlsx_sheets(xlsx_path)
num_col <- function(rows, k) vapply(rows, function(r)
  suppressWarnings(as.numeric(xcell(r, k))), numeric(1))
chr_col <- function(rows, k) vapply(rows, function(r) xcell(r, k), character(1))

# --- A1: Table 1 --------------------------------------------------------------
tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
t1 <- Filter(function(t) grepl("Frequencies of SNP classes", xml2::xml_text(
  xml2::xml_find_first(t, ".//*[local-name()='caption']"))), tw)
if (length(t1) != 1L) stop("Table 1 not found")
cells <- lapply(xml2::xml_find_all(t1[[1L]], ".//*[local-name()='tr']"), function(r)
  trimws(xml2::xml_text(xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']"))))
data_rows <- Filter(function(v) length(v) == 7L && grepl("^[0-9.]+%$", v[[3L]]), cells)
if (length(data_rows) != 3L) stop("Expected 3 dataset rows in Table 1, got ", length(data_rows))
snp_class <- tibble::tibble(
  dataset = vapply(data_rows, `[[`, character(1), 1L),
  genus   = vapply(data_rows, `[[`, character(1), 2L),
  pct_class_I   = as.numeric(sub("%$", "", vapply(data_rows, `[[`, character(1), 3L))),
  pct_class_II  = as.numeric(sub("%$", "", vapply(data_rows, `[[`, character(1), 4L))),
  pct_class_III = as.numeric(sub("%$", "", vapply(data_rows, `[[`, character(1), 5L))),
  pct_class_IV  = as.numeric(sub("%$", "", vapply(data_rows, `[[`, character(1), 6L))),
  total_snps = as.numeric(gsub(",", "", vapply(data_rows, `[[`, character(1), 7L))))
# The printed percentages are rounded and, for AF/TF, do not sum to 100 in the
# source itself (AF 104.5, TF 94.9). That is a property of the paper; assert the
# exact printed values so a silent re-parse cannot drift.
pct_sum <- with(snp_class, pct_class_I + pct_class_II + pct_class_III + pct_class_IV)
if (!isTRUE(all.equal(round(pct_sum, 2), c(104.5, 94.9, 100.0), tolerance = 1e-9)))
  stop("Table 1 percentages sum to ", paste(round(pct_sum, 2), collapse = "/"),
       " (source prints 104.5/94.9/100)")
if (!isTRUE(all.equal(snp_class$total_snps, c(8372, 6142, 4138822))))
  stop("Table 1 totals are ", paste(snp_class$total_snps, collapse = "/"))
if (!any(grepl("8372 loci for AF SNPs and 6142 loci for TF SNPs", paras, fixed = TRUE)))
  stop("Methods do not restate the 8372 / 6142 totals")

# --- A2: funnel ---------------------------------------------------------------
f <- c(exp1_designed = 48L, exp1_working = 38L, exp2_pairs = 32L, exp2_working = 24L,
       exp3_designed = 24L, exp3_polymorphic = 18L, exp4_genotypes = 672L,
       exp4_discordant = 10L, exp4_unknown = 9L, exp4_identical = 653L)
if (!any(grepl("38 (79.2%) primer pairs", paras, fixed = TRUE))) stop("38 (79.2%) not found")
if (!any(grepl("comprised NNNs-HRM markers", paras, fixed = TRUE)) ||
    !any(grepl("24 (75.0%) showed single amplification and polymorphism", paras,
               fixed = TRUE)))
  stop("24 (75.0%) not found")
if (!any(grepl("The remaining 18 primers (75%) showed co-dominant segregation", paras,
               fixed = TRUE))) stop("18 (75%) not found")
if (!any(grepl("3 showed no polymorphism between parents and another 3 showed a dominant segregation pattern",
               paras, fixed = TRUE))) stop("Exp. 3 drop counts not found")
funnel <- tibble::tibble(
  experiment = c("Exp. 1", "Exp. 2", "Exp. 3"),
  rule = c("48 primer pairs on AF SNPs, <= 80 bp amplicon; success = single band AND polymorphic",
           "16 NNNs-HRM + 16 non-NNNs-HRM pairs at the same 16 SNP loci",
           "24 NNNs-HRM pairs, |dG-deg3bp| > 1.43 kcal/mol, 54-100 bp products"),
  n_designed = c(f[["exp1_designed"]], f[["exp2_pairs"]], f[["exp3_designed"]]),
  n_success = c(f[["exp1_working"]], f[["exp2_working"]], f[["exp3_polymorphic"]]),
  pct_reported = c(79.2, 75.0, 75.0))
if (any(abs(100 * funnel$n_success / funnel$n_designed - funnel$pct_reported) > 0.05))
  stop("funnel percentages do not reproduce from the counts")

# --- A3: Supplemental Table 2 -------------------------------------------------
t2_all <- sh[["Supplemental_Table2"]][-(1:2)]
# The last row of the sheet is the footnote, not a marker.
if (length(t2_all) != 39L || grepl("^a)", chr_col(t2_all[39L], "A")) == FALSE)
  stop("Supplemental Table 2 layout changed: expected 38 marker rows + 1 footnote")
t2_rows <- t2_all[1:38]
m2 <- tibble::tibble(
  marker_id = as.integer(chr_col(t2_rows, "A")), chromosome = chr_col(t2_rows, "B"),
  product_bp = num_col(t2_rows, "D"), ref_triplet = chr_col(t2_rows, "E"),
  snp_triplet = chr_col(t2_rows, "F"), tm_aokimame = num_col(t2_rows, "G"),
  tm_fukuyutaka = num_col(t2_rows, "H"), dtm_obs = num_col(t2_rows, "I"),
  dg3bp = num_col(t2_rows, "J"), dtm_pred = num_col(t2_rows, "K"))
if (nrow(m2) != 38L) stop("Table 2 should have 38 marker rows, got ", nrow(m2))
if (any(!m2$marker_id %in% 1:48) || anyDuplicated(m2$marker_id))
  stop("Table 2 marker IDs are not a unique subset of 1:48")
if (length(setdiff(1:48, m2$marker_id)) != 48L - 38L)
  stop("Table 2 does not account for 48 - 38 failed markers")
bad <- abs(m2$dtm_obs - abs(m2$tm_aokimame - m2$tm_fukuyutaka)) > 1e-9
if (any(bad, na.rm = TRUE))
  stop("Table 2: dTmobs != |Tm_Aokimame - Tm_Fukuyutaka| in row(s) ", paste(which(bad), collapse = ", "))
r_dg <- stats::cor(m2$dtm_obs, m2$dg3bp, use = "complete.obs")
if (abs(r_dg - 0.647) > 1e-3)
  stop("Table 2 gives r(dTmobs, |dG3bp|) = ", round(r_dg, 4), ", prose says 0.647")
m2$tm_order <- ifelse(m2$tm_aokimame > m2$tm_fukuyutaka, "Aokimame_higher", "Fukuyutaka_higher")

# --- A4: Supplemental Table 5 -------------------------------------------------
t5_rows <- sh[["Supplemental_Table5"]][3:20]
m5 <- tibble::tibble(
  marker_id = chr_col(t5_rows, "A"), chromosome = chr_col(t5_rows, "B"),
  product_bp = num_col(t5_rows, "F"), nnn_position = chr_col(t5_rows, "G"),
  nnn_combination = chr_col(t5_rows, "H"), dg3bp_nnns = num_col(t5_rows, "I"),
  dg3bp_nonnns = num_col(t5_rows, "J"), tm_fukuyutaka = num_col(t5_rows, "K"),
  tm_toyoshirome = num_col(t5_rows, "L"), dtm_obs = num_col(t5_rows, "M"))
if (nrow(m5) != 18L) stop("Table 5 should have 18 marker rows, got ", nrow(m5))
if (any(abs(m5$dtm_obs - abs(m5$tm_fukuyutaka - m5$tm_toyoshirome)) > 1e-9))
  stop("Table 5: dTmobs != |Tm_F - Tm_T|")
if (any(m5$product_bp < 54 | m5$product_bp > 100)) stop("Table 5 sizes outside 54..100 bp")
if (any(m5$dg3bp_nnns < 1.43 - 1e-9)) stop("Table 5 has |dG3bp| below the 1.43 kcal/mol rule")
if (!any(grepl("54\u2013100 bp", paras, fixed = TRUE))) stop("prose 54-100 bp not found")
if (!any(grepl("mean, 0.74°C", paras, fixed = TRUE))) stop("prose mean 0.74 C not found")
if (abs(mean(m5$dtm_obs) - 0.74) > 0.005)
  stop("Table 5 mean dTmobs is ", round(mean(m5$dtm_obs), 4), ", prose says 0.74")
r_len <- stats::cor(m5$product_bp, m5$dtm_obs)
if (!any(grepl("r = \u22120.75", paras, fixed = TRUE)) || r_len > -0.7)
  stop("dTmobs is not negatively correlated with product length as stated")

# --- A5: Supplemental Table 4 -------------------------------------------------
t4_rows <- sh[["Supplemental_Table4"]][5:17]
m4 <- tibble::tibble(
  marker_nonnns = chr_col(t4_rows, "A"), marker_nnns = chr_col(t4_rows, "B"),
  nnn_position = chr_col(t4_rows, "G"), nnn_combo = chr_col(t4_rows, "H"),
  product_bp = num_col(t4_rows, "I"), dg3bp_after = num_col(t4_rows, "J"),
  dg3bp_before = num_col(t4_rows, "K"), dtm_nonnns = num_col(t4_rows, "T"),
  dtm_nnns = num_col(t4_rows, "U"), effect = num_col(t4_rows, "V"))
if (nrow(m4) != 13L) stop("Table 4 should have 13 marker pairs, got ", nrow(m4))
if (any(abs(m4$effect - (m4$dtm_nnns - m4$dtm_nonnns)) > 1e-9, na.rm = TRUE))
  stop("Table 4: effect != NNNs - non-NNNs")
usable <- !is.na(m4$dtm_nonnns) & !is.na(m4$dtm_nnns)
if (sum(usable) != 10L) stop("Table 4: expected 10 pairs with both dTm values, got ", sum(usable))
tt4 <- stats::t.test(m4$dtm_nnns[usable], m4$dtm_nonnns[usable], paired = TRUE)
m4_stats <- c(n_pairs = nrow(m4), n_both_values = sum(usable),
              mean_nonnns = mean(m4$dtm_nonnns[usable]), mean_nnns = mean(m4$dtm_nnns[usable]),
              t_stat = unname(tt4$statistic), df = unname(tt4$parameter))
if (abs(m4_stats[["mean_nonnns"]] - 0.43) > 0.005)
  stop("Table 4 mean non-NNNs dTm is ", round(m4_stats[["mean_nonnns"]], 4), ", prose says 0.43")
m4_match <- c(mean_nnns = abs(m4_stats[["mean_nnns"]] - 0.63) <= 0.005,
              df = identical(m4_stats[["df"]], 10), t = abs(m4_stats[["t_stat"]] - 7.6) < 0.05)

# --- A6: Supplemental Table 6 -------------------------------------------------
t6_rows <- sh[["Supplemental_Talbe6"]][4:7]
lev <- c("Homozygous for Fukyutaka allele", "Homozygous for Aokimame allele",
         "Heterozygous", "Unknown/Negative")
conc <- do.call(rbind, lapply(t6_rows, function(r) {
  idx <- match(xcell(r, "A"), lev)
  v <- vapply(c("B", "C", "D", "E"), function(k) {
    z <- xcell(r, k)
    if (is.na(z) || !grepl("^[0-9]+", z)) 0 else as.numeric(sub("^([0-9]+).*$", "\\1", z))
  }, numeric(1))
  c(diag = v[idx], total = sum(v))
}))
# Table rows sum to 674; the paper reports 672 (footnote 10/672, prose 653/672).
# The discrepancy is in the source, so assert it rather than hide it.
if (sum(conc[, "total"]) != 674L)
  stop("Table 6 totals ", sum(conc[, "total"]), " genotypes, source rows say 674")
cat("jsbbs2018 NOTE: Table 6 rows sum to 674 while prose/footnotes use 672\n")
# 10/672 and 9/672 are footnotes in the supplement sheet, not in the article XML.
fn6 <- grep("672", vapply(sh[["Supplemental_Talbe6"]], function(r) xcell(r, "A"), character(1)),
            value = TRUE, fixed = TRUE)
if (!any(grepl("10/672, 1.49%", fn6, fixed = TRUE)) ||
    !any(grepl("9/672, 1.34%", fn6, fixed = TRUE)))
  stop("Table 6 footnotes 10/672 (1.49%) and 9/672 (1.34%) not found in the supplement")
if (f[["exp4_genotypes"]] - f[["exp4_discordant"]] - f[["exp4_unknown"]] != f[["exp4_identical"]])
  stop("672 - 10 - 9 != 653")
if (!any(grepl("653 out of 672 genotypes (97.2%)", paras, fixed = TRUE)))
  stop("prose 653/672 (97.2%) not found")
if (abs(100 * f[["exp4_discordant"]] / f[["exp4_genotypes"]] - 1.49) > 0.005 ||
    abs(100 * f[["exp4_unknown"]] / f[["exp4_genotypes"]] - 1.34) > 0.005)
  stop("Table 6 footnote percentages do not reproduce from 10/672 and 9/672")

# --- write --------------------------------------------------------------------
TIER <- "measured"
readr::write_csv(
  dplyr::mutate(m2, doi = DOI, evidence_tier = TIER,
                source_file = "supplements/68_18048_2.xlsx", sheet = "Supplemental_Table2",
                outcome_definition = paste(
                  "Outcome = dTmobs (C), the absolute difference between the maximum -dF/dt",
                  "peaks of the two parental PCR fragments under HRM, for ONE primer pair at",
                  "ONE SNP (Exp. 1, AF SNPs, <= 80 bp). Homozygote-versus-homozygote",
                  "discrimination in plants, not heterozygote detection and not mutation",
                  "scanning. tm_order records which parent melted higher (the absolute",
                  "value hides the sign flips).")),
  fs::path(out_dir, "hrm_jsbbs2018_dtm_exp1.csv"))
readr::write_csv(
  dplyr::mutate(m5, doi = DOI, evidence_tier = TIER,
                source_file = "supplements/68_18048_2.xlsx", sheet = "Supplemental_Table5",
                outcome_definition = paste(
                  "Outcome = dTmobs (C) between the two parents for the 18 NNNs-HRM markers",
                  "of Exp. 3 that segregated co-dominantly. Pre-selected on |dG-deg3bp| >",
                  "1.43 kcal/mol and 54-100 bp products, so the distribution is CONDITIONED",
                  "on the paper's selection rule and understates unselected markers.")),
  fs::path(out_dir, "hrm_jsbbs2018_dtm_exp3.csv"))
readr::write_csv(
  dplyr::mutate(m4, doi = DOI, evidence_tier = TIER,
                source_file = "supplements/68_18048_2.xlsx", sheet = "Supplemental_Table4",
                outcome_definition = paste(
                  "Outcome = dTmobs (C) for a matched pair of markers at the same SNP,",
                  "with and without the NNN substitution in the primer. effect =",
                  "dTmobs(NNNs) - dTmobs(non-NNNs) is a DESIGN gain, not a melting",
                  "constant. NA where the table prints 'No polymorphism', 'PCR failure'",
                  "or 'Not caluclated' [sic].")),
  fs::path(out_dir, "hrm_jsbbs2018_nnn_effect.csv"))
readr::write_csv(
  dplyr::mutate(funnel, doi = DOI, evidence_tier = TIER, source_file = "fulltext.xml",
                outcome_definition = paste(
                  "Outcome = a primer pair that passed that experiment's own stated test",
                  "(single band AND parental polymorphism; co-dominant segregation in",
                  "Exp. 3). Denominators are the pairs designed per experiment; the three",
                  "experiments are NOT nested and must not be pooled into a design rate.")),
  fs::path(out_dir, "hrm_jsbbs2018_funnel.csv"))
readr::write_csv(
  dplyr::mutate(snp_class, doi = DOI, evidence_tier = TIER, source_file = "fulltext.xml",
                sheet = "Table1",
                outcome_definition = paste(
                  "Outcome = percentage of SNPs in each of the four classes in a whole-genome",
                  "variant set (I = C/T and G/A transitions, II = other transversions,",
                  "III = G/C, IV = A/T). A population composition, not an HRM outcome; it",
                  "bounds how often a class III/IV SNP - the hardest for HRM - appears.")),
  fs::path(out_dir, "hrm_jsbbs2018_snp_class_frequencies.csv"))
readr::write_csv(
  tibble::tibble(
    criterion = c("maximise |dG-deg3bp|", "|dG-deg3bp| selection threshold",
                  "amplicon length", "C/T SNP optimal NNN", "G/A SNP optimal NNN",
                  "C/A SNP optimal NNN", "G/T SNP optimal NNN", "class III/IV SNP NNN"),
    value = c("mutate the 5' or 3' nearest neighbour to the base that maximises |dG-deg3bp|",
              "1.43", "< 80", "G (5') or A (3')", "T (5') or C (3')",
              "G (either side)", "C (either side)", "no clear trend"),
    units = c(NA, "kcal/mol", "bp", NA, NA, NA, NA, NA),
    basis = c("r = 0.647 with dTmobs, P < 0.01, 39 markers", "Exp. 3 marker set",
              "paper's primer-design cap", "Fig. 3", "Fig. 3", "Fig. 3", "Fig. 3",
              "Fig. 3"),
    outcome_definition = paste(
      "Rule, not an outcome. |dG-deg3bp| is the change in Gibbs free energy over the",
      "3 bp (i-NNN, i, i+NNN) using SantaLucia (1998) unified NN parameters at 37 C;",
      "it is a free-energy proxy for dTm, NOT a melting-temperature threshold. The",
      "1.43 kcal/mol cut is the authors' choice for Exp. 3; they state only that",
      "> 1.0 kcal/mol with < 80 bp 'would generally guarantee dTm > 0.5 C'."),
    evidence_tier = "measured", doi = DOI, source_file = "fulltext.xml"),
  fs::path(out_dir, "hrm_jsbbs2018_selection_rule.csv"))
readr::write_csv(
  tibble::tibble(
    quantity = c("n_markers_exp1", "dtm_exp1_min_c", "dtm_exp1_median_c", "dtm_exp1_max_c",
                 "n_markers_exp3", "dtm_exp3_mean_c", "r_dtmobs_dg3bp_exp1",
                 "r_dtmobs_length_exp3", "n_pairs_table4", "n_pairs_table4_both_values",
                 "mean_dtm_nonnns_table4_c", "mean_dtm_nnns_table4_c",
                 "t_table4_recomputed", "df_table4_recomputed", "prose_mean_dtm_nonnns_c",
                 "prose_mean_dtm_nnns_c", "prose_t_table4", "prose_df_table4"),
    value = c(nrow(m2), min(m2$dtm_obs), stats::median(m2$dtm_obs), max(m2$dtm_obs),
              nrow(m5), mean(m5$dtm_obs), r_dg, r_len, m4_stats[["n_pairs"]],
              m4_stats[["n_both_values"]], m4_stats[["mean_nonnns"]], m4_stats[["mean_nnns"]],
              m4_stats[["t_stat"]], m4_stats[["df"]], 0.43, 0.63, 7.6, 10),
    recomputed_from_file = c(rep(NA, 12), TRUE, TRUE, TRUE, unname(m4_match)),
    outcome_definition = paste(
      "Outcome = descriptives of the observed dTmobs distribution. dtm_exp1_min_c is the",
      "smallest dTm the authors still called polymorphic; it is NOT a validated detection",
      "floor (no false-negative rate is reported at it). The Table 4 rows record a known",
      "non-reproducibility: 0.63 C, df = 10 and t = 7.6 in the prose cannot be recomputed",
      "from the 13-row supplementary table, so they are carried with",
      "recomputed_from_file = FALSE rather than adjusted."),
    evidence_tier = TIER, doi = DOI, source_file = "supplements/68_18048_2.xlsx"),
  fs::path(out_dir, "hrm_jsbbs2018_summary.csv"))

message("jsbbs2018: ", f[["exp1_working"]], "/48, ", f[["exp2_working"]], "/32, ",
        f[["exp3_polymorphic"]], "/24; dTm Exp1 ", round(min(m2$dtm_obs), 3), "-",
        round(max(m2$dtm_obs), 3), " C; r=", round(r_dg, 3), "; Table4 ",
        m4_stats[["n_both_values"]], "/13 usable, t=", round(m4_stats[["t_stat"]], 2),
        " df=", m4_stats[["df"]], " prose-match ", sum(m4_match), "/3 -- assertions passed")
