# =============================================================================
#  data-raw/meng2026_soybean_outcomes.R  --  Z8.2: per-assay soybean KASP outcomes
# =============================================================================
#  Meng X, Fu Z, Zhao W, Wu X, Ma C, Feng Y, Du S, Zhao X, Wang Y, Han Y.
#  "Association Mapping of Seedling Resistance to Fusarium graminearum Root Rot
#  and Development of KASP Assays in Soybean", Plants 2026, 15:2479,
#  10.3390/plants15162479, PMC13516154, PMID 42654881.
#
#  This is the only soybean source found that reaches a per-assay (not a
#  per-paper) reading of KASP conversion, and it is a small one: of the six
#  exonic SNPs it prioritised (Table 1), two were submitted to KASP and both
#  gave usable allelic discrimination.  Four SNPs were never submitted, so their
#  outcome is `not_attempted`, NOT `failed`; this file does not manufacture a
#  failure denominator out of them.
#
#  Numbers below are printed in the paper's Results (Table 1 and Section 3.5):
#    * 336 accessions = 17 HR + 119 R + 191 S + 9 HS (HR 0<DSI<=20, R 20<DSI<=40,
#      S 40<DSI<=80, HS >80)
#    * six exonic SNPs in six candidate genes, one favourable allele each
#    * KASP-S17_32244510 (C/T, Glyma.17g202500) and KASP-S18_55105706 (A/G,
#      Glyma.18g266700) were designed and both resolved homozygous genotype
#      clusters ("successfully converted into KASP assays")
#    * screening efficiency 26 -> 19 accessions = 73.08 %;
#      35 -> 26 accessions = 74.29 % (printed as 73.08 % and 74.29 %)
#    * 308 of 336 accessions had unambiguous homozygous calls at both loci
#
#  What it does NOT report: a failure list (no assay in this paper failed), the
#  primer sequences for the four unsubmitted SNPs, and any reason a KASP assay
#  would be discarded.  Supplementary Table S7 (primer sequences) was not
#  readable (MDPI returns 403 to a scripted fetch), so primer sequences are not
#  cited here; only the outcomes and the counts above are.
#
#  Run:  Rscript data-raw/meng2026_soybean_outcomes.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

doi <- "10.3390/plants15162479"
pmcid <- "PMC13516154"

pheno <- c(HR = 17L, R = 119L, S = 191L, HS = 9L)
n_panel <- 336L
n_called_both <- 308L

# Table 1: the six exonic SNPs the paper prioritised, with the favourable
# allele and the mean DSI of each allele class.
snps <- data.frame(
  snp = c(
    "S09_39290818", "S12_37955440", "S17_32166476",
    "S17_32244510", "S18_2618132", "S18_55105706"
  ),
  gene = c(
    "Glyma.09g168300", "Glyma.12g220200", "Glyma.17g202000",
    "Glyma.17g202500", "Glyma.18g033900", "Glyma.18g266700"
  ),
  ref_allele = c("T", "A", "T", "C", "A", "A"),
  alt_allele = c("C", "G", "A", "T", "G", "G"),
  ref_mean_dsi = c(45.30, 47.58, 45.96, 45.97, 46.78, 46.40),
  alt_mean_dsi = c(38.64, 38.82, 37.15, 33.19, 39.00, 32.33),
  delta_dsi = c(-6.66, -8.77, -8.82, -12.78, -7.78, -14.07),
  stringsAsFactors = FALSE
)
# Table 1 prints ΔDSI as a negative number sized REF-ALT (its note says
# ALT-REF, so the sign is the paper's, the magnitude is REF-ALT).

# Section 3.5: the two SNPs that were submitted to KASP.
kasp <- data.frame(
  snp = c("S17_32244510", "S18_55105706"),
  assay_id = c("KASP-S17_32244510", "KASP-S18_55105706"),
  favourable_genotype = c("TT", "GG"),
  n_favourable = c(26L, 35L),
  n_favourable_resistant = c(19L, 26L),
  efficiency_printed = c(73.08, 74.29),
  stringsAsFactors = FALSE
)

stopifnot(
  # panel arithmetic as the paper states it
  n_panel == sum(pheno),
  unname(pheno[["HR"]] + pheno[["R"]]) == 136L,
  n_called_both < n_panel,
  # Table 1 is self-consistent
  nrow(snps) == 6L,
  !anyDuplicated(snps$snp),
  # Table 1 is self-consistent to the printed precision: |ΔDSI| = REF - ALT,
  # within 0.01 at two loci (the paper's means are rounded, so an exact match
  # happens to hold at four of six only).
  all(abs(abs(snps$delta_dsi) - (snps$ref_mean_dsi - snps$alt_mean_dsi)) <= 0.011),
  # the ALT allele is the favourable one at every locus
  all(snps$alt_mean_dsi < snps$ref_mean_dsi),
  # the two KASP loci are in Table 1 and take its favourable allele
  all(kasp$snp %in% snps$snp),
  # the favourable genotype is the ALT allele repeated (homozygous)
  identical(
    snps$alt_allele[match(kasp$snp, snps$snp)],
    substring(kasp$favourable_genotype, 1L, 1L)
  ),
  identical(substring(kasp$favourable_genotype, 1L, 1L), substring(kasp$favourable_genotype, 2L, 2L)),
  # screening efficiency is reproducible from the counts it is computed from
  all(round(100 * kasp$n_favourable_resistant / kasp$n_favourable, 2) ==
    kasp$efficiency_printed),
  all(kasp$n_favourable_resistant <= kasp$n_favourable),
  # selection, not failure: the four untouched SNPs stay out of the outcome rows
  nrow(kasp) == 2L,
  nrow(snps) - nrow(kasp) == 4L
)

basis_worked <- paste(
  "Table 1 lists the SNP; Section 3.5 reports the assay designed and resolved",
  "homozygous genotype clusters on the 336-accession panel"
)

outcomes <- data.frame(
  source = rep("meng2026", nrow(snps)),
  assay_id = ifelse(
    snps$snp %in% kasp$snp,
    kasp$assay_id[match(snps$snp, kasp$snp)],
    paste0("KASP-", snps$snp)
  ),
  target = paste0(snps$snp, " (", snps$gene, ", ", snps$ref_allele, "/",
    snps$alt_allele, ")"),
  outcome = ifelse(snps$snp %in% kasp$snp, "worked", "not_attempted"),
  basis = ifelse(snps$snp %in% kasp$snp, basis_worked,
    "Table 1 only; Section 3.5 selected two of six SNPs for KASP, this one was not submitted"),
  n_genotyped = ifelse(snps$snp %in% kasp$snp, n_panel, NA_integer_),
  screening_eff_n = ifelse(snps$snp %in% kasp$snp,
    kasp$n_favourable[match(snps$snp, kasp$snp)], NA_integer_),
  screening_eff_worked = ifelse(snps$snp %in% kasp$snp,
    kasp$n_favourable_resistant[match(snps$snp, kasp$snp)], NA_integer_),
  screening_efficiency_pct = ifelse(snps$snp %in% kasp$snp,
    kasp$efficiency_printed[match(snps$snp, kasp$snp)], NA_real_),
  doi = doi,
  pmcid = pmcid,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(outcomes) == 6L,
  identical(sum(outcomes$outcome == "worked"), 2L),
  identical(sum(outcomes$outcome == "not_attempted"), 4L),
  !any(outcomes$outcome %in% c("failed", "no_amplification", "no_polymorphism", "ambiguous"))
)

out <- "inst/extdata/meng2026_soybean_per_assay.csv"
utils::write.csv(outcomes, out, row.names = FALSE)
cat(
  "meng2026_soybean_outcomes:", nrow(outcomes), "SNPs;",
  sum(outcomes$outcome == "worked"), "worked KASP,",
  sum(outcomes$outcome == "not_attempted"), "not attempted;",
  "efficiency", paste0(outcomes$screening_efficiency_pct[!is.na(outcomes$screening_efficiency_pct)], "%", collapse = " / "),
  "\n"
)
