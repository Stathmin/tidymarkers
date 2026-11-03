# =============================================================================
#  data-raw/snapshot_multiplex_design.R  --  ROADMAP §U: SNaPshot priors
# =============================================================================
#  The SNaPshot / single-base-extension scheme in R/marker_snapshot.R carries
#  design rules whose only sources were two application papers (Fanis 2014,
#  Namdar Aligoodarzi 2023), so every criterion would land as `rule_of_thumb`.
#  This generator puts the rules that a named source STATES into a table with
#  provenance.  Nothing is wired into R/marker_scheme_registry.R.
#
#  A rule is included only when the full text states it.  A number printed in a
#  paper's own Table is `measured` (an observation of that panel); a design
#  recommendation the authors make is `stated`; a value we derive by reading a
#  table counts is `derived`.  Rules with no source stay out of the table and
#  are named in the `absent` block at the bottom of this script, which is also
#  what data-raw/proposed/u-snapshot-priors.md lists as rule_of_thumb.
#
#  Evidence tiers: "moderate" = a stated rule plus the panel that follows it in
#  one named source; "weak" = a stated rule in a single source, no independent
#  confirmation; "none" = we hold the statement but no quantitative outcome.
#
#  Verbatim sources checked (Europe PMC fullTextXML / PMC HTML):
#    Atanasovska 2012 PLoS One 7:e48167 10.1371/journal.pone.0048167 PMC3482202
#    Di Cristofaro 2010 J Mol Diagn 12:453 10.2353/jmoldx.2010.090222 PMC2893629
#    Fanis 2014 BMC Genomics 15:108 10.1186/1471-2164-15-108 PMC3922441
#    Coutinho 2014 PLoS One 9:e93292 10.1371/journal.pone.0093292 PMC3966882
#    Mrazkova 2021 Diagnostics 11:301 10.3390/diagnostics11020301 PMC7918147
#    Namdar Aligoodarzi 2023 Iran Biomed J 27:46 10.52547/ibj.3856 PMC9971712
#    Zhao 2026 BMC Genomics 27 10.1186/s12864-026-13173-0 PMC13644189
#    Bouvet 2020 Int J Mol Sci 21:9650 10.3390/ijms21249650 PMC7766892
#    Riordan 2022 Plant Methods 18 10.1186/s13007-022-00923-w PMC9250251
#  Cited by those papers but READ ONLY AS TITLE BY US (no full-text access):
#    Pastinen 1997 Genome Res 7:606 10.1101/gr.7.6.606 (no PMCID in EPMC)
#    Lindblad-Toh 2000 Nat Genet 24:381 (paywalled; not in EPMC)
#    Norton 2002 Hum Mutat 19:628 (not in EPMC)
#    Sanchez 2003 (ref [35] of Fanis 2014); its numbers are only seen second-hand.
#
#  Run:  Rscript data-raw/snapshot_multiplex_design.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()

rules <- data.frame(
  rule_id = c(
    "sbe_3prime_adjacent",
    "sbe_3prime_no_flexibility",
    "sbe_core_length",
    "sbe_spacing_min",
    "sbe_spacing_used_4bp",
    "sbe_spacing_used_ct5bp",
    "sbe_tail_5prime",
    "sbe_tail_keeps_tm",
    "sbe_amp_primer_tm",
    "sbe_sbe_tm_gc_window",
    "sbe_inter_primer_tm_delta",
    "sbe_primer_ratio_reported",
    "sbe_peak_height_floor",
    "sbe_het_balance_floor",
    "sbe_dye_yield_asymmetry",
    "sbe_structure_check",
    "sbe_nonhomologous_window",
    "sbe_no_template_control",
    "sbe_apply_format"
  ),
  parameter = c(
    "extension primer 3' position", "extension primer 3' position",
    "core primer length", "same-dye fragment spacing", "same-dye fragment spacing",
    "same-dye fragment spacing", "tail", "tail composition",
    "amplification primer Tm", "SBE primer Tm & GC", "inter-primer Tm spread",
    "peak height vs input", "peak height floor", "heterozygote peak ratio",
    "allele peak height asymmetry", "secondary structure", "template window",
    "template-free extension", "format"
  ),
  rule_text = c(
    "the 3' base of the minisequencing primer is the base immediately 5' (adjacent) of the variant, so the first incorporated ddNTP is the allele",
    "the primer position is fixed by the variant; only length, tail and strand are free design variables",
    "hybridise the core primer over 15-25 nt of template (also stated as 18-25, 17-25); the 3' base stays the base before the SNP",
    "oligos carrying the same fluorescent label must differ by at least 3 nt (a 4 nt spacing was designed into one panel)",
    "panels keep neighbouring extension products at least 4 nt apart (AmericaPlex26; also 4, 8 and 9 nt in the PAH panel)",
    "poly(CT) tails are used to force at least 5 nt between products in the same dye channel",
    "a non-homologous neutral sequence (or poly-CT / poly-T / poly-dC) is added at the 5' end, not the 3' end",
    "the tail is added only to reach the wanted size; Tm is computed on the core annealing region with the tail excluded",
    "multiplex amplification primers were designed to Tm 58-63 C",
    "SBE primers: Tm 53-60 C, GC 40-55 %",
    "difference in Tm between SBE primers in one panel <= 3 C",
    "peak height follows template input and primer concentration; panels tune both",
    "electropherogram peaks kept within 700-4000 rfu",
    "minor/major peak height ratio < 0.1 marks allelic imbalance",
    "the two dye terminators of a heterozygote can differ strongly in height (A vs C reported as ~1:3.5)",
    "extension primers are screened for hairpins, self-dimers and heterodimers; 3'-end structures can prime template-free extension",
    "amplicon length is set by the template window that the primer sits in, not by a fixed ceiling",
    "a multiplex SBE with the primer mix and no template is run to detect template-free products",
    "SNaPshot MultiReady kit: one premix carries the polymerase and the four fluorescent ddNTPs, so the four dNTPs are never user-mixed"
  ),
  value = c(
    "3' base at variant-1", "n/a", "17-25", ">= 3", ">= 4", ">= 5", "5' end", "excluded",
    "58-63", "53-60", "<= 3", "not quantified", "700-4000", "< 0.1", "1:3.5 (A vs C)",
    "checked", "not quantified", "required", "premix"
  ),
  range = c(
    NA, NA, "15-25", "3-Inf", "4-Inf", "5-Inf", NA, NA, "58-63", "53-60", NA, NA,
    "700-4000", "0-0.1", NA, NA, NA, NA, NA
  ),
  unit = c(
    "position", "n/a", "nt", "nt", "nt", "nt", "end", "n/a", "degC", "degC", "degC",
    "n/a", "rfu", "ratio", "fold", "n/a", "n/a", "n/a", "n/a"
  ),
  outcome_definition = c(
    "the allele is read as the identity of the single ddNTP added 3' of the primer",
    "primer length and strand shift the peak; the 3' base cannot move",
    "core length that the papers' tables print for working primers",
    "resolved (non-overlapping) peaks for equal-labelled products",
    "resolved peaks in the same dye channel",
    "resolved peaks in the same dye channel",
    "separation and balanced base composition of the extension products",
    "Tm of the core only, so the tail does not change annealing",
    "mPCR primers anneal in one multiplex",
    "SBE primers anneal at one extension temperature (50 C in the kit cycle)",
    "one SBE temperature works for the whole panel",
    "uniform peak heights across loci",
    "peaks above the analysis threshold and below saturation",
    "a heterozygote is called, not flagged as imbalance",
    "heterozygote gives two resolved peaks, not one dominant peak",
    "no background peaks from hairpin/dimer extension",
    "the amplicon spans the variant with room for both PCR primers",
    "no peaks in the template-free reaction",
    "the dNTP constraint is kit-internal, not a design free parameter"
  ),
  evidence_tier = c(
    "moderate", "moderate", "moderate", "moderate", "weak", "weak", "moderate",
    "moderate", "weak", "weak", "weak", "none", "weak", "weak", "weak",
    "moderate", "weak", "moderate", "weak"
  ),
  source_file = c(
    paste("PMC3482202 :: Results/Discussion; PMC3922441 :: Methods;",
      "PMC9971712 :: Table 3 caption; PMC13644189 :: Methods"),
    "PMC3482202 :: Discussion",
    "PMC9971712 Table 3; PMC3482202 Table 2; PMC3922441 Table 3",
    paste("PMC3482202 :: Methods; PMC3922441 :: Discussion ('minimum of 3 nt')"),
    "PMC3966882 :: Methods (AmericaPlex26); PMC9971712 :: Discussion",
    "PMC13644189 :: Methods",
    paste("PMC9971712 :: Discussion (40 nt neutral); PMC3922441 :: Methods (poly-C);",
      "PMC3966882 :: Table 1 (poly-CT); PMC2893629 :: Table 2 (poly-T)"),
    "PMC3482202 :: Methods (dC tail)",
    "PMC9971712 :: Methods",
    "PMC13644189 :: Methods",
    "PMC13644189 :: Methods",
    "PMC2893629 :: Discussion; PMC9971712 :: Discussion; PMC7766892 :: Results",
    "PMC2893629 :: Discussion",
    "PMC13644189 :: Methods (GeneMapper ID-X settings)",
    "PMC2893629 :: Discussion",
    paste("PMC3482202 :: Discussion (3'-end structures); PMC3922441 :: Methods;",
      "PMC13644189 :: Methods (AutoDimer)"),
    "PMC3922441 :: Methods (amplicons 60-649 bp); PMC9971712 Methods (143-882 bp)",
    "PMC7918147 :: Methods; PMC3482202 :: Table S2",
    "PMC2893629 :: Methods; PMC7918147 :: Methods"
  ),
  doi = c(
    rep("10.1371/journal.pone.0048167", 2), "10.52547/ibj.3856",
    "10.1371/journal.pone.0048167", "10.1371/journal.pone.0093292",
    "10.1186/s12864-026-13173-0", "10.52547/ibj.3856", "10.1371/journal.pone.0048167",
    "10.52547/ibj.3856", "10.1186/s12864-026-13173-0", "10.1186/s12864-026-13173-0",
    "10.2353/jmoldx.2010.090222", "10.2353/jmoldx.2010.090222",
    "10.1186/s12864-026-13173-0", "10.2353/jmoldx.2010.090222",
    "10.1371/journal.pone.0048167", "10.1186/1471-2164-15-108",
    "10.3390/diagnostics11020301", "10.2353/jmoldx.2010.090222"
  ),
  pmcid = c(
    rep("PMC3482202", 2), "PMC9971712", "PMC3482202", "PMC3966882", "PMC13644189",
    "PMC9971712", "PMC3482202", "PMC9971712", "PMC13644189", "PMC13644189",
    "PMC2893629", "PMC2893629", "PMC13644189", "PMC2893629", "PMC3482202",
    "PMC3922441", "PMC7918147", "PMC2893629"
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  !anyDuplicated(rules$rule_id),
  all(c("rule_id", "parameter", "rule_text", "value", "range", "unit",
    "outcome_definition", "evidence_tier", "source_file", "doi", "pmcid") %in%
    names(rules)),
  all(rules$evidence_tier %in% c("moderate", "weak", "none")),
  all(grepl("^PMC[0-9]+$", rules$pmcid)),
  all(grepl("^10\\.", rules$doi))
)

## `value` may carry a qualifier (">= 3", "< 0.1"), a range ("700-4000") or a
## colon ratio ("1:3.5"); keep the first number and, for a ratio, the second
num <- function(id) {
  v <- gsub("<", "", rules$value[match(id, rules$rule_id)])
  as.numeric(sub("^[^0-9]*-?([0-9]+(\\.[0-9]+)?).*$", "\\1", v))
}
den <- function(id) {
  v <- gsub("<", "", rules$value[match(id, rules$rule_id)])
  as.numeric(sub(".*:([0-9]+(\\.[0-9]+)?).*$", "\\1", v))
}
## the spacing ladder is strictly ordered: same dye, same product, more spacing
stopifnot(
  num("sbe_spacing_min") == 3,
  num("sbe_spacing_used_4bp") == 4,
  num("sbe_spacing_used_ct5bp") == 5,
  num("sbe_spacing_used_ct5bp") > num("sbe_spacing_used_4bp"),
  num("sbe_spacing_used_4bp") > num("sbe_spacing_min")
)
## Tm window and spread are the Zhao 2026 design targets, verbatim
tm <- regmatches(rules$range[match("sbe_sbe_tm_gc_window", rules$rule_id)],
  regexpr("[0-9]+-[0-9]+", rules$range[match("sbe_sbe_tm_gc_window", rules$rule_id)])
)[1]
tm <- as.numeric(strsplit(tm, "-")[[1L]])
stopifnot(
  identical(tm, c(53, 60)),
  num("sbe_inter_primer_tm_delta") == 3,
  tm[2L] - tm[1L] > num("sbe_inter_primer_tm_delta"),
  num("sbe_peak_height_floor") == 700,
  as.numeric(sub(".*-([0-9]+).*$", "\\1",
    rules$value[match("sbe_peak_height_floor", rules$rule_id)])) == 4000,
  num("sbe_het_balance_floor") == 0.1,
  is.na(suppressWarnings(num("sbe_primer_ratio_reported"))),
  num("sbe_dye_yield_asymmetry") == 1,
  den("sbe_dye_yield_asymmetry") == 3.5
)
## every rule with a numeric Tm window clears the kit's 50 C extension step
stopifnot(min(tm) > 50)

out <- file.path(root, "inst", "extdata", "snapshot_multiplex_design_rules.csv")
utils::write.csv(rules, out, row.names = FALSE)

absent <- data.frame(
  rule_id = c(
    "sbe_extension_temperature", "sbe_extension_time_cycles",
    "sbe_dnaddition_dntp_purity", "sbe_primer_conc_absolute",
    "sbe_copy_number_ceiling", "sbe_amplicon_max",
    "sbe_polymerase_fidelity", "sbe_hetero_height_symmetric"
  ),
  why_absent = c(
    "papers quote their own reaction temperature (50 C, kit cycle) but none states it as a design rule for the primer",
    "the kit cycle (30 x 96 C/10 s, 50 C/5 s, 60 C/30 s) is protocol, not a primer-design rule",
    "the SNaPshot premix supplies the ddNTPs; no paper states a user-side dNTP-mix constraint",
    "concentrations are reported per panel (0.01-0.4 uM in the PAH panel) but never as a rule",
    "no source gives a template copy-number or input ceiling for SBE",
    "published amplicons run 60-882 bp; no source states a maximum length",
    "no source reports a polymerase requirement beyond the kit premix",
    "measured asymmetry is dye-dependent (1:3.5 for A vs C); no source sets a symmetric height rule"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(absent, file.path(root, "data-raw", "proposed", "u-snapshot-absent.csv"),
  row.names = FALSE)

cat("snapshot_multiplex_design: ", nrow(rules), " rules, ",
  sum(rules$evidence_tier == "moderate"), " moderate / ",
  sum(rules$evidence_tier == "weak"), " weak / ",
  sum(rules$evidence_tier == "none"), " none; ", nrow(absent),
  " rule_of_thumb items left out\n", sep = "")
