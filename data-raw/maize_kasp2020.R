# =============================================================================
#  data-raw/maize_kasp2020.R  --  Jagtap, Vikal & Johal 2020, genome-wide maize
# =============================================================================
#  Jagtap AB, Vikal Y, Johal GS.  "Genome-Wide Development and Validation of
#  Cost-Effective KASP Marker Assays for Genetic Dissection of Heat Stress
#  Tolerance in Maize."  Int J Mol Sci 2020;21:7386.
#  DOI 10.3390/ijms21197386, PMC7582619.
#
#  This paper is the CLEANEST statement of the distinction in the family: it
#  defines the two rates itself and reports both for the same 100 SNVs --
#
#    "The SNP conversion rate is composed of two components--the design success
#     rate and the work success rate. ... In our study, we developed KASP
#     markers with an assay design success rate of 100% and a work success rate
#     of 92%."
#
#  and then splits its own 29% of failures into named count categories:
#
#    71 (71%) polymorphic -> usable
#    21 (21%) "monomorphic and/or heterozygous type" -> unusable
#     8 ( 8%) "failed to generate a useful amplification signal" -> unusable
#
#  71 + 21 + 8 = 100, and 71 + 21 = 92 is the paper's own work success rate.
#  Those two arithmetic identities are asserted here.  The 92% is the fraction
#  of DESIGNED SNVs that produced a usable two-allele call on the parent pair,
#  NOT the fraction that amplified -- amplification is 92 + 21 = 93%.
#
#  The supplement (supporting data S1/S2) lists the 100 SNV sites and their
#  primer triples, but it carries NO per-marker outcome column.  So the extract
#  is at STUDY level: per-marker outcomes are deliberately not invented.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)
dir.create(file.path(root, ".tmp", "lit", "maize_kasp2020"), recursive = TRUE,
           showWarnings = FALSE)
source(file.path(root, "data-raw", "kasp_lib.R"))

DATE_VERIFIED <- "2026-10-21"
DOI <- "10.3390/ijms21197386"
OUT <- fs::path(root, "inst", "extdata")

s1 <- kasp_input(root, "maize_kasp2020",
                 "supplements/s001/SupplementaryData/supportingdataS1.xlsx",
                 "maize_kasp2020_sequences.xlsx", "maize_kasp2020 supporting data S1")
s2 <- kasp_input(root, "maize_kasp2020",
                 "supplements/s001/SupplementaryData/supportingdataS2.xlsx",
                 "maize_kasp2020_primers.xlsx", "maize_kasp2020 supporting data S2")
xml_path <- kasp_input(root, "maize_kasp2020", "fulltext.xml",
                       "maize_kasp2020_fulltext.xml", "maize_kasp2020 full text")

## --- supporting data S1: the 100 SNV sites -----------------------------------
stopifnot(file.size(s1) > 10000L)
rows1 <- lapply(kasp_sheets(s1)[[1L]], kasp_row)
hdr1 <- rows1[[2L]]
if (!identical(unname(hdr1[c("B", "C", "D", "F")]),
               c("Chromosome", "Marker name (KASP assay) ", "Position", "Polymorphism"))) {
  stop("supporting data S1's header is not where it was: ",
       paste(hdr1[1:7], collapse = " | "))
}
d1 <- Filter(function(r) grepl("^[0-9]+$", r[["A"]]) && grepl("^MKAM_", r[["C"]]), rows1)
sites <- do.call(rbind, lapply(d1, function(r) {
  data.frame(no = as.integer(r[["A"]]), chromosome = r[["B"]], marker = trimws(r[["C"]]),
             position = r[["D"]], left_flank = r[["E"]], polymorphism = r[["F"]],
             right_flank = r[["G"]], stringsAsFactors = FALSE)
}))
if (nrow(sites) != 100L || !identical(sort(sites$no), 1:100)) {
  stop("supporting data S1 does not hold 100 numbered sites (parsed ",
       nrow(sites), ")")
}
# the paper's own description: "10 SNPs per chromosome" over all 10 chromosomes
per_chr <- table(sites$chromosome)
if (length(per_chr) != 10L || !all(per_chr == 10L)) {
  stop("the 100 sites are not 10 per chromosome as the text says: ",
       paste(names(per_chr), per_chr, collapse = " "))
}

## --- supporting data S2: the primer triples ----------------------------------
rows2 <- lapply(kasp_sheets(s2)[[1L]], kasp_row)
hdr2 <- rows2[[2L]]
if (!identical(unname(hdr2[c("A", "B", "C", "D", "E", "F")]),
               c("KASP assay marker name", "Primer_AlleleX  (5'-3')",
                 "Primer_AlleleY (5'-3')", "Primer_Common (5'-3')",
                 "AlleleX (LM11)", "AlleleY (CML25)"))) {
  stop("supporting data S2's header is not where it was")
}
d2 <- Filter(function(r) grepl("^MKAM_", r[["A"]]) && nzchar(r[["D"]]), rows2)
primers <- do.call(rbind, lapply(d2, function(r) {
  data.frame(marker = r[["A"]], primer_allele_x = r[["B"]], primer_allele_y = r[["C"]],
             primer_common = r[["D"]], allele_x = r[["E"]], allele_y = r[["F"]],
             stringsAsFactors = FALSE)
}))
# The paper says "All selected 100 SNP sites (100%) were converted"; supporting
# data S2 prints primer triples for only 92 of them.  The 8 absences are NOT
# read as design failures: the paper's 8 named failures are "failed to
# generate a useful amplification signal", i.e. a WORK failure, and the
# supplement does not say the 8 are the same assays.  The 92 is therefore kept
# as its own observed count, and the design rate is asserted from the paper's
# text rather than from the file.
n_s2 <- nrow(primers)
if (n_s2 != 92L) {
  stop("supporting data S2 holds ", n_s2, " primer triples; the study-level ",
       "rows and the report assume 92 -- re-read the supplement if this changes")
}
if (!all(primers$marker %in% sites$marker)) {
  stop("supporting data S2 names a marker that supporting data S1 does not")
}

# right join: every one of the 100 sites stays in the table, and the 8 without
# a printed primer triple get NA primer columns rather than being dropped
joined <- merge(sites, primers, by = "marker", all.x = TRUE)
joined <- joined[order(joined$no), ]
if (nrow(joined) != 100L) stop("the site/primer join lost rows")
# the SNP printed in S1 must sit at the two alleles printed in S2
has_prim <- !is.na(joined$allele_x)
allele_ok <- vapply(seq_len(nrow(joined)), function(i) {
  if (!has_prim[i]) return(TRUE)
  p <- strsplit(joined$polymorphism[i], "/")[[1L]]
  identical(sort(p), sort(c(joined$allele_x[i], joined$allele_y[i])))
}, logical(1))
if (!all(allele_ok)) {
  stop(sum(!allele_ok), " of the 100 sites print SNP alleles that disagree with ",
       "the primer alleles in S2: ", paste(joined$marker[!allele_ok], collapse = ", "))
}
# the allele-specific primers must differ at their 3' base
as_x <- substr(joined$primer_allele_x, nchar(joined$primer_allele_x),
               nchar(joined$primer_allele_x))
as_y <- substr(joined$primer_allele_y, nchar(joined$primer_allele_y),
               nchar(joined$primer_allele_y))
if (any(has_prim & as_x == as_y)) {
  stop(sum(has_prim & as_x == as_y), " primer pairs do not differ at the 3' base")
}
# The X/Y labels are NOT tied to the allele letters: on many rows the 3' base
# of primer X is the complement of AlleleX (the synthesis may carry the other
# strand), and on at least one row X and Y are swapped.  What must hold is the
# UNORDERED pair, because that is what makes the assay two-allele.
comp_base <- c(A = "T", T = "A", C = "G", G = "C")
pair_ok <- vapply(seq_len(nrow(joined)), function(i) {
  if (!has_prim[i]) return(TRUE)
  got <- sort(c(as_x[i], as_y[i]))
  printed <- sort(c(joined$allele_x[i], joined$allele_y[i]))
  compd <- sort(c(comp_base[[as_x[i]]], comp_base[[as_y[i]]]))
  identical(got, printed) || identical(compd, printed)
}, logical(1))
if (!all(pair_ok)) {
  stop(sum(!pair_ok), " primer pairs do not carry the site's two alleles in some ",
       "orientation: ", paste(joined$marker[!pair_ok], collapse = ", "))
}
joined$allele_labels_note <- ifelse(
  has_prim & (as_x == joined$allele_x & as_y == joined$allele_y),
  "",
  "3'-base labels are complemented and/or X/Y swapped relative to AlleleX/AlleleY; the unordered pair matches")
as_x[!has_prim] <- NA_character_
as_y[!has_prim] <- NA_character_

## --- the full text: the two rates and the three categories -------------------
stopifnot(file.size(xml_path) > 20000L)
txt <- paste(readLines(xml_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
if (!grepl("<article", txt, fixed = TRUE)) stop("not a JATS article")
flat <- kasp_norm(gsub("<[^>]+>", " ", txt))

# the paper's own definitions
if (!grepl("The SNP conversion rate is composed of two components", flat, fixed = TRUE)) {
  stop("the design/work definition sentence is absent; this source's whole ",
       "value is that sentence")
}
if (!grepl("an assay design success rate of 100% and a work success rate of 92%",
           flat, fixed = TRUE)) {
  stop("the '100% design / 92% work' sentence is absent or has changed")
}

# the counted split: 71 polymorphic, 21 monomorphic/heterozygous, 8 no signal
a_poly <- kasp_find(flat, "of these MKAMs were found to be polymorphic")
site_poly <- substr(flat, a_poly$start - 40L, a_poly$start + 260L)
if (!grepl("71 [(]71%[)] .*polymorphic", site_poly, perl = TRUE)) {
  stop("the 71-polymorphic sentence has changed: ", site_poly)
}
n_poly <- 71L
if (!grepl("21 [(]21%[)] .*monomorphic and/or heterozygous", site_poly, perl = TRUE)) {
  stop("the 21-monomorphic/heterozygous sentence has changed: ", site_poly)
}
n_mono <- 21L
if (!grepl("8 [(]8%[)] failed to generate a useful amplification signal",
           site_poly, perl = TRUE)) {
  stop("the 8-no-amplification sentence has changed: ", site_poly)
}
n_noamp <- 8L

n_attempted <- nrow(sites)   # 100 SNV sites selected, all 100 designed
n_design <- 100L  # the paper: "an assay design success rate of 100%"
if (n_poly + n_mono + n_noamp != n_attempted) {
  stop("the three outcome categories do not sum to the 100 sites: ",
       n_poly, " + ", n_mono, " + ", n_noamp, " = ", n_poly + n_mono + n_noamp)
}
if (!grepl("All selected 100 SNP sites (100%) were converted to KASP assay markers",
           flat, fixed = TRUE)) {
  stop("the 'all 100 sites were converted' sentence is absent")
}
if (n_design != n_attempted) {
  stop("design success (", n_design, ") is not 100% of the ",
       n_attempted, " sites")
}
# The paper's own work-success number is (polymorphic + monomorphic/het) /
# designed = 92/100.  The fraction that produced a USABLE two-allele call is
# 71/100.  Both are computed and both are written; neither is renamed.
work_rate <- (n_poly + n_mono) / n_design
if (round(work_rate * 100) != 92L) {
  stop("(71 + 21)/100 is not the paper's 92% work success rate (", work_rate, ")")
}
usable_rate <- n_poly / n_design
# and amplification = usable + monomorphic = 92 / 100, so the 8 no-signal assays
# are the ONLY ones that did not amplify
amp_rate <- (n_poly + n_mono) / n_design
if (round(amp_rate * 100) != 92L) {
  stop("the amplified fraction is not 92/100 (", amp_rate, ")")
}

## --- write: one row per designed assay (site + primers, NO outcome) ---------
out <- data.frame(
  marker = joined$marker,
  chromosome = joined$chromosome,
  position = joined$position,
  polymorphism = joined$polymorphism,
  left_flank_50bp = joined$left_flank,
  right_flank_50bp = joined$right_flank,
  primer_allele_x = joined$primer_allele_x,
  primer_allele_y = joined$primer_allele_y,
  primer_common = joined$primer_common,
  allele_x_lm11 = joined$allele_x,
  allele_y_cml25 = joined$allele_y,
  as_3p_x = as_x,
  as_3p_y = as_y,
  allele_labels_note = joined$allele_labels_note,
  primer_triple_in_s2 = has_prim,
  design_success = TRUE,
  assay_outcome = NA_character_,
  stringsAsFactors = FALSE
)
out$assay_outcome_note <- paste(
  "NOT IN THE FILE: supporting data S1/S2 carry the site, the flanks and the",
  "primer triple only. The paper reports outcomes as three study-level counts",
  "(71/21/8), so per-assay outcomes are left NA rather than assigned."
)
out$doi <- DOI
out$evidence_tier <- "measured"
out$source_file <- basename(s1)
out$source_table <- "supporting data S1 (sites) joined to S2 (primers) on MKAM name"
out$date_verified <- DATE_VERIFIED
out$conversion_definition <- paste(
  "DESIGN: 100/100 = 1.00. The paper: 'All selected 100 SNP sites (100%) were",
  "converted to KASP assay markers' and 'an assay design success rate of 100%'.",
  "Supporting data S2 prints a primer triple for 92 of the 100 sites; the ",
  "other 8 have no entry there, which the paper does not explain."
)
out$outcome_definition <- paste(
  "WORK (paper's own definition): 'the number of SNP sites that can generate",
  "genotype calls via primers to the number of SNP sites with successfully",
  "designed primers'. Reported for the PARENT PAIR LM11 x CML25 and, later, the",
  "90-plant F2: 71/100 = 71%, which the paper calls a 'work success rate of 92%'",
  "only after counting the 21 monomorphic/heterozygous assays as amplifying.",
  "Both numbers are kept apart in maize_kasp2020_outcomes.csv."
)
kasp_write_checked(out, fs::path(OUT, "maize_kasp2020_assays.csv"))

## --- write: the study-level outcomes -----------------------------------------
dfn_design <- paste(
  "DESIGN (paper's definition): 'the ratio of the SNP sites that can be used to",
  "design primers to the total number of SNP sites'. Here 100/100."
)
dfn_work <- paste(
  "WORK (paper's definition): 'the number of SNP sites that can generate",
  "genotype calls via primers to the number of SNP sites with successfully",
  "designed primers'. The paper's 92% is 92/100 = (71 polymorphic + 21",
  "monomorphic-or-heterozygous)/100; that is the fraction that produced ANY",
  "call, not the fraction that produced a TWO-ALLELE (usable) call, which is",
  "71%. This extract records both so neither can be quoted as the other."
)
dfn_amp <- paste(
  "AMPLIFICATION: 100 - 8 = 92 of 100 assays 'generated a useful amplification",
  "signal'. Identical to the paper's 92% work rate by its own arithmetic",
  "(71 + 21), which is why the paper and this extract must be read carefully:",
  "amplifying is not the same as being a usable marker."
)
dfn_mono <- paste(
  "USABLE (paper's category): '71 (71%) of these MKAMs were found to be",
  "polymorphic'. This is the only count in the paper that means 'a working,",
  "informative assay on this material'."
)
outs <- data.frame(
  quantity = c("snv_sites_attempted", "assays_designed", "design_success_rate",
               "amplified", "amplification_rate", "work_success_rate_paper",
               "usable_polymorphic", "unusable_monomorphic_or_het",
               "no_amplification_signal", "primer_triples_in_s2"),
  n_numerator = c(100L, n_design, n_design, n_poly + n_mono, n_poly + n_mono,
                  n_poly + n_mono, n_poly, n_mono, n_noamp, n_s2),
  n_denominator = c(100L, 100L, 100L, 100L, 100L, 100L, 100L, 100L, 100L, 100L),
  value_printed = c("100", "100", "100%", "92", "92%", "92%", "71 (71%)",
                    "21 (21%)", "8 (8%)", "92 of 100"),
  rate_recomputed = c(1, 1, 1, (n_poly + n_mono) / 100, (n_poly + n_mono) / 100,
                      (n_poly + n_mono) / 100, n_poly / 100, n_mono / 100,
                      n_noamp / 100, n_s2 / 100),
  outcome_definition = c(
    "The selected RNA-seq SNV sites this study starts from ('a total of 100 SNPs site spanning all 10 chromosomes').",
    dfn_design, dfn_design, dfn_work, dfn_amp, dfn_work, dfn_mono,
    paste("UNUSABLE (paper's category): '21 (21%) were monomorphic and/or",
          "heterozygous type and were deemed unusable as markers'. Amplified but",
          "not informative on LM11 x CML25."),
    paste("NO SIGNAL (paper's category): 'The remaining 8 (8%) failed to",
          "generate a useful amplification signal and thus were unusable'."),
    paste("SUPPLEMENT COVERAGE, not an outcome: supporting data S2 prints a",
          "primer triple for 92 of the 100 sites. This is what the file shows,",
          "NOT a design rate -- the paper claims 100% design success and the",
          "missing 8 are not identified anywhere.")),
  note = c(
    "stated as 'a total of 100 SNPs site spanning all 10 chromosomes'",
    "paper's claim; supporting data S2 prints only 92 primer triples",
    "the paper's 'assay design success rate of 100%'",
    "derived: 100 - 8 (the only category that did not amplify)",
    "derived: identical to the paper's work rate by its own arithmetic",
    "'with an assay design success rate of 100% and a work success rate of 92%'",
    "the paper's own count; this is the number a breeder would call 'working'",
    "the paper's own count",
    "the paper's own count",
    "counted from supporting data S2"),
  doi = DOI, evidence_tier = "measured", source_file = basename(xml_path),
  source_table = "Results 2.3 and Discussion",
  date_verified = DATE_VERIFIED, stringsAsFactors = FALSE
)
outs <- outs[c("quantity", "n_numerator", "n_denominator", "value_printed",
               "rate_recomputed", "outcome_definition", "note", "doi",
               "evidence_tier", "source_file", "source_table", "date_verified")]
kasp_write_checked(outs, fs::path(OUT, "maize_kasp2020_outcomes.csv"))

message(
  "maize_kasp2020: ", nrow(out), " sites, ", n_s2, " with primer triples; ",
  "design 100/100; categories ", n_poly, "+", n_mono, "+", n_noamp, " = 100; ",
  "work ", n_poly + n_mono, "/100 = ", work_rate, " (paper 92%), usable ",
  n_poly, "/100"
)
