# =============================================================================
#  data-raw/lefever2019-supplement.R  --  what a SINGLE 3' mismatch is worth
# =============================================================================
#  Lefever S, Rihani A, Van der Meulen J, et al. "Cost-effective and robust
#  genotyping using double-mismatch allele-specific quantitative PCR."
#  Sci Rep 2019;9:2835. doi:10.1038/s41598-019-38581-z  PMC6377641 (CC-BY)
#
#  WHY THIS FILE EXISTS
#  Z7 proposes deliberate second mismatches for weak transition SNPs. The claim
#  "a weak transition SNP is strengthened" cannot be tested inside our own
#  ddG model (circular). Lefever gives the external measurement: allele-specific
#  qPCR against SYNTHETIC templates of known genotype, so dCq is clean.
#  Two results matter:
#    (a) a single 3'-terminal mismatch discriminates badly for most SNPs;
#    (b) an artificial mismatch at position 3 (4th nucleotide from the 3' end)
#        maximises discrimination, with no base-type effect at that position.
#
#  THE ASSERTION
#  Supplemental Table 1 reports, for 6 SNPs, dCq in three genotype groups and
#  ddCq = |dCq(group) - dCq(het)| for each homozygote. We recompute ddCq from
#  the means and require agreement to 0.01. That is the paper's own arithmetic:
#  if our transcription or the supplement's rounding were wrong, it fails.
#
#  WHAT IS NOT HERE
#  The per-position violin data (Fig 2) is a raster figure, not a table, so the
#  "position 3 is optimal" statement stays a CITED RULE, not a number we own.
#  It is recorded as evidence='primary literature' in mismatch_discrimination
#  and must not be promoted to 'fitted here'.
#
#  Run:  Rscript data-raw/lefever2019-supplement.R [--dir .tmp/lit/lefever2019]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!exists("root")) {
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else getwd()
  setwd(root)
}
suppressPackageStartupMessages(library(xml2))

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "lefever2019")
supp <- file.path(lit, "supplements", "41598_2019_38581_MOESM1_ESM.pdf")
ft <- file.path(lit, "fulltext.xml")
for (f in c(supp, ft)) if (!file.exists(f)) {
  stop("missing ", f, "\n  fetch with: Rscript data-raw/fetch-sources.R")
}

DOI <- "10.1038/s41598-019-38581-z"
PDB <- Sys.which("pdftotext")
if (!nzchar(PDB)) stop("pdftotext is required to read the supplemental PDF")

txt_path <- tempfile(fileext = ".txt")
system2(PDB, c("-layout", shQuote(supp), shQuote(txt_path)), stdout = FALSE)
lines <- readLines(txt_path, warn = FALSE)
unlink(txt_path)

## --- Supplemental Table 1: single-mismatch dCq by genotype ------------------
# Six rows of the shape
#   rs6734469  1.89 (± 0.00)  6.22 (± 4.83)  -12.83 (± 0.00)  4.33  14.72
# Columns: SNP | dCq A/B | dCq A/A | dCq B/B | ddCq AA-AB | ddCq BB-AB
num <- "(-?[0-9]+\\.[0-9]+)"
# The +/- symbol is not reliably ISO-8859-1 in pdftotext output, so match any
# short run of bytes between parenthesis and number rather than "\u00b1".
# NB: R's default regex engine is TRE and does NOT support \s; use [:space:].
sp <- "[[:space:]]+"
sd <- function(n) sprintf("\\([^)]*?%s\\)", n)
row_re <- paste0("^[[:space:]]*(rs[0-9]+)", sp, num, " ", sd(num), sp,
                 num, " ", sd(num), sp, num, " ", sd(num), sp, num, sp, num,
                 "[[:space:]]*$")
hit <- regmatches(lines, regexec(row_re, lines))
# full match + rsID + 8 numbers = 10 elements; the other rsID rows in the PDF
# (figure axis ticks) fail the numeric tail and are dropped here
hit <- Filter(function(h) length(h) == 10L, hit)
if (length(hit) != 6L) {
  stop("expected 6 SNP rows in Supplemental Table 1, parsed ", length(hit))
}
tab1 <- do.call(rbind, lapply(hit, function(h) {
  data.frame(
    snp = h[2],
    dCq_AB = as.numeric(h[3]), sd_AB = as.numeric(h[4]),
    dCq_AA = as.numeric(h[5]), sd_AA = as.numeric(h[6]),
    dCq_BB = as.numeric(h[7]), sd_BB = as.numeric(h[8]),
    ddCq_AA_AB_reported = as.numeric(h[9]),
    ddCq_BB_AB_reported = as.numeric(h[10]),
    stringsAsFactors = FALSE)
}))

# ddCq is defined in the legend as the absolute difference between the mean
# dCq of the heterozygous samples and that of the homozygous samples. The
# supplement rounds its INPUTS to 2 dp, so recomputing from the printed means
# can land one hundredth away from the printed ddCq (rs2273953: 1.62-0.68=0.94
# printed as 0.95). Consistency is therefore asserted at 0.01 + rounding, i.e.
# 0.011, NOT at 0.01; the 3-dp tail we keep is the reported value itself.
tab1$ddCq_AA_AB_recomputed <- round(abs(tab1$dCq_AA - tab1$dCq_AB), 3)
tab1$ddCq_BB_AB_recomputed <- round(abs(tab1$dCq_BB - tab1$dCq_AB), 3)
bad <- abs(tab1$ddCq_AA_AB_recomputed - tab1$ddCq_AA_AB_reported) > 0.0101 |
  abs(tab1$ddCq_BB_AB_recomputed - tab1$ddCq_BB_AB_reported) > 0.0101
if (any(bad)) {
  print(tab1[bad, ])
  stop("recomputed ddCq disagrees with the supplement's printed ddCq")
}

# The claim this table supports: a single mismatch separates genotypes only
# erratically. Quantify it on the paper's own definition (ddCq), not on a cut
# we invent: report how many SNPs clear 2 ddCq cycles on BOTH homozygotes.
tab1$min_ddCq <- pmin(tab1$ddCq_AA_AB_reported, tab1$ddCq_BB_AB_reported)
tab1$separates_at_2cycles <- tab1$min_ddCq >= 2

## --- full-text guard: quote the position-3 rule verbatim --------------------
body <- xml_text(xml_find_all(read_xml(ft), "//body//p"))
rule <- grep("mismatch on position three appears to be the optimal", body, value = TRUE)
base_effect <- grep("no differential effects when comparing different base pair",
                    body, value = TRUE)
if (!length(rule) || !length(base_effect)) {
  stop("full text no longer contains the position-3 and base-type statements")
}

out <- data.frame(
  snp = tab1$snp,
  dCq_het = tab1$dCq_AB, sd_het = tab1$sd_AB,
  dCq_hom_1 = tab1$dCq_AA, sd_hom_1 = tab1$sd_AA,
  dCq_hom_2 = tab1$dCq_BB, sd_hom_2 = tab1$sd_BB,
  ddCq_hom_1 = tab1$ddCq_AA_AB_reported,
  ddCq_hom_2 = tab1$ddCq_BB_AB_reported,
  min_ddCq = tab1$min_ddCq,
  separates_at_2cycles = tab1$separates_at_2cycles,
  design = "single 3'-terminal mismatch, no artificial mismatch",
  assay = "allele-specific qPCR (SYBR Green), two parallel AS reactions",
  template = "synthetic 60 bp templates, 10 000 molecules, known genotype",
  outcome_def = paste(
    "dCq = mean Cq(perfect-match AS primer) - Cq(mismatch AS primer), by",
    "genotype group; ddCq = |dCq(homozygote) - dCq(heterozygote)|, the paper's",
    "measure of discriminating power. 'separates_at_2cycles' is OUR summary",
    "(min ddCq >= 2), not a claim by the authors."),
  evidence = "primary literature",
  calibrated = FALSE,
  doi = DOI,
  stringsAsFactors = FALSE
)
write.csv(out, file.path(root, "inst", "extdata", "lefever_single_mismatch.csv"),
          row.names = FALSE, na = "")

# The design rule as cited prose, so downstream code cannot silently retype it.
write.csv(data.frame(
  rule = c("artificial mismatch position",
           "base-type effect at position 3",
           "single 3'-terminal mismatch"),
  statement = c(sub("^.*?Data for six SNPs", "Data for six SNPs", rule),
                base_effect,
                paste("the discrimination power of a single 3'-terminal mismatch",
                      "is not sufficient for reliable genotype calling (Suppl. Table 1)")),
  evidence = "primary literature",
  verified_against = "fulltext.xml (verbatim grep)",
  doi = DOI
), file.path(root, "inst", "extdata", "lefever_design_rules.csv"),
row.names = FALSE, na = "")

message(
  "lefever2019: ", nrow(out), " SNPs, ddCq recomputed and matched; ",
  sum(out$separates_at_2cycles), "/", nrow(out),
  " separate the homozygotes by >= 2 cycles on both alleles with a single",
  " 3'-terminal mismatch"
)
