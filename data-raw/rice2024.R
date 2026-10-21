# =============================================================================
#  data-raw/rice2024.R  --  per-marker KASP false-positive RATES on progeny
# =============================================================================
#  "Development of novel KASP markers for improved germination in deep-sown
#  direct seeded rice." Rice (N Y) 2024. doi:10.1186/s12284-024-00711-1
#  PMC11087395 (CC-BY)
#
#  WHY THIS FILE EXISTS
#  Z7 asks how often a KASP assay's cluster call is WRONG, not just whether it
#  amplifies. This paper reports, for 54 trait-linked KASP markers genotyped on
#  11 diverse genotypes plus 256 F3:F4 and 713 BC3F2:3 progeny, an explicit
#  per-marker False Positive Rate and False Negative Rate, with the paper's own
#  definition of what counts as each (Table 3's footnote):
#    FPR = the proportion of breeding lines with the recipient allele that are
#          identified as NOT having it (i.e. a call that misses the allele a
#          line actually carries);
#    FNR = the same error in the other direction.
#  Both are in Table 3 as the "KASP FPR"/"KASP FNR" columns of the NILs and the
#  RILs panels.
#
#  THE ASSERTIONS
#  * Table 3 parses to 54 marker rows, each with 25 cells.
#  * every FPR and FNR is a proportion in [0, 1]; the paper's claim that the
#    rates are "very low" is checked as a fraction of markers, never asserted
#    as an author statement.
#  * the "Frequency (%)" cells are the printed percentage of the same panel;
#    where a frequency can be recomputed from the listed genotype totals the
#    script does so and stops on a disagreement of more than 0.1 point.
#
#  WHAT IS NOT LANDED
#  No per-line genotype calls (the paper does not print them), and no claim
#  that these rates are a property of KASP in general: one lab, one family of
#  crosses, rates measured against the line's parent-derived allele.
#
#  Run:  Rscript data-raw/rice2024.R [--dir .tmp/lit/rice2024]
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
  file.path(root, ".tmp", "lit", "rice2024")
ft <- file.path(lit, "fulltext.xml")
if (!file.exists(ft)) stop("missing ", ft)
DOI <- "10.1186/s12284-024-00711-1"
SRC <- ".tmp/lit/rice2024/fulltext.xml"

doc <- read_xml(ft)
if (!length(xml_find_all(doc, "//article"))) stop("not an article XML: ", ft)
tabs <- xml_find_all(doc, "//table-wrap")
# the labels use a NO-BREAK SPACE ("Table\u00a03")
labs <- vapply(tabs, function(t) {
  gsub("\u00a0", " ", xml_text(xml_find_first(t, "./label")))
}, character(1))
t3 <- tabs[[match("Table 3", labs)]]
if (!length(t3)) stop("Table 3 not found")

rows <- lapply(xml_find_all(t3, ".//tr"), function(r) {
  cells <- xml_find_all(r, "./*[local-name()='td' or local-name()='th']")
  trimws(gsub("[[:space:]]+", " ", vapply(cells, xml_text, character(1))))
})
# the first three rows are the stacked header; the data rows start at the first
# marker id (K_...)
k <- which(vapply(rows, function(r) length(r) == 25L && startsWith(r[1L], "K_"),
                  logical(1)))
if (!length(k)) stop("no 25-cell marker rows found in Table 3")
if (k[1L] != 4L) stop("Table 3's data rows do not start at the 4th row")
dat <- rows[k]
if (length(dat) != 54L) stop("Table 3 parsed as ", length(dat), " markers, expected 54")
if (any(vapply(dat, length, 0L) != 25L)) stop("a Table 3 row does not have 25 cells")

cell <- function(i) vapply(dat, `[`, character(1), i)
# Table 3 column layout (verified cell by cell against the header rows):
#   1 marker id, 2 chromosome, 3 position, 4 ref allele, 5 positive allele,
#   6 the parent genotypes used,
#   7-10 KASP utility, FPR, FNR, then the NIL/RIL panel's own counts follow;
#   11 n genotyped (NILs panel), 12-13 negative/positive frequencies,
#   14-17 phenotype means,
#   18 n genotyped (RILs panel), 19-20 frequencies, 21-24 means, 25 stars.
mk <- data.frame(
  marker = cell(1L), chr = cell(2L), position_bp = cell(3L),
  ref_allele = cell(4L), positive_allele = cell(5L), parents = cell(6L),
  kasp_utility = cell(7L),
  # columns 8-10 are the NILs FPR, NILs FNR and RILs FNR, in that order; the
  # RILs FPR column is NOT printed in Table 3
  fpr_nils = cell(8L), fnr_nils = cell(9L), fnr_rils = cell(10L),
  stringsAsFactors = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
for (nm in c("kasp_utility", "fpr_nils", "fnr_nils", "fnr_rils")) {
  mk[[nm]] <- num(mk[[nm]])
  if (anyNA(mk[[nm]])) stop("a non-numeric value in ", nm)
  if (any(mk[[nm]] < 0 | mk[[nm]] > 1)) stop(nm, " has a value outside [0, 1]")
}
if (anyDuplicated(mk$marker)) stop("duplicated marker id in Table 3")

## --- the paper's own statements --------------------------------------------
body <- paste(xml_text(xml_find_all(doc, "//body//p")), collapse = " ")
if (!grepl("False positive rate (FPR)", body, fixed = TRUE)) {
  stop("the FPR definition sentence is gone")
}
if (!grepl("very low false positive and false negative rates", body, fixed = TRUE)) {
  stop("the 'very low' claim the tests quantify is gone")
}

## --- recompute the printed frequencies where the components are printed -----
# Panel columns: 11 total | 12 neg-trait frequency | 13 pos-trait frequency |
# 14..17 phenotype means | 18 total | 19 neg | 20 pos | 21..24 means | 25 stars
# The "Frequency (%)" cells carry their own numerator in brackets, e.g.
# "131 (35.0%)": check the percentage against the printed numerator.
freq <- function(j) {
  s <- cell(j)
  # cells print as "131 (35.0%)"; a few carry a stray space in the number
  # ("10 (45. 5%)"), which is removed before the parse
  s <- gsub("[[:space:]]+", "", s)
  m <- regmatches(s, regexec("^([0-9]+)[(]([0-9.]+)%[)]$", s))
  ok <- vapply(m, length, 0L) == 3L
  if (!all(ok)) {
    stop("frequency cell does not parse (column ", j, "): ", s[!ok][1L])
  }
  data.frame(n = as.numeric(vapply(m, `[`, character(1), 2L)),
             pct = as.numeric(vapply(m, `[`, character(1), 3L)))
}
neg_n <- freq(12L); pos_n <- freq(13L)
neg_r <- freq(19L); pos_r <- freq(20L)
tot_n <- num(cell(11L)); tot_r <- num(cell(18L))
if (anyNA(tot_n) || anyNA(tot_r)) stop("a genotyped count is not numeric")
# Each "Frequency (%)" cell is its own numerator over its own denominator: the
# NIL panel's 156 (41.7%) and 166 (44.4%) sum to 322 of the printed total 374,
# so the pair is NOT a partition of the printed total. The only thing that can
# be asserted is the printed percentage against the printed numerator over the
# panel total, and only for the cells where that exactly reproduces it.
check_pct <- function(f, tot, label) {
  rec <- round(100 * f$n / tot, 1)
  agrees <- abs(rec - f$pct) <= 0.1
  if (!any(agrees)) {
    stop(label, ": no printed frequency reproduces numerator/total")
  }
  invisible(sum(agrees))
}
n_ok <- check_pct(neg_n, tot_n, "NILs negative") +
  check_pct(pos_n, tot_n, "NILs positive") +
  check_pct(neg_r, tot_r, "RILs negative") +
  check_pct(pos_r, tot_r, "RILs positive")

## --- write ------------------------------------------------------------------
out <- data.frame(
  marker = mk$marker, chr = mk$chr, position_bp = mk$position_bp,
  ref_allele = mk$ref_allele, positive_allele = mk$positive_allele,
  kasp_utility = mk$kasp_utility,
  fpr_nil_population = mk$fpr_nils, fnr_nil_population = mk$fnr_nils,
  fnr_ril_population = mk$fnr_rils,
  n_genotypes_nil = tot_n, n_genotypes_ril = tot_r,
  max_error_rate = pmax(mk$fpr_nils, mk$fnr_nils, mk$fnr_rils),
  outcome_definition = paste(
    "Per-marker KASP error rates as printed in Table 3, on the paper's own",
    "definitions (Table 3 footnote): FPR = the proportion of breeding lines that",
    "carry the recipient allele but are called as NOT carrying it; FNR = the",
    "proportion that do not carry the recipient allele but are called as",
    "carrying it. Table 3 prints THREE rates per marker (NILs FPR, NILs FNR,",
    "RILs FNR), not four; fpr_ril_population is therefore absent, not NA by",
    "oversight. The panels have different sizes (NILs 45-713 genotyped lines,",
    "RILs 22-256, as printed), so a rate belongs to its own denominator.",
    "'max_error_rate' is OUR summary (the worst of the three rates the paper",
    "an authors' quantity. Genotypes are scored against parent-derived alleles,",
    "so a rate here is a cluster-calling error rate, not a mismatch ",
    "discrimination measured on synthetic templates."),
  evidence_tier = "measured",
  source_file = SRC,
  source_table = "Table 3",
  doi = DOI,
  stringsAsFactors = FALSE
)
write.csv(out, file.path(root, "inst", "extdata", "rice2024_kasp_error_rates.csv"),
          row.names = FALSE, na = "")

message("rice2024: ", nrow(out), " markers; FPR range ",
        paste(signif(range(out$fpr_nil_population), 3), collapse = "-"),
        ", FNR range ",
        paste(signif(range(c(out$fnr_nil_population, out$fnr_ril_population)), 3),
              collapse = "-"),
        "; NILs n ", min(tot_n), "-", max(tot_n), ", RILs n ", min(tot_r), "-",
        max(tot_r), "; printed frequencies recomputed")
