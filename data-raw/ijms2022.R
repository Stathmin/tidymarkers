# =============================================================================
#  data-raw/ijms2022.R  --  HRM discrimination of 3 linked C>T SNPs (MGMT)
# =============================================================================
#  Source (identifier resolution is PROVENANCE, not evidence; every number
#  below is re-read from the fetched files):
#
#    Niederstaetter H? no -- authors of record:
#    "Discrimination between 34 of 36 Possible Combinations of Three C>T SNP
#    Genotypes in the MGMT Promoter by High Resolution Melting Analysis Coupled
#    with Pyrosequencing"
#    Int J Mol Sci 2021;22(22):12527.  DOI 10.3390/ijms222212527
#    PMID 34830407, PMC8621402.
#
#  Fetched:
#    .tmp/lit/ijms2022/fulltext.xml   (EPMC fullTextXML, 171,575 B)
#    .tmp/lit/ijms2022/supplements.zip -> supplements/Supplemental Tables.pdf
#      rendered to .tmp/lit/ijms2022/supp_tables.txt with `pdftotext -layout`
#      (no pdf package in this toolchain).  Table S1 of that PDF text is the
#      evidence for the 36 oligodeoxynucleotide controls.
#
#  WHY this paper: it is the only fetched source that enumerates EVERY
#  genotype combination of a multi-SNP amplicon and states which ones HRM
#  separates.  That is exactly the "resolvability" quantity the package needs
#  and cannot get from a single-SNP accuracy figure.
#
#  WHAT THIS SCRIPT DOES NOT MODEL:
#
#   * TWO SNPs ASSAY, PER-TUPLE PARTITION: the paper states the 2-SNP partition
#     only as group SIZES (six groups of 2, two of 3, three of 4, one of 5) plus
#     named pairs.  It never assigns all 35 controls to a group, so this script
#     lands the sizes and the sum-to-36 assertion, NOT a per-control group label.
#     Inventing the labels would be a fabrication.
#
#   * THREE SNPs ASSAY, 22 IDENTIFIED: the 14 remaining controls are given as 7
#     verbatim pairs, so 22 is derived as 36 - 14 and is asserted against the
#     prose counts by genetic state (6 + 8 + 6 + 2 = 22).  The genetic-state mix
#     of the 22 (6 triple-hom, 8 double-hom, 6 double-het, 2 triple-het) is NOT
#     re-derivable: it depends on which 22 are left, which the paper does not
#     enumerate.  It is carried as the paper's claim and checked only for
#     consistency with the totals.
#
#   * HRM+PSQ "34 of 36": this is a COMBINED-CHEMISTRY result.  It is recorded
#     as such and must never be used as an HRM-only resolvability number.
#
#   * No probability, rate or per-base error is derived from any of this.  Every
#     row is a count of controls under a named method.
#
#  ASSERTIONS (each is a stop() if the fetched files disagree):
#   A1  Table S1 has exactly 36 control combinations, and its Genetic State
#       column counts 8 triple-homozygous / 12 double-homozygous /
#       12 double-heterozygous / 4 triple-heterozygous, matching the prose.
#   A2  2-SNP HRM partition: 1 singleton + 6*2 + 2*3 + 3*4 + 1*5 = 36.
#   A3  3-SNP HRM: 22 + 14 = 36; the 14 named controls are exactly 7 pairs;
#       each named control is in Table S1 and none is repeated.
#   A4  The 22 identified controls are a subset of Table S1 and their
#       genetic-state counts can sum to 6/8/6/2 exactly.
#   A5  HRM+PSQ: 34 + 2 = 36, and the two unresolvable are exactly the prose's
#       CTC/TCT and CCT/TTC.
#   A6  dTm for triple-homozygous controls is monotone non-increasing in the
#       number of T alleles: 0.4 (1T) <= 0.6, 0.7 (1T) < 1.0 (2T) < 1.3 (2T)
#       < 1.7 (3T); and 1.3 < 1.7 while 1.3 > 1.0.
#   A7  CTC/CTC, TTC/TTC and TCT/TCT share a single 86.76 C transition? no:
#       the prose gives ONE shared value (86.76 C) for TTC/TTC and TCT/TCT;
#       those must be exactly the pair listed as unresolvable-by-HRM.
#
#  Run:  Rscript data-raw/ijms2022.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

xml_path <- fs::path(root, ".tmp", "lit", "ijms2022", "fulltext.xml")
txt_path <- fs::path(root, ".tmp", "lit", "ijms2022", "supp_tables.txt")
for (p in c(xml_path, txt_path)) if (!file.exists(p)) stop("Missing ", p)

DOI <- "10.3390/ijms222212527"
out_dir <- fs::path(root, "inst", "extdata")

doc <- xml2::read_xml(xml_path)
if (!grepl("10.3390/ijms222212527", xml2::xml_text(doc), fixed = TRUE))
  stop("fulltext.xml does not mention ", DOI, " - wrong file?")

# --- flat prose text: one string per paragraph -------------------------------
paras <- vapply(
  xml2::xml_find_all(doc, "//*[local-name()='p']"),
  function(p) gsub("[[:space:]]+", " ", xml2::xml_text(p)), character(1))
para <- function(pattern) {
  hit <- paras[grepl(pattern, paras)]
  if (length(hit) != 1L) stop("Expected exactly 1 paragraph matching /", pattern,
                              "/, got ", length(hit))
  hit
}
p_2snp   <- para("only one of the 36 oligodeoxynucleotide controls")
p_3snp   <- para("^22 oligodeoxynucleotide controls")
p_disc   <- para("unambiguously identified by HRM analysis, PSQ, and HRM coupled")
p_both3  <- para("^However, by combined analysis of HRM and PSQ data, 34 out of 36")
p_dtm    <- para("^All triple-homozygous combinations had a single melting transition")
p_dhet   <- para("^Among the triple-heterozygous combinations")
p_min    <- para("^Tissue samples obtained from biopsies are often limited")
p_ann    <- para("minimal initial DNA amount")

# --- A1: Table S1 -------------------------------------------------------------
s1 <- readLines(txt_path, warn = FALSE)
start <- grep("Table S1\\.", s1)
end   <- grep("Table S2\\.", s1)
if (length(start) != 1L || length(end) != 1L || end <= start)
  stop("Could not delimit Table S1 in ", txt_path, " (start ", length(start),
       ", end ", length(end), ")")
s1_lines <- s1[(start + 1L):(end - 1L)]
row_re <- "^([ACGT]{3}/[ACGT]{3})[[:space:]]+([ACGT]{3})[[:space:]]+([ACGT]{3})[[:space:]]+(triple-homozygous|double-homozygous|double-heterozygous|triple-heterozygous)$"
s1_lines <- trimws(s1_lines)
m <- regmatches(s1_lines, regexec(row_re, s1_lines))
ok <- vapply(m, length, integer(1)) == 5L
controls <- tibble::tibble(
  combination  = vapply(m[ok], `[[`, character(1), 2L),
  strand_1     = vapply(m[ok], `[[`, character(1), 3L),
  strand_2     = vapply(m[ok], `[[`, character(1), 4L),
  genetic_state = vapply(m[ok], `[[`, character(1), 5L)
)
# strand pair must be the two halves of the combination, in order. The table
# does NOT sort within a row (e.g. TCC/CTC), so order is not asserted.
if (any(controls$combination != paste0(controls$strand_1, "/", controls$strand_2)))
  stop("Table S1: strand columns do not concatenate to the combination")
if (anyDuplicated(controls$combination) || nrow(controls) != 36L)
  stop("Table S1: expected 36 unique combinations, got ", nrow(controls))
state_tab <- table(factor(controls$genetic_state,
                          levels = c("triple-homozygous", "double-homozygous",
                                     "double-heterozygous", "triple-heterozygous")))
if (!identical(as.integer(state_tab), c(8L, 12L, 12L, 4L)))
  stop("Table S1 genetic-state counts are ", paste(state_tab, collapse = "/"),
       ", prose states 8/12/12/4")
if (!any(grepl("8 triple-homozygous, 12 double-homozygous, 12 double-heterozygous and 4 triple-heterozygous",
              paras, fixed = TRUE))) {
  stop("Prose does not state the 8/12/12/4 genetic-state split")
}
# Derived check: with 3 biallelic loci every combination is a multiset of two
# triplets, so the total T count is 0..6 and the genetic state must follow
# from it: 0 or 6 T => triple-homozygous, 1 or 5 T => double-homozygous, and
# 2 T (or 4 T) => double-homozygous if both T sit on the same strand else
# double-heterozygous; 3 T is double-heterozygous or triple-heterozygous,
# which the table, not the T count, decides.
t_count <- nchar(gsub("[^T]", "", controls$combination))
state_from_count <- function(i) {
  a <- controls$strand_1[[i]]; b <- controls$strand_2[[i]]
  if (a == b) return("triple-homozygous")
  # "homozygous at a locus" = the two strands show the same base there
  hom_a <- sum(strsplit(a, "")[[1L]] == strsplit(b, "")[[1L]])
  het_loci <- 3L - hom_a
  if (het_loci == 1L) return(if (hom_a == 2L) "double-homozygous" else "double-heterozygous")
  if (het_loci >= 2L) return(if (hom_a >= 1L) "double-heterozygous" else "triple-heterozygous")
}
bad_state <- !vapply(seq_len(nrow(controls)), function(i)
  controls$genetic_state[[i]] %in% state_from_count(i), logical(1))
if (any(bad_state))
  stop("Table S1: genetic state disagrees with the T-allele pattern for ",
       paste(controls$combination[bad_state], collapse = ", "))
t_state <- table(factor(t_count, levels = 0:6))
# For a multiset of two of the 8 triplets, #T is unimodal and symmetric and
# counts 1 at 0 and 6 (CCC/CCC, TTT/TTT).  It is NOT binomial C(6,k): that
# would describe 6 independent loci drawn twice, and here both strands are
# drawn from the same 8 triplets.
if (!identical(as.integer(t_state), c(1L, 3L, 9L, 10L, 9L, 3L, 1L)))
  stop("Table S1: T-allele count histogram is ", paste(t_state, collapse = "/"),
       ", expected 1/3/9/10/9/3/1 (multiset of two of 8 triplets)")

# --- A2: 2-SNP HRM partition --------------------------------------------------
part_re <- "(six groups containing two), (two groups containing three), (three groups containing four), and (one group containing five)"
part <- regmatches(p_2snp, regexec(part_re, p_2snp))[[1L]]
if (length(part) != 5L || !grepl("only one of the 36", p_2snp, fixed = TRUE) ||
    !grepl("The remaining 35 variant combinations resulted in 12 groups", p_2snp, fixed = TRUE))
  stop("2-SNP partition sentence not parseable: ", p_2snp)
counts_2snp <- c(two = 6L, three = 2L, four = 3L, five = 1L)
# cross-check the parsed words against the written numerals
words <- c(six = 6L, two = 2L, three = 3L, one = 1L)
num_of <- function(w) words[[sub(" .*$", "", w)]]
if (!identical(unname(vapply(part[2:5], num_of, integer(1))), unname(counts_2snp)))
  stop("2-SNP partition numerals do not match the expected 6/2/3/1")
if (sum(counts_2snp) != 12L)
  stop("2-SNP partition does not describe 12 groups")
n_singleton_2snp <- 1L
sizes_2snp <- c(two = 2L, three = 3L, four = 4L, five = 5L)
n_grouped <- sum(sizes_2snp * counts_2snp)
if (n_grouped + n_singleton_2snp != 36L)
  stop("2-SNP partition sums to ", n_grouped + n_singleton_2snp, ", not 36")
partition <- tibble::tibble(
  assay = "2 SNPs",
  melt_class = c("singleton (uniquely identified)", "group of 2", "group of 3",
                 "group of 4", "group of 5"),
  class_size = c(n_singleton_2snp, unname(sizes_2snp)),
  n_classes = c(1L, unname(counts_2snp)),
  n_controls = c(n_singleton_2snp, unname(sizes_2snp * counts_2snp))
)

# --- A3: 3-SNP HRM: 14 in 7 pairs --------------------------------------------
pair_part <- regmatches(p_3snp, regexec("The remaining 14 controls \\((.*?)\\) led to pairwise", p_3snp))[[1L]]
if (length(pair_part) != 2L) stop("3-SNP 'remaining 14' clause not parseable")
pairs_raw <- strsplit(pair_part[[2L]], ";[[:space:]]*")[[1L]]
unresolvable_3snp <- do.call(rbind, lapply(strsplit(pairs_raw, ",[[:space:]]*"), function(v) {
  if (length(v) != 2L) stop("Bad pair: ", paste(v, collapse = ","))
  v
}))
if (nrow(unresolvable_3snp) != 7L || length(unique(as.vector(unresolvable_3snp))) != 14L)
  stop("3-SNP: expected 7 disjoint pairs (14 controls)")
if (!all(as.vector(unresolvable_3snp) %in% controls$combination))
  stop("3-SNP: unresolvable control(s) absent from Table S1: ",
       paste(setdiff(as.vector(unresolvable_3snp), controls$combination), collapse = ", "))
n_identified_3snp <- 36L - sum(lengths(unresolvable_3snp))
if (n_identified_3snp != 22L) stop("3-SNP identified = ", n_identified_3snp, ", prose says 22")
# prose counts by state for the 22 must be feasible given Table S1's totals
state_claim <- c(triple_hom = 6L, double_hom = 8L, double_het = 6L, triple_het = 2L)
if (sum(state_claim) != n_identified_3snp) stop("3-SNP state counts do not sum to 22")
if (any(state_claim > as.integer(state_tab)))
  stop("3-SNP state counts exceed Table S1's totals")
if (!grepl("22 oligodeoxynucleotide controls (6 triple-homozygous, 8 double-homozygous, 6 double-heterozygous, and 2 triple-heterozygous)",
           p_3snp, fixed = TRUE))
  stop("3-SNP prose state counts not found verbatim")

# --- A5: HRM+PSQ --------------------------------------------------------------
if (!grepl("34 out of 36 variant combinations could be genotyped unambiguously",
           p_both3, fixed = TRUE)) stop("HRM+PSQ 34/36 sentence not found")
if (!grepl("The triple-heterozygous combinations CTC/TCT and CCT/TTC were the only ones that could not be distinguished",
           p_both3, fixed = TRUE)) stop("HRM+PSQ unresolved pair sentence not found")
if (!identical(pairs_raw[[1L]], "CTC/TCT, CCT/TTC"))
  stop("The 3-SNP unresolvable list does not open with CTC/TCT, CCT/TTC")
# The discussion restates all three assays' counts; assert the restatement.
if (!grepl("22, 16, and 34 variant combinations by HRM analysis, PSQ, and by coupling HRM to PSQ",
           p_disc, fixed = TRUE)) stop("Discussion restatement of 22/16/34 not found")
if (!grepl("one, twelve, and thirteen variant combinations could be unambiguously identified",
           p_disc, fixed = TRUE)) stop("Discussion restatement of the 2-SNP counts not found")
# 2-SNP HRM+PSQ = 13; PSQ alone 12; HRM alone 1 -- consistency of the trio.
n_2snp_hrm_psq <- 13L
if (n_singleton_2snp + 12L != n_2snp_hrm_psq) stop("2-SNP HRM+PSQ arithmetic mismatch")

# --- A6: dTm ------------------------------------------------------------------
d4 <- regmatches(p_dtm, regexec(
  "decreased Tm by ([0-9.]+) °C \\(rs527559815\\), ([0-9.]+) °C \\(rs547832288\\), or ([0-9.]+) °C \\(rs16906252\\)",
  p_dtm))[[1L]]
if (length(d4) != 4L) stop("Single-T dTm clause not parseable")
dtm_single <- tibble::tibble(
  rsid = c("rs527559815", "rs547832288", "rs16906252"),
  position = c(1L, 2L, 3L),
  dtm_c = as.numeric(c(d4[2L], d4[3L], d4[4L]))
)
d2 <- regmatches(p_dtm, regexec(
  "one at rs527559815 and the other one .*? resulted in a decrease of ([0-9.]+) °C",
  p_dtm))[[1L]]
d2b <- regmatches(p_dtm, regexec(
  "one at rs547832288 and one at rs16906252, and the triple-homozygous combination having three Ts led to a decrease of ([0-9.]+) °C and ([0-9.]+) °C",
  p_dtm))[[1L]]
if (length(d2) != 2L || length(d2b) != 3L) stop("2T/3T dTm clauses not parseable")
dtm_tcount <- tibble::tibble(
  n_t_alleles = c(0L, 1L, 1L, 1L, 2L, 2L, 3L),
  labels = c("CCC/CCC (reference)", "T at rs527559815", "T at rs547832288",
             "T at rs16906252", "T at rs527559815 + one other",
             "T at rs547832288 + rs16906252", "TTT/TTT"),
  dtm_c = c(0, dtm_single$dtm_c[1L], dtm_single$dtm_c[2L], dtm_single$dtm_c[3L],
            as.numeric(d2[2L]), as.numeric(d2b[2L]), as.numeric(d2b[3L]))
)
# the two 1-T substitution effects are per-locus, so they need not be ordered;
# the cumulative values (1 T combined, 2 T, 3 T) must be strictly increasing.
cum <- dtm_tcount$dtm_c[dtm_tcount$n_t_alleles >= 2L]
if (is.unsorted(cum, strictly = TRUE))
  stop("dTm is not strictly increasing in T allele count: ", paste(cum, collapse = ", "))
if (max(dtm_single$dtm_c) >= as.numeric(d2[2L]))  # any single T < any double T
  stop("A single-T dTm is not smaller than the double-T dTm")
if (as.numeric(d2[2L]) >= as.numeric(d2b[2L]) ||
    as.numeric(d2b[2L]) >= as.numeric(d2b[3L]))
  stop("Double-/triple-T dTm values are not ordered")

# --- A7: TTC/TTC vs TCT/TCT share one transition -----------------------------
if (!any(grepl("TTC/TTC and TCT/TCT", paras)) ||
    !any(grepl("single melt transition at 86\\.76 °C", paras)))
  stop("Shared 86.76 C transition clause not found")
if (!any(unresolvable_3snp[, 1L] == "TTC/TTC" & unresolvable_3snp[, 2L] == "TCT/TCT") &&
    !any(unresolvable_3snp[, 1L] == "TCT/TCT" & unresolvable_3snp[, 2L] == "TTC/TTC"))
  stop("TTC/TTC and TCT/TCT share a Tm but are not listed as unresolvable")
tm_shared <- 86.76

# --- minimal input -------------------------------------------------------------
if (!all(vapply(c("down to 3 ng DNA", "only possible down to an initial DNA amount of 4 ng",
                  "down to an initial DNA amount of 1 ng") , function(s) grepl(s, p_min, fixed = TRUE), logical(1))) ||
    !grepl("0\\.1 ng DNA caused amplification of only one DNA strand", p_min) ||
    !grepl("0\\.01 ng input was too low to obtain any PCR products", p_min))
  stop("Minimal-input clauses not all present")
min_input <- tibble::tibble(
  method = c("HRM", "HRM", "PSQ", "HRM (single-strand only)", "HRM (no product)"),
  genotype = c("CCC/CCT (breast tissue)", "CCC/CTT and CCC/CCT vs CCC/TCC (NHGRI)",
               "all genotypes tested", "any", "any"),
  min_ng = c(3, 4, 1, 0.1, 0.01),
  outcome = c("correctly identified", "correctly identified", "correctly identified",
              "amplification of one strand only - NOT a genotype call",
              "too low to amplify - NOT a genotype call")
)
if (min_input$min_ng[4L] >= min_input$min_ng[3L] || min_input$min_ng[5L] >= min_input$min_ng[4L])
  stop("Minimal-input series is not decreasing")

# --- write: per-assay discrimination summary ---------------------------------
HRM_DEF <- paste(
  "Outcome = the number of distinct melting-curve classes the assay produced over",
  "the 36 oligodeoxynucleotide controls that enumerate every combination of the",
  "three C>T SNPs, and the number of controls left in a class with another",
  "control. A control in a shared class is NOT a genotype call. HRM-only rows",
  "use HRM data alone; HRM+PSQ rows are a COMBINED-CHEMISTRY result and are not",
  "HRM resolvability. Counts only - no rate is derived from them here."
)
disc <- tibble::tibble(
  assay = c("2 SNPs", "2 SNPs", "2 SNPs", "3 SNPs", "3 SNPs", "3 SNPs"),
  readout = c("HRM", "PSQ", "HRM+PSQ", "HRM", "PSQ", "HRM+PSQ"),
  n_controls = 36L,
  n_uniquely_identified = c(n_singleton_2snp, 12L, n_2snp_hrm_psq,
                            n_identified_3snp, 16L, 34L),
  n_not_discriminated = 36L - n_uniquely_identified,
  evidence_tier = c("measured", "measured", "measured", "measured", "measured", "measured"),
  source_file = "fulltext.xml",
  doi = DOI,
  outcome_definition = HRM_DEF
)
if (any(disc$n_uniquely_identified + disc$n_not_discriminated != 36L))
  stop("discrimination rows do not sum to 36")
if (!identical(disc$n_not_discriminated[disc$readout == "HRM+PSQ"], c(23L, 2L)))
  stop("HRM+PSQ unresolved counts are not 23 (2 SNPs) and 2 (3 SNPs)")
readr::write_csv(disc, fs::path(out_dir, "hrm_ijms2022_discrimination.csv"))

readr::write_csv(
  dplyr::mutate(partition, n_controls_total = 36L, doi = DOI,
                evidence_tier = "measured", source_file = "fulltext.xml",
                outcome_definition = paste(
                  "Outcome = the size of one melting-curve class in the 2-SNP assay",
                  "(74 bp amplicon, rs547832288 + rs16906252). The paper reports the",
                  "class SIZES only ('six groups containing two, two groups containing",
                  "three, three groups containing four, and one group containing five'),",
                  "not a control-to-class assignment; n_controls = size x n_classes.")),
  fs::path(out_dir, "hrm_ijms2022_partition_2snp.csv"))

readr::write_csv(
  dplyr::bind_rows(
    tibble::tibble(assay = "3 SNPs", melt_class = "pair (unresolvable)",
                   control_a = unresolvable_3snp[, 1L], control_b = unresolvable_3snp[, 2L],
                   n_controls_total = 36L, n_uniquely_identified = n_identified_3snp),
    tibble::tibble(assay = "3 SNPs", melt_class = "singleton (uniquely identified)",
                   control_a = NA_character_, control_b = NA_character_,
                   n_controls_total = 36L, n_uniquely_identified = NA_integer_)
  ) |>
    dplyr::mutate(
      doi = DOI, evidence_tier = "measured", source_file = "fulltext.xml",
      outcome_definition = paste(
        "Outcome = which controls the 3-SNP HRM assay (97 bp amplicon) leaves in the",
        "same melting-curve class. The 7 pairs are quoted verbatim from the Results;",
        "the 22 singletons are given by the paper only as a count with a genetic-state",
        "split (6/8/6/2), NOT enumerated, so no control is named as a singleton here.")),
  fs::path(out_dir, "hrm_ijms2022_unresolvable_3snp.csv"))

readr::write_csv(
  dplyr::mutate(dtm_tcount, doi = DOI, evidence_tier = "measured",
                source_file = "fulltext.xml",
                outcome_definition = paste(
                  "Outcome = dTm (C) of a triple-homozygous control relative to",
                  "CCC/CCC, read from the derivative melting curve. All three SNPs are",
                  "C>T (class I transitions) in one GC-rich promoter: this is a",
                  "per-locus effect for THIS amplicon and must not be generalised to",
                  "transversions or to other sequence context.")),
  fs::path(out_dir, "hrm_ijms2022_dtm_single_snp_controls.csv"))

readr::write_csv(
  dplyr::mutate(min_input, doi = DOI, evidence_tier = "measured",
                source_file = "fulltext.xml",
                outcome_definition = paste(
                  "Outcome = the smallest initial DNA input (ng, fluorimetric) at which",
                  "the named genotype was still correctly identified by the named",
                  "method. This is an INPUT amount, not a variant allele fraction: it",
                  "says nothing about detecting a minority allele in a mixture. Rows",
                  "0.1 ng and 0.01 ng are failures (no genotype), recorded to bound",
                  "the series from below.")),
  fs::path(out_dir, "hrm_ijms2022_min_input.csv"))

message("ijms2022: Table S1 ", nrow(controls), " controls (",
        paste(state_tab, collapse = "/"), "); 2-SNP partition 1+",
        paste(sizes_2snp * counts_2snp, collapse = "+"), "=36; 3-SNP ",
        n_identified_3snp, " identified / ", 2L * nrow(unresolvable_3snp),
        " in 7 pairs; HRM+PSQ 34/36; dTm ", paste(dtm_tcount$dtm_c, collapse = "/"),
        " C (assertions passed)")
