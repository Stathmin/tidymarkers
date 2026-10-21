# =============================================================================
#  data-raw/sharplab2021.R  --  the Chinese maize RNA-seq KASP set (Chen et al.)
# =============================================================================
#  Chen Z, Tang D, Ni J, Li P, Wang L, Zhou J, Li C, Lan H, Li L, Liu J.
#  "Development of genic KASP SNP markers from RNA-Seq data for map-based
#  cloning and marker-assisted selection in maize."
#  BMC Plant Biol 2021;21:157.  DOI 10.1186/s12870-021-02932-8, PMC8004444.
#
#  This source is the counter-example to the usual conflation.  It reports THREE
#  different denominators for the same 50 markers:
#
#    design       50 of 50 got a synthesised primer triple (Additional file 4
#                 lists all 50 with FA/FD/R primers)              -> 100%
#    amplification  "The success ratio of PCR amplification for 50 KASP SNP
#                 markers in the 20 inbred lines was 98.5%"      -> a RATIO, no
#                 count is printed, so no count is landed
#    polymorphism  46 of 50 "showed polymorphism between the maize inbred lines
#                 1212 and B73 using PCR"                          -> 92%
#
#  The 92% is NOT a design rate and NOT an amplification rate: it is the
#  fraction of the 50 designed assays whose PCR product was polymorphic between
#  two specific inbreds.  Each column of the extract says which one it is.
#
#  The brief for this worker expected a "cross-crop conversion survey" in this
#  paper.  It is not there: the paper compares its average PIC (0.13) with
#  sesame (0.38), radish (0.32), grapevine (0.25) and sweet orange (0.2) --
#  a PIC survey, not a KASP conversion survey.  Nothing cross-crop is landed.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)
dir.create(file.path(root, ".tmp", "lit", "sharplab2021"), recursive = TRUE,
           showWarnings = FALSE)
source(file.path(root, "data-raw", "kasp_lib.R"))

DATE_VERIFIED <- "2026-10-21"
DOI <- "10.1186/s12870-021-02932-8"
OUT <- fs::path(root, "inst", "extdata")

xlsx <- kasp_input(
  root, "sharplab2021", "supplements/12870_2021_2932_MOESM4_ESM.xlsx",
  "sharplab2021_validation.xlsx", "sharplab2021 Additional file 4"
)
xml_path <- kasp_input(
  root, "sharplab2021", "fulltext.xml",
  "sharplab2021_fulltext.xml", "sharplab2021 full text"
)

## --- Additional file 4: the 50 validated markers -----------------------------
stopifnot(file.size(xlsx) > 10000L)
sheets <- kasp_sheets(xlsx)
sh <- if (!is.null(sheets[["3"]])) sheets[["3"]] else sheets[[1L]]
rows <- lapply(sh, kasp_row)

# header: row 8 names the columns, row 9 carries the primer sub-headers
h8 <- rows[[8L]]
if (!identical(unname(h8[c("B", "C", "D", "E", "F", "G", "H")]),
               c("Primer_Name", "Chromosome", "Locus", "SNP", "PIC(in 368 RNA-seq)",
                 "polymorphism between 1212 and B73", "PIC(in 20 maize inbred lines)"))) {
  stop("Additional file 4's header is not where it was: ", paste(h8[1:8], collapse = " | "))
}
h9 <- rows[[9L]]
if (!identical(unname(h9[c("I", "J", "K")]), c("FA", "FD", "R"))) {
  stop("Additional file 4's primer sub-header (FA/FD/R) is not where it was")
}

# a data row: numbered 1..50 in A, an id like Chr1-4122559 in B, a primer in I
dr <- Filter(function(r) {
  grepl("^[0-9]+$", r[["A"]]) && nzchar(r[["B"]]) && grepl("^Chr", r[["B"]]) &&
    nzchar(r[["K"]])
}, rows)
markers <- do.call(rbind, lapply(dr, function(r) {
  data.frame(
    no = as.integer(r[["A"]]), marker = r[["B"]], chromosome = r[["C"]],
    locus = r[["D"]], snp = r[["E"]],
    pic_rnaseq = suppressWarnings(as.numeric(r[["F"]])),
    polymorphism_1212_vs_b73 = r[["G"]],
    pic_20_lines = suppressWarnings(as.numeric(r[["H"]])),
    primer_fa = r[["I"]], primer_fd = r[["J"]], primer_r = r[["K"]],
    stringsAsFactors = FALSE
  )
}))

if (nrow(markers) != 50L || !identical(sort(markers$no), 1:50)) {
  stop("Additional file 4 does not hold markers numbered 1..50 (parsed ",
       nrow(markers), " rows)")
}
if (anyDuplicated(markers$marker)) stop("the marker ids are not unique")

poly <- markers$polymorphism_1212_vs_b73
if (!all(poly %in% c("\u221a", "\u00d7"))) {
  stop("the polymorphism column is not the tick/cross vocabulary: ",
       paste(unique(poly), collapse = " "))
}
n_poly <- sum(poly == "\u221a")
n_nopoly <- sum(poly == "\u00d7")
# the crosses have no PIC in the 20 lines, which is the authors' own marker of
# "not polymorphic"; assert that rather than tolerate a silent NA
has_no_pic <- is.na(markers$pic_20_lines)
if (!identical(has_no_pic, poly == "\u00d7")) {
  stop("the '--' PIC cells do not coincide with the crosses: ",
       paste(markers$no[has_no_pic != (poly == "\u00d7")], collapse = ", "))
}

## --- the full text: the three denominators ----------------------------------
stopifnot(file.size(xml_path) > 20000L)
txt <- paste(readLines(xml_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
if (!grepl("<article", txt, fixed = TRUE)) stop("not a JATS article")
flat <- kasp_norm(gsub("<[^>]+>", " ", txt))

# 46 (92.00%) out of 50 ... "showed polymorphism between the maize inbred lines
# 1212 and B73 using PCR amplification"
a46 <- kasp_find(flat, "out of 50 KASP SNP markers")
# "... Following validation, 46 (92.00%) out of 50 KASP SNP markers showed
# polymorphism ...".  The engine indexes inconsistently here, so the sentence
# is cut by its two stable ends and the numbers are read from the cut text.
pre46 <- substr(flat, a46$start - 60L, a46$start + 60L)
if (!grepl("Following validation, 46 [(]92.00%[)]", pre46, perl = TRUE)) {
  stop("the 92.00% sentence has changed; the text at its site is: ", pre46)
}
num46 <- as.integer(sub("^.*validation, ([0-9]+) [(].*$", "\\1", pre46))
if (is.na(num46)) stop("cannot read the 46 from: ", pre46)
if (!grepl("(92.00%)", pre46, fixed = TRUE)) {
  stop("the printed rate after 'out of 50 KASP SNP markers' is no longer 92.00%")
}
if (num46 != n_poly) {
  stop("the paper says ", num46, " of 50 are polymorphic but Additional file 4 ",
       "has ", n_poly, " ticks")
}
if (!grepl("92.00", flat, fixed = TRUE)) {
  stop("the paper no longer prints 92.00%; the rate being claimed has changed")
}
if (n_poly / 50 != 0.92) stop("46/50 is not 0.92")

# "The success ratio of PCR amplification for 50 KASP SNP markers in the 20
# inbred lines was 98.5%" -- a ratio with no count, so no count is landed
a985 <- kasp_find(flat, "success ratio of PCR amplification")
site985 <- kasp_site(a985, flat, 150L)
if (!grepl("in the 20 inbred lines was 98.5%", site985, fixed = TRUE)) {
  stop("the 98.5% sentence is not where it was: ", site985)
}
pct985 <- 98.5

# "all 46 SNPs showing polymorphism in the 20 maize inbred lines"
if (!grepl("all 46 SNPs showing polymorphism", flat, fixed = TRUE)) {
  stop("the 'all 46 SNPs ... in the 20 maize inbred lines' sentence is absent")
}
n_in_20 <- 46L

# in-silico conversion for the whole genome
a71 <- kasp_find(flat, "were identified in 2,948,985")
site71 <- substr(flat, a71$start - 30L, a71$start + 50L)
if (!grepl("71,311", site71, fixed = TRUE)) {
  stop("the 71,311 in-silico KASP count is not where it was: ", site71)
}
n_insilico <- 71311L
n_snp <- 2948985L
if (!grepl("2,948,985 SNP", flat, fixed = TRUE)) {
  stop("the 2,948,985 SNP denominator is not where it was")
}
if (!grepl("2.42%", flat, fixed = TRUE)) {
  stop("the paper no longer prints 2.42% for the in-silico conversion")
}
rate_insilico <- round(n_insilico / n_snp, 4)
if (rate_insilico != 0.0242) {
  stop("71,311 / 2,948,985 = ", rate_insilico, ", not the printed 2.42%")
}

# PIC ranges: the abstract says 0.14-0.5 in 368 lines, the results say
# 0.10-0.50 in 20 lines.  Both are checked against the file, and the file wins
# where they conflict (its minimum is recorded verbatim).
pic_rnaseq <- range(markers$pic_rnaseq, na.rm = TRUE)
pic_20 <- range(markers$pic_20_lines, na.rm = TRUE)
if (abs(pic_rnaseq[1L] - 0.14) > 1e-9 || abs(pic_rnaseq[2L] - 0.50) > 1e-9) {
  stop("the 368-line PIC range in the file is ", paste(pic_rnaseq, collapse = "-"),
       ", not the printed 0.14-0.50")
}
if (abs(pic_20[1L] - 0.10) > 1e-9) {
  stop("the 20-line PIC minimum in the file is ", pic_20[1L],
       ", not the printed 0.10")
}
pic_discrepancy <- paste0(
  "PAPER: the abstract says the 20-line PIC ranged 0.11-0.50; the results say ",
  "0.10-0.50; Additional file 4's minimum is ", format(pic_20[1L], nsmall = 2L),
  ". The file is recorded and the abstract's 0.11 is not used."
)

## --- write: one row per marker ------------------------------------------------
markers$design_success <- TRUE   # a primer triple exists in Additional file 4
markers$assay_worked <- poly == "\u221a"
markers$outcome_polymorphism <- ifelse(poly == "\u221a",
                                       "polymorphic_1212_vs_b73",
                                       "not_polymorphic_1212_vs_b73")
markers$pic_20_printed_cell <- ifelse(is.na(markers$pic_20_lines), "--",
                                      format(markers$pic_20_lines, nsmall = 2L))
markers$doi <- DOI
markers$evidence_tier <- "measured"
markers$source_file <- basename(xlsx)
markers$source_table <- "Additional file 4 (sheet '3'), rows 10-59"
markers$date_verified <- DATE_VERIFIED
markers$conversion_definition <- paste(
  "DESIGN: 50/50 = 1.00 -- every one of the 50 selected markers carries a",
  "synthesised FA/FD/R primer triple in Additional file 4. This is the design",
  "rate and ONLY the design rate."
)
markers$outcome_definition <- paste(
  "WORK (paper's own words): '46 (92.00%) out of 50 KASP SNP markers showed",
  "polymorphism between the maize inbred lines 1212 and B73 using PCR",
  "amplification'. So `assay_worked` = the PCR product was polymorphic between",
  "ONE pair of inbreds, 1212 and B73 -- it is not 'the assay amplified' (the",
  "paper reports a separate 98.5% amplification ratio) and it is not 'the assay",
  "genotyped the 20-line panel' (all 46 polymorphic markers did). A marker that",
  "amplified but was monomorphic between 1212 and B73 is counted here as not",
  "worked, which is the paper's own classification."
)
markers <- markers[c("no", "marker", "chromosome", "locus", "snp", "pic_rnaseq",
                     "polymorphism_1212_vs_b73", "pic_20_printed_cell",
                     "pic_20_lines", "primer_fa", "primer_fd", "primer_r",
                     "design_success", "assay_worked", "outcome_polymorphism",
                     "conversion_definition", "outcome_definition", "doi",
                     "evidence_tier", "source_file", "source_table",
                     "date_verified")]
kasp_write_checked(markers, fs::path(OUT, "sharplab2021_assays.csv"))

## --- write: the study-level denominators --------------------------------------
def_design <- "DESIGN (a synthesised primer triple exists): 50/50 from Additional file 4."
def_poly <- paste(
  "POLYMORPHISM between two specific inbreds by PCR: 46/50 = 92.00% (the",
  "paper's number, recomputed from the tick column of Additional file 4). Not",
  "a design rate and not an amplification rate."
)
def_amp <- paste(
  "AMPLIFICATION ratio stated for the 20-inbred panel: 98.5%. The paper prints",
  "a RATIO and no count, so this extract records the printed percentage only",
  "and derives no count from it."
)
summary_out <- data.frame(
  quantity = c("design_success", "polymorphic_1212_vs_b73", "amplified_20_lines",
               "polymorphic_in_20_lines", "insilico_kasp_markers",
               "rna_seq_snps", "insilico_conversion"),
  n_numerator = c(50L, n_poly, NA_integer_, n_in_20, n_insilico, n_snp, n_insilico),
  n_denominator = c(50L, 50L, NA_integer_, n_poly, n_snp, NA_integer_, n_snp),
  value_printed = c("50/50 (Additional file 4)", "46 (92.00%)", "98.5%",
                    "all 46", "71,311", "2,948,985", "2.42%"),
  outcome_definition = c(def_design, def_poly, def_amp,
    paste("POLYMORPHISM in the 20-line panel: all 46 polymorphic markers also",
          "polymorphic in the 20 lines. Not a conversion rate."),
    paste("IN-SILICO (no wet lab): markers passing the authors' sequence",
          "criteria (unique region, bi-allelic, PIC >= 0.4, conserved primer",
          "sequence) out of all RNA-seq SNPs. This is a sequence-criterion",
          "yield, not an assay outcome and not a wet-lab conversion rate."),
    "The RNA-seq SNP call set this source starts from.",
    paste("IN-SILICO conversion: 71,311 / 2,948,985 = 2.42%, the paper's own",
          "number. Comparable to nothing in conversion_rates.csv.")),
  note = c("every marker in the file",
           paste("the 92% denominator: 50 markers that received a primer triple",
                 "and were tested on 1212 x B73"),
           "printed as a ratio, no count",
           "the paper's phrase is 'all 46 SNPs'",
           "printed as '71,311 SNP markers were identified'",
           "printed as 2,948,985 SNP loci",
           "recomputed from the two printed counts"),
  pic_note = c(NA_character_, NA_character_, NA_character_, NA_character_,
               NA_character_, NA_character_, pic_discrepancy),
  doi = DOI, evidence_tier = "measured", source_file = basename(xml_path),
  source_table = "Results / Methods text",
  date_verified = DATE_VERIFIED,
  stringsAsFactors = FALSE
)
summary_out <- summary_out[c("quantity", "n_numerator", "n_denominator",
                             "value_printed", "outcome_definition", "note",
                             "pic_note", "doi", "evidence_tier", "source_file",
                             "source_table", "date_verified")]
kasp_write_checked(summary_out, fs::path(OUT, "sharplab2021_outcomes.csv"))

# land the primer triples too: they are the Y3.2 artifact for maize and are
# identical in structure to makhoul_kasp_primers.csv
prim <- markers[c("marker", "primer_fa", "primer_fd", "primer_r", "snp")]
prim$as_3p_differs <- substr(prim$primer_fa, nchar(prim$primer_fa), nchar(prim$primer_fa)) !=
  substr(prim$primer_fd, nchar(prim$primer_fd), nchar(prim$primer_fd))
prim$alleles_on_primers <- paste0(
  substr(prim$primer_fa, nchar(prim$primer_fa), nchar(prim$primer_fa)), "/",
  substr(prim$primer_fd, nchar(prim$primer_fd), nchar(prim$primer_fd)))
prim$alleles_printed <- prim$snp
prim$alleles_match <- vapply(seq_len(nrow(prim)), function(i) {
  printed <- sort(strsplit(prim$alleles_printed[i], "/")[[1L]])
  onp <- sort(strsplit(prim$alleles_on_primers[i], "/")[[1L]])
  identical(printed, onp)
}, logical(1))
prim$doi <- DOI
prim$evidence_tier <- "measured"
prim$source_file <- basename(xlsx)
prim$source_table <- "Additional file 4, columns I-K"
prim <- prim[c("marker", "snp", "primer_fa", "primer_fd", "primer_r",
               "as_3p_differs", "alleles_on_primers", "alleles_printed",
               "alleles_match", "doi", "evidence_tier", "source_file",
               "source_table")]
kasp_write_checked(prim, fs::path(OUT, "sharplab2021_primers.csv"))

message(
  "sharplab2021: ", nrow(markers), " markers (paper 50); design 50/50; ",
  "polymorphic ", n_poly, "/50 = ", n_poly / 50, " (paper 92.00%); ",
  "primer 3' bases differ on ", sum(prim$as_3p_differs), "/", nrow(prim),
  "; printed alleles match the 3' bases on ", sum(prim$alleles_match), "/",
  nrow(prim)
)
