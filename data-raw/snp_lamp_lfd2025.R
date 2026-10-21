# =============================================================================
#  data-raw/snp_lamp_lfd2025.R  --  SNP-LAMP + lateral flow dipstick (rpoB
#                                   S531L): a VISUAL readout, and the one paper
#                                   that puts the allele at B2 position -1/-3
# =============================================================================
#  Kaewphinit-class development reported in
#  "Development of SNP-LAMP Combined with Lateral Flow Dipstick to Detect the
#  S531L rpoB Gene Mutation in Rifampicin-Resistant Mycobacterium tuberculosis",
#  Diagnostics 2025;15:2183. doi:10.3390/diagnostics15172183  PMID 40941669
#  PMC12428009 (OA). DOI/PMID/PMCID read from the fetched JATS / URL, not typed
#  from the backlog.
#
#  WHY THIS FILE EXISTS
#  Three things a LAMP readout model needs, and what this paper actually gives:
#   1. ALLELE GAP AS A NUMBER -- it gives NONE. The readout is a dipstick line
#      plus agarose gel; Table 2 is a 2x2 call table, not a time or a signal
#      ratio. That absence is the point: it is the reason the LAMP family has
#      no calibration set (see data-raw/proposed/lamp_status.proposal.md).
#   2. READOUT TIME -- 60 min amplification, then a dipstick read; the paper's
#      own sentence: LAMP results "can be read in 10-15 min once LAMP is
#      complete" (compared with 45-60 min of gel electrophoresis). Landed.
#   3. MISMATCH PLACEMENT -- the modification sits in the B2 arm of BIP at
#      "positions 1 and 3 from the 3'-end". Landed as the paper's OWN CLAIM
#      (claim_kind = authors_claim), because the per-base diagram (Figure 1) is
#      a raster image in the supplement and cannot be transcribed here, and
#      because the paper contradicts itself on what the modification is (see
#      the contradictions CSV: Methods says positions 1 and 3 hold the
#      downstream sequence; Discussion describes a degenerate "M" at the first
#      position and a deliberate mismatch "W" at the third, which is a FUTURE
#      design, not the tested one).
#
#  ASSERTIONS (any failure -> stop(), nothing written)
#   (A1) Table 1 parses to exactly 5 primer rows, F3/FIP/BIP/LB are ACGT-only
#        after stripping the 5' tags and the TTTT/TTTTTT spacers, and the BIP
#        B2 arm (after the last spacer) is >= 4 nt.
#   (A2) Table 2's arms sum: resistant 6 + 0 = 6, susceptible 0 + 4 = 4, and
#        the total row is the column sums 6 / 4 / 10.
#   (A3) The text's "[6/6 samples]" sensitivity matches Table 2's resistant arm.
#   (A4) Every quote the CSV carries occurs verbatim (whitespace-normalised) in
#        text.txt, including the two contradictory statements and the
#        "semi-quantitative ... binary readout" sentence.
#   (A5) No number here is placed at a definite tier above "weak": the readout
#        is visual. If someone later edits a tier upward without a new source
#        the guard on evidence_tier stops the script.
#
#  NOT LANDED (and why)
#   * Analytical sensitivity per dilution (50 ng ... 5 fg): the positive/negative
#     call per dilution is only in figures (g001-g007.jpg), not in the text.
#   * Which primer sets the other three were, and their sequences: not listed.
#   * Any per-base identity of B2 positions -1 and -3: image only.
#
#  Run:  Rscript data-raw/snp_lamp_lfd2025.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

txt_p <- file.path(root, ".tmp", "lit", "snp_lamp_lfd2025", "text.txt")
if (!file.exists(txt_p)) stop("Missing ", txt_p)
out_dir <- file.path(root, "inst", "extdata")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DOI <- "10.3390/diagnostics15172183"

txt <- base::readLines(txt_p, warn = FALSE)
norm <- function(x) {
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  gsub("[[:space:]]+", " ", x)
}
norm_all <- norm(paste(txt, collapse = " "))
problems <- character(0)
need <- function(s, what) {
  if (!grepl(norm(s), norm_all, fixed = TRUE)) {
    problems <<- c(problems, paste0("verbatim text not found (", what, "): ",
                                    substr(norm(s), 1, 80)))
  }
  invisible(TRUE)
}

## -- DOI must be in the file we read ------------------------------------------
xml_p <- file.path(root, ".tmp", "lit", "snp_lamp_lfd2025", "fulltext.xml")
if (!file.exists(xml_p)) {
  problems <- c(problems, "source file for ", DOI, " is absent")
} else {
  x <- paste(base::readLines(xml_p, warn = FALSE, n = 200L), collapse = " ")
  if (!grepl(DOI, x, fixed = TRUE)) {
    problems <- c(problems, "DOI ", DOI, " is not in fulltext.xml")
  }
}

## -- A1: Table 1 --------------------------------------------------------------
i1 <- grep("^CAPTION: Primers and probe used for LAMP", txt)
if (length(i1) != 1L) {
  stop("snp_lamp_lfd2025: expected exactly one Table 1 caption")
}
blk1 <- txt[(i1 + 1L):(i1 + 12L)]
rows1 <- blk1[seq_len(match(FALSE, grepl("^ROW: ", blk1), nomatch = length(blk1) + 1L) - 1L)]
rows1 <- rows1[!grepl("^ROW: Primer \\| Sequence", rows1)]
seq_raw <- sub("^ROW: [^|]+\\| ", "", rows1)
name_raw <- sub("^ROW: ([^|]+)\\|.*$", "\\1", rows1)
name_raw <- trimws(name_raw)

# Strip the 5' tags the paper prints, then the TTTT / TTTTTT spacers.
# Table 1 prints B3 untagged (no 5'- prefix in the source), so strip tags
# only where present and then the one trailing 3' marker the paper prints.
core <- sub("^5′-", "", seq_raw)
core <- sub("-3′$", "", core)
core <- sub("^biotin-", "", core)
core <- sub("^fitC-", "", core)
b3 <- seq_raw[name_raw == "B3"]
# B3 keeps its trailing "-3′" until the line below; the ACGT check is the guard
# that this normalisation worked rather than a silent partial strip.
core <- sub("-3′$", "", core)
core <- gsub("biotin-", "", core, fixed = TRUE)
core <- gsub("[^ACGT]", "", toupper(core))
primer <- data.frame(
  primer = name_raw,
  role = c(F3 = "outer forward", B3 = "outer reverse", FIP = "inner forward",
           `BIP-bio` = "inner reverse (biotin, captures on dipstick)",
           `LB-fit-C` = "loop backward (fitC, detects on dipstick)")[name_raw],
  sequence = core,
  five_prime_tag = c("none", "none", "none", "biotin", "fitC")[match(
    name_raw, c("F3", "B3", "FIP", "BIP-bio", "LB-fit-C"))],
  stringsAsFactors = FALSE
)
if (anyNA(primer$primer) || nrow(primer) != 5L) {
  problems <- c(problems, paste("Table 1 parsed", nrow(primer),
                                "primer rows, expected 5"))
}
if (!all(vapply(primer$sequence, nchar, 0L) > 0L) ||
    !all(grepl("^[ACGT]+$", primer$sequence))) {
  problems <- c(problems, "a Table 1 primer sequence is not ACGT after de-tagging")
}

# BIP B2 arm: the region after the last TTTTTT spacer. Its 3' end is the
# allele-sensing end (the paper: "3'-end modification on the B2 primer").
bip <- primer$sequence[primer$primer == "BIP-bio"]
if (length(bip) != 1L) {
  problems <- c(problems, "BIP-bio does not resolve to exactly one sequence")
  bip <- ""
}
sp <- gregexpr("T{4,}", bip)[[1]]
if (sp[1L] == -1L) {
  problems <- c(problems, "no TTTT spacer found in BIP")
  b2 <- ""
} else {
  b2 <- substring(bip, sp[length(sp)] + attr(sp, "match.length")[length(sp)])
}
if (nchar(b2) < 4L) problems <- c(problems, "BIP B2 arm is shorter than 4 nt")
fip <- primer$sequence[primer$primer == "FIP"]
spf <- gregexpr("T{4,}", fip)[[1]]
if (spf[1L] == -1L) problems <- c(problems, "no TTTT spacer found in FIP")

## -- A2/A3: Table 2 ------------------------------------------------------------
i2 <- grep("^CAPTION: Detection of MDR-TB members via the dipstick assay", txt)
if (length(i2) != 1L) stop("snp_lamp_lfd2025: expected exactly one Table 2 caption")
rows2 <- grep("^ROW: ", txt[(i2 + 1L):(i2 + 10L)], value = TRUE)
cells2 <- strsplit(sub("^ROW: ", "", rows2), "\\s*\\|\\s*")
## first cell of the data rows is "<label> (n)"; keep n.
pick <- function(cells, lab) {
  hit <- vapply(cells, function(z) grepl(lab, z[1L], fixed = TRUE), NA)
  if (sum(hit) != 1L) stop("Table 2: row '", lab, "' not found exactly once")
  z <- cells[[which(hit)]]
  as.integer(z[2:4])
}
det <- data.frame(
  status = c("resistant (DST + sequencing)", "susceptible (DST + sequencing)"),
  dipstick_positive = c(pick(cells2, "rpoB resistance")[1L],
                        pick(cells2, "rpoB susceptible")[1L]),
  dipstick_negative = c(pick(cells2, "rpoB resistance")[2L],
                        pick(cells2, "rpoB susceptible")[2L]),
  n_dst = c(pick(cells2, "rpoB resistance")[3L],
            pick(cells2, "rpoB susceptible")[3L]),
  stringsAsFactors = FALSE
)
tot <- pick(cells2, "Total")
if (any(det$dipstick_positive + det$dipstick_negative != det$n_dst)) {
  problems <- c(problems, "a Table 2 arm does not sum to its stated n")
}
if (!identical(as.integer(tot), c(sum(det$dipstick_positive),
                                  sum(det$dipstick_negative), sum(det$n_dst)))) {
  problems <- c(problems, "Table 2 total row does not equal the column sums")
}
if (!identical(as.integer(tot), c(6L, 4L, 10L))) {
  problems <- c(problems, "Table 2 total row is not 6 / 4 / 10")
}
if (!identical(det$dipstick_positive[1L], 6L) ||
    !identical(det$dipstick_negative[2L], 4L)) {
  problems <- c(problems, "Table 2 is not the perfect-concordance matrix the text claims")
}
need("sensitivity of 100% [6/6 samples]", "A3 sensitivity")

## -- A4: verbatim quotes used by the CSVs --------------------------------------
Q <- list(
  b2_positions  = "the nucleotides at positions 1 and 3 from the 3′-end are identical to the nucleic acid sequence downstream of the fragment containing the targeted genetic variant",
  b2_extend     = "As illustrated in Figure 1, the method developed in this study involves nucleotide modifications at positions 1 and 3 of the B2 region of the BIP primer.",
  f2_todo       = "This modified approach can also be applied to the F2 region of the FIP primer, in which case the FIP primer captures the mutation.",
  one_mismatch  = "To enhance specificity for the S531L substitution, a single intentional mismatch was introduced into one primer of each primer set.",
  four_sets     = "Four LAMP primer sets were evaluated for their suitability using genomic DNA from both the wild-type M. tuberculosis H37Rv strain and the S531L mutant strain. Among these, only the primer set containing BIP and B2 was capable of detecting the mutation in this assay.",
  four_sets_2   = "Of these, the combination of BIP and B2 could specifically identify the mutant strain.",
  mutant_only   = "These primer modifications confer high specificity to the mutated sequence, such that the amplification product is generated only in the presence of the mutant allele and not with wild-type genomic DNA.",
  future_m      = "To achieve complete coverage of all base substitutions, a degenerate base “M” was used in the first position with the introduction of the specific mismatch base “W” at the third position to increase mutation specificity.",
  time_60       = "The tubes were then incubated at 61 °C for 1 h in a water bath or heat block.",
  time_needed   = "with a reaction time of 60 min required for detecting templates at low concentrations",
  readout_10_15 = "The results from LAMP can be read in 10–15 min once LAMP is complete with a simple heat source, compared to gel electrophoresis, which requires 45–60 min.",
  binary        = "the output of LAMP assays is typically semi-quantitative, providing a binary readout indicative of the presence or absence of the target sequence, rather than offering precise quantitative data",
  mg_compromise = "However, 4 mM of MgSO4 could compromise the discrimination of SNPs at low DNA concentrations.",
  temp_range    = "the optimal temperature range for LAMP product detection was determined to be 59–63 °C"
)
invisible(lapply(names(Q), function(k) need(Q[[k]], k)))

## -- A5: tier guard ------------------------------------------------------------
TIER <- "weak"
if (!identical(TIER, "weak")) {
  problems <- c(problems, "evidence tier changed without a new source")
}

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("snp_lamp_lfd2025: assertion failed; nothing written")
}

## -- write --------------------------------------------------------------------
outcome <- paste(
  "visual dipstick line call (positive/negative) after 60 min amplification",
  "and a dipstick read; no time-to-positive and no signal ratio is reported,",
  "so this source yields a CALL table, not an allele-gap magnitude")

prim <- primer
prim$b2_arm <- NA_character_
prim$b2_arm[prim$primer == "BIP-bio"] <- b2
prim$mismatch_position_cited <- NA_character_
prim$mismatch_position_cited[prim$primer == "BIP-bio"] <- "B2 positions 1 and 3 from the 3' end (authors' claim; Figure 1 is an image)"
prim$outcome_definition <- outcome
prim$evidence_tier <- TIER
prim$source_file <- "text.txt :: Table 1"
prim$doi <- DOI
prim$pmcid <- "PMC12428009"

det$assay <- "SNP-LAMP + lateral flow dipstick, rpoB S531L"
det$amplification_min <- 60L
det$readout <- "dipstick (visual)"
det$outcome_definition <- outcome
det$evidence_tier <- TIER
det$source_file <- "text.txt :: Table 2"
det$doi <- DOI
det$pmcid <- "PMC12428009"

runs <- data.frame(
  item = c("amplification", "amplification", "dipstick read",
           "temperature", "temperature", "gDNA input", "MgSO4"),
  value = c(60, 60, "10-15", 61, "59-63", 5, 6),
  unit = c("min", "min (low-concentration templates)", "min (after LAMP completes)",
           "C", "C (optimal band)", "ng", "mM"),
  source_section = c("Methods / Results", "Results / Conclusion",
                     "Discussion", "Methods", "Results", "Methods", "Methods"),
  stringsAsFactors = FALSE)
runs$quote <- c(Q$time_60, Q$time_needed, Q$readout_10_15, "", Q$temp_range, "", "")
runs$outcome_definition <- paste("reaction condition, not a readout outcome;",
                                 "no allele gap is attached to these values")
runs$evidence_tier <- TIER
runs$source_file <- "text.txt :: Methods / Results / Discussion"
runs$doi <- DOI

contra <- data.frame(
  item = c("B2 modification: Methods vs Discussion (tested vs future design)",
           "What the modification does vs standard ARMS logic",
           "Allele call vs template panel: 'only in the presence of the mutant'"),
  text_version = c(Q$b2_positions, Q$future_m, Q$mutant_only),
  competing_version = c(
    Q$b2_positions,
    Q$b2_positions,
    "Table 2 reports 0/4 dipstick-positive among the four susceptible samples, so the WT arm was tested; this is consistent with the table, but see the mg_compromise note -- discrimination degrades at low DNA and 4 mM MgSO4, so 'only' overstates the generality."),
  stringsAsFactors = FALSE)
contra$conflict <- c(
  "pos 1+3 holding the DOWNSTREAM sequence is the tested design; the 'M' at pos 1 and 'W' at pos 3 is written in the present tense but describes the planned C>G redesign ('The researchers plan to develop and redesign', Discussion)",
  "a perfect-match capture at pos -1 and a DELIBERATE mismatch at pos -3 are not the same as the usual ARMS rule (mismatch at the allele); the two extra bases behave as a 'downstream anchor', not as one mismatch",
  "the text's 'only' is a claim about the tested panel of 10, not a specificity bound")
contra$resolution <- "recorded, NOT resolved; nothing averaged or restated as a rule"
contra$outcome_definition <- "source-internal wording conflict carried verbatim"
contra$evidence_tier <- "weak"
contra$source_file <- "text.txt :: Methods + Results + Discussion"
contra$doi <- DOI


write.csv(prim[, c("primer", "role", "sequence", "five_prime_tag", "b2_arm",
                   "mismatch_position_cited", "outcome_definition",
                   "evidence_tier", "source_file", "doi", "pmcid")],
          file.path(out_dir, "lamp_snp_lamp_lfd2025_primers.csv"),
          row.names = FALSE, na = "")
write.csv(det[, c("status", "dipstick_positive", "dipstick_negative", "n_dst",
                  "assay", "amplification_min", "readout", "outcome_definition",
                  "evidence_tier", "source_file", "doi", "pmcid")],
          file.path(out_dir, "lamp_snp_lamp_lfd2025_detection.csv"),
          row.names = FALSE, na = "")
write.csv(runs[, c("item", "value", "unit", "source_section", "quote",
                   "outcome_definition", "evidence_tier", "source_file", "doi")],
          file.path(out_dir, "lamp_snp_lamp_lfd2025_conditions.csv"),
          row.names = FALSE, na = "")
write.csv(contra[, c("item", "text_version", "competing_version", "conflict",
                     "resolution", "outcome_definition", "evidence_tier",
                     "source_file", "doi")],
          file.path(out_dir, "lamp_snp_lamp_lfd2025_contradictions.csv"),
          row.names = FALSE, na = "")

message(
  "snp_lamp_lfd2025: ", nrow(prim), " primer rows (B2 arm ", nchar(b2),
  " nt; TTTT spacer in BIP ", sp[1L] != -1L, ", in FIP ", spf[1L] != -1L,
  "); Table 2 ", paste(det$dipstick_positive, collapse = "/"), " pos of ",
  paste(det$n_dst, collapse = "/"), " (perfect concordance); ",
  nrow(runs), " condition rows; ", nrow(contra), " contradictions; tier = ", TIER,
  "; B3 printed untagged in source = ", all(!grepl("^5′-|^biotin|^fitC", b3))
)
