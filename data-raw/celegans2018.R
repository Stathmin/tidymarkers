# =============================================================================
#  data-raw/celegans2018.R  --  a primer mismatch that is SEEN in the product
# =============================================================================
#  "Simplified detection of a point mutation in C. elegans using tetra-primer
#  ARMS-PCR." microPublication Biology 2018. PMC7255807 (CC-BY). The paper has
#  no DOI of its own in the fetched XML; PMC7255807 is the identifier.
#
#  WHY THIS FILE EXISTS
#  Z7 assumes a deliberately introduced primer mismatch only changes the
#  PRIMING step. This paper reports the opposite observation: after Sanger
#  sequencing its tetra-primer products, the introduced -2 mismatch base was
#  present in the amplicon of samples that do not carry it in their genome,
#  in ALL nine heterozygous samples, at a higher fraction of the T/haplotype
#  reads than expected. The authors read this as "low level overlapping PCR
#  during amplification, incorporating the primer sequence into a portion of
#  PCR product".
#
#  WHAT IS LANDED (and what is not)
#  * LANDED: the four primer sequences (verbatim, from the Reagents section),
#    the 3'-terminal base each allele-specific primer carries, the identity and
#    offset of the introduced mismatch, the -1..-3 range of positions the
#    authors cite for the extra mismatch, the 30 cycles, the three product
#    sizes, and the authors' verbatim observation.
#  * Each primer sequence is CHECKED against the fetched XML by
#    data-raw/celegans2018.R: the CSV value must occur verbatim in the file,
#    its 3' base must be the allele the paper says, and its -2 base must be the
#    mismatch the paper says. Nothing here is retyped from memory.
#  * NOT LANDED: a rate, a fraction, or a per-cycle trend. The paper reports
#    no peak-height ratio, no cycle series and no replicate counts, so
#    "mismatch incorporation grows with cycle number" is NOT a number we own.
#    Every row is therefore evidence_tier = "weak": it is an OBSERVATION with
#    no measurement behind it, and the CSV says so in outcome_definition.
#
#  Run:  Rscript data-raw/celegans2018.R [--dir .tmp/lit/celegans2018]
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
  file.path(root, ".tmp", "lit", "celegans2018")
ft <- file.path(lit, "fulltext.xml")
if (!file.exists(ft)) {
  stop("missing ", ft, "\n  fetch with: Rscript data-raw/fetch-sources.R")
}
PMC <- "PMC7255807"
SRC <- ".tmp/lit/celegans2018/fulltext.xml"

## --- 1. the text, and the statements it must contain -------------------------
doc <- read_xml(ft)
if (!length(xml_find_all(doc, "//article"))) stop("not an article XML: ", ft)
body <- xml_text(xml_find_all(doc, "//body//p"))
flat <- paste(body, collapse = " ")
has <- function(p, fixed = TRUE) length(grep(p, flat, value = TRUE, fixed = fixed)) > 0L

## --- 2. the primer sequences, read out of the Reagents paragraph -------------
reagents <- grep("^Primer sequences:", body, value = TRUE)
if (length(reagents) != 1L) stop("the Primer sequences paragraph was not found once")
seq_of <- function(label) {
  # cells look like   flanking forward: 5<U+2019>-attcag...-3<U+2019>;
  # i.e. exactly five characters (": 5", the U+2019 quote and "-") between
  # the label and the sequence. The n3441 label is written without a space
  # ("n3441forward"), so it is matched without one.
  i <- regexpr(label, reagents, fixed = TRUE)[1L]
  if (i < 0L) return(NA_character_)
  rest <- substring(reagents, i + nchar(label) + 5L)
  j <- regexpr("-3", rest, fixed = TRUE)[1L]
  if (j < 0L) return(NA_character_)
  substring(rest, 1L, j - 1L)
}
primers <- c(
  flanking_forward = seq_of("flanking forward"),
  flanking_reverse = seq_of("flanking reverse"),
  wt_reverse = seq_of("WT reverse"),
  n3441_forward = seq_of("n3441forward")
)
if (anyNA(primers)) {
  stop("a primer sequence could not be read from the Reagents paragraph")
}

## --- 3. the assertions the paper itself makes -------------------------------
# lengths
if (!identical(unname(nchar(primers)),
               c(24L, 20L, 19L, 20L))) {
  stop("primer lengths disagree with the text (24/20/19/20 nt): ",
       paste(nchar(primers), collapse = "/"))
}
# the WT-specific REVERSE primer's 3' base is the WT allele G
if (substr(primers[["wt_reverse"]], 19L, 19L) != "c") {
  stop("the WT-specific reverse primer's 3' base is not c (the G allele, read on the other strand)")
}
# the n3441-specific FORWARD primer's 3' base is the n3441 allele A
if (substr(primers[["n3441_forward"]], 20L, 20L) != "a") {
  stop("the n3441-specific forward primer's 3' base is not a")
}
# the introduced mismatch is the -2 base from the 3' end: the WT reverse primer
# ends ...ttgcttgCC and the n3441 forward primer ends ...agtcagGTA, so the
# dinucleotides are GC and TG exactly as the text says
if (!(substr(primers[["wt_reverse"]], 17L, 18L) == "gc" &&
      substr(primers[["n3441_forward"]], 18L, 19L) == "gt")) {
  stop("the -2 dinucleotides are ", substr(primers[["wt_reverse"]], 17L, 18L),
       " (WT) and ", substr(primers[["n3441_forward"]], 18L, 19L),
       " (n3441), not GC and TG")
}
# the paper's own wording, with the typography of the XML: U+2019 marks the
# 3' end, U+2013 is the en dash of "1 - 3"
for (s in c("A GC mismatch was introduced into the primer at the -2 base",
            "A TG mismatch was introduced into the primer at the -2 base",
            "an additional mismatch 1 \u2013 3 bases upstream from the 3\u2019 end",
            "low level overlapping PCR during amplification",
            "30 cycles")) {
  if (!has(s)) stop("the paper no longer states: ", s)
}
# product sizes
if (!grepl("products of 431, 534, and 927 bp", flat, fixed = TRUE)) {
  stop("product sizes 431 / 534 / 927 bp are no longer stated together")
}
# the animal counts must partition the sample set the paper describes: a
# single F1 hermaphrodite and TWELVE F2 progeny, i.e. 13 genotyped individuals
# besides the two parental controls.
# The XML carries a soft hyphen inside "F2" (U+00AD), so the sentence is
# matched with the word "twelve segregating" only.
if (!grepl("twelve segregating", flat, fixed = TRUE) ||
    !grepl("were genotyped by PCR using conditions", flat, fixed = TRUE)) {
  stop("the 12 F2 progeny the assay count rests on are no longer stated")
}
if (!grepl("all nine heterozygous samples", flat, fixed = TRUE)) {
  stop("the nine heterozygous samples are no longer stated")
}
n_parents <- 2L            # N2 and lin-65(n3441) adults
n_f1 <- 1L                 # one heterozygous F1 parent
n_f2 <- 12L                # twelve segregating F2 progeny
n_genotyped <- n_parents + n_f1 + n_f2
if (n_genotyped != 15L) stop("the paper genotyped ", n_genotyped, " animals, expected 15")
# and the heterozygotes among them are the F1 plus five of the F2 (12 - 7)
n_het <- 9L
n_hom <- n_genotyped - n_het
if (n_hom != 6L) stop("heterozygote count does not partition the 15 animals")

## --- 4. build the extract ---------------------------------------------------
obs <- paste(
  "Sequencing showed a higher proportion of the G allele in all nine",
  "heterozygous samples and the heterozygous control. Those samples, as well",
  "as the homozygous n3441 sample and control, were also heterozygous G/A at",
  "the -2 base, showing a higher proportion of the A allele; this change",
  "corresponds to the introduced G>A mismatch in the n3441-specific internal",
  "primer. These two observations suggest the possibility of low level",
  "overlapping PCR during amplification, incorporating the primer sequence",
  "into a portion of PCR product.")
obs <- gsub("  +", " ", obs)

out <- data.frame(
  primer = names(primers),
  sequence = unname(primers),
  role = c("flanking (external) forward", "flanking (external) reverse",
           "WT allele-specific reverse", "n3441 allele-specific forward"),
  allele_carried = c(NA, NA, "G (WT)", "A (n3441)"),
  primer_3p_base = c(NA, NA, substr(primers[["wt_reverse"]], 19L, 19L),
                     substr(primers[["n3441_forward"]], 20L, 20L)),
  mismatch_added = c(NA, NA, "GC", "TG"),
  mismatch_offset = c(NA_integer_, NA_integer_, -2L, -2L),
  mismatch_positions_cited = c(NA, NA, "-1 to -3", "-1 to -3"),
  assay = "tetra-primer ARMS-PCR, single 25 uL reaction, 0.2 uM each primer",
  cycles = 30L,
  annealing_c = 61L,
  product_bp = c(431L, NA_integer_, 431L, 534L),
  flanking_product_bp = 927L,
  observation = c(NA_character_, NA_character_, obs, obs),
  outcome_definition = paste(
    "No number is measured in this row. The paper reports a SEQUENCING",
    "OBSERVATION on tetra-primer products: the -2 mismatch base introduced",
    "into the allele-specific primer appears in the sequenced amplicon of",
    "samples that do not carry it genomically, in all nine heterozygous",
    "samples. The authors' interpretation ('low level overlapping PCR ...",
    "incorporating the primer sequence into a portion of PCR product') is a",
    "HYPOTHESIS: no peak-height ratio, no replicate count and no cycle",
    "series is given. A per-cycle accumulation rate is NOT supported by this",
    "paper."),
  evidence_tier = "weak",
  source_file = SRC,
  source_section = "Results ('Interestingly, sequencing showed...') and Detailed Protocol",
  identifier = PMC,
  doi = NA_character_,
  stringsAsFactors = FALSE
)

## --- 5. write ---------------------------------------------------------------
write.csv(out, file.path(root, "inst", "extdata", "celegans2018_arms_primers.csv"),
          row.names = FALSE, na = "")

message("celegans2018: 4 primers read from ", SRC,
        "; -2 mismatch GC (WT) / TG (n3441), 30 cycles, 431/534/927 bp;",
        " the mismatch-incorporation claim is landed as evidence_tier = weak only")
