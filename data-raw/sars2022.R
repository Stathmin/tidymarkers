# =============================================================================
#  data-raw/sars2022.R  --  the double mismatch, measured as a Ct gap
# =============================================================================
#  Hale R, et al. "Highly sensitive and specific detection of the SARS-CoV-2
#  Delta variant by double-mismatch allele-specific real time reverse
#  transcription PCR." J Clin Virol 2022;146:105049.
#  doi:10.1016/j.jcv.2021.105049  PMC8628602 (CC-BY)
#
#  WHY THIS FILE EXISTS
#  Z7's core design claim is that a SECOND deliberate mismatch, four bases from
#  the 3' end, turns a weak allele-specific primer into a reliable one. Lefever
#  2019 measured that on synthetic templates with SYBR Green; this paper is the
#  test on real samples: two parallel reactions whose allele-specific primers
#  differ ONLY in the 3'-terminal base, each also carrying a deliberate
#  internal mismatch at the 4th base from the 3' end (marked by <underline> in
#  Table 1). The readout is a Ct gap between the matching and the non-matching
#  reaction among 42 clinical samples.
#
#  WHAT IS LANDED
#  * The four allele-specific primer sequences and the common probe, with the
#    artificial mismatch base and the 3'-terminal base parsed from the MARKUP
#    (<underline> = the artificial mismatch, <bold><italic> = the 3'-terminal
#    base), never transcribed by eye.
#  * Per-sample Ct for the two reactions (Table 3), and the derived gap
#    Ct(non-matching reaction) - Ct(matching reaction) for every sample where
#    BOTH reactions gave a Ct. "Negative" is not a Ct and is never imputed.
#  * The paper's own summary of that gap ("typically around 17 cycles
#    (range 14.3 to 18.8 cycles)"), which is ASSERTED against the values
#    recomputed here. It holds only for the three samples in which the
#    non-matching reaction gave a Ct at all.
#
#  WHAT IS DELIBERATELY NOT LANDED
#  * The 17-cycle headline is NOT a property of the double mismatch in general:
#    in 38 of the 42 samples the non-matching reaction produced NO Ct, i.e.
#    the gap is unbounded, and the authors' 14.3-18.8 range is measured on the
#    4 samples with a weak cross-reaction. The CSV says which is which; the
#    tests keep the two apart.
#  * No reaction efficiency, no mismatch Tm penalty and no cycle series.
#
#  Run:  Rscript data-raw/sars2022.R [--dir .tmp/lit/sars2022]
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
  file.path(root, ".tmp", "lit", "sars2022")
ft <- file.path(lit, "fulltext.xml")
if (!file.exists(ft)) stop("missing ", ft)
DOI <- "10.1016/j.jcv.2021.105049"
SRC <- ".tmp/lit/sars2022/fulltext.xml"

doc <- read_xml(ft)
if (!length(xml_find_all(doc, "//article"))) stop("not an article XML: ", ft)
tabs <- xml_find_all(doc, "//table-wrap")
labs <- vapply(tabs, function(t) xml_text(xml_find_first(t, "./label")), character(1))
if (!all(c("Table 1", "Table 2", "Table 3") %in% labs)) {
  stop("the three tables of the paper were not found")
}
tab <- function(label) tabs[[match(label, labs)]]
rows_of <- function(t) {
  lapply(xml_find_all(t, ".//tr"), function(r) {
    cells <- xml_find_all(r, "./*[local-name()='td' or local-name()='th']")
    trimws(gsub("[[:space:]]+", " ", vapply(cells, xml_text, character(1))))
  })
}

## --- 1. Table 1: sequences, with the markup read as the legend says ---------
# legend: underlined = the deliberate internal mismatch; bold italic = the
# 3'-terminal base. Reads one <td> and returns the sequence plus both marks.
read_seq <- function(tab, name) {
  tds <- xml_find_all(tab, ".//*[local-name()='td']")
  hit <- which(vapply(tds, xml_text, character(1)) == name)
  if (length(hit) != 1L) stop("oligo row not found: ", name)
  seq_cell <- tds[[hit + 1L]]
  seq_txt <- gsub("[[:space:]]+", "", xml_text(seq_cell))
  seq_txt <- sub("^5.\u2032?", "", seq_txt)   # the 5' mark
  seq_txt <- sub("^FAM-", "", seq_txt)
  seq_txt <- sub("-BHQ.*$", "", seq_txt, useBytes = FALSE)
  seq_txt <- gsub("\u2021", "", seq_txt)  # the quencher's footnote dagger
  u <- xml_find_all(seq_cell, ".//underline")
  b <- xml_find_all(seq_cell, ".//bold")
  if (length(u) > 1L) stop(name, ": more than one underlined base")
  term <- if (length(b)) xml_text(b) else NA_character_
  # the probe's only bold node is the footnote dagger of the quencher
  if (!is.na(term) && !(term %in% c("A", "C", "G", "T"))) term <- NA_character_
  list(seq = seq_txt,
       mismatch_base = if (length(u)) xml_text(u) else NA_character_,
       terminal_base = term,
       len = nchar(seq_txt))
}
t1 <- tab("Table 1")
seqs <- lapply(c("Non-Delta_DMAS_F", "Non-Delta_DMAS_R", "Delta_DMAS_F",
                 "Delta_DMAS_R", "Common_A&B_Prb"), function(nm) read_seq(t1, nm))
# the probe carries no deliberate mismatch (its superscript marks the quencher
# footnote), so its two marks are recorded as NA and excluded below
names(seqs) <- c("non_delta_fwd", "non_delta_rev", "delta_fwd", "delta_rev", "probe")

# the two forward primers differ only at the 3'-terminal base
if (substr(seqs$non_delta_fwd$seq, 24L, 24L) == substr(seqs$delta_fwd$seq, 24L, 24L)) {
  stop("the two forward primers do not differ at the 3' base")
}
if (substr(seqs$non_delta_fwd$seq, 1L, 23L) != substr(seqs$delta_fwd$seq, 1L, 23L)) {
  stop("the two forward primers differ somewhere other than the 3' base")
}
# the artificial mismatch is 4 bases from the 3' end on every allele-specific primer
for (nm in c("non_delta_fwd", "non_delta_rev", "delta_fwd", "delta_rev")) {
  s <- seqs[[nm]]
  if (substr(s$seq, s$len - 3L, s$len - 3L) != s$mismatch_base) {
    stop(nm, ": the underlined base is not the 4th from the 3' end")
  }
  if (substr(s$seq, s$len, s$len) != s$terminal_base) {
    stop(nm, ": the bold base is not the 3'-terminal base")
  }
  if (s$mismatch_base != "C") {
    stop(nm, ": the paper chooses C as the artificial mismatch, got ", s$mismatch_base)
  }
  if (!(s$len %in% c(24L, 22L))) stop(nm, ": unexpected primer length ", s$len)
}

## --- 2. Table 3: per-sample Ct, and the gap between the two reactions --------
rows <- rows_of(tab("Table 3"))
# drop header rows and empty cells
rows <- rows[!grepl("^Sample No", vapply(rows, `[`, character(1), 1L), ignore.case = TRUE)]
keep <- vapply(rows, function(r) length(r) >= 6L, logical(1))
rows <- rows[keep]
vals <- do.call(rbind, lapply(rows, function(r) {
  data.frame(sample = r[1L], genotype_sequencing = r[2L], e_gene_ct = r[3L],
             ct_non_delta = r[4L], ct_delta = r[5L], call = r[6L],
             stringsAsFactors = FALSE)
}))
num <- function(x) suppressWarnings(as.numeric(x))
is_neg <- function(x) x == "Negative"
if (nrow(vals) != 42L) stop("Table 3 parsed as ", nrow(vals), " samples, expected 42")
if (!all(vals$call %in% c("Delta variant", "Not Delta"))) {
  stop("unexpected genotype call: ", paste(unique(vals$call), collapse = " / "))
}
vals$ct_non_delta_num <- num(vals$ct_non_delta)
vals$ct_delta_num <- num(vals$ct_delta)
# "Negative" must be the only non-numeric value, never a number in disguise
bad <- !is_neg(vals$ct_non_delta) & is.na(vals$ct_non_delta_num)
if (any(bad)) stop("unparsed non-Delta Ct: ", vals$ct_non_delta[bad])
bad <- !is_neg(vals$ct_delta) & is.na(vals$ct_delta_num)
if (any(bad)) stop("unparsed Delta Ct: ", vals$ct_delta[bad])

# which reaction matches the sample's genotype: Delta -> reaction B, else A
vals$delta <- vals$call == "Delta variant"
vals$ct_matching <- ifelse(vals$delta, vals$ct_delta_num, vals$ct_non_delta_num)
vals$ct_nonmatching <- ifelse(vals$delta, vals$ct_non_delta_num, vals$ct_delta_num)
vals$nonmatching_detected <- !is_neg(ifelse(vals$delta, vals$ct_non_delta,
                                            vals$ct_delta))
vals$ct_gap <- ifelse(vals$nonmatching_detected,
                      round(vals$ct_nonmatching - vals$ct_matching, 1), NA_real_)
# samples whose own genotype could not be established ("Insufficient virus") are
# kept but carry no gap: the paper calls them by the reaction that amplified
# "Insufficient virus" and "Not done" are the paper's own labels: for those
# samples the genotype comes from the reaction call alone, so the gap there is
# a consistency check, not an orthogonal validation.

## --- 3. the paper's own statements about the gap ----------------------------
body <- xml_text(xml_find_all(doc, "//body//p"))
flat <- paste(body, collapse = " ")
has <- function(p) grepl(p, flat, fixed = TRUE)
if (!has("around 17 cycles (range 14.3 to 18.8 cycles)")) {
  stop("the '17 cycles (range 14.3 to 18.8)' sentence is gone")
}
if (!has("four samples (1, 2, 20 and 37)")) {
  stop("the four cross-reacting samples are no longer named")
}
# FOUR samples cross-reacted (1, 2, 20 and 37), all with E-gene Ct <= 22. The
# paper quotes no per-sample gap; read from Table 3 they are 17.4 (sample 1,
# 37.4 - 20.0), 17.3 (2), 14.3 (20) and 18.8 (37), which is exactly the
# "range 14.3 to 18.8" the text states, with min and max taken over all four.
gap_ref <- c("1" = 17.4, "2" = 17.3, "20" = 14.3, "37" = 18.8)
gap4 <- vals$ct_gap[match(names(gap_ref), vals$sample)]
if (!isTRUE(all.equal(unname(gap4), unname(gap_ref)))) {
  stop("recomputed gaps for samples 1, 2, 20, 37 are ",
       paste(gap4, collapse = "/"), ", not 17.4/17.3/14.3/18.8")
}
gaps <- sort(stats::na.omit(vals$ct_gap))
if (length(gaps) != 4L || min(gaps) != 14.3 || max(gaps) != 18.8) {
  stop("the recomputed Ct-gap range is ", min(gaps), " to ", max(gaps),
       ", the paper states 14.3 to 18.8")
}
# sample 1 is the one the paper describes by hand: B = 20.0, A = 37.4
if (!identical(as.numeric(vals$ct_gap[vals$sample == "1"]), 17.4)) {
  stop("sample 1's gap is not 37.4 - 20.0 = 17.4")
}
# and every sample with a cross-reaction is one of the four high-viral-load ones
if (!setequal(vals$sample[!is.na(vals$ct_gap)], c("1", "2", "20", "37"))) {
  stop("a Ct gap exists outside samples 1, 2, 20 and 37")
}

## --- 4. write ---------------------------------------------------------------
out <- data.frame(
  sample = vals$sample,
  genotype_by_sequencing = vals$genotype_sequencing,
  genotype_call = vals$call,
  e_gene_ct = vals$e_gene_ct,
  ct_non_delta_reaction = vals$ct_non_delta,
  ct_delta_reaction = vals$ct_delta,
  ct_matching = vals$ct_matching,
  ct_nonmatching = ifelse(vals$nonmatching_detected,
                          as.character(vals$ct_nonmatching), "Negative"),
  nonmatching_reaction_detected = vals$nonmatching_detected,
  ct_gap = vals$ct_gap,
  design = "double-mismatch allele-specific RT-PCR: 3'-terminal allele base plus a deliberate C at the 4th base from the 3' end, one common hydrolysis probe",
  outcome_definition = paste(
    "Ct values as printed in Table 3 of the paper. ct_gap = Ct(non-matching",
    "reaction) - Ct(matching reaction), computed here; NA where the non-matching",
    "reaction gave no Ct ('Negative'), in which case the gap is unbounded and is",
    "NOT set to an arbitrary number. The paper's summary 'around 17 cycles (range",
    "14.3 to 18.8 cycles)' is the min and max of the FOUR samples where the",
    "non-matching reaction gave any Ct at all (1, 2, 20, 37); the other 38",
    "samples have no gap, so for them the method's margin is unbounded, not 17."),
  evidence_tier = "measured",
  source_file = SRC,
  source_table = "Table 3",
  doi = DOI,
  stringsAsFactors = FALSE
)
write.csv(out, file.path(root, "inst", "extdata", "sars2022_dmas_ct.csv"),
          row.names = FALSE, na = "")

primers <- data.frame(
  oligo = names(seqs),
  role = c("allele-specific forward (reaction A, non-Delta)",
           "allele-specific reverse (reaction A, non-Delta)",
           "allele-specific forward (reaction B, Delta)",
           "allele-specific reverse (reaction B, Delta)",
           "common hydrolysis probe (antisense strand)"),
  sequence = vapply(seqs, `[[`, character(1), "seq"),
  length_nt = vapply(seqs, `[[`, integer(1), "len"),
  primer_3p_base = ifelse(is.na(vapply(seqs, `[[`, character(1), "terminal_base")),
                          NA_character_,
                          vapply(seqs, `[[`, character(1), "terminal_base")),
  artificial_mismatch_base = vapply(seqs, `[[`, character(1), "mismatch_base"),
  artificial_mismatch_offset = ifelse(is.na(vapply(seqs, `[[`, character(1),
                                          "mismatch_base")), NA_integer_, -3L),
  outcome_definition = paste(
    "Read from the MARKUP of Table 1: <underline> is the deliberate internal",
    "mismatch and <bold><italic> the 3'-terminal base, as the table's own",
    "legend states. artificial_mismatch_offset = -3 is OUR offset convention",
    "(-1 = the 3'-terminal base), i.e. the 4th base from the 3' end, which is",
    "the position the text describes. The paper picks C as the artificial",
    "mismatch because it is the least stable choice."),
  evidence_tier = "measured",
  source_file = SRC,
  source_table = "Table 1",
  doi = DOI,
  stringsAsFactors = FALSE
)
write.csv(primers, file.path(root, "inst", "extdata", "sars2022_dmas_primers.csv"),
          row.names = FALSE, na = "")

message("sars2022: ", nrow(out), " samples from Table 3, ",
        sum(out$nonmatching_reaction_detected), " with a detectable non-matching",
        " reaction (gaps ", paste(gaps, collapse = "/"),
        "); 5 oligos from Table 1 (artificial C at -3, 3' base is the allele)")
