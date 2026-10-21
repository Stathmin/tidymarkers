# =============================================================================
#  data-raw/sorghum_arms2021.R  --  a FIFTH independent mismatch-position probe,
#                                   and the first one whose position is read off
#                                   the primer sequence instead of prose
# =============================================================================
#  Zhu X, Wu M, Deng R, et al. "Amplification refractory mutation system (ARMS)-PCR
#  for waxy sorghum authentication with single-nucleotide resolution."
#  Foods 2021;10:2218.  doi:10.3390/foods10092218  PMID 34574328  PMC8467681  (OA)
#
#  WHY THIS PAPER MATTERS TO Z7
#  The package carries a deliberate-mismatch POSITION contradiction (Ye 2001 says
#  -2, Liu 2012 measures -2, Lefever 2019 measures -3). Zhu 2021 is a crop paper
#  that independently scanned the position, and unlike the others its scan is
#  written down as primer SEQUENCES (Table S1), so the position is not a matter
#  of interpretation: we align each allele-specific primer to the paper's own
#  perfectly-matching primer and read the differing columns off the alignment.
#
#  WHAT THE PAPER SAYS, AND WHAT THE SEQUENCES SAY
#  Text: "when the mismatched base was placed fourth from the 3' end of wxt F4
#  (substitution of C for T), aARMS-PCR yielded the highest discrimination
#  ability".  Decoded by alignment, wxt F4 carries TWO mismatches, at -1 and -4,
#  not one at -4. So the paper's prose and its own sequences disagree, and the
#  CSV lands the SERQUENCE positions with the paper's prose kept verbatim
#  beside them. This is deliberate: it is a measured fact about the paper.
#
#  WHAT IS NOT LANDED
#  The F0-F3 screening per-primer Cq values exist only in Figure S2 (an image in
#  a PDF). They are the actual evidence that "F0-F3 could not specifically
#  identify", but they are not text, so no number is extracted from them. The
#  CSV marks those rows outcome_available = FALSE.
#
#  THE ASSERTIONS (each one stops)
#  (1) Table S1 parses to exactly 10 allele-specific primers, all 20 nt.
#  (2) Aligning each primer to wxc_F0 reproduces the EXACT mismatch map we
#      intend to land (-1, -2, -3, -4 positions per primer). If the PDF layout
#      changes and a sequence is truncated, the map moves and we stop.
#  (3) The CT quartet in the main text is internally consistent: for each
#      template the matched primer's CT is lower than the mismatched primer's,
#      and the discrimination (difference) is recomputed, not copied.
#  (4) The matched/mismatched CT values are exactly the four printed numbers.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!exists("root")) {
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else {
    getwd()
  }
  setwd(root)
}

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "sorghum_arms2021")
pdf_s1 <- file.path(lit, "supplements", "s001", "foods-1353494-supplementary.pdf")
ft <- file.path(lit, "fulltext.xml")
for (f in c(pdf_s1, ft)) {
  if (!file.exists(f)) {
    stop("missing ", f, "\n  fetch with:\n",
         "  curl -o ", ft, " https://www.ebi.ac.uk/europepmc/webservices/rest/",
         "PMC8467681/fullTextXML\n",
         "  curl -o ", file.path(lit, "sups.zip"),
         " https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8467681/",
         "supplementaryFiles")
  }
}

DOI <- "10.3390/foods10092218"
out_dir <- file.path(root, "inst", "extdata")

## --- 1. Table S1 from the supplementary PDF ---------------------------------
txt_file <- tempfile(fileext = ".txt")
status <- system2("pdftotext", c("-layout", shQuote(pdf_s1), shQuote(txt_file)),
                  stdout = TRUE, stderr = TRUE)
if (!file.exists(txt_file)) {
  stop("pdftotext failed on ", pdf_s1, ": ", paste(status, collapse = " "))
}
lines <- readLines(txt_file, warn = FALSE)
unlink(txt_file)

# One primer per line: a name field then a 20-mer. Anchored on the 20-mer so a
# caption sentence cannot be mistaken for a primer row. R1 (the common reverse)
# is ALSO all-ACGT and would match, so it is excluded explicitly: this extract is
# about the allele-specific forward primers, where the deliberate mismatch sits.
pat <- "^[[:space:]]*([Ww]x[ct][[:space:]]*F[0-4])[[:space:]]+([ACGT]{20})[[:space:]]*$"
hit <- regmatches(lines, regexec(pat, lines))
hit <- Filter(length, hit)
if (!length(hit)) stop("no primer rows in Table S1; PDF layout changed")
prim <- data.frame(
  primer_raw = vapply(hit, `[[`, character(1), 2L),
  sequence   = vapply(hit, `[[`, character(1), 3L),
  stringsAsFactors = FALSE
)
prim$primer <- gsub("[[:space:]]+", "_", tolower(prim$primer_raw))

if (nrow(prim) != 10L) {
  stop("Table S1 yielded ", nrow(prim), " primers, expected 10 allele-specific ",
       "forward primers (wxc/wxt F0-F4)")
}
if (!all(nchar(prim$sequence) == 20L)) {
  stop("not all Table S1 primers are 20 nt: ",
       paste(prim$primer[nchar(prim$sequence) != 20L], collapse = ", "))
}

## --- 2. read the mismatch positions off the alignment -----------------------
# wxc_F0 (no deliberate mismatch) is the reference: it differs from every other
# primer only where the authors put a base on purpose.
ref_seq <- prim$sequence[prim$primer == "wxc_f0"]
if (length(ref_seq) != 1L) stop("wxc_F0 is not unique in Table S1")

positions_from_3p <- function(seq, ref) {
  # primer 1-based index i sits at -(n - i + 1) from the 3' end (last base = -1).
  # as.integer() is load-bearing: -() on an empty index returns double(0), which
  # would make identical() fail on the primer that has no deliberate mismatch.
  idx <- which(strsplit(seq, "")[[1L]] != strsplit(ref, "")[[1L]])
  # Ordered 5'->3' along the primer, i.e. ascending (most distal from the 3' end
  # first), so a two-mismatch primer reads (-4, -1).
  as.integer(sort(-(nchar(seq) - idx + 1L), decreasing = FALSE))
}
prim$mismatch_positions <- stats::setNames(
  lapply(prim$sequence, positions_from_3p, ref = ref_seq), prim$primer)
prim$base_at_mismatch <- vapply(seq_len(nrow(prim)), function(i) {
  seq <- prim$sequence[i]
  idx <- which(strsplit(seq, "")[[1L]] != strsplit(ref_seq, "")[[1L]])
  if (!length(idx)) return("")
  paste0(strsplit(seq, "")[[1L]][idx], collapse = "/")
}, character(1))

# Expected map, read by hand from the Table S1 sequences and asserted here.
# Two-mismatch primers list positions 5'->3' along the primer (see the sorter).
expected <- list(
  wxc_f0 = integer(0), wxc_f1 = -2L, wxc_f2 = -2L, wxc_f3 = -3L, wxc_f4 = -2L,
  wxt_f0 = -1L, wxt_f1 = c(-2L, -1L), wxt_f2 = c(-3L, -1L),
  wxt_f3 = c(-2L, -1L), wxt_f4 = c(-4L, -1L)
)
for (nm in names(expected)) {
  got <- as.integer(prim$mismatch_positions[[nm]])
  if (!identical(got, expected[[nm]])) {
    stop("mismatch map for ", nm, " is (", paste(got, collapse = ", "),
         "), expected (", paste(expected[[nm]], collapse = ", "), "); ",
         "Table S1 layout or a sequence changed")
  }
}

# The textual claims we are testing against. If either sentence is absent the
# CSV would be quoting something the paper does not say.
ft_text <- paste(xml2::xml_text(xml2::xml_find_all(xml2::read_xml(ft), "//text()[not(ancestor::ref-list)]")),
                collapse = " ")
ft_text <- gsub("[[:space:]]+", " ", ft_text)
quote_fourth <- "when the mismatched base was placed fourth from the 3"
quote_second <- "introduced second from the 3"
if (!grepl(quote_fourth, ft_text, fixed = TRUE)) {
  stop("the 'fourth from the 3-prime end' sentence is not in the full text")
}
if (!grepl(quote_second, ft_text, fixed = TRUE)) {
  stop("the 'second from the 3-prime end' sentence is not in the full text")
}

## --- 3. the CT quartet, recomputed ------------------------------------------
# Figure 1 caption: (C) amplification curves, (D) CT values, both "with
# allele-specific wxt F4 and wxc F4"; (E) gel at 50 ng. Text states the quartet.
ct_text <- regmatches(ft_text, regexpr(
  "C T values of .Hongyingzi. amplified by wxc F4 and wxt F4 were [0-9.]+ and [0-9.]+, respectively",
  ft_text))
ct_numbers <- as.numeric(regmatches(ct_text, gregexpr("[0-9]+\\.[0-9]+", ct_text))[[1L]])
if (length(ct_numbers) != 2L) stop("could not read the Hongyingzi CT pair")
ct2_text <- regmatches(ft_text, regexpr(
  "C T of .Jinnuoliang No.6. amplified by wxt F4 and wxc F4 were [0-9.]+ and [0-9.]+",
  ft_text))
if (!length(ct2_text) || grepl("Hongyingzi", ct2_text)) {
  stop("could not isolate the Jinnuoliang CT pair")
}
ct2_numbers <- as.numeric(regmatches(ct2_text, gregexpr("[0-9]+\\.[0-9]+", ct2_text))[[1L]])
if (length(ct2_numbers) != 2L) stop("could not read the Jinnuoliang CT pair")

cts <- data.frame(
  template        = c("Hongyingzi", "Jinnuoliang No.6"),
  primer          = c("wxc_f4", "wxt_f4"),
  role            = c("matched", "matched"),
  ct              = c(ct_numbers[1L], ct2_numbers[1L]),
  mismatch_primer = c("wxt_F4", "wxc_F4"),
  ct_mismatch     = c(ct_numbers[2L], ct2_numbers[2L]),
  stringsAsFactors = FALSE
)
cts$delta_ct <- cts$ct_mismatch - cts$ct
# Direction assertion: a wrong pairing makes the matched primer the slower one.
if (any(cts$delta_ct <= 0)) {
  stop("matched primer is not faster than the mismatched one: ",
       paste(cts$template, cts$delta_ct, sep = "=", collapse = "; "))
}

## --- 4. write the outcome table ---------------------------------------------
outcome_definition <- paste(
  "CT = qPCR cycle threshold (Cq) at 400 ng template DNA.",
  "matched = primer fully complementary to that template;",
  "mismatch = the other allele-specific primer.",
  "delta_ct = CT(mismatched primer) - CT(matched primer): larger = more",
  "discrimination. Design screening of F0-F3 is reported only in Figures S1/S2",
  "and is NOT extracted (outcome_available = FALSE)."
)

out <- data.frame(
  primer_group   = prim$primer,
  sequence       = prim$sequence,
  n_mismatch     = vapply(prim$mismatch_positions, length, integer(1)),
  mismatch_positions = vapply(prim$mismatch_positions, function(x) {
    if (!length(x)) return("0") else paste(x, collapse = ";")
  }, character(1)),
  base_at_mismatch = prim$base_at_mismatch,
  paper_claim    = vapply(prim$primer, function(p) {
    switch(p,
      wxt_f4 = "mismatched base placed fourth from the 3' end of wxt F4 gave the highest discrimination",
      wxc_f4 = "mismatch introduced second from the 3' end of wxc F4 (A for C) made it specific",
      paste0("F0-F3 screening: primers ", p, " did not specifically identify the target"))
  }, character(1)),
  outcome_available = FALSE,
  outcome_definition = outcome_definition,
  evidence_tier  = "primary literature",
  source_file    = "sorghum_arms2021/supplements/s001/foods-1353494-supplementary.pdf (Table S1); fulltext.xml Fig. 1 text",
  doi            = DOI,
  stringsAsFactors = FALSE
)
# Attach the two measured rows' CT outcomes to the two scored primers.
out$ct_matched <- NA_real_
out$ct_mismatched <- NA_real_
out$delta_ct <- NA_real_
for (i in seq_len(nrow(cts))) {
  j <- which(out$primer_group == cts$primer[i])
  if (length(j) != 1L) stop("no unique row for scored primer ", cts$primer[i])
  out$ct_matched[j] <- cts$ct[i]
  out$ct_mismatched[j] <- cts$ct_mismatch[i]
  out$delta_ct[j] <- cts$delta_ct[i]
  out$outcome_available[j] <- TRUE
}

# The paper says only F4 (one per allele) was scored: exactly two rows carry a
# CT outcome, and they are the two the text names.
if (sum(out$outcome_available) != 2L ||
    !setequal(out$primer_group[out$outcome_available], c("wxc_f4", "wxt_f4"))) {
  stop("scored primers are not exactly the wxc_F4/wxt_F4 pair the text names")
}

outfile <- file.path(out_dir, "arms_sorghum2021_outcomes.csv")
utils::write.csv(out, outfile, row.names = FALSE, na = "")
message("wrote ", outfile, " (", nrow(out), " primers; ",
        sum(out$outcome_available), " with a measured CT outcome)")
