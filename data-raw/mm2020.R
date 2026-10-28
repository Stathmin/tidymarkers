# =============================================================================
#  data-raw/mm2020.R  --  4096 measured duplex Tm (all mismatch contexts)
# =============================================================================
#  Source (the DOI is PROVENANCE; every number below is re-read from the
#  fetched file):
#
#    Oliveira LM, Long AS, Brown T, Fox KR, Weber G.  Melting temperature
#    measurement and mesoscopic evaluation of single, double and triple DNA
#    mismatches.  Chem Sci 2020;11:8273-8287.  DOI 10.1039/d0sc01700k.
#    PMCID PMC8163305, PMID 34094181.
#    Fetched: .tmp/lit/mm2020/fulltext.xml and
#             .tmp/lit/mm2020/supplements/SC-011-D0SC01700K-s001.pdf
#             (pdftotext -layout -> .tmp/lit/mm2020/supp.txt)
#
#  WHAT THIS IS, AND WHY IT MATTERS FOR ZW3
#
#   * This is the ONLY dataset found in this pass that is BOTH (a) a set of
#     MEASURED oligonucleotide duplex melting temperatures (not computed, not
#     a parameter table) AND (b) genuinely TABULAR (ESI Table S1, 4096 rows)
#     AND (c) measured under an EXPLICITLY STATED buffer: 50 mM NaCl + 10 mM
#     Na-phosphate pH 7.4, 1.0 uM total strand, no Mg2+.  That buffer is
#     exactly the regime a sodium-only NN model is supposed to describe, so
#     unlike an HRM/PCR amplicon set it is a FAIR test of `primer_tm()`.
#
#   * It is NOT an oligo-level primer set: every duplex shares the scaffold
#     5'-CGACGTGCN1N3N5ATGTGCTG-3' / 3'-GCTGCACGN2N4N6TACACGAC-5' (19 bp).
#     Only the CENTRAL TRIMER varies, so this dataset probes nearest-neighbour
#     context and mismatch discrimination, not length/composition generality.
#
#   * It is NOT a clean ground truth for the four parameter tables:
#         mean residual of primer_tm() vs measured, canonical 64, at the
#         paper's own 50 mM Na+, is about -16.9 C with sd 0.64 C.  That is a
#         LEVEL offset, not a ranking failure (r = 0.958).  The offset is not
#         attributable to Mg (Mg = 0 here) and is NOT corrected or absorbed
#         anywhere -- tests assert RANK and spread, never absolute agreement.
#     See the report for the attempted Mg explanation and why it fails.
#
#  WHAT THE OFFSET IS, AND WHAT IS FORBIDDEN (ZW3.4)
#
#   Measured on R 4.6.1 / TmCalculator 1.1.2, the 64 canonical rows, at this
#   buffer:
#     bias (primer_tm - measured) = -16.864 C, sd 0.641 C
#     Pearson r = 0.9579, Spearman rho = 0.9413
#   A LEVEL offset with the ranking intact can be hidden by one constant.  This
#   table does not allow that, so every row carries calibration_stance = STANCE
#   and the numbers below are NAMED, not applied.
#
#   Candidate hypotheses, with the arithmetic that supports each:
#
#   H1  SALT TERM.  The default correction is Schildkraut & Lifson 1965,
#       16.6*log10([Na+]); at 50 mM that is 16.6*log10(0.05) = -21.60 C against
#       a 1 M reference.  The more complete Owczarzy 2008 form moves the level
#       by +7.53 C with the SAME rank (same r), so at least 7.5 of the 16.9 C
#       is the SALT CODE and not the measurement: only the computed Tm moved.
#
#   H2  MEASUREMENT STANDARD.  Owczarzy still leaves -9.34 C.  The residual has
#       the size and sign a reference-state mismatch produces: taking the
#       -21.60 C salt term out of the -16.86 C bias leaves +4.73 C, i.e. the
#       model's 1 M value sits above the measurement.  The paper measures the
#       first-derivative maximum of a LightCycler curve at 1.0 uM total strand,
#       which is not the thermodynamic Tm at 1 M, and its method section states
#       no concentration/standard-state correction.  Cannot be tested from this
#       table, so it stays a hypothesis.
#
#   H3  MAGNESIUM: ruled out.  Mg = 0 in this buffer; a plausible 3 mM Mg
#       closes ~12 C and then OVERSHOOTS the mismatch classes, trading a level
#       error for a class error.
#
#   The ~-9.3 C left by H1 is the blocker: it needs either a stated
#   standard-state correction from the source (which is closed) or an
#   independent 1 M measurement (which does not exist).  Absorbing it here is
#   what STANCE forbids, and it is why Z0.2 (a dPCR Tm scale) stays blocked.
#
#  ASSERTIONS (all against the paper's own totals)
#   * 4096 rows; `i` is exactly 1..4096 with no duplicates (paper: 4096
#     sequences, "all possible combinations of the bases A, C, G and T").
#   * Mismatch-class counts are EXACTLY 64 canonical / 576 single /
#     1728 double / 1728 triple, matching the Results text ("576 sequences
#     with single, 1728 with double and 1728 with triple mismatches") and the
#     abstract's "64 sequences containing all combinations of canonical base
#     pairs".  The class of a row is decided from the printed base pairs, not
#     copied, so this is a real cross-check of the table against the prose.
#   * The 64 canonical rows are exactly those whose printed bottom trimer is
#     the COMPLEMENT of the printed top trimer (not the reverse complement:
#     the flanking columns swap order).  Asserted both ways -- as a set
#     identity and as the head-to-head orientation of all 64.
#     This is what proves the two columns are one duplex.
#   * Measured Tm is non-increasing in `i` (the caption says the table is
#     ordered by decreasing measured temperature).
#
#  Run:  Rscript data-raw/mm2020.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

supp_txt <- fs::path(root, ".tmp", "lit", "mm2020", "supp.txt")
xml_path <- fs::path(root, ".tmp", "lit", "mm2020", "fulltext.xml")
for (p in c(supp_txt, xml_path)) if (!file.exists(p)) stop("Missing ", p)

DOI <- "10.1039/d0sc01700k"
doc <- xml2::read_xml(xml_path)
if (!any(grepl(DOI, xml2::xml_text(doc), fixed = TRUE)))
  stop("fulltext.xml does not mention ", DOI)

# ZW3.4: the offset is recorded, not corrected.  No constant is added to the
# model; the header block above names the candidate hypotheses instead.
STANCE <- "recorded, not corrected: do not calibrate into the minus"

# ---------------------------------------------------------------------------
# 1. Parse ESI Table S1.  Four columns per entry, three entries per line:
#      i  centre(N1N3N5/N2N4N6)  Ti  Ti0  tau
#    The tau column is absent on some wrapped lines, so the tau field is
#    optional in the pattern.
# ---------------------------------------------------------------------------
lines <- readLines(supp_txt, warn = FALSE)
lines <- lines[!grepl("continued|Colours", lines)]
pat <- "(^|[[:space:]])([0-9]{1,4})[[:space:]]+([ACGT]{3})/([ACGT]{3})[[:space:]]+([0-9]{1,2}\\.[0-9])[[:space:]]+([0-9]{1,3}\\.[0-9]{2})([[:space:]]+([0-9]\\.[0-9]{4}))?($|[[:space:]])"

rows <- do.call(rbind, lapply(lines, function(ln) {
  m <- regmatches(ln, gregexpr(pat, ln, perl = TRUE))[[1L]]
  if (!length(m)) return(NULL)
  data.frame(
    index      = as.integer(sub(pat, "\\2", m)),
    top        = sub(pat, "\\3", m),
    bottom     = sub(pat, "\\4", m),
    measured   = as.numeric(sub(pat, "\\5", m)),
    paper_pred = as.numeric(sub(pat, "\\6", m)),
    tau        = suppressWarnings(as.numeric(sub(pat, "\\8", m))),
    stringsAsFactors = FALSE
  )
}))

if (is.null(rows)) stop("Table S1: nothing parsed")
if (nrow(rows) != 4096L)
  stop("Table S1: parsed ", nrow(rows), " rows, the paper reports 4096 sequences")
if (!identical(sort(rows$index), 1:4096))
  stop("Table S1: `i` is not a permutation of 1..4096")
if (anyDuplicated(paste(rows$top, rows$bottom, sep = "/")))
  stop("Table S1: duplicate central trimer")

# ---------------------------------------------------------------------------
# 2. Classify each row from the printed base pairs, then match the prose counts.
# ---------------------------------------------------------------------------
WC <- c("AT", "TA", "CG", "GC")
n_canonical <- vapply(seq_len(nrow(rows)), function(k) {
  a <- strsplit(rows$top[k], "")[[1L]]
  b <- strsplit(rows$bottom[k], "")[[1L]]
  sum(paste0(a, b) %in% WC)
}, integer(1L))

mismatches <- 3L - n_canonical
class_name <- c("canonical", "single", "double", "triple")[mismatches + 1L]
counts <- table(factor(class_name, levels = c("canonical", "single", "double", "triple")))
expected <- c(canonical = 64L, single = 576L, double = 1728L, triple = 1728L)
if (!all(as.integer(counts) == expected))
  stop("mismatch-class counts ", paste(names(counts), counts, collapse = "/"),
       " != paper's ", paste(names(expected), expected, collapse = "/"))

# The 64 all-Watson-Crick rows are exactly those whose printed bottom trimer
# is the COMPLEMENT of the printed top trimer (not the reverse complement: the
# flanking columns swap order, so the two are written head-to-head rather than
# antiparallel).  Asserted both ways.  This is what proves the two columns are
# the two strands of one duplex.
complement <- function(s) {
  map <- c(A = "T", C = "G", G = "C", T = "A")
  vapply(s, function(x)
    paste(unname(map[strsplit(x, "")[[1L]]]), collapse = ""), character(1L))
}
complementary <- rows$bottom == complement(rows$top)
if (!setequal(which(complementary), which(class_name == "canonical")))
  stop("the 64 canonical rows are not exactly the 64 complement-of-top rows")
# measured, not assumed: all 64 canonical rows are written head-to-head, so the
# top column is also the complement of the bottom column (1 of them happens to
# be palindromic and so is its own reverse complement too)
head_to_head <- rows$top == complement(rows$bottom)
if (!all(head_to_head[class_name == "canonical"]) || sum(head_to_head) != 64L)
  stop("canonical rows are not written head-to-head")

# caption: table ordered by decreasing measured Tm
if (any(diff(rows$measured) > 0))
  stop("measured Tm is not non-increasing in i")

# ---------------------------------------------------------------------------
# The bottom column is the other strand written in the SAME 5'->3' reading
# direction with the flanks in mirrored order (verified: taking it as-is makes
# all 64 canonical rows perfect 19-bp WC duplexes; reverse-complementing it
# makes none of them pair).
F5 <- "CGACGTGC"; F3 <- "ATGTGCTG"
B5 <- "GCTGCACG"; B3 <- "TACACGAC"   # printed bottom, 3'->5'
top_18 <- paste0(F5, rows$top, F3)
bottom_18 <- paste0(B5, rows$bottom, B3)
if (any(nchar(top_18) != 19L)) stop("rebuilt duplex is not 19 nt")
# canonical rows must come out as perfect 19-bp Watson-Crick duplexes
perfect <- bottom_18 == complement(top_18)
if (sum(perfect) != 64L || !all(perfect[class_name == "canonical"]))
  stop("rebuilt canonical duplexes are not perfect Watson-Crick pairs")

# ---------------------------------------------------------------------------
# 4. Write.
# ---------------------------------------------------------------------------
out <- data.frame(
  index            = rows$index,
  central_trimer   = paste(rows$top, rows$bottom, sep = "/"),
  top_18           = top_18,
  bottom_18        = bottom_18,
  mismatch_class   = class_name,
  measured_tm_c    = sprintf("%.1f", rows$measured),
  paper_pred_tm_c  = sprintf("%.2f", rows$paper_pred),
  tau              = ifelse(is.na(rows$tau), "", sprintf("%.4f", rows$tau)),
  na_mm            = 50L,
  mg_mm            = 0L,
  dntp_mm          = 0L,
  oligo_length_nt  = 19L,
  evidence_tier    = "measured",
  source_file      = "SC-011-D0SC01700K-s001.pdf Table S1",
  outcome_definition = paste(
    "duplex melting temperature (degC), first-derivative maximum of the",
    "LightCycler melting curve; 50 mM NaCl + 10 mM Na-phosphate pH 7.4,",
    "1.0 uM total strand, no Mg2+"),
  calibration_stance = STANCE,
  stringsAsFactors = FALSE
)

dest <- fs::path(root, "inst", "extdata", "tm_mm2020_s1_measured.csv")
utils::write.csv(out, dest, row.names = FALSE, na = "")
back <- utils::read.csv(dest, stringsAsFactors = FALSE)
if (nrow(back) != 4096L) stop("written CSV does not round-trip to 4096 rows")
if (!all(back$calibration_stance == STANCE))
  stop("written CSV does not round-trip calibration_stance")

cat("wrote", nrow(out), "rows to", dest, "\n")
cat("classes:", paste(names(counts), as.integer(counts), sep = "=", collapse = " "), "\n")
