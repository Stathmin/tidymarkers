# =============================================================================
#  data-raw/snap2caps-supplement.R  --  SNP2CAPS predicted-vs-observed CAPS
# =============================================================================
#  Source (identifier resolution is PROVENANCE, not evidence; every number
#  below is re-read from the fetched full text):
#
#    Thiel T, Kota R, Grosse I, Stein N, Graner A.  SNP2CAPS: a SNP and INDEL
#    analysis and CAPS-marker development toolbox.  Nucleic Acids Research 2004,
#    32(Web Server issue):W205-W209.  DOI 10.1093/nar/gnh006  (PMC373308)
#
#  Fetched: Europe PMC `fullTextXML` for PMC373308 returns HTTP 500; the OUP
#  full text is 403.  The working artifact is the PMC HTML page, stored as
#  .tmp/lit/snap2caps/fulltext.html, plus its text as fulltext.txt (the failed
#  EPMC response is kept verbatim as fulltext.epmc-500.json).  The script
#  parses fulltext.txt, which is where Table 2 lives.
#
#  WHY this paper: it is the only source in the registry that reports, for
#  every tested marker, the COMPUTED cut site and the OBSERVED fragment sizes
#  side by side, with the accessions that were expected to be cut and those
#  expected not to be.  That makes the prediction falsifiable: a wrong site
#  call shows up as enumerated fragments that do not sum to the amplicon.
#
#  WHAT THIS SCRIPT DOES NOT MODEL  (please read before citing a rate):
#
#   * The paper's own claim is "all 14 EST-restriction enzyme pairs ... revealed
#     the predicted restriction pattern" (Results; Table 2 footnote: "In all
#     cases the predictions of SNP2CAPS could be confirmed").  So the outcome
#     here is PREDICTION CONFIRMATION for a marker, i.e. 14/14 - NOT "14 SNPs
#     converted into genotyping assays".  Table 2 shows many accessions
#     predicted to stay uncut, because the SNP sits in only part of the panel;
#     those uncut lanes are part of the confirmation and are NOT failures.
#   * It is not a conversion rate for arbitrary SNPs: the 14 were "randomly
#     selected" only AFTER screening alignments with a set of common enzymes,
#     and SNP2CAPS itself warns that it does not judge whether a predicted
#     pattern is resolvable on a gel.
#   * The accessions that ARE cut is an accession-level polymorphism tally,
#     recorded but unasserted: the paper verified by gel, and no independent
#     second count of those lanes exists in the text to check against.
#   * This validation predates and is separate from the 2021 VCF2CAPS paper
#     (different program, different crop) - the two must not be pooled.
#   * No rate is converted into a probability for internal use.
#
#  Run:  Rscript data-raw/snap2caps-supplement.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

txt_path <- fs::path(root, ".tmp", "lit", "snap2caps", "fulltext.txt")
if (!file.exists(txt_path)) stop("Missing ", txt_path)
lines <- readLines(txt_path, encoding = "UTF-8", warn = FALSE)
out_dir <- fs::path(root, "inst", "extdata")
DOI <- "10.1093/nar/gnh006"

# Table 2 sits between its caption and the "Open in a new tab" that follows it.
start <- grep("^[[:space:]]*Table 2[.] Experimental validation of 14", lines)
if (length(start) != 1L) stop("Table 2 caption not found exactly once")
rest_idx <- (start + 1L):length(lines)
rest_idx <- rest_idx[grepl("^[[:space:]]*Open in a new tab", lines[rest_idx])]
if (!length(rest_idx)) stop("Table 2 has no closing 'Open in a new tab' marker")
end <- rest_idx[1L]

# The full text came from the PMC HTML, so every TABLE CELL is its own line
# (pdftotext -layout is NOT involved here; fulltext.txt was extracted from the
# HTML).  After dropping blank/NBSP lines the body is a sequence of cells:
# a marker id, then the enzyme, then 4-cell rows.  Blank rows in the HTML are
# already gone, so rows are consecutive quads - except that a row may carry two
# predicted sites, which the HTML puts in ONE cell ("88, 194").
tokens <- trimws(lines[(start + 1L):(end - 1L)])
tokens <- tokens[nzchar(tokens) & tokens != "\u00a0"]

header <- which(tokens == "Accessions")
if (length(header) == 1L) tokens <- tokens[seq_along(tokens) > header]

is_marker <- grepl("^GBS[0-9]+$", tokens)
if (sum(is_marker) != 14L) {
  stop("Expected 14 GBS marker tokens in Table 2, found ", sum(is_marker))
}

blocks <- list()
i <- 1L
while (i <= length(tokens)) {
  if (!is_marker[i]) { i <- i + 1L; next }
  j <- i + 1L
  while (j <= length(tokens) && !is_marker[j]) j <- j + 1L
  blocks[[length(blocks) + 1L]] <- c(tokens[i], tokens[(i + 1L):(j - 1L)])
  i <- j
}

# Body = enzyme, then rows of <product> <site(s)> <confirmed fragments> <accessions>
parse_block <- function(b) {
  marker <- b[1L]
  enzyme <- b[2L]
  v <- b[-(1:2)]
  if (length(v) %% 4L != 0L) {
    stop("Block ", marker, " has ", length(v), " cells, not a whole number of ",
         "<product> <site> <fragments> <accessions> rows: ", paste(v, collapse = " | "))
  }
  rows <- lapply(split(v, rep(seq_len(length(v) / 4L), each = 4L)), function(r) {
    sites <- if (r[2L] %in% c("\u2013", "-")) NA_integer_ else
      as.integer(strsplit(r[2L], ", ")[[1L]][1L])
    cut <- !is.na(sites)
    frags <- as.integer(strsplit(r[3L], ", ")[[1L]])
    if ((cut && is.na(sites)) || anyNA(frags)) {
      stop("Unparsable Table 2 row in ", marker, ": ", paste(r, collapse = " | "))
    }
    tibble::tibble(
      marker = marker, enzyme = enzyme, pcr_product_bp = as.integer(r[1L]),
      predicted_sites_bp = if (cut) r[2L] else NA_character_,
      predicted_cut = cut,
      observed_fragments_bp = r[3L],
      fragment_sum_bp = sum(frags),
      accessions = r[4L]
    )
  })
  dplyr::bind_rows(rows)
}
groups <- dplyr::bind_rows(lapply(blocks, parse_block))
groups$marker <- factor(groups$marker, levels = unique(groups$marker))

markers <- groups |>
  dplyr::group_by(marker, enzyme) |>
  dplyr::summarise(
    n_accession_groups = dplyr::n(),
    pcr_product_bp = paste(sort(unique(pcr_product_bp)), collapse = ","),
    n_cut_groups = sum(predicted_cut),
    n_uncut_groups = sum(!predicted_cut),
    .groups = "drop"
  )

## --- the assertion that makes this a *verified* extract --------------------
# Two checks against the paper's own statements:
#  (a) 14 marker blocks, of which GBS0560 appears twice with different enzymes
#      (HaeIII and PstI) -> 14 distinct marker x enzyme tests, as the caption
#      says ("14 randomly chosen EST-restriction enzyme pairs");
#  (b) the enumerated confirmed fragments of every group SUM to that group's
#      PCR product size (e.g. 204 + 131 = 335).  This is a real cross-check of
#      the transcription: it catches a dropped or mis-read fragment number.
#      It is also exactly what the paper claims to have observed, so a bad
#      parse surfaces as a sum failure rather than a silent wrong number.
n_tests <- nrow(unique(groups[, c("marker", "enzyme")]))
stopifnot(length(n_tests) == 1L, !is.na(n_tests))
sum_bad <- groups[!is.na(groups$pcr_product_bp) &
                    groups$fragment_sum_bp != groups$pcr_product_bp, ]
if (n_tests != 14L || length(blocks) != 14L || nrow(sum_bad) > 0L) {
  stop(
    "Table 2 transcription failed: ", length(blocks), " blocks, ", n_tests,
    " marker x enzyme tests (expected 14); ",
    nrow(sum_bad), " group(s) whose fragments do not sum to the PCR product",
    if (nrow(sum_bad)) paste0(" (e.g. ", sum_bad$marker[1], ": ",
                              sum_bad$fragment_sum_bp[1], " vs ",
                              sum_bad$pcr_product_bp[1], " bp)") else ""
  )
}

OUTCOME_DEF <- paste(
  "SNP2CAPS 'success' = the gel showed the restriction pattern the program",
  "predicted for that accession panel: fragments of the predicted sizes where",
  "the SNP leaves the site intact, and the uncut PCR product where it does not",
  "(Table 2 / Fig 3 footnote: 'In all cases the predictions of SNP2CAPS could",
  "be confirmed'). It is a prediction-confirmation outcome on one gel per",
  "marker x enzyme pair, NOT a conversion rate over arbitrary SNPs."
)

aug <- function(x) dplyr::mutate(x, outcome = "predicted pattern confirmed on gel",
                                 outcome_definition = OUTCOME_DEF, doi = DOI)
readr::write_csv(aug(markers), fs::path(out_dir, "snap2caps_markers.csv"))
readr::write_csv(aug(dplyr::ungroup(groups)),
                 fs::path(out_dir, "snap2caps_accession_groups.csv"))

message("snap2caps: ", length(blocks), " blocks / ", n_tests,
        " marker x enzyme pairs, ", nrow(groups),
        " accession-group rows, all fragment sums match (assertion passed)")
