# =============================================================================
#  data-raw/belv2024_assay_concordance.R  --  eleven assays, 55 pairwise
#                                             kappa / agreement cells
# =============================================================================
#  "Inter-laboratory comparison of eleven quantitative or digital PCR assays
#  for detection of proviral bovine leukemia virus", BMC Vet Res 2024,
#  10.1186/s12917-024-04228-z, PMC11346035.
#
#  WHY THIS FILE EXISTS
#  `inst/extdata/mismatch_family_*.csv` and the dPCR Tm tables each hold ONE
#  platform's readout.  This paper is the only source in the corpus that
#  prints a FULL pairwise concordance matrix across ten qPCR assays and one
#  ddPCR assay on the same samples (Table 2).  It is the real-world answer to
#  "how much does the assay choice move the call", which the design DAGs
#  currently only assert in prose.
#
#  TABLE SHAPE (Table 2, "Assessment of qPCR results by pairwise comparison")
#  Header rows: one explanatory row, then the 11 assay names.
#  Each data row: row label + 11 cells.  The matrix is triangular:
#    * cell (i, j), j > i  -> Cohen's kappa of assay i vs assay j
#    * cell (i, j), j < i  -> percent agreement of assay i vs assay j
#    * diagonal            -> "-"
#  The upper and lower triangles are DIFFERENT quantities for the same pair,
#  so both are read from the same cell coordinates and joined by (i, j).
#
#  TWO TRAPS RECORDED IN THE PAPER AS-SHIPPED
#   * One kappa prints with a DECIMAL COMMA: the paper's `0,075` is the
#     qPCR3 x qPCR5 pair.  Naive `as.numeric()` turns it into NA silently.
#   * The header's sixth assay is `ddPCR6`, not a qPCR; a positional read
#     that assumes "qPCR6" shifts every later column.  The header guard below
#     stops on that.
#
#  FETCH (verified route, 2026-11)
#    curl -sS -L -m 120 -o .tmp/lit/belv2024/full.xml \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11346035/fullTextXML'
#
#  Run:  Rscript data-raw/belv2024_assay_concordance.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

xml_path <- Sys.getenv(
  "BELV_XML",
  file.path(root, ".tmp", "lit", "belv2024", "full.xml")
)
if (!file.exists(xml_path)) {
  stop("missing ", xml_path, "; see the FETCH block at the top of this script")
}

assays <- c(
  "qPCR1", "qPCR2", "qPCR3", "qPCR4", "qPCR5", "ddPCR6",
  "qPCR7", "qPCR8", "qPCR9", "qPCR10", "qPCR11"
)
n_assay <- length(assays)

doc <- xml2::read_xml(xml_path)
tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
rows_of <- function(t) {
  lapply(
    xml2::xml_find_all(t, ".//*[local-name()='tr']"),
    function(r) {
      trimws(xml2::xml_text(
        xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']")
      ))
    }
  )
}
hit <- Filter(
  function(t) {
    grepl(
      "Assessment of qPCR results by pairwise comparison",
      xml2::xml_text(xml2::xml_find_first(t, ".//*[local-name()='caption']")),
      fixed = TRUE
    )
  },
  tw
)
if (length(hit) != 1L) {
  stop("Table 2 not found exactly once; got ", length(hit))
}
rows <- rows_of(hit[[1L]])

# the header row carries a leading empty cell (row-span corner), so drop blanks
header <- Filter(function(v) identical(v[v != ""], assays), rows)
if (length(header) != 1L) {
  stop("assay header row not found; saw: ",
       paste(vapply(rows, function(v) paste(v, collapse = "|"), character(1)), collapse = " / "))
}
data_rows <- Filter(function(v) length(v) == n_assay + 1L && v[[1L]] != "", rows)
if (length(data_rows) != n_assay) {
  stop("expected ", n_assay, " data rows, got ", length(data_rows))
}
labels <- vapply(data_rows, `[[`, character(1), 1L)
if (!identical(labels, assays)) {
  stop("row labels are ", paste(labels, collapse = ","), ", not the header order")
}

# strip the row label; the remaining 11 cells follow the header order
cells <- lapply(data_rows, function(v) v[-1L])
num_of <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", sub("%$", "", trimws(x)), fixed = TRUE)))
}
m <- matrix(NA_real_, n_assay, n_assay)
for (i in seq_len(n_assay)) {
  for (j in seq_len(n_assay)) {
    m[i, j] <- num_of(cells[[i]][[j]])
  }
}

pairs <- do.call(rbind, lapply(seq_len(n_assay - 1L), function(i) {
  data.frame(i = i, j = seq.int(i + 1L, n_assay))
}))
kappa <- m[cbind(pairs$i, pairs$j)]
agreement <- m[cbind(pairs$j, pairs$i)]

stopifnot(
  !anyNA(kappa),
  !anyNA(agreement),
  # the decimal-comma kappa, and no other comma cell in the matrix
  which(abs(kappa - 0.075) < 1e-9) == which(pairs$i == 3L & pairs$j == 5L),
  sum(grepl(",", unlist(cells), fixed = TRUE)) == 1L,
  # kappa and agreement are the two different triangles, never the same cell
  !isTRUE(all.equal(kappa, agreement))
)

out <- data.frame(
  assay_a = assays[pairs$i],
  assay_b = assays[pairs$j],
  # keep kappa as a fixed 3-decimal character so "1.000" and "0.070" survive
  kappa = sprintf("%.3f", kappa),
  agreement_pct = agreement,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(out) == 55L,
  !anyDuplicated(paste(out$assay_a, out$assay_b)),
  # measured splits; a shifted column breaks these before it reaches the CSV
  sum(kappa < 0.4) == 30L,
  sum(kappa >= 0.8) == 3L,
  identical(sort(kappa[kappa >= 0.8]), c(0.845, 0.927, 1.000)),
  min(kappa) == 0.070,
  max(kappa) == 1.000,
  min(agreement) == 66.8,
  max(agreement) == 100
)

readr::write_csv(out, "inst/extdata/belv2024_assay_concordance.csv")
cat(
  "belv2024_assay_concordance:", nrow(out), "pairs on", n_assay, "assays;",
  sprintf("kappa %.3f-%.3f (%d below 0.4, %d at/above 0.8);",
          min(kappa), max(kappa), sum(kappa < 0.4), sum(kappa >= 0.8)),
  sprintf("agreement %.1f-%.1f %%\n", min(agreement), max(agreement))
)
