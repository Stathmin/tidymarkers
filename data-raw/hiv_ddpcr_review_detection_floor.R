# =============================================================================
#  data-raw/hiv_ddpcr_review_detection_floor.R  --  per-patient DRM calls at
#                                                  the bulk/ddPCR detection floor
# =============================================================================
#  "Low-Frequency Drug Resistance Mutations in 16 Patients Sensitive to Sanger
#  Sequencing", PMC12188199, Table 2.
#
#  WHY THIS FILE EXISTS
#  Every design DAG in this package reasons about a minority allele that the
#  reference workflow cannot see.  This table is the corpus's measured version
#  of that claim: 16 patients Sanger-negative, each called by bulk PCR and by
#  ddPCR, with the variant allele frequency printed per mutation.  It supplies
#  the numbers behind "bulk misses below ~20 %, ddPCR reports down to 0.01".
#
#  TABLE SHAPE (Table 2)
#    row 0  6-cell span header ("Next-Generation Sequencing")
#    row 1  Case ID | ART at Time of DRM Test | RNA, Copies/mL | cDNA, ng/uL
#           | Bulk PCR | ddPCR | Current | Outcomes
#    rows 2..17  16 patients
#
#  DERIVATION RULES
#   * a DRM cell reads `M46I (0.05), I50L (0.02)`; the parenthesised value is a
#     VAF FRACTION (0.05 = 5 %), read as printed.
#   * `...` in both assay cells -> patient 1758: nothing reported, all derived
#     columns NA.  That is NOT the same as zero discordant variants.
#   * `none` in a `*_only_drms` cell means "the set difference is empty";
#     NA means "the paper reports no calls at all".
#   * `M46IV` is verbatim from the paper (no such substitution exists; it is
#     almost certainly M46I/M46V).  Kept as printed so the extract stays
#     checkable against the source.
#   * `RNA, Copies/mL` is printed with digit-group spaces ("3 000 000").
#
#  FETCH (verified route, 2026-11)
#    curl -sS -L -m 120 -o .tmp/lit/hiv2025/full.xml \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12188199/fullTextXML'
#
#  Run:  Rscript data-raw/hiv_ddpcr_review_detection_floor.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

xml_path <- Sys.getenv(
  "HIV_XML",
  file.path(root, ".tmp", "lit", "hiv2025", "full.xml")
)
if (!file.exists(xml_path)) {
  stop("missing ", xml_path, "; see the FETCH block at the top of this script")
}

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
      "16 Patients Sensitive to Sanger Sequencing",
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

expected_header <- c(
  "Case ID", "ART at Time of DRM Test", "RNA, Copies/mL", "cDNA, ng/uL",
  "Bulk PCR", "ddPCR", "Current", "Outcomes"
)
# the micro sign and the Greek mu both appear in the wild; compare on "u"
norm <- function(x) gsub("[\u00b5\u03bc]", "u", x)
hdr <- Filter(function(v) identical(norm(v), expected_header), rows)
if (length(hdr) != 1L) {
  stop("8-column header row not found; saw:\n",
       paste(vapply(rows, function(v) paste(v, collapse = " | "), character(1)),
             collapse = "\n"))
}
data_rows <- Filter(
  function(v) length(v) == 8L && grepl("^[0-9]+$", v[[1L]]),
  rows
)
if (length(data_rows) != 16L) {
  stop("expected 16 patient rows, got ", length(data_rows))
}

# `NAME (vaf)` tokens; the VAF is a printed fraction
token_re <- "([A-Z][0-9]+[A-Z]*)\\s*\\(\\s*([0-9.]+)\\s*\\)"
parse_drms <- function(cell) {
  if (identical(cell, "\u2026") || identical(cell, "...")) {
    return(NULL) # not reported in the paper
  }
  m <- regmatches(cell, gregexpr(token_re, cell, perl = TRUE))[[1L]]
  if (!length(m)) {
    stop("no DRM token in cell: ", cell)
  }
  list(
    name = sub(token_re, "\\1", m, perl = TRUE),
    vaf = as.numeric(sub(token_re, "\\2", m, perl = TRUE))
  )
}

# digit-group spaces: "3 000 000" -> 3000000 (any space flavour)
copies <- function(x) as.numeric(gsub("[[:space:]]", "", x))

out <- do.call(rbind, lapply(data_rows, function(v) {
  bulk <- parse_drms(v[[5L]])
  digital <- parse_drms(v[[6L]])
  if (is.null(bulk) || is.null(digital)) {
    return(data.frame(
      case_id = as.integer(v[[1L]]),
      rna_copies_per_ml = copies(v[[3L]]),
      n_bulk = NA_integer_, n_ddpcr = NA_integer_,
      min_vaf = NA_real_, max_vaf = NA_real_,
      ddpcr_only_drms = NA_character_, bulk_only_drms = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  set_or_none <- function(x) if (length(x)) paste(x, collapse = ";") else "none"
  data.frame(
    case_id = as.integer(v[[1L]]),
    rna_copies_per_ml = copies(v[[3L]]),
    n_bulk = length(bulk$name),
    n_ddpcr = length(digital$name),
    min_vaf = min(bulk$vaf, digital$vaf),
    max_vaf = max(bulk$vaf, digital$vaf),
    # keep each assay's own reading order, not alphabetical
    ddpcr_only_drms = set_or_none(digital$name[!digital$name %in% bulk$name]),
    bulk_only_drms = set_or_none(bulk$name[!bulk$name %in% digital$name]),
    stringsAsFactors = FALSE
  )
}))

stopifnot(
  nrow(out) == 16L,
  !anyDuplicated(out$case_id),
  # measured totals: ddPCR calls more DRMs than bulk on this cohort
  sum(out$n_bulk, na.rm = TRUE) == 31L,
  sum(out$n_ddpcr, na.rm = TRUE) == 39L,
  # the detection floor the table demonstrates
  min(out$min_vaf, na.rm = TRUE) == 0.01,
  max(out$max_vaf, na.rm = TRUE) == 0.10,
  # of 16 patients: 6 carry a ddPCR-only DRM, 2 a bulk-only DRM, 1 is
  # unreported, so 9 and 13 columns legitimately read `none`
  sum(out$ddpcr_only_drms %in% "none" & !is.na(out$ddpcr_only_drms)) == 9L,
  sum(!out$ddpcr_only_drms %in% "none" & !is.na(out$ddpcr_only_drms)) == 6L,
  sum(!out$bulk_only_drms %in% "none" & !is.na(out$bulk_only_drms)) == 2L,
  # patient 1758 is the paper's `...`: not zero, unknown.  The RNA copy count
  # IS reported, so only the derived columns are all-NA
  all(is.na(unlist(out[out$case_id == 1758L,
    c("n_bulk", "n_ddpcr", "min_vaf", "max_vaf",
      "ddpcr_only_drms", "bulk_only_drms")]))),
  # the two bulk-only rows are real, not a parse artefact
  out$bulk_only_drms[out$case_id == 1931L] == "M46I",
  out$bulk_only_drms[out$case_id == 2443L] == "M46IV"
)

readr::write_csv(out, "inst/extdata/hiv_ddpcr_review_detection_floor.csv")
cat(
  "hiv_ddpcr_review_detection_floor:", nrow(out), "patients;",
  sprintf("bulk %d DRMs, ddPCR %d DRMs;", sum(out$n_bulk, na.rm = TRUE),
          sum(out$n_ddpcr, na.rm = TRUE)),
  sprintf("VAF %.2f-%.2f; %d ddPCR-only, %d bulk-only, 1 unreported\n",
          min(out$min_vaf, na.rm = TRUE), max(out$max_vaf, na.rm = TRUE),
          sum(!out$ddpcr_only_drms %in% "none" & !is.na(out$ddpcr_only_drms)),
          sum(!out$bulk_only_drms %in% "none" & !is.na(out$bulk_only_drms)))
)
