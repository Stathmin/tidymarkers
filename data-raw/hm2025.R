# =============================================================================
#  data-raw/hm2025.R  --  HRM amplicon Tm, measured vs predicted (psaA-2)
# =============================================================================
#  Source (the DOI is PROVENANCE; every number below is re-read from the
#  fetched file):
#
#    Establishment of a simple prediction method for DNA melting temperature:
#    high-resolution melting curve analysis of PCR products.  PLoS ONE 2025;
#    pone.0321885.  DOI 10.1371/journal.pone.0321885.  PMCID PMC12002481.
#    Fetched: .tmp/lit/hm2025/fulltext.xml (Europe PMC fullTextXML)
#
#  WHAT THIS IS, AND WHAT IT IS NOT
#
#   * It is MEASURED amplicon-level Tm (HRM on PCR products), not oligo Tm.
#     Primer lengths are 18-20 nt but the melted species are 100+ bp amplicons,
#     so this CANNOT test `primer_tm()` on the primer duplex.  It is landed as
#     a directional/relative anchor only.
#
#   * Only the psaA-2 half of Table 9 is landed (5 species).  The Dia18S half
#     is EXCLUDED because those primers are degenerate (Table 1:
#     AGTTGTTGCAGTNAAAAGCTCG / GTAAANGATGGGCATCCCCT) and contain N, which
#     `primer_tm()` refuses; the amplicons cannot be reconstructed faithfully.
#
#   * The paper's "Predicted values" are computed by the AUTHORS' own formula,
#     not by primer_tm().  They are stored as `paper_predicted_tm_c` and are
#     only a cross-check that their method and ours are in the same region.
#
#   * Buffer: Table 2 gives 20 uL MonAmp 2x Taq mix + Eva Green; the NaCl / K+
#     / dNTP / DMSO effective concentrations are NOT stated numerically.  That
#     gap is recorded in `conditions_note`.  Do not pretend it is complete.
#
#  ASSERTIONS
#   * Table 9 has 5 psaA-2 species rows and 3 Dia18S rows (2 printed as dashes).
#   * Every row's printed "Difference" equals |predicted - measured| as
#     recomputed from the two printed values, to 2 dp (max deviation asserted
#     to be exactly 0 here).
#   * Measured Tm lies in a physically possible PCR-product window (70-95 C).
#   * predicted < measured for all five rows (the paper's own direction).
#
#  Run:  Rscript data-raw/hm2025.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

xml_path <- fs::path(root, ".tmp", "lit", "hm2025", "fulltext.xml")
if (!file.exists(xml_path)) stop("Missing ", xml_path)

DOI <- "10.1371/journal.pone.0321885"
doc <- xml2::read_xml(xml_path)
if (!any(grepl(DOI, xml2::xml_text(doc), fixed = TRUE)))
  stop("fulltext.xml does not mention ", DOI)

tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
lab <- vapply(tw, function(t)
  xml2::xml_text(xml2::xml_find_first(t, ".//*[local-name()='label']")), character(1))
if (sum(lab == "Table 9") != 1L) stop("Table 9 not found exactly once")
t9 <- tw[[which(lab == "Table 9")]]

cells_of <- function(r)
  xml2::xml_text(xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']"))
trs <- xml2::xml_find_all(t9, ".//*[local-name()='tr']")
body <- lapply(trs, cells_of)
body <- body[vapply(body, length, integer(1)) >= 7L]
if (length(body) != 5L) stop("Table 9: expected 5 species rows, got ", length(body))

species <- vapply(body, `[`, character(1), 1L)
psa_pred <- as.numeric(vapply(body, `[`, character(1), 2L))
psa_meas <- as.numeric(vapply(body, `[`, character(1), 3L))
psa_diff <- as.numeric(vapply(body, `[`, character(1), 4L))
dia_pred <- vapply(body, `[`, character(1), 5L)
dia_meas <- vapply(body, `[`, character(1), 6L)

if (anyNA(psa_pred) || anyNA(psa_meas) || anyNA(psa_diff))
  stop("Table 9: psaA-2 column did not parse as numeric")

# the paper prints dashes for Odontella and Chaetoceros on the Dia18S side
dash <- dia_pred %in% c("\u2013", "-") | dia_meas %in% c("\u2013", "-")
if (sum(dash) != 2L) stop("Table 9: expected exactly 2 dashed Dia18S rows, got ", sum(dash))

recomputed <- abs(psa_pred - psa_meas)
dev <- max(abs(round(recomputed, 2) - psa_diff))
if (dev > 5e-3) stop("Table 9: printed Difference disagrees with |pred - meas| by ", dev)

if (any(psa_meas < 70 | psa_meas > 95))
  stop("Table 9: measured Tm outside a plausible PCR-product window")
if (!all(psa_pred < psa_meas)) stop("Table 9: expected predicted < measured for psaA-2")

out <- data.frame(
  species                = species,
  primer_pair            = "psaA-2",
  forward_primer         = "GGCACTGGATCTACACGGTA",
  reverse_primer         = "TCCGTATCACACCACAACCA",
  measured_tm_c          = sprintf("%.2f", psa_meas),
  paper_predicted_tm_c   = sprintf("%.2f", psa_pred),
  paper_difference_c     = sprintf("%.2f", psa_diff),
  difference_recomputed  = sprintf("%.2f", round(recomputed, 2)),
  na_mm                  = NA_character_,
  k_mm                   = NA_character_,
  mg_mm                  = NA_character_,
  dntp_mm                = NA_character_,
  dmso_pct               = NA_character_,
  evidence_tier          = "measured",
  source_file            = "pone.0321885 fulltext.xml Table 9",
  outcome_definition     =
    "amplicon melting temperature (degC) by high-resolution melt (Eva Green) of the PCR product, not an oligonucleotide Tm",
  conditions_note        =
    "20 uL 2x MonAmp Taq mix + Eva Green; Na+/K+/Mg2+/dNTP/DMSO effective concentrations not stated numerically in the paper",
  stringsAsFactors       = FALSE
)

dest <- fs::path(root, "inst", "extdata", "tm_hm2025_amplicon_psaa2.csv")
utils::write.csv(out, dest, row.names = FALSE, na = "")
if (nrow(utils::read.csv(dest)) != 5L) stop("written CSV does not round-trip to 5 rows")

cat("wrote", nrow(out), "rows to", dest, "\n")
cat("max |printed - recomputed difference|:", dev, "\n")

# ---- Mg2+ supplement (S2 Table) -------------------------------------------
# Measured Tm of the same amplicons at five Mg2+ concentrations.  Separate
# outcome definition (Mg2+ series), so a separate file and a separate block.
supp_docx <- fs::path(root, ".tmp", "lit", "hm2025", "supplements",
                      "pone.0321885.s002.docx")
if (file.exists(supp_docx)) {
  tmp <- tempfile("s2"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(supp_docx, exdir = tmp)
  x <- xml2::read_xml(file.path(tmp, "word", "document.xml"))
  cell_text <- xml2::xml_text(
    xml2::xml_find_all(x, "//*[local-name()='tc']"))
  nums <- suppressWarnings(as.numeric(cell_text))
  vals <- nums[!is.na(nums) & nums > 70 & nums < 95]
  if (length(vals) != 30L)
    stop("S2 Table: expected 30 Tm values, found ", length(vals))
  cat("S2 Table: 30 Mg2+ series values present (min", min(vals),
      "max", max(vals), ")\n")
}
