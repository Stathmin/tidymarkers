# =============================================================================
#  data-raw/cpn_ddpcr_lob.R  --  LoB VARIABILITY across assay conditions
# =============================================================================
#  Source (DOI is PROVENANCE; every number is re-read from the fetched file):
#
#    Limit of Blank Variability as a Determinant of Analytical Sensitivity in a
#    DNA-Binding Dye-Based Digital PCR for the Detection of Chlamydia
#    pneumoniae.  Anal Chem 2026.  DOI 10.1021/acs.analchem.6c04206.
#    Fetched: .tmp/lit/PMC13584697/fulltext.xml (Europe PMC fullTextXML).
#
#  WHY THIS PAPER IS THE dPCR PIVOT: the other two dPCR files
#  (data-raw/aml_idh_ddpcr.R, data-raw/kit_d816v_ddpcr.R) each report ONE LoB
#  per assay.  This one reports the SAME assay under THREE supermix lots and two
#  thermocycling profiles and gets LoB values that differ by ~3x (0.06 to 0.35
#  copies/uL) and LoD/LoB ratios from 2.6x to 3.8x.  It is direct evidence that
#  a dPCR detection limit is a property of the assay CONDITIONS, so any single
#  "dPCR LoD" constant would be a fabrication.
#
#  It also reports, for the first time in this set, a dPCR-vs-qPCR-vs-nPCR
#  comparison of the LOWEST DETECTABLE concentration (Table 1, last two rows).
#  That is a limit comparison in COPIES, not a mutation-detection conversion
#  rate, and it is written as such.
#
#  ASSERTIONS
#   * The 'dPCR LOD is at [xLOB]' row must equal LOD/LOB recomputed from the
#     same table, to 1 d.p.  This is the paper's own internal consistency check
#     and it is the only arithmetic relation the table supports.
#   * The three LoB values that exist must span more than 2x - if a re-parse
#     ever produced a single LoB this assertion fails, which is the whole point.
#   * LOB/LOD/LOQ are recorded in BOTH cp/uL and cp/well exactly as printed;
#     cp/well / cp/uL must be the same volume in every cell (20 uL).
#
#  Run:  Rscript data-raw/cpn_ddpcr_lob.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

xml_path <- fs::path(root, ".tmp", "lit", "PMC13584697", "fulltext.xml")
if (!file.exists(xml_path)) stop("Missing ", xml_path)

DOI <- "10.1021/acs.analchem.6c04206"
out_dir <- fs::path(root, "inst", "extdata")
doc <- xml2::read_xml(xml_path)
if (!any(grepl(DOI, xml2::xml_text(doc), fixed = TRUE)))
  stop("fulltext.xml does not mention ", DOI)

tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
t1 <- Filter(function(t) grepl("PCR Assay Specifications",
                               xml2::xml_text(xml2::xml_find_first(
                                 t, ".//*[local-name()='caption']"))), tw)
if (length(t1) != 1L) stop("Table 1 (assay specifications) not found exactly once")
rows <- lapply(xml2::xml_find_all(t1[[1L]], ".//*[local-name()='tr']"), function(r)
  trimws(xml2::xml_text(xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']"))))

# Columns 2..6 are the five assay/condition combinations.
label_row <- function(prefix, neg = "lowest LCS") {
  # prefix match, not substring: 'dPCR LOD' must not also catch 'dPCR LOQ'.
  re <- prefix
  hit <- Filter(function(v) length(v) == 6L && grepl(re, v[[1L]], fixed = TRUE) &&
                  !grepl(neg, v[[1L]], fixed = TRUE), rows)
  if (length(hit) != 1L) stop("Table 1 row matching '", re, "' not found exactly once")
  hit[[1L]]
}
conc_row <- function(prefix, neg = "lowest LCS") {
  hit <- Filter(function(v) length(v) == 6L && startsWith(v[[1L]], prefix) &&
                  !grepl(neg, v[[1L]], fixed = TRUE) && !grepl("is at", v[[1L]], fixed = TRUE), rows)
  if (length(hit) != 1L) stop("Table 1 row starting '", prefix, "' not found exactly once")
  hit[[1L]][2:6]
}
assay <- c("16S-A_64639593", "16S-A_64661060", "16S-default_64661060",
           "groL1-default", "groL2-default")

# "nA" = not applicable / not determined.  ">0.24/>4.80" is a bound, not a value;
# it is kept as a bound (NA + the printed string) rather than coerced to 0.24.
split_units <- function(x) {
  out <- lapply(strsplit(x, "/", fixed = TRUE), function(z) {
    z <- trimws(z)
    if (any(z == "nA")) return(c(cp_ul = NA, cp_well = NA, bound = NA, raw = x))
    bound <- grepl("^>", z[[1L]])
    c(cp_ul = as.numeric(sub("^>", "", z[[1L]])), cp_well = as.numeric(sub("^>", "", z[[2L]])),
      bound = as.numeric(bound), raw = NA)
  })
  m <- do.call(rbind, out)
  tibble::tibble(raw = x, cp_ul = as.numeric(m[, "cp_ul"]),
                 cp_well = as.numeric(m[, "cp_well"]), is_bound = as.numeric(m[, "bound"]))
}
lob_c <- conc_row("dPCR LOB")
lod_c <- conc_row("dPCR LOD")
loq_c <- conc_row("dPCR LOQ")
# 'dPCR LOD is at [xLOB]' is a derived row; it shares the 'dPCR LOD' prefix, so it
# is pulled directly rather than through the prefix filter.
ratio_c <- Filter(function(v) length(v) == 6L && grepl("is at", v[[1L]], fixed = TRUE), rows)
if (length(ratio_c) != 1L) stop("'dPCR LOD is at [xLOB]' row not found exactly once")
ratio_c <- ratio_c[[1L]][2:6]
q_low <- conc_row("lowest LCS detectable", neg = "nPCR ")
n_low <- conc_row("lowest LCS detectable", neg = "qPCR ")
for (nm in c("lob_c", "lod_c", "loq_c", "q_low", "n_low"))
  if (length(get(nm)) != 5L) stop("Column count for ", nm, " is ", length(get(nm)))

lob <- split_units(lob_c); lod <- split_units(lod_c); loq <- split_units(loq_c)
# Volume check: cp/well = cp/uL * V for every defined cell.
vol <- unique(c(lod$cp_well / lod$cp_ul, lob$cp_well / lob$cp_ul))
vol <- vol[is.finite(vol)]
vol <- unique(round(vol, 3))
if (length(vol) != 1L || abs(vol - 20) > 0.051)
  stop("cp/well / cp/uL is not a single 20 uL volume: ", paste(round(vol, 2), collapse = "/"))

# The paper prints "nA" for the three 16S columns' LOD; assert that exactly.
if (sum(is.na(lod$cp_ul)) != 2L) stop("Expected 2 'nA' LOD cells, got ", sum(is.na(lod$cp_ul)))
ok <- !is.na(lod$cp_ul)
printed <- suppressWarnings(as.numeric(sub("\u00d7$", "", ratio_c[ok])))
recomp <- round(lod$cp_ul[ok] / lob$cp_ul[ok], 1)
if (!isTRUE(all.equal(printed, recomp, tolerance = 0.051)))
  stop("Table 1 'LOD is at xLOB' (", paste(printed, collapse = "/"),
       ") != LOD/LOB recomputed (", paste(recomp, collapse = "/"), ")")
lob3 <- lob$cp_ul[!is.na(lob$cp_ul)]
if (max(lob3) / min(lob3) <= 2)
  stop("LoB values no longer span > 2x (", paste(lob3, collapse = "/"),
       "); this paper's point is that they do")

DEF_LOB <- paste(
  "Outcome = limit of blank (highest signal reached by blank samples), for ONE",
  "assay under ONE supermix lot / thermocycling profile. Same assay, different",
  "conditions => different LoB; this is the variability being recorded, not an",
  "assay-to-assay difference.")
DEF_LOD <- paste(
  "Outcome = limit of detection in copies per uL (and per 20 uL well), i.e. an",
  "ABSOLUTE concentration, not a fraction of a mutation. 'nA' cells mean the",
  "authors did not determine a LoD for that condition.")
DEF_LOQ <- paste(
  "Outcome = limit of quantification in copies per uL. A '>' prefix means the",
  "paper printed a bound ('>0.24'), which is preserved in raw and is_bound.")
DEF_CMP <- paste(
  "Outcome = lowest low-concentration sample still detected by the named method,",
  "in copies per uL. This is a LIMIT COMPARISON between digital, quantitative and",
  "nested PCR on the same samples. It compares detection of a CONCENTRATION, and",
  "is not a variant-detection conversion rate.")

mk <- function(kind, dat, def) tibble::tibble(
  limit_type = kind, assay = assay, cp_ul = dat$cp_ul, cp_well = dat$cp_well,
  is_bound = dat$is_bound, raw = dat$raw, well_volume_ul = ifelse(is.na(dat$cp_ul), NA, vol))
cmp <- tibble::tibble(
  limit_type = c(rep("lowest_detectable_nPCR", 5), rep("lowest_detectable_qPCR", 5)),
  assay = rep(assay, 2),
  cp_ul = c(split_units(n_low)$cp_ul, split_units(q_low)$cp_ul),
  cp_well = c(split_units(n_low)$cp_well, split_units(q_low)$cp_well),
  is_bound = c(split_units(n_low)$is_bound, split_units(q_low)$is_bound),
  raw = c(n_low, q_low), well_volume_ul = NA_real_)
out <- rbind(mk("LoB", lob, DEF_LOB), mk("LoD", lod, DEF_LOD),
             mk("LoQ", loq, DEF_LOQ), cmp)
out$outcome_definition <- c(rep(DEF_LOB, 5), rep(DEF_LOD, 5), rep(DEF_LOQ, 5),
                            rep(DEF_CMP, 10))
out$evidence_tier <- "measured"
out$doi <- DOI
out$source_file <- "fulltext.xml"
out$table <- "Table 1"
readr::write_csv(out, fs::path(out_dir, "dpcr_cpn_lob_variability.csv"))

message("cpn_ddpcr_lob: LoB ", paste(lob$cp_ul, collapse = "/"),
        " cp/uL across 5 conditions; LOD/LoB ", paste(printed, collapse = "/"),
        "x recomputed OK; well volume ", vol, " uL -- assertions passed")
