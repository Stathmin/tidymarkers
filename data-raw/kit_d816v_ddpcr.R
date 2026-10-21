# =============================================================================
#  data-raw/kit_d816v_ddpcr.R  --  LoB / LoD / LoQ for a multiplex KIT D816V assay
# =============================================================================
#  Source (DOI is PROVENANCE; every number is re-read from the fetched file):
#
#    Development of a High-Sensitivity Multiplex mRNA-Based KIT D816V Droplet
#    Digital PCR Assay and Correlation with Tumor Load in Systemic Mastocytosis.
#    Int J Mol Sci 2026.  DOI 10.3390/ijms27146314.
#    Fetched: .tmp/lit/PMC13410330/fulltext.xml (Europe PMC fullTextXML).
#
#  WHY: the AML IDH file (data-raw/aml_idh_ddpcr.R) reports LoB and LoD only.
#  This paper reports LoQ as well, on THREE DIFFERENT NORMALISATION BASES
#  (genomic DNA VAF %, and two cDNA ratios, KIT/ABL1 and KIT/EAB), which is the
#  clearest available demonstration that a dPCR "limit" is an ASSAY PARAMETER
#  carried on a specific normalisation, not a portable sensitivity.
#
#  WHAT THIS IS NOT: there is no qPCR or Sanger head-to-head in Table 1, so no
#  conversion rate is derived.  The CV table (Table 2) is kept only as the
#  context in which LoQ was set; no CV is turned into an error probability.
#
#  ASSERTION: LoB <= LoD <= LoQ in every one of the three columns.  That is the
#  only relation the paper's own numbers must satisfy, and it is asserted here
#  so a re-parse that mixes columns (the three bases have similar magnitudes)
#  fails loudly.  Two of the three columns have LoQ == LoD exactly, which is
#  recorded verbatim rather than "corrected".
#
#  Run:  Rscript data-raw/kit_d816v_ddpcr.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

xml_path <- fs::path(root, ".tmp", "lit", "PMC13410330", "fulltext.xml")
if (!file.exists(xml_path)) stop("Missing ", xml_path)

DOI <- "10.3390/ijms27146314"
out_dir <- fs::path(root, "inst", "extdata")
doc <- xml2::read_xml(xml_path)
if (!any(grepl(DOI, xml2::xml_text(doc), fixed = TRUE)))
  stop("fulltext.xml does not mention ", DOI)

tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
t1 <- Filter(function(t) grepl("Analytical sensitivity and specificity",
                               xml2::xml_text(xml2::xml_find_first(
                                 t, ".//*[local-name()='caption']"))), tw)
if (length(t1) != 1L) stop("Table 1 (analytical sensitivity) not found exactly once")
rows <- lapply(xml2::xml_find_all(t1[[1L]], ".//*[local-name()='tr']"), function(r)
  trimws(xml2::xml_text(xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']"))))

grab <- function(prefix) {
  hit <- Filter(function(v) length(v) == 4L && grepl(paste0("^", prefix), v[[1L]]), rows)
  if (length(hit) != 1L) stop("Row '", prefix, "' not found exactly once")
  num <- function(x) as.numeric(sub("\\s.*$", "", x))
  c(gDNA_VAF = num(hit[[1L]][[2L]]), cDNA_KIT_ABL1 = num(hit[[1L]][[3L]]),
    cDNA_KIT_EAB = num(hit[[1L]][[4L]]))
}
lob <- grab("Limit of blank")
lod <- grab("Limit of detection")
loq <- grab("Limit of quantification")
if (!identical(names(lob), names(lod)) || !identical(names(lod), names(loq)))
  stop("Table 1 columns did not parse consistently")
if (any(lob > lod) || any(lod > loq))
  stop("Table 1 violates LoB <= LoD <= LoQ: LoB ", paste(lob, collapse = "/"),
       "; LoD ", paste(lod, collapse = "/"), "; LoQ ", paste(loq, collapse = "/"))

basis <- c(gDNA_VAF = "genomic DNA, variant allele fraction (%)",
           cDNA_KIT_ABL1 = "cDNA, KIT / ABL1 ratio (%)",
           cDNA_KIT_EAB = "cDNA, KIT / EAB ratio (%)")
DEF <- c(
  LoB = paste("Limit of blank: the signal level blank (HMC-1.2-free) samples reach.",
              "One number per normalisation basis; not a per-droplet probability."),
  LoD = paste("Limit of detection: CLSI-style lowest level detected above LoB.",
              "Assay- and normalisation-specific; not transferable to another primer",
              "set, amplicon or matrix."),
  LoQ = paste("Limit of quantification: lowest level with acceptable replicate CV.",
              "Distinct from LoD even where numerically equal; a value equal to LoD",
              "means the paper could not demonstrate an LoQ above LoD."))
out <- do.call(rbind, lapply(c("LoB", "LoD", "LoQ"), function(k) {
  v <- get(tolower(k))
  tibble::tibble(limit_type = k, normalisation = unname(basis[names(v)]),
                 value_pct = unname(v), assay = "KIT D816V multiplex mRNA assay",
                 equal_to_lod = if (k == "LoQ") unname(v == lod) else NA,
                 outcome_definition = DEF[[k]], evidence_tier = "measured",
                 doi = DOI, source_file = "fulltext.xml", table = "Table 1")
}))
readr::write_csv(out, fs::path(out_dir, "dpcr_kit_d816v_limits.csv"))

message("kit_d816v_ddpcr: LoB ", paste(lob, collapse = "/"), "; LoD ",
        paste(lod, collapse = "/"), "; LoQ ", paste(loq, collapse = "/"), " (",
        paste(names(basis), collapse = "/"), "); LoQ==LoD in ",
        sum(loq == lod), "/3 columns -- assertions passed")
