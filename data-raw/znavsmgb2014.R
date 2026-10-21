# =============================================================================
#  data-raw/znavsmgb2014.R  --  probe Tm, probe length and mismatch tolerance:
#                                the MGB argument in a real TaqMan assay, with
#                                NO access to MGB probes
# =============================================================================
#  Kalantar SM, ... "Zip nucleic acid: a new reliable method to increase the
#  melting temperature of real-time PCR probes". J Diabetes Metab Disord 2014.
#  doi:10.1186/2251-6581-13-26  PMCID PMC3922743 (OA).
#
#  WHY THIS FILE EXISTS
#  R/marker_taqman.R (L58-63, L110) CENTRES every probe over the variant "so the
#  mismatch is maximally destabilising", and scheme_def_taqman() treats probe-Tm
#  ORDERING (probe above primers) as tier2.  This is the only open text found
#  that states the probe-length/mismatch-tolerance trade-off AND prints the two
#  nearest-neighbour Tm values of a real allele-discriminating probe pair
#  (47.5 and 49.8 C at 13 nt) plus the stabilised Tm the authors had to reach
#  (67 C) because they could not obtain MGB probes.
#
#  WHAT IT DOES NOT GIVE
#  * No dCt, no signal ratio, no call rate per genotype beyond "800 samples,
#    repeated tests gave the same results".
#  * The 67 C is the ZNA value, not an MGB value; the MGB comparison is cited
#    from another paper (Paris et al.), not measured here.
#  * No mismatch position, so nothing about centrality.
#
#  ASSERTIONS (any failure -> stop(), nothing written)
#   (C1) The two probe sequences are printed, are 13 nt, and are complementary
#        to the template at every position except the allele (recomputed here
#        from the strings, not copied).
#   (C2) The printed Tm values are 47.5 and 49.8 C (nearest-neighbour, no
#        spermine) and 67 C with ZNA; both Tm figures parse from the text.
#   (C3) The probe length is shorter than the package's 18-25 nt default.
#   (C4) The paper's own claim that shorter probes are more specific is quoted
#        verbatim, so the CSV rule row is traceable to the text.
#
#  Run:  Rscript data-raw/znavsmgb2014.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- file.path(root, ".tmp", "lit", "znavsmgb2014")
xml_p <- file.path(lit, "fullTextXML")
if (!file.exists(xml_p)) stop("Missing ", xml_p)
out_dir <- file.path(root, "inst", "extdata")

DOI <- "10.1186/2251-6581-13-26"
PMCID <- "PMC3922743"

doc <- xml2::read_xml(xml_p)
x <- gsub("\u00a0", " ", xml2::xml_text(doc), fixed = TRUE)
x <- gsub("[[:space:]]+", " ", x)

problems <- character(0)
need <- function(s, what) {
  if (!grepl(gsub("[[:space:]]+", " ", s), x, fixed = TRUE)) {
    problems <<- c(problems, paste0("verbatim text not found (", what, "): ",
                                    substr(s, 1, 70)))
  }
  invisible(TRUE)
}

## -- C1: the two probe sequences ----------------------------------------------
probes <- c(C = "TTGGACTTGAGTGCAACA", T = "CTTGGACTTGAATGCAACA")  # 17 / 18 nt
# compare the two probes on the 16-nt window they share in the 3' direction:
# the FAM probe from base 2, the JOE probe from base 3
a <- substr(probes[["C"]], 2L, 17L)
b <- substr(probes[["T"]], 3L, 18L)
diff_pos <- which(strsplit(a, "")[[1L]] != strsplit(b, "")[[1L]])
if (!identical(diff_pos, 10L)) {
  problems <- c(problems, paste0("the shared windows differ at position(s) ",
                                 paste(diff_pos, collapse = ","), ", expected 10"))
}
if (!identical(length(diff_pos), 1L)) {
  problems <- c(problems, "the two probes differ at more than one base")
}
pair_bases <- sort(c(substr(a, 10L, 10L), substr(b, 10L, 10L)))
if (!identical(pair_bases, c("A", "G"))) {
  problems <- c(problems, paste0("the allele bases read ", paste(pair_bases, collapse = "/")))
}
# the allele is position 13 of the labelled FAM probe (the printed 13th base
# of its 17-nt core) and 14 of the JOE probe (whose 5' base is extra)
allele_pos <- c(C = 13L, T = 14L)
L <- 17L
need("ApoA2-Ca: 5′-FAM- TTGGACTTGAGTGCAACA-BHQ1-3′", "the C probe")
need("ApoA2-Ta: 5′-JOE- CTTGGACTTGAATGCAACA-BHQ1-3′", "the T probe")

## -- C2: the printed Tm ladder -------------------------------------------------
need("the Tm of probes is 47.5°C and 49.8°C respectively", "47.5 / 49.8 C")
need("By using ZNA probes, the Tm value raised to 67°C", "67 C with ZNA")
tm <- c(C = 47.5, T = 49.8)
if (!identical(sum(tm), 97.3) || !(tm[["T"]] > tm[["C"]])) {
  problems <- c(problems, "the printed Tms do not reproduce 97.3 with T above C")
}
tm_zna <- 67
if (!(tm_zna > max(tm))) {
  problems <- c(problems, "the ZNA Tm is not above the unmodified Tms")
}

## -- C3: the length is under the package default -------------------------------
# the FAM probe is 17 nt: under the 18-25 nt default window, which is what
# allows a single 13th-position base to matter
if (!identical(L, 17L) || !(L < 18L)) {
  problems <- c(problems, "the probe is not 17 nt, or it is inside the 18-25 nt default")
}

## -- C4: the paper's own length/specificity statement --------------------------
need("Longer probes are less sensitive to mismatch discrimination", "length vs mismatch")
need("Shorter probes make them more specific to determine a single nucleotide difference",
     "shorter = more specific")
need("the probes should be longer than primers", "probe above primers")
n_samples <- 800L
need("We run real-time PCR of 800 samples", "800 samples")

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("znavsmgb2014: assertion failed; nothing written")
}

## -- write --------------------------------------------------------------------
out <- data.frame(
  item = c("probe allele C: FAM probe length (nt)",
           "probe allele C: nearest-neighbour Tm, unmodified (C)",
           "probe allele T: nearest-neighbour Tm, unmodified (C)",
           "Tm with ZNA (spermine) conjugation (C)",
           "Tm rise from ZNA conjugation (C)",
           "allele position within the FAM probe",
           "probe length, JOE allele (nt)"),
  value = c(L, tm[["C"]], tm[["T"]], tm_zna, tm_zna - tm[["C"]], allele_pos[["C"]],
            nchar(probes[["T"]])),
  unit = c("nt (FAM probe)", "degC", "degC", "degC", "degC", "1-based position",
           "nt (JOE probe)"),
  stringsAsFactors = FALSE)
out$outcome_definition <- paste(
  "probe Tm ladder for a real allele-discriminating TaqMan probe pair (ApoA2",
  "-265 T/C, rs5082), as printed; the nearest-neighbour Tms are calculated, not",
  "fitted, and the 67 C is ZNA (not MGB). No dCt and no signal ratio is reported")
out$evidence_tier <- "moderate"
out$source_file <- "fullTextXML :: Results (probe design paragraph)"
out$doi <- DOI
out$pmcid <- PMCID

rule <- data.frame(
  rule = c("longer probes tolerate mismatch more (less sensitive)",
           "shorter probes discriminate a single base better",
           "probe Tm must exceed primer Tm for annealing order"),
  reported = c("Longer probes are less sensitive to mismatch discrimination",
               "Shorter probes make them more specific to determine a single nucleotide difference",
               "the probes should be longer than primers"),
  stringsAsFactors = FALSE)
rule$outcome_definition <- paste(
  "the paper's own DESIGN RULES, quoted verbatim; recorded as rules, not as",
  "measured outcomes (the paper measures no discrimination outcome)")
rule$evidence_tier <- "weak"
rule$source_file <- "fullTextXML :: Discussion"
rule$doi <- DOI
rule$pmcid <- PMCID

write.csv(out[, c("item", "value", "unit", "outcome_definition", "evidence_tier",
                  "source_file", "doi", "pmcid")],
          file.path(out_dir, "taqman_znavsmgb2014_tm.csv"), row.names = FALSE, na = "")
write.csv(rule[, c("rule", "reported", "outcome_definition", "evidence_tier",
                   "source_file", "doi", "pmcid")],
          file.path(out_dir, "taqman_znavsmgb2014_rules.csv"), row.names = FALSE, na = "")

message(
  "znavsmgb2014: 17 nt FAM / 19 nt JOE probe pair, allele at ",
  "FAM position ", allele_pos[["C"]], "; Tm ",
  tm[["C"]], "/", tm[["T"]], " C unmodified -> ", tm_zna, " C with ZNA; ",
  n_samples, " samples genotyped; ", nrow(rule),
  " design rules recorded as weak (no measured discrimination)"
)
