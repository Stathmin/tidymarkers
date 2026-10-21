# =============================================================================
#  data-raw/lamp_comp2018.R  --  competitive SNP-LAMP probes: an ENDPOINT
#                                FLUORESCENCE separation, and per-experiment
#                                LAMP durations read from the supplement PDF
# =============================================================================
#  Wu/… "Competitive SNP-LAMP probes for rapid and robust SNP detection",
#  Cell Rep Methods 2022. doi:10.1016/j.crmeth.2022.100242  PMID 35880021
#  PMC9308130 (OA). DOI/PMID read from the fetched JATS.
#
#  WHY THIS FILE EXISTS
#  This is the only source in the LAMP family that reports an allele gap as a
#  continuous SIGNAL, not a call and not a time. It matters for the proposal
#  because it proves the gap is real and large (2.7x) -- and because that gap
#  is measured on a STRAND-DISPLACEMENT PROBE DUPLEX after amplification, i.e.
#  it belongs to a case that ALSO carried six of the sequenced samples.
#
#  WHAT IT DOES NOT GIVE
#  * No time-to-positive and no Ct. The discrimination is endpoint fluorescence
#    after a fixed LAMP duration (Table S3: 60 or 75 min, one experiment 30).
#    So this paper CANNOT calibrate a kinetics model; it calibrates a
#    probe-endpoint model, which is a different object.
#  * The 2.7x is from synthetic DNA oligos (Figure 1B), not genomic samples.
#
#  ASSERTIONS (any failure -> stop(), nothing written)
#   (C1) Table S3 in mmc1.txt parses to exactly 7 experiment rows; LAMP
#        durations are a subset of {30, 60, 75} min and the 'Anneal only' rows
#        carry no LAMP duration.
#   (C2) Every headline number below occurs verbatim in text.txt:
#        2.7 times, 15%, p = 1.6e-7, p = 0.46, p = 8.3e-5, 10,000, 0.45%,
#        313.2, 47.6, 85.5%, 7 GA generations, population 128, and the
#        ~1.1 signal:background limitation.
#   (C3) The one-pot inhibition statement is present, and the two-step
#        preference for sensitivity is present.
#
#  Run:  Rscript data-raw/lamp_comp2018.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- file.path(root, ".tmp", "lit", "lamp_comp2018")
txt_p <- file.path(lit, "text.txt")
s3_p <- file.path(lit, "sup", "mmc1.txt")
for (p in c(txt_p, s3_p)) if (!file.exists(p)) stop("Missing ", p)
out_dir <- file.path(root, "inst", "extdata")

DOI <- "10.1016/j.crmeth.2022.100242"
PMCID <- "PMC9308130"

txt <- base::readLines(txt_p, warn = FALSE)
norm <- function(x) {
  # Fold non-breaking spaces FIRST: the JATS text uses them inside
  # "60–75 min" and TRE's [[:space:]] does not match U+00A0.
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  gsub("[[:space:]]+", " ", x)
}
norm_all <- norm(paste(txt, collapse = " "))
problems <- character(0)
need <- function(s, what) {
  # `<<-`, not `<-`: a local assignment here would silently swallow every
  # verbatim check (the function's own scope, not the script's).
  if (!grepl(norm(s), norm_all, fixed = TRUE)) {
    problems <<- c(problems, paste0("verbatim text not found (", what, "): ",
                                    substr(norm(s), 1, 70)))
  }
  invisible(TRUE)
}
need("DOI: 10.1016/j.crmeth.2022.100242", "doi")   # the dump's own header line
# The dump's own PMCID line reads 'PMCID: NA' (EPMC does not emit it in the
# JATS); PMC9308130 was taken from the fetch URL and is not asserted from text.
xml_p <- file.path(lit, "fulltext.xml")
if (file.exists(xml_p)) {
  x <- paste(base::readLines(xml_p, warn = FALSE, n = 200L), collapse = " ")
  if (!grepl(DOI, x, fixed = TRUE)) problems <- c(problems, "DOI not in fulltext.xml")
}

## -- C1: Table S3 --------------------------------------------------------------
s3 <- base::readLines(s3_p, warn = FALSE)
i <- grep("Table S3", s3)
if (length(i) != 1L) stop("lamp_comp2018: expected exactly one Table S3 in mmc1.txt")
# The table is column-aligned text; the experiment label is the text left of the
# first run of 2+ spaces after the heading, and the duration is the last field.
# Rows are taken until the footnote line starting with '*'.
end <- grep("^\\*Standard LAMP primer", s3)[1L]
if (is.na(end)) stop("lamp_comp2018: Table S3 footnote not found")
body <- s3[(i + 1L):(end - 1L)]
body <- body[nzchar(trimws(body))]
body <- body[!grepl("^\\s*LAMP\\s*$|^Experiment|^\\s*Duration\\s*$", body)]
exp_lab <- c("ACTB Probe", "MT-CO2 Probe", "NOTCH1 Probe", "TP53 Probe",
             "NOTCH1 V.A.F.", "TP53 V.A.F.",
             "Spike D614G", "Spike 614D and")
fmt <- vapply(body, function(l) {
  if (grepl("One-pot", l)) "One-pot" else if (grepl("Two-step", l)) "Two-step" else
    if (grepl("Anneal only", l)) "Anneal only" else NA_character_
}, "")
keep <- !is.na(fmt)   # data rows only; the 18-line block contains wrapped labels
fmt <- fmt[keep]
rowtxt <- trimws(body[keep])
dur_v <- vapply(rowtxt, function(l) {
  m <- regmatches(l, regexpr("[0-9]+ minutes", l))
  if (length(m)) as.integer(sub(" minutes", "", m)) else NA_integer_
}, 0L)
s3_tab <- data.frame(raw = rowtxt, format = unname(fmt), minutes = unname(dur_v),
                     stringsAsFactors = FALSE)
n_exp <- nrow(s3_tab)
if (!identical(n_exp, 8L)) {
  problems <- c(problems, paste("Table S3 parsed", n_exp, "data rows, expected 8",
                                "(ACTB + 3 one-pot + 2 V.A.F. + 2 SARS-CoV-2)"))
}
if (any(!is.na(s3_tab$minutes) & !s3_tab$minutes %in% c(30L, 60L, 75L))) {
  problems <- c(problems, "a Table S3 LAMP duration is not 30/60/75 min")
}
if (any(s3_tab$format == "Anneal only" & !is.na(s3_tab$minutes))) {
  problems <- c(problems, "an 'Anneal only' Table S3 row carries a LAMP duration")
}
if (any(s3_tab$format == "One-pot" & is.na(s3_tab$minutes))) {
  problems <- c(problems, "an 'One-pot' Table S3 row has no LAMP duration")
}
# The paper's own prose: 60-75 min -> 30 min when the faster primer set was
# found; the two-step SARS-CoV-2 experiment is the 30-min one.
need("This reduced the LAMP reaction time required from 60–75 min to 30 min.", "time reduction")

## -- C2: headline numbers ------------------------------------------------------
need("The fluorescence signal of the SNP sample was 2.7 times higher than the non-SNP sample.", "2.7x")
need("the non-SNP signal was only 15% higher than the non-template control", "15%")
need("p = 1.6 × 10−7", "p 1.6e-7")
need("no significant difference in signal between the two (p = 0.46)", "p 0.46")
need("Adding a sink complex greatly reduced the signal from a representative WT oligo target (p = 8.3 × 10−5)", "p 8.3e-5")
need("we simulated 10,000 random probe sets", "10000")
need("only 0.45% of these random designs displaying a signal greater than 70% and specificity more than 100-fold", "0.45%")
need("a median of 313.2 and a standard deviation (SD) of 47.6", "313.2/47.6")
need("produced 85.5% of its maximum possible fluorescent signal", "85.5%")
need("With only 7 GA generations and an initial population size of 128", "GA size")
need("the signal:background ratio is only ∼1.1", "1.1 signal:background")
need("only include up to 100 nM of each probe and sink strand", "100 nM cap")

## -- C3: one-pot inhibition ----------------------------------------------------
need("we have previously observed slight LAMP inhibition from the probe in one-pot reactions", "one-pot inhibition")
need("Ultimately, we found two-step SNP-LAMP assays to be more powerful, since 10-fold higher probe concentrations can be used.", "two-step better")

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("lamp_comp2018: assertion failed; nothing written")
}

## -- write --------------------------------------------------------------------
outcome <- paste(
  "endpoint fluorescence ratio after a fixed-duration LAMP reaction;",
  "NO time-to-positive and no Ct are reported, so this is a probe-endpoint",
  "separation, not a kinetics gap")

sep <- data.frame(
  comparison = c("SNP vs non-SNP (synthetic oligos, ACTB c.776A>C)",
                 "non-SNP vs no-template control",
                 "SNP vs non-SNP, WITHOUT the sink complex",
                 "SNP signal after sink vs before sink"),
  reported = c("2.7 times higher", "15% higher", "no significant difference",
               "no significant change"),
  magnitude = c(2.7, 1.15, NA, NA),
  p_value = c("1.6e-7", NA, "0.46", "0.65"),
  signal_type = "endpoint fluorescence",
  stringsAsFactors = FALSE)
sep$location <- "text.txt :: Results, Figure 1B"
sep$outcome_definition <- outcome
sep$evidence_tier <- "moderate"
sep$source_file <- "text.txt :: Results (Figure 1B)"
sep$doi <- DOI
sep$pmcid <- PMCID

design <- data.frame(
  item = c("random probe sets simulated",
           "fraction reaching >70% signal and >100x specificity",
           "median SNP/WT specificity over 10 GA runs",
           "SD of that specificity",
           "fraction of maximum fluorescence achieved",
           "GA generations",
           "GA initial population",
           "one-pot probe+sink concentration cap",
           "signal:background at the cap"),
  value = c(10000, 0.45, 313.2, 47.6, 85.5, 7, 128, 100, 1.1),
  unit = c("sets", "percent", "fold", "fold", "percent", "generations", "individuals",
           "nM each strand", "ratio"),
  stringsAsFactors = FALSE)
design$outcome_definition <- paste(
  "design-space and optimisation statistics for a PROBE-based SNP-LAMP;",
  "none of these is an allele gap in time")
design$evidence_tier <- "moderate"
design$source_file <- "text.txt :: Results + STAR Methods"
design$doi <- DOI
design$pmcid <- PMCID

dur <- s3_tab
if (!identical(nrow(dur), length(exp_lab))) {
  problems <- c(problems, "the Table S3 row count does not match the 8 known experiments")
}
dur$experiment <- c(exp_lab, rep(NA_character_, nrow(dur) - length(exp_lab)))[
  seq_len(nrow(dur))]
if (any(vapply(exp_lab, function(e) !any(startsWith(dur$raw, e)), NA))) {
  problems <- c(problems, "a known Table S3 experiment label was not found in its row")
}
dur$outcome_definition <- paste(
  "fixed LAMP duration per experiment (Table S3), from the supplement PDF's",
  "text layer; NOT a time-to-positive")
dur$evidence_tier <- "moderate"
dur$source_file <- "sup/mmc1.txt (pdftotext -layout) :: Table S3"
dur$doi <- DOI
dur$pmcid <- PMCID

write.csv(sep[, c("comparison", "reported", "magnitude", "p_value", "signal_type",
                  "location", "outcome_definition", "evidence_tier", "source_file",
                  "doi", "pmcid")],
          file.path(out_dir, "lamp_comp2018_separation.csv"), row.names = FALSE, na = "")
write.csv(design[, c("item", "value", "unit", "outcome_definition", "evidence_tier",
                     "source_file", "doi", "pmcid")],
          file.path(out_dir, "lamp_comp2018_design_space.csv"), row.names = FALSE, na = "")
write.csv(dur[, c("raw", "format", "minutes", "experiment", "outcome_definition",
                  "evidence_tier", "source_file", "doi", "pmcid")],
          file.path(out_dir, "lamp_comp2018_durations.csv"), row.names = FALSE, na = "")

message(
  "lamp_comp2018: ", nrow(sep), " separation rows (2.7x endpoint fluorescence,",
  " no kinetics); ", nrow(design), " design-space rows (0.45% of 10,000; ",
  "specificity 313.2 +/- 47.6); Table S3 ", nrow(dur), " experiment rows, ",
  "durations ", paste(sort(unique(stats::na.omit(dur$minutes))), collapse = "/"),
  " min; one-pot inhibition recorded"
)
