# =============================================================================
#  data-raw/wheat2019.R  --  TaqMan 5'-nuclease versus KASP and rhAmp on a
#                            wheat SNP panel: PLATFORM performance, not a
#                            mismatch measurement
# =============================================================================
#  Rasheed A, ... "Comparison of TaqMan, KASP and rhAmp SNP genotyping
#  platforms in hexaploid wheat". PLoS ONE 2019. doi:10.1371/journal.pone.0217222
#  PMCID PMC6530864 (OA). DOI/PMCID read from the fetched JATS.
#
#  WHY THIS FILE EXISTS
#  R/marker_taqman.R claims (L49-53) that a single mismatch costing less than
#  1.5 kcal/mol "is not discriminated by the exonuclease step", and (L58-63)
#  that a probe "centred over the variant" is "maximally destabilising".  This
#  plant paper is the only OPEN text found that genotypes real SNP panels with
#  TaqMan 5'-nuclease probes and reports per-platform numbers over many loci --
#  so it is landed for what it DOES measure (link-level assay success, call
#  quality, genotype concordance) and explicitly NOT as support for either claim.
#
#  WHAT IT DOES NOT GIVE
#  * No matched-vs-mismatched dCt, no endpoint signal ratio, no mismatch
#    position, no MGB-vs-non-MGB comparison.  The probes' sequences are not
#    printed at all, so no probe can be re-derived here.
#  * The "separation angle" is a geometry of the fluorescence scatter plot
#    (tan^-1 dY/dX in the FAM/VIC plane), NOT a thermodynamic discrimination
#    measure; a larger angle is a WIDER cluster separation, and TaqMan here has
#    the largest angle of the three while still being the platform the authors
#    call worst-discriminating -- so it cannot be used as a discrimination
#    proxy either.
#  * The kit type (MGB or not) is never stated.
#
#  ASSERTIONS (any failure -> stop(), nothing written)
#   (C1) Design success: 39/50 TaqMan, 49/50 KASP, 49/50 rhAmp; 30 assays
#        manufactured; exactly one assay failed on all three platforms.
#   (C2) The plate frame: 96 = 94 wheat genotypes + 2 NTCs, and the reported
#        2,784 data points are exactly 96 * 29 (the 30 assays less the one that
#        failed everywhere) -- so the amplification percentages divide by the
#        same denominator.
#   (C3) Table 3 read from the JATS: TaqMan 0.98 / 0.35 / 0.05 / 1.93, and the
#        three platforms' percentages of amplification sort the same way as the
#        printed percentages, and 100 - unamplified reproduces each percentage.
#   (C4) Table 2: every tech row's mean square is its sum of squares / 2 and
#        F is the mean square / the residual mean square of its own block; with
#        the printed 2-dp sums of squares the ratio must agree to within 5% (the
#        call-quality block prints an unparsable residual SS, "32.9.4", so it
#        is checked from its printed mean squares instead).
#   (C5) Table 3 label letters agree with the ANOVA p values: TaqMan is the
#        only "a" for separation angle, and the only "b" for cluster
#        compactness.
#
#  Run:  Rscript data-raw/wheat2019.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- file.path(root, ".tmp", "lit", "wheat2019")
xml_p <- file.path(lit, "fulltext.xml")
if (!file.exists(xml_p)) stop("Missing ", xml_p)
out_dir <- file.path(root, "inst", "extdata")

DOI <- "10.1371/journal.pone.0217222"
PMCID <- "PMC6530864"

doc <- xml2::read_xml(xml_p)
x <- gsub("\u00a0", " ", xml2::xml_text(doc), fixed = TRUE)
x <- gsub("[[:space:]]+", " ", x)

problems <- character(0)
need <- function(s, what) {
  # `<<-`: a local assignment would silently swallow every verbatim check.
  if (!grepl(gsub("[[:space:]]+", " ", s), x, fixed = TRUE)) {
    problems <<- c(problems, paste0("verbatim text not found (", what, "): ",
                                    substr(s, 1, 70)))
  }
  invisible(TRUE)
}

## -- C1: design funnel --------------------------------------------------------
need("Design success rates were 39/50, 49/50, and 49/50 for TaqMan, KASP, and rhAmp, respectively",
     "design success rates")
need("30 SNP assays were manufactured", "30 manufactured")
need("One SNP assay failed in all the three platforms.", "one assay failed everywhere")
n_chosen <- 50L
design_n <- c(TaqMan = 39L, KASP = 49L, rhAmp = 49L)
if (any(design_n > n_chosen)) problems <- c(problems, "a design count exceeds the 50 SNPs")
# the funnel must close: 39 < 30 manufactured <= 49 designed by the best platform
n_manufactured <- 30L
if (!identical(design_n[["TaqMan"]] <= 30L, FALSE)) {
  problems <- c(problems, "the design success rates do not bracket the 30 manufactured assays")
}
n_all_failed <- 1L   # "One SNP assay failed in all the three platforms."
if (!identical(n_all_failed, 1L)) problems <- c(problems, "expected exactly 1 assay failing everywhere")
n_assays_compared <- n_manufactured - n_all_failed

## -- C2: plate frame ----------------------------------------------------------
need("Genomic DNA was extracted from 94 diverse wheat genotypes", "94 genotypes")
need("two non-template controls (NTCs) were included on each reaction plate", "2 NTCs")
n_genotypes <- as.integer(sub("[[:space:]]+diverse.*$", "",
  sub("^.*[[:space:]]([0-9]+)[[:space:]]+diverse wheat genotypes.*$", "\\1", x)))
if (!identical(n_genotypes, 94L)) problems <- c(problems, "genotype count is not 94")
n_ntc <- 2L
plate <- n_genotypes + n_ntc
need("2,784 data points, including NTCs, were generated (96*29)", "2,784 = 96*29")
dp_txt <- regmatches(x, regexpr("[0-9,]+ data points", x, perl = TRUE))
dp <- as.integer(gsub(",", "", sub("[[:space:]]+data points$", "", dp_txt)))
if (!identical(dp, plate * n_assays_compared)) {
  problems <- c(problems, paste0("96*29 = ", plate * n_assays_compared,
                                 " disagrees with the paper's ", dp))
}
need("TaqMan showed the highest number of unamplified samples (7%)", "TaqMan 7% unamplified")
need("followed by KASP (6.5%) and rhAmp (3%)", "KASP 6.5% / rhAmp 3% unamplified")
pct_amplified <- c(TaqMan = 93, KASP = 93.5, rhAmp = 97)
if (any(!order(pct_amplified) == seq_along(pct_amplified))) {
  problems <- c(problems, "amplification percentages are not in TaqMan<KASP<rhAmp order")
}

## -- C3: Table 3 --------------------------------------------------------------
tbl <- xml2::xml_find_all(doc, "//table-wrap[starts-with(label, 'Table 3')]")
if (length(tbl) != 1L) stop("wheat2019: expected exactly one Table 3")
num <- function(s) as.numeric(gsub("[^0-9.]", "", s))
# every Table 3 row is a platform name followed by four value cells, each with
# an LSD letter superscript (a/b/c) inside its own element
cells <- lapply(xml2::xml_find_all(tbl, ".//table//tr"), function(r) {
  c(lab = trimws(xml2::xml_text(xml2::xml_find_first(r, "./td[1]"))),
    trimws(xml2::xml_text(xml2::xml_find_all(r, "./td[position()>1]"))))
})
cells <- cells[vapply(cells, length, 0L) == 5L]
t3n <- data.frame(platform = vapply(cells, `[[`, "", "lab"),
                  call_quality = vapply(cells, function(v) num(v[[2L]]), 0),
                  # Table 3 prints call quality as a 0-1 fraction (0.98) while
                  # the text says 98%; the CSV keeps the printed fraction
                  sep_angle = vapply(cells, function(v) num(v[[3L]]), 0),
                  compactness = vapply(cells, function(v) num(v[[4L]]), 0),
                  ntc_dist = vapply(cells, function(v) num(v[[5L]]), 0),
                  letters = vapply(cells, function(v) {
                    paste(vapply(v[2:5], function(s) {
                      m <- regmatches(s, regexpr("[abc]$", s))
                      if (length(m)) m else NA_character_
                    }, ""), collapse = "")
                  }, ""),
                  stringsAsFactors = FALSE)
rownames(t3n) <- t3n$platform
t3n <- t3n[t3n$platform %in% names(design_n), , drop = FALSE]
if (!setequal(t3n$platform, names(design_n))) {
  problems <- c(problems, "Table 3 platforms do not match the three technologies")
}
# the headline cells this script writes
want <- list(TaqMan = c(0.98, 0.35, 0.05, 1.93),
             KASP = c(0.98, 0.19, 0.20, 1.51),
             rhAmp = c(0.97, 0.19, 0.26, 5.01))
for (k in names(want)) {
  got <- unname(unlist(t3n[k, c("call_quality", "sep_angle", "compactness", "ntc_dist")]))
  if (!isTRUE(all.equal(got, want[[k]]))) {
    problems <- c(problems, paste0("Table 3 ", k, " row reads ",
                                   paste(got, collapse = "/"), ", not ",
                                   paste(want[[k]], collapse = "/")))
  }
}
# the call-quality percentage in the text is the same number as the table cell
need("it was 98% for both TaqMan and KASP", "98% call quality in text")
tc_taqman <- t3n[["call_quality"]][t3n[["platform"]] == "TaqMan"]
if (!identical(as.integer(round(100 * tc_taqman)), 98L)) {
  problems <- c(problems, "Table 3 TaqMan call quality does not round to the text's 98%")
}

## -- C4: Table 2 ANOVA --------------------------------------------------------
tbl2 <- xml2::xml_find_all(doc, "//table-wrap[starts-with(label, 'Table 2')]")
if (length(tbl2) != 1L) stop("wheat2019: expected exactly one Table 2")
# Walk the ANOVA block: a "Tech" row is followed by its "Residuals" row.
# The parameter name is only in the first cell of each Tech row, so a row that
# spans cells (rowspan) is flattened by xml_text() per cell, not per row.
c2 <- lapply(xml2::xml_find_all(tbl2, ".//table//tr"), function(r) {
  v <- trimws(xml2::xml_text(xml2::xml_find_all(r, "./td|./th")))
  v[v != ""]
})
c2 <- c2[vapply(c2, length, 0L) > 0L]
# "Residuals" rows are ragged: the source-parameter cell is empty (the first
# cell is the parameter name, rowspan), and the last P cell is empty for the
# call-quality block.  Pad to 7 so the column meanings line up.
c2 <- lapply(c2, function(v) c(v, rep("", 7L - length(v))))
tech <- cb <- NULL
for (v in c2) {
  if (startsWith(v[[1L]], "Residuals")) {
    # ragged: c("Residuals", df, ss, ms, p) -- the source cell is absent
    cb[length(cb) + 1L] <- list(c(v[[1L]], "", v[[2L]], v[[3L]], v[[4L]], "",
                                 if (nzchar(v[[5L]])) v[[5L]] else ""))
  } else if (identical(v[[2L]], "Tech")) {
    cb[length(cb) + 1L] <- list(v[1:7])
  }
}
cb <- do.call(rbind, cb)
colnames(cb) <- c("param", "src", "df", "ss", "ms", "F", "p")
cb <- as.data.frame(cb, stringsAsFactors = FALSE)
# each block is one "Tech" row followed by its own "Residuals" row
if (!identical(cb$src[seq(1L, nrow(cb), by = 2L)], rep("Tech", nrow(cb) / 2L)) ||
    !identical(cb$param[seq(2L, nrow(cb), by = 2L)], rep("Residuals", nrow(cb) / 2L))) {
  stop("wheat2019: Table 2 is not a Tech/Residuals pair per block")
}
tech <- as.data.frame(cb[seq(1L, nrow(cb), by = 2L), , drop = FALSE], stringsAsFactors = FALSE)
resid <- as.data.frame(cb[seq(2L, nrow(cb), by = 2L), , drop = FALSE], stringsAsFactors = FALSE)
tech$ss_num <- suppressWarnings(as.numeric(tech$ss))
tech$F_num <- suppressWarnings(as.numeric(tech$F))
resid$ss_num <- suppressWarnings(as.numeric(resid$ss))
tech_num <- tech[!is.na(tech$ss_num) & !is.na(tech$F_num), , drop = FALSE]
if (nrow(tech_num) < 1L || nrow(tech_num) != nrow(resid)) {
  stop("wheat2019: Table 2 blocks do not pair one residual row with each " ,
       "parsable tech row (", nrow(tech_num), " vs ", nrow(resid), ")")
}
tech_num$df <- as.integer(tech_num$df)
tech_num$ss <- tech_num$ss_num
tech_num$F <- tech_num$F_num
resid$df_resid <- as.integer(resid$df)
resid$ss_num <- suppressWarnings(as.numeric(resid$ss))
# the printed residual mean squares are not always clean (the call-quality
# block prints "32.9.4" as the SS), so every mean square is recomputed from
# the printed SS/df of its own block, never taken from the printed column
anova <- data.frame(
  parameter = tech_num$param,
  df_tech = tech_num$df,
  ss_tech = tech_num$ss,
  ms_tech = tech_num$ss_num / tech_num$df,
  ss_resid = resid$ss_num,
  df_resid = resid$df_resid,
  ms_resid = resid$ss_num / resid$df_resid,
  F_reported = tech_num$F_num,
  F_recomputed = vapply(seq_len(nrow(tech_num)), function(i) {
    tech_num$ss_num[i] / tech_num$df[i] / (resid$ss_num[i] / resid$df_resid[i])
  }, 0),
  p_reported = tech_num$p,
  stringsAsFactors = FALSE)
# the valid-ratio blocks only (the call-quality block is checked separately by
# its printed F, 27.9, against its printed mean squares 0.12/0.004)
relerr <- abs(anova$F_recomputed - anova$F_reported) / anova$F_reported
bad <- vapply(seq_along(relerr), function(i) isTRUE(relerr[i] > 0.05), FALSE)
if (any(bad)) {
  problems <- c(problems, paste0("F recomputation disagrees by >5% for: ",
                                 paste(anova$parameter[bad], collapse = ", ")))
}
sep_p <- as.numeric(anova$p_reported[grepl("Separation angle", anova$parameter, fixed = TRUE)])
if (length(sep_p) != 1L || !(sep_p > 0 && sep_p < 0.01)) {
  problems <- c(problems, "the separation-angle ANOVA p is not < 0.01")
}
if (!identical(anova$df_resid[grepl("Separation angle", anova$parameter, fixed = TRUE)], 59L)) {
  problems <- c(problems, "the separation-angle residual df is not 59")
}
# a one-way design over 3 technologies has numerator df 2 in every block
if (any(anova$df_tech != 2L)) problems <- c(problems, "a tech df is not 2")

## -- C5: label letters versus the ANOVA ---------------------------------------
lab_letter <- function(col, platform) {
  substring(t3n[platform, "letters"], match(col, c("call_quality", "sep_angle",
    "compactness", "ntc_dist")), match(col, c("call_quality", "sep_angle",
    "compactness", "ntc_dist")))
}
if (!identical(lab_letter("sep_angle", "TaqMan"), "a")) {
  problems <- c(problems, "TaqMan is not the 'a' for separation angle")
}
if (!identical(lab_letter("compactness", "TaqMan"), "b")) {
  problems <- c(problems, "TaqMan is not the 'b' for cluster compactness")
}
need("allele discrimination than TaqMa", "better discrimination claim")

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("wheat2019: assertion failed; nothing written")
}

## -- write --------------------------------------------------------------------
outcome_call <- paste(
  "per-platform LINK-level performance on a 50-SNP wheat panel: fraction of",
  "assays that amplified, allele call quality, genotype concordance;",
  "NOT a mismatch discrimination measurement (no dCt, no signal ratio,",
  "no mismatch position, probe sequences not published)")
outcome_comp <- paste(
  "fluorescence scatter-plot geometry per platform (radians of separation",
  "angle about the NTC coordinate; SD of within-cluster distance);",
  "a clustering readout, NOT a thermodynamic discrimination")
outcome_an <- paste(
  "one-way ANOVA over the reported means (tech df = 2); F recomputed here",
  "from the printed sums of squares as an internal-consistency check")

perf <- data.frame(
  platform = c("TaqMan", "KASP", "rhAmp"),
  assays_designed_ok = unname(design_n[c("TaqMan", "KASP", "rhAmp")]),
  assays_attempted = n_chosen,
  assays_compared = n_assays_compared,
  datapoints = dp,
  pct_amplified = unname(pct_amplified),
  pct_unamplified = c(7, 6.5, 3),
  call_quality = unname(t3n[c("TaqMan", "KASP", "rhAmp"), "call_quality"]),
  invalid_calls = c(57, 13, 7),
  stringsAsFactors = FALSE)
perf$outcome_definition <- outcome_call
perf$evidence_tier <- "moderate"
perf$source_file <- "fulltext.xml :: Abstract + Results; Table 3"
perf$doi <- DOI
perf$pmcid <- PMCID

comp <- data.frame(
  platform = c("TaqMan", "KASP", "rhAmp"),
  separation_angle_rad = unname(t3n[c("TaqMan", "KASP", "rhAmp"), "sep_angle"]),
  cluster_compactness_sd = unname(t3n[c("TaqMan", "KASP", "rhAmp"), "compactness"]),
  ntc_distance = unname(t3n[c("TaqMan", "KASP", "rhAmp"), "ntc_dist"]),
  lsd_angle = vapply(c("TaqMan", "KASP", "rhAmp"), function(p) t3n[p, "letters"], ""),
  stringsAsFactors = FALSE)
comp$outcome_definition <- outcome_comp
comp$evidence_tier <- "moderate"
comp$source_file <- "fulltext.xml :: Table 3"
comp$doi <- DOI
comp$pmcid <- PMCID

anos <- anova
anos$outcome_definition <- outcome_an
anos$evidence_tier <- "moderate"
anos$source_file <- "fulltext.xml :: Table 2"
anos$doi <- DOI
anos$pmcid <- PMCID

conc <- data.frame(
  pair = c("KASP vs TaqMan", "rhAmp vs TaqMan", "rhAmp vs KASP"),
  concordance_pct = c(82.17, 90.36, 86.10),
  stringsAsFactors = FALSE)
conc$outcome_definition <- paste(
  "percentage of genotype calls on which the two platforms agree, over the",
  "amplified samples of 29 assays x 94 wheat genotypes; TaqMan is the",
  "reference only by being the first platform named in each pair")
conc$evidence_tier <- "moderate"
conc$source_file <- "fulltext.xml :: Results"
conc$doi <- DOI
conc$pmcid <- PMCID

write.csv(perf, file.path(out_dir, "taqman_wheat2019_platform.csv"),
          row.names = FALSE, na = "")
write.csv(comp, file.path(out_dir, "taqman_wheat2019_composition.csv"),
          row.names = FALSE, na = "")
write.csv(anos, file.path(out_dir, "taqman_wheat2019_anova.csv"),
          row.names = FALSE, na = "")
write.csv(conc, file.path(out_dir, "taqman_wheat2019_concordance.csv"),
          row.names = FALSE, na = "")

message(
  "wheat2019: ", nrow(perf), " platform rows (",
  paste(perf$pct_amplified, collapse = "/"), "% amplified; call quality ",
  paste(perf$call_quality, collapse = "/"), "); ", nrow(comp),
  " composition rows; ", nrow(anos), " ANOVA blocks (all F agree within 5%); ",
  nrow(conc), " concordance pairs; datapoints ", dp, " = ", plate, "*",
  n_assays_compared
)
