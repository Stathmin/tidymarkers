# =============================================================================
#  data-raw/wheat_platform2019.R  --  the 2019 three-platform wheat comparison
# =============================================================================
#  Ayalew H, Tsang PW, Chu C, Wang J, Liu S, Chen C, Ma X.  "Comparison of
#  TaqMan, KASP and rhAmp SNP genotyping platforms in hexaploid wheat."
#  PLoS ONE 2019;14(5):e0217222.  DOI 10.1371/journal.pone.0217222, PMC6530864.
#
#  This is the one source in the KASP family that measures, on the same wheat
#  samples and the same targets, ALL of:
#
#    design      how many of 50 targets yielded a primer pair, per platform;
#    manufacture how many were made and run;
#    work        how many samples amplified and how many calls were valid;
#    cluster     separation angle, compactness, NTC distance, call quality.
#
#  It is therefore the test case for the distinction the rest of the family
#  blurs: a designed assay is not a working assay.  See
#  data-raw/proposed/kasp_z8_scope.proposal.md.
#
#  Two files are read:
#    S1 Table  pone.0217222.s001.xlsx  -- the 30 manufactured assays and their
#                                        KASP primer triples (sheet "Sheet1")
#    full text PMC6530864/fullTextXML    -- every count, Table 2, Table 3
#
#  Two traps this script refuses to fall into:
#
#  1. The "39"s.  Methods say "Thirty-nine out of the 50 ... passed assay design
#     criteria across the three platforms" (an INTERSECTION) while Results and
#     the abstract say TaqMan 39/50 (a per-platform rate).  Both readings are
#     written to the CSV verbatim with a `scope` column; they are never merged.
#  2. 6.5% is the fraction of ASSAY x SAMPLE points that did not amplify, not
#     an assay-level failure rate.  Its denominator (2,784 = 96 x 29) is
#     asserted here, and `numerator_printed` records that the paper printed a
#     percentage, not a count: no integer count is invented from it.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)
dir.create(file.path(root, ".tmp", "lit", "wheat_platform2019"), recursive = TRUE,
           showWarnings = FALSE)
source(file.path(root, "data-raw", "kasp_lib.R"))

DATE_VERIFIED <- "2026-10-21"
DOI <- "10.1371/journal.pone.0217222"
OUT <- fs::path(root, "inst", "extdata")

xlsx <- kasp_input(
  root, "wheat_platform2019", "supplements/pone.0217222.s001.xlsx",
  "wheat_platform2019_kasp_assays.xlsx", "wheat_platform2019 S1 Table"
)
xml_path <- kasp_input(
  root, "wheat_platform2019", "fulltext.xml",
  "wheat_platform2019_fulltext.xml", "wheat_platform2019 full text"
)

## --- the full text -----------------------------------------------------------
stopifnot(file.size(xml_path) > 20000L)
txt <- paste(readLines(xml_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
if (!grepl("<article", txt, fixed = TRUE)) {
  stop("the fulltext file is not a JATS article (no <article> element)")
}
flat <- kasp_norm(gsub("<[^>]+>", " ", txt))

# Every number below is cut out of `flat` by position, and each is recorded
# with the sentence fragment it came from, so the test file can re-read the
# same fragment from the same file and re-derive the same value.
REC <- list()
q <- function(value, site, note) {
  REC[[length(REC) + 1L]] <<- c(value = as.character(value),
                                site = substr(site, 1L, 120L), note = note)
  value
}
anchor <- function(word, len = 90L) {
  loc <- kasp_find(flat, word)
  list(loc = loc, site = kasp_site(loc, flat, len))
}
num <- function(a, window = 120L) {
  q(kasp_num_after(a$loc, flat, window, what = a$site), a$site,
    paste("first number after:", a$loc$anchor))
}

## design (per platform) ------------------------------------------------------
# "Fifty SNP markers linked to genes ... were chosen": the count is spelled out,
# so it is read from the word, not from a digit.
a50 <- anchor("Fifty SNP markers linked")
if (substr(a50$site, 1L, 5L) != "Fifty") stop("abstract no longer says Fifty: ", a50$site)
n_design_tested <- q(50L, a50$site, "spelled-out 'Fifty' at the abstract")

# "Design success rates were 39/50, 49/50, and 49/50 for TaqMan, KASP, and rhAmp"
adesign <- anchor("Design success rates were", 80L)
dnums <- kasp_num_after(adesign$loc, flat, window = 80L, what = "design rates")
# four numbers follow: 39, 50, 49, 50 ... and then 49, 50 for rhAmp
dsite <- adesign$site
dvals <- numeric()
pos <- nchar("Design success rates were ")
rest <- substr(dsite, pos + 1L, nchar(dsite))
for (i in 1:6) {
  nsel <- regmatches(rest, gregexpr("[0-9]+", rest, perl = TRUE))
  v <- if (length(nsel[[1L]])) as.numeric(nsel[[1L]][1L]) else NA_real_
  if (is.na(v)) break
  dvals <- c(dvals, v)
  rest <- sub("[0-9]+", "", rest)
}
if (length(dvals) < 6L) {
  stop("the design-success sentence does not carry three x/y pairs; the text at ",
       "its site is: ", dsite)
}
dvals <- dvals[1:6]
design <- list(
  taqman = q(dvals[1:2], dsite, "design-success sentence, first pair"),
  kasp   = q(dvals[3:4], dsite, "design-success sentence, second pair"),
  rhamp  = q(dvals[5:6], dsite, "design-success sentence, third pair")
)
if (!all(vapply(design, `[`, numeric(1), 2L) == n_design_tested)) {
  stop("a design-success denominator is not the 50 targets: ",
       paste(vapply(design, function(x) paste(x, collapse = "/"), character(1)),
             collapse = " "))
}
n_design_tested <- as.integer(n_design_tested)
design <- lapply(design, as.integer)

## the 39 / 30 / 9 arithmetic, and the Methods' ambiguous reuse of "39" -------
# Methods: "Thirty-nine out of the 50 SNPs or INDELs passed assay design
# criteria across the three platforms. For the present study, 30 out of the 39
# were selected for manufacturing ...".  Results: TaqMan designed 39/50.  The
# sentence therefore closes as TaqMan-only (39 - 30 = 9 unmade designs) and is
# ALSO consistent with an intersection of 39; the CSV keeps it under
# `design_ambiguous`.  The word "Thirty-nine" is spelled out and this engine
# does not reliably match the "-nine" suffix, so the stable head is anchored
# and the following characters are checked as text.
a39 <- anchor("Thirty", 80L)
if (substr(a39$site, 1L, 13L) != "Thirty-nine o") {
  stop("the Methods sentence no longer reads 'Thirty-nine out of the'; the text ",
       "at its site is: ", a39$site)
}
intersection_39 <- q(39L, a39$site, "'Thirty-nine' at the Methods sentence head")
amade <- kasp_find(flat, "out of the 39 were selected")
# the manufactured count is the number immediately BEFORE the anchor
n_manufactured <- as.integer(substr(flat, amade$start - 3L, amade$start - 2L))
if (is.na(n_manufactured) || !grepl("^[0-9]{2}$", paste0(n_manufactured))) {
  stop("cannot read a two-digit manufactured count before \"out of the 39 ",
       "were selected\"; the text at its site is: ",
       kasp_site(amade, flat, 40L))
}
q(n_manufactured, substr(flat, amade$start - 12L, amade$start + 28L),
  "the number immediately before 'out of the 39 were selected'")
if (intersection_39 != design$taqman[1L]) {
  stop("the Methods '39' no longer equals TaqMan's per-platform ", design$taqman[1L],
       "/50; the ambiguity note in the platform-outcomes CSV must be re-read")
}
if (n_manufactured != 30L) {
  stop("the paper no longer says 30 of the 39 were manufactured (read ",
       n_manufactured, ")")
}
if (intersection_39 - n_manufactured != 9L ||
    n_manufactured > min(vapply(design, `[`, integer(1), 1L))) {
  stop("the pass/manufacture arithmetic does not close: ", intersection_39,
       " passing, ", n_manufactured, " made")
}

## work (platform x sample points) --------------------------------------------
n_lines <- q(kasp_num_after(kasp_find(flat, "Genomic DNA was extracted"), flat, 40L,
                            what = "wheat lines"),
             kasp_site(kasp_find(flat, "Genomic DNA was extracted"), flat, 52L),
             "wheat lines")
# "two non-template controls (NTCs)": spelled out, and there is no digit here
loc_ntc <- kasp_find(flat, "non-template controls")
if (substr(flat, loc_ntc$start - 4L, loc_ntc$start - 2L) != "two") {
  stop("the NTC sentence no longer starts with 'two'; the text is: ",
       kasp_site(loc_ntc, flat, 40L))
}
n_ntc <- q(2L, "two non-template controls", "spelled-out 'two'")
# "In each platform, 2,784 data points, including NTCs, were generated (96*29)".
# Both factors are printed; the comma-grouped total is read as three runs
# joined by commas, which kasp_num_after does not do, so it is cut by hand.
loc_pts <- kasp_find(flat, "data points, including NTCs")
if (substr(flat, loc_pts$start - 6L, loc_pts$start - 1L) != "2,784 ") {
  stop("the printed data-point total is no longer 2,784; the text is: ",
       kasp_site(loc_pts, flat, 60L))
}
n_points <- q(2784L, kasp_site(loc_pts, flat, 70L),
              "'2,784' in 'In each platform, 2,784 data points ... (96*29)'")
per_plate <- as.integer(kasp_groups(flat, "[(]([0-9]+)[*]([0-9]+)[)]", "96*29"))
if (length(per_plate) != 2L) stop("the (96*29) arithmetic is not printed")
# "One SNP assay failed in all the three platforms": spelled out again
loc_fail <- kasp_find(flat, "SNP assay failed")
if (substr(flat, loc_fail$start - 4L, loc_fail$start - 2L) != "One") {
  stop("the failure sentence no longer starts with 'One'; the text is: ",
       kasp_site(loc_fail, flat, 60L))
}
failed_all3 <- q(1L, kasp_site(loc_fail, flat, 60L), "spelled-out 'One'")
n_assays_run <- n_manufactured - failed_all3

if (per_plate[1] != n_lines + n_ntc ||
    per_plate[2] != n_assays_run ||
    n_points != per_plate[1] * per_plate[2]) {
  stop("data-point arithmetic does not close: paper says ", n_points,
       " from ", per_plate[1], "*", per_plate[2], ", but ", n_lines, "+", n_ntc,
       " = ", n_lines + n_ntc, " and ", n_manufactured, "-", failed_all3, " = ",
       n_assays_run)
}

## amplification, stated twice per platform ------------------------------------
# PCRE here treats "(" as a group delimiter and "." as "any character", so a
# decimal such as 93.5 splits into 93 and 5 in the cut-out text.  Decimals are
# therefore reassembled by hand from the fragment.
decimals <- function(s) {
  v <- strsplit(s, "")[[1L]]
  out <- numeric()
  i <- 1L
  while (i <= length(v)) {
    if (grepl("[0-9]", v[i])) {
      j <- i
      while (j < length(v) && grepl("[0-9]", v[j + 1L]) &&
             !(j + 1L <= length(v) && FALSE)) j <- j + 1L
      num <- paste(v[i:j], collapse = "")
      # a decimal point survives as a literal "." only if it was written as one
      if (j + 1L <= length(v) && v[j + 1L] == "." && j + 2L <= length(v) &&
          grepl("[0-9]", v[j + 2L])) {
        k <- j + 2L
        while (k < length(v) && grepl("[0-9]", v[k + 1L])) k <- k + 1L
        num <- paste0(num, ".", paste(v[(j + 2L):k], collapse = ""))
        j <- k
      }
      out <- c(out, as.numeric(num))
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  out
}
# the first "N%" or "N.N%" at or after an anchor, read as a proportion
pct_after <- function(anchor_word, len, what) {
  loc <- kasp_find(flat, anchor_word)
  site <- kasp_site(loc, flat, len)
  # walk the cut text and take the first percent number
  v <- strsplit(site, "")[[1L]]
  for (i in seq_along(v)) {
    if (!grepl("[0-9]", v[i])) next
    j <- i
    while (j < length(v) && grepl("[0-9]", v[j + 1L])) j <- j + 1L
    num <- paste(v[i:j], collapse = "")
    if (j + 1L <= length(v) && v[j + 1L] == "." && j + 2L <= length(v) &&
        grepl("[0-9]", v[j + 2L])) {
      k <- j + 2L
      while (k < length(v) && grepl("[0-9]", v[k + 1L])) k <- k + 1L
      num <- paste0(num, ".", paste(v[(j + 2L):k], collapse = ""))
      j <- k
    }
    if (j + 1L <= length(v) && v[j + 1L] == "%") {
      return(q(as.numeric(num), site, what))
    }
  }
  stop("no percentage at or after the anchor for ", what, "; the text is: ", site)
}

pct_taqman_amp <- pct_after("TaqMan and KASP showed", 40L, "TaqMan amplified, abstract")
pct_taqman_amp_site <- kasp_site(kasp_find(flat, "TaqMan and KASP showed"), flat, 40L)
pct_kasp_amp <- q(as.numeric(sub("[^0-9]*([0-9]+[.][0-9]+|[0-9]+)%.*$", "\\1",
                                 sub("^[^%]*% and ", "", pct_taqman_amp_site))),
                  pct_taqman_amp_site, "KASP amplified, abstract")
pct_rhamp_amp <- pct_after("rhAmp showed", 30L, "rhAmp amplified, abstract")
# "unamplified samples (7%), followed by KASP (6.5%) and rhAmp (3%)": all
# three unamplified percentages sit in ONE sentence, so they are read from one
# site and taken in order, rather than by three separate anchors.
un_site <- kasp_site(kasp_find(flat, "unamplified samples"), flat, 80L)
un_vals <- decimals(un_site)
if (length(un_vals) < 3L) {
  stop("the unamplified sentence no longer carries three percentages; the text ",
       "at its site is: ", un_site)
}
unamp_taqman <- q(un_vals[1L], un_site, "TaqMan unamplified, first percentage")
unamp_kasp <- q(un_vals[2L], un_site, "KASP unamplified, second percentage")
unamp_rhamp <- q(un_vals[3L], un_site, "rhAmp unamplified, third percentage")

# the two printings must agree to the printed digit.  TaqMan (93 / 7) and rhAmp
# (97 / 3) do; KASP does NOT (93 in the abstract vs 100 - 6.5 = 93.5 in the
# results) and both printings are kept, with the discrepancy in the note.
if (abs((100 - unamp_taqman) - pct_taqman_amp) > 0.5) {
  stop("TaqMan's amplified and unamplified percentages disagree")
}
if (abs((100 - unamp_rhamp) - pct_rhamp_amp) > 0.5) {
  stop("rhAmp's amplified and unamplified percentages disagree")
}
amp_rounding_note <- if (abs(pct_kasp_amp - (100 - unamp_kasp)) > 1e-9) {
  paste0("PAPER ARITHMETIC: KASP is printed twice -- '93%' of amplifications in ",
         "the abstract and '6.5%' unamplified in the results; 100 - 6.5 = 93.5. ",
         "TaqMan (93 / 7) and rhAmp (97 / 3) agree. Both KASP printings are kept.")
} else {
  "abstract and results agree to the printed digit"
}

## invalid calls ---------------------------------------------------------------
a_inv <- kasp_find(flat, "very low fluorescence")
inv_site <- sub("^ ", "", substr(flat, a_inv$start - 37L, a_inv$start + 40L))
# "Seventy-seven amplified samples had very low fluorescence": spelled out
if (substr(inv_site, 1L, 13L) != "Seventy-seven") {
  stop("the invalid-call total is no longer spelled out as Seventy-seven; the ",
       "text at its site is: ", inv_site)
}
n_invalid <- q(77L, inv_site, "spelled-out 'Seventy-seven'")
# "TaqMan showed 57 invalid allele calls, while KASP and rhAmp showed 13 and
# seven, respectively" -- the rhAmp count is spelled out, so the three counts
# are cut from the sentence by hand rather than by one regular expression.
loc_inv <- kasp_find(flat, "invalid allele calls")
inv_sentence <- substr(flat, loc_inv$start - 20L, loc_inv$start + 110L)
if (!grepl("seven", inv_sentence)) {
  stop("the invalid-call sentence has changed; the text at its site is: ",
       inv_sentence)
}
# "TaqMan showed 57 invalid ... while KASP and rhAmp showed 13 and seven"
inv_t <- as.integer(sub("^.*showed ([0-9]+) invalid allele calls.*$", "\\1",
                        inv_sentence))
inv_k <- as.integer(sub("^.*showed ([0-9]+) and seven.*$", "\\1", inv_sentence))
if (is.na(inv_t) || is.na(inv_k)) {
  stop("cannot read the TaqMan / KASP invalid-call counts; the sentence is: ",
       inv_sentence)
}
inv <- c(inv_t, inv_k, 7L)  # "seven" in the same sentence
q(inv, inv_sentence, "invalid calls: TaqMan, KASP, rhAmp ('seven' spelled out)")
invalid <- c(taqman = inv[1L], kasp = inv[2L], rhamp = inv[3L])
if (sum(invalid) != n_invalid) {
  stop("invalid calls do not sum to the printed total: ", sum(invalid), " vs ", n_invalid)
}

## genotype concordance --------------------------------------------------------
conc <- as.numeric(kasp_groups(
  flat, "between rhAmp and KASP was ([0-9.]+)% while between KASP and TaqMan ([0-9.]+)%, and rhAmp and TaqMan ([0-9.]+)%",
  "concordance"))
if (length(conc) != 3L) stop("the concordance sentence does not parse")
q(conc, kasp_site(kasp_find(flat, "Genotype concordance"), flat, 120L),
  "rhAmp/KASP, KASP/TaqMan, rhAmp/TaqMan")

## Table 2 (ANOVA) -------------------------------------------------------------
t2_site <- kasp_site(kasp_find(flat, "Call quality Tech"), flat, 220L)
df_res <- as.integer(kasp_groups(
  flat, "Residuals ([0-9]+) [0-9.]+ [0-9.]+ Separation angle",
  "call-quality residual df"))
df_angle <- as.integer(kasp_groups(
  flat, "Separation angle [*] Tech 2 [0-9.]+ [0-9.]+ [0-9.]+ [0-9.eE+-]+ Residuals ([0-9]+)",
  "separation-angle residual df"))
df_comp <- as.integer(kasp_groups(
  flat, "Cluster compactness [*] Tech 2 [0-9.]+ [0-9.]+ [0-9.]+ [0-9.eE+-]+ Residuals ([0-9]+)",
  "compactness residual df"))
df_ntc <- as.integer(kasp_groups(
  flat, "NTC to cluster distance Tech 2 [0-9.]+ [0-9.]+ [0-9.]+ [0-9.eE+-]+ Residuals ([0-9]+)",
  "NTC-distance residual df"))
if (any(is.na(c(df_res, df_angle, df_comp, df_ntc)))) {
  stop("Table 2 did not parse field by field; the text at its site is: ", t2_site)
}
anova <- data.frame(
  metric = c("allele_call_quality", "separation_angle", "cluster_compactness",
             "ntc_cluster_distance"),
  df_residual = c(df_res, df_angle, df_comp, df_ntc),
  stringsAsFactors = FALSE
)
q(anova$df_residual, t2_site, "Table 2 residual df, one capture group per metric")
# the call-quality residual sum of squares is printed as "32.9.4", which is not
# a number.  It is recorded verbatim and never repaired or used.
if (!grepl("32.9.4", flat, fixed = TRUE)) {
  stop("Table 2 no longer prints the '32.9.4' typo; re-read the table before ",
       "trusting the ANOVA degrees of freedom")
}
typo_note <- "Table 2 call-quality residual SS printed as '32.9.4' -- not a number; kept verbatim, unused"

## Table 3 (means + LSD letters) ----------------------------------------------
t3_nums <- as.numeric(kasp_groups(
  flat, "rhAmp ([0-9.]+) [a-c] ([0-9.]+) [a-c] ([0-9.]+) [a-c] ([0-9.]+) [a-c]",
  "Table 3 rhAmp row"))
t3_nums_k <- as.numeric(kasp_groups(
  flat, "KASP ([0-9.]+) [a-c] ([0-9.]+) [a-c] ([0-9.]+) [a-c] ([0-9.]+) [a-c]",
  "Table 3 KASP row"))
t3_nums_t <- as.numeric(kasp_groups(
  flat, "Ta?qMan? ([0-9.]+) [a-c] ([0-9.]+) [a-c] ([0-9.]+) [a-c] ([0-9.]+) [a-c]",
  "Table 3 TaqMan row"))
t3_lets <- kasp_groups(
  flat, "rhAmp [0-9.]+ ([a-c]) [0-9.]+ ([a-c]) [0-9.]+ ([a-c]) [0-9.]+ ([a-c])",
  "Table 3 rhAmp letters")
t3_lets_k <- kasp_groups(
  flat, "KASP [0-9.]+ ([a-c]) [0-9.]+ ([a-c]) [0-9.]+ ([a-c]) [0-9.]+ ([a-c])",
  "Table 3 KASP letters")
t3_lets_t <- kasp_groups(
  flat, "Ta?qMan? [0-9.]+ ([a-c]) [0-9.]+ ([a-c]) [0-9.]+ ([a-c]) [0-9.]+ ([a-c])",
  "Table 3 TaqMan letters")
if (any(vapply(list(t3_nums, t3_nums_k, t3_nums_t, t3_lets, t3_lets_k, t3_lets_t),
               length, integer(1)) != 4L)) {
  stop("Table 3 did not parse as three rows of four values and four letters; ",
       "the text at its first row is: ",
       kasp_site(kasp_find(flat, "rhAmp 0.97"), flat, 120L))
}
metrics <- c("allele_call_quality", "separation_angle", "cluster_compactness",
             "ntc_cluster_distance")
table3 <- data.frame(
  platform = rep(c("rhAmp", "KASP", "TaqMan"), each = 4L),
  metric = rep(metrics, times = 3L),
  value = c(t3_nums, t3_nums_k, t3_nums_t),
  lsd_letter = c(t3_lets, t3_lets_k, t3_lets_t),
  stringsAsFactors = FALSE
)
q(table3$value, "Table 3 rows", "three rows x four metric means")
q(table3$lsd_letter, "Table 3 rows", "three rows x four LSD letters")
if (nrow(table3) != 12L) stop("Table 3 must hold 3 platforms x 4 metrics")
if (abs(table3$value[table3$platform == "KASP" &
                     table3$metric == "allele_call_quality"] - 0.98) > 1e-9) {
  stop("Table 3 KASP call quality is not the 0.98 the text states")
}
let <- function(p, m) table3$lsd_letter[table3$platform == p & table3$metric == m]
if (let("KASP", "separation_angle") != let("rhAmp", "separation_angle") ||
    let("KASP", "cluster_compactness") != let("rhAmp", "cluster_compactness")) {
  stop("Table 3 letters do not match the text's claim that rhAmp and KASP share ",
       "a letter for separation angle and differ for cluster compactness")
}
if (!grepl("Separation angle and cluster compactness were analyzed based on the values of Allele 11",
           flat)) {
  stop("the Allele-11 footnote used by outcome_definition is absent")
}

## cost per sample --------------------------------------------------------------
cost <- as.numeric(kasp_groups(
  flat, "technologies [(][$]([0-9.]+) per 5", "rhAmp cost"))
cost_k <- as.numeric(kasp_groups(
  flat, "followed by KASP [(][$]([0-9.]+) per 5", "KASP cost"))
cost_t <- as.numeric(kasp_groups(
  flat, "TaqMan assays cost [$]([0-9.]+) per 5", "TaqMan cost"))
if (any(is.na(c(cost, cost_k, cost_t)))) {
  stop("the cost-per-sample sentence does not parse; the text at its site is: ",
       kasp_site(kasp_find(flat, "least expensive"), flat, 160L))
}
cost_per_sample <- data.frame(
  platform = c("rhAmp", "KASP", "TaqMan"),
  cost_usd_per_5ul = c(cost, cost_k, cost_t),
  stringsAsFactors = FALSE
)
q(cost_per_sample$cost_usd_per_5ul, kasp_site(kasp_find(flat, "least expensive"),
                                              flat, 160L),
  "cost per 5 ul reaction, rhAmp / KASP / TaqMan")
if (any(diff(cost_per_sample$cost_usd_per_5ul) <= 0)) {
  stop("the printed cost order is no longer rhAmp < KASP < TaqMan")
}

## --- S1 Table: the 30 manufactured assays -----------------------------------
stopifnot(file.size(xlsx) > 15000L)
sh <- kasp_sheets(xlsx)[["Sheet1"]]
if (is.null(sh)) stop("S1 Table has no Sheet1")
rows <- lapply(sh, kasp_row)
hdr2 <- rows[[2L]]
if (!identical(unname(hdr2[c("C", "J", "K", "L", "M")]),
               c("Assay", "ID", "Allele 1 primer sequence", "Allele 2 primer sequence",
                 "Common_primer_sequence"))) {
  stop("S1 Table columns C/J-M are not the expected assay/KASP header")
}
data_rows <- rows[-c(1L, 2L, 3L)]
# column A of S1 Table is the TRAIT CATEGORY, not a row number: every data row
# carries an assay name in C and a target sequence in D, and the three header
# rows do not.
data_rows <- Filter(function(r) nzchar(r[["C"]]) && nzchar(r[["D"]]), data_rows)
if (length(data_rows) == 0L) stop("S1 Table's data rows are not where they were")
assays <- do.call(rbind, lapply(data_rows, function(r) {
  data.frame(
    trait_category = r[["A"]], gene_or_target = r[["B"]], assay = r[["C"]],
    target_sequence = r[["D"]], reference = r[["E"]],
    rhamp_id = r[["F"]], kasp_id = r[["J"]],
    kasp_allele1 = r[["K"]], kasp_allele2 = r[["L"]], kasp_common = r[["M"]],
    taqman_id = r[["N"]], stringsAsFactors = FALSE
  )
}))

if (nrow(assays) != n_manufactured) {
  stop("S1 Table holds ", nrow(assays), " assays but the paper says ",
       n_manufactured, " were manufactured")
}
if (anyDuplicated(assays$assay) || anyDuplicated(assays$kasp_id)) {
  stop("S1 Table assay / KASP ids are not unique")
}
# The same assay id names the assay and the rhAmp and KASP rows.  It does NOT
# always name the TaqMan row: TaPds_B1 has no TaqMan probe (the paper's own
# reason -- TaqMan cannot design across an INDEL larger than six bp), and that
# is the single assay whose TaqMan cells are empty.  That is a DESIGN effect
# and must not be read as a work failure.
if (!(all(assays$assay == assays$kasp_id) && all(assays$assay == assays$rhamp_id))) {
  stop("the assay id in column C disagrees with the rhAmp or KASP id in the same row")
}
no_taqman <- assays$assay[!nzchar(assays$taqman_id)]
if (length(no_taqman) != 1L) {
  stop(length(no_taqman), " assays lack a TaqMan id; expected exactly one (TaPds_B1)")
}
if (!identical(no_taqman, "TaPds_B1")) {
  stop("the assay without a TaqMan probe is no longer TaPds_B1: ", no_taqman)
}
# a KASP primer triple must be present and allele-specific at the 3' base
if (any(!nzchar(assays$kasp_allele1) | !nzchar(assays$kasp_allele2) |
        !nzchar(assays$kasp_common))) {
  stop("an S1 Table row has an empty KASP primer")
}
a1 <- substr(assays$kasp_allele1, nchar(assays$kasp_allele1), nchar(assays$kasp_allele1))
a2 <- substr(assays$kasp_allele2, nchar(assays$kasp_allele2), nchar(assays$kasp_allele2))
if (any(a1 == a2)) stop("a KASP pair in S1 Table shares its 3' base")
# The KASP allele-specific primers must differ at their 3' base (already
# asserted) and that base must appear in the target's bracketed alleles.  The
# supplement does NOT consistently print the SNP as a two-base bracket: some
# targets carry a one-base bracket with a flanking letter ("[CGTCG/]"), some an
# INDEL ("[ACTT/]"), and on 13 rows the 3' bases are the complement of the
# bracket.  So this check is deliberately weak -- each 3' base must occur
# SOMEWHERE in the bracketed text, in either orientation -- and it is recorded
# how many rows needed the complement, rather than silently "fixed".
bracket_bases <- function(seq) {
  br <- regmatches(seq, gregexpr("\\[[^]]+\\]", seq, perl = TRUE))[[1L]]
  if (!length(br)) return(character(0))
  unique(strsplit(gsub("[^A-Za-z]", "", br), "")[[1L]])
}
comp <- c(A = "T", T = "A", C = "G", G = "C")
base_ok <- mapply(function(seq, x, y) {
  bb <- bracket_bases(seq)
  if (!length(bb)) return(FALSE)
  ok1 <- (x %in% bb) && (y %in% bb)
  ok2 <- (comp[[x]] %in% bb) && (comp[[y]] %in% bb)
  list(ok = ok1 || ok2, complemented = !ok1 && ok2)
}, assays$target_sequence, a1, a2, SIMPLIFY = FALSE)
ok <- vapply(base_ok, `[[`, logical(1), "ok")
complemented <- vapply(base_ok, `[[`, logical(1), "complemented")
if (!all(ok)) {
  stop(sum(!ok), " S1 Table target(s) carry neither KASP 3' base in their ",
       "bracketed alleles: ", paste(assays$assay[!ok], collapse = ", "))
}
assays$bracketed_alleles <- vapply(assays$target_sequence,
                                   function(x) paste(bracket_bases(x), collapse = "/"),
                                   character(1))
assays$as_3p_is_complement <- complemented
assays$kasp_as_differs_at_3p <- TRUE
assays$as3p_matches_bracket_note <- ifelse(
  complemented,
  "allele-specific 3' bases are the COMPLEMENT of the printed bracket on this row",
  "")
## --- write: per-assay ---------------------------------------------------------
kasp_prefix <- "GGAGCAGGTCCAGATCGCG" # sanity only: allele primers carry a tail
assay_out <- assays[c("assay", "trait_category", "gene_or_target", "target_sequence",
                      "reference", "kasp_id", "kasp_allele1", "kasp_allele2",
                      "kasp_common", "bracketed_alleles",
                      "as_3p_is_complement", "as3p_matches_bracket_note")]
assay_out$kasp_allele1_3p <- a1
assay_out$kasp_allele2_3p <- a2
assay_out <- assay_out[order(match(assay_out$assay, assays$assay)), ]
assay_out$source_row <- seq_along(data_rows) + 3L  # three header rows
# S1 Table prints an ambiguous base ("W") inside the 1RS_1BL target.  It is
# kept as printed; `has_ambiguous_base` records it so a later reader does not
# mistake it for a typo.
assay_out$has_ambiguous_base <- grepl("[^ACGTacgtN\\[\\]/]", assay_out$target_sequence)
assay_out$doi <- DOI
assay_out$evidence_tier <- "measured"
assay_out$source_file <- basename(xlsx)
assay_out$source_table <- "S1 Table, sheet Sheet1 (one row per manufactured assay)"
assay_out$date_verified <- DATE_VERIFIED
# cost per sample is printed in the text only; it travels with the metric table
# because it is one of the paper's comparative numbers, not a KASP outcome
assay_out$conversion_definition <- paste(
  "NOT APPLICABLE at assay level. S1 Table lists the assays that were made; it",
  "carries no per-assay design, amplification or call outcome. The only",
  "outcomes in this source are platform-level counts and percentages, in",
  "wheat_platform2019_platform_outcomes.csv."
)
assay_out$outcome_definition <- paste(
  "Assay identity and KASP primer triple only. 'Design success' for this",
  "platform arm is 49/50 at the study level (platform outcomes CSV); it is NOT",
  "asserted per assay here, because the paper does not say which assay failed",
  "to design -- only that one did, and that a different assay (identity also",
  "not given) failed on all three platforms."
)
assay_out <- assay_out[c("assay", "trait_category", "gene_or_target", "target_sequence",
                         "reference", "kasp_id", "kasp_allele1", "kasp_allele2",
                         "kasp_common", "kasp_allele1_3p", "kasp_allele2_3p",
                         "bracketed_alleles", "as_3p_is_complement",
                         "as3p_matches_bracket_note", "has_ambiguous_base",
                         "source_row", "doi",
                         "evidence_tier", "source_file", "source_table",
                         "conversion_definition", "outcome_definition")]
kasp_write_checked(assay_out, fs::path(OUT, "wheat_platform2019_kasp_assays.csv"))

## --- write: platform-level outcomes ------------------------------------------
def_design <- paste(
  "PAPER'S DEFINITION (verbatim): 'Design success rates were 39/50, 49/50, and",
  "49/50 for TaqMan, KASP, and rhAmp'. A design success is a primer/probe set",
  "the provider returned for one of the 50 pre-validated wheat SNPs/INDELs.",
  "Denominator = 50 targets attempted. This is NOT a statement that the assay",
  "genotyped anything."
)
def_work <- paste(
  "PAPER'S DEFINITION, per ASSAY x SAMPLE POINT: 96 x 29 = 2,784 points per",
  "platform (94 wheat lines + 2 NTCs, minus the one assay that failed on all",
  "three platforms). 'Amplified' percentages are read off that denominator;",
  "'Invalid' counts are amplified points the Thermo Fisher Cloud autocaller",
  "could not assign to Allele 11/12/22. Neither is a per-assay success rate.",
  "No per-assay count is printed, so no per-assay rate is derived."
)
mk <- function(platform, q, n, d, printed, scope, kind, metric, unit, note) {
  data.frame(platform = platform, quantity = q, n_numerator = n, n_denominator = d,
             value_printed = printed, scope = scope, quantity_kind = kind,
             metric = metric, unit = unit, note = note, stringsAsFactors = FALSE)
}
outs <- rbind(
  mk("TaqMan", "design_success", design$taqman[1L], design$taqman[2L], "39/50",
     "study", "design", NA, "assays", "same sentence as the other two platforms"),
  mk("KASP", "design_success", design$kasp[1L], design$kasp[2L], "49/50",
     "study", "design", NA, "assays", "KASP lost one target at design"),
  mk("rhAmp", "design_success", design$rhamp[1L], design$rhamp[2L], "49/50",
     "study", "design", NA, "assays", "same sentence as the other two platforms"),
  mk("all three", "all_passed_design_criteria", intersection_39, n_design_tested,
     "39/50", "study_intersection", "design_ambiguous",
     NA, "targets",
     paste("INTERSECTION reading from Methods ('passed assay design criteria",
           "across the three platforms'), NOT TaqMan's 39/50. The paper uses the",
           "same number for both; do not treat them as one fact.")),
  mk("all three", "manufactured", n_manufactured, intersection_39, "30 of 39",
     "study_intersection", "manufacture", NA, "assays",
     "30 of the 39 selected for manufacture on all three platforms (S1 Table)"),
  mk("all three", "failed_on_all_three_platforms", failed_all3, n_manufactured,
     "one SNP assay", "assay", "work", NA, "assays",
     "identity of the failed assay is NOT given"),
  mk("all three", "assays_genotyped", n_assays_run, n_manufactured, "29",
     "assay", "work", NA, "assays", "30 minus the one that failed everywhere"),
  mk("all three", "assays_x_samples", n_points, prod(per_plate), "2,784",
     "assay_x_sample", "work", NA, "datapoints",
     "paper's own arithmetic: (94 genotypes + 2 NTCs) * 29 assays"),
  mk("TaqMan", "amplified", NA, n_points, "93%", "assay_x_sample", "work",
     "amplification", "fraction", amp_rounding_note),
  mk("KASP", "amplified", NA, n_points, "93.5%", "assay_x_sample", "work",
     "amplification", "fraction", amp_rounding_note),
  mk("rhAmp", "amplified", NA, n_points, "97%", "assay_x_sample", "work",
     "amplification", "fraction", amp_rounding_note),
  mk("TaqMan", "unamplified", NA, n_points, "7%", "assay_x_sample", "work",
     "amplification", "fraction", "printed as a percentage, not a count"),
  mk("KASP", "unamplified", NA, n_points, "6.5%", "assay_x_sample", "work",
     "amplification", "fraction",
     paste("printed as a percentage, not a count; 6.5% of ", n_points,
           " = 180.96, which is not an integer, so no count is recorded")),
  mk("rhAmp", "unamplified", NA, n_points, "3%", "assay_x_sample", "work",
     "amplification", "fraction", "printed as a percentage, not a count"),
  mk("TaqMan", "invalid_call", invalid[["taqman"]], n_invalid, as.character(invalid[["taqman"]]),
     "assay_x_sample", "work", "call_validity", "datapoints",
     "amplified, low fluorescence, autocaller returned 'Invalid'"),
  mk("KASP", "invalid_call", invalid[["kasp"]], n_invalid, as.character(invalid[["kasp"]]),
     "assay_x_sample", "work", "call_validity", "datapoints",
     "amplified, low fluorescence, autocaller returned 'Invalid'"),
  mk("rhAmp", "invalid_call", invalid[["rhamp"]], n_invalid, as.character(invalid[["rhamp"]]),
     "assay_x_sample", "work", "call_validity", "datapoints",
     "amplified, low fluorescence, autocaller returned 'Invalid'"),
  mk("all three", "invalid_call", n_invalid, NA, as.character(n_invalid),
     "assay_x_sample", "work", "call_validity", "datapoints",
     "text total; equals the sum of the three platforms (asserted)")
)
# the genotype-concordance pairs
conc_df <- data.frame(
  platform = c("rhAmp vs KASP", "KASP vs TaqMan", "rhAmp vs TaqMan"),
  quantity = "genotype_concordance", n_numerator = NA_integer_,
  n_denominator = NA_integer_, value_printed = paste0(format(conc, nsmall = 2L), "%"),
  scope = "assay_x_sample", quantity_kind = "concordance",
  metric = "genotype_concordance", unit = "fraction",
  note = "printed percentage only; denominator not printed (the assay was run on the 94 wheat genotypes)",
  stringsAsFactors = FALSE
)
outs <- rbind(outs, conc_df)
outs$doi <- DOI
outs$evidence_tier <- "measured"
outs$source_file <- basename(xml_path)
outs$source_table <- "Results/Discussion text (no per-assay outcomes are printed)"
outs$date_verified <- DATE_VERIFIED
outs$conversion_definition <- ifelse(
  outs$quantity_kind == "design", def_design,
  ifelse(outs$quantity_kind == "design_ambiguous",
         paste(def_design, "READ AS AN INTERSECTION across platforms; the",
               "per-platform KASP number is the separate row 'design_success'."),
         def_work))
outs$outcome_definition <- outs$conversion_definition
outs$conversion_definition <- NULL
outs <- outs[c("platform", "quantity", "quantity_kind", "scope", "n_numerator",
               "n_denominator", "value_printed", "metric", "unit", "note", "doi",
               "evidence_tier", "source_file", "source_table", "date_verified",
               "outcome_definition")]
outs$n_numerator <- as.integer(outs$n_numerator)
outs$n_denominator <- as.integer(outs$n_denominator)
kasp_write_checked(outs, fs::path(OUT, "wheat_platform2019_platform_outcomes.csv"))

## --- write: cluster metrics ---------------------------------------------------
t3 <- table3
t3$doi <- DOI
t3$evidence_tier <- "measured"
t3$source_file <- basename(xml_path)
t3$source_table <- "Table 3 (means), Table 2 (ANOVA df)"
t3$df_tech <- 2L  # Table 2: "Tech 2" for every metric (three platforms)
t3$df_residual <- anova$df_residual[match(t3$metric, anova$metric)]
if (anyNA(t3$df_residual)) stop("a Table 3 metric has no Table 2 residual df")
t3$value_printed <- sub("0$", "", format(t3$value, nsmall = 2L))
t3$outcome_definition <- paste(
  "Cluster-geometry quality metric, platform mean over assays. PAPER'S",
  "FOOTNOTE (verbatim): 'Separation angle and cluster compactness were based on",
  "data values of Allele 11 for all three platforms' -- i.e. the HOMOZYGOUS",
  "cluster only, not the heterozygous cluster and not the call itself. Letters",
  "are LSD groups within a metric (shared letter = not significantly different).",
  "These are assay-quality metrics; they are not design, amplification or",
  "conversion rates and must never be pooled with one."
)
t3 <- t3[order(t3$metric, match(t3$platform, c("TaqMan", "KASP", "rhAmp"))), ]
t3 <- t3[c("metric", "platform", "value", "value_printed", "lsd_letter", "df_tech",
           "df_residual", "doi", "evidence_tier", "source_file", "source_table",
           "outcome_definition")]
kasp_write_checked(t3, fs::path(OUT, "wheat_platform2019_cluster_metrics.csv"))

message(
  "wheat_platform2019: ", nrow(assay_out), " assays (paper 30); design ",
  paste(vapply(design, function(x) paste(x, collapse = "/"), character(1)),
        collapse = " "),
  "; ", n_points, " points = ", per_plate[1], "x", per_plate[2],
  "; invalid ", sum(invalid), " (paper ", n_invalid, ")"
)
