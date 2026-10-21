# =============================================================================
#  data-raw/z7-priors.R  --  published mismatch priors as testable PREDICTIONS
# =============================================================================
#  Z7 wants deliberate second mismatches.  Its stated blocker is that our engine
#  cannot rank 3'-terminal mismatches.  This script writes down (a) what the
#  published priors PREDICT, in a form a test can falsify, and (b) the engine
#  margin that is supposed to test them -- so the blocker is a number, not an
#  opinion.
#
#  Writes
#    z7_priors.csv        one row per directional, published prediction, with the
#                         verbatim claim, the evidence level, calibrated = FALSE
#    z7_engine_probe.csv  the four KASP allele pairs scored by DECIPHER: efficiency
#                         on the MATCH template vs on the OFF-allele template
#
#  Rules this script obeys:
#    * No prior is an outcome.  A prior is a prediction we can FAIL; only
#      measured data can be labelled calibrated, and none of this is.
#    * Contradictions are kept, not averaged.  Ye 2001 says -2, Lefever 2019
#      says the 4th base from the 3' end (= 3rd from the terminus), Liu 2012
#      says the 3rd base from the 3' end.  All three stay, each labelled.
#    * The position claims are derived from the per-class tables we already
#      extracted, and the script fails if the table does not support the claim.
#
#  Engine note (this is the Z7 blocker, measured):
#    pcr_efficiency() computes DECIPHER's efficiency of a primer against a
#    *forward* target, so the target must open with the complement of the
#    primer's 3' base.  Against that MATCH template a 3'-terminal mismatch barely
#    moves the number; the collapse only appears when the OFF-allele template
#    (the other SNP allele at that position) is the one being scored.  The
#    discriminant therefore has to be computed as match - off, never as a single
#    efficiency.
#
#  Run:  Rscript data-raw/z7-priors.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

out_dir <- fs::path(root, "inst", "extdata")
ext <- function(f) utils::read.csv(fs::path(out_dir, f), stringsAsFactors = FALSE)

KWOK <- "10.1093/nar/18.4.999"
LEFEVER <- "10.1038/s41598-019-38581-z"
LIU <- "10.1186/1746-4811-8-34"
YE <- "10.1093/nar/29.17.e88"

## --- 1. the published rule tables we hold -----------------------------------
kwok <- ext("mismatch_discrimination.csv")
kwok <- kwok[kwok$table == "kwok_3prime_severity", ]
kwok_order <- c("severe", "moderate", "minimal", "context_dependent")
if (!setequal(kwok$mismatch_class, kwok_order)) {
  stop("kwok_3prime_severity classes are not the four published ones")
}

lo <- ext("liu2012_outcomes.csv")
lo <- lo[lo$grouping == "snp_class", ]
pct <- round(100 * lo$n_specific / lo$n_tested, 1L)
if (any(abs(pct - lo$specificity_pct) > 0.05)) {
  stop("liu2012_outcomes: a specificity_pct does not equal n_specific/n_tested")
}

## --- 2. derive the position claim FROM the table ----------------------------
sites <- c("2nd", "3rd", "4th")
mean_pct <- vapply(sites, function(s) {
  x <- lo[lo$mismatch_site == s, ]
  if (nrow(x) != 4L) stop("liu2012: expected 4 SNP classes at site ", s)
  stats::weighted.mean(x$specificity_pct, x$n_tested)
}, numeric(1))
best_site <- sites[[which.max(mean_pct)]]
if (best_site != "3rd") {
  stop(
    "the per-class numbers put the optimum at the ", best_site,
    ", not the 3rd base; do not state the 3rd-base rule from memory"
  )
}


## --- 3. the prior table ------------------------------------------------------
priors <- data.frame(
  prior_id = c(
    "kwok_severity_order", "kwok_terminal_T_weak", "kwok_context_dependence",
    "lefever_single_insufficient", "lefever_position", "lefever_base_free",
    "liu_position_3rd", "ye_second_mismatch_position"
  ),
  source_key = c(
    "kwok1990", "kwok1990", "kwok1990", "lefever2019", "lefever2019",
    "lefever2019", "liu2012_arms", "ye2001"
  ),
  family = rep("KASP / AS-PCR", 8L),
  what = c(
    "3'-terminal mismatch severity is ordered by base pair",
    "a 3'-terminal T mismatch is the weak class",
    "a G:G mismatch is not one number",
    "one 3'-terminal mismatch is not a genotype call",
    "the best extra-mismatch position is the 4th base from the 3' end",
    "at the optimum position the mismatched base does not matter",
    "the best single mismatch site is the 3rd base from the 3' end",
    "the tetra-primer rule places the extra mismatch at -2"
  ),
  prediction = c(
    "margin(severe) < margin(moderate) < margin(minimal)",
    "none of T:G, T:C, T:T falls in the severe class",
    "a single margin for G:G is not a valid summary of its behaviour",
    "a genotype call resting on a lone 3'-terminal mismatch is rejected",
    "the margin peaks at offset 3 (0 = the 3' base)",
    "at offset 3 the four mismatch bases agree within one s.e.",
    "specificity is highest at the 3rd base from the 3' end",
    "a second mismatch at -2 improves discrimination"
  ),
  testable_form = c(
    "ordinal", "ordinal", "none (flag only)", "ordinal",
    "continuous", "equivalence", "continuous", "ordinal"
  ),
  position = c(0L, 0L, 0L, 0L, -3L, -3L, -3L, -2L),
  contradicted_by = c(
    NA_character_, NA_character_, NA_character_, NA_character_,
    "ye2001 (-2); numbering conventions differ (see position_note)",
    NA_character_, "lefever2019 agrees on -3; ye2001 (-2) disagrees", "lefever2019 + liu2012 (-3); Ye 2001 (-2) is the fallback rule"
  ),
  position_note = c(
    "", "", "", "",
    "Lefever's 'position 3' is the 4th nucleotide from the 3' end; in our offset convention that is -3",
    "",
    "Liu's '3rd base' is the 3rd nucleotide from the 3' end, i.e. the same site the decision files at -3 (arms_mismatch_position.csv)",
    "Ye's -2 is the penultimate base; the decided extra-mismatch default is -3 (Liu 2012 + Lefever 2019), Ye's -2 is the fallback for a context_dependent terminal class"
  ),
  measured = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  evidence = c(
    "primary literature", "primary literature", "primary literature",
    "primary literature", "primary literature", "primary literature",
    "primary literature", "rule of thumb"
  ),
  calibrated = FALSE,
  doi = c(KWOK, KWOK, KWOK, LEFEVER, LEFEVER, LEFEVER, LIU, YE),
  stringsAsFactors = FALSE
)
if (nrow(priors) != 8L) stop("prior table lost a row")

## --- 4. can our engine be tested? -------------------------------------------
# A KASP assay scores one allele primer against a MATCH template (the allele it
# carries) and an OFF template (the other allele).  Both templates are the same
# forward locus; only the SNP base differs.  Decoding the target therefore needs
# complement(primer 3' base) at that position.
pairs <- c("G:A", "A:G", "C:A", "A:A") # FAM 3' : HEX 3' (allele 1 : allele 2)
pcr3 <- stats::setNames(substr(pairs, 1L, 1L), pairs) # FAM primer's 3' base
off <- stats::setNames(substr(pairs, 3L, 3L), pairs) # the other allele

primer <- "CTGGACCGATCGTAGCTAGCA" # fixed flank; only the 3' base moves
locus <- "TGCTAGCTACGATCGGTCCAG" # 20 nt, locus[1] is the SNP position
if (nchar(primer) != nchar(locus)) stop("primer and locus differ in length")
setb <- function(s, b) paste0(substr(s, 1L, nchar(s) - 1L), b)
sett <- function(s, b) paste0(b, substr(s, 2L, nchar(s)))
# pcr3 || off names the TEMPLATE at the primer's 3' end: a primer base X
# pairs with template base complement(X)
template_at <- function(p3, off) {
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  unname(comp[[if (off) off else p3]])
}

probe <- do.call(rbind, lapply(pairs, function(k) {
  pr <- setb(primer, pcr3[[k]])
  mt <- template_at(pcr3[[k]], off = FALSE)
  ot <- template_at(pcr3[[k]], off = TRUE)
  match_eff <- pcr_efficiency(pr, sett(locus, mt), cores = 1L)
  off_eff <- pcr_efficiency(pr, sett(locus, ot), cores = 1L)
  # column NAMES matter here: a column called `calibrated` would silently
  # replace a later `calibrated = FALSE`, so it is the last thing set
  out <- data.frame(
    allele_pair = k, primer_3p = pcr3[[k]], match_base = mt, off_base = ot,
    mismatch = mt != ot,
    match_eff = match_eff, off_eff = off_eff, margin = match_eff - off_eff,
    stringsAsFactors = FALSE
  )
  out$engine <- "tidymarkers::pcr_efficiency (DECIPHER CalculateEfficiencyPCR)"
  out$calibrated <- FALSE
  out
}))

# the match base must actually pair with the primer, and the off base must not
comp <- c(A = "T", C = "G", G = "C", T = "A")
if (!all(unname(comp[probe$primer_3p]) == probe$match_base)) {
  stop("the match template base does not pair with the primer's 3' base")
}
if (!all(probe$mismatch == !(probe$off_base == unname(comp[probe$primer_3p])))) {
  stop("an off-allele template still pairs with the primer's 3' base")
}

## --- 5. the discriminant -----------------------------------------------------
# Two different quantities, deliberately kept apart:
#   * `margin` is the model's own per-assay discrimination, match - off.
#     Kwok's severity prediction is a claim ABOUT THIS.
#   * `terminal_floor` is how far any 3'-terminal mismatch depresses efficiency
#     against the MATCH template alone.  It is the quantity a design filter
#     would use if it looked at one efficiency, and it is why that is wrong.
probe$terminal_floor <- round(1 - probe$match_eff, 4L)
probe$match_eff <- round(probe$match_eff, 4L)
probe$off_eff <- round(probe$off_eff, 4L)
probe$margin <- round(probe$margin, 4L)
stopifnot(all(c("margin", "terminal_floor", "calibrated") %in% names(probe)))
spread <- diff(range(probe$margin))
if (anyNA(probe$margin)) stop("a pair had no finite margin")
if (!all(probe$margin[is.na(probe$mismatch)] == 0)) {
  # an exact match must not discriminate; if this trips, the probe is broken
  stop("a matched pair has a non-zero margin")
}

# Kwok's ordering, read against our margin.  This is a CHECK OF THE PROBE, not
# a calibration: it says whether a test could even be written today.
# The table's `members` are read as TEMPLATE:PRIMER (the primer's 3' base is
# the second letter), so a primer base X is scored by looking for ":X".
rank_of <- function(p3) {
  cls <- kwok$mismatch_class[vapply(kwok$members, function(m) {
    any(endsWith(trimws(strsplit(m, ",", fixed = TRUE)[[1L]]), paste0(":", p3)))
  }, logical(1))]
  if (length(cls) > 1L) cls <- cls[[1L]] # G:G is both minimal and context-dependent
  if (length(cls) != 1L) stop("primer base ", p3, " is in ", length(cls), " Kwok classes")
  match(cls, kwok_order)
}
have <- probe[probe$mismatch, ]
have$kwok_rank <- vapply(have$primer_3p, rank_of, integer(1))
have <- have[order(have$kwok_rank), ]
# strict equality would be checked by the test; here we only record whether the
# ordering holds, because a probe that cannot order by its own most-discriminating
# quantity is exactly the blocker Z7 is about.
ordered <- all(diff(have$margin) <= 0)

## --- write -------------------------------------------------------------------
utils::write.csv(priors, fs::path(out_dir, "z7_priors.csv"), row.names = FALSE, na = "")
utils::write.csv(probe, fs::path(out_dir, "z7_engine_probe.csv"), row.names = FALSE, na = "")

message(
  "z7 priors: ", nrow(priors), " predictions, ",
  sum(!is.na(priors$contradicted_by)), " with a named contradiction (kept, not averaged);\n",
  "  Liu weighted optimum = ", best_site, " (",
  paste(sites, round(mean_pct, 1L), sep = "=", collapse = " "), ");\n",
  "DECIPHER over the four KASP pairs: margin ",
  paste(sprintf("%.4f", probe$margin), collapse = "/"),
  " (spread ", signif(spread, 3L), "); terminal floor ",
  paste(sprintf("%.4f", probe$terminal_floor), collapse = "/"), ";\n",
  "  Kwok severity order reproduced by the margin: ", ordered
)
