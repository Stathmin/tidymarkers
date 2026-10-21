# =============================================================================
#  data-raw/desilva2009-supplement.R  --  HRM amplicon length vs mutation-scan
#                                          sensitivity/specificity (TP53)
# =============================================================================
#  Garritano S, Gemignani F, Voegele C, Nguyen-Dumont T, Le Calvez-Kelm F,
#  De Silva D, Lesueur F, Landi S, Tavtigian SV. "Determining the effectiveness
#  of High Resolution Melting analysis for SNP genotyping and mutation scanning
#  at the TP53 locus." BMC Genet 2009;10:5.
#  doi:10.1186/1471-2156-10-5  PMID 19173742  PMC2648999  (OA)
#
#  PROVENANCE CORRECTION -- READ THIS BEFORE CITING
#  The backlog row `desilva2009` carried doi:10.1186/1471-2156-10-70, which
#  resolves to a DIFFERENT paper (Hallander & Waldmann, optimisation of
#  selection contribution in monoecious tree breeding). The TP53 HRM paper is
#  10.1186/1471-2156-10-5, and De Silva is a co-author in the middle of the
#  author list, not the first author. The registry error was the kind that
#  produces a citation that "resolves" but supports nothing; it is corrected in
#  scientific_backlog.csv and here.
#
#  WHY THIS FILE EXISTS
#  Z7 wants an external anchor for allele-discrimination design rules. This is
#  the only source we hold that measures the cost of a LONGER amplicon for
#  variant detection: the same 47 samples scanned twice, split by amplicon
#  length at 400 bp, with a true/false positive/negative breakdown per group.
#  It is a mutation-scanning readout, not an allele-specific primer rule, but it
#  is a real measurement of the regime in which discrimination degrades, and it
#  is the outcome data behind the amplicon-composition guidance in Table 2.
#
#  WHAT IS FETCHED
#  * fulltext.xml   -- Europe PMC REST, PMC2648999 (JATS)
#                      Table 1: 22 amplicons (primers, location, size)
#                      Table 2: 5 "genetic situation -> HRM solution" rows
#                      Results paragraphs: the per-group TP/FP/TN/FN counts
#  * supplements/1471-2156-10-5-S1.doc -- Additional file 1, a legacy Word file.
#                      Read with `strings`, NOT a Word library: the text streams
#                      are plain ASCII and appear in exactly the table order
#                      (rs id, forward, reverse, optional length). Rows whose
#                      length is absent are real blanks in the file, and rows
#                      naming several rs ids are genuine multi-SNP amplicons.
#  The BMC and PMC PDF endpoints are bot-walled (an HTML interstitial is served
#  in place of the PDF), so the PDF was not used and is not needed: every
#  number below is in fulltext.xml, which we hold.
#
#  THE ASSERTION (what makes this verified, not transcribed)
#  Table 1's 22 amplicon sizes are re-derived into the paper's own two groups
#  and must reproduce its printed averages:
#     (1) 9 amplicons < 400 nt with mean exactly 286 bp  -- the paper prints
#         "Nine ... with an average length of 286 bp".
#     (2) 13 amplicons >= 400 nt, and the paper's group counts 9 + 13 = 22
#         cover all 22 Table-1 rows.
#     (3) the totals 38+0 and 295+60 sum to 47 (the samples), and 105+23 and
#         339+69 also sum to 47.
#  NOTE (do NOT "fix" this): the >400 bp group's PRINTED mean is 544; our
#  sizes give 545 for the same 13 rows. The size split is unambiguous (9/13 and
#  286 both reproduce exactly), so this is a contradiction between two printed
#  numbers of one paper and is recorded, not resolved. It is 0.2%.
#  If (1)-(3) fail the script stops and writes nothing.
#
#  NOT MODELLED
#  * The per-amplicon TP/FP/TN/FN calls. They exist in the paper ONLY as the
#    two pooled rows (38/0/295/60 and 105/23/339/69); the paper does not print
#    a per-amplicon contingency table, so no per-amplicon outcome is available
#    and none is invented. The 22 primer rows therefore carry primers and size
#    only.
#  * The 45 amplicons of the 67 that are not in Table 1: sizes unknown.
#  * The pooled counts are kept as two rows. They are NOT decomposed in
#    desilva2009_amplicons.csv -- doing so would mean inventing a split.
#  * The 30-SNP x 270-sample genotyping arm reports concordance qualitatively
#    (no counts), so it yields no numeric row.
#  * A 47 x 22 pass/dropout count model: Table 1 does not say which samples
#    had no call, so the 0.83/0.84 specificities cannot be cross-checked.
#
#  CALIBRATION: nothing here is calibrated. The 400 bp split is the authors'
#  choice of grouping, not a threshold we endorse or have fitted.
#
#  Run:  Rscript data-raw/desilva2009-supplement.R [--dir .tmp/lit/desilva2009]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!exists("root")) {
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else {
    getwd()
  }
  setwd(root)
}
suppressPackageStartupMessages({
  library(xml2)
})

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "desilva2009")
ft <- file.path(lit, "fulltext.xml")
s1 <- file.path(lit, "supplements", "1471-2156-10-5-S1.doc")
for (f in c(ft, s1)) if (!file.exists(f)) {
  stop("missing ", f, " (PMC2648999; fetch by hand, see SOURCES-arms.md)")
}

DOI <- "10.1186/1471-2156-10-5"
PMID <- "19173742"
PMCID <- "PMC2648999"
SRC <- "Garritano S, et al. (De Silva D co-author). BMC Genet 2009;10:5"
out_dir <- file.path(root, "inst", "extdata")

doc <- read_xml(ft)
body_p <- xml_text(xml_find_all(doc, "//body//p"))

## --- 1. Table 1: amplicons ---------------------------------------------------
# Column order: Amplicon | Forward | Reverse | Location | size. NOTE: 7 of the
# 22 rows carry a size ONLY in Table 1 -- Additional file 1 leaves their size
# cell blank, so Table 1 is the sole source for those and there is nothing to
# cross-check them against.
tab <- function(ft, label) {
  tw <- xml_find_all(ft, "//table-wrap")
  labs <- vapply(tw, function(t) {
    x <- xml_find_first(t, "./*[local-name()='label']")
    if (inherits(x, "xml_missing")) NA_character_ else
      trimws(gsub("\u00a0", " ", xml_text(x)))
  }, character(1))
  hit <- which(labs == label)
  if (!length(hit)) stop("table not found: ", label)
  lapply(xml_find_all(tw[[hit]], ".//tr"), function(r) {
    cells <- xml_find_all(r, "./*[local-name()='td' or local-name()='th']")
    trimws(gsub("[[:space:]]+", " ", xml_text(cells)))
  })
}

t1 <- tab(doc, "Table 1")
t1 <- Filter(function(r) length(r) == 5L, t1)
# The first row is the header (Amplicon / Forward / Reverse / Location / size).
if (t1[[1]][1] != "Amplicon") stop("Table 1: unexpected first row")
t1 <- t1[-1]
if (any(!nzchar(vapply(t1, `[`, character(1), 1L)))) {
  stop("Table 1: blank amplicon id")
}
amp <- data.frame(
  amplicon = vapply(t1, `[`, character(1), 1L),
  forward = vapply(t1, `[`, character(1), 2L),
  reverse = vapply(t1, `[`, character(1), 3L),
  location = vapply(t1, `[`, character(1), 4L),
  size_bp = suppressWarnings(as.integer(vapply(t1, `[`, character(1), 5L))),
  stringsAsFactors = FALSE)
if (nrow(amp) != 22L || anyNA(amp$size_bp)) {
  stop("Table 1 did not parse into 22 amplicons with sizes")
}
if (any(!grepl("^[ACGT]+$", amp$forward)) || any(!grepl("^[ACGT]+$", amp$reverse))) {
  stop("Table 1: a primer field is not a bare ACGT sequence")
}

## --- 2. Additional file 1: per-SNP primer pairs -------------------------------
# `strings` yields the text streams in table order. The column header is split
# over three lines ("Forward sequence 5" / " > 3" / "Reverse sequence 5"), so
# the body starts at the first rs id. The LAST row (rs17886760) has no length
# cell, so the end cannot be anchored on a length: anchor on the last rs id and
# take the two sequence lines that follow it.
lines <- system2("strings", c("-n", "3", shQuote(s1)), stdout = TRUE)
is_rs_all <- grepl("^rs[0-9]", lines)
first <- which(is_rs_all)[1]
last_row <- tail(which(is_rs_all), 1L)
stopifnot(is.finite(first), is.finite(last_row), last_row > first)
last_len <- last_row + 2L
if (last_len > length(lines) || !all(grepl("^[ACGT]{15,}$", lines[(last_row + 1L):last_len]))) {
  stop("Additional file 1: the last rs row is not followed by two sequences")
}
body <- lines[first:last_len]
end <- length(body)
is_rs <- grepl("^rs[0-9]", body)
is_seq <- grepl("^[ACGT]{15,}$", body)
is_len <- grepl("^[0-9]{2,4}$", body)
if (any(!(is_rs | is_seq | is_len))) {
  stop("Additional file 1: unrecognised line(s): ",
       paste(head(body[!(is_rs | is_seq | is_len)], 3L), collapse = " | "))
}

# Column order in the file is rs | forward | reverse | size, so for each rs
# line take the next two sequence lines, then a length if one immediately
# follows. Rows with no length are a real blank in the file.
rows <- lapply(which(is_rs), function(i) {
  ## forward and reverse are the first two sequence lines after the rs id
  rel <- which(is_seq[i:end])
  if (length(rel) < 2L) {
    stop("Additional file 1: rs line without two sequences: ", body[i])
  }
  j <- i + rel[1:2] - 1L
  k <- j[2] + 1L
  len <- if (k <= end && is_len[k]) as.integer(body[k]) else NA_integer_
  data.frame(snp_ids = body[i], forward = body[j[1]],
             reverse = body[j[2]], size_bp = len, stringsAsFactors = FALSE)
})
snp <- do.call(rbind, rows)
snp$n_snps <- vapply(strsplit(snp$snp_ids, ","), length, 0L)
snp$multi_snp_amplicon <- snp$n_snps > 1L

## --- 3. body: the two groups' printed averages and counts --------------------
num1 <- function(pattern) {
  hit <- grep(pattern, body_p, value = TRUE)
  if (length(hit) != 1L) stop("expected one paragraph matching: ", pattern)
  trimws(gsub("[[:space:]]+", " ", hit))
}
p_small <- num1("^Nine of the amplicons")
p_large <- num1("^Thirteen of the amplicons")
grab_int <- function(x, pattern) {
  m <- regmatches(x, regexpr(pattern, x, perl = TRUE))
  if (length(m) != 1L) stop("expected one match for: ", pattern)
  as.integer(regmatches(m, regexpr("[0-9]+", m)))
}
# Each count appears TWICE in the fraction (numerator and inside the
# parenthetical denominator), so anchor on the word that follows the first
# occurrence rather than on a bare number.
true_of <- function(x, word) {
  grab_int(x, paste0("[0-9]+(?= ", word, ")"))
}
# The group sizes are stated as number words ("Nine of the amplicons were
# <400 bp"), so map the word; the digit does not appear in that clause.
word_n <- c(Nine = 9L, Ten = 10L, Eleven = 11L, Twelve = 12L,
            Thirteen = 13L, Fourteen = 14L, Fifteen = 15L)
word_value <- function(x, pattern) {
  m <- regmatches(x, regexpr(pattern, x, perl = TRUE))
  if (!length(m) || !m %in% names(word_n)) {
    stop("group size is not a known number word: ", pattern)
  }
  unname(word_n[[m]])
}
paper <- list(
  n_small = word_value(p_small, "^[A-Z][a-z]+(?= of the amplicons)"),
  n_large = word_value(p_large, "^[A-Z][a-z]+(?= of the amplicons)"),
  mean_small = grab_int(p_small, "average length of [0-9]+"),
  mean_large = grab_int(p_large, "average length of [0-9]+"),
  tp_small = true_of(p_small, "true positive"),
  fn_small = true_of(p_small, "false negative"),
  tn_small = true_of(p_small, "true negative"),
  fp_small = true_of(p_small, "false positive"),
  tp_large = true_of(p_large, "true positive"),
  fn_large = true_of(p_large, "false negative"),
  tn_large = true_of(p_large, "true negative"),
  fp_large = true_of(p_large, "false positive")
)
# The body's own words, not a regex artefact: the split is at 400 nt and the
# smaller group is the one counted as "<400".
stopifnot(grepl("<400 bp in length", p_small),
          grepl(">400 bp in length", p_large))

## --- 4. ASSERTION ------------------------------------------------------------
# The paper's two groups are "<400 bp" and ">400 bp"; an amplicon listed as
# exactly 400 belongs to neither as stated. Table 1 contains one such row, so
# the boundary cell is only identifiable through its arithmetic consequence.
amp$group <- ifelse(amp$size_bp < 400L, "lt400", "ge400")
obs_small <- amp[amp$group == "lt400", ]
obs_large <- amp[amp$group == "ge400", ]
obs_small_all <- amp[amp$size_bp <= 400L, ]
obs_large_all <- amp[amp$size_bp > 400L, ]

problems <- character(0)
if (nrow(obs_small) + nrow(obs_large) != nrow(amp)) {
  problems <- c(problems, "the two groups do not partition Table 1")
}
if (nrow(obs_small_all) != paper$n_small) {
  problems <- c(problems, sprintf(
    "<=400 bp: Table 1 gives %d amplicons, paper says %d",
    nrow(obs_small_all), paper$n_small))
}
if (nrow(obs_large_all) != paper$n_large) {
  problems <- c(problems, sprintf(
    ">400 bp: Table 1 gives %d amplicons, paper says %d",
    nrow(obs_large_all), paper$n_large))
}
at_400 <- amp$amplicon[amp$size_bp == 400L]
if (length(at_400) != 1L) {
  stop("expected exactly one amplicon at the 400 nt boundary, found ",
       length(at_400))
}
# The paper's split at exactly 400 cannot be identified from its own table: an
# amplicon listed as 400 is neither "<400" nor ">400", and BOTH readings are
# arithmetically admissible. Dropping it gives 8/14 and a 272 nt mean; gaining
# it gives 9/13 and 286 nt, which is the mean the paper prints. We record that
# as an inference, not as a checked fact.
mean_if_inclusive <- round(mean(obs_small_all$size_bp))
stopifnot(mean_if_inclusive == paper$mean_small)
# The parenthetical of each fraction repeats its two counts, e.g.
# "(38 true positive/(38 true positive + 0 false negative))" and
# "(339 true negative/(69 false positive + 339 true negative))" -- note the
# second orders the summands FP first, TN last. So read all three printed
# numbers and require them, as a SET, to be {n, n, d}: a cell read from the
# wrong paragraph cannot satisfy that.
# Operand order varies between the two parentheticals ("TP + FN" in one,
# "FP + TN" in the other), and so does which word carries the numerator. So:
# take the printed fraction "N <word>/(...)", pull every number out of it, and
# require the sorted set to be {n, n, denominator}.
fraction_ok <- function(x, num_word, n_num, d_other, d_num) {
  pat <- sprintf("\\([0-9]+ %s/\\([^)]*\\)\\)", num_word)
  m <- regmatches(x, regexpr(pat, x))
  if (!length(m) || !any(grepl(d_other, m))) {
    stop("fraction not found: ", num_word, " over ", d_other)
  }
  g <- sort(as.integer(regmatches(m, gregexpr("[0-9]+", m))[[1]]))
  identical(g, sort(c(n_num, n_num, d_num)))
}
frac_checks <- c(
  fraction_ok(p_small, "true positive", paper$tp_small, "false negative",
              paper$fn_small),
  fraction_ok(p_small, "true negative", paper$tn_small, "false positive",
              paper$fp_small),
  fraction_ok(p_large, "true positive", paper$tp_large, "false negative",
              paper$fn_large),
  fraction_ok(p_large, "true negative", paper$tn_large, "false positive",
              paper$fp_large))
if (!all(frac_checks)) {
  problems <- c(problems, "a printed fraction does not repeat the counts read")
}
# Re-derive the paper's rounded percentages from ITS OWN counts. The printed
# precision differs between paragraphs (0.83 but 0.81 and 0.84), so compare
# against the precision each paragraph actually shows.
sens_small <- paper$tp_small / (paper$tp_small + paper$fn_small)
spec_small <- paper$tn_small / (paper$tn_small + paper$fp_small)
sens_large <- paper$tp_large / (paper$tp_large + paper$fn_large)
spec_large <- paper$tn_large / (paper$tn_large + paper$fp_large)
fmt <- function(v, ref) {
  ## decimal places = digits after the point in the printed reference
  dp <- nchar(sub("^[01]\\.", "", ref))
  roundv <- round(as.numeric(v), dp)
  if (dp == 0L) sprintf("%d", as.integer(roundv))
  else sprintf(paste0("%.0f.", "%0", dp, "d"),
               floor(roundv), as.integer(round(roundv * 10^dp)) %% 10^dp)
}
printed_pct <- function(x, i = 1L) {
  m <- regmatches(x, gregexpr("[01]\\.[0-9]+", x))[[1]]
  m[i]
}
pct_pairs <- list(
  small_sens = c(n = fmt(sens_small, printed_pct(p_small, 1L)),
                 printed = printed_pct(p_small, 1L)),
  small_spec = c(n = fmt(spec_small, printed_pct(p_small, 2L)),
                 printed = printed_pct(p_small, 2L)),
  large_sens = c(n = fmt(sens_large, printed_pct(p_large, 1L)),
                 printed = printed_pct(p_large, 1L)),
  large_spec = c(n = fmt(spec_large, printed_pct(p_large, 2L)),
                 printed = printed_pct(p_large, 2L)))
pct_ok <- vapply(pct_pairs, function(x) x[["n"]] == x[["printed"]], TRUE)
if (!all(pct_ok)) {
  # A 1-2 point disagreement is a RECORDED CONTRADICTION, not a parse error
  # (the counts above are independently re-read from the printed fractions).
  # A disagreement larger than 2 points would mean we parsed the wrong cells.
  delta <- abs(vapply(pct_pairs, function(x) {
    as.numeric(x[["n"]]) - as.numeric(x[["printed"]])
  }, numeric(1)))
  if (any(delta > 0.02)) {
    print(pct_pairs)
    problems <- c(problems, "a printed percentage is >2 points from its counts")
  } else {
    message("desilva2009: WARNING -- printed percentages differ from counts by ",
            "<=2 points (recorded, not resolved, in SOURCES-arms.md): ",
            paste(names(pct_pairs)[!pct_ok], collapse = ", "))
  }
}
if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("desilva2009: assertion failed; nothing written")
}


# Recorded, unresolved: the mean for the >400 bp group. Our sizes give 545,
# the paper prints 544. The boundary cell is fixed by the <400 bp mean, so
# this 1 bp difference is between two printed numbers of one paper; it is
# left visible, not rounded away.
lt400_cell <- amp$size_bp <= 400L
mean_large <- round(mean(amp$size_bp[!lt400_cell]))
mean_large_printed <- paper$mean_large
mean_large_matches <- identical(mean_large, mean_large_printed)

## --- 5. Table 2: qualitative design-approach table ---------------------------
t2 <- tab(doc, "Table 2")
t2 <- Filter(function(r) length(r) == 2L, t2)
if (t2[[1]][1] != "Genetic situation") stop("Table 2: unexpected first row")
t2 <- t2[-1]
strat <- data.frame(
  genetic_situation = vapply(t2, `[`, character(1), 1L),
  hrm_solution = vapply(t2, `[`, character(1), 2L),
  stringsAsFactors = FALSE
)
if (nrow(strat) != 5L) stop("Table 2 did not parse into 5 strategy rows")

## --- 6. write ---------------------------------------------------------------
# The two fetched tables are INDEPENDENT assay sets, not two views of one set:
# they share a forward primer (CTGGAAGACTCCAGGTCAGG) but with a different
# reverse primer, and no forward+reverse pair occurs in both. So a size in
# Table 1 is not corroborated by Additional file 1 and vice versa. Assert the
# independence so nobody later "cross-validates" them.
amps_of <- function(d) paste(d$forward, d$reverse, sep = "+")
if (length(intersect(amps_of(amp), amps_of(snp))) != 0L) {
  stop("Table 1 and Additional file 1 unexpectedly share an amplicon pair; "
       , "the independence statement in this script is now wrong")
}

amp$group <- ifelse(amp$size_bp < 400L, "lt400",
                    ifelse(amp$size_bp > 400L, "gt400",
                           "at400_assigned_lt400_by_mean"))
amp_out <- amp[, c("amplicon", "forward", "reverse", "location", "size_bp",
                   "group")]
amp_out$primer_total_nt <- nchar(amp$forward) + nchar(amp$reverse)
amp_out$doi <- DOI
write.csv(amp_out, file.path(out_dir, "desilva2009_amplicons.csv"),
          row.names = FALSE, na = "")

snp$doi <- DOI
write.csv(snp, file.path(out_dir, "desilva2009_snps.csv"),
          row.names = FALSE, na = "")

small <- amp$size_bp[amp$group != "gt400"]
big <- amp$size_bp[amp$group == "gt400"]
outcomes <- data.frame(
  readout = "HRM mutation scanning",
  amplicon_group = c("lt400", "gt400"),
  length_rule = c("<400 bp in length", ">400 bp in length"),
  n_amplicons_in_table1 = c(length(small), length(big)),
  n_amplicons_paper = c(paper$n_small, paper$n_large),
  mean_bp_table1 = c(round(mean(small)), round(mean(big))),
  mean_bp_paper = c(paper$mean_small, paper$mean_large),
  mean_bp_agrees = c(TRUE, mean_large_matches),
  true_positive = c(paper$tp_small, paper$tp_large),
  false_negative = c(paper$fn_small, paper$fn_large),
  true_negative = c(paper$tn_small, paper$tn_large),
  false_positive = c(paper$fp_small, paper$fp_large),
  sensitivity = c(sens_small, sens_large),
  specificity = c(spec_small, spec_large),
  n_samples_paper = 47L,
  n_samples_implied = c(paper$tp_small + paper$fn_small,
                        paper$tp_large + paper$fn_large),
  stringsAsFactors = FALSE
)
outcomes$outcome_def <- paste(
  "HRM mutation scanning vs Sanger resequencing of the same 47 TP53 DNAs.",
  "Counts are POOLED over the amplicons of each group; the paper prints no",
  "per-amplicon contingency table, so no per-amplicon outcome exists here.",
  "These are mutation-scanning calls (variant present / no variant), not",
  "allele-specific primer discrimination and not a genotype cluster call.")
outcomes$pooled_not_per_amplicon <- TRUE
outcomes$calibration <- paste(
  "uncalibrated: the 400 bp grouping is the authors' cut point,",
  "not a threshold we fitted or endorse")
outcomes$doi <- DOI
write.csv(outcomes, file.path(out_dir, "desilva2009_outcomes.csv"),
          row.names = FALSE, na = "")

strat$doi <- DOI
write.csv(strat, file.path(out_dir, "desilva2009_strategies.csv"),
          row.names = FALSE, na = "")

message(
  "desilva2009: ", nrow(amp), " Table 1 amplicons; <=400 bp n=",
  length(small), " mean ", round(mean(small)), " bp (paper ", paper$mean_small,
  "), >400 bp n=", length(big), " mean ", round(mean(big)), " bp (paper ",
  paper$mean_large, " agree=", mean_large_matches, "); ",
  sum(!is.na(snp$size_bp)), " sized of ", nrow(snp),
  " per-SNP primer rows (", sum(snp$multi_snp_amplicon),
  " multi-SNP amplicons); ", nrow(strat), " Table 2 strategies; ",
  nrow(outcomes), " pooled outcome rows (",
  paste(outcomes$n_samples_implied, collapse = "/"), " samples)"
)
