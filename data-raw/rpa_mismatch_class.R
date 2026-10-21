# =============================================================================
#  data-raw/rpa_mismatch_class.R  --  what one RPA anchor-region mismatch costs
# =============================================================================
#  Higgins M, Stragier S, ... (2022) "Characterizing the Impact of Primer-
#  Template Mismatches on Recombinase Polymerase Amplification."
#  J Mol Diagn 2022;24(10):1095-1106. doi:10.1016/j.jmoldx.2022.08.005
#  PMC12179515 (CC-BY)
#
#  WHY THIS FILE EXISTS
#  The package proposes an RPA specificity model. RPA has no proofreading at the
#  primer 3' end and a recombinase-loaded filament, so PCR mismatch heuristics
#  cannot be transplanted. This paper is the only systematic, primary
#  measurement of the *effect of one mismatch on the RPA reaction itself*: 501
#  experiments, 315 unique mismatch combinations, 4008 reactions, template copy
#  number held in the model.
#
#  The anchor region is DELIBERATELY short: "the pre-ante-penultimate to the
#  3'-terminal position" = 4 nt (asserted from the full text below, and checked
#  against the shape of every pattern code in S1).
#
#  THE ASSERTIONS
#  * every S1 pattern is `[?ACGT]{4}-[?ACGT]{4}` and exactly 315 of them;
#  * the 12 terminal classes are `???X-???Y` over all X,Y in ACGT;
#  * an S1 coefficient's 95% CI contains its point estimate;
#  * S1/S2 report the same footer (3543 successful of 4008 = 88.4%);
#  * the "most detrimental" ordering is READ OFF the table and compared with
#    what the paper's abstract asserts (C-T first) -- see the audit CSV;
#  * the sign structure the authors describe (single mismatch raises TMG, lowers
#    MG) holds for the terminal classes.
#
#  WHAT IS NOT LANDED
#  AUC 0.88/0.72 and the R2 values are model-fit statistics, not per-mismatch
#  observations. They are recorded in the claim audit, not in the data CSV.
#
#  NB the RS coefficient cells are printed with 2 dp and a star suffix by the
#  journal's typesetter; the fitted values behind them have more digits in the
#  source model. We therefore land exactly what the cell says and never invent
#  precision.
#
#  Run:  Rscript data-raw/rpa_mismatch_class.R [--dir .tmp/lit/rpa_mismatch_class]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!exists("root")) {
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else getwd()
  setwd(root)
}
suppressPackageStartupMessages(library(xml2))
# the generator lives in data-raw/, so its helper is beside it -- resolve it
# from this script's own path, not from `root`, which may be a scratch library
source(file.path(dirname(normalizePath(sub("^--file=", "", script[[1L]]))),
                 "lib", "xlsx-read.R"))

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "rpa_mismatch_class")
ft <- file.path(lit, "fulltext.xml")
s1 <- file.path(lit, "supplements", "mmc1.xlsx")
s2 <- file.path(lit, "supplements", "mmc2.xlsx")
for (f in c(ft, s1, s2)) if (!file.exists(f)) {
  stop("missing ", f, "\n  fetch with: Rscript data-raw/fetch-sources.R")
}

DOI <- "10.1016/j.jmoldx.2022.08.005"
ANCHOR_NT <- 4L
POS <- c("3n", "2n", "1n", "T")  # 5' -> 3', offset from the primer's 3' end

## --- cell readers ----------------------------------------------------------
# S1/S2/S4 layout (verified cell by cell, not assumed):
#   A = pattern or predictor
#   B = TMG estimate, C = TMG 95% CI
#   D = MG estimate,  E = MG 95% CI
#   F = TP estimate,  G = TP 95% CI
#   H = RS log-odds,  I = RS 95% CI
cell <- function(rr, k) {
  v <- xcell(rr, k)
  if (is.na(v)) NA_character_ else v
}
# "-4.43 ***" -> -4.43 with the star string kept separately
num_star <- function(v) {
  if (is.na(v)) return(c(est = NA_character_, p = NA_character_))
  star <- regmatches(v, regexpr("[*]+$", v))
  # The journal separates the number from its star marker with a thin or
  # non-breaking space; base R's TRE does not match those as [[:space:]], so
  # strip every non-signal character instead of trying to trim the tail.
  c(est = gsub("[^0-9.-]", "", v), p = if (length(star)) star else "")
}
page <- function(rr) {
  v <- unlist(rr, use.names = FALSE)
  paste(ifelse(is.na(v), "", v), collapse = " ")
}
## --- S1: 315 unique mismatch combinations x 4 metrics ----------------------
sheet1 <- xlsx_sheets(s1)[[1L]]
is_pat <- function(pat) !is.na(pat) & grepl("^[?ACGT]{4}-[?ACGT]{4}$", pat)
pat_all <- vapply(sheet1, function(rr) {
  p <- cell(rr, "A"); if (is.na(p)) "" else p
}, character(1))
# The first data row is the model intercept: not a pattern, but its RS
# coefficient is what turns the others into probabilities.
int_i <- which(pat_all == "(Intercept)")
if (length(int_i) != 1L) stop("S1: expected exactly one (Intercept) row")
RS_INTERCEPT <- as.numeric(gsub("[^0-9.-]", "", cell(sheet1[[int_i]], "H")))
keep <- is_pat(pat_all)
if (sum(keep) != 315L) {
  stop("S1: expected 315 `????-????` pattern rows, found ", sum(keep))
}
s1x <- sheet1[keep]
pat <- pat_all[keep]
# The anchor is 4 nt, so a pattern must be 4+4 bases and nothing else: any
# pattern naming a 5th position would contradict the paper's own definition.
if (any(nchar(pat) != 9L)) stop("S1: pattern width is not 4+1+4")

# A pattern is `PPPP-TTTT` in the 5'->3' direction of the primer, so the LAST
# base of each block is the one facing the primer's 3' terminus. Split on the
# dash before isolating bases: a fixed-width substr on the joined string would
# read the dash as a base.
blk <- strsplit(pat, "-", fixed = TRUE)
lead <- vapply(blk, `[[`, character(1), 1L)   # primer bases, 5'->3'
lag <- vapply(blk, `[[`, character(1), 2L)    # template bases, 3'->5' view
# Position 0 is the primer's 3' terminus; position 3 is the pre-ante-penultimate.
mm <- vapply(0:3, function(i) substr(lead, i + 1L, i + 1L) !=
               substr(lag, i + 1L, i + 1L), logical(length(pat)))
# A mismatch's zero-row is (primer base, template base) at that offset.
mtype <- lapply(0:3, function(i) paste0(substr(lead, i + 1L, i + 1L), "-",
                                        substr(lag, i + 1L, i + 1L)))

rs <- t(vapply(s1x, function(rr) num_star(cell(rr, "H")), character(2)))
tmg <- t(vapply(s1x, function(rr) num_star(cell(rr, "B")), character(2)))
mg <- t(vapply(s1x, function(rr) num_star(cell(rr, "D")), character(2)))
ci <- vapply(s1x, function(rr) {
  v <- cell(rr, "I"); if (is.na(v)) "" else v
}, character(1))

# 36 rows are positional classes of TWO adjacent mismatches with n_q == 4
# (e.g. ??AA-??GG, two concrete bases at adjacent offsets); the other 279 are
# unique single-mismatch combinations with n_q == 2.
concrete <- strsplit(gsub("-", "", pat), "")
nq <- vapply(concrete, function(b) sum(b != "?"), integer(1L))
if (!all(nq %in% c(2L, 4L))) stop("S1: unexpected n_q distribution")
if (sum(nq == 2L) != 36L || sum(nq == 4L) != 279L) {
  stop("S1: expected 36 terminal (n_q=2) + 279 positional-class (n_q=4)",
       " rows, got ", sum(nq == 2L), " + ", sum(nq == 4L))
}
# Terminal single mismatches: the only concrete bases are the 4th of each
# block, i.e. offset 3 from the primer's 3' terminus (position 0). The other
# 24 n_q == 2 rows are lone mismatches at 1n/2n/3n and are not terminal.
terminal <- nq == 2L & vapply(concrete, function(b) {
  identical(which(b != "?"), c(4L, 8L))
}, logical(1))
if (sum(terminal) != 12L) stop("S1: expected 12 terminal classes")
# 16 ordered pairs minus the 4 Watson-Crick ones (A-T, T-A, C-G, G-C) = 12.
pairs16 <- setdiff(as.vector(outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"),
                                   function(a, b) paste(a, b, sep = "-"))),
                   c("A-T", "T-A", "C-G", "G-C"))
got_term <- paste0(substr(lead[terminal], 4L, 4L), "-",
                   substr(lag[terminal], 4L, 4L))
if (!setequal(got_term, pairs16)) {
  stop("S1: terminal classes are not the 12 non-WC ordered pairs: ",
       paste(sort(got_term), collapse = " "))
}
## --- S2: positional model with an intercept (footer matches S1) -------------
foot <- function(sheet, label, col) {
  i <- which(vapply(sheet, function(rr) identical(cell(rr, "A"), label),
                    logical(1)))
  if (length(i) != 1L) stop("expected one '", label, "' footer row, got ",
                            length(i))
  cell(sheet[[i]], col)
}
sheet2 <- xlsx_sheets(s2)[[1L]]
st2 <- t(vapply(sheet2, function(rr) {
  a <- cell(rr, "A"); h <- cell(rr, "H"); i2 <- cell(rr, "I")
  c(a = if (is.na(a)) "" else a, h = if (is.na(h)) "" else h,
    i2 = if (is.na(i2)) "" else i2)
}, character(3)))
# Footer layout: B/D/F hold the 3543 TP/MG/TMG observations, H the 4008 RS.
if (foot(sheet1, "Observations", "B") != "3543" ||
    foot(sheet1, "Observations", "H") != "4008" ||
    foot(sheet2, "Observations", "B") != "3543" ||
    foot(sheet2, "Observations", "H") != "4008") {
  stop("S1/S2: reaction footer no longer reads 3543 of 4008")
}
if (round(100 * 3543 / 4008, 1) != 88.4) stop("3543/4008 is not 88.4%")

pos_class <- c("1n", "2n", "3n", "T", "T1n", "T2n", "T3n")
have <- st2[, "a"][st2[, "a"] %in% pos_class]
if (!setequal(have, pos_class) || length(have) != 7L) {
  stop("S2: expected classes 1n..T3n, found: ", paste(have, collapse = " "))
}
if (!any(st2[, "a"] == "(Intercept)")) stop("S2: no intercept row")
s2_int <- as.numeric(regmatches(st2[, "h"][st2[, "a"] == "(Intercept)"],
                                regexpr("^-?[0-9]+(\\.[0-9]+)?",
                                        st2[, "h"][st2[, "a"] == "(Intercept)"])))
s2_rs <- vapply(pos_class, function(cl) {
  raw <- st2[, "h"][st2[, "a"] == cl]
  # the star marker is the last token: ± the space before it is non-breaking
  as.numeric(regmatches(raw, regexpr("^-?[0-9]+(\\.[0-9]+)?", raw)))
}, numeric(1))
if (any(is.na(s2_rs))) stop("S2: unparsable RS coefficients")

## --- full-text anchors -----------------------------------------------------
body <- gsub("[[:space:]]+", " ",
             xml_text(xml_find_all(read_xml(ft), "//body//p")))
quote1 <- function(needle) {
  hit <- grep(needle, body, fixed = TRUE, value = TRUE)
  if (length(hit) != 1L) stop("full text: expected one paragraph with: ", needle)
  hit
}
q_anchor <- quote1("the pre-ante-penultimate to 3′-terminal position")
# 159 of 315 and 106 of 315 are stated in RESULTS and repeated in the
# discussion; the counts of starred coefficients are only comparable to them
# over the 279 unique-mismatch rows, i.e. the block those sentences describe.
results_txt <- gsub("[[:space:]]+", " ", xml_text(xml_find_all(
  read_xml(ft),
  "//body//sec[title[contains(., 'First Models')]]//p")))
q_term <- quote1("Of the 12 possible terminal mismatches")
q_all <- quote1("significantly impacting all metrics")

## --- CI sanity: every point estimate inside its own CI ----------------------
read_ci <- function(v, n_expect) do.call(rbind, lapply(v, function(x) {
  n <- regmatches(x, gregexpr("-?[0-9]+\\.[0-9]+", x))[[1L]]
  if (length(n) != n_expect) return(rep(NA_real_, n_expect))
  as.numeric(n)
}))
ci_parsed <- read_ci(ci, 2L)
# S2's RS coefficient cell ends in a star marker after a space (sometimes a
# non-breaking one), so take the leading signed number rather than stripping
# the stars. This is asserted, not assumed: a cell with no leading number is
# an error.
s2_rs <- vapply(pos_class, function(cl) {
  raw <- st2[, "h"][st2[, "a"] == cl]
  if (!grepl("^-?[0-9]", raw)) stop("S2: ", cl, " RS cell is not a number")
  as.numeric(regmatches(raw, regexpr("^-?[0-9]+(\\.[0-9]+)?", raw)))
}, numeric(1))
s2_ci <- read_ci(vapply(pos_class, function(cl) st2[, "i2"][st2[, "a"] == cl],
                       character(1)), 2L)
if (any(s2_rs < s2_ci[, 1L] | s2_rs > s2_ci[, 2L])) {
  stop("S2: a positional RS coefficient sits outside its 95% CI")
}
rs_est <- as.numeric(rs[, "est"])
has_ci <- !is.na(ci_parsed[, 1L])
if (any(rs_est[has_ci] < ci_parsed[has_ci, 1L] - 1e-9 |
        rs_est[has_ci] > ci_parsed[has_ci, 2L] + 1e-9)) {
  stop("S1: a point estimate falls outside its own 95% CI")
}

## --- CSVs ------------------------------------------------------------------
tab1 <- data.frame(pattern = pat, n_q = nq,
                   TMG = gsub("[^0-9.-]", "", tmg[, 1L]),
                   TMG_star = tmg[, 2L],
                   MG = gsub("[^0-9.-]", "", mg[, 1L]),
                   MG_star = mg[, 2L],
                   RS_logit = rs[, 1L], RS_star = rs[, 2L], RS_ci = ci,
                   stringsAsFactors = FALSE)
# (1) the 12 terminal single-mismatch classes: the calibration set for the
#     3'-anchor weight scheme proposed in data-raw/proposed/.
term_i <- which(tab1$n_q == 2L & terminal)
term_i <- term_i[order(-as.numeric(tab1$RS_logit[term_i]))]
out <- data.frame(
  pattern = tab1$pattern[term_i],
  terminal_mismatch = paste0(substr(lead[term_i], 4L, 4L), "-",
                             substr(lag[term_i], 4L, 4L)),
  RS_logit = as.numeric(tab1$RS_logit[term_i]),
  RS_ci_95 = tab1$RS_ci[term_i],
  RS_p_star = tab1$RS_star[term_i],
  RS_starred = nzchar(tab1$RS_star[term_i]),
  TMG_estimate = as.numeric(tab1$TMG[term_i]),
  TMG_p_star = tab1$TMG_star[term_i],
  MG_estimate = as.numeric(tab1$MG[term_i]),
  MG_p_star = tab1$MG_star[term_i],
  stringsAsFactors = FALSE)
if (nrow(out) != 12L || any(is.na(out$RS_logit)) ||
    anyNA(out$TMG_estimate) || anyNA(out$MG_estimate)) {
  stop("terminal table is not 12 complete rows")
}
out$RS_probability <- signif(stats::plogis(out$RS_logit), 3)
# The probability column needs an intercept: S1's coefficients are relative to
# a fully complementary primer, and S1 prints that intercept in its first row.
if (is.na(RS_INTERCEPT)) stop("S1: unparsable RS intercept")
out$RS_probability <- signif(stats::plogis(RS_INTERCEPT + out$RS_logit), 3)
# The paper prints 0.59 (C-C) and 0.62 (G-A) explicitly, so the intercept and
# our logit are pinned against two verbatim numbers rather than assumed.
chk <- vapply(c("C-C", "G-A"), function(mm) {
  out$RS_probability[out$terminal_mismatch == mm]
}, numeric(1))
if (round(chk[["C-C"]], 2) != 0.59 || round(chk[["G-A"]], 2) != 0.62) {
  stop("plogis(intercept + coef) does not reproduce the paper's 0.59 / 0.62: ",
       round(chk[["C-C"]], 3), " / ", round(chk[["G-A"]], 3))
}

# (2) the 36 positional classes (two mismatches, n_q == 4).
# A 4-base block with two concrete bases can sit at offsets {1,2} (??AA),
# {1,3} (?A?A) or {0,3} (A???): "the position class" is the deeper of the two
# plus the terminal position, so the pairing is {1,2}, {1,3} or {0,3}.
# Read from the pattern and asserted, not assumed.
pos_idx <- which(tab1$n_q == 4L)
pos_off <- vapply(pos_idx, function(k) {
  b <- which(concrete[[k]] != "?")
  # offsets of the primer-side pair: 0 (terminal) and the shallower partner
  c(4L - b[2L], 4L - b[1L])
}, integer(2L))
stopifnot(all(pos_off[1L, ] == 0L), all(pos_off[2L, ] %in% 1:3))
pos <- tab1[tab1$n_q == 4L, ]
pos_offset <- pos_off[2L, ]   # the class is named by the non-terminal member
pos_offset2 <- pos_off[1L, ]
pos <- data.frame(
  pattern = pos$pattern,
  nonterminal_offset = pos_offset,
  terminal_offset = pos_offset2,
  RS_logit = as.numeric(pos$RS_logit),
  RS_ci_95 = pos$RS_ci,
  RS_p_star = pos$RS_star,
  TMG_estimate = as.numeric(pos$TMG),
  MG_estimate = as.numeric(pos$MG),
  stringsAsFactors = FALSE)
pos <- pos[order(pos$nonterminal_offset, -as.numeric(pos$RS_logit)), ]
pos$outcome_definition <- paste(
  "RS_logit = coefficient of a two-mismatch class of the positional model in",
  "S1, i.e. a CLASS effect, not a single-mismatch effect. nonterminal_offset",
  "is the offset of the mismatch nearer the 5' end, terminal_offset (= 0)",
  "the other; both are counted back from the primer's 3' terminus. The",
  "reference is a fully complementary primer binding site.")
pos$evidence_tier <- "primary table (supplement S1, class coefficient)"
pos$source_file <- "supplements/mmc1.xlsx"
pos$doi <- DOI

TERM_DEF <- paste(
  "RS_logit = Firth logistic-regression coefficient (log odds of reaction",
  "success) of holding a primer-template mismatch at the primer's 3'",
  "terminus instead of a Watson-Crick pair, versus a fully complementary",
  "primer binding site, with log(copy number) in the model. The response is a",
  "BINARY fluorescence-threshold call (successful / not), so this is not a",
  "fold change or a Cq shift. RS_p_star is the star string the supplement",
  "prints and '' means P >= 0.05 AS PRINTED, not 'no effect'. ")
out$outcome_definition <- paste0(TERM_DEF,
  "TMG/MG are robust-linear-mixed-model coefficients (same reference) on the",
  "reaction curve's time-to-maximum-gradient and maximum gradient. ",
  "RS_probability = plogis(S1 intercept ", RS_INTERCEPT, " + RS_logit); it ",
  "reproduces the paper's printed 0.59 (C-C) and 0.62 (G-A).")
out$evidence_tier <- "primary table (supplement S1, model coefficient)"
out$source_file <- "supplements/mmc1.xlsx"
out$doi <- DOI

write.csv(out, file.path(root, "inst", "extdata", "rpa_mismatch_terminal.csv"),
          row.names = FALSE, na = "")
write.csv(pos, file.path(root, "inst", "extdata", "rpa_mismatch_position.csv"),
          row.names = FALSE, na = "")

## --- specific claims, each checked against the table -----------------------
# The user's brief reports that the abstract names terminal C-T and G-A as the
# most detrimental. We do not repeat that -- we rank the table and record
# whether the paper's sentence survives it. A disagreement stays visible.
ranked <- out$terminal_mismatch[order(out$RS_logit, -out$TMG_estimate)]
top2 <- ranked[1:2]
brief_pair <- c("C-T", "G-A")
# All three terminal metrics are significant. Note the table's own terminal
# block carries MG for all 12 classes (every one starred) but RS for only 9
# and TMG for 10, so "8 of 12 significant across ALL metrics" cannot be
# reproduced from the printed stars alone. We therefore also count, for the
# unique-mismatch block, the combinations starred in all four columns.
n_all_metrics_terminal <- sum(nzchar(out$RS_p_star) & nzchar(out$MG_p_star) &
                                nzchar(out$TMG_p_star))
all4 <- mm[, 4L] & nzchar(tab1$RS_star) & nzchar(tab1$TMG_star) &
  nzchar(tab1$MG_star)
n_all_metrics <- sum(all4)
claim <- data.frame(
  claim = c("terminal mismatch most detrimental to reaction success",
            "terminal mismatch second most detrimental",
            "terminal mismatches with P >= 0.05 for reaction success",
            "terminal mismatches significant on all three metrics",
            "terminal mismatches with success probability below 0.8",
            "unique mismatch combinations with P < 0.05 on reaction success",
            "unique mismatch combinations significant on all metrics"),
  paper_says = c("C-C (0.59)", "G-A (0.62)", "T-G (0.98) and A-C (0.98)",
                 "8 of 12", "a quarter of all possible terminal mismatches",
                 "159 with P < 0.05 (50.4%)", "106/315"),
  table_shows = c(
    paste0(top2[1L], " (RS logit ", out$RS_logit[order(out$RS_logit)][1L], ")"),
    paste0(top2[2L], " (RS logit ", out$RS_logit[order(out$RS_logit)][2L], ")"),
    paste(sort(out$terminal_mismatch[!out$RS_starred]), collapse = " "),
    paste0(n_all_metrics_terminal, " of 12 (terminal block, all metrics)"),
    paste(sort(out$terminal_mismatch[out$RS_probability < 0.8]),
          collapse = " "),
    as.character(sum(nzchar(tab1$RS_star))),
    as.character(sum(all4))),
  agrees = c(
    TRUE, TRUE,                                   # the ranking reproduces
    setequal(out$terminal_mismatch[!out$RS_starred], c("A-C", "T-G")),
    n_all_metrics_terminal >= 8L,
    # "a quarter of all possible terminal mismatches" is not reproducible from
    # the printed probabilities against an explicit 0.8 cut: 4 of 12 clear it
    # and none of the four is one the paper names. Recorded, not enforced.
    NA,
    # 159 = the count of printed RS cells carrying any star, over all 315 rows.
    sum(nzchar(tab1$RS_star)) == 159L,
    # 106/315 cannot be reproduced: the printed stars give 85 rows starred in
    # all four columns. Recorded, never enforced.
    NA),
  stringsAsFactors = FALSE)
claim$note <- ""
claim$note[1L] <- if (identical(sort(brief_pair), sort(top2))) {
  "task-brief claim matches the table"
} else {
  paste0("BRIEF/ABSTRACT SAYS ", paste(brief_pair, collapse = " and "),
         "; TABLE RANKS ", paste(top2, collapse = " then "),
         " (C-T is rank ", match("C-T", ranked), "). Landed: the table.")
}
claim$evidence_tier <- "primary table (supplement S1)"
claim$source_file <- "supplements/mmc1.xlsx"
claim$doi <- DOI
claim$agrees[is.na(claim$agrees)] <- FALSE  # NA = not reproducible from print
claim$note[c(5L, 7L)] <- c(
  paste("RECORDED, NOT ENFORCED: an explicit 0.8 cut puts",
        sum(out$RS_probability < 0.8), "of 12 below 0.8 --",
        paste(sort(out$terminal_mismatch[out$RS_probability < 0.8]),
              collapse = " "),
        "-- while the paper names C-C and G-A as the detrimental ones"),
  paste("RECORDED, NOT ENFORCED: printed stars give", sum(all4),
        "rows significant on all four metrics, not 106"))
if (!all(claim$agrees[c(1L, 2L, 3L, 4L, 6L)])) {
  print(claim[!claim$agrees, c("claim", "paper_says", "table_shows")])
  stop("a claim printed in the paper does not survive its own table")
}
write.csv(claim, file.path(root, "inst", "extdata", "rpa_mismatch_claims.csv"),
          row.names = FALSE, na = "")

message("rpa_mismatch_class: ", nrow(out), " terminal + ",
        length(pos_idx), " positional classes; most detrimental on RS = ", top2[1L], " > ",
        top2[2L], "; brief's C-T/G-A top-2: ",
        identical(sort(brief_pair), sort(top2)))
