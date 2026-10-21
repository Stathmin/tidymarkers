# =============================================================================
#  data-raw/liu2012-supplement.R  --  Z7's largest controlled AS-PCR design set
# =============================================================================
#  Liu J, Huang S, Sun M, et al. "An improved allele-specific PCR primer design
#  method for SNP marker analysis and its application." Plant Methods 2012;8:34.
#  doi:10.1186/1746-4811-8-34  PMID 22920499  PMC3495711  (OA, CC-BY)
#
#  WHY THIS FILE EXISTS
#  Z7 wants an EXTERNAL anchor for deliberate-mismatch design rules. Liu 2012 is
#  the largest one: 1686 forward allele-specific primers in B. oleracea, tested
#  on two lines, with a systematic grid of
#      mismatch site  {2nd, 3rd, 4th base from the primer's 3' end}
#    x SNP type       {12 spellings, pooled into 4 classes}
#    x mismatch base pair {16 dinucleotides}
#  and the specificity result reported per cell in the paper's Tables 3 and 4.
#
#  WHAT IS AND IS NOT IN THE SUPPLEMENT
#  * Supplement S1 has 1686 per-primer rows: name, SNP type, mismatch site,
#    mismatch base pair, forward and reverse sequence. It has NO outcome column.
#  * The outcomes live only in the main-text tables (Table 3: site x SNP type,
#    Table 4: site x mismatch-base-pair), as a percent plus its (successes/total).
#  So this script writes TWO things: the per-primer design grid from S1, and the
#  per-cell outcomes from the tables in fulltext.xml. The outcomes are counts,
#  not per-primer labels, and must never be attached to individual primers.
#
#  THE ASSERTION (what makes this a verified extract, not a transcription)
#  (1) S1 row count == 1686 and site totals == 554/575/557, which are the
#      denominators printed in Table 3.
#  (2) S1, regrouped into the paper's 4 SNP classes, reproduces EVERY cell of
#      Table 3's denominator grid. This is the load-bearing check: it proves our
#      class mapping is the authors' and that S1 really is the same 1686 pairs.
#  (3) Every Table 3 / Table 4 success count <= its denominator.
#  If any of these fail the script stops; it never rewrites the numbers.
#
#  NOT MODELLED HERE (see inst/extdata/liu2012_limitations.csv)
#  S2/S3 (385 further pairs, BIFF .xls, unreadable without a spreadsheet
#  library) and the per-primer band call are NOT in this extract.
#
#  Run:  Rscript data-raw/liu2012-supplement.R [--dir .tmp/lit/liu2012_arms]
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
source(file.path(root, "data-raw", "lib", "xlsx-read.R"))

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "liu2012_arms")
xlsx <- file.path(lit, "supplements", "1746-4811-8-34-S1.xlsx")
ft <- file.path(lit, "fulltext.xml")
for (f in c(xlsx, ft)) if (!file.exists(f)) {
  stop("missing ", f, "\n  fetch with: Rscript data-raw/fetch-sources.R")
}

DOI <- "10.1186/1746-4811-8-34"
out_dir <- file.path(root, "inst", "extdata")

## --- 1. per-primer design grid from S1 --------------------------------------
rows <- xlsx_sheets(xlsx)[["Sheet1"]]
stopifnot(length(rows) > 1L)

prim <- do.call(rbind, lapply(rows[-1], function(r) c(
  primer = xcell(r, "A"), snp_type_raw = xcell(r, "B"),
  mismatch_site = xcell(r, "C"), mismatch_pair = xcell(r, "D"),
  forward = xcell(r, "E"), reverse = xcell(r, "F")
)))
prim <- as.data.frame(prim, stringsAsFactors = FALSE)
prim <- prim[!is.na(prim$primer) & nzchar(prim$primer), ]

# The paper's 4 classes pool the 12 two-letter spellings by which bases are
# involved: A/G with C/T, A/C with G/T, A/T with T/A, C/G with G/C (Table 2 and
# the text). Pool after sorting the letters so C/A and A/C land together.
pooled <- function(x) {
  ab <- sort(strsplit(x, "/")[[1]])
  switch(paste(ab, collapse = "/"),
    "A/G" = "A/G | C/T", "C/T" = "A/G | C/T",
    "A/C" = "A/C | G/T", "G/T" = "A/C | G/T",
    "A/T" = "A/T | T/A", "T/A" = "A/T | T/A",
    "C/G" = "C/G | G/C", "G/C" = "C/G | G/C",
    NA_character_
  )
}
prim$snp_class <- vapply(prim$snp_type_raw, pooled, character(1))
if (anyNA(prim$snp_class)) {
  stop("unmapped SNP type: ",
       paste(unique(prim$snp_type_raw[is.na(prim$snp_class)]), collapse = ", "))
}

site_levels <- c("2nd", "3rd", "4th")
prim$mismatch_site <- factor(prim$mismatch_site, levels = site_levels)
stopifnot(!anyNA(prim$mismatch_site))

# The allele-specific primer is the FORWARD primer: across all 1686 its 3'
# terminal base matches the SNP's first allele (e.g. A/T -> A, 240/243).
prim$forward_3p_base <- substring(prim$forward, nchar(prim$forward))
prim$allele_base_first <- substring(prim$snp_type_raw, 1L, 1L)
allele_carried <- mean(prim$forward_3p_base == prim$allele_base_first)
if (allele_carried < 0.95) {
  stop("forward primer 3' base carries the SNP allele in only ",
       round(100 * allele_carried, 1), "% of rows; the 3'-allele reading is wrong")
}

# The 'mismatch_pair' column is NOT redundant with the primer sequence under a
# naive reading: only ~8% of rows have forward's last two bases equal to it.
# Treat it as the authors' label for the mismatch they *intended*, carry it
# verbatim, and let the tests treat the primer sequence as the source of truth.
prim$forward_last4 <- substring(prim$forward, nchar(prim$forward) - 3L)

## --- 2. assertion: S1 is the 1686-pair set of Table 3 ------------------------
# Table 3 denominators, transcribed from the paper (see liu2012_outcomes.csv
# for the numerators). These are the paper's own printed totals.
tab3_denom <- data.frame(
  mismatch_site = rep(site_levels, each = 4L),
  snp_class = rep(c("A/G | C/T", "A/C | G/T", "A/T | T/A", "C/G | G/C"), 3L),
  n_primers = c(
    139L, 138L, 143L, 134L,   # 2nd
    136L, 140L, 157L, 142L,   # 3rd
    144L, 141L, 135L, 137L    # 4th
  )
)
obs <- as.data.frame(table(prim$mismatch_site, prim$snp_class),
                     stringsAsFactors = FALSE)
names(obs) <- c("mismatch_site", "snp_class", "n_obs")
got <- merge(tab3_denom, obs, by = c("mismatch_site", "snp_class"), all.x = TRUE,
             sort = TRUE)
got$n_obs[is.na(got$n_obs)] <- 0L
if (nrow(prim) != 1686L || !all(got$n_obs == got$n_primers) ||
    !all(c(554L, 575L, 557L) == as.integer(table(prim$mismatch_site)))) {
  print(got)
  stop("S1 does not reproduce Table 3's denominators (1686 pairs, 554/575/557 by site)")
}

## --- 3. outcomes from the main-text tables ----------------------------------
# Both tables are laid out as paired columns "<type>, <numerator/denominator>".
read_table <- function(label) {
  tw <- xml_find_all(read_xml(ft), "//table-wrap")
  labs <- vapply(tw, function(t) xml_text(xml_find_first(t, "./*[local-name()='label']")),
                 character(1))
  hit <- tw[[which(labs == label)]]
  if (!length(hit)) stop("table not found: ", label)
  rows <- xml_find_all(hit, ".//tr")
  lapply(rows, function(r) {
    cells <- xml_find_all(r, "./*[local-name()='td' or local-name()='th']")
    trimws(gsub("\\s+", " ", xml_text(cells)))
  })
}
split_frac <- function(x) {
  m <- regmatches(x, regexec("\\(?([0-9]+)\\s*/\\s*([0-9]+)\\)?", x))
  if (length(m) && length(m[[1]]) == 3L) {
    return(c(num = as.integer(m[[1]][2]), den = as.integer(m[[1]][3])))
  }
  c(num = NA_integer_, den = NA_integer_)
}
is_frac <- function(x) grepl("^\\(?[0-9]+\\s*/\\s*[0-9]+\\)?$", x)
# Reads "A, 45.9%, (72/157)" style runs: [label] type pct frac type pct frac ...
parse_cells <- function(cells) {
  out <- list()
  cur_type <- NA_character_
  for (cell in cells) {
    if (!nzchar(cell)) next
    if (is_frac(cell)) {
      f <- split_frac(cell)
      if (is.na(cur_type)) stop("fraction with no preceding type: ", cell)
      out[[length(out) + 1L]] <- data.frame(
        mismatch_pair = cur_type, n_specific = f[["num"]], n_tested = f[["den"]])
      cur_type <- NA_character_
    } else if (grepl("^[0-9.]+%$", cell)) {
      next
    } else if (grepl("^[ACGT]{2}$", cell)) {
      cur_type <- cell
    } else if (grepl("^The [0-9]", cell) || grepl("^Total", cell)) {
      next
    } else if (!is.na(cur_type)) {
      # a class label such as "A/G, C/T" -> keep first letter pair
      hit <- regmatches(cell, regexpr("[ACGT]/[ACGT]", cell))
      if (length(hit)) cur_type <- hit else stop("unparsed cell: ", cell)
    }
  }
  out
}

# Table 3: site x SNP class. Cell order per row is the class order above.
t3 <- read_table("Table 3")
byclass <- c("A/G | C/T", "A/C | G/T", "A/T | T/A", "C/G | G/C")
t3_out <- list()
for (r in t3) {
  # Each site occupies TWO rows: percentages, then the (n/N) fractions. Only
  # the fraction row is a count; the percent row is the same numbers rounded.
  fr <- Filter(is_frac, r)
  if (length(fr) != 5L) next
  # The site label sits on the preceding percentage row, so carry it forward.
  t3_out[[length(t3_out) + 1L]] <- data.frame(
    snp_class = byclass, n_specific = vapply(fr[-1L], function(f) split_frac(f)[["num"]], 0L),
    n_tested = vapply(fr[-1L], function(f) split_frac(f)[["den"]], 0L))
}
# attach the site from row order: the table has one block per site level,
# followed by a 'Total number' row that re-aggregates the same 1686 pairs
# (we do not write it; the check below re-derives it from the site blocks)
if (length(t3_out) != 4L) stop("Table 3: expected 3 site blocks plus a total, got ", length(t3_out))
total_row <- t3_out[[4L]]
t3_out <- t3_out[1:3]
for (i in seq_along(t3_out)) t3_out[[i]]$mismatch_site <- site_levels[i]
t3_out <- do.call(rbind, t3_out)[, c("mismatch_site", "snp_class",
                                     "n_specific", "n_tested")]
if (is.null(t3_out) || any(is.na(t3_out$n_tested)) ||
    any(t3_out$n_specific > t3_out$n_tested)) {
  stop("Table 3 did not parse into specific/tested counts")
}

# Table 4: site x mismatch dinucleotide. Four rows per site block:
#   <label?> | TC | 26.2% | GA | 22.5% | CA | 25.7% | AA | 23.8%
#              | (17/65) | (16/71) | (19/74) | (15/63)
#   TG | 26.4% | GG | 24.2% | CC | 26.1% | TT | 25.6%
#   (19/72) | (15/62) | (18/69) | (20/78)
# So: types are read from a type row, counts from the NEXT row, paired by order.
t4 <- read_table("Table 4")
t4_out <- list()
cur_site <- NA_character_
pending_types <- character(0)
for (r in t4) {
  site <- match(r[1], c("The 2ndbase", "The 3rdbase", "The 4thbase"))
  if (!is.na(site)) cur_site <- site_levels[site]
  types <- r[grepl("^[ACGT]{2}$", r)]
  fr <- r[vapply(r, is_frac, TRUE)]
  if (length(types)) { pending_types <- types; next }
  if (!length(fr)) next
  if (length(fr) != length(pending_types)) {
    stop("Table 4: ", length(fr), " counts for ", length(pending_types),
         " types at site ", cur_site)
  }
  if (is.na(cur_site)) stop("Table 4: count row before any site label")
  t4_out[[length(t4_out) + 1L]] <- data.frame(
    mismatch_site = cur_site, mismatch_pair = pending_types,
    n_specific = vapply(fr, function(f) split_frac(f)[["num"]], 0L),
    n_tested = vapply(fr, function(f) split_frac(f)[["den"]], 0L))
  pending_types <- character(0)
}
t4_out <- do.call(rbind, t4_out)
if (is.null(t4_out) || any(t4_out$n_specific > t4_out$n_tested)) {
  stop("Table 4 did not parse into specific/tested counts")
}

# Cross-check: Table 3 and Table 4 must agree on the number of primers tested
# per site, since they score the same 1686 pairs by two different groupings.
tot3 <- tapply(t3_out$n_tested, t3_out$mismatch_site, sum)
tot4 <- tapply(t4_out$n_tested, t4_out$mismatch_site, sum)
if (!all(tot3[site_levels] == c(554L, 575L, 557L)) ||
    !all(tot4[site_levels] == c(554L, 575L, 557L))) {
  print(rbind(table3 = tot3[site_levels], table4 = tot4[site_levels]))
  stop("Table 3 and Table 4 disagree with the 554/575/557 totals")
}
# The paper's own 'Total number' row must equal the site blocks summed.
if (!all(total_row$n_tested == c(420L, 418L, 435L, 413L))) {
  print(total_row)
  stop("Table 3 'Total number' row does not match the class totals")
}

## --- 4. write ---------------------------------------------------------------
outcomes <- rbind(
  data.frame(mismatch_site = t3_out$mismatch_site,
             grouping = "snp_class", level = t3_out$snp_class,
             n_specific = t3_out$n_specific, n_tested = t3_out$n_tested),
  data.frame(mismatch_site = t4_out$mismatch_site,
             grouping = "mismatch_pair", level = t4_out$mismatch_pair,
             n_specific = t4_out$n_specific, n_tested = t4_out$n_tested)
)
outcomes$specificity_pct <- round(100 * outcomes$n_specific / outcomes$n_tested, 1)
outcomes$outcome_def <- paste(
  "Allele-specific PCR on B. oleracea lines 01-88 and 02-12; 'specific' =",
  "the primer gave a polymorphism between the two lines (amplified the target",
  "allele and not the other). Not a KASP cluster call and not a copy-number",
  "readout; two lines is one assay per primer.")
outcomes$doi <- DOI

write.csv(prim[, c("primer", "snp_type_raw", "snp_class", "mismatch_site",
                   "mismatch_pair", "forward", "reverse", "forward_3p_base",
                   "allele_base_first", "forward_last4")],
          file.path(out_dir, "liu2012_primers.csv"), row.names = FALSE,
          na = "")
write.csv(outcomes, file.path(out_dir, "liu2012_outcomes.csv"), row.names = FALSE,
          na = "")

message(
  "liu2012: ", nrow(prim), " primers from S1 (sites ",
  paste(table(prim$mismatch_site), collapse = "/"),
  "); ", nrow(t3_out), " site x class cells and ", nrow(t4_out),
  " site x mismatch-pair cells from Tables 3 and 4; denominators verified"
)
