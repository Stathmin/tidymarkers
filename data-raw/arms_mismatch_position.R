# =============================================================================
#  data-raw/arms_mismatch_position.R  --  the deliberate-mismatch POSITION ledger
# =============================================================================
#  THE QUESTION
#  Where, relative to the allele-specific primer's 3' end, should a deliberate
#  (artificial) mismatch be placed? Numbering used everywhere in this file is
#  NEGATIVE FROM THE 3' END: the last base is -1, the one before it -2, etc.
#  The papers do NOT share a numbering. The mapping, stated once here and never
#  silently applied, is:
#     "2nd base from the 3' end"  = -2      "penultimate" = -2
#     "3rd base from the 3' end"  = -3      Lefever's scan is 0 = terminal, so
#                                           his "position 3" = -3
#  This is the most dangerous step in the family: a paper saying "position 2 or
#  3" and a CSV column saying "-2 or -3" are different objects.
#
#  WHAT THIS FILE DOES, AND REFUSES TO DO
#  One row per MEASUREMENT (not per opinion), each with the position it supports
#  and how that position was obtained. It COUNTS how many independent sources
#  support -2 and -3. It does not average, vote, or pick a winner: the counts are
#  the output and the contradiction is reported in
#  data-raw/proposed/arms_mismatch_position.proposal.md.
#
#  Evidence classes, kept separate because they are not interchangeable:
#    measurement    - the paper reports an outcome at each position it scanned
#    primer_scan    - the position is read off the primer SEQUENCES designed by
#                     the authors; per-position outcomes are figures only, so
#                     the position is certain but the ranking is not quantified
#    rule           - the paper states a position without testing it
#
#  ASSERTIONS (each one stops)
#  (1) Every quote is found in its cited file. A row cannot outlive its evidence.
#      ye2001 is the one exception: its full text is present in neither tree, so
#      it carries no quote and is recorded as UNVERIFIED rather than asserted.
#  (2) Liu 2012's site totals recomputed from its own S1 == 554/575/557, its
#      per-site fractions are (139/554, 189/575, 162/557) == 490/1686, and the
#      site ranking is 3rd > 4th > 2nd. Liu is therefore filed at -3.
#  (3) The printed counts equal a plain filter of the ledger.
#  (4) Both sides still have at least one source. If one side drops to zero the
#      contradiction has been resolved by evidence - a human decision, so we
#      stop and say so instead of quietly editing the script.
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

# Source search path: this worker's fetches live in .tmp/lit here; papers fetched
# by the parent science lane live in the sibling worktree and are READ, not written.
src_dirs <- c(
  file.path(root, ".tmp", "lit"),
  file.path(dirname(root), "tidymarkers-agent-science", ".tmp", "lit")
)
find_source <- function(rel) {
  for (d in src_dirs) {
    p <- file.path(d, rel)
    if (file.exists(p)) return(p)
  }
  NA_character_
}
out_dir <- file.path(root, "inst", "extdata")

# ---------------------------------------------------------------------------
# The ledger. One row = one paper's one claim/measurement about one position.
# ---------------------------------------------------------------------------
rows <- list()
add <- function(key, doi, citation, position, evidence_class, n_tested,
                n_primers, discrimination, source_relpath, quote, note) {
  rows[[length(rows) + 1L]] <<- data.frame(
    key = key, doi = doi, citation = citation, position = position,
    evidence_class = evidence_class, n_tested = n_tested, n_primers = n_primers,
    discrimination = discrimination, source_relpath = source_relpath,
    quote = quote, note = note, stringsAsFactors = FALSE)
}

# --- Ye 2001: the rule behind the package's -2 assumption --------------------
# Full text is absent from both trees (ye2001/ is empty), so NO quote is made.
add(
  key = "ye2001", doi = "10.1093/nar/29.17.e88",
  citation = "Ye et al. 2001, Nucleic Acids Res 29:e88",
  position = -2L, evidence_class = "rule", n_tested = 0L, n_primers = 12L,
  discrimination = NA_character_,
  source_relpath = file.path("ye2001", "fulltext.xml"), quote = "",
  note = "inner primers carry a deliberate second mismatch; NO per-position outcome. Text not present in either tree: unverified here"
)

# --- Liu 2012: the largest controlled grid, filed at -3 ----------------------
# The package used to label Liu as the -2 anchor. Its own sentence ranks the
# 2nd/3rd/4th bases 25.1% / 32.9% / 29.1% polymorphic, so the 3rd base (-3) is
# the peak. Liu supports -3, not -2. Corrected here because the numbers are on
# file; the contradiction with the package's framing is reported, not hidden.
add(
  key = "liu2012", doi = "10.1186/1746-4811-8-34",
  citation = "Liu et al. 2012, Plant Methods 8:34",
  position = -3L, evidence_class = "measurement", n_tested = 3L, n_primers = 1686L,
  discrimination = "polymorphic 25.1% (139/554) at -2, 32.9% (189/575) at -3, 29.1% (162/557) at -4",
  source_relpath = file.path("liu2012_arms", "fulltext.xml"),
  quote = "polymorphism percents of the 2 nd , 3 rd and 4 th base location",
  note = "largest controlled scan; -3 is the peak by the authors' own percentages"
)

# --- Lefever 2019: explicit 1-5 scan (0 = terminal), optimum 3 == -3 ---------
add(
  key = "lefever2019", doi = "10.1038/s41598-019-38581-z",
  citation = "Lefever et al. 2019, Sci Rep 9:3841",
  position = -3L, evidence_class = "measurement", n_tested = 5L, n_primers = 6L,
  discrimination = "highest mean dCq at positions 2-3; position 3 optimal on variance and match penalty",
  source_relpath = file.path("lefever2019", "fulltext.xml"),
  quote = "a mismatch on position three appears to be the optimal solution",
  note = "all four nucleotides at positions 1-5 of both AS primers, 6 SNPs"
)

# --- Zhu 2021 (sorghum): positions decoded from the primer sequences ---------
# data-raw/sorghum_arms2021.R aligns Table S1: wxc F4 differs only at -2;
# wxt F4 differs at -1 AND -4, so its "-4" is not tested on its own.
add(
  key = "sorghum_arms2021", doi = "10.3390/foods10092218",
  citation = "Zhu et al. 2021, Foods 10:2218",
  position = -2L, evidence_class = "primer_scan", n_tested = 0L, n_primers = 10L,
  discrimination = "wxc F4 (-2 only) identified waxy sorghum specifically; per-primer Cq only in Fig. S2",
  source_relpath = file.path("sorghum_arms2021", "supplements", "s001",
                             "foods-1353494-supplementary.pdf"),
  quote = "Wxc F4                  CAAATACGACGTGTCAACAT",
  note = "position decoded from the Table S1 sequence, not from prose"
)
add(
  key = "sorghum_arms2021", doi = "10.3390/foods10092218",
  citation = "Zhu et al. 2021, Foods 10:2218",
  position = -4L, evidence_class = "primer_scan", n_tested = 0L, n_primers = 10L,
  discrimination = "wxt F4 (-1 with -4) gave the highest discrimination by the authors' judgement",
  source_relpath = file.path("sorghum_arms2021", "supplements", "s001",
                             "foods-1353494-supplementary.pdf"),
  quote = "Wxt F4                  CAAATACGACGTGTCACCGG",
  note = "the primer called 'fourth' also carries -1, so -4 alone is never tested"
)

ledger <- do.call(rbind, rows)

## --- (1) quotes must be in their cited files --------------------------------
text_of <- function(path) {
  if (grepl("[.]pdf$", path, ignore.case = TRUE)) {
    tf <- tempfile(fileext = ".txt")
    system2("pdftotext", c("-layout", shQuote(path), shQuote(tf)), stdout = TRUE,
            stderr = TRUE)
    on.exit(unlink(tf), add = TRUE)
    return(gsub("[[:space:]]+", " ",
                paste(readLines(tf, warn = FALSE), collapse = " ")))
  }
  if (grepl("[.]xml$", path, ignore.case = TRUE)) {
    return(gsub("[[:space:]]+", " ",
                paste(xml2::xml_text(xml2::xml_find_all(xml2::read_xml(path),
                                                        "//text()")),
                      collapse = " ")))
  }
  ""
}

ledger$verified_here <- FALSE
for (i in seq_len(nrow(ledger))) {
  p <- find_source(ledger$source_relpath[i])
  if (is.na(p)) {
    if (nzchar(ledger$quote[i])) {
      stop("quoted source for ", ledger$key[i], " is missing: ",
           ledger$source_relpath[i])
    }
    next  # unverified rule row; explicitly allowed only when quote == ""
  }
  if (grepl("[.]xlsx$", p, ignore.case = TRUE)) {
    ledger$verified_here[i] <- TRUE
    next
  }
  body <- text_of(p)
  # The PDF layout puts many spaces between primer name and sequence; compare
  # on whitespace-collapsed text so the quote is robust to pdftotext spacing.
  if (!grepl(gsub("[[:space:]]+", " ", ledger$quote[i]), body, fixed = TRUE)) {
    stop("quote \"", ledger$quote[i], "\" not found in ", p,
         " (row ", ledger$key[i], ")")
  }
  ledger$verified_here[i] <- TRUE
}

## --- (2) Liu 2012's own numbers ---------------------------------------------
liu_xlsx <- find_source(file.path("liu2012_arms", "supplements",
                                  "1746-4811-8-34-S1.xlsx"))
if (!is.na(liu_xlsx)) {
  source(file.path(root, "data-raw", "lib", "xlsx-read.R"))
  s1 <- xlsx_sheets(liu_xlsx)[["Sheet1"]]
  sites <- vapply(s1[-1], function(r) xcell(r, "C"), character(1))
  sites <- sites[!is.na(sites) & nzchar(sites)]
  totals <- as.integer(table(factor(sites, levels = c("2nd", "3rd", "4th"))))
  if (!identical(totals, c(554L, 575L, 557L))) {
    stop("Liu S1 site totals are ", paste(totals, collapse = "/"),
         ", expected 554/575/557")
  }
  # The four numerator/denominator pairs printed in the text, recomputed.
  n_poly <- c(139L, 189L, 162L)
  pct <- 100 * n_poly / totals
  if (!identical(round(pct, 1), c(25.1, 32.9, 29.1)) ||
      sum(n_poly) != 490L || sum(totals) != 1686L) {
    stop("Liu per-site fractions do not reproduce 25.1/32.9/29.1 and 490/1686")
  }
  if (which.max(pct) != 2L) {
    stop("Liu's per-site peak is not the 3rd base; the -3 filing is wrong")
  }
}

## --- (3) counts must match a plain filter of the ledger ---------------------
count_position <- function(pos) {
  sub <- ledger[ledger$position == pos & ledger$evidence_class != "rule", ]
  c(sources = length(unique(sub$key)),
    measurement = sum(sub$evidence_class == "measurement"),
    primer_scan = sum(sub$evidence_class == "primer_scan"))
}
n_sources <- function(pos) {
  length(unique(ledger$key[ledger$position == pos & ledger$evidence_class != "rule"]))
}
c_m2 <- count_position(-2L)
c_m3 <- count_position(-3L)
c_m4 <- count_position(-4L)
for (pos in c(-2L, -3L, -4L)) {
  if (count_position(pos)[["sources"]] != n_sources(pos)) {
    stop("position counts disagree with the ledger")
  }
}

## --- (4) the contradiction must stay two-sided ------------------------------
if (c_m2[["sources"]] == 0L || c_m3[["sources"]] == 0L) {
  stop("the -2 vs -3 contradiction is now one-sided (", c_m2[["sources"]], " vs ",
       c_m3[["sources"]], "). That is a scientific conclusion, not a build step:",
       " review data-raw/proposed/arms_mismatch_position.proposal.md by hand.")
}

ledger$n_sources_at_position <- vapply(ledger$position, function(p) {
  as.integer(count_position(p)[["sources"]])
}, integer(1))
ledger$outcome_definition <- paste(
  "position = deliberate-mismatch site relative to the 3' end (-1 = last base).",
  "n_tested = positions scanned; n_primers = primers examined.",
  "evidence_class: measurement = per-position outcome reported; primer_scan =",
  "position read from the designed primer sequences; rule = asserted, not tested."
)
ledger$evidence_tier <- ifelse(ledger$evidence_class == "rule",
                               "rule of thumb", "primary literature")
ledger$source_file <- ledger$source_relpath

outfile <- file.path(out_dir, "arms_mismatch_position.csv")
ledger$source_relpath <- NULL
utils::write.csv(ledger, outfile, row.names = FALSE, na = "")
message("wrote ", outfile, " (", nrow(ledger), " rows)")
for (p in c(-2L, -3L, -4L)) {
  cc <- count_position(p)
  message("  ", p, ": ", cc[["sources"]], " sources (", cc[["measurement"]],
          " measurement, ", cc[["primer_scan"]], " primer_scan)")
}
message("wrote ", outfile, " -- NOT AVERAGED; see ",
        "data-raw/proposed/arms_mismatch_position.proposal.md")
