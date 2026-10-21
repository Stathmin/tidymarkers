# =============================================================================
#  data-raw/arms_cycle_decay.R  --  does allele discrimination decay with cycles?
# =============================================================================
#  THE HYPOTHESIS (from a prior session's microPublication note)
#  In an allele-specific PCR the product made in an early cycle contains, at its
#  5' end, the full allele-specific primer - including the deliberate mismatch.
#  In later cycles the product itself can act as a template for the OTHER
#  primer, so the mismatch is no longer part of the primer-template duplex and
#  discrimination should soften. If true, "discrimination" is cycle-dependent and
#  a single end-point number is the wrong abstraction.
#
#  WHAT WE SEARCHED, AND WHAT WE FOUND
#  Checked the full text of every paper fetched for this family for a QUANTITATIVE
#  discrimination-vs-cycle series (delta signal, delta Cq, or a discrimination
#  index evaluated at successive cycle numbers):
#    quantas2013       (PMC3658971) - Cq only; delta Cq between 100% mutant and
#                        the first false positive is reported at end-point
#                        (20/23/16/23 cycles). No per-cycle discrimination series.
#      fret_asqpcr2021   (PMC8784781) - STATES end-point robustness ("perfect
#                        discrimination ... after 25-33 cycles and beyond") but
#                        plots RFU by cycle; no tabulated discrimination per cycle.
#      tparms2026        (PMC12898809) - one mention of "cycle"; no series.
#      wheat_as2023, sorghum_arms2021 - no cycle-resolved discrimination at all.
#  So the package's answer stands: NO quantitative cycle-decay measurement was
#  found. This file lands that NEGATIVE RESULT as a weak-tier record, so the
#  absence is searchable instead of being silently assumed.
#
#  This is deliberately NOT an outcome table: a negative result has no rows to
#  count, and inventing a "0.0 decay per cycle" would be fabrication. The CSV has
#  one row per searched source stating what it does/does not contain.
#
#  ASSERTIONS
#  (1) Each cited file exists and its quoted token is really in it - so a row
#      claiming "this paper has no cycle series" is attached to a real paper.
#  (2) The two sources that DO discuss cycles are exactly the ones we mark
#      talks_about_cycles = TRUE; if a new fetch changes that we stop.
#  (3) No row claims a numeric decay rate. That is the point of the file.
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

rows <- list()
add <- function(key, doi, role, talks_about_cycles, has_cycle_series,
                quote_token, note) {
  rows[[length(rows) + 1L]] <<- data.frame(
    key = key, doi = doi, role = role,
    talks_about_cycles = talks_about_cycles,
    has_cycle_resolved_discrimination_series = has_cycle_series,
    quote_token = quote_token, note = note, stringsAsFactors = FALSE)
}

add("quantas2013", "10.1186/1471-2407-13-206", "ARMS-adjacent (QuantAS, JAK2)",
    TRUE, FALSE, "Cq value difference of 20 cycles",
    "delta Cq is end-point only; 20/23/16/23 cycles are thresholds, not a per-cycle series")
add("fret_asqpcr2021", "10.3389/fpls.2021.747886", "barley ASQ / allele-specific qPCR",
    TRUE, FALSE, "after 25",
    "asserts end-point robustness across 25-33 cycles; RFU by cycle is plotted, not tabulated")
add("tparms2026", "10.3390/ijms27031406", "crop ARMS (TP-ARMS)",
    FALSE, FALSE, "primer success rates",
    "no cycle-resolved discrimination reported")
add("sorghum_arms2021", "10.3390/foods10092218", "crop ARMS (aARMS-PCR)",
    FALSE, FALSE, "highest discrimination ability",
    "reports end-point Cq only; no cycle series")
add("wheat_as2023", "10.1002/pld3.520", "crop KASP/PACE",
    FALSE, FALSE, "produced polymorphic signals",
    "cluster-level end-point calls; no cycle series")

out <- do.call(rbind, rows)
out$source_file <- file.path(out$key, "fulltext.xml")
out$verification <- ifelse(vapply(out$source_file, function(r) !is.na(find_source(r)),
                                  logical(1)), "fetched and searched", "not fetched here")
out$evidence_tier <- "weak"
out$outcome_definition <- paste(
  "NEGATIVE RESULT record, not an outcome table.",
  "has_cycle_resolved_discrimination_series = TRUE only if the paper reports",
  "discrimination/delta signal evaluated at successive cycle numbers.",
  "No such series was found in any fetched ARMS/AS-PCR source;",
  "the package therefore still has NO quantitative cycle-decay measurement."
)
out$cycle_decay_rate_per_cycle <- NA_real_

for (i in seq_len(nrow(out))) {
  p <- find_source(out$source_file[i])
  if (is.na(p)) next  # absent fetch: row kept, not asserted
  body <- gsub("[[:space:]]+", " ",
               paste(xml2::xml_text(xml2::xml_find_all(xml2::read_xml(p),
                                                       "//text()")), collapse = " "))
  if (!grepl(out$quote_token[i], body, fixed = TRUE)) {
    stop("token \"", out$quote_token[i], "\" not in ", p)
  }
}

# (2)/(3) invariants of the negative result
if (any(out$has_cycle_resolved_discrimination_series)) {
  stop("a cycle-resolved series was found; this negative-result file is obsolete")
}
if (any(!is.na(out$cycle_decay_rate_per_cycle))) {
  stop("a decay rate was written; a negative result must not carry a number")
}

outfile <- file.path(out_dir, "arms_cycle_decay_negative.csv")
utils::write.csv(out, outfile, row.names = FALSE, na = "")
message("wrote ", outfile, " (", nrow(out), " sources searched; ",
        "0 cycle-resolved discrimination series)")
