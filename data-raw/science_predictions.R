# =============================================================================
#  data-raw/science_predictions.R  --  every source must predict something
# =============================================================================
#  A registry of papers that nothing uses is a reading list, not science.  This
#  table is the join between the literature registry (scientific_sources.csv),
#  the extracted data (inst/extdata/*.csv), and the tests that consume it:
#
#      source  ->  extract file  ->  test id  ->  verdict column
#
#  Every row must (a) name a source that is in the registry, (b) name an extract
#  that exists, (c) name a test file that exists, and (d) say which column of
#  which evaluation it would move.  tests/testthat/test-science-ledger.R checks
#  that each named test actually mentions the extract it is paired with, so a
#  row cannot claim a consumer that has since stopped reading the data.
#
#  `use` is one of
#    test      the extract is an oracle an existing test asserts against
#    artefact  the extract is shown to the user (a design artifact column)
#    rule      the extract is a design rule the constructor may apply
#    blocked   read and extracted, but nothing consumes it yet -- with the reason
#  A `blocked` row is honest, and must carry `reason`.
#
#  Run:  Rscript data-raw/science_predictions.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
out_dir <- fs::path(root, "inst", "extdata")

sources <- utils::read.csv(fs::path(out_dir, "scientific_sources.csv"), stringsAsFactors = FALSE)
# keys whose registry identifier is a URL rather than a DOI
source_doi_fallback <- "10.1093/bioinformatics/btv069"
# The extract scripts' own DOI constants are the identifier of record for a
# source row, so the ledger can name a source whose registry row lives on an
# agent branch without guessing a DOI here.
source_doi <- function(key) {
  # a software source is registered with its repository as the identifier and
  # the paper DOI in a second column; the DOI is the one extracts cite
  if (key %in% sources$key) {
    d <- sources$doi[sources$key == key][[1L]]
    if (!grepl("^10\\.", d)) d <- source_doi_fallback
    return(d)
  }
  # a script may be named for the backlog key (liu2012_arms) or for the paper
  # (liu2012); accept either, but require exactly one
  cand <- fs::path(root, "data-raw", paste0(c(key, sub("_arms$", "", key)), "-supplement.R"))
  cand <- cand[file.exists(cand)]
  if (length(cand) != 1L) {
    stop("no registry row and no single extract script for ", key)
  }
  # the constant is at the top of every extract script; the grep is bounded to
  # the first 80 lines so a DOI inside a quoted paper sentence is not picked up
  txt <- base::readLines(cand, warn = FALSE, n = 80L)
  hit <- grep('^DOI *<- *"', txt, value = TRUE)
  if (length(hit) != 1L) {
    stop("extract script for ", key, " has ", length(hit), " DOI constants near the top")
  }
  sub('^DOI *<- *"([^"]+)".*', "\\1", hit)
}

# one row per (source, extract, consumer).  Rows are written one at a time
# so counts cannot drift out of alignment, which is what bit-rotted the
# column form of this table three times.
row <- function(source_key, extract, use, test_file, verdict_moved, reason = NA_character_) {
  data.frame(source_key = source_key, extract = extract, use = use,
             test_file = test_file, verdict_moved = verdict_moved,
             reason = reason, stringsAsFactors = FALSE)
}
# test_file is a fact about the test THAT exists today: the Z7/Z8/Y3.2 tests
# land in the same wave as this ledger entry.
ledger <- do.call(rbind, list(
  row("kwok1990", "mismatch_discrimination.csv", "test", "test-z7-priors.R", "3'-terminal severity class", NA_character_),
  row("kwok1990", "z7_priors.csv", "test", "test-z7-priors.R", "kwok_severity_order", NA_character_),
  row("lefever2019", "lefever_single_mismatch.csv", "test", "test-z7-priors.R", "offset of the artificial mismatch", NA_character_),
  row("lefever2019", "z7_priors.csv", "test", "test-z7-priors.R", "lefever_position", NA_character_),
  row("liu2012_arms", "liu2012_outcomes.csv", "test", "test-z7-priors.R", "mismatch site (2nd/3rd/4th)", NA_character_),
  row("ye2001", "mismatch_discrimination.csv", "test", "test-z7-priors.R", "ye_second_mismatch_position", NA_character_),
  row("makhoul2020", "makhoul_design_strategies.csv", "test", "test-z8-outcomes.R", "design_outcome(assay): KASP usable / not usable", NA_character_),
  row("makhoul2020", "makhoul_kasp_primers.csv", "artefact", "test-y32-homeolog.R", "artifact$common_3p (subgenome axis)", NA_character_),
  row("grewal2020", "grewal_kasp_axes.csv", "test", "test-y32-homeolog.R", "artifact$as_fam_3p / artifact$as_hex_3p (allele axis)", NA_character_),
  row("grewal2020", "kasp_homeolog_axes.csv", "test", "test-y32-homeolog.R", "the two axes are separate columns, never one", NA_character_),
  row("grewal2020", "conversion_rates.csv", "test", "test-z8-outcomes.R", "conversion_rates: gate-pass -> usable counts", NA_character_),
  row("ramirezgonzalez2015", "conversion_rates.csv", "test", "test-z8-outcomes.R", "conversion_rates: outcome_scale stays ordinal", NA_character_),
  row("ramirezgonzalez2014", "conversion_rates.csv", "test", "test-z8-outcomes.R", "conversion_rates: outcome_scale stays ordinal", NA_character_),
  row("grewal2022_amblyopyrum", "conversion_rates.csv", "test", "test-z8-outcomes.R", "conversion_rates: failure classes are not pooled", NA_character_),
  row("polymarker_source", "polymarker_design_rules.csv", "rule", "test-y32-homeolog.R", "constraint: the common primer 3' base is the subgenome anchor", NA_character_),
  row("polymarker_source", "polymarker_design_rules.csv", "artefact", "test-y32-homeolog.R", "artifact$subgenome_anchor (which primer carries the genome)", NA_character_),
  row("newton1989", "mismatch_discrimination.csv", "test", "test-z7-priors.R", "mechanism prior only: a mismatched 3' base is not a primer", NA_character_),
  row("ye2001_primer_set", "ye2001_primer_set.csv (deliberately not written)", "blocked", NA_character_, "none yet (provenance only)", "the full text was read but no algebraic assertion is possible: the primer set carries no per-primer outcome to reproduce"),
  row("desilva2009", "desilva2009_outcomes.csv", "test", "test-z8-outcomes.R", "hrm: outcome counts with their own denominators", NA_character_),
  row("ace1r_aslamp2015", "ace1r_aslamp2015_counts.csv", "test", "test-z8-outcomes.R", "lamp: time-to-positive window, not a binary call", NA_character_),
  row("liu2016_caps", "liu2016_caps_funnel.csv", "test", "test-z8-outcomes.R", "caps: restriction-site conversion rate", NA_character_),
  row("soybean_est_caps", "soybean_est_caps_funnel.csv", "test", "test-z8-outcomes.R", "caps: restriction-site conversion rate", NA_character_),
  row("snap2caps", "snap2caps_markers.csv", "test", "test-z8-outcomes.R", "caps: dCAPS site creation", NA_character_),
  row("jak2v617f_ddpcr", "jak2v617f_ddpcr_limits.csv", "test", "test-z8-outcomes.R", "dpcr: LoB / LoD of a partition readout", NA_character_)
))

## --- checks -------------------------------------------------------------------
# Every registry row finishes with its own `DOI <- "..."` constant, so the
# identifier of record is the script's, not this file's.
for (i in seq_len(nrow(ledger))) {
  ext <- ledger$extract[[i]]
  if (ledger$use[[i]] == "blocked") next
  path <- fs::path(out_dir, ext)
  if (!file.exists(path)) stop("ledger names an extract that does not exist: ", ext)
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("doi" %in% names(d)) {
    want <- source_doi(ledger$source_key[[i]])
    got <- trimws(as.character(stats::na.omit(d$doi)))
    # some agents write an unquoted DOI such as 10.1055/s-0044-1785537, which a
    # CSV parser reads as a number; a file may cite several sources (its rows
    # say which), but the row's own source must be among them
    if (!length(got)) {
      row <- base::readLines(path, n = 2L)[[2L]]
      stop(ext, " has a doi column, but row 1 parses without one: ", substr(row, 1L, 120L))
    }
    if (!want %in% got && !(ext %in% c("conversion_rates.csv", "kasp_homeolog_axes.csv", "mismatch_discrimination.csv"))) {
      stop(ext, " does not cite ", want, " (it cites: ",
           paste(unique(got), collapse = ", "), ")")
    }
  }
}

bad_use <- setdiff(ledger$use, c("test", "artefact", "rule", "blocked"))
if (length(bad_use)) stop("unknown use: ", paste(bad_use, collapse = ", "))

# a test/artefact/rule row must name a test file, and that file must exist
needs_test <- ledger$use != "blocked"
if (any(needs_test & is.na(ledger$test_file))) {
  stop("rows with use != blocked must name a test_file: ",
       paste(ledger$extract[needs_test & is.na(ledger$test_file)], collapse = ", "))
}
for (f in unique(stats::na.omit(ledger$test_file))) {
  if (!file.exists(fs::path(root, "tests", "testthat", f))) {
    stop("ledger names a missing test file: ", f)
  }
}
# every blocked row carries its reason, and every non-blocked row does not
if (any(ledger$use == "blocked" & is.na(ledger$reason))) stop("a blocked row has no reason")
if (any(ledger$use != "blocked" & !is.na(ledger$reason))) stop("a non-blocked row carries a reason")

utils::write.csv(ledger, fs::path(out_dir, "science_ledger.csv"), row.names = FALSE, na = "")

message(
  "science ledger: ", nrow(ledger), " source->extract->test rows over ",
  length(unique(ledger$source_key)), " sources and ",
  length(unique(ledger$extract)), " extracts; ",
  sum(ledger$use == "test"), " tested, ", sum(ledger$use == "artefact"), " shown, ",
  sum(ledger$use == "rule"), " applied as a rule, ", sum(ledger$use == "blocked"), " blocked"
)
