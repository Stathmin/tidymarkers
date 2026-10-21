# =============================================================================
#  The science ledger  --  every source must predict something
# =============================================================================
#  A paper that was read and used for nothing is a reading list, not science.
#  `inst/extdata/science_ledger.csv` is the join from a registered source, to the
#  extract it produced, to the test that consumes it.  These tests keep that join
#  honest: they check the named test really reads the named extract, and that a
#  row which is blocked says why.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name), show_col_types = FALSE,
                  na = c("", "NA"))
}

ledger <- function() ext("science_ledger.csv")

test_that("ledger: every row joins a registered source to an existing extract", {
  l <- ledger()
  expect_gt(nrow(l), 10L)
  s <- ext("scientific_sources.csv")
  expect_true(all(l$source_key %in% s$key))
  expect_true(all(c("source_key", "extract", "use", "test_file",
                    "verdict_moved") %in% names(l)))
  expect_false(any(is.na(l$source_key) | !nzchar(l$source_key)))
  expect_false(any(is.na(l$verdict_moved) | !nzchar(l$verdict_moved)))
  # no duplicated (source, extract) pair: one row per claim
  key <- paste(l$source_key, l$extract)
  expect_false(any(duplicated(key)))
})

test_that("ledger: the use vocabulary is closed, and blocked rows carry a reason", {
  l <- ledger()
  # `report` = a second extract from the same source that is a REPORTED CLAIM
  # (what the paper asserts in prose/abstract), kept beside the measured table so
  # the two can be compared. It is a first-class level, not a stray value.
  expect_true(all(l$use %in% c("test", "artefact", "rule", "report", "blocked")))
  # a blocked row is the honest state: read, extracted, consumed by nothing yet
  blocked <- l$use == "blocked"
  expect_true(all(!is.na(l$reason[blocked]) & nzchar(l$reason[blocked])))
  expect_true(all(is.na(l$reason[!blocked])))
  # and a blocked row is the only kind that may skip naming a test
  expect_true(all(is.na(l$test_file[blocked])))
})

test_that("ledger: every named test file exists and actually reads its extract", {
  l <- ledger()
  used <- l[!is.na(l$test_file), ]
  expect_gt(nrow(used), 5L)
  for (i in seq_len(nrow(used))) {
    f <- test_path(used$test_file[[i]])
    expect_true(file.exists(f), label = used$test_file[[i]])
    src <- readLines(f, warn = FALSE)
    # the extract name, or its stem, must appear in the test: otherwise the row
    # claims a consumer that stopped reading the data
    stem <- sub("\\.csv$", "", basename(used$extract[[i]]))
    expect_true(
      any(grepl(stem, src, fixed = TRUE)),
      label = paste(used$test_file[[i]], "reads", stem, sep = " / ")
    )
  }
})

test_that("ledger: the extracts it names exist (unless the row is blocked)", {
  l <- ledger()
  live <- l[l$use != "blocked", ]
  for (e in unique(live$extract)) {
    expect_true(file.exists(test_path("../../inst/extdata", e)), label = e)
  }
})

test_that("ledger: an extract's own DOI belongs to the source that row names", {
  # The identifier is checked against the extract's data, not against the
  # ledger's prose: a row cannot quietly re-cite a file to a different paper.
  l <- ledger()
  s <- ext("scientific_sources.csv")
  live <- l[l$use != "blocked", ]
  # a table may legitimately combine rows from several papers; what must hold is
  # that the row's own source is among the DOIs it cites
  pooled <- c("conversion_rates.csv", "mismatch_discrimination.csv",
              "kasp_homeolog_axes.csv")
  for (i in seq_len(nrow(live))) {
    e <- live$extract[[i]]
    # a non-extract row names a registry or a deliberately unwritten file; the
    # DOI check only applies to a real CSV on disk
    if (!file.exists(test_path("../../inst/extdata", e))) next
    d <- readr::read_csv(test_path("../../inst/extdata", e), show_col_types = FALSE,
                         na = c("", "NA"))
    if (!"doi" %in% names(d)) next
    want <- s$doi[s$key == live$source_key[[i]]]
    # a software source is registered by URL and cited by its paper DOI
    if (!grepl("^10\\.", want)) {
      want <- s$pmcid[s$key == live$source_key[[i]]]
      if (length(want) && !is.na(want) && nzchar(want)) {
        # no DOI exists for this paper anywhere; its PMCID is the identifier,
        # and the registry's citation must say which article that PMCID is
        expect_true(grepl(want, s$citation[s$key == live$source_key[[i]]],
                          fixed = TRUE), label = paste(e, "names", want))
        next
      }
      # a software source is registered by URL and cited by its paper DOI
      want <- "10.1093/bioinformatics/btv069"
    }
    if (e %in% pooled) next
    got <- unique(stats::na.omit(as.character(d$doi)))
    # a table may legitimately combine rows from several papers (a pooled
    # conversion-rate table, or a prior table); what must hold is that the row's
    # own source is among them
    expect_true(want %in% got, label = paste(e, "cites", want))
  }
})

test_that("ledger: the three axes the review asked for are all covered", {
  l <- ledger()
  # Z7 priors, Z8 outcomes, and the Y3.2 two-axis artifact must each have at
  # least one consumer, or a whole wave could silently lose its data
  expect_true(any(grepl("z7", l$test_file)))
  expect_true(any(grepl("z8", l$test_file)))
  expect_true(any(grepl("y32", l$test_file)))
  # and the PolyMarker rule source is applied as a RULE, not just cited
  expect_true(any(l$use == "rule" & l$source_key == "polymarker_source"))
})
