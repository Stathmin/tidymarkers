# =============================================================================
#  probes/rice2018  --  the KASP validation RATE, and what a failure is
# =============================================================================
#  Rice 2018 (Steele et al., Mol Breed 38:38, doi:10.1007/s11032-018-0777-2)
#  is the first number in the Z7 family that is a rate rather than a ranking:
#  how often a KASP assay built on an allele-specific primer actually genotypes
#  segregating progeny. Table S5 of the supplement prints the summary, Table S6
#  lists every marker x cross combination it is counted from.
#
#  These tests guard the transcription and the arithmetic. They do NOT claim
#  the rate transfers to another crop, another polymerase or synthetic
#  templates: the paper's own outcome definition (in the CSV) says the failures
#  mix primer non-discrimination, locus non-polymorphism and readout failure.

ext_probe <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name), show_col_types = FALSE,
                  na = c("", "NA"), progress = FALSE)
}

test_that("probes/rice2018: Table S5 is 3 levels x 2 strata x 3 marker types", {
  o <- ext_probe("rice2018_kasp_validation.csv")
  expect_equal(nrow(o), 18L)
  expect_setequal(unique(o$level), c("MARKER--LEVEL", "CROSS--LEVEL", "ASSAY--LEVEL"))
  expect_setequal(unique(o$filter_criteria), c("met", "not_met"))
  expect_setequal(unique(o$marker_type), c("Existing", "Novel", "Total"))
  expect_equal(nrow(unique(o[, c("level", "filter_criteria", "marker_type")])), 18L)
  expect_true(all(o$n_tested == o$n_validated + o$n_not_validated))
  expect_true(all(o$n_validated >= 0L & o$n_tested >= o$n_validated))
  expect_true(all(o$evidence_tier == "measured"))
})

test_that("probes/rice2018: Total really is Existing + Novel in every stratum", {
  o <- ext_probe("rice2018_kasp_validation.csv")
  for (lv in unique(o$level)) for (fc in unique(o$filter_criteria)) {
    s <- o[o$level == lv & o$filter_criteria == fc, ]
    stopifnot(nrow(s) == 3L)
    tot <- s[s$marker_type == "Total", ]
    ex <- s[s$marker_type == "Existing", ]
    nov <- s[s$marker_type == "Novel", ]
    expect_equal(tot$n_validated, ex$n_validated + nov$n_validated,
                 info = paste(lv, fc))
    expect_equal(tot$n_not_validated, ex$n_not_validated + nov$n_not_validated,
                 info = paste(lv, fc))
  }
})

test_that("probes/rice2018: the printed percentages are the printed counts", {
  o <- ext_probe("rice2018_kasp_validation.csv")
  # one printed cell (assay-level Existing/not_met) is a rounding artefact of
  # the authors' unrounded counts: 86.62 recomputes, 86.63 is printed. The
  # generator documents this; the tolerance is one rounding step, never more.
  expect_true(all(abs(o$validated_pct_recomputed - o$validated_pct_printed) <= 0.0101))
  expect_lt(sum(abs(o$validated_pct_recomputed - o$validated_pct_printed) > 0.005), 2L)
  expect_true(all(o$validated_pct_printed >= 0 & o$validated_pct_printed <= 100))
})

test_that("probes/rice2018: the cross-level rate is the one the S6 rows give", {
  s6 <- ext_probe("rice2018_marker_cross.csv")
  expect_equal(nrow(s6), 647L)
  # 46 new designs with ids novel_01..novel_46, plus 75 existing
  expect_equal(sum(substr(s6$assay_id, 1L, 6L) == "novel_"), 122L)
  nov_ids <- unique(s6$assay_id[substr(s6$assay_id, 1L, 6L) == "novel_"])
  expect_equal(length(nov_ids), 46L)
  ex_ids <- setdiff(unique(s6$assay_id), nov_ids)
  expect_equal(length(ex_ids), 75L)
  # the prose totals: 412 met combinations, 394 validated
  expect_equal(sum(s6$met_criteria), 412L)
  expect_equal(sum(s6$validated & s6$met_criteria), 394L)
  # the paper's own arithmetic slips by 3 on the not-met side
  expect_equal(sum(!s6$met_criteria), 235L)
  expect_equal(sum(s6$validated & !s6$met_criteria), 201L)
})

test_that("probes/rice2018: a failed assay never produced a result column", {
  s6 <- ext_probe("rice2018_marker_cross.csv")
  # 'validated' is the paper's marker-cross call, which is independent of the
  # success/failure counts (a cross with 34 failures can still be validated),
  # so the only invariant is that both count columns are non-negative and that
  # a validated row is not one with zero successful results.
  expect_true(all(s6$n_successful_assays >= 0L & s6$n_failed_assays >= 0L))
  v <- s6[s6$validated, ]
  expect_true(all(v$n_successful_assays > 0L))
})

test_that("probes/rice2018: the extract still matches the fetched supplement", {
  lit <- test_path("../../.tmp/lit/rice2018")
  supp <- file.path(lit, "supplements", "11032_2018_777_MOESM3_ESM.pdf")
  skip_if_not(file.exists(supp), "rice2018 supplement not fetched")
  skip_if_not(nzchar(Sys.which("pdftotext")), "pdftotext not available")
  txt <- tempfile(fileext = ".txt")
  system2("pdftotext", c("-layout", supp, txt), stdout = FALSE)
  lay <- readLines(txt, warn = FALSE)
  unlink(txt)
  # the six count rows are printed verbatim; a hand-edit of the CSV would break
  # at least one of them
  printed <- c("48", "30", "5", "78", "335", "394", "18", "158", "201", "31",
               "7124", "12637", "1610", "3465", "5177", "1036")
  joined <- paste(lay, collapse = "\n")
  expect_true(all(vapply(printed, function(p) {
    grepl(paste0("criteria[[:space:]]+(?:[0-9.]+[[:space:]]+)*", p, "([[:space:]]|$)"),
          joined, perl = TRUE)
  }, logical(1))))
  o <- ext_probe("rice2018_kasp_validation.csv")
  expect_equal(o$n_validated[o$level == "CROSS--LEVEL" & o$filter_criteria == "met" &
                               o$marker_type == "Total"], 394L)
})
