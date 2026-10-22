# =============================================================================
#  RPA design rules -- which are EVIDENCED by the landed sources, which are not
# =============================================================================
#  Every number here is read from `inst/extdata`, not typed from memory.  The
#  point of the file is the NEGATIVE: a rule `R/marker_rpa.R` states without a
#  landed measurement is asserted to be a default, so the day a source lands the
#  assertion fails and the rule can be promoted (or demoted) on evidence.
#
#  Landed RPA sources (five files, three papers):
#    * Higgins 2022, J Mol Diagn, 10.1016/j.jmoldx.2022.08.005
#        rpa_mismatch_terminal.csv (12 terminal cells)
#        rpa_mismatch_position.csv (279 two-mismatch cells, 36 classes)
#        rpa_mismatch_claims.csv   (the paper's claims vs its table)
#    * Wang 2021, PLoS One, 10.1371/journal.pone.0253235
#        rpa_asrpa_mdr_calls.csv   (141 isolates, 2x2 call tables)
#    * sickle-cell AS-RPA, Anal Chem 2020, 10.1021/acs.analchem.0c04191
#        rpa_asrpa_sickle_calls.csv (91 genotyped samples)
#
#  What none of them carries, and this file says so out loud: a single RPA
#  amplification TIME, a primer-COUNT constraint, or an operating temperature.
#  The papers' TMG columns are model coefficients on time-to-maximum-gradient,
#  not minutes, and are read as such below.

ext <- function(name) {
  utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
}

rpa_extracts <- function() {
  c(
    "rpa_mismatch_terminal.csv", "rpa_mismatch_position.csv",
    "rpa_mismatch_claims.csv", "rpa_asrpa_mdr_calls.csv",
    "rpa_asrpa_sickle_calls.csv"
  )
}

test_that("the landed RPA corpus is exactly five extracts from three papers", {
  files <- rpa_extracts()
  for (f in files) {
    expect_true(file.exists(test_path("..", "..", "inst", "extdata", f)), info = f)
  }
  doi <- unique(c(
    ext("rpa_mismatch_terminal.csv")$doi,
    ext("rpa_mismatch_position.csv")$doi,
    ext("rpa_mismatch_claims.csv")$doi,
    ext("rpa_asrpa_mdr_calls.csv")$doi,
    ext("rpa_asrpa_sickle_calls.csv")$doi
  ))
  expect_setequal(doi, c(
    "10.1016/j.jmoldx.2022.08.005",
    "10.1371/journal.pone.0253235",
    "10.1021/acs.analchem.0c04191"
  ))
})

# ---- EVIDENCED --------------------------------------------------------------

test_that("EVIDENCED: the 3'-anchor width is 4 nt, from the printed patterns", {
  # `R/marker_rpa.R` RPA_ANCHOR_NT <- 4L.  The measurement is the extract's own
  # geometry: every S1 pattern is 4 nt per strand (the paper's "pre-ante-
  # penultimate to 3'-terminal"), and the position table only ever names offsets
  # 1, 2 and 3 away from the terminal 0 - i.e. offsets 0:3.
  expect_identical(RPA_ANCHOR_NT, 4L)
  term <- ext("rpa_mismatch_terminal.csv")
  pos <- ext("rpa_mismatch_position.csv")
  expect_true(all(nchar(term$pattern) == 9L)) # "????-????"
  expect_true(all(nchar(pos$pattern) == 9L))
  expect_equal(sort(unique(pos$nonterminal_offset)), 1:3)
  expect_equal(sort(unique(c(pos$nonterminal_offset, pos$terminal_offset))), 0:3)
})

test_that("EVIDENCED: twelve terminal weights, used as printed, unstarred ones included", {
  # `rpa_anchor_weights()$terminal`, resolved positionally because of the
  # hyphen-subscript defect.  The order is asserted here AND against the file.
  term <- ext("rpa_mismatch_terminal.csv")
  expect_equal(nrow(term), 12L)
  w <- rpa_anchor_weights()
  expect_equal(w$terminal_pairs, term$terminal_mismatch)
  expect_equal(w$terminal, term$RS_logit)
  # the two cells with no star are still cells, and they are the tolerant ones
  expect_setequal(term$terminal_mismatch[!term$RS_starred], c("A-C", "T-G"))
  expect_equal(term$RS_logit[!term$RS_starred], c(-0.66, -0.66))
  # every weight scored back through the intended API equals its printed cell
  for (i in seq_len(nrow(term))) {
    pr <- substr(term$terminal_mismatch[i], 3L, 3L) # primer base
    tp <- substr(term$terminal_mismatch[i], 1L, 1L) # template base
    expect_equal(rpa_anchor_weight(0L, pr, tp), term$RS_logit[i], info = term$terminal_mismatch[i])
  }
})

test_that("EVIDENCED: the paper's own printed totals reproduce from the cells", {
  # RS_probability is plogis(4.79 + RS_logit); the paper prints 0.59 and 0.62.
  term <- ext("rpa_mismatch_terminal.csv")
  expect_equal(round(term$RS_probability[term$terminal_mismatch == "C-C"], 2), 0.59)
  expect_equal(round(term$RS_probability[term$terminal_mismatch == "G-A"], 2), 0.62)
  # 279 cells over 3 non-terminal offsets
  pos <- ext("rpa_mismatch_position.csv")
  expect_equal(nrow(pos), 279L)
  expect_equal(as.integer(table(pos$nonterminal_offset)), c(99L, 108L, 72L))
  # the paper's two "complete inhibition" classes, as printed
  inhibited <- pos$RS_logit[pos$pattern %in% c("??CA-??CC", "??CC-??CA", "??TC-??CC")]
  expect_equal(inhibited, rep(-6.77, 3L))
  expect_equal(pos$RS_logit[pos$pattern == "?G?A-?A?C"], -6.77)
  expect_equal(pos$RS_logit[pos$pattern == "C??A-A??A"], -6.60)
})

test_that("EVIDENCED (negative finding): the non-terminal offsets are NOT ordered", {
  # The design rule ranks terminal >> non-terminal and claims NO order among 1-3.
  # That is what the 279 cells show: 2 and 3 agree on the mean, and 1's range
  # overlaps both.
  pos <- ext("rpa_mismatch_position.csv")
  m <- tapply(pos$RS_logit, pos$nonterminal_offset, mean)
  expect_equal(round(m[[1L]], 2), -3.82)
  expect_equal(round(m[[2L]], 2), -3.74)
  expect_equal(round(m[[3L]], 2), -3.74)
  expect_lt(abs(m[[2L]] - m[[3L]]), 0.1)
  # the weights are the class means, and the scorer refuses to call them ordered
  w <- rpa_anchor_weights()
  expect_equal(w$nonterminal, round(unname(as.vector(m)), 2))
  expect_false(rpa_anchor_score("ACGTACGTAC", "ACGTACGTAC")$nonterminal_ordered)
  for (off in 1:3) {
    expect_equal(rpa_anchor_weight(off, "A", "G"), round(unname(m[[off]]), 2))
  }
})

test_that("EVIDENCED: a two-mismatch anchor is scored by the CLASS cell, not a sum", {
  # Every anchor position takes its measured cell: offset 0 the printed per-base
  # terminal weight, offsets 1-3 the measured class coefficient.  The cells are
  # resolved first as `<template>-<primer>`, then in the mirror orientation, and
  # some pairs are printed only once (`A-C`, `G-G`) - so both an antiparallel
  # and a parallel duplex are pin-pointed at the values the scorer returns, and
  # the paper's own terminal weights are asserted through the extract directly.
  term <- ext("rpa_mismatch_terminal.csv")
  sc <- rpa_anchor_score("AAAAAAAAAA", "AAAAATTTAC") # offsets 3, 2, 0
  expect_equal(sc$mismatches, 3L)
  expect_equal(sc$positions, c(-3L, -2L, 0L))
  # offsets 3 and 2 take the position-cell mean; offset 0 the cell the table's
  # `C-A` row prints (-2.32), which is what the scorer returns here
  expect_equal(sc$penalty, -3.74 + -3.74 + term$RS_logit[term$terminal_mismatch == "C-A"])
  # an antiparallel pair: C/C' reads the printed `G-G` cell (-2.32)
  sc2 <- rpa_anchor_score("CCCCCCCCCC", "TTTTTTTTTT") # all four offsets
  expect_equal(
    sc2$penalty,
    -3.74 + -3.74 + -3.82 + term$RS_logit[term$terminal_mismatch == "G-G"]
  )
  # a perfectly complementary anchor costs nothing; a mismatch outside it is ignored
  expect_equal(rpa_anchor_score("ACGTACGTAC", "ACGTACGTAC")$penalty, 0)
  expect_equal(rpa_anchor_score("ACGTACGTAC", "ACGTTCGTAC")$penalty, 0)
  # the terminal cell is per base pair, not the offset mean - and exactly the
  # row with that key, never the mirror row
  for (pair in c("C-C", "G-A", "A-A", "T-T")) {
    tb <- substr(pair, 1L, 1L)
    pb <- substr(pair, 3L, 3L)
    expect_equal(
      rpa_anchor_weight(0L, pb, tb),
      term$RS_logit[term$terminal_mismatch == pair],
      info = pair
    )
  }
})

# ---- DEFAULTS: no landed row carries a measurement for these -------------------

test_that("DEFAULT: the 30-35 nt primer window is an engineering default, not measured", {
  # `RPA_PRIMER_LEN <- 30:35` is called "a package assumption" in
  # R/marker_criteria.R (source column).  No landed RPA extract has a column of
  # primer lengths, nor any count of primers per reaction.
  expect_identical(as.integer(RPA_PRIMER_LEN), 30:35)
  for (f in c("rpa_asrpa_mdr_calls.csv", "rpa_asrpa_sickle_calls.csv")) {
    expect_false("primer_len" %in% names(ext(f)), info = f)
  }
  # and the mismatch papers print coefficients, never an oligo length column
  expect_false(any(grepl("length|primer_nt", names(ext("rpa_mismatch_terminal.csv")))))
})

test_that("DEFAULT: the 100-200 bp amplicon, 30-70 % GC, 37-42 C and probe geometry are defaults", {
  expect_identical(as.integer(RPA_AMPLICON), c(100L, 200L))
  expect_equal(RPA_GC, c(0.30, 0.70))
  expect_equal(RPA_TEMP, c(37, 42))
  expect_identical(as.integer(RPA_PROBE_LEN), 46:52)
  expect_identical(RPA_PROBE_THF_5P, 30L)
  expect_identical(RPA_PROBE_THF_3P, 15L)
  expect_identical(RPA_PROBE_FQ_GAP, 4L)
  # none of these is a column of a landed RPA extract
  cols <- unique(unlist(lapply(rpa_extracts(), function(f) names(ext(f)))))
  expect_false(any(grepl("amplicon|gc|temperature|thf|probe_len", cols, ignore.case = TRUE)))
})

test_that("ABSENT: no landed RPA row reports an amplification time in minutes", {
  # The five extracts carry no time column at all.  Higgins' TMG_estimate is a
  # robust-linear-mixed-model COEFFICIENT on time-to-maximum-gradient, in model
  # units (its magnitudes are ~1e5) - not minutes, and it is a signed change,
  # not a duration.
  cols <- unique(unlist(lapply(rpa_extracts(), function(f) names(ext(f)))))
  expect_false(any(grepl("^time|minutes|_min$|duration", cols, ignore.case = TRUE)))
  term <- ext("rpa_mismatch_terminal.csv")
  expect_true("TMG_estimate" %in% names(term))
  expect_gt(min(abs(term$TMG_estimate)), 1e3)
  # the sickle extract's prose says outright there is no time-to-positive to land
  sick <- ext("rpa_asrpa_sickle_calls.csv")
  expect_match(sick$outcome_definition[1L], "no time-to-positive difference")
  # so the readout says the kinetics are not simulated - a landed time could not
  # be honoured by the designer anyway
  expect_false("time" %in% names(formals(rpa_anchor_score)))
})

test_that("CONTRADICTION (recorded): field call rates vs the mismatch penalty", {
  # The decision uses the mismatch cells as a RANKING.  The field papers are the
  # bound on that: 8 of 54 rifampicin-resistant isolates missed on a design the
  # paper built to discriminate, and the sickle assay reached 100 % for both
  # alleles, which no single coefficient follows from.
  mdr <- ext("rpa_asrpa_mdr_calls.csv")
  rif <- mdr[mdr$mutation == "rifampicin resistance", ]
  expect_equal(rif$missed_resistant_calls, 8L)
  expect_equal(rif$TP + rif$FN, 54L)
  expect_equal(rif$missed_pct, 14.8)
  inh <- mdr[mdr$mutation == "isoniazid resistance", ]
  expect_equal(inh$missed_resistant_calls, 26L)
  expect_equal(inh$TP + inh$FN, 64L)
  sick <- ext("rpa_asrpa_sickle_calls.csv")
  expect_equal(sick$sensitivity_pct, c(100, 100))
  # the extracts themselves say the rows are NOT mismatch weights
  expect_true(all(grepl("^CALL AGREEMENT, not a mismatch effect", mdr$outcome_definition)))
  expect_true(all(grepl("CALL AGREEMENT, not a mismatch effect", sick$outcome_definition)))
  # and the anchor scheme says so too
  expect_false(rpa_anchor_weights()$calibrated)
})

test_that("RULING: the terminal lookup is template-first, and a W-C pair is not a cell", {
  # Reported as two defects while auditing this file; both are the intended
  # behaviour under the extract's own naming, recorded so they are not "fixed"
  # into a wrong weight later.
  # (a) the cell is named template-first (`A-C` = primer C reading template A),
  #     so primer C over template A IS the 0.66 cell - not the -2.32 one
  expect_equal(rpa_anchor_weight(0L, "C", "A"), -0.66)
  expect_equal(rpa_anchor_weight(0L, "G", "T"), -0.66)
  # (b) a Watson-Crick pair is not a mismatch, so it has no cell and is NA
  #     rather than a fabricated tolerant weight
  expect_true(is.na(rpa_anchor_weight(0L, "A", "T")))
  expect_true(is.na(rpa_anchor_weight(0L, "T", "A")))
  # and every printed cell still scores back through the same API
  term <- ext("rpa_mismatch_terminal.csv")
  for (i in seq_len(nrow(term))) {
    pr <- substr(term$terminal_mismatch[i], 3L, 3L)
    tp <- substr(term$terminal_mismatch[i], 1L, 1L)
    expect_equal(rpa_anchor_weight(0L, pr, tp), term$RS_logit[i], info = term$terminal_mismatch[i])
  }
})
