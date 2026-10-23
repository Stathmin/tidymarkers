# RPA 3'-anchor mismatch scheme
#
# Decision: data-raw/proposed/rpa_anchor_rule.decision.md.  The weights RANK
# placements; they do not predict a call rate, and the paper's own position
# cells do not order the non-terminal offsets.

ext <- function(name) {
  utils::read.csv(
    tm_extdata(name),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
}

test_that("the terminal weights are the extract's twelve cells, unstarred ones included", {
  t <- ext("rpa_mismatch_terminal.csv")
  expect_equal(nrow(t), 12L)
  w <- rpa_anchor_weights()
  # the weights are held in the extract's row ORDER, so the order is asserted
  expect_equal(w$terminal_pairs, t$terminal_mismatch)
  expect_equal(w$terminal, t$RS_logit, tolerance = 1e-9)
  # the two unstarred cells are a FINDING, not missing data: RPA tolerates them
  unstarred <- t$terminal_mismatch[!t$RS_starred]
  expect_setequal(unstarred, c("A-C", "T-G"))
  expect_true(all(w$terminal[w$terminal_pairs %in% unstarred] > -1))
  expect_true(all(abs(t$RS_probability[!t$RS_starred] - 0.984) < 1e-9))
})

test_that("the anchor is 4 nt and the non-terminal offsets are NOT ordered", {
  expect_equal(RPA_ANCHOR_NT, 4L)
  w <- rpa_anchor_weights()
  expect_length(w$nonterminal, 3L)
  # 2 and 3 are indistinguishable on the paper's own cells, so the weights are
  # equal - a strict order would be invented
  expect_equal(w$nonterminal[[2L]], w$nonterminal[[3L]])
  expect_false(w$calibrated)
})

test_that("rpa_anchor_score reads a primer's 3' anchor", {
  # a perfect anchor costs nothing
  s <- rpa_anchor_score("ACGTACGTAC", "ACGTACGTAC")
  expect_equal(s$penalty, 0)
  expect_equal(s$mismatches, 0L)
  expect_equal(s$basis, "exact")
  # a single terminal mismatch is exactly its printed cell: primer T over
  # template C is the extract's `C-T` row (-3.40), and the extract names each
  # cell template-first
  s <- rpa_anchor_score("ACGTACGTAT", "ACGTACGTAC")
  expect_equal(s$mismatches, 1L)
  expect_equal(s$positions, 0L)
  expect_equal(s$penalty, -3.40)
  # and one of the paper's TOLERANT, unstarred cells (A-C) costs far less
  s <- rpa_anchor_score("ACGTACGTAC", "ACGTACGTAA")
  expect_equal(s$penalty, -0.66)
  # a mismatch OUTSIDE the anchor is not scored
  s <- rpa_anchor_score("ACGTACGTAC", "ACGTTCGTAC")
  expect_equal(s$mismatches, 0L)
  expect_equal(s$penalty, 0)
  # two inside the anchor: one class cell per offset, and no order is claimed
  # between the non-terminal offsets
  s <- rpa_anchor_score("ACGTACGTAC", "ACGTACGTTA")
  expect_equal(s$mismatches, 2L)
  expect_false(s$nonterminal_ordered)
})

test_that("the paper's own 279 cells do not order the non-terminal offsets", {
  # this is the decision's central finding: the design rule says 3 is worst, the
  # cells say 2 and 3 are indistinguishable.  If the cells ever DO order them,
  # the decision needs revisiting - so the test fails then.
  p <- ext("rpa_mismatch_position.csv")
  expect_equal(nrow(p), 279L)
  expect_setequal(unique(p$nonterminal_offset), 1:3)
  mean_by <- tapply(p$RS_logit, p$nonterminal_offset, mean)
  expect_lt(abs(mean_by[["2"]] - mean_by[["3"]]), 0.1)
  # and the spread within each offset is far larger than the gap between means
  expect_gt(diff(range(mean_by)), 0)
  expect_lt(diff(range(mean_by)), diff(range(p$RS_logit)) / 10)
})
