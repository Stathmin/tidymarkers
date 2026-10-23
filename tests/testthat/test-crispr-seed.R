# Cas12a seed: a POSITION gate, with the measured magnitude kept as evidence
#
# Decision: data-raw/proposed/crispr_seed_rule.decision.md.  The point of these
# tests is that the gate never turns "inside the seed" into a magnitude claim,
# and that the evidence claims are the extract's own numbers.

ext <- function(name) {
  utils::read.csv(
    tm_extdata(name),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
}

test_that("the seed gate is a position rule for both effectors", {
  # Cas12a: the variant's offset in the spacer must fall in the PAM-proximal seed
  expect_true(crispr_seed_gate("cas12a", 20L, 0L))
  expect_true(crispr_seed_gate("cas12a", 20L, 7L))
  expect_false(crispr_seed_gate("cas12a", 20L, 8L))
  expect_false(crispr_seed_gate("cas12a", 23L, 12L))
  # Cas13: the variant must fall in the central critical region
  expect_true(crispr_seed_gate("cas13a", 28L, 14L))
  expect_false(crispr_seed_gate("cas13a", 28L, 0L))
  # an unknown offset is not a pass
  expect_false(crispr_seed_gate("cas12a", 20L, NA_integer_))
})

test_that("the gate alone says nothing about the magnitude, and the evidence says why", {
  ev <- crispr_seed_evidence()
  expect_s3_class(ev, "marker_crispr_seed_evidence")
  expect_match(ev$gate, "position only", fixed = TRUE)
  expect_equal(ev$magnitude, "unknown - no model predicts the ratio for a target")
  expect_false(ev$calibrated)
  expect_equal(ev$source, "10.1038/s42003-024-07173-7")
})

test_that("the measured seed evidence is the extract's own numbers", {
  d <- ext("communbiol2024_crrna_ratio.csv")
  expect_true(all(c("crrna", "n_mismatch", "mismatch_position", "mean_fc_vs_wt", "adjacency") %in% names(d)))

  # clause 1: an unmodified crRNA does not discriminate at any position
  wt <- d[d$crrna == "WT crRNA" & d$n_mismatch == 1L, ]
  expect_equal(nrow(wt), 20L)
  expect_equal(range(wt$mean_fc_vs_wt), c(0.766701, 1.393679), tolerance = 1e-5)

  # clause 2: a single seed mismatch on the optimised crRNA keeps most of the
  # signal - the range is the decision's, not a round number of ours
  fit <- d[d$crrna == "M39" & d$n_mismatch == 1L & d$mismatch_position %in% 1:8, ]
  expect_equal(nrow(fit), 8L)
  expect_equal(range(fit$mean_fc_vs_wt), c(0.3728637, 0.9679837), tolerance = 1e-6)
  ev <- crispr_seed_evidence()
  expect_equal(unname(ev$ratio_range), range(fit$mean_fc_vs_wt), tolerance = 1e-6)
  expect_equal(unname(ev$unmodified_crrna_range), range(wt$mean_fc_vs_wt), tolerance = 1e-5)

  # clause 3: adjacent doubles are usually but not always stronger - a TREND
  adj <- d[d$crrna == "M39" & d$n_mismatch == 2L & d$adjacency == "adjacent", ]
  singles <- d[d$crrna == "M39" & d$n_mismatch == 1L, ]
  expect_equal(nrow(adj), 19L)
  expect_equal(nrow(singles), 20L)
  expect_equal(sum(adj$mean_fc_vs_wt < 0.5), 17L)
  expect_equal(sum(singles$mean_fc_vs_wt < 0.5), 8L)
  # and it is not a rule: at least one adjacent pair is NOT better than both of
  # its constituent singles, which is why it is not encoded as a filter
  expect_lt(sum(adj$mean_fc_vs_wt < 0.5), nrow(adj))
})
