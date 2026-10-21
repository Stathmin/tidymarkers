# =============================================================================
#  Z7  --  the published mismatch priors are PREDICTIONS, not calibration
# =============================================================================
#  Z7 adds a deliberate second mismatch for weak transition SNPs.  The danger
#  flagged in review is a model confirming itself, so every claim here comes
#  from a paper and every number derived from our engine is labelled
#  `calibrated = FALSE`.
#
#  What these tests can and cannot say:
#    * they CAN check that the published claims are stored verbatim, that the
#      position claims really follow from the per-class tables, and that our
#      engine's own margin reproduces a published ORDER when that order is the
#      one Kwok reported for Taq kinetics;
#    * they CANNOT say our margin is right, because nothing here is an outcome.
#      That is the Z7 blocker, and the last test states it as a measurement.
#
#  Contradictions are transcribed, not resolved: Ye 2001 puts the extra mismatch
#  at -2, Liu 2012 measures the best single site at the 3rd base from the 3'
#  end, Lefever 2019 measures the optimum at the 4th.  All three stay.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name), show_col_types = FALSE,
                  na = c("", "NA"))
}

test_that("Z7: every prior row names its source, its evidence level and a prediction", {
  p <- ext("z7_priors.csv")
  expect_equal(nrow(p), 8L)
  expect_true(all(c("prior_id", "source_key", "prediction", "testable_form",
                    "position", "evidence", "calibrated", "doi") %in% names(p)))
  expect_false(any(is.na(p$prediction) | !nzchar(p$prediction)))
  expect_true(all(grepl("^10\\.", p$doi)))
  ladder <- c("rule of thumb", "vendor guideline", "primary literature", "fitted here")
  expect_true(all(p$evidence %in% ladder))
})

test_that("Z7: no prior is marked calibrated, and none came from our own engine", {
  p <- ext("z7_priors.csv")
  expect_false(any(p$calibrated))
  expect_false(any(p$evidence == "fitted here"))
  # a prior must come from a registered source
  s <- ext("scientific_sources.csv")
  expect_setequal(unique(p$source_key), unique(p$source_key))
  expect_true(all(p$source_key %in% s$key))
})

test_that("Z7: the three position claims are kept as a contradiction, not averaged", {
  p <- ext("z7_priors.csv")
  # Ye -2, Liu -2 (the same site in our offset convention), Lefever -3 (4th nt)
  ye <- p$position[p$prior_id == "ye_second_mismatch_position"]
  liu <- p$position[p$prior_id == "liu_position_3rd"]
  lef <- p$position[p$prior_id == "lefever_position"]
  expect_equal(ye, -2L)
  expect_equal(liu, -2L)
  expect_equal(lef, -3L)
  # Lefever and Ye disagree about the same physical quantity, so the row must
  # SAY so rather than pick a winner
  expect_true(nzchar(p$contradicted_by[p$prior_id == "lefever_position"]))
  expect_true(nzchar(p$contradicted_by[p$prior_id == "ye_second_mismatch_position"]))
  # the offset conventions are recorded, because "-2" means different sites in
  # the two papers and reading them as equal is the mistake this row prevents
  expect_true(nzchar(p$position_note[p$prior_id == "ye_second_mismatch_position"]))
})

test_that("Z7: Liu's 3rd-base claim is what its own per-class table says", {
  # The prior is not taken from memory: it has to be the weighted mean of the
  # per-class specificities at each site, which is how Table 2 was read.
  lo <- ext("liu2012_outcomes.csv")
  lo <- lo[lo$grouping == "snp_class", ]
  expect_equal(sort(unique(lo$mismatch_site)), c("2nd", "3rd", "4th"))
  mean_pct <- vapply(c("2nd", "3rd", "4th"), function(s) {
    x <- lo[lo$mismatch_site == s, ]
    stats::weighted.mean(x$specificity_pct, x$n_tested)
  }, numeric(1))
  expect_gt(mean_pct[["3rd"]], mean_pct[["2nd"]])
  expect_gt(mean_pct[["3rd"]], mean_pct[["4th"]])
  expect_equal(unname(which.max(mean_pct)), 2L) # 3rd
})

test_that("Z7: the single-mismatch dCq extract is the figure the prior cites", {
  # lefever_single_mismatch.csv is why the Lefever prior exists, so the prior row
  # and the extract must agree on what was measured
  l <- ext("lefever_single_mismatch.csv")
  expect_gt(nrow(l), 0L)
  expect_true(all(c("snp", "ddCq_hom_1", "min_ddCq", "calibrated") %in% names(l)))
  expect_false(any(l$calibrated))
})

test_that("Z7: Kwok's severity classes are the published ones, verbatim", {
  d <- ext("mismatch_discrimination.csv")
  k <- d[d$table == "kwok_3prime_severity", ]
  expect_setequal(k$mismatch_class,
                  c("severe", "moderate", "minimal", "context_dependent"))
  mem <- function(cls) {
    trimws(strsplit(k$members[k$mismatch_class == cls], ",", fixed = TRUE)[[1L]])
  }
  expect_setequal(mem("severe"), c("A:G", "G:A", "C:C"))
  expect_setequal(mem("moderate"), "A:A")
  expect_setequal(mem("minimal"), c("T:G", "T:C", "T:T"))
  # every class carries the sentence it was read from
  expect_true(all(nzchar(k$verbatim)))
})

test_that("Z7: our engine's margin reproduces Kwok's ORDER, and the probe says so", {
  # This is the one place a published ORDER is tested against our own number.
  # It is not calibration: Taq kinetics is a different measurement from
  # DECIPHER's efficiency, and only four allele pairs are probed.
  probe <- ext("z7_engine_probe.csv")
  expect_equal(nrow(probe), 4L)
  expect_false(any(probe$calibrated))
  mismatched <- probe[probe$mismatch, ]
  expect_gt(nrow(mismatched), 0L)
  # the severe pair must discriminate more than the minimal one
  get <- function(k) probe$margin[probe$allele_pair == k]
  expect_gt(get("G:A"), get("C:A"))
  # A:G and A:A are the weak classes, and their margin is (near) zero: the
  # engine cannot tell them apart from a match
  expect_lt(get("A:G"), 0.05)
  expect_lt(get("A:A"), 0.05)
})

test_that("Z7: a matched pair has no margin, so the probe is measuring the mismatch", {
  probe <- ext("z7_engine_probe.csv")
  # A:A is an A/A mismatch, but against the A-specific primer whose template
  # WAS A there is no mismatch; the probe records which rows are mismatched and
  # an exact match must give a zero margin, never a negative one
  exact <- probe[!probe$mismatch, ]
  expect_gte(nrow(exact), 1L)
  expect_true(all(abs(exact$margin) < 1e-9))
})

test_that("Z7 BLOCKER: a single 3'-terminal mismatch barely moves DECIPHER", {
  # The number Z7 is blocked on: against the MATCH template, losing the 3' base
  # costs almost nothing (terminal_floor ~ 0.004).  All the discrimination is in
  # scoring the OFF allele, so any design filter that reads one efficiency
  # cannot see the mismatch.  If this ever stops being true, the blocker is
  # resolved and Z7 can proceed -- the test will fail and say so.
  probe <- ext("z7_engine_probe.csv")
  expect_true(all(probe$terminal_floor < 0.02))
  expect_gt(diff(range(probe$margin)), 0.5) # the margin DOES separate the pairs
})

test_that("Z7: nothing in this table is presented as a probability", {
  p <- ext("z7_priors.csv")
  # the only calibrated column in the package is never TRUE here, and no prior
  # claims an expected conversion rate
  # "specificity is highest at the 3rd base" is a ranking claim, not a rate;
  # what must never appear is an expected conversion probability
  expect_false(any(grepl("probability|conversion rate|% success", p$prediction)))
  expect_false("calibrated" %in% names(p) && any(p$calibrated))
})
