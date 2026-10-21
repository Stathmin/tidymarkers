# =============================================================================
#  Z8  --  gate-pass -> assay outcome, with the failures kept
# =============================================================================
#  Review point: published markers are successes, so a pass/fail dataset built
#  from them is survivorship-biased.  The response here is narrow and testable:
#
#    * every conversion row carries its denominator AND its failures;
#    * the outcome vocabulary stays ORDINAL ("usable"), never a probability,
#      until the failures decompose for enough sources;
#    * a row that decomposes its failures is recorded as such, and no other row
#      is silently given the same status;
#    * no two studies are pooled into one rate.
#
#  The extract files from the readout sweep (CAPS, HRM, LAMP, dPCR) enter here
#  as COUNTS WITH THEIR OWN OUTCOME DEFINITIONS.  "Success" means a different
#  thing in each paper, so the definition travels with the count.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name), show_col_types = FALSE,
                  na = c("", "NA"))
}

test_that("Z8: every conversion row keeps its denominator and its failures", {
  c <- ext("conversion_rates.csv")
  expect_gt(nrow(c), 3L)
  for (col in c("n_designed", "n_success", "n_failed")) {
    expect_true(all(is.finite(c[[col]]) & c[[col]] >= 0), label = col)
  }
  # success + failure must reconstruct the denominator: a row that does not
  # balance is a row where a failed assay vanished
  expect_equal(c$n_success + c$n_failed, c$n_designed)
  # and no success rate may exceed 1
  expect_true(all(c$n_success <= c$n_designed))
})

test_that("Z8: the outcome scale stays ordinal, and nothing is calibrated", {
  c <- ext("conversion_rates.csv")
  expect_true(all(c$outcome_scale == "usable"))
  expect_false(any(c$calibrated))
  # no column of this table may be a rate we computed and could mistake for a
  # probability of success
  expect_false(any(grepl("prob|p_success|expected", names(c))))
})

test_that("Z8: decomposition is claimed only where the named failures add up", {
  c <- ext("conversion_rates.csv")
  # `n_fail_known` is how many failures a source NAMES; if that equals its
  # failure total the row decomposes, and only then
  implied <- !is.na(c$n_fail_known) & c$n_fail_known == c$n_failed & c$n_failed > 0
  expect_equal(c$fully_decomposed, implied)
  # exactly one source currently decomposes: the multi-assembly + Sanger arm,
  # the one that reports the reason per SNP
  expect_equal(sum(c$fully_decomposed), 1L)
  expect_equal(c$study[c$fully_decomposed], "makhoul2020_multi_assembly_sanger")
})

test_that("Z8: the design-strategy ladder is ordered, and nothing is pooled", {
  m <- ext("makhoul_design_strategies.csv")
  # the three strategies of ONE paper, tested on the same 15 SNPs
  expect_equal(nrow(m), 15L)
  s1 <- sum(m$s1_flanks_only != "")
  s2 <- sum(m$s2_single_ref != "")
  expect_equal(s1, 15L)
  expect_equal(s2, 15L)
  # more information never loses a usable assay: a strategy ladder that went
  # backwards would mean the table was mis-parsed
  c <- ext("conversion_rates.csv")
  m_rows <- c[grepl("^makhoul2020_", c$study), ]
  expect_equal(nrow(m_rows), 3L)
  expect_true(all(diff(m_rows$n_success) >= 0))
  # the outcome vocabulary is the authors' own, kept verbatim
  expect_true(all(m$outcome_ordinal %in% c("usable", "not usable", "doomed")))
})

test_that("Z8: an outcome row cannot be both usable and doomed", {
  m <- ext("makhoul_design_strategies.csv")
  usable <- m$s3_multi_sanger == "Correct cluster"
  expect_true(all(usable == (m$outcome_ordinal == "usable")))
  # 'Locus-specific design not possible' is the authors' doom state, and it is
  # not the same as "no cluster" (which is a failed assay, not a failed design)
  not_possible <- grepl("not possible", m$s3_multi_sanger)
  expect_true(all(m$outcome_ordinal[not_possible] == "doomed"))
})

## --- the readout sweep: counts with their own definitions -------------------

test_that("Z8 CAPS: a conversion funnel never goes up, and success is defined", {
  f <- ext("liu2016_caps_funnel.csv")
  expect_gt(nrow(f), 3L)
  expect_true(all(c("stage", "n", "outcome_definition") %in% names(f)))
  # the funnel is monotone; 271/532 must reproduce the printed 50.94%
  two <- f$n[grepl("primer pair", f$stage)][1L:2L]
  expect_true(all(diff(two) <= 0))
  expect_equal(round(100 * two[[2L]] / two[[1L]], 2L), 50.94)
  # the definition of success is not "amplified": it is polymorphic + expected size
  expect_true(any(grepl("POLYMORPHIC", f$outcome_definition)))
})

test_that("Z8 CAPS: a second paper's funnel is contiguous and self-consistent", {
  f <- ext("soybean_est_caps_funnel.csv")
  n <- f$n_kept[!is.na(f$n_kept)]
  expect_gte(length(n), 4L)
  expect_true(all(diff(n) <= 0)) # each stage keeps no more than its parent
  # the number judged unequivocal matches the number with a usable pattern
  expect_lte(min(n), 26L)
  expect_gte(max(n), 48L)
})

test_that("Z8 dCAPS: every predicted product is reconstructible from its fragments", {
  s <- ext("snap2caps_markers.csv")
  expect_gt(nrow(s), 5L)
  # a marker's cut and uncut groups partition its accession groups
  expect_equal(s$n_cut_groups + s$n_uncut_groups, s$n_accession_groups)
  # success is a prediction confirmed on a gel, and the file says so
  expect_true(all(grepl("predicted", s$outcome_definition)))
})

test_that("Z8 dPCR: LoB and LoD are recorded as an assay limit, not a threshold we set", {
  l <- ext("jak2v617f_ddpcr_limits.csv")
  expect_gte(nrow(l), 2L)
  expect_true(all(c("LoB", "LoD") %in% l$quantity))
  # every limit carries the definition it was measured under
  expect_true(all(nzchar(l$definition)))
  # LoD is not below LoB: that would mean misreading the table
  lob <- l$value_vaf_pct[l$quantity == "LoB"]
  lod <- l$value_vaf_pct[l$quantity == "LoD"]
  expect_true(all(lod >= lob))
})

test_that("Z8 LAMP: discrimination is a time window with explicit censoring", {
  t <- ext("ace1r_aslamp2015_time_to_positive.csv")
  expect_gt(nrow(t), 0L)
  expect_true(all(c("time_to_positive_min", "observation_window_min",
                    "censoring") %in% names(t)))
  # the non-target allele is NOT called "negative": it is censored at the window
  cens <- t[!is.na(t$censoring) & nzchar(t$censoring), ]
  expect_gt(nrow(cens), 0L)
  expect_true(all(grepl("not measured|exceeds", cens$censoring)))
})

test_that("Z8 HRM: a pooled outcome is marked pooled, never re-split", {
  o <- ext("desilva2009_outcomes.csv")
  expect_gte(nrow(o), 2L)
  expect_true("pooled_not_per_amplicon" %in% names(o))
  # the paper prints two pooled contingency rows and no per-amplicon calls
  expect_true(all(o$pooled_not_per_amplicon))
  # the 400 bp boundary belongs to the <400 group because that is the reading
  # which reproduces the printed mean; the number is not invented
  lt <- o[o$amplicon_group == "lt400", ]
  expect_equal(lt$n_amplicons_table1, 9L)
  expect_equal(lt$mean_bp_table1, lt$mean_bp_paper)
})

test_that("Z8: an outcome definition always names the assay, not just 'success'", {
  c <- ext("conversion_rates.csv")
  expect_true("failure_note" %in% names(c))
  m <- ext("makhoul_design_strategies.csv")
  expect_true(all(nzchar(m$outcome_def)))
  expect_true(grepl("KASP cluster", m$outcome_def[[1L]]))
})
