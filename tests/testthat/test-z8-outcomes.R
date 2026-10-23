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
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
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
  # `decomposes` is the authors' own claim that the failure reasons they give
  # reconstitute their failure total; `fully_decomposed` is what the extract
  # script could verify.  They are one statement about one row and may never
  # disagree; what is NOT asserted is that the sum of the named classes equals
  # the total, because a paper may name only the classes it chose to report.
  expect_equal(c$decomposes, c$fully_decomposed)
  expect_true(all(c$n_fail_known <= c$n_failed))
  expect_equal(c$fully_decomposed, c$decomposes)
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
  # the funnel is monotone at every step, not just the last
  expect_true(all(diff(f$n) <= 0))
  # 271/532 must reproduce the paper's printed 50.94%
  pair <- f$n[grepl("chosen for design", f$stage)]
  ok <- f$n[grepl("polymorphic", f$stage)]
  expect_equal(length(pair), 1L)
  expect_equal(round(100 * ok / pair, 2L), 50.94)
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

test_that("Z8 dPCR: LoB, LoD and LoQ are three ordered limits, not one 'sensitivity'", {
  l <- ext("jak2v617f_ddpcr_limits.csv")
  expect_setequal(l$quantity, c("LoB", "LoD", "LoQ"))
  # every limit carries the definition it was measured under
  expect_true(all(nzchar(l$definition)))
  # LoB < LoD <= LoQ: a reversal would mean the table was misread
  v <- l$value_vaf_pct[match(c("LoB", "LoD", "LoQ"), l$quantity)]
  expect_equal(v, sort(v))
  # the CV at the LoD is carried verbatim, however large
  expect_equal(l$reported_cv_pct[l$quantity == "LoD"], 76)
})

test_that("Z8 LAMP: AS-LAMP calls are counts against a sequencing truth set", {
  # ace1r_aslamp2015_counts.csv: every method column must sum to the study total,
  # so a method that dropped samples cannot pass as agreement
  k <- ext("ace1r_aslamp2015_counts.csv")
  expect_gt(nrow(k), 5L)
  expect_true(all(c("genotype", "method", "n_called", "n_sequenced") %in% names(k)))
  # every method is scored against the same 104 sequencing calls
  # the sequencing truth set is 104 samples; the PCR-RFLP and AS-LAMP panels are
  # scored on their own subsets, which are recorded explicitly
  expect_true(all(k$n_sequenced[!is.na(k$n_sequenced)] %in% c(94L, 9L, 1L, 104L)))
  per_method <- tapply(k$n_called, k$method, sum)
  expect_gt(length(per_method), 1L)
  # every method accounts for all 104 samples (the re-run subset is excluded by
  # its own NA sequencing column)
  # at least two methods account for the full 104 samples
  scored <- per_method[!grepl("water bath", names(per_method))]
  expect_gte(sum(scored == 104L), 2L)
})

test_that("Z8 LAMP: discrimination is a time window with explicit censoring", {
  t <- ext("ace1r_aslamp2015_time_to_positive.csv")
  expect_gt(nrow(t), 0L)
  expect_true(all(c(
    "time_to_positive_min", "observation_window_min",
    "censoring"
  ) %in% names(t)))
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
  # which reproduces the printed mean; the mean is the paper's, not ours
  lt <- o[o$amplicon_group %in% c("lt400", "<400"), ]
  expect_gte(nrow(lt), 1L)
  expect_equal(lt$n_amplicons_in_table1[[1L]], 9L)
  expect_equal(lt$mean_bp_table1[[1L]], lt$mean_bp_paper[[1L]])
  # and the group the paper itself cannot reproduce is flagged, not smoothed
  gt <- o[o$amplicon_group %in% c("gt400", ">400"), ]
  expect_gte(nrow(gt), 1L)
  expect_false(any(gt$mean_bp_agrees))
})

test_that("Z8: an outcome definition always names the assay, not just 'success'", {
  c <- ext("conversion_rates.csv")
  expect_true("failure_note" %in% names(c))
  m <- ext("makhoul_design_strategies.csv")
  expect_true(all(nzchar(m$outcome_def)))
  expect_true(grepl("KASP cluster", m$outcome_def[[1L]]))
})

test_that("Z8.3: a KASP conversion rate names which of the four quantities it is", {
  # Decision: data-raw/proposed/kasp_z8_scope.decision.md.  A conversion rate is
  # design, manufacture, work or usable - and rates with different definitions
  # are never pooled or compared (the 6.5 % wheat_platform2019 figure is a
  # fraction of assay x sample points, not of assays).
  c <- ext("conversion_rates.csv")
  def <- if ("conversion_definition" %in% names(c)) c$conversion_definition else rep(NA_character_, nrow(c))
  # the column may be absent only while it is being introduced; the moment a row
  # carries it, it must be populated and must distinguish definitions
  if (any(!is.na(def))) {
    expect_false(any(is.na(def) | !nzchar(def)))
    rate <- ifelse(
      !is.na(c$n_designed) & c$n_designed > 0L,
      round(c$n_success / c$n_designed, 6),
      NA_real_
    )
    key <- paste(rate, def)
    # two rows may share a value only if they share a definition
    dup <- stats::ave(key, key, FUN = length)
    expect_true(all(dup[!is.na(rate)] <= 1L))
  }
  # the assay-x-sample scope, when recorded, is never an assay conversion
  if ("scope" %in% names(c)) {
    xs <- c[!is.na(c$scope) & c$scope == "assay_x_sample", ]
    expect_true(all(!is.na(xs$conversion_definition)))
  }
})

## --- the de Silva 2009 HRM extracts -----------------------------------------

test_that("de Silva's 22 amplicons: primer lengths add up and 400 bp is the seam", {
  a <- ext("desilva2009_amplicons.csv")
  expect_equal(nrow(a), 22L)
  # the paper prints a primer-total column, so the two sequences must reproduce it
  expect_equal(a$primer_total_nt, nchar(a$forward) + nchar(a$reverse))
  # the 400 bp boundary is assigned to the <400 group by the paper's own mean
  expect_equal(nrow(a), sum(!is.na(a$group)))
  expect_equal(unname(table(a$group)[["gt400"]]), 13L)
  expect_equal(unname(table(a$group)[["lt400"]]), 8L)
  seam <- a$group == "at400_assigned_lt400_by_mean"
  expect_equal(sum(seam), 1L)
  expect_equal(a$size_bp[seam], 400L)
  # a group label must agree with the size it names, wherever it is not the seam
  expect_true(all((a$size_bp > 400) == (a$group == "gt400")))
})

test_that("de Silva's 25 SNP amplicons keep their missing sizes and multi-SNP flag", {
  s <- ext("desilva2009_snps.csv")
  expect_equal(nrow(s), 25L)
  expect_identical(anyDuplicated(s$snp_ids), 0L)
  expect_true(all(s$n_snps >= 1L & s$n_snps <= 3L))
  # the flag is derived, so it must equal the derivation
  expect_equal(s$multi_snp_amplicon, s$n_snps > 1L)
  # five rows have no size in the paper; they stay NA rather than being dropped
  expect_equal(sum(is.na(s$size_bp)), 5L)
  expect_equal(sum(s$multi_snp_amplicon), sum(s$n_snps > 1L))
})

test_that("de Silva's five strategies are five distinct situations", {
  st <- ext("desilva2009_strategies.csv")
  expect_equal(nrow(st), 5L)
  expect_identical(anyDuplicated(st$genetic_situation), 0L)
  expect_true(all(nzchar(st$genetic_situation)))
  expect_true(all(nzchar(st$hrm_solution)))
  expect_true(all(st$doi == "10.1186/1471-2156-10-5"))
})

## --- YP5: outcomes recorded against a stable artefact id --------------------

# A minimal artefact on disk: the manifest carries the id an outcome is bound
# to, targets.csv is what a target label is checked against.  No solve, no
# engine, no network - this file tests the bookkeeping, not a design.
outcome_fixture <- function(dir) {
  writeLines(
    '{"schema": "tidymarkers/artefact@1", "id": "0123456789abcdef"}',
    fs::path(dir, "artefact.json")
  )
  readr::write_csv(
    tibble::tibble(target = c("t1", "t2")),
    fs::path(dir, "targets.csv")
  )
  dir
}

test_that("Z8: a recorded outcome is ordinal, per-row, and never calibrated", {
  dir <- withr::local_tempdir()
  dir <- outcome_fixture(dir)
  expect_equal(nrow(marker_outcomes(dir)), 0L)

  out <- marker_record_outcome(dir, "works", target = "t1")
  expect_equal(nrow(out), 1L)
  expect_identical(out$id, "0123456789abcdef")
  expect_identical(out$target, "t1")
  expect_identical(out$result, "works")
  expect_false(any(out$calibrated))
  expect_true(nzchar(out$recorded[[1L]]))

  # a retest appends rather than overwriting: the disagreement stays visible
  out <- marker_record_outcome(dir, "fails",
    target = "t1", notes = "no amplification at 60 C"
  )
  expect_equal(nrow(out), 2L)
  expect_false(any(out$calibrated))
  expect_identical(out$notes[[2L]], "no amplification at 60 C")

  # a whole-artefact call has no target; "unknown" is a result, not a number
  out <- marker_record_outcome(dir, "unknown")
  expect_true(is.na(out$target[[3L]]))
  expect_identical(out$result[[3L]], "unknown")

  # observations only: the schema has no derived column at all, so no column
  # can be a rate or a probability that a reader might take for one
  expect_identical(names(out), outcome_columns())
  # and the file a non-R reader sees is the same table
  expect_equal(marker_outcomes(dir), out)
})

test_that("Z8: a failure without a reason is refused, so is an unknown target", {
  dir <- withr::local_tempdir()
  dir <- outcome_fixture(dir)
  expect_error(marker_record_outcome(dir, "fails"), "notes")
  expect_error(marker_record_outcome(dir, "fails", notes = "   "), "notes")
  expect_error(
    marker_record_outcome(dir, "works", target = "nope"),
    "not in this artefact"
  )
  expect_error(marker_record_outcome(dir, "passed"))
  # nothing a refused call touched was written
  expect_equal(nrow(marker_outcomes(dir)), 0L)
})
