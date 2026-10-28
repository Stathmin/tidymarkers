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
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("Z7: every prior row names its source, its evidence level and a prediction", {
  p <- ext("z7_priors.csv")
  expect_equal(nrow(p), 8L)
  expect_true(all(c(
    "prior_id", "source_key", "prediction", "testable_form",
    "position", "evidence", "calibrated", "doi"
  ) %in% names(p)))
  expect_false(any(is.na(p$prediction) | !nzchar(p$prediction)))
  expect_true(all(grepl("^10\\.", p$doi)))
  ladder <- c("rule of thumb", "vendor guideline", "primary literature", "fitted here")
  expect_true(all(p$evidence %in% ladder))
})

test_that("Z7.2: every anchor names the chemistry it was measured on", {
  p <- ext("z7_priors.csv")
  expect_true(all(c("anchor_chemistry", "transfer_assumption") %in% names(p)))
  expect_false(any(is.na(p$anchor_chemistry) | !nzchar(p$anchor_chemistry)))
  # one readout per source: the anchor does not change between its predictions
  per_source <- tapply(p$anchor_chemistry, p$source_key, function(x) {
    length(unique(x))
  })
  expect_true(all(per_source == 1L))
  # four anchors, four readouts, and none of them is KASP fluorescence
  expect_equal(length(unique(p$anchor_chemistry)), 4L)
  expect_false(any(grepl("KASP", p$anchor_chemistry, ignore.case = TRUE)))
  chem <- function(key) unique(p$anchor_chemistry[p$source_key == key])
  expect_match(chem("kwok1990"), "Taq")
  expect_match(chem("lefever2019"), "real-time")
  expect_match(chem("liu2012_arms"), "allele-specific")
  expect_match(chem("ye2001"), "ARMS")
  # the transfer is recorded as an assumption, identical on every row
  expect_length(unique(p$transfer_assumption), 1L)
  expect_match(p$transfer_assumption[[1L]], "assumption")
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
  # Liu (3rd base from the 3' end) and Lefever (4th nt) both sit at -3 in our
  # offset convention; only Ye 2001 puts the extra mismatch at -2, as a rule
  ye <- p$position[p$prior_id == "ye_second_mismatch_position"]
  liu <- p$position[p$prior_id == "liu_position_3rd"]
  lef <- p$position[p$prior_id == "lefever_position"]
  expect_equal(ye, -2L)
  # Liu's own per-class table files it at -3 (see the ledger row); the prior row
  # used to say -2 with a note, which put one paper on both sides of the split
  expect_equal(abs(liu), 3L)
  expect_equal(abs(lef), 3L)
  # the contradiction is recorded rather than resolved away
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
  expect_setequal(
    k$mismatch_class,
    c("severe", "moderate", "minimal", "context_dependent")
  )
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

## --- Z7.5: the added mismatch must not weaken another gate -------------------

#  The roadmap's worry is a closed gate paying for another: a design that gains
#  the deliberate mismatch must still be graded by the same engine, the same
#  checks and the same terminal-base gate.  Two tests, one per layer:
#    * the primer layer - one live pair, with and without the mismatch, at the
#      SAME geometry, so the only difference is the one added base;
#    * the design layer - the full ARMS pipeline on one synthetic locus, run both
#      ways, so the gate vocabulary is compared across real designs.
#  Statuses are NOT required to be equal: a single base can move Tm, GC or a
#  dimer.  What must be equal is the gate SET (ids, order, stage, level, limit)
#  and its applicability, plus the terminal-base gate the mismatch sits behind.

test_that("Z7.5: an added mismatch is one 3'-distal base and drops no gate", {
  skip_if_no_oligoarrayaux()
  gen <- syn_genome(seed = 11L, n = 2000L)
  w <- Biostrings::DNAString(gen$seq)
  pos <- 1000L
  ref <- substr(gen$seq, pos, pos)
  # A:G is one of Kwok's classes, so the terminal gate carries a real value
  # (`severe`) instead of `unknown` - "did not move" then means something.
  alt <- "G"
  d <- design_as_primers(
    w, pos,
    ref = ref, alt = alt, lengths = 20L, strands = "+",
    stab_ats = c(NA, abs(ARMS_EXTRA_MISMATCH_POSITION))
  )
  mm_rows <- d[!is.na(d$stab_at), , drop = FALSE]
  expect_gt(nrow(mm_rows), 0L)
  mm <- mm_rows[1L, ]
  # the SAME geometry without the deliberate mismatch
  plain <- d[
    is.na(d$stab_at) & d$as_length == mm$as_length &
      d$tmpl_start == mm$tmpl_start & d$strand == mm$strand, ,
    drop = FALSE
  ]
  expect_equal(nrow(plain), 1L)
  plain <- plain[1L, ]
  # the added mismatch is a single base, at the decided distance, in BOTH primers
  # of the pair; the allele-specific 3' base itself is untouched
  chars <- function(x) strsplit(x, "", fixed = TRUE)[[1L]]
  at <- mm$as_length - abs(ARMS_EXTRA_MISMATCH_POSITION) + 1L
  expect_equal(which(chars(mm$ref_primer) != chars(plain$ref_primer)), at)
  expect_equal(which(chars(mm$alt_primer) != chars(plain$alt_primer)), at)
  expect_equal(
    substr(mm$ref_primer, mm$as_length, mm$as_length),
    substr(plain$ref_primer, mm$as_length, mm$as_length)
  )

  # the same pipeline before and after: identical gate ids, order, stage, level
  # and limit, and no applicable gate may become inapplicable
  com <- "ACGTACGTACGTACGTACGT"
  a0 <- make_assay(plain$ref_primer, plain$alt_primer, com, ref = ref, alt = alt)
  a1 <- make_assay(mm$ref_primer, mm$alt_primer, com, ref = ref, alt = alt)
  e0 <- marker_design_checks(a0)
  e1 <- marker_design_checks(a1)
  expect_identical(e0$criterion, e1$criterion)
  expect_identical(e0$stage, e1$stage)
  expect_identical(e0$level, e1$level)
  expect_identical(e0$limit, e1$limit)
  expect_false(any(e1$status == "not_applicable" & e0$status != "not_applicable"))
  # the mismatch reaches the gates: AS Tm is computed from the mismatched pair,
  # and it moved - it cannot be a stale, un-mismatched read
  as_tm <- e1$observed[e1$criterion == "AS Tm (C)"]
  expect_equal(as_tm, mean(primer_tm(c(mm$ref_primer, mm$alt_primer))[1:2]))
  expect_gt(abs(as_tm - e0$observed[e0$criterion == "AS Tm (C)"]), 0.01)
  # the terminal-base gate is the one the mismatch is layered ON, and it must not
  # move: the 3' base is the same in both primers
  term <- "AS 3' mismatch base pair"
  expect_equal(as.character(e0$status[e0$criterion == term]), "pass")
  expect_equal(as.character(e1$status[e1$criterion == term]), "pass")
  # ZW0 M2 still holds with the mismatch present: swapping REF/ALT moves neither
  # the verdict nor the gate set
  sw <- a1
  sw$snp$ref <- a1$snp$alt
  sw$snp$alt <- a1$snp$ref
  sw$assay$as_ref <- a1$assay$as_alt
  sw$assay$as_alt <- a1$assay$as_ref
  sw$thermo$tm <- c(a1$thermo$tm[2L], a1$thermo$tm[1L], a1$thermo$tm[3L])
  e_sw <- marker_design_checks(sw)
  expect_equal(evaluation_verdict(e1), evaluation_verdict(e_sw))
  expect_identical(e1$criterion, e_sw$criterion)
})

test_that("Z7.5: the full ARMS pipeline gates a mismatched design the same way", {
  skip_if_no_oligoarrayaux()
  # The same synthetic locus, designed twice through the real pipeline: once with
  # no deliberate mismatch and once with the decided one.  Both are live markers
  # (outer pair, inner pair, band sizes); the gate vocabulary must not shrink.
  gen <- syn_genome(seed = 11L, n = 2000L)
  g <- marker_genome(gen$path)
  pos <- 1000L
  ref <- substr(gen$seq, pos, pos)
  alt <- if (ref %in% c("A", "G")) "C" else "T"
  snp <- tibble::tibble(chrom = "Chr1", pos = as.integer(pos), ref = ref, alt = alt)
  plain <- suppressWarnings(design_arms_marker(g, snp, flank = 300L, stab_ats = NA))
  mism <- suppressWarnings(design_arms_marker(
    g, snp,
    flank = 300L, stab_ats = abs(ARMS_EXTRA_MISMATCH_POSITION)
  ))
  testthat::skip_if(
    is.null(plain) || is.null(mism),
    "the synthetic locus did not design both ways"
  )
  # the mismatch is recorded on the winner, and it is the policy's default row
  expect_identical(plain$inner$ref$stab_policy, "none")
  expect_identical(mism$inner$ref$stab_policy, "default")
  expect_equal(mism$inner$ref$stab_at, abs(ARMS_EXTRA_MISMATCH_POSITION))
  # the mismatched inner primer is among the oligos the evaluator scores
  expect_true(mism$inner$ref$primer %in% marker_primers(mism))
  # same pipeline -> identical gate ids, stage, level and limit, and no gate
  # becomes inapplicable because the mismatch is present
  e0 <- marker_design_checks(plain)
  e1 <- marker_design_checks(mism)
  expect_identical(e0$criterion, e1$criterion)
  expect_identical(e0$stage, e1$stage)
  expect_identical(e0$level, e1$level)
  expect_identical(e0$limit, e1$limit)
  expect_false(any(e1$status == "not_applicable" & e0$status != "not_applicable"))
})
