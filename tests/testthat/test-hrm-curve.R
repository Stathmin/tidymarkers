# HRM melting-curve realism: what the landed dTm rows do and do not support.
#
# The decision (data-raw/proposed/hrm_resolvability.decision.md) grades an HRM
# design by the measured dTm of its bimorphism class against a 0.5 C bar. This
# file asks the narrower question the decision leaves open: does any landed row
# tie a dTm to a *stated instrument or buffer resolution*, i.e. is 0.5 C an
# instrument-scoped default or a demonstrated floor? Along the way it pins each
# table's own printed totals and its smallest called dTm, so a later extract
# cannot quietly widen the claim.

hrm_csv <- function(f) {
  utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", f),
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
}

test_that("jsbbs2018 dTm distribution is a screening set, not a detection floor", {
  e1 <- hrm_csv("hrm_jsbbs2018_dtm_exp1.csv")
  s <- hrm_csv("hrm_jsbbs2018_summary.csv")
  val <- function(q) s$value[s$quantity == q]

  expect_equal(nrow(e1), 38L)
  expect_equal(val("n_markers_exp1"), 38)
  # the summary's own descriptives are these 38 rows, exactly
  expect_equal(min(e1$dtm_obs), val("dtm_exp1_min_c"), tolerance = 1e-9)
  expect_equal(stats::median(e1$dtm_obs), val("dtm_exp1_median_c"),
    tolerance = 1e-9
  )
  expect_equal(max(e1$dtm_obs), val("dtm_exp1_max_c"), tolerance = 1e-9)
  # the smallest dTm the authors still called polymorphic, and how many of the
  # 38 screened markers sit under the 0.5 C bar the decision applies
  expect_equal(min(e1$dtm_obs), 0.1976, tolerance = 1e-3)
  expect_equal(sum(e1$dtm_obs < HRM_DTM_RESOLVABLE), 8L)
  # a screened set is not a detection limit, and the extract must say so in
  # words; otherwise the 0.198 C minimum would read as a validated floor
  min_row <- s[s$quantity == "dtm_exp1_min_c", ]
  expect_match(min_row$outcome_definition, "NOT a validated detection floor",
    fixed = TRUE
  )
  expect_match(min_row$outcome_definition, "no false-negative rate",
    fixed = TRUE
  )
})

test_that("the 0.5 C bar is the paper's screening outcome, above its own low call", {
  e1 <- hrm_csv("hrm_jsbbs2018_dtm_exp1.csv")
  e3 <- hrm_csv("hrm_jsbbs2018_dtm_exp3.csv")
  bm <- hrm_csv("hrm_jsbbs2018_dtm_by_bimorphism.csv")

  # the 5 class rows are pooled from the two experiments and from nothing else
  expect_equal(sum(bm$n), 56L)
  expect_equal(sum(bm$n), nrow(e1) + nrow(e3))
  expect_equal(sort(unique(bm$source_file)), unique(bm$source_file))
  expect_true(all(bm$evidence_tier == "measured"))

  # Exp. 3 is the pre-filtered co-dominant set: every one of its 18 markers
  # clears the bar, so the pooled minimum can only come from Exp. 1
  expect_equal(nrow(e3), 18L)
  expect_equal(mean(e3$dtm_obs), 0.7393, tolerance = 1e-3)
  expect_equal(min(e3$dtm_obs), 0.5871, tolerance = 1e-3)
  expect_equal(sum(e3$dtm_obs < HRM_DTM_RESOLVABLE), 0L)
  expect_equal(
    min(bm$dtm_min),
    min(e1$dtm_obs),
    tolerance = 1e-9
  )
  # ... and that pooled minimum is BELOW the bar, so a bar is not enforced by
  # the source that states it: it is a rule the authors applied afterwards
  expect_lt(min(bm$dtm_min), HRM_DTM_RESOLVABLE)
})

test_that("a class can contain sub-bar markers and still be graded resolvable", {
  bm <- hrm_csv("hrm_jsbbs2018_dtm_by_bimorphism.csv")
  # four of the five landed classes have at least one member under the bar
  expect_equal(sum(bm$dtm_min < HRM_DTM_RESOLVABLE), 4L)
  for (i in seq_len(nrow(bm))) {
    parts <- strsplit(bm$bimorphism[[i]], "/", fixed = TRUE)[[1L]]
    r <- hrm_resolvability(parts[[1L]], parts[[2L]])
    expect_equal(r$dtm_min, bm$dtm_min[[i]], tolerance = 1e-3)
    expect_match(r$evidence, "measured", fixed = TRUE)
  }
  # only C/G, whose single measured member is under the bar, is condemned
  expect_equal(hrm_resolvability("C", "G")$status, "below_threshold")
  # G/T's mean clears the bar by 0.019 C: the thinnest margin in the extract
  gt <- hrm_resolvability("G", "T")
  expect_equal(gt$status, "resolvable")
  expect_equal(gt$dtm_mean - HRM_DTM_RESOLVABLE, 0.0187, tolerance = 1e-3)
})

test_that("gundry2008: a 0.02 C inter-class gap was called at 64.3 percent", {
  a <- hrm_csv("hrm_gundry2008_amplicon_tm.csv")
  h <- hrm_csv("hrm_gundry2008_homozygote_accuracy.csv")

  expect_equal(length(unique(a$amplicon)), 8L)
  expect_equal(nrow(a), 16L) # two homozygote rows per amplicon
  # the paper's own totals: six loci, 1147 samples
  expect_equal(nrow(h), 6L)
  expect_equal(sum(h$n), 1147)
  expect_equal(sort(h$n), c(42L, 49L, 184L, 248L, 296L, 328L))

  # this is the counterexample to a dTm floor: the smallest reported dTm
  # (BX647987, 0.02 C after calibration) is ALSO the locus with the worst
  # accuracy (64.3%), while 0.13 C at CPS1 is called at 100%
  expect_equal(min(a$dtm_before), 0.09, tolerance = 1e-9)
  expect_equal(min(a$dtm_after), 0.02, tolerance = 1e-9)
  expect_equal(range(a$dtm_before), c(0.09, 0.62), tolerance = 1e-9)
  expect_equal(range(a$dtm_after), c(0.02, 0.62), tolerance = 1e-9)
  bx <- h$gene == "BX647987"
  cps <- h$gene == "CPS1"
  expect_equal(h$accuracy_before_pct[bx], 69.0)
  expect_equal(h$accuracy_after_pct[bx], 64.3)
  expect_lt(h$accuracy_after_pct[bx], h$accuracy_after_pct[cps])
  expect_equal(h$accuracy_after_pct[cps], 100)
  # so a smaller dTm did not guarantee a worse call: BX647987 at 0.02 C was
  # called worse than CPS1 at 0.15 C, and CHK2 at 0.62 C was called perfectly,
  # which is what "not monotone in dTm" means here
  expect_equal(h$accuracy_after_pct[h$gene == "CHK2"], 100)
})

test_that("the landings disagree on the smallest dTm called, so no floor is demonstrated", {
  bm <- hrm_csv("hrm_jsbbs2018_dtm_by_bimorphism.csv")
  g <- hrm_csv("hrm_gundry2008_amplicon_tm.csv")
  ij <- hrm_csv("hrm_ijms2022_dtm_single_snp_controls.csv")

  # smallest reported NON-REFERENCE dTm per fetched source, in C: 0.198
  # (jsbbs2018, probe Tm), 0.02 (Gundry 2008, amplicon Tm), 0.4 (ijms2022, the
  # lowest non-zero step of the 3-locus C>T panel, whose reference row is 0)
  expect_equal(min(bm$dtm_min), 0.1976, tolerance = 1e-3)
  expect_equal(min(g$dtm_after), 0.02, tolerance = 1e-9)
  expect_equal(min(ij$dtm_c[ij$n_t_alleles > 0]), 0.4, tolerance = 1e-9)
  # they span a 10-fold (0.02-0.20 C) range of "smallest called", which is what
  # an instrument-dependent floor looks like; a physical floor would be one
  # number. So the bar must stay a class-level default, not a property of dTm.
  expect_gt(min(bm$dtm_min) / min(g$dtm_after), 9)

  # and each table's OWN header must scope it to its platform, or the
  # cross-source comparison above would be reading apples against oranges
  expect_match(g$outcome_definition[[1L]], "LightScanner", fixed = TRUE)
  expect_match(ij$outcome_definition[[1L]], "THIS amplicon", fixed = TRUE)
})

test_that("ijms2022 gives partitions and dTm steps, never a resolution figure", {
  d <- hrm_csv("hrm_ijms2022_discrimination.csv")
  p <- hrm_csv("hrm_ijms2022_partition_2snp.csv")
  u <- hrm_csv("hrm_ijms2022_unresolvable_3snp.csv")
  ij <- hrm_csv("hrm_ijms2022_dtm_single_snp_controls.csv")

  expect_equal(d$n_controls, rep(36L, 6L))
  expect_equal(
    d$n_uniquely_identified[d$assay == "2 SNPs" & d$readout == "HRM"],
    1L
  )
  expect_equal(
    d$n_uniquely_identified[d$assay == "3 SNPs" & d$readout == "HRM"],
    22L
  )
  # the 2-SNP classes are a partition of all 36 controls, both ways
  expect_equal(sum(p$n_controls), 36L)
  expect_equal(sum(p$class_size * p$n_classes), 36L)
  expect_equal(sum(p$n_classes), 13L) # HRM alone found one class per 13 groups
  # the 3-SNP misses are stated as pairs: 7 pairs hold 14 of the 36 controls
  expect_equal(sum(!is.na(u$control_a)), 7L)
  expect_equal(sum(!is.na(u$control_a)) * 2, 36L - 22L)

  # dTm here is a step between controls, not a resolution: the reference is
  # 0 and the ladder is these six steps
  expect_equal(sort(ij$dtm_c[ij$n_t_alleles > 0]), c(0.4, 0.6, 0.7, 1.0, 1.3, 1.7))
  expect_equal(ij$dtm_c[ij$n_t_alleles == 0], 0)
})

test_that("NEGATIVE: no landed HRM row states an instrument or buffer resolution", {
  files <- list.files(
    testthat::test_path("..", "..", "inst", "extdata"),
    pattern = "^hrm_.*[.]csv$",
    full.names = TRUE
  )
  expect_gt(length(files), 10L)

  # (a) no column is a resolution, an instrument or a buffer field; the source
  # file and sheet are recorded, but not what the rig could resolve
  cols <- unique(unlist(lapply(files, function(f) names(hrm_csv(basename(f))))))
  expect_false(any(grepl(
    "resolution|instrument|buffer|platform", cols,
    ignore.case = TRUE
  )))

  # (b) the platform name sits ONLY in prose or in code metadata, never in a
  # data field. Exactly one HRM CSV mentions an instrument at all (Gundry
  # 2008's LightScanner, inside outcome_definition, because its Tm column was
  # read off that rig), and no row states what dTm that rig could resolve.
  platform <- "LightCycler|LightScanner|Rotor-Gene"
  naming <- function(f) {
    d <- hrm_csv(basename(f))
    any(vapply(
      d, function(x) is.character(x) && any(grepl(platform, x)),
      logical(1L)
    ))
  }
  named <- vapply(files, naming, logical(1L))
  expect_equal(sum(named), 1L)
  expect_identical(basename(files[named]), "hrm_gundry2008_amplicon_tm.csv")
  # and that mention is prose, not a resolution column value
  gdef <- hrm_csv("hrm_gundry2008_amplicon_tm.csv")$outcome_definition
  expect_true(all(grepl("LightScanner", gdef, fixed = TRUE)))

  # jsbbs2018 scopes its 0.5 C rule to an instrument in PROSE the landed files
  # do not carry; the code carries it as a string, so the bar cannot be read
  # back out of the data as if the data had measured a resolution
  expect_match(hrm_resolvability("C", "T")$evidence, "LightCycler",
    fixed = TRUE
  )

  # (c) the one source that could have supplied a cross-source table is closed
  # access and stays secondary_cited; it must never be promoted to measured
  lw <- hrm_csv("hrm_liew2004_restated_bands.csv")
  expect_equal(nrow(lw), 6L)
  expect_true(all(lw$evidence_tier == "secondary_cited"))
  expect_false(any(lw$doi == "10.1373/clinchem.2004.032136"))
  expect_true(all(lw$doi == "10.1093/nar/gkn204")) # read via Gundry 2008 prose

  # (d) therefore the decision text stands unchanged: the bar is a
  # platform-scoped default, and no test may turn it into a measured floor
  expect_equal(HRM_DTM_RESOLVABLE, 0.5)
  r <- hrm_resolvability("C", "G")
  expect_equal(r$status, "below_threshold")
  # unmeasured pairs are still unknown, never failed: nothing landed resolves them
  for (pair in c("A/T", "A/A", "T/T")) {
    parts <- strsplit(pair, "/", fixed = TRUE)[[1L]]
    expect_equal(hrm_resolvability(parts[[1L]], parts[[2L]])$status, "unknown")
  }
})
