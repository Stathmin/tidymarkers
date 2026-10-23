# Z7 — where the deliberate (extra) mismatch goes
#
# The decision is recorded in data-raw/proposed/arms_mismatch_position.decision.md:
# -3 is the measured default (Liu 2012 + Lefever 2019), -2 is a fallback for a
# `context_dependent` terminal class (Ye 2001, a rule).  These tests hold the
# decision to the evidence it cites, so neither the code nor the ledger can drift.

ext <- function(name) {
  utils::read.csv(
    testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
}

test_that("the decided default is -3 and the fallback is -2", {
  expect_equal(ARMS_EXTRA_MISMATCH_POSITION, -3L)
  expect_equal(ARMS_EXTRA_MISMATCH_FALLBACK, -2L)
})

test_that("the fallback is spent on the weak terminal class only", {
  expect_equal(arms_extra_mismatch_position("severe"), -3L)
  expect_equal(arms_extra_mismatch_position("moderate"), -3L)
  expect_equal(arms_extra_mismatch_position("minimal"), -2L)
  expect_equal(arms_extra_mismatch_position("context_dependent"), -3L)
  # no severity class -> no position policy: `NA`, never a guess
  expect_true(is.na(arms_extra_mismatch_position(NA_character_)))
  expect_true(is.na(arms_extra_mismatch_position("")))
  expect_true(is.na(arms_extra_mismatch_position("made-up")))
})

test_that("the ledger still puts the measurements at -3 and none at -2", {
  # the decision's whole basis: two controlled scans at -3, nothing controlled
  # at -2.  If the ledger ever gains a -2 measurement, the default is revisitable.
  led <- ext("arms_mismatch_position.csv")
  expect_true(all(c("position", "evidence_class") %in% names(led)))
  meas <- led$position[led$evidence_class == "measurement"]
  expect_equal(sort(unique(meas)), -3L)
  expect_equal(sum(meas == -3L), 2L)
  expect_false(any(led$position[led$evidence_class == "measurement"] == -2L))
  # the -2 rows are the rule and the primer scan, and they are NOT measurements
  expect_setequal(
    led$evidence_class[led$position == -2L],
    c("rule", "primer_scan")
  )
})

test_that("z7_priors puts the measured optimum at -3, and that row is Liu's", {
  p <- ext("z7_priors.csv")
  liu <- p[p$prior_id == "liu_position_3rd", ]
  lef <- p[p$prior_id == "lefever_position", ]
  ye <- p[p$prior_id == "ye_second_mismatch_position", ]
  # Liu's "3rd base from the 3' end" is -3 in our offset convention; the ledger
  # row (arms_mismatch_position.csv) already files it there.  The prior row used
  # to say -2 with a note, which put the same paper on both sides.
  expect_equal(liu$position, -3L)
  expect_equal(liu$evidence, "primary literature")
  expect_true(liu$measured)
  expect_equal(lef$position, -3L)
  expect_equal(ye$position, -2L)
  expect_equal(ye$evidence, "rule of thumb")
  expect_false(ye$measured)
  # and the contradiction is still recorded rather than resolved away
  expect_true(nzchar(ye$contradicted_by))
})

test_that("the AS engine ranks -3 ahead of -2 at equal discrimination", {
  # a synthetic window long enough for a 20-mer AS primer on both strands
  w <- Biostrings::DNAString(paste(rep("ACGT", 30L), collapse = ""))
  out <- design_as_primers(
    w, 40L,
    ref = substr(w, 40L, 40L), alt = "A", lengths = 20L,
    stab_ats = c(NA, 2L, 3L), strands = c("+", "-")
  )
  expect_s3_class(out, "tbl_df")
  expect_true("stab_policy" %in% names(out))
  rank <- match(out$stab_policy, c("default", "fallback", "none", "other"))
  # the order is the policy, not discrimination: default first, then fallback
  expect_equal(rank, sort(rank))
  # every -3 row is the default; a -2 row is the fallback only where the
  # terminal base is a weak one
  expect_true(all(out$stab_policy[!is.na(out$stab_at) & out$stab_at == -3L] == "default"))
  expect_true(all(!out$stab_policy[!is.na(out$stab_at) & out$stab_at == -3L] == "fallback"))
  mm2 <- out[!is.na(out$stab_at) & out$stab_at == -2L, ]
  expect_true(all(mm2$stab_policy == ifelse(mm2$weak_terminal, "fallback", "other")))
  # the policy reorders, it does not delete: both positions are still returned
  expect_setequal(unique(stats::na.omit(out$stab_at)), c(2L, 3L))
})

test_that("a caller who asks for another position gets it, marked as not-policy", {
  w <- Biostrings::DNAString(paste(rep("ACGT", 30L), collapse = ""))
  out <- design_as_primers(
    w, 40L,
    ref = substr(w, 40L, 40L), alt = "A", lengths = 20L,
    stab_ats = 4L, strands = "+"
  )
  expect_true(all(out$stab_policy == "other"))
  expect_true(all(out$stab_at == 4L))
})

test_that("the signed policy constant is not a stab_ats distance", {
  # ARMS_EXTRA_MISMATCH_POSITION is a signed 3'-relative offset (-3); the
  # designer's `stab_ats` are positive distances from the 3' end (3).  Passing
  # the constant straight in used to place the mismatch past the 5' end in
  # silence, so a technology composed with it did not design at all.
  expect_error(
    design_as_primers("ACGTACGTACGTACGTACGT", 10L, "A", "C",
      strands = "+", stab_ats = ARMS_EXTRA_MISMATCH_POSITION
    ),
    "must be positive"
  )
  # the distance form is accepted and its magnitude is the policy's position
  expect_identical(abs(ARMS_EXTRA_MISMATCH_POSITION), 3L)
})

test_that("the chosen position survives into the returned marker (Z7.1)", {
  skip_if_no_oligoarrayaux()
  g <- ath_genome()
  snp <- pick_target_snp(ensembl_variants("4", 1500000L, 1560000L), g, ensembl_gene("AT4G00650"))
  m <- suppressWarnings(design_arms_marker(g, snp, flank = 200L))
  testthat::skip_if(is.null(m), "no ARMS design at this locus")
  # the winner names the policy row that produced it, on both inner primers
  for (side in c("ref", "alt")) {
    p <- m$inner[[side]]
    expect_true("stab_at" %in% names(p))
    expect_true("stab_policy" %in% names(p))
    expect_true(p$stab_policy %in% c("none", "default", "fallback", "other"))
  }
  # a position that is reported is the one the policy chose
  if (!is.na(m$inner$ref$stab_at)) {
    expect_true(m$inner$ref$stab_policy %in% c("default", "fallback"))
    expect_equal(abs(m$inner$ref$stab_at), abs(ARMS_EXTRA_MISMATCH_POSITION))
  }
})

test_that("a technology carries the decided position through marker_technology()", {
  skip_if_no_oligoarrayaux()
  g <- ath_genome()
  snp <- pick_target_snp(ensembl_variants("4", 1500000L, 1560000L), g, ensembl_gene("AT4G00650"))
  # the distance form the designer takes
  marker_technology("arms_z71_dist",
    base = "arms",
    params = list(stab_ats = abs(ARMS_EXTRA_MISMATCH_POSITION)), overwrite = TRUE
  )
  m <- suppressWarnings(marker_design_at(g, snp, scheme = "arms_z71_dist", flank = 200L))
  testthat::skip_if(is.null(m), "no ARMS design at this locus")
  expect_identical(m$inner$ref$stab_at, abs(ARMS_EXTRA_MISMATCH_POSITION))
  expect_identical(m$inner$ref$stab_policy, "default")
})

## --- the Liu 2012 tetra-primer extract ---------------------------------------

# 1686 published ASP primers, tabulated with the position of the deliberate
# mismatch. The table is checked against its own sequences, not against the
# paper's prose: a wrong 3' base would invalidate every position claim.
test_that("liu2012 primers: the 3' base and the last four are the sequence's own", {
  d <- ext("liu2012_primers.csv")
  expect_equal(nrow(d), 1686L)
  n <- nchar(d$forward)
  expect_equal(d$forward_3p_base, substr(d$forward, n, n))
  expect_equal(d$forward_last4, substr(d$forward, n - 3L, n))
  # the deliberate mismatch sits at a named position, never "1st" (that is the
  # allele-specific base) and never unstated
  expect_setequal(unique(d$mismatch_site), c("2nd", "3rd", "4th"))
  # four canonical SNP classes, each written as a pair of pairs
  expect_setequal(
    unique(d$snp_class),
    c("A/C | G/T", "A/G | C/T", "A/T | T/A", "C/G | G/C")
  )
  expect_true(all(nchar(d$mismatch_pair) == 2L))
  expect_identical(anyDuplicated(d$primer), 0L)
})
