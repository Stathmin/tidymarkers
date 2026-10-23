## Commun Biol 2024 (Wang et al., doi:10.1038/s42003-024-07173-7) Cas12a
## single / adjacent-double mismatch data. The extracts are produced by
## data-raw/communbiol2024_crrna.R, whose own assertions run against the
## fetched workbook. These tests check the LANDED artefacts only, so they are
## fast and need no network.

wd <- tm_extdata()
f_ratio <- file.path(wd, "communbiol2024_crrna_ratio.csv")
f_design <- file.path(wd, "communbiol2024_crrna_activators.csv")
skip_if_not(file.exists(f_ratio) && file.exists(f_design), "extracts not built")

ratio <- utils::read.csv(f_ratio, stringsAsFactors = FALSE)
design <- utils::read.csv(f_design, stringsAsFactors = FALSE)
CR <- c("WT crRNA", "M37", "M39", "M42+M43", "M45", "M48", "M63", "M65")

test_that("the counts are the paper's own: 20 single, 19 adjacent double", {
  ## the paper: "Mut1-20 and mut21-39 contained one mismatch and two adjacent
  ## mismatches relative to the on-target site"
  fig <- design[design$measured_in_fig1, ]
  expect_equal(nrow(fig), 39L)
  expect_equal(sum(fig$n_mismatch == 1L), 20L)
  expect_equal(sum(fig$adjacency == "adjacent"), 19L)
  ## and they are contiguous id blocks, not an arbitrary subset
  expect_setequal(fig$activator[fig$n_mismatch == 1L], sprintf("mut%d", 1:20))
  expect_setequal(fig$activator[fig$adjacency == "adjacent"], sprintf("mut%d", 21:39))
  ## every "adjacent" double is exactly 1 nt apart
  expect_true(all(fig$mismatch_spacing[fig$adjacency == "adjacent"] == 1L))
})

test_that("every mismatched position 1..20 of the spacer is covered", {
  singles <- design[design$n_mismatch == 1L, ]
  pos <- as.integer(singles$mismatch_position)
  expect_setequal(pos[!is.na(pos)], 1:20)
  expect_false(anyNA(pos))
})

test_that("the measurement matrix is 40 activators x 8 crRNAs and complete", {
  expect_equal(nrow(ratio), 320L)
  expect_setequal(unique(ratio$activator), c("CPSIT_0429 -WT", sprintf("mut%d", 1:39)))
  expect_setequal(unique(ratio$crrna), CR)
  expect_true(all(table(ratio$activator) == 8L))
  expect_false(anyNA(ratio$mean_fc_vs_wt))
  ## the ratio is computed from the triplicates it sits next to, not read
  expect_equal(ratio$mean_fluorescence,
    rowMeans(ratio[, c("replicate_1", "replicate_2", "replicate_3")]),
    tolerance = 1e-9
  )
  expect_equal(ratio$mean_fc_vs_wt,
    ratio$mean_fluorescence / ratio$wt_mean_fluorescence_same_crrna,
    tolerance = 1e-9
  )
})

test_that("the WT row is the reference: ratio exactly 1 in every crRNA", {
  wt <- ratio[ratio$activator == "CPSIT_0429 -WT", ]
  expect_equal(nrow(wt), 8L)
  expect_equal(wt$mean_fc_vs_wt, rep(1, 8), tolerance = 1e-12)
})

test_that("a single mismatch in the seed leaves a large mutant signal", {
  ## This is the number that justifies correcting "seed < 8 nt is intolerant".
  ## On the paper's own recommended crRNA (M39) the WORST single mismatch in
  ## positions 1-8 still transmits 37% of the WT signal (A6G, 0.373).
  s <- ratio[ratio$crrna == "M39" & ratio$n_mismatch == 1L, ]
  seed <- s[as.integer(s$mismatch_position) <= 8L, ]
  expect_equal(nrow(seed), 8L)
  expect_gt(max(seed$mean_fc_vs_wt), 0.30)
  ## and on the unmodified crRNA no single mismatch anywhere in the spacer
  ## drops below ~0.75: the unmodified guide does not discriminate at all
  unmod <- ratio[ratio$crrna == "WT crRNA" & ratio$n_mismatch == 1L, ]
  expect_gt(min(unmod$mean_fc_vs_wt), 0.70)
  expect_lt(max(unmod$mean_fc_vs_wt), 1.50)
})

test_that("two adjacent mismatches are usually, but NOT always, stronger than one", {
  ## The ordering holds as a trend and is NOT a rule. Landed numbers:
  ##   - on M39, 17 of 19 adjacent doubles fall below 0.5 while only 8 of 20
  ##     singles do;
  ##   - the double is below BOTH of its constituent single mismatches on all
  ##     8 crRNAs for 9 of 19 position pairs, and only 4 of 8 for T18C&A19G
  ##     and A19G&C20A (spacer-end pairs).
  ## The test asserts the trend plus the weakest observed double, so a future
  ## re-extraction that changes either is caught.
  wide <- reshape(ratio[ratio$n_mismatch > 0, c("crrna", "mutation", "mean_fc_vs_wt")],
    idvar = "mutation", timevar = "crrna", direction = "wide"
  )
  names(wide) <- sub("^mean_fc_vs_wt\\.", "", names(wide))
  dbl <- wide[grepl("&", wide$mutation), setdiff(names(wide), "mutation")]
  sgl <- wide[!grepl("&", wide$mutation), setdiff(names(wide), "mutation")]

  expect_equal(nrow(dbl), 19L)
  expect_equal(nrow(sgl), 20L)
  ## trend, on the recommended crRNA
  expect_lt(median(dbl[["M39"]]), median(sgl[["M39"]]))
  expect_equal(sum(dbl[["M39"]] < 0.5), 17L)
  expect_equal(sum(sgl[["M39"]] < 0.5), 8L)
  ## the strongest adjacent double in the whole set (C12A&T13C on M65, 0.0422)
  ## -- a real number, so the "doubles are stronger" claim has a magnitude
  expect_lt(min(unlist(dbl)), 0.05)
  ## the counter-example is kept visible: the weakest adjacent double still
  ## transmits more than the WT itself (A19G&C20A on M65, 1.191)
  expect_gt(max(unlist(dbl)), 1.0)
})

test_that("every row carries its outcome definition, tier and source file", {
  for (d in list(ratio, design)) {
    expect_true(all(nzchar(d$outcome_definition)))
    expect_true(all(d$evidence_tier == "measured"))
    expect_true(all(grepl("MOESM3_ESM\\.xlsx", d$source_file, fixed = FALSE)))
    expect_true(all(d$doi == "10.1038/s42003-024-07173-7"))
    expect_true(all(d$pmcid == "PMC11541961"))
  }
})
