# What a TaqMan design may claim about mismatch
# (data-raw/proposed/taqman_status.decision.md). No model change: the constants
# keep their values and their status changes to "engineering default".

tq_read <- function(name) {
  utils::read.csv(tm_extdata(name),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
}

tq_source <- function() {
  f <- testthat::test_path("..", "..", "R", "marker_taqman.R")
  # the R source is not shipped with the installed package, so this assertion
  # only exists in a checkout
  skip_if_not(file.exists(f), "package source not present (installed check)")
  readLines(f)
}

test_that("the mismatch threshold is declared an engineering default, not a finding", {
  src <- tq_source()
  i <- grep("^TAQMAN_MISMATCH_DG", src)
  expect_length(i, 1L)
  # the comment block above the constant must say both things: it is a rule the
  # package chose, and no open paper measures it
  head_block <- paste(src[max(1L, i - 20L):i], collapse = " ")
  expect_match(head_block, "PACKAGE RULE, NOT A MEASUREMENT")
  expect_match(head_block, "No open paper")
  # and it must name the limit: the exonuclease step is not modelled
  expect_match(head_block, "exonuclease")
  # the value itself is unchanged by this decision
  expect_identical(eval(parse(text = src[i])[[1L]][[3L]]), 1.5)
})

test_that("the TaqMan readout is not simulated, so no discrimination is claimed", {
  expect_match(readout_assumption("taqman"), "not simulated|not modelled")
  m <- structure(list(scheme = "taqman"), class = "marker_taqman")
  r <- marker_readout(m)
  expect_false(r$simulated)
  expect_null(r$table)
})

test_that("centrality is a preference: real assays place the allele off centre", {
  p <- tq_read("taqman_citrus2020_probes.csv")
  expect_true("centrality" %in% names(p))
  cen <- p$centrality[is.finite(p$centrality)]
  expect_gt(length(cen), 5L)
  # decision clause 3: no landed assay is central (1.0), so centrality must not
  # be graded as a hard requirement
  expect_true(all(cen < 1))
  expect_lt(min(cen), 0.8)
  # five of eight FAM probes are below the 18:25 default range
  expect_true(any(is.finite(p$probe_len_fam) & p$probe_len_fam < 18L))
})

test_that("the one evidenced TaqMan failure is a footprint gap, not a mismatch", {
  a <- tq_read("taqman_apol1_2019_failure.csv")
  # decision clause 4: the 6-bp deletion (G2) removes bases under the probe, so
  # G2/G2 gives no amplification and G1GM/G2 is mis-called
  noamp <- a[a$predesigned_call == "No amplification", ]
  mis <- a[a$predesigned_call != a$custom_call & a$predesigned_call != "No amplification", ]
  expect_equal(noamp$freq_pct, 1.4)
  expect_equal(mis$freq_pct, 5.9)
  expect_true(all(a$n_samples == 1110))
  expect_true(all(grepl("deletion", a$outcome_definition)))
  # the artefact says it is a failure MODE, not mismatch discrimination
  expect_false(any(grepl("delta_ct|dct", names(a), ignore.case = TRUE)))
})
