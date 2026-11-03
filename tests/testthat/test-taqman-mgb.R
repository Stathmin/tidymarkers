# TaqMan probe chemistry modifiers (Z9.3): MGB as a sourced rule-of-thumb window
# shift, LNA honestly unsupported.  The offset is the same-sequence 3'-MGB Tm gain
# Kutyavin 2000 measured (inst/extdata/taqman_kutyavin2000_mgb_tm.csv).

test_that("probe_modification = 'none' is byte-identical to the default design", {
  f <- snp_fixture()
  m_default <- design_taqman_marker(f$genome, f$snp, flank = 150L)
  m_none <- design_taqman_marker(f$genome, f$snp, flank = 150L, probe_modification = "none")
  expect_identical(m_none, m_default)
  # the unmodified path records no offset or modifier
  expect_null(m_default$probe_tm_offset)
  expect_null(m_default$probe_modification)
})

test_that("the MGB offset is read from the Kutyavin ledger, not hardcoded", {
  csv <- utils::read.csv(tm_extdata("taqman_kutyavin2000_mgb_tm.csv"), stringsAsFactors = FALSE)
  expected <- csv$value[csv$item == "ODN1 15mer, Tm gain from MGB"]
  expect_length(expected, 1L)
  expect_gt(expected, 0)
  expect_equal(taqman_mgb_tm_offset(), expected)
})

test_that("mgb shifts the accepted probe Tm window down and records the rule of thumb", {
  f <- snp_fixture()
  offset <- taqman_mgb_tm_offset()
  # short MGB probes (12-20 nt) anneal like longer unmodified probes: the
  # effective window is the declared one shifted down by the sourced gain
  m <- design_taqman_marker(
    f$genome, f$snp,
    flank = 150L,
    probe_lengths = 12:20, probe_modification = "mgb"
  )
  expect_s3_class(m, "marker_taqman")
  # the winner sits in the shifted window, below the unshifted floor of 56 C
  expect_lt(m$probes$tm[1L], 56)
  expect_true(m$selection$gates[["gate_taqman_probe_tm"]])
  expect_identical(m$probe_modification, "mgb")
  expect_equal(m$probe_tm_offset, offset)
  expect_identical(m$probe_tm_evidence, "rule_of_thumb")
})

test_that("lna is honestly unsupported (no LNA arm in the source)", {
  f <- snp_fixture()
  expect_error(
    design_taqman_marker(f$genome, f$snp, flank = 150L, probe_modification = "lna"),
    "no LNA arm"
  )
})
