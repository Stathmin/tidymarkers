# =============================================================================
#  ZW3 external oracle  --  hm2025: HRM amplicon Tm, measured vs predicted
# =============================================================================
#  PLoS ONE 2025, DOI 10.1371/journal.pone.0321885, Table 9 (psaA-2 half).
#
#  This is an AMPLICON-level measurement, not an oligo Tm, so it is NOT a test
#  of primer_tm() on the primer duplex.  It is carried for two reasons:
#    * it is a real measured tabular Tm set with an internal arithmetic check
#      (printed Difference == |predicted - measured|), so a hand-edit breaks it;
#    * it records the ORDERING claim those authors make (predicted below
#      measured for every amplicon) against their own model, so a future
#      external comparison has a reference point.
#  Dia18S rows are excluded: degenerate primers (contain N).

hm2025_path <- function() {
  tm_extdata("tm_hm2025_amplicon_psaa2.csv")
}

test_that("ZW3/hm2025: amplicon Tm table is arithmetically self-consistent", {
  skip_if_not(file.exists(hm2025_path()), "tm_hm2025 CSV not present")
  d <- utils::read.csv(hm2025_path(), stringsAsFactors = FALSE)

  expect_equal(nrow(d), 5L)
  expect_true(all(d$primer_pair == "psaA-2"))

  # the printed difference must equal |predicted - measured|
  expect_equal(
    round(abs(d$paper_predicted_tm_c - d$measured_tm_c), 2),
    round(d$paper_difference_c, 2)
  )

  # physically possible PCR-product window
  expect_true(all(d$measured_tm_c > 70 & d$measured_tm_c < 95))

  # the authors' own directional claim
  expect_true(all(d$paper_predicted_tm_c < d$measured_tm_c))

  # no degenerate (N-bearing) primer slipped in
  expect_false(any(grepl("[^ACGT]", paste0(d$forward_primer, d$reverse_primer))))
})
