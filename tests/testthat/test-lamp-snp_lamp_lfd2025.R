# SNP-LAMP + dipstick (rpoB S531L): what the landed CSVs must keep true.
# Data-only: reads the CSVs from inst/extdata, so it is fast and never needs
# .tmp/lit. The generator data-raw/snp_lamp_lfd2025.R is the re-parse check.

lfd_read <- function(name) {
  utils::read.csv(testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE
  )
}
lfd_files <- c(
  "lamp_snp_lamp_lfd2025_primers.csv",
  "lamp_snp_lamp_lfd2025_detection.csv",
  "lamp_snp_lamp_lfd2025_conditions.csv",
  "lamp_snp_lamp_lfd2025_contradictions.csv"
)

test_that("the LFD extract files exist and carry the required columns", {
  for (f in lfd_files) {
    d <- lfd_read(f)
    expect_gt(nrow(d), 0L)
    expect_true(all(c("outcome_definition", "evidence_tier", "source_file") %in% names(d)))
    expect_true(all(d$evidence_tier == "weak"))
    expect_true(all(nzchar(d$source_file)))
  }
})

test_that("primers are ACGT, the tags are recorded, and B2 is the allele arm", {
  p <- lfd_read(lfd_files[1])
  expect_identical(p$primer, c("F3", "B3", "FIP", "BIP-bio", "LB-fit-C"))
  expect_true(all(grepl("^[ACGT]+$", p$sequence)))
  expect_identical(p$five_prime_tag, c("none", "none", "none", "biotin", "fitC"))
  b2 <- p$b2_arm[p$primer == "BIP-bio"]
  expect_true(length(b2) == 1L || all(is.na(b2)))
  expect_gte(nchar(b2), 4L)
  # the B2 arm is the 3' end of BIP: it must end where BIP's sequence ends
  if (length(b2) == 1L && !is.na(b2)) {
    expect_true(endsWith(p$sequence[p$primer == "BIP-bio"], b2))
  }
  # no per-base mismatch identity is claimed (it is in an image)
  expect_match(p$mismatch_position_cited[p$primer == "BIP-bio"], "image")
})

test_that("the dipstick table is 6/6 and 0/4 with consistent totals", {
  d <- lfd_read(lfd_files[2])
  expect_identical(d$dipstick_positive, c(6L, 0L))
  expect_identical(d$dipstick_negative, c(0L, 4L))
  expect_identical(d$n_dst, c(6L, 4L))
  expect_identical(sum(d$n_dst), 10L)
  expect_true(all(d$dipstick_positive + d$dipstick_negative == d$n_dst))
})

test_that("conditions are the paper's, and no allele gap is claimed", {
  cnd <- lfd_read(lfd_files[3])
  amp <- cnd[cnd$item == "amplification", ]
  expect_true(all(amp$value == 60))
  expect_true(any(cnd$item == "dipstick read"))
  # the point of this source: a CALL, not a magnitude
  d <- lfd_read(lfd_files[2])
  expect_match(unique(d$outcome_definition), "no time-to-positive")
  expect_false(any(c("time_to_positive_min", "signal_ratio") %in% names(d)))
})

test_that("the contradictions are recorded, not resolved", {
  k <- lfd_read(lfd_files[4])
  expect_gte(nrow(k), 3L)
  expect_true(all(grepl("recorded", k$resolution)))
  expect_true(any(grepl("positions 1 and 3", k$text_version)))
})
