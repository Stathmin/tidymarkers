# ARMS/AS-PCR family, batch 2: the blocked barley source and the cycle-decay
# negative result. Both files exist to record an ABSENCE, so their tests check
# that the absence is still real and still documented.

blocked_csv <- function() {
  path <- tm_extdata("arms_barley_arms2004_blocked.csv")
  skip_if_not(file.exists(path), "barley blocked record not built")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

cycle_csv <- function() {
  path <- tm_extdata("arms_cycle_decay_negative.csv")
  skip_if_not(file.exists(path), "cycle-decay record not built")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

test_that("barley_arms2004 is recorded as blocked with no numbers", {
  b <- blocked_csv()
  expect_equal(nrow(b), 1L)
  expect_equal(b$status, "blocked")
  expect_equal(b$evidence_tier, "blocked")
  expect_true(grepl("barley", b$citation, ignore.case = TRUE))
  expect_true(is.na(b$source_file))
  # a blocked record must not smuggle in an outcome
  expect_true(grepl("NONE", b$outcome_definition, fixed = TRUE))
})

test_that("barley_arms2004 still has no fetched full text here", {
  lit <- test_path("../../.tmp/lit/barley_arms2004")
  if (!dir.exists(lit)) skip("no fetch directory")
  files <- list.files(lit, recursive = TRUE)
  big <- files[file.size(file.path(lit, files)) > 20000]
  expect_length(big, 0L)
})

test_that("the cycle-decay record is a genuine negative result", {
  c <- cycle_csv()
  expect_gt(nrow(c), 0L)
  expect_true(all(!c$has_cycle_resolved_discrimination_series))
  expect_true(all(is.na(c$cycle_decay_rate_per_cycle)))
  expect_true(all(c$evidence_tier == "weak"))
})

test_that("every searched source states whether it discusses cycles", {
  c <- cycle_csv()
  expect_true(all(!is.na(c$talks_about_cycles)))
  expect_type(c$talks_about_cycles, "logical")
  # at least one source must actually engage with cycle number, otherwise the
  # search itself was too thin to support a negative result
  expect_true(any(c$talks_about_cycles))
  expect_true(all(nzchar(c$note)))
})
