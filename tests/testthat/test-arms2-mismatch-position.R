# ARMS/AS-PCR family, batch 2: the mismatch-position ledger.
#
# The point of these tests is that a hand-edit to arms_mismatch_position.csv
# breaks something. They assert the invariants the generator established, and
# - critically - that the file does NOT collapse the -2 vs -3 contradiction
# into a single preferred position.

pos_csv <- function() {
  path <- tm_extdata("arms_mismatch_position.csv")
  skip_if_not(file.exists(path), "ledger CSV not built")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

test_that("the ledger keeps both sides of the mismatch-position contradiction", {
  led <- pos_csv()
  measured <- led[led$evidence_class != "rule", ]
  expect_true(-2L %in% measured$position)
  expect_true(-3L %in% measured$position)
  # A single averaged/preferred position would collapse this to one level.
  expect_gte(length(unique(measured$position)), 2L)
})

test_that("only measurements (not rules) are counted as support", {
  led <- pos_csv()
  rules <- led[led$evidence_class == "rule", ]
  expect_true(all(is.na(rules$n_tested) | rules$n_tested == 0L))
  # every counted row must still carry its per-position evidence size
  expect_false(any(is.na(led$n_sources_at_position)))
  expect_true(all(led$n_sources_at_position >= 1L))
})

test_that("positions are negative-from-3'-end and inside the primer body", {
  led <- pos_csv()
  expect_true(all(led$position %in% c(-1L, -2L, -3L, -4L, -5L)))
  expect_true(all(led$position < 0L))
  # a deliberate mismatch cannot sit beyond a 20-mer
  expect_true(all(abs(led$position) <= 20L))
})

test_that("Liu 2012 is filed at the position its own percentages support", {
  led <- pos_csv()
  liu <- led[led$key == "liu2012", ]
  expect_equal(nrow(liu), 1L)
  expect_equal(liu$position, -3L)
  # the sentence printed in `discrimination` must rank -3 above -2 and -4
  pcts <- as.numeric(regmatches(
    liu$discrimination,
    gregexpr("[0-9]+\\.[0-9]+", liu$discrimination)
  )[[1L]])
  expect_equal(which.max(pcts), 2L) # -2, -3, -4 in that order
})

test_that("no row was invented without a source file or a quote", {
  led <- pos_csv()
  expect_true(all(nzchar(led$source_file)))
  expect_true(all(led$outcome_definition == led$outcome_definition[1L]))
  # ye2001's text is not in either tree, so it must not claim verification
  ye <- led[led$key == "ye2001", ]
  expect_false(ye$verified_here)
  unverified <- led[!led$verified_here, ]
  expect_true(all(!nzchar(unverified$quote)))
})
