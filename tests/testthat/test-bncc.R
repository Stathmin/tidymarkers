# BNCC (ZW4.6 phase C): band priority, weighted Chebyshev regret, Pareto front.

test_that("every declared band partition is a permutation of the stage's tiers", {
  for (stage in names(BNCC_BANDS)) {
    pos <- unlist(BNCC_BANDS[[stage]], use.names = FALSE)
    expect_identical(sort(pos), seq_along(pos), info = stage)
  }
  # an undeclared stage keeps the identity partition (plain lexicographic)
  tiers <- c("a", "b", "c")
  expect_identical(unname(unlist(bncc_bands("not_a_stage", tiers))), tiers)
  # a partition that does not cover the tiers fails loudly, not silently
  expect_error(bncc_bands("as", c("a", "b", "c")))
})

test_that("regret is bounded, monotone, and keeps unknown margins unknown", {
  x <- c(-2, 0, 1, 3)
  r <- bncc_regret(x, ideal = 3, nadir = -2)
  expect_equal(r, c(1, 0.6, 0.4, 0))
  # monotone: a larger margin never has more regret
  expect_true(all(diff(r) <= 0))
  # beyond the ideal / nadir clamps
  expect_equal(bncc_regret(c(10, -10), ideal = 3, nadir = -2), c(0, 1))
  # no spread at all is "equally good", not an error
  expect_equal(bncc_regret(c(5, 5)), c(0, 0))
  # an unknown margin is not evidence
  expect_true(is.na(bncc_regret(c(1, NA, 3))[2]))
})

test_that("a band is scored by its bottleneck, but keeps a Pareto front", {
  d <- tibble::tibble(
    id = 1:2,
    tier1 = c(3, 1), # regret 0 vs 1 (band A)
    tier2 = c(3, 1),
    tier3 = c(1, 3) # the first is worse in band B
  )
  ranked <- bncc_rank(d, NULL, c("tier1", "tier2", "tier3"),
    bands = list(A = c("tier1", "tier2"), B = "tier3")
  )
  # band A is lexicographic first: the row with the better worst-case Band A wins
  expect_equal(ranked$data$id[1], 1L)
  expect_true(ranked$data$.pareto[1])
  # the two rows trade A against B, so both sit on the front
  expect_equal(sum(ranked$data$.pareto), 2L)
})

test_that("a candidate that wins band A beats one that only wins band B", {
  d <- tibble::tibble(
    id = 1:2,
    tier1 = c(5, 0), tier2 = c(5, 0), # band A: row 1 far better
    tier3 = c(0, 5) # band B: row 2 better, but it cannot compensate
  )
  ranked <- bncc_rank(d, NULL, c("tier1", "tier2", "tier3"),
    bands = list(A = c("tier1", "tier2"), B = "tier3")
  )
  expect_equal(ranked$data$id[1], 1L)
})

test_that("the Chebyshev score sits between the bottleneck and the mean", {
  d <- tibble::tibble(tier1 = 2, tier2 = 8) # regrets 0.2 and 0.8 on a 0..10 range
  s <- bncc_scores(d, c("tier1", "tier2"), list(A = c("tier1", "tier2")),
    ideal = c(tier1 = 10, tier2 = 10), nadir = c(tier1 = 0, tier2 = 0)
  )
  expect_equal(s$.bncc_A, 0.8) # the worst component, not the mean 0.5
  expect_gte(s$.bncc_A, mean(c(0.2, 0.8)))
})

test_that("an empty tier is skipped, not dropped, inside a band", {
  d <- tibble::tibble(tier1 = c(NA, 2), tier2 = c(3, NA))
  s <- bncc_scores(d, c("tier1", "tier2"), list(A = c("tier1", "tier2")))
  expect_equal(s$.bncc_A[1], 0) # the known tier carries the row
  expect_equal(s$.bncc_A[2], 0)
  # a band where every tier is unknown has no opinion
  d2 <- tibble::tibble(tier1 = NA_real_, tier2 = NA_real_)
  s2 <- bncc_scores(d2, c("tier1", "tier2"), list(A = c("tier1", "tier2")))
  expect_true(is.na(s2$.bncc_A))
  # and an all-NA band sorts last
  d3 <- tibble::tibble(id = 1:2, tier1 = c(NA, 1))
  ranked <- bncc_rank(d3, "pair", "tier1", bands = list(A = "tier1"))
  expect_equal(ranked$data$id[1], 2L)
})

test_that("the diagonal case reproduces the lexicographic order exactly", {
  set.seed(1)
  d <- tibble::tibble(
    id = 1:20,
    tier1 = runif(20), tier2 = runif(20), tier3 = runif(20)
  )
  tiers <- c("tier1", "tier2", "tier3")
  ident <- as.list(tiers)
  names(ident) <- LETTERS[seq_along(ident)]
  by_bncc <- bncc_rank(d, NULL, tiers, bands = ident)$data
  by_lex <- order_by_margins(d, tiers, character(0))
  # one tier per band with equal weights is the current order, not a new one
  expect_equal(
    match(by_lex$id, d$id),
    match(by_bncc$id, d$id)
  )
})

test_that("select_with_ladder keeps lexicographic by default and can opt in", {
  d <- tibble::tibble(
    tier1 = c(1, 5), tier2 = c(5, 1),
    g = c(TRUE, TRUE), m_g = c(1, 1), id = 1:2
  )
  args <- list(
    data = d, gates = "g", relax_order = character(0), tiers = c("tier1", "tier2"),
    ties = "id", gate_margin = c(g = "m_g"), gate_floor = c(g = 0)
  )
  lex <- do.call(select_with_ladder, args)
  expect_equal(lex$best$id, 2L)
  bncc <- do.call(select_with_ladder, c(args, list(
    bncc = list(stage = NULL, bands = list(A = "tier1"))
  )))
  expect_equal(bncc$best$id, 2L)
  expect_true(".pareto" %in% names(bncc$best))
})

test_that("the calibration hook is inert: no criterion row is fitted here", {
  expect_false(BNCC_CALIBRATION$fitted)
  expect_length(BNCC_CALIBRATION$weights, 0L)
  spec <- tidymarkers:::criterion_spec(KASP_LIMITS)
  expect_identical(any(as.character(spec$evidence) == "fitted_here"), FALSE)
})
