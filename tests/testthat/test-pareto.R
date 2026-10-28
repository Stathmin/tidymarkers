# Pareto front over marker candidates: minimise cost and linkage risk, maximise
# quality.  A dominated point never enters the front; the recommended point is
# always one of the front points.

pareto_points <- function() {
  tibble::tibble(
    label = c("A", "B", "C"),
    cost = c(1, 2, 0.5),
    quality = c(0.9, 0.8, 0.5),
    linkage_risk = c(0.1, 0.2, 0.5)
  )
}

test_that("a point dominated on every axis is not on the front", {
  r <- marker_pareto(pareto_points())
  # A is cheaper, better and lower-risk than B: B is dominated
  expect_false(r$on_front[r$label == "B"])
  expect_true(r$on_front[r$label == "A"])
  expect_true(r$on_front[r$label == "C"])
  expect_equal(sum(r$on_front), 2L)
})

test_that("being worse on one axis does not remove a point that wins on another", {
  pts <- tibble::tibble(
    label = c("cheap", "good"),
    cost = c(0.5, 2),
    quality = c(0.4, 0.95),
    linkage_risk = c(0.3, 0.3)
  )
  r <- marker_pareto(pts)
  expect_true(all(r$on_front))
})

test_that("the recommended point is always on the front", {
  set.seed(11L)
  for (k in seq_len(6L)) {
    pts <- tibble::tibble(
      label = paste0("p", seq_len(8L)),
      cost = runif(8L, 0.5, 10),
      quality = runif(8L, 0, 1),
      linkage_risk = runif(8L, 0, 0.5)
    )
    r <- marker_pareto(pts)
    expect_equal(sum(r$recommended), 1L)
    expect_true(any(r$recommended[r$on_front]))
  }
})

test_that("a lower-risk point at the same cost and quality is preferred", {
  pts <- tibble::tibble(
    label = c("risky", "safe"),
    cost = c(1, 1),
    quality = c(0.8, 0.8),
    linkage_risk = c(0.4, 0.01)
  )
  r <- marker_pareto(pts)
  expect_true(r$recommended[r$label == "safe"])
  expect_false(r$on_front[r$label == "risky"])
})

test_that("a marker_solved supplies its report costs", {
  fake <- structure(
    list(
      report = tibble::tibble(
        label = c("t1", "t2"), chrom = "Chr1", pos = c(100L, 300L),
        cost = c(1, 4)
      ),
      plan = list(region = NULL)
    ),
    class = "marker_solved"
  )
  r <- marker_pareto(fake, quality = c(0.9, 0.4), linkage_risk = c(0, 0))
  expect_equal(r$cost, c(1, 4))
  expect_equal(r$quality, c(0.9, 0.4))
  expect_equal(r$label, c("t1", "t2"))
})

test_that("linkage risk grows with distance and is a rule of thumb", {
  fake <- structure(
    list(
      report = tibble::tibble(
        label = c("near", "far"), chrom = "Chr1", pos = c(1000L, 101000L),
        cost = c(1, 1)
      ),
      plan = list(region = NULL)
    ),
    class = "marker_solved"
  )
  r <- marker_pareto(fake,
    quality = c(0.8, 0.8), anchor = 1000L,
    recombination = 1e-8
  )
  expect_equal(r$linkage_risk[r$label == "near"], 0)
  expect_gt(r$linkage_risk[r$label == "far"], 0)
})

test_that("malformed candidates are rejected, not guessed", {
  expect_error(marker_pareto(tibble::tibble(cost = 1)), "quality")
  expect_error(marker_pareto(pareto_points(), quality = 1), "quality")
  expect_error(
    marker_pareto(tibble::tibble(cost = -1, quality = 1, linkage_risk = 0)),
    "cost"
  )
  expect_error(
    marker_pareto(tibble::tibble(cost = 1, quality = 2, linkage_risk = 0)),
    "quality"
  )
})
