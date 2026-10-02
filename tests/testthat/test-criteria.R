# The criterion registry is the single vocabulary for selection + acceptance.

crit_mock <- function() {
  p <- c(AS_ref = "CATGGACCTAGCATTGACGTA", AS_alt = "CATGGACCTAGCATTGACGTT", COM = "GCTTAGCCATGGTCAATGCTA")
  list(
    assay = tibble::tibble(
      amplicon = 90L, as_ref = p[[1]], as_alt = p[[2]], com = p[[3]], com_tm = 60,
      as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L, as_dropout = 0, com_dropout = 0
    ),
    thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), hairpin_dg = c(-1, -1, -1)),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3), b = rep(names(p), times = 3),
      dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
    ),
    amplicons = tibble::tibble(primer = c("AS_ref", "AS_alt"), n_products = c(1L, 0L)),
    efficiency = tibble::tibble(discrimination = 0.9)
  )
}

test_that("kasp_criteria declares gates, tiers and checks with metadata", {
  crit <- kasp_criteria()
  expect_setequal(unique(crit$role), c("gate", "tier", "check"))
  expect_true(all(c("stage", "id", "label", "level", "hard", "limit", "order") %in% names(crit)))
  expect_true(all(crit$stage %in% c("site", "as", "pair", "assay")))
  expect_true(any(crit$hard))
  # the gates are the constants, and the hard ones are exactly those off the ladder
  pgate <- crit$id[crit$role == "gate" & crit$stage == "pair"]
  expect_setequal(pgate, paste0("pair:gate:", KASP_GATES_PAIR))
  expect_setequal(
    crit$id[crit$role == "gate" & crit$hard],
    c(
      paste0("as:gate:", setdiff(KASP_GATES_AS, KASP_RELAX_AS)),
      paste0("pair:gate:", setdiff(KASP_GATES_PAIR, KASP_RELAX_PAIR))
    )
  )
})

test_that("kasp_design_checks returns exactly the registry's checks, in order", {
  crit <- kasp_criteria()
  checks <- crit[crit$role == "check", ]
  ch <- kasp_design_checks(crit_mock())
  expect_equal(ch$criterion, checks$label)
  expect_equal(ch$limit, checks$limit)
  expect_true(all(ch$pass)) # a clean mock passes every acceptance check
})
