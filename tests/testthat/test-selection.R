# The selection logic: gate the grave sins, then rank on tier margins.

# A minimal candidate table with the columns grade_assay_pairs()/select need.
pair_table <- function(discrimination, com_homodimer = -1, com_hairpin = -0.5,
                       com_tm = 60, amplicon = 80, n_var = 0L, as_n_3p = 0L,
                       run = 0L, ssr = 0L) {
  tibble::tibble(
    discrimination = discrimination,
    com_homodimer = com_homodimer, com_hairpin = com_hairpin,
    com_tm = com_tm, amplicon = amplicon, n_var = as.integer(n_var),
    as_n_3p = as.integer(as_n_3p), run = as.integer(run), ssr = as.integer(ssr),
    com_start = seq_along(discrimination)
  )
}

test_that("grade_assay_pairs builds gates and bottleneck tier margins", {
  d <- pair_table(c(0.9, 0.9), com_homodimer = c(-4, -8), com_hairpin = c(-1, -6))
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_true(all(c("tier1", "tier2", "m_disc") %in% names(g$data)))
  # tier2 is the bottleneck of the two structure margins
  expect_equal(
    g$data$tier2,
    pmin(c(-4, -8) - KASP_LIMITS$homodimer_dg, c(-1, -6) - KASP_LIMITS$hairpin_dg)
  )
  expect_true(all(g$data$gate_discrimination))
})

test_that("a clean pairing beats a prettier one with a COM footprint variant", {
  d <- pair_table(
    discrimination = c(0.99, 0.60),
    com_tm = c(59.95, 61.5),
    n_var = c(2L, 0L) # the first is otherwise perfect but sits on known SNPs
  )
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  pick <- select_with_ladder(g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties)
  expect_equal(pick$best$n_var, 0L)
  expect_length(pick$relaxed, 0L)
  expect_length(pick$violated, 0L)
})

test_that("the ladder relaxes only as far as needed, least grave first", {
  # every candidate sits on known SNPs: the COM-footprint gate must give way,
  # and only after the milder gates have been tried.  The two rows are otherwise
  # equal, so the least-laden pairing wins on the footprint tie-break.
  d <- pair_table(discrimination = c(0.9, 0.9), n_var = c(3L, 1L))
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  pick <- select_with_ladder(g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties)
  expect_true("gate_com_footprint" %in% pick$relaxed)
  expect_false("gate_as_window" %in% pick$relaxed) # the gravest is never reached
  expect_true("gate_com_footprint" %in% pick$violated)
  expect_equal(pick$best$n_var, 1L) # still the least-laden pairing
})

test_that("margins are compared within a tier, not summed", {
  # same tier2 bottleneck (min = 1) but very different totals -> still tied
  d <- pair_table(
    discrimination = c(0.9, 0.9),
    com_homodimer = c(-8, -8), # margin 1 both
    com_hairpin = c(-0.1, -2.0) # margin 2.9 vs 1.0; bottleneck is 1 for both
  )
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_equal(g$data$tier2[1], g$data$tier2[2])
})

test_that("a variant in the AS 3' window is graver than everything else", {
  d <- pair_table(discrimination = c(0.95, 0.60), as_n_3p = c(1L, 0L))
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  pick <- select_with_ladder(g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties)
  expect_equal(pick$best$as_n_3p, 0L)
  expect_length(pick$relaxed, 0L)
})

test_that("weighted load ranks one common variant below three rare ones", {
  d <- pair_table(discrimination = c(0.9, 0.9), n_var = c(1L, 3L))
  d$w_var <- c(0.40, 0.03) # one common SNP is the graver load here
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_equal(g$data$m_footprint, -c(0.40, 0.03))
  # the gate is still the plain count: both rows have variants under the footprint
  expect_false(any(g$data$gate_com_footprint))
})
