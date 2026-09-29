# KASP reaction recipe: the AS:AS:COM = 2:2:5 ratio and volume bookkeeping.

test_that("the primer mix keeps the AS:AS:COM = 2:2:5 ratio", {
  r <- kasp_mix_recipe(mix_uL = 100)
  pm <- dplyr::filter(r, .data$stage == "primer mix")
  fam <- pm$final_uM[pm$component == "AS1_FAM"]
  hex <- pm$final_uM[pm$component == "AS2_HEX"]
  com <- pm$final_uM[pm$component == "COM"]
  expect_equal(fam, hex)
  expect_equal(com / fam, 2.5)                       # 5 / 2
  # stock volumes are in the same 2:2:5 proportion
  mix <- dplyr::filter(pm, .data$component != "nuclease-free water")
  expect_equal(mix$total_uL, c(12, 12, 30))          # uL for a 100 uL mix
  expect_equal(sum(pm$total_uL), 100)                # stocks + water fill the mix
})

test_that("the reaction volumes add up to the requested total", {
  r <- kasp_mix_recipe(n_reactions = 96, reaction_uL = 10)
  rx <- dplyr::filter(r, .data$stage == "reaction")
  expect_equal(sum(rx$per_rxn_uL), 10)
  total_rxn <- ceiling(96 * 1.1)
  expect_equal(rx$total_uL, rx$per_rxn_uL * total_rxn)
})

test_that("a too-concentrated working mix is rejected", {
  expect_error(kasp_mix_recipe(stock_uM = 40, as_uM = 12, com_uM = 30), "stock")
})
