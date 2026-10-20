# A composed technology (roadmap Y1): a base chemistry + designer params +
# a specificity policy, first-class wherever a scheme is.

test_that("a composed technology is discovered with a valid DAG and its policy", {
  marker_technology("giga_snp",
    base = "snp", params = list(as_lengths = 24:30),
    specificity = "unique", label = "giga SNP KASP"
  )
  on.exit(rm("giga_snp", envir = .marker_tech_store), add = TRUE)

  expect_true("giga_snp" %in% marker_schemes()$scheme)
  expect_identical(marker_scheme("giga_snp")$family, "kasp")
  # the specificity policy is carried (Y3) and visible to the DAG contract
  expect_identical(marker_specificity("giga_snp"), "unique")
  expect_dag_invariants("giga_snp")

  t <- marker_technologies()
  expect_true("giga_snp" %in% t$scheme)
  expect_identical(t$base[t$scheme == "giga_snp"], "snp")
  expect_identical(t$specificity[t$scheme == "giga_snp"], "unique")
})

test_that("a composed technology designs via its base designer", {
  skip_if_no_oligoarrayaux()
  marker_technology("giga_test", base = "snp", params = list())
  on.exit(rm("giga_test", envir = .marker_tech_store), add = TRUE)
  f <- snp_fixture()
  genome <- f$genome
  snp <- f$snp
  a <- suppressWarnings(marker_design_at(genome, snp, "giga_test"))
  expect_false(is.null(a))
})

test_that("marker_technology validates its arguments", {
  expect_error(marker_technology("", base = "snp"), "non-empty")
  expect_error(marker_technology("x", base = "nope"), "Unknown marker scheme")
  expect_error(marker_technology("snp"), "already exists")
  expect_error(marker_technology("y", base = "snp", specificity = "wat"), "Unknown specificity")
  expect_error(marker_technology("z", base = "snp", params = list(1, 2)), "named list")
})

test_that("marker_filter_blocks names the composable space", {
  b <- marker_filter_blocks()
  expect_true(all(c("base", "specificity", "params") %in% b$block))
  expect_setequal(b$values[b$block == "specificity"], paste(KASP_SPECIFICITY_MODES, collapse = " | "))
})
