# =============================================================================
#  SSR stutter: the Woerner 2017 slope model (FORENSIC HUMAN STR ONLY)
# =============================================================================
#  Source: PMC5704242 = Genes 2017;8(11):329, DOI 10.3390/genes8110329.
#  NOT the closed Vilsen 2018 FSI Genetics paper -- see the generator header.
#  These tests hold the scope flag, and hold the numbers to what the paper
#  prints in Tables 3 and 4.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("stutter: the file is a per-locus slope model with both measures", {
  d <- ext("stutter_woerner2017.csv")
  expect_equal(nrow(d), 20L) # 10 loci x 2 measures
  expect_equal(length(unique(d$locus)), 10L)
  expect_setequal(unique(d$measure), c("PAL", "LUS"))
  # the two measures cover the same loci, so a slope is comparable
  expect_setequal(d$locus[d$measure == "PAL"], d$locus[d$measure == "LUS"])
})

test_that("stutter: slopes and intercepts stay in a physically possible range", {
  d <- ext("stutter_woerner2017.csv")
  # stutter ratio gained per extra repeat unit: positive, and well under 1
  expect_true(all(d$slope_hap1 > 0 & d$slope_hap1 < 1))
  # stutter does not start at a large fraction of the true allele
  expect_true(all(abs(d$intercept_hap1) < 1))
  # p-values are probabilities.  p_same_slope is printed as "--" for every
  # locus except D5S818 (the paper selects the non-interaction model there),
  # so NA is the correct encoding of "no test was run" and is excluded.
  for (col in c("p_same_intercept", "p_correct")) {
    expect_true(all(d[[col]] >= 0 & d[[col]] <= 1), label = col)
  }
  expect_true(all(is.na(d$p_same_slope) == (d$locus != "D5S818")))
  ps <- d$p_same_slope[!is.na(d$p_same_slope)]
  expect_true(all(ps >= 0 & ps <= 1))
  # R2 is a proportion
  expect_true(all(d$r2 >= 0 & d$r2 <= 1))
})

test_that("stutter: D5S818 is the only locus the paper gives two slopes", {
  d <- ext("stutter_woerner2017.csv")
  two <- unique(d$locus[!is.na(d$slope_hap2)])
  expect_equal(two, "D5S818")
  ds <- d[d$locus == "D5S818", ]
  expect_equal(nrow(ds), 2L)
  expect_true(all(ds$slope_hap2 == 0.009))
})

test_that("stutter: the paper's own worked example reproduces", {
  d <- ext("stutter_woerner2017.csv")
  i <- which(d$locus == "D9S1122" & d$measure == "PAL")
  expect_length(i, 1L)
  gap <- abs(d$intercept_hap1[i] - d$intercept_hap2[i])
  expect_equal(round(gap, 3L), 0.04)
  # "roughly equals the effect of having an additional three repeat units"
  expect_lt(abs(3 * d$slope_hap1[i] - gap), 0.002)
})

test_that("stutter: the scope flag says HUMAN STR, not plant SSR", {
  d <- ext("stutter_woerner2017.csv")
  expect_true(all(grepl("HUMAN STR", d$scope_flag, fixed = TRUE)))
  expect_true(all(grepl("not plant SSR", d$scope_flag, fixed = TRUE)))
  expect_false(any(grepl("plant SSR behaviour by measurement", d$scope_flag, fixed = TRUE)))
  # the file must not silently claim a plant measurement anywhere
  expect_false(any(grepl("\\bplant\\b.*measured|measured.*\\bplant\\b", d$outcome_definition)))
})
