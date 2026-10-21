# Readout-physics extract tables: HRM (ijms2022, jsbbs2018), ddPCR limits,
# Sanger detection floor.
#
# These tests read the CSVs that data-raw/*.R writes. They deliberately do not
# re-parse .tmp/lit (that is the generator's job, and .tmp is untracked); when
# the source tree IS present the last test re-runs every generator and requires
# it to succeed, which is where the "stop() on disagreement" guarantee is
# exercised end to end.

extdata <- function(...) testthat::test_path("..", "..", "inst", "extdata", ...)
read_ext <- function(name) utils::read.csv(extdata(name), stringsAsFactors = FALSE)

test_that("ijms2022: every control is accounted for in each partition", {
  part <- read_ext("hrm_ijms2022_partition_2snp.csv")
  # class_size x n_classes must reproduce the 36 controls exactly
  expect_equal(sum(part$class_size * part$n_classes), 36L)
  expect_equal(unique(part$n_controls_total), 36L)
  # the number of classes is the sum over classes present
  expect_equal(sum(part$n_classes), 1L + 6L + 2L + 3L + 1L)
  expect_equal(sum(part$class_size == 1L), 1L)   # exactly one singleton

  unr <- read_ext("hrm_ijms2022_unresolvable_3snp.csv")
  pairs <- unr[unr$melt_class == "pair (unresolvable)", ]
  expect_equal(nrow(pairs), 7L)
  expect_equal(2L * nrow(pairs), 14L)
  # 36 - 14 in pairs = 22 uniquely identified in the 3-SNP assay
  expect_setequal(na.omit(unique(unr$n_uniquely_identified)), 36L - 14L)
})

test_that("ijms2022: discrimination counts rise with panel size and never fall", {
  d <- read_ext("hrm_ijms2022_discrimination.csv")
  for (rd in unique(d$readout)) {
    v <- d[d$readout == rd, ]
    v <- v[order(v$assay), ]
    expect_true(all(v$n_uniquely_identified >= 1L))
    expect_true(all(v$n_uniquely_identified + v$n_not_discriminated == v$n_controls))
    expect_gte(v$n_uniquely_identified[[2L]], v$n_uniquely_identified[[1L]])
  }
  # HRM+PSQ >= each readout alone, and <= the number of controls
  expect_equal(d$n_uniquely_identified[d$readout == "HRM+PSQ"], c(13L, 34L))
})

test_that("ijms2022: dTm is physical and rises with the number of T alleles", {
  dtm <- read_ext("hrm_ijms2022_dtm_single_snp_controls.csv")
  expect_true(all(dtm$dtm_c >= 0 & dtm$dtm_c < 5))
  expect_equal(dtm$dtm_c[dtm$n_t_alleles == 0L], 0)
  # ordering between loci with the SAME T count is not monotone (0.4/0.7/0.6
  # for one T at three different loci) - the neighbours differ, so only the
  # between-count relation is asserted.
  by_n <- tapply(dtm$dtm_c, dtm$n_t_alleles, range)
  expect_true(all(vapply(by_n[-1L], function(z) z[[1L]], numeric(1)) >=
                    vapply(by_n[-length(by_n)], function(z) z[[2L]], numeric(1))))
})

test_that("jsbbs2018: dTm distributions are physical and conditioned", {
  e1 <- read_ext("hrm_jsbbs2018_dtm_exp1.csv")
  e3 <- read_ext("hrm_jsbbs2018_dtm_exp3.csv")
  fun <- read_ext("hrm_jsbbs2018_funnel.csv")
  rule <- read_ext("hrm_jsbbs2018_selection_rule.csv")
  expect_equal(nrow(e1), 38L)
  expect_equal(nrow(e3), 18L)
  expect_true(all(e1$dtm_obs > 0 & e1$dtm_obs < 5))
  expect_equal(e1$dtm_obs, abs(e1$tm_aokimame - e1$tm_fukuyutaka), tolerance = 1e-9)
  expect_true(all(e3$product_bp >= 54 & e3$product_bp <= 100))
  # exactly one rule is quoted in kcal/mol: the free-energy proxy. A dTm
  # threshold in degrees C would be a different (and unsourced) claim.
  expect_equal(sum(rule$units == "kcal/mol", na.rm = TRUE), 1L)
  expect_false(any(rule$units == "C", na.rm = TRUE))
  expect_true(all(fun$n_success <= fun$n_designed))
  expect_equal(fun$n_designed, c(48L, 32L, 24L))
})

test_that("ddPCR: LoB <= LoD <= LoQ holds within every assay", {
  kit <- read_ext("dpcr_kit_d816v_limits.csv")
  for (nm in unique(kit$normalisation)) {
    v <- kit[kit$normalisation == nm, ]
    expect_lt(v$value_pct[v$limit_type == "LoB"], v$value_pct[v$limit_type == "LoD"] + 1e-12)
    expect_lte(v$value_pct[v$limit_type == "LoD"], v$value_pct[v$limit_type == "LoQ"] + 1e-12)
  }
  aml <- read_ext("dpcr_aml_idh_limits.csv")
  expect_setequal(unique(aml$limit_type), c("LoB", "LoD"))
  expect_equal(sum(aml$limit_type == "LoB"), 4L)
  expect_true(all(aml$lod_fa <= 1, na.rm = TRUE))
})

test_that("ddPCR: LoB varies across conditions and no conversion rate ships", {
  cpn <- read_ext("dpcr_cpn_lob_variability.csv")
  lob <- cpn$cp_ul[cpn$limit_type == "LoB"]
  expect_equal(length(lob), 5L)
  expect_gt(max(lob) / min(lob), 2)
  for (f in c("dpcr_cpn_lob_variability.csv", "dpcr_aml_idh_limits.csv",
              "dpcr_kit_d816v_limits.csv"))
    expect_false(any(grepl("conversion|sensitivity_rate", names(read_ext(f)))))
})

test_that("Sanger floor: sourced, bounded, and not usable as calibration", {
  sf <- read_ext("sanger_detection_floor.csv")
  expect_gt(nrow(sf), 0L)
  expect_true(all(sf$value_pct >= 5 & sf$value_pct <= 30))
  expect_true(all(sf$usable_as_calibration == (sf$evidence_tier == "measured")))
  expect_false(any(sf$usable_as_calibration))
  expect_true(all(nzchar(sf$source_file)))
})

test_that("generators reproduce the shipped extracts when sources are present", {
  skip_if(!dir.exists(testthat::test_path("..", "..", ".tmp", "lit")),
          "source tree .tmp/lit not fetched in this checkout")
  gens <- c("ijms2022.R", "jsbbs2018.R", "aml_idh_ddpcr.R", "kit_d816v_ddpcr.R",
            "cpn_ddpcr_lob.R", "sanger_floor.R")
  repo <- testthat::test_path("..", "..")
  for (g in gens) {
    out <- suppressMessages(system2(
      file.path(R.home("bin"), "Rscript"), c("--vanilla", file.path(repo, "data-raw", g)),
      stdout = TRUE, stderr = TRUE))
    expect_null(attr(out, "status"), info = paste(g, paste(out, collapse = "\n")))
  }
})
