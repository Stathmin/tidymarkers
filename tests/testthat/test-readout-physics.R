# Readout-physics extract tables: HRM (ijms2022, jsbbs2018), ddPCR limits,
# Sanger detection floor.
#
# These tests read the CSVs that data-raw/*.R writes. They deliberately do not
# re-parse .tmp/lit (that is the generator's job, and .tmp is untracked); when
# the source tree IS present the last test re-runs every generator and requires
# it to succeed, which is where the "stop() on disagreement" guarantee is
# exercised end to end.

extdata <- function(...) tm_extdata(...)
read_ext <- function(name) utils::read.csv(extdata(name), stringsAsFactors = FALSE)

test_that("ijms2022: every control is accounted for in each partition", {
  part <- read_ext("hrm_ijms2022_partition_2snp.csv")
  # class_size x n_classes must reproduce the 36 controls exactly
  expect_equal(sum(part$class_size * part$n_classes), 36L)
  expect_equal(unique(part$n_controls_total), 36L)
  # the number of classes is the sum over classes present
  expect_equal(sum(part$n_classes), 1L + 6L + 2L + 3L + 1L)
  expect_equal(sum(part$class_size == 1L), 1L) # exactly one singleton

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
  expect_true(all(
    vapply(by_n[-1L], function(z) z[[1L]], numeric(1)) >=
      vapply(by_n[-length(by_n)], function(z) z[[2L]], numeric(1))
  ))
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
  for (f in c(
    "dpcr_cpn_lob_variability.csv", "dpcr_aml_idh_limits.csv",
    "dpcr_kit_d816v_limits.csv"
  )) {
    expect_false(any(grepl("conversion|sensitivity_rate", names(read_ext(f)))))
  }
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
  # .tmp/lit may exist for other projects, so the guard must name the exact
  # fetched files these six generators read, not just the directory.
  lit <- testthat::test_path("..", "..", ".tmp", "lit")
  required <- c(
    "ijms2022/fulltext.xml", "ijms2022/supp_tables.txt",
    "jsbbs2018/fulltext.xml", "jsbbs2018/supplements/68_18048_2.xlsx",
    "PMC12610676/fulltext.xml", "PMC13410330/fulltext.xml",
    "PMC13584697/fulltext.xml", "PMC12188199/fulltext.xml",
    "PMC13410318/fulltext.xml", "mosha2014/fulltext.xml"
  )
  skip_if(
    !all(file.exists(file.path(lit, required))),
    "source tree .tmp/lit not fetched in this checkout"
  )
  gens <- c(
    "ijms2022.R", "jsbbs2018.R", "aml_idh_ddpcr.R", "kit_d816v_ddpcr.R",
    "cpn_ddpcr_lob.R", "sanger_floor.R"
  )
  repo <- testthat::test_path("..", "..")
  for (g in gens) {
    out <- suppressMessages(system2(
      file.path(R.home("bin"), "Rscript"),
      c("--vanilla", file.path(repo, "data-raw", g)),
      stdout = TRUE, stderr = TRUE
    ))
    expect_null(attr(out, "status"), info = paste(g, paste(out, collapse = "\n")))
  }
})

## --- three readout tables that were on disk with no consumer -----------------

test_that("JAK2 V617F precision: the CV reproduces above 1 % and does NOT at 0.1 %", {
  d <- read_ext("jak2v617f_ddpcr_precision.csv")
  expect_equal(nrow(d), 16L)
  expect_equal(length(unique(d$target_vaf_pct)), 4L)
  expect_equal(length(unique(d$condition)), 4L)
  # CV = sd / mean, in percent: recomputed from the printed mean and sd
  d$cv_calc <- 100 * d$sd_pct / d$mean_vaf_pct
  hi <- d$target_vaf_pct >= 1
  expect_true(all(abs(d$cv_calc[hi] - d$cv_pct[hi]) < 0.05))
  # at 0.1 % the paper's cells are rounded hard enough that the identity fails.
  # The extract keeps the printed values, so the test asserts the FAILURE: it is
  # what stops anyone calibrating a LoD on a recomputed number.
  lo <- d$target_vaf_pct == 0.1
  expect_true(all(abs(d$cv_calc[lo] - d$cv_pct[lo]) > 0.05))
})

test_that("ijms2022 minimum input: the smallest amount is not a genotype call", {
  d <- read_ext("hrm_ijms2022_min_input.csv")
  expect_equal(nrow(d), 5L)
  expect_equal(d$min_ng, c(3, 4, 1, 0.1, 0.01))
  expect_true(all(d$evidence_tier == "measured"))
  # the two smallest rows are the limits of the method, not genotypes
  expect_true(all(grepl("NOT a genotype call", d$outcome[4:5])))
  expect_true(all(grepl("correctly identified", d$outcome[1:3])))
  expect_setequal(unique(d$method), c("HRM", "PSQ", "HRM (single-strand only)", "HRM (no product)"))
})

## --- HRM melt: the diatom Tm set, one row per amplicon ----------------------

test_that("HRM diatom melt: 28 amplicons, one run, and the paper's own slips", {
  d <- read_ext("readout_hrm_diatom2025.csv")
  expect_equal(nrow(d), 28L)
  expect_setequal(
    d$primer_pair,
    c("UPA99", "18S-V7", "ND", "Primer1", "psaA-2", "Dia18S")
  )
  n_by <- table(d$primer_pair)
  expect_equal(as.integer(n_by[["Dia18S"]]), 3L)
  expect_true(all(as.integer(n_by[names(n_by) != "Dia18S"]) == 5L))
  # constant columns are provenance, not per-row choices
  expect_setequal(unique(d$chemistry), "hrm_melt")
  expect_equal(length(unique(d$assay)), 1L)
  expect_true(all(d$evidence_tier == "measured"))
  expect_true(all(d$doi == "10.1371/journal.pone.0321885"))
  expect_true(all(d$pmcid == "PMC12002481"))
  expect_true(all(nzchar(d$source_file)))
  # the CSV must keep saying what it is: absolute Tm of one amplicon, not a dTm
  expect_true(all(grepl("NOT a dTm", d$outcome_definition)))
})

test_that("HRM diatom melt: thermodynamics scale with length and stay negative", {
  d <- read_ext("readout_hrm_diatom2025.csv")
  expect_true(all(d$tm_measured_c > 70 & d$tm_measured_c < 90))
  expect_true(all(d$dh_kcal_mol < 0))
  expect_true(all(d$ds_cal_k_mol < 0))
  expect_true(all(d$dh_kcal_mol / d$amplicon_bp > -8.5))
  expect_true(all(d$dh_kcal_mol / d$amplicon_bp < -7.8))
  expect_true(all(d$ds_cal_k_mol / d$amplicon_bp > -22.5))
  expect_true(all(d$ds_cal_k_mol / d$amplicon_bp < -21.0))
  # the printed extremes, which are the paper's own (its prose claims 60-152 bp)
  expect_equal(range(d$amplicon_bp), c(59L, 154L))
  expect_equal(range(d$gc_pct), c(34.62, 55.24))
})

test_that("HRM diatom melt: Table 7 vs Table 8 disagreements are kept, not patched", {
  d <- read_ext("readout_hrm_diatom2025.csv")
  slip <- d$gc_printed_vs_recomputed_slip
  expect_equal(sum(slip), 3L)
  # the paper's three arithmetic slips: printed GC% HIGHER than recomputed
  expect_true(all(d$gc_pct[slip] > d$gc_from_seq[slip]))
  expect_true(all(abs(d$gc_pct[slip] - d$gc_from_seq[slip]) >= 0.06))
  ok <- !slip & !is.na(d$gc_from_seq)
  expect_true(all(abs(d$gc_pct[ok] - d$gc_from_seq[ok]) < 0.06))
  # no sequence is printed for the three Dia18S cells, so no GC is recomputed
  expect_true(all(is.na(d$gc_from_seq) == (d$primer_pair == "Dia18S")))
  # the +1 nt between Table 7 length and Table 8 bp, kept as a column
  expect_equal(d$seq_len - d$amplicon_bp, d$seq_len_vs_table_bp)
  expect_equal(sum(d$seq_len_vs_table_bp == 1L, na.rm = TRUE), 23L)
  expect_equal(sum(d$seq_len_vs_table_bp == 0L, na.rm = TRUE), 2L)
  expect_equal(sum(is.na(d$seq_len_vs_table_bp)), 3L)
})

test_that("HRM diatom melt: no salt, so the set cannot constrain a Tm model alone", {
  d <- read_ext("readout_hrm_diatom2025.csv")
  # Mg2+ sits inside a proprietary 2x master mix; the paper prints neither it
  # nor the DNA input, so a salt term cannot be fitted. This is why
  # R/marker_readout.R keeps HRM at simulated = FALSE.
  expect_true(all(is.na(d$mg_mm)))
  expect_true(all(is.na(d$dna_ng_ul)))
  # Absolute Tm of DIFFERENT amplicons: within each primer pair the products
  # span < 3 C, and the six pair means span > 5 C. There is no dTm of one
  # locus across genotypes in this set.
  spread <- tapply(d$tm_measured_c, d$primer_pair, function(z) diff(range(z)))
  expect_true(all(spread < 3))
  means <- sort(tapply(d$tm_measured_c, d$primer_pair, mean))
  expect_gt(means[[length(means)]] - means[[1L]], 5)
})

test_that("Sanger floor: a textbook number may not calibrate anything", {
  d <- read_ext("sanger_detection_floor.csv")
  expect_equal(nrow(d), 2L)
  expect_true(all(d$metric == "detection_floor_pct"))
  expect_true(all(d$value_pct == 20))
  # both rows are restatements, so neither is a measurement by these authors
  expect_true(all(d$evidence_tier == "textbook"))
  expect_true(all(!d$usable_as_calibration))
  expect_true(all(nzchar(d$reason)))
})
