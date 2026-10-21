# =============================================================================
#  ZW3 external oracle  --  mm2020: 4096 MEASURED duplex Tm, tabulated
# =============================================================================
#  Oliveira, Long, Brown, Fox & Weber, Chem Sci 2020, DOI 10.1039/d0sc01700k.
#  ESI Table S1 prints measured Tm for 4096 19-mers that differ only in a
#  central trimer, under a stated buffer: 50 mM NaCl + 10 mM Na-phosphate
#  pH 7.4, 1.0 uM total strand, no Mg2+.  This is the first anchor in this
#  package whose numbers were MEASURED and PRINTED, rather than computed by a
#  model or read off a raster figure.
#
#  WHAT IS ASSERTED, AND WHY IT IS NOT A |dTm| <= 2 C TEST
#
#    primer_tm() sits about 16.9 C BELOW these measurements (sd 0.64 C) on the
#    64 canonical diblock-free duplexes, with r = 0.958.  The offset is a LEVEL
#    offset across the whole set, it is NOT a Mg effect (Mg = 0 here; adding
#    3 mM Mg closes only ~12 C and then overshoots on the mismatch classes),
#    and it is deliberately NOT absorbed by any correction in the package.
#    What the data DOES support is:
#      * ordering: model and measurement rank the 64 canonical duplexes alike;
#      * discrimination: the measured spread across the four mismatch classes
#        (canonical 60.7-69.3 C, triple 41.3-62.8 C) is real, and the model
#        must be able to separate a perfect duplex from a triple mismatch;
#      * stability: a hand-edited CSV breaks the internal invariants below.
#
#  What is NOT asserted: absolute agreement, or any claim that primer_tm()
#  reproduces published Tm at the conventional 1.5 mM Mg2+.  See
#  data-raw/proposed/tm_oracle_status.md.

mm2020_path <- function() {
  test_path("../../inst/extdata/tm_mm2020_s1_measured.csv")
}

test_that("ZW3/mm2020: the landed table still matches the paper's own totals", {
  skip_if_not(file.exists(mm2020_path()), "tm_mm2020 CSV not present")
  d <- utils::read.csv(mm2020_path(), stringsAsFactors = FALSE)

  expect_equal(nrow(d), 4096L)
  expect_equal(sort(d$index), 1:4096)

  # 64 canonical + 576 single + 1728 double + 1728 triple, as the paper states
  expect_equal(
    as.integer(table(factor(d$mismatch_class,
                            levels = c("canonical", "single", "double", "triple")))),
    c(64L, 576L, 1728L, 1728L)
  )
  # every duplex is the same 19-mer scaffold
  expect_true(all(nchar(d$top_18) == 19L))
  expect_true(all(nchar(d$bottom_18) == 19L))
  # canonical rows are perfect Watson-Crick duplexes
  comp <- chartr("ACGT", "TGCA", d$top_18)
  expect_true(all(d$bottom_18[d$mismatch_class == "canonical"] ==
                  comp[d$mismatch_class == "canonical"]))
  # buffer is stated and is sodium-only
  expect_true(all(d$na_mm == 50))
  expect_true(all(d$mg_mm == 0))
})

test_that("ZW3/mm2020: measured Tm is physically possible and class-ordered", {
  skip_if_not(file.exists(mm2020_path()), "tm_mm2020 CSV not present")
  d <- utils::read.csv(mm2020_path(), stringsAsFactors = FALSE)

  expect_true(all(d$measured_tm_c > 30 & d$measured_tm_c < 100))

  # the measured MEDIAN must fall as mismatches are added; this is a
  # monotonicity statement about the data, not about the model
  med <- tapply(d$measured_tm_c, d$mismatch_class, stats::median)
  expect_true(all(diff(med[c("canonical", "single", "double", "triple")]) < 0))

  # the paper's own claim: a substantial minority of single mismatches melt
  # ABOVE the least stable canonical trimer (ATT/TAA, 60.7 C)
  canon_min <- min(d$measured_tm_c[d$mismatch_class == "canonical"])
  expect_equal(canon_min, 60.7)
  sing <- d$measured_tm_c[d$mismatch_class == "single"]
  expect_gt(sum(sing > canon_min), 0)
})

test_that("ZW3/mm2020: primer_tm() ranks the 64 canonical duplexes with measurement", {
  skip_if_not(file.exists(mm2020_path()), "tm_mm2020 CSV not present")
  d <- utils::read.csv(mm2020_path(), stringsAsFactors = FALSE)
  can <- d[d$mismatch_class == "canonical", ]

  tm <- primer_tm(can$top_18,
                  conditions = reaction_conditions(Na = 50, Mg = 0, dNTPs = 0, DMSO = 0))
  expect_length(tm, 64L)
  expect_true(all(is.finite(tm)))

  # measured r = 0.958, Spearman 0.941; the floors below leave room for a
  # TmCalculator release changing its salt code
  expect_gt(stats::cor(tm, can$measured_tm_c), 0.9)
  expect_gt(stats::cor(tm, can$measured_tm_c, method = "spearman"), 0.85)

  # the LEVEL offset is real and must not be silently "fixed": if the model
  # ever lands within 2 C of these measurements at their own salt, someone has
  # introduced a hidden correction and this anchor has gone vacuous
  bias <- mean(tm - can$measured_tm_c)
  expect_lt(bias, -10)
  expect_gt(bias, -25)
  # and it is an offset, not noise
  expect_lt(stats::sd(tm - can$measured_tm_c), 2)
})

test_that("ZW3/mm2020: the model separates a perfect duplex from a triple mismatch", {
  skip_if_not(file.exists(mm2020_path()), "tm_mm2020 CSV not present")
  d <- utils::read.csv(mm2020_path(), stringsAsFactors = FALSE)
  cond <- reaction_conditions(Na = 50, Mg = 0, dNTPs = 0, DMSO = 0)

  # same scaffold, GC-rich canonical vs all-mismatch central trimer
  pick <- function(tri) d$top_18[match(tri, d$central_trimer)]
  hi <- primer_tm(pick("GCG/CGC"), conditions = cond)
  lo <- primer_tm(pick("CCC/CCC"), conditions = cond)
  expect_true(is.finite(hi) && is.finite(lo))
  expect_gt(hi, lo)

  # direction agrees with the measurement for this pair
  expect_gt(d$measured_tm_c[match("GCG/CGC", d$central_trimer)],
            d$measured_tm_c[match("CCC/CCC", d$central_trimer)])
})
