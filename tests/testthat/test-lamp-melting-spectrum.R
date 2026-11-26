# =============================================================================
#  Z8/ZY17  --  the quenching-probe LAMP melting table of Imai et al. 2025
# =============================================================================
#  Source key lamp_melting_allele_spectrum2025, DOI 10.1128/spectrum.01565-25
#  (PMC12671098), is the only source in the corpus with a per-sample LAMP
#  melting readout.
#  This test is the consumer gate: it pins the three probe partitions, the
#  allele -> call mapping, the samples the paper could not resolve, and the
#  reaction chemistry that gives `SCHEME_UNMODELLED$lamp = "betaine"` a number.
#
#  `na = character()` is deliberate: the supplement writes the string "NA" for
#  an unassembled call and "-" for an untested dipstick, and both must survive
#  as text.  Every "NA" call still carries a numeric Tm.

calls <- function() {
  readr::read_csv(tm_extdata("lamp_tp0705_melting_calls.csv"),
    show_col_types = FALSE,
    na = character()
  )
}

conditions <- function() {
  utils::read.csv(tm_extdata("lamp_tp0705_reaction_conditions.csv"),
    stringsAsFactors = FALSE
  )
}

test_that("lamp melting: the mismatch family table states what the endpoint carries", {
  m <- readr::read_csv(tm_extdata("mismatch_family_lamp.csv"),
    show_col_types = FALSE,
    na = c("", "NA")
  )
  expect_equal(nrow(m), 3L)
  expect_setequal(m$family, "lamp")
  expect_setequal(m$doi, "10.1128/spectrum.01565-25")
  expect_setequal(m$pmcid, "PMC12671098")
  # the two lineage-specific peaks are 1.25 C apart, just outside the paper's
  # own +/-1 C assignment margin
  v <- function(metric) m$value[m$metric == metric]
  expect_equal(v("melting_tm_c"), 92.46)
  expect_equal(v("delta_melting_tm_c"), 1.25)
  expect_equal(v("tm_shift_direction"), -1)
  # the table carries a SIGN where the paper prints no per-mismatch magnitude,
  # so n_mismatches is genuinely unset rather than zero
  expect_true(all(is.na(m$n_mismatches)))
  expect_true(all(c("melting_tm_c", "delta_melting_tm_c", "tm_shift_direction") %in% m$metric))
})

test_that("lamp melting: the sample table keeps all 133 typed records", {
  d <- calls()
  expect_equal(nrow(d), 133L)
  expect_equal(ncol(d), 13L)
  expect_false(anyDuplicated(d$sample_no) > 0L)
  expect_setequal(
    names(d),
    c(
      "sample_no", "sample_type", "strain", "lineage", "tp0705_allele",
      "tm_506", "call_506", "tm_625", "call_625", "tm_708", "call_708",
      "tp0548_nichols_lamp", "tp0548_ss14_lamp"
    )
  )
  expect_equal(sum(d$sample_type == "Swab"), 105L)
  expect_equal(sum(d$sample_type == "Saliva"), 28L)
  expect_equal(sum(d$lineage == "SS14-East Asia"), 110L)
  expect_equal(sum(d$lineage == "Nichols"), 9L)
  expect_equal(sum(d$lineage == "SS14\u03a9-A"), 7L)
  expect_equal(sum(d$lineage == "SS14\u03a9-B 1A"), 2L)
})

test_that("lamp melting: each probe splits into the paper's alleles", {
  d <- calls()
  expect_equal(
    c(sum(d$call_506 == "V"), sum(d$call_506 == "T"), sum(d$call_506 == "A")),
    c(106L, 11L, 12L)
  )
  expect_equal(c(sum(d$call_625 == "M"), sum(d$call_625 == "V")), c(123L, 7L))
  expect_equal(c(sum(d$call_708 == "G"), sum(d$call_708 == "S")), c(117L, 12L))

  # the peaks are disjoint: that is the whole readout
  expect_equal(range(d$tm_506[d$call_506 == "V"]), c(60.08, 61.52))
  expect_equal(range(d$tm_506[d$call_506 == "T"]), c(63.64, 64.03))
  expect_equal(range(d$tm_506[d$call_506 == "A"]), c(71.52, 71.91))
  expect_equal(range(d$tm_625[d$call_625 == "M"]), c(63.09, 63.86))
  expect_equal(range(d$tm_625[d$call_625 == "V"]), c(69.96, 70.94))
  expect_equal(range(d$tm_708[d$call_708 == "G"]), c(72.94, 74.74))
  expect_equal(range(d$tm_708[d$call_708 == "S"]), c(67.38, 68.00))

  # every gap a caller would have to resolve is wider than the paper's own
  # +/-1 C decision margin
  expect_gt(min(d$tm_506[d$call_506 == "T"]) - max(d$tm_506[d$call_506 == "V"]), 2)
  expect_gt(min(d$tm_506[d$call_506 == "A"]) - max(d$tm_506[d$call_506 == "T"]), 2)
  expect_gt(min(d$tm_625[d$call_625 == "V"]) - max(d$tm_625[d$call_625 == "M"]), 2)
  expect_gt(max(d$tm_708[d$call_708 == "G"]) - min(d$tm_708[d$call_708 == "S"]), 2)
})

test_that("lamp melting: the call follows the tp0705 allele", {
  d <- calls()
  typed <- d[d$call_506 != "NA", ]
  expect_true(all(typed$call_506[typed$tp0705_allele == "8"] == "V"))
  expect_true(all(typed$call_506[typed$tp0705_allele == "3"] == "T"))
  expect_true(all(typed$call_506[typed$tp0705_allele %in% c("1", "9")] == "A"))
  expect_equal(sum(d$tp0705_allele == "8"), 110L)
  expect_equal(sum(d$tp0705_allele == "3"), 11L)
  expect_equal(sum(d$tp0705_allele == "1"), 7L)
  expect_equal(sum(d$tp0705_allele == "9"), 5L)
})

test_that("lamp melting: unresolved calls stay verbatim and keep a Tm", {
  d <- calls()
  expect_equal(sum(d$call_506 == "NA"), 4L)
  expect_equal(sum(d$call_625 == "NA"), 3L)
  expect_equal(sum(d$call_708 == "NA"), 4L)
  expect_equal(d$sample_no[d$call_506 == "NA"], c(112L, 113L, 123L, 128L))
  expect_equal(d$sample_no[d$call_625 == "NA"], c(112L, 123L, 128L))
  expect_equal(d$sample_no[d$call_708 == "NA"], c(112L, 113L, 123L, 128L))
  expect_false(anyNA(d$tm_506))
  expect_false(anyNA(d$tm_625))
  expect_false(anyNA(d$tm_708))
  # the four 506 "NA" samples are all saliva, all allele 8 -- the paper's note
  expect_true(all(d$sample_type[d$call_506 == "NA"] == "Saliva"))
  expect_true(all(d$tp0705_allele[d$call_506 == "NA"] == "8"))
})

test_that("lamp melting: the tp0548 dipstick table keeps its partition", {
  d <- calls()
  expect_equal(sum(d$tp0548_nichols_lamp == "-" & d$tp0548_ss14_lamp == "-"), 122L)
  expect_equal(
    sum(d$tp0548_nichols_lamp == "Positive" & d$tp0548_ss14_lamp == "Negative"),
    9L
  )
  expect_equal(
    sum(d$tp0548_nichols_lamp == "Negative" & d$tp0548_ss14_lamp == "Positive"),
    2L
  )
})

test_that("lamp melting: the conditions file names the betaine magnitude", {
  cnd <- conditions()
  expect_equal(nrow(cnd), 21L)
  expect_setequal(
    names(cnd),
    c(
      "item", "value", "unit", "source_section", "quote",
      "outcome_definition", "evidence_tier", "source_file", "doi"
    )
  )
  expect_true(all(cnd$evidence_tier == "weak"))
  expect_true(all(nzchar(cnd$quote)))
  expect_true(all(cnd$doi == "10.1128/spectrum.01565-25"))

  value_of <- function(item) cnd$value[cnd$item == item]
  expect_equal(value_of("betaine"), "0.2")
  expect_equal(cnd$unit[cnd$item == "betaine"], "M")
  expect_equal(value_of("MgSO4"), "6")
  expect_equal(value_of("amplification temperature"), "68")
  expect_equal(value_of("QProbe"), "80")
  expect_equal(value_of("melt ramp"), "0.2")
  expect_equal(value_of("Tm decision margin"), "+/-1")
})
