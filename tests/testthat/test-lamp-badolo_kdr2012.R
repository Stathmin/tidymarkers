# AS-LAMP kdr-w (Malar J 2012;11:227): the second time window, and the source
# that pins the B2 terminal dinucleotide. Data-only, from inst/extdata.

bd_read <- function(name) {
  utils::read.csv(testthat::test_path("..", "..", "inst", "extdata", name),
                  stringsAsFactors = FALSE)
}

test_that("the five badolo_kdr2012 files exist with required columns", {
  for (f in c("lamp_badolo_kdr2012_primers.csv",
              "lamp_badolo_kdr2012_time_to_positive.csv",
              "lamp_badolo_kdr2012_counts.csv",
              "lamp_badolo_kdr2012_metrics.csv",
              "lamp_badolo_kdr2012_contradictions.csv")) {
    d <- bd_read(f)
    expect_gt(nrow(d), 0L)
    expect_true(all(c("outcome_definition", "evidence_tier", "source_file") %in% names(d)))
  }
})

test_that("the DOI landed is the AS-LAMP paper, not the monitoring survey", {
  d <- bd_read("lamp_badolo_kdr2012_time_to_positive.csv")
  expect_identical(unique(d$doi), "10.1186/1475-2875-11-227")
  expect_false(any(d$doi == "10.1186/1475-2875-11-232"))
})

test_that("the two BIP arms differ only at the B2 3' terminal base", {
  p <- bd_read("lamp_badolo_kdr2012_primers.csv")
  expect_identical(p$primer, c("F3", "B3", "FIP", "BIP-wild", "BIP-kdr-w"))
  w <- p$sequence_as_printed[p$primer == "BIP-wild"]
  k <- p$sequence_as_printed[p$primer == "BIP-kdr-w"]
  expect_identical(nchar(w), nchar(k))
  expect_identical(sub("[a-zA-Z]$", "", w), sub("[a-zA-Z]$", "", k))
  expect_identical(p$terminal_bases[p$primer == "BIP-wild"], "gT")
  expect_identical(p$terminal_bases[p$primer == "BIP-kdr-w"], "gA")
  # F3/B3/FIP are shared between the sets, so discrimination rests on BIP
  expect_true(all(grepl("^[ACGT]+$", p$sequence_as_printed[1:2])))
})

test_that("Table 1 sums to 120 for every method", {
  c <- bd_read("lamp_badolo_kdr2012_counts.csv")
  for (m in unique(c$method)) {
    d <- c[c$method == m, ]
    expect_identical(sum(d$n_called), 120L)
  }
  seq_row <- c[c$method == "sequencing (truth)" & c$genotype != "unidentified", ]
  expect_identical(sort(seq_row$n_called), c(20L, 26L, 70L))
})

test_that("printed metrics are verbatim and flagged when not reproducible", {
  m <- bd_read("lamp_badolo_kdr2012_metrics.csv")
  expect_identical(nrow(m), 6L)
  expect_true(all(m$printed_sensitivity >= 0 & m$printed_sensitivity <= 1))
  expect_true(all(m$printed_specificity >= 0 & m$printed_specificity <= 1))
  expect_true(all(is.logical(m$sens_reproduced)))
  # AS-LAMP called 23 kdr-w/kdr-w of 20 sequenced: the printed 0.92 cannot be
  # called/sequenced, and the row says so rather than hiding it
  r <- m[m$method == "AS-LAMP" & m$genotype == "kdr-w/kdr-w", ]
  expect_identical(r$printed_sensitivity, 0.92)
  expect_gt(r$recomputed_sens_called_over_truth, 1)
  expect_false(r$sens_reproduced)
})

test_that("non-target arms are censored at the 75-min window, never timed", {
  tp <- bd_read("lamp_badolo_kdr2012_time_to_positive.csv")
  expect_identical(nrow(tp), 10L)
  expect_true(all(tp$observation_window_min == 75L))
  expect_true(all(is.na(tp$time_to_positive_min[!tp$is_target])))
  expect_true(all(grepl("not amplified within 75 min",
                        tp$censoring[!tp$is_target])))
  expect_true(all(tp$time_to_positive_min <= tp$observation_window_min, na.rm = TRUE))
  # the heterozygote is the slowest DETECTED arm of both sets: it sits ON the cut-off
  het <- tp$time_to_positive_min[tp$template == "gDNA kdr-w/+"]
  expect_equal(het, c(75, 65))
})

test_that("the penultimate vs 3'-end wording conflict is recorded", {
  k <- bd_read("lamp_badolo_kdr2012_contradictions.csv")
  expect_true(any(grepl("PENULTIMATE", k$a_version)))
  expect_true(any(grepl("3' end of the BIP", k$b_version)))
  expect_true(all(grepl("recorded", k$resolution)))
})
