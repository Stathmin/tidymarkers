# What the LAMP-SNP lane may claim (data-raw/proposed/lamp_status.decision.md).
# Data-only: these assert the boundary of the claim, not a LAMP model.

ls_read <- function(name) {
  utils::read.csv(tm_extdata(name),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
}

test_that("LAMP has no simulated readout, so it reports no predicted gap", {
  expect_match(readout_assumption("lamp"), "not modelled")
  # decision clause 1: no landed row is a per-template Tt under two alleles, so
  # a LAMP design reports constructs and no signal at all
  m <- structure(list(scheme = "lamp"), class = "marker_lamp")
  r <- marker_readout(m)
  expect_identical(r$scheme, "lamp")
  expect_false(r$simulated)
  expect_null(r$table)
  expect_true(is.na(r$model))
})

test_that("LAMP discrimination is a window with censoring, never a binary call", {
  t <- ls_read("lamp_badolo_kdr2012_time_to_positive.csv")
  expect_true("censoring" %in% names(t))
  # the non-target result is a sentence, not a token: a censored arm must not
  # carry a point value, and every censored arm must be the non-target allele
  cen <- t[t$censoring != "measured", ]
  expect_gt(nrow(cen), 0L)
  expect_true(all(cen$is_target %in% FALSE))
  expect_true(all(is.na(cen$time_to_positive_min)))
  expect_true(all(cen$observation_window_min == 75))
  expect_match(cen$censoring[1L], "exceeds the window")
  meas <- t[t$censoring == "measured", ]
  expect_true(all(meas$is_target %in% TRUE))
  expect_true(all(is.finite(meas$time_to_positive_min)))
  # the censored arms are above the window, so nothing was measured at >75 min
  expect_true(all(is.na(cen$time_to_positive_min)))
})

test_that("the comp2018 2.7x ratio is an endpoint ratio between oligos, not specificity", {
  s <- ls_read("lamp_comp2018_separation.csv")
  row <- s[!is.na(s$magnitude) & s$magnitude == 2.7, ]
  expect_equal(nrow(row), 1L)
  expect_match(row$signal_type, "endpoint fluorescence")
  expect_match(row$outcome_definition, "NO time-to-positive")
  # decision clause 3: never quoted beside a Tt-based window
  expect_false(any(grepl("specificity", names(s), ignore.case = TRUE)))
  expect_false(any(c("time_to_positive_min", "ct", "tt") %in% tolower(names(s))))
})

test_that("the two Badolo 2012 articles are keyed to their own identifiers", {
  src <- ls_read("scientific_sources.csv")
  mon <- src[src$key == "badolo2012", ]
  lamp <- src[src$key == "badolo_kdr2012", ]
  # the monitoring survey (PMC3489511) and the AS-LAMP paper (PMC3407793) are
  # different articles; the source record used to carry ...-11-231, which
  # resolves to neither
  expect_identical(mon$doi, "10.1186/1475-2875-11-232")
  expect_identical(mon$pmcid, "PMC3489511")
  expect_identical(mon$pmid, 22799568L)
  expect_identical(lamp$doi, "10.1186/1475-2875-11-227")
  expect_identical(lamp$pmcid, "PMC3407793")
  expect_identical(lamp$pmid, 22770418L)

  bl <- ls_read("scientific_backlog.csv")
  b5 <- bl[bl$key == "badolo2012", ]
  expect_identical(b5$doi, "10.1186/1475-2875-11-232")
  expect_match(b5$note, "11-227")
})

test_that("the FIP/F2 3' end is the evidenced 3' mismatch position", {
  p <- ls_read("lamp_mismatch_rtlamp2019_position_rules.csv")
  # two FIP rows: the 5' end of FIP (F1c) and the 3' end of FIP (= F2)
  fip2 <- p[grepl("F2", p$item), ]
  fip5 <- p[grepl("FIP", p$item) & !grepl("F2", p$item), ]
  expect_equal(nrow(fip2), 1L)
  expect_equal(nrow(fip5), 1L)
  # decision clause 4: F2 3' 13.5-17.3 min vs F1c 5' 6.9-13.4 and F3 1.1-3.5
  expect_equal(fip2$conventional_tt_over_wt_min, 13.5)
  expect_equal(fip2$conventional_tt_over_wt_max_min, 17.3)
  expect_gt(fip2$conventional_tt_over_wt_min, fip5$conventional_tt_over_wt_min)
  f3 <- p[grepl("F3", p$item), ]
  expect_lt(f3$conventional_tt_over_wt_max_min, fip2$conventional_tt_over_wt_min)
  # and the artefact says it is a tolerance measure, not discrimination
  expect_true(all(grepl("TOLERANCE|tolerance", p$outcome_definition)))
})

# The backlog is the other half of the registry: what was read, what could not be
# read, and what was never fetched.  A status outside that vocabulary would make
# "how much of the literature is covered" unanswerable.
test_that("scientific_backlog: three states, and every extract has a source", {
  b <- ls_read("scientific_backlog.csv")
  expect_equal(nrow(b), 45L)
  expect_identical(anyDuplicated(b$key), 0L)
  expect_setequal(unique(b$status), c("extracted", "not_extractable", "unfetched"))
  # a family label is what makes a gap visible per assay type
  expect_true(all(nzchar(b$family)))
  expect_gte(length(unique(b$family)), 10L)
  # the anchored claim is why the row exists at all
  expect_true(all(nzchar(b$what_it_would_anchor)))
  s <- ls_read("scientific_sources.csv")
  extracted <- b$key[b$status == "extracted"]
  expect_true(all(extracted %in% s$key))
})
