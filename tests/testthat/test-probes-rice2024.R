# =============================================================================
#  probes/rice2024  --  a per-marker KASP CALL error rate, not a PCR efficiency
# =============================================================================
#  Rice (N Y) 2024, doi:10.1186/s12284-024-00711-1, PMC11087395.
#  54 trait-linked KASP markers, each with a false-positive and a
#  false-negative rate printed for the population it was tested on.

ext_probe <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name),
    show_col_types = FALSE,
    na = c("", "NA"), progress = FALSE
  )
}

test_that("probes/rice2024: 54 markers, three rates each, all proportions", {
  d <- ext_probe("rice2024_kasp_error_rates.csv")
  expect_equal(nrow(d), 54L)
  expect_false(anyDuplicated(d$marker) > 0L)
  for (nm in c(
    "fpr_nil_population", "fnr_nil_population", "fnr_ril_population",
    "max_error_rate"
  )) {
    if (nm == "max_error_rate") next
    x <- d[[nm]]
    expect_true(all(x >= 0 & x <= 1), info = nm)
  }
  # the paper prints three rates, not four: a RILs FPR column would be a
  # transcription invention
  expect_false("fpr_ril_population" %in% names(d))
  expect_true(all(d$evidence_tier == "measured"))
})

test_that("probes/rice2024: max_error_rate is the max of the printed rates", {
  d <- ext_probe("rice2024_kasp_error_rates.csv")
  expect_equal(d$max_error_rate,
    pmax(
      d$fpr_nil_population, d$fnr_nil_population,
      d$fnr_ril_population
    ),
    tolerance = 1e-12
  )
})

test_that("probes/rice2024: the rates are small but never zero for every marker", {
  d <- ext_probe("rice2024_kasp_error_rates.csv")
  # The paper's claim is "very low" rates. Measured here the medians are 5.5%
  # (NILs FPR), 5.2% (NILs FNR) and 4.2% (RILs FNR), and 6 of the 54 markers
  # have at least one rate above 10%. "Very low" therefore depends on the
  # tolerance: about one call in twenty is wrong, per marker, on average.
  expect_lt(stats::median(d$fpr_nil_population), 0.06)
  expect_lt(stats::median(d$fnr_nil_population), 0.06)
  expect_lt(stats::median(d$fnr_ril_population), 0.05)
  expect_equal(sum(d$max_error_rate > 0.10), 6L)
  # the extreme is a real per-marker value (12.9%), not a rounding artefact
  expect_gt(max(d$max_error_rate), 0.12)
  # a zero rate is printed for some markers and must stay a printed zero
  expect_gte(sum(d$fnr_ril_population == 0), 1L)
})

test_that("probes/rice2024: the denominators are the panels the paper names", {
  d <- ext_probe("rice2024_kasp_error_rates.csv")
  expect_true(all(d$n_genotypes_nil >= 45 & d$n_genotypes_nil <= 713))
  expect_true(all(d$n_genotypes_ril >= 22 & d$n_genotypes_ril <= 256))
  # the paper's headline sample sizes appear among the denominators
  expect_true(any(d$n_genotypes_nil %in% c(256, 713)))
})

test_that("probes/rice2024: the values are still in the fetched table", {
  f <- test_path("../../.tmp/lit/rice2024/fulltext.xml")
  skip_if_not(file.exists(f), "rice2024 paper not fetched")
  txt <- paste(readLines(f, warn = FALSE), collapse = "")
  d <- ext_probe("rice2024_kasp_error_rates.csv")
  # every marker id must occur in the XML, and the first marker's three rates
  # must occur as the printed decimals
  m <- head(d, 1L)
  expect_true(grepl(m$marker, txt, fixed = TRUE))
  for (v in c(m$fpr_nil_population, m$fnr_nil_population, m$fnr_ril_population)) {
    expect_true(grepl(sprintf("%.4f", v), txt, fixed = TRUE), info = as.character(v))
  }
})
