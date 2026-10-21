# =============================================================================
#  VCF2CAPS: the "no wet-lab outcomes" premise was false
# =============================================================================
#  Source: PMC8186816 (DOI 10.1371/journal.pcbi.1008980).  10 markers x 30
#  plants, outcomes in prose.  These tests hold the headline counts and the
#  absence of a per-plant table as a stated fact rather than a silence.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name), show_col_types = FALSE,
                  na = c("", "NA"))
}

test_that("VCF2CAPS: the headline is a readable key/value table", {
  h <- ext("vcf2caps2021_headline.csv")
  expect_setequal(names(h), c("quantity", "value"))
  expect_false(anyNA(h$quantity))
  expect_false(anyNA(h$value))
  get <- function(k) h$value[h$quantity == k]
  expect_equal(get("n_markers_verified"), 10)
  expect_equal(get("n_plants_each"), 30)
  expect_equal(get("n_plant_x_marker_observations"), 300)
  # the design was 10 non-cleavable homs + 10 hets + 10 cleavable homs
  expect_equal(get("n_plants_each"), 3 * 10)
  # the per-marker observations reconcile with the counts above
  expect_equal(get("n_plant_x_marker_observations"),
               get("n_markers_verified") * get("n_plants_each"))
})

test_that("VCF2CAPS: the wording groups partition the ten markers", {
  h <- ext("vcf2caps2021_headline.csv")
  groups <- c("exact_match", "weaker_cleavable_allele", "close_single_het",
              "five_noncleavable_reclassified", "cleavable_hom_read_as_het_banding",
              "het_read_as_noncleavable")
  g <- h$value[h$quantity %in% groups]
  expect_length(g, length(groups))
  expect_equal(sum(g), 10)
  # the one class-level correction cannot exceed its own denominator
  expect_lte(h$value[h$quantity == "marker_3_363_reclassified_as_het"], 10)
  expect_gt(h$value[h$quantity == "marker_3_363_reclassified_as_het"], 0)
})

test_that("VCF2CAPS: per-marker rows never claim a count the design forbids", {
  v <- ext("vcf2caps2021_verification.csv")
  expect_equal(nrow(v), 10L)
  expect_true(all(v$n_plants_per_marker == 30L))
  # the design is recorded as the paper states it, per class: "10/10/10"
  expect_true(all(v$n_plants_per_class_design == "10/10/10"))
  per_class <- as.integer(strsplit(v$n_plants_per_class_design[1L], "/", fixed = TRUE)[[1L]])
  expect_equal(sum(per_class), v$n_plants_per_marker[1L])
  expect_true(all(v$n_plant_x_marker_observations == 300L))
  # only ONE marker has a stated reclassification count; the other nine are NA
  # because the prose does not count them.  An NA here is the honest value, so
  # the test asserts the shape of the silence rather than filling it in.
  stated <- !is.na(v$observed_reclassified_as_het)
  expect_equal(sum(stated), 1L)
  expect_equal(v$marker_abbrev[stated], "3_363")
  # where a count IS stated it must fit inside the class it was predicted into
  expect_true(all(v$observed_reclassified_as_het[stated] <= max(per_class)))
  expect_true(all(v$observed_supports_prediction[stated] >= 0L))
  # and the unstated rows must not carry a number in the support column
  expect_true(all(is.na(v$observed_supports_prediction[!stated])))
  # every marker keeps the paper's own wording, so the row can be audited
  expect_true(all(nzchar(v$reported_verdict_verbatim)))
  expect_true(all(grepl("prose", v$evidence_basis, ignore.case = TRUE)))
  # the absence of a call table is stated, not left implicit
  expect_true(all(grepl("no per-plant call table", v$evidence_basis, fixed = TRUE)))
})

test_that("VCF2CAPS: the wrong dispatch DOI is not the one recorded", {
  v <- ext("vcf2caps2021_verification.csv")
  expect_true(all(v$doi == "10.1371/journal.pcbi.1008980"))
  # 10.1371/journal.pcbi.1009001 resolves to an unrelated paper
  expect_false(any(v$doi == "10.1371/journal.pcbi.1009001"))
  expect_true(all(v$pmcid == "PMC8186816"))
})
