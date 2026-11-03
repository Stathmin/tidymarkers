# Whole-amplicon/target folding (Z10.1): the raw engine answer plus the limits a
# LAMP/RPA/Cas13-RNA or GC-rich HRM design has to live with, or "unknown".

test_that("fold_amplicon() folds a GC-rich target and returns dG and structure", {
  skip_if_no_oligoarrayaux()
  skip_if_not(
    fs::file_exists(fs::path(tidymarkers:::oligoarrayaux_paths()$bin, "hybrid-ss-min")),
    "hybrid-ss-min not installed"
  )
  # a two-helix GC-rich target, as HRM/LAMP amplicons are
  res <- fold_amplicon(paste0("GCGCGCGCGCGGGGCCCCCAATT", "GCGCGCGCGCGGGGCCCCCAATT"))
  expect_identical(res$status, "ok")
  expect_true(is.finite(res$mfe))
  expect_lt(res$mfe, 0) # a stem that long is stabilising
  expect_identical(nchar(res$structure), res$length)
  expect_false(grepl("[^().]", res$structure))
  expect_true(any(strsplit(res$structure, "")[[1L]] != "."))
  # the whole sequence was folded, not a primer prefix
  expect_identical(res$length, 46L)
})

test_that("fold_amplicon() answers unknown, never a number, without the engine", {
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  expect_identical(capabilities(marker_engines()$thermo, "fold"), FALSE)
  res <- fold_amplicon("GCGCGCGCGCGGGGCCCCCAATTGCGCGCGCGCGGGGCCCCCAATT")
  expect_identical(res$status, "unknown")
  expect_true(is.na(res$mfe))
  expect_true(is.na(res$structure))
  expect_match(res$note, "unknown")
})

test_that("an RNA target is labelled RNA, never silently folded as DNA", {
  skip_if_no_oligoarrayaux()
  res <- fold_amplicon("ACGUACGUACGUGGCAUGCUAGC", kind = "RNA")
  # this OligoArrayAux build has no RNA energy tables: the honest answer is
  # unknown, and it is reported as an RNA result rather than a DNA one
  expect_identical(res$kind, "RNA")
  expect_false(isTRUE(res$salt_modelled))
  if (identical(res$status, "ok")) {
    expect_true(is.finite(res$mfe))
    # folded as RNA: the U-alphabet target is accepted (a DNA-only engine rejects it)
    expect_identical(res$length, 23L)
  } else {
    expect_identical(res$status, "unknown")
    expect_true(is.na(res$mfe))
    expect_match(res$note, "not folded as DNA", fixed = TRUE)
  }
  # even when unknown, the label says which alphabet was asked for
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  expect_identical(fold_amplicon("ACGUACGU", kind = "RNA")$kind, "RNA")
})

test_that("fold_amplicon() rejects a bad alphabet and an empty sequence", {
  expect_error(fold_amplicon("ACGT", kind = "protein"), "arg")
  expect_error(fold_amplicon("  "), "sequence")
})

# Z10.1: the amplicon fold check.  The whole amplicon is folded, not a primer.
AMP <- paste0(
  "ACGTTGCAAGGCTTACGATCGATCGGATCCGTTACGGATCCTAGCATCGATGGCATCG",
  "GATCGTAGCTAGGCATCGGATCCAAGCTTAGCGATCGATCGGATCCGATGCTAGCAT"
)

test_that("amplicon_fold_check() is unknown, never a number, without the engine", {
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  expect_identical(capabilities(marker_engines()$thermo, "fold"), FALSE)
  res <- amplicon_fold_check(AMP)
  expect_identical(nrow(res), 1L)
  expect_identical(res$status, "unknown")
  expect_true(is.na(res$mfe))
  expect_true(is.na(res$single_stranded))
  expect_match(res$note, "unknown")
})

test_that("amplicon_fold_check() folds the whole amplicon, not the primer", {
  skip_if_no_oligoarrayaux()
  skip_if_not(isTRUE(capabilities(marker_engines()$thermo, "fold")), "fold engine absent")
  res <- amplicon_fold_check(AMP)
  expect_identical(res$status %in% c("pass", "warn"), TRUE)
  expect_identical(res$length, nchar(AMP))
  # the oligo engine would fold the first 20 nt only: a different question
  expect_gt(res$length, 20L)
  expect_true(is.finite(res$mfe))
  expect_true(is.finite(res$single_stranded))
  expect_true(res$single_stranded >= 0 && res$single_stranded <= 1)
})

test_that("an RNA target is unknown, not a DNA folding number", {
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  res <- amplicon_fold_check("ACGUACGUACGUGGCAUGCUAGC", kind = "RNA")
  expect_identical(nrow(res), 1L)
  expect_identical(res$kind, "RNA")
  expect_identical(res$status, "unknown")
  expect_true(is.na(res$mfe))
})

test_that("amplicon_fold_check() is unknown past max_length and empty input errors", {
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  res <- amplicon_fold_check(c(short = AMP, long = paste(rep("ACGT", 100L), collapse = "")),
    max_length = 200L
  )
  expect_identical(res$length, c(nchar(AMP), 400L))
  expect_identical(res$status, c("unknown", "unknown"))
  expect_error(amplicon_fold_check(c("ACGT", "")), "sequence")
  expect_error(amplicon_fold_check(character(0)), "sequence")
})
