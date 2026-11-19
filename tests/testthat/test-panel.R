# GWAS-to-panel end to end (roadmap ZY1): one row per target, a design or its
# cause chain, the one-pot estimate, and an abort when nothing designs.
#
# The fixture is committed and hand-written: a 600 bp Chr1 (no homopolymer run
# longer than 4 bp, with the target at base 300) and a 300 bp Chr2 homopolymer,
# plus three phased-haplotype rows (H1 carries the Chr1:300 ALT at 40 %, both
# haplotypes carry Chr1:100 C>G).  Nothing here is random, so the output is
# reproducible.

fixture_genome <- function() {
  Biostrings::readDNAStringSet(test_path("fixtures", "panel", "genome.fa"))
}

fixture_haplotypes <- function() {
  readr::read_csv(test_path("fixtures", "panel", "haplotypes.csv"), show_col_types = FALSE)
}

fixture_targets <- function() {
  tibble::tibble(
    chrom = c("Chr1", "Chr2"),
    pos = c(300L, 150L),
    ref = c("A", "A"),
    alt = c("C", "C")
  )
}

design_panel <- function(dir, targets = fixture_targets(), chemistries = c("snp", "hrm"),
                         variants = NULL, haplotypes = fixture_haplotypes()) {
  marker_design_panel(
    fixture_genome(),
    variants = variants,
    haplotypes = haplotypes,
    chemistries = chemistries,
    targets = targets,
    dir = dir
  )
}

test_that("every target yields a design or its full cause chain, never a drop", {
  skip_if_no_oligoarrayaux()
  dir <- withr::local_tempdir()
  panel <- design_panel(dir, chemistries = c("snp", "hrm", "caps"))

  expect_s3_class(panel, "marker_panel")
  expect_true(panel$experimental)

  t <- panel$targets
  expect_equal(nrow(t), 2L)
  expect_setequal(t$target, c("Chr1:300 A>C", "Chr2:150 A>C"))
  # each row is exactly one of: designed, or failed with a cause chain
  expect_true(all(ifelse(t$ok, lengths(t$cause) == 0L, lengths(t$cause) > 0L)))

  good <- t[t$target == "Chr1:300 A>C", , drop = FALSE]
  expect_true(good$ok)
  expect_identical(good$scheme, "snp")
  expect_s3_class(good$assay[[1L]], "marker_assay")
  expect_true(is.finite(good$footprint_dropout))
  expect_gte(good$footprint_dropout, 0)
  expect_lte(good$footprint_dropout, KASP_LIMITS$footprint_dropout)

  bad <- t[t$target == "Chr2:150 A>C", , drop = FALSE]
  expect_false(bad$ok)
  expect_true(is.na(bad$scheme))
  expect_null(bad$assay[[1L]])
  expect_true(is.na(bad$footprint_dropout))
  causes <- bad$cause[[1L]]
  expect_gte(length(causes), 1L)
  expect_true(all(grepl("^[a-z]+: .+", causes)))
  expect_setequal(sub(":.*$", "", causes), c("snp", "hrm", "caps"))

  # the one-pot estimate holds every designed assay, placed or explained
  pot <- panel$one_pot
  expect_true(pot$experimental)
  expect_equal(nrow(pot$groups) + nrow(pot$unassigned), sum(t$ok))

  # the method itself (the controller's roxygenise adds the S3 registration)
  out <- capture.output(tidymarkers:::print.marker_panel(panel))
  expect_true(any(grepl("<marker_panel> 2 target(s), 1 designed, 1 failed", out, fixed = TRUE)))
  expect_true(any(grepl("footprint_dropout", out, fixed = TRUE)))
  expect_true(any(grepl("Chr2:150 A>C: snp: ", out, fixed = TRUE)))
  expect_true(any(grepl("experimental", out, fixed = TRUE)))
})

test_that("no chemistry working aborts naming every target and every cause", {
  skip_if_no_oligoarrayaux()
  bad <- tibble::tibble(chrom = "Chr2", pos = 150L, ref = "A", alt = "C")
  err <- tryCatch(
    design_panel(withr::local_tempdir(), targets = bad, chemistries = c("snp", "hrm")),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_match(msg, "Chr2:150 A>C", fixed = TRUE)
  expect_match(msg, "No AS primer passed the GC/Tm filter", fixed = TRUE)
  expect_match(msg, "returned no marker object", fixed = TRUE)
})

test_that("variant context and a haplotype path feed the plan, not the panel", {
  skip_if_no_oligoarrayaux()
  # Chr1:100 C>G sits in the haplotype fixture and away from the target; the
  # design of Chr1:300 A>C stays the same, so the context path changes nothing
  # the panel claims beyond the population check it feeds.  The haplotypes come
  # from the CSV path here, which is the other accepted input shape.
  variants <- tibble::tibble(chrom = "Chr1", pos = 100L, ref = "C", alt = "G", pop_freq = 0.6)
  panel <- design_panel(
    withr::local_tempdir(),
    variants = variants,
    haplotypes = test_path("fixtures", "panel", "haplotypes.csv")
  )

  good <- panel$targets[panel$targets$target == "Chr1:300 A>C", , drop = FALSE]
  expect_true(good$ok)
  expect_identical(good$scheme, "snp")
  expect_identical(good$footprint_dropout, 0)
})

test_that("the same call twice gives the same targets and one-pot estimate", {
  skip_if_no_oligoarrayaux()
  p1 <- design_panel(withr::local_tempdir())
  p2 <- design_panel(withr::local_tempdir())

  # marker_assay$decision carries an igraph pointer and a random id, so the
  # panel's determinism is the design, the tolerance and the pot grouping
  expect_identical(p1$targets$target, p2$targets$target)
  expect_identical(p1$targets$scheme, p2$targets$scheme)
  expect_identical(p1$targets$footprint_dropout, p2$targets$footprint_dropout)
  expect_identical(p1$targets$cause, p2$targets$cause)
  expect_identical(p1$targets$ok, p2$targets$ok)
  expect_identical(p1$targets$assay[[1L]]$assay, p2$targets$assay[[1L]]$assay)
  expect_identical(p1$one_pot, p2$one_pot)
})
