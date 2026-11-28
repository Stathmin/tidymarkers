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
                         variants = NULL, haplotypes = fixture_haplotypes(),
                         ld = NULL, costs = NULL, created = NULL) {
  marker_design_panel(
    fixture_genome(),
    variants = variants,
    haplotypes = haplotypes,
    chemistries = chemistries,
    targets = targets,
    ld = ld,
    costs = costs,
    created = created,
    dir = dir
  )
}

test_that("every target yields a design or its full cause chain, never a drop", {
  skip_if_no_oligoarrayaux()
  dir <- withr::local_tempdir()
  panel <- design_panel(dir, chemistries = c("snp", "hrm", "caps"))

  expect_s3_class(panel, "marker_panel")
  expect_false(panel$experimental)

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

# -- the seams (T4/T5) -------------------------------------------------------

fixture_panel_file <- function(name) {
  readr::read_csv(test_path("fixtures", "panel", name), show_col_types = FALSE)
}

ld_targets <- function() {
  tidymarkers:::normalise_targets(tibble::tibble(
    chrom = c("Chr1", "Chr1"), pos = c(300L, 500L),
    ref = c("A", "T"), alt = c("C", "G")
  ))
}

test_that("the haplotype contract rejects unphased input and a frequency sum above one", {
  targets <- tidymarkers:::normalise_targets(
    tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "C")
  )
  unphased <- fixture_panel_file("haplotypes_unphased.csv")
  err <- tryCatch(
    tidymarkers:::panel_haplotype_contract(unphased, targets),
    error = function(e) e
  )
  expect_s3_class(err, "tidymarkers_haplotype_error")
  expect_match(conditionMessage(err), "H1", fixed = TRUE)
  expect_match(conditionMessage(err), "Chr1:300", fixed = TRUE)

  over <- tibble::tibble(
    haplotype = c("H1", "H2"), freq = c(0.7, 0.6), chrom = "Chr1",
    pos = c(300L, 100L), ref = c("A", "C"), alt = c("C", "G")
  )
  err <- tryCatch(
    tidymarkers:::panel_haplotype_contract(over, targets),
    error = function(e) e
  )
  expect_s3_class(err, "tidymarkers_haplotype_error")
  expect_match(conditionMessage(err), "Chr1", fixed = TRUE)

  # absence at a target is reference and is counted, not silently invented
  h <- fixture_panel_file("haplotypes.csv")
  ct <- tidymarkers:::panel_haplotype_contract(h, targets)
  expect_identical(attr(ct, "reference_calls"), 1L)

  # a non-target row must exist in the supplied variant table
  variants <- tibble::tibble(chrom = "Chr1", pos = 100L, ref = "C", alt = "A")
  err <- tryCatch(
    tidymarkers:::panel_haplotype_contract(h, targets, variants = variants),
    error = function(e) e
  )
  expect_s3_class(err, "tidymarkers_haplotype_error")
  expect_match(conditionMessage(err), "Chr1:100", fixed = TRUE)
})

test_that("a duplicate target aborts before any chemistry runs", {
  duplicate <- rbind(fixture_targets()[1L, ], fixture_targets()[1L, ])
  err <- tryCatch(
    marker_design_panel(
      fixture_genome(),
      chemistries = "snp",
      targets = duplicate, haplotypes = fixture_haplotypes()
    ),
    error = function(e) e
  )
  expect_s3_class(err, "tidymarkers_duplicate_target")
  expect_match(conditionMessage(err), "Duplicate target", fixed = TRUE)
})

test_that("the LD seam collapses a perfect pair only when a threshold is given", {
  h <- fixture_panel_file("haplotypes_ld.csv")
  targets <- ld_targets()

  not_applied <- tidymarkers:::panel_ld(targets, h, NULL)
  expect_equal(nrow(not_applied$targets), 2L)
  expect_identical(not_applied$record$method, "not_applied")
  expect_null(not_applied$record$dropped)

  pruned <- tidymarkers:::panel_ld(targets, h, 0.8)
  expect_equal(nrow(pruned$targets), 1L)
  expect_equal(pruned$targets$pos, 300L)
  expect_identical(pruned$record$method, "threshold")
  expect_equal(pruned$record$threshold, 0.8)
  expect_identical(pruned$record$dropped, "Chr1:500 T>G")
  expect_identical(pruned$record$reason, "ld")

  # an explicit linkage table needs no threshold
  table_ld <- tibble::tibble(target1 = "Chr1:300 A>C", target2 = "Chr1:500 T>G")
  expect_equal(nrow(tidymarkers:::panel_ld(targets, h, table_ld)$targets), 1L)
})

test_that("the panel applies LD and writes the record into the artefact", {
  skip_if_no_oligoarrayaux()
  h <- fixture_panel_file("haplotypes_ld.csv")
  targets <- ld_targets()

  both <- design_panel(withr::local_tempdir(), targets = targets, haplotypes = h)
  expect_equal(nrow(both$targets), 2L)
  expect_identical(both$ld$method, "not_applied")
  expect_identical(both$artefact$manifest$plan$ld$method, "not_applied")

  one <- design_panel(withr::local_tempdir(), targets = targets, haplotypes = h, ld = 0.8)
  expect_equal(nrow(one$targets), 1L)
  expect_identical(one$targets$target, "Chr1:300 A>C")
  expect_identical(one$artefact$manifest$plan$ld$method, "threshold")
  expect_equal(one$artefact$manifest$plan$ld$threshold, 0.8)
  expect_identical(one$artefact$manifest$plan$ld$dropped, "Chr1:500 T>G")
  expect_false(identical(both$artefact$id, one$artefact$id))

  pip <- one$artefact$manifest$plan$panel
  expect_identical(
    pip$order,
    c("gwas", "ld", "haplotypes", "chemistry", "artefact")
  )
  expect_identical(pip$gwas$n_signals, 2L)
  expect_identical(pip$ld$method, "threshold")
  expect_identical(pip$haplotypes$reference_calls, 1L)
  expect_identical(pip$chemistry$chosen[["Chr1:300 A>C"]], "snp")
  expect_identical(pip$one_pot$groups, 1L)
})

test_that("LD runs before the haplotype contract", {
  skip_if_no_oligoarrayaux()
  targets <- ld_targets()
  # H1 is unphased at Chr1:500, the target ld = 0.8 drops; the same row on a
  # surviving target must still abort
  unphased <- tibble::tibble(
    haplotype = c("H1", "H1", "H1", "H2"),
    freq = c(0.4, 0.4, 0.4, 0.6),
    chrom = "Chr1", pos = c(300L, 500L, 500L, 100L),
    ref = c("A", "T", "T", "C"), alt = c("C", "G", "A", "G")
  )
  pruned <- design_panel(
    withr::local_tempdir(),
    targets = targets, haplotypes = unphased, ld = 0.8
  )
  expect_identical(pruned$targets$target, "Chr1:300 A>C")

  err <- tryCatch(
    design_panel(
      withr::local_tempdir(),
      targets = targets, haplotypes = unphased
    ),
    error = function(e) e
  )
  expect_s3_class(err, "tidymarkers_haplotype_error")
  expect_match(conditionMessage(err), "Chr1:500", fixed = TRUE)
})

test_that("reversing the chemistry costs changes the chosen scheme", {
  skip_if_no_oligoarrayaux()
  target <- fixture_targets()[1L, ]
  cheap_snp <- design_panel(withr::local_tempdir(), targets = target)
  expect_identical(cheap_snp$targets$scheme, "snp")
  cheap_hrm <- design_panel(
    withr::local_tempdir(),
    targets = target,
    costs = marker_costs(snp = 5, hrm = 1)
  )
  expect_identical(cheap_hrm$targets$scheme, "hrm")
})

test_that("a fixed created time makes the text artefact byte-identical", {
  expect_identical(
    tidymarkers:::artefact_created("2026-01-01T00:00:00Z"),
    "2026-01-01T00:00:00Z"
  )
  expect_identical(tidymarkers:::artefact_created(0), "1970-01-01T00:00:00Z")
  withr::with_envvar(c(SOURCE_DATE_EPOCH = "86400"), {
    expect_identical(tidymarkers:::artefact_created(), "1970-01-02T00:00:00Z")
  })

  skip_if_no_oligoarrayaux()
  target <- fixture_targets()[1L, ]
  created <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  first <- design_panel(withr::local_tempdir(), targets = target, created = created)
  second <- design_panel(withr::local_tempdir(), targets = target, created = created)

  files <- c("artefact.json", "decisions.csv", "checks.csv", "why.csv", "targets.csv")
  for (file in files) {
    expect_identical(
      readBin(file.path(first$artefact$dir, file), "raw", file.size(file.path(first$artefact$dir, file))),
      readBin(file.path(second$artefact$dir, file), "raw", file.size(file.path(second$artefact$dir, file))),
      info = file
    )
  }
  # excluded by name: solved.rds serialises an igraph pointer and a random id,
  # so a byte comparison would test the serializer, not the design
  expect_true(file.exists(file.path(first$artefact$dir, "solved.rds")))
})

test_that("the panel is a typed input to pareto, artefact and export", {
  skip_if_no_oligoarrayaux()
  dir <- withr::local_tempdir()
  panel <- design_panel(dir, targets = fixture_targets()[1L, ], created = "2026-01-01T00:00:00Z")

  solved <- as_marker_solved(panel)
  expect_s3_class(solved, "marker_solved")
  expect_identical(solved$report, panel$solved$report)
  expect_identical(attr(solved, "panel"), panel$artefact$manifest$plan$panel)

  front <- marker_pareto(panel)
  expect_s3_class(front, "tbl_df")
  expect_true(all(c("label", "cost", "on_front") %in% names(front)))

  export_dir <- withr::local_tempdir()
  expect_identical(marker_export(panel, dir = export_dir, artefact = FALSE), export_dir)
  expect_true(all(file.exists(file.path(export_dir, c("markers.csv", "attempts.csv", "primers.csv")))))

  art <- marker_artefact(panel, dir = dir, created = "2026-01-01T00:00:00Z")
  expect_s3_class(art, "marker_artefact")
  expect_true(file.exists(file.path(art$dir, "artefact.json")))
})

test_that("a panel variant table position is whole in its spelling (ADR 0076)", {
  bad <- tibble::tibble(chrom = "Chr1", pos = 100.9, ref = "C", alt = "A")
  expect_error(tidymarkers:::panel_variants(bad), "positive whole `pos`")
  bad$pos <- "100.9"
  expect_error(tidymarkers:::panel_variants(bad), "positive whole `pos`")
  ok <- tibble::tibble(chrom = "Chr1", pos = 100L, ref = "C", alt = "A")
  expect_identical(tidymarkers:::panel_variants(ok), ok)
})
