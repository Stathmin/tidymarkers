# J2: a VCF batch report continues as a typed object, then checks and CSV.
# ADR 0008: `as_marker_solved()` is the only bridge from the report.

test_that("J2: VCF to a solved batch with checks and exported tables", {
  skip_if_no_oligoarrayaux()
  genome <- read_genome_fasta(scenario_file("j1_genome.fa"))
  variants <- scenario_step(
    "J2", "ingest", "variants",
    read_variants_vcf(scenario_file("j2_loci.vcf"))
  )
  expect_scenario_class(variants, "variants")
  expect_equal(nrow(variants), 2L)

  report <- scenario_step(
    "J2", "design", "report",
    suppressWarnings(design_kasp_markers(genome, variants, flank = 150L))
  )
  expect_scenario_class(report, "marker_batch")
  expect_identical(
    names(report),
    c(
      "chrom", "pos", "ref", "alt", "type", "strand", "amplicon",
      "discrimination", "n_fail", "verdict", "note"
    )
  )
  expect_true(all(report$verdict %in% c("pass", "warn")))

  solved <- scenario_step("J2", "design", "solved", as_marker_solved(report))
  expect_scenario_class(solved, "marker_solved")
  expect_equal(nrow(solved$report), 2L)
  expect_true(all(solved$report$ok))
  expect_equal(length(solved$report$marker), 2L)
  expect_null(solved$plan)

  checks <- scenario_step("J2", "checks", "solved", marker_design_checks(solved))
  expect_equal(sort(unique(checks$label)), sort(solved$report$label))
  expect_true(all(c("label", "criterion", "pass") %in% names(checks)))
  expect_gt(nrow(checks), 0L)

  dir <- withr::local_tempdir(pattern = "scenario-j2")
  scenario_step("J2", "artefact", "tables", marker_export(solved, dir, artefact = FALSE))
  expect_true(all(c("markers.csv", "attempts.csv", "primers.csv", "solved.rds") %in% list.files(dir)))
})
