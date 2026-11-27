# J3: GWAS targets -> panel -> pareto, artefact and export.  The committed
# genome and haplotype fixtures are the same ones the panel tests use, so the
# scenario proves the typed chain, not a new design.

test_that("J3: a panel is a typed step to pareto, artefact and export", {
  skip_if_no_oligoarrayaux()
  scenario_timer("J3")
  genome_path <- scenario_file("j1_genome.fa")
  targets_path <- scenario_file("j1_target.csv")
  haplotypes_path <- scenario_file("j3_haplotypes.csv")
  created <- "2026-01-01T00:00:00Z"
  dir <- withr::local_tempdir()

  genome <- scenario_step("j3", "ingest", "genome", read_genome_fasta(genome_path))
  expect_scenario_class(genome, "genome")
  targets <- scenario_step("j3", "ingest", "targets", utils::read.csv(targets_path))
  haplotypes <- scenario_step("j3", "ingest", "haplotypes", utils::read.csv(haplotypes_path))

  panel <- scenario_step(
    "j3", "design", "panel",
    marker_design_panel(
      genome,
      haplotypes = haplotypes,
      chemistries = c("snp", "hrm"),
      targets = targets,
      created = created,
      dir = dir
    )
  )
  expect_scenario_class(panel, "marker_panel")
  expect_identical(panel$targets$target, "Chr1:300 A>G")
  expect_true(panel$targets$ok[[1L]])
  expect_identical(panel$ld$method, "not_applied")
  expect_identical(
    panel$artefact$manifest$plan$panel$order,
    c("gwas", "ld", "haplotypes", "chemistry", "artefact")
  )

  solved <- scenario_step("j3", "design", "solved", as_marker_solved(panel))
  expect_scenario_class(solved, "marker_solved")
  expect_identical(solved$report$label, "Chr1:300 A>G")

  pareto <- scenario_step("j3", "design", "pareto", marker_pareto(solved))
  expect_scenario_class(pareto, "tbl_df")
  expect_true(all(c("label", "cost", "on_front") %in% names(pareto)))

  artefact <- scenario_step(
    "j3", "artefact", "artefact",
    marker_artefact(solved, dir = dir, created = created, scenario = "J3")
  )
  expect_scenario_class(artefact, "marker_artefact")
  expect_true(file.exists(file.path(dir, "artefact.json")))

  scenario_step(
    "j3", "artefact", "export",
    marker_export(solved, dir = dir, artefact = FALSE)
  )
  expect_true(all(file.exists(file.path(
    dir, c("markers.csv", "attempts.csv", "primers.csv")
  ))))

  # the same inputs and timestamp must reproduce the manifest byte for byte
  dir2 <- withr::local_tempdir()
  panel2 <- marker_design_panel(
    genome,
    haplotypes = haplotypes,
    chemistries = c("snp", "hrm"),
    targets = targets,
    created = created,
    dir = dir2
  )
  expect_identical(
    readLines(file.path(dir, "artefact.json")),
    readLines(file.path(dir2, "artefact.json"))
  )
  expect_identical(panel2$artefact$id, artefact$id)
})
