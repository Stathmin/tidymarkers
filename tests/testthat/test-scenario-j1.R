# J1: a FASTA plus one target goes to checks, a readable reason and CSV files.
# ADR 0008: between typed steps nothing is pulled out of the object by hand.

test_that("J1: genome to checked, explained and exported marker", {
  skip_if_no_oligoarrayaux()
  scenario_timer("J1")
  genome <- scenario_step("J1", "ingest", "genome", read_genome_fasta(scenario_file("j1_genome.fa")))
  expect_scenario_class(genome, "genome")

  target <- utils::read.csv(scenario_file("j1_target.csv"), colClasses = "character")
  target$pos <- as.integer(target$pos)
  assay <- scenario_step(
    "J1", "design", "assay",
    assemble_kasp_assay(genome, target, variants = target, flank = 150L, cores = 1L)
  )
  expect_scenario_class(assay, "marker_assay")
  expect_true(is_target_designable(assay))

  checks <- scenario_step("J1", "checks", "assay", marker_design_checks(assay))
  expect_s3_class(checks, "data.frame")
  expect_true(all(checks$pass))

  explain <- scenario_step("J1", "checks", "explain", marker_explain(checks, scheme = "snp"))
  expect_scenario_class(explain, "marker_explain")

  dir <- withr::local_tempdir(pattern = "scenario-j1")
  files <- scenario_step(
    "J1", "artefact", "results",
    write_marker_results(assay, list(symbol = "FIXTURE1"), dir = dir)
  )
  expect_true(all(file.exists(files)))
  expect_true(all(c("marker_design.csv", "marker_primer_qc.csv") %in% basename(files)))

  primers <- scenario_step(
    "J1", "artefact", "primers",
    export_primers(assay, list(symbol = "FIXTURE1"), dir = dir)
  )
  expect_s3_class(primers, "data.frame")
  expect_gt(nrow(primers), 0L)
  expect_true(any(grepl("_idt[.]csv$", list.files(dir))))
})
