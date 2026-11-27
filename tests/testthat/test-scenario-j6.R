# J6: a released artefact carries its contract, exports, and replays.
#
# The scenario starts at `marker_artefact()` (ADR 0011): the record is the input
# to every release export, and the replay reads it rather than a private copy.

test_that("J6: the artefact records its contract and replays to the same verdicts", {
  skip_if_no_oligoarrayaux()
  genome <- read_genome_fasta(scenario_file("j1_genome.fa"))
  targets <- utils::read.csv(scenario_file("j1_target.csv"))
  plan <- marker_plan(genome) +
    marker_targets(targets) +
    marker_tech(schemes = "snp")
  solved <- marker_solve(plan, progress = FALSE)

  dir <- withr::local_tempdir("j6-artefact-")
  art <- marker_artefact(
    solved, dir,
    created = "2026-01-01T00:00:00Z", scenario = "J6"
  )
  m <- art$manifest
  expect_identical(m$schema, "tidymarkers/artefact@2")
  expect_identical(m$contract$version, "tidymarkers/contract@1")
  expect_identical(m$contract$scenario, "J6")
  expect_true(any(vapply(m$inputs, function(i) identical(i$role, "targets"), logical(1))))
  expect_identical(m$files$solved$name, "solved.rds")
  expect_match(m$files$solved$excluded, "igraph")
  expect_true(file.exists(fs::path(dir, "solved.rds")))

  # the same plan under another scenario is a different record (the contract is
  # part of the id), while the timestamp alone is not
  other <- marker_artefact(
    solved, withr::local_tempdir("j6-other-"),
    created = "2026-01-01T00:00:00Z", scenario = "J5"
  )
  expect_false(identical(other$manifest$id, m$id))
  expect_identical(
    marker_artefact(solved, withr::local_tempdir("j6-again-"), created = m$created, scenario = "J6")$manifest$id,
    m$id
  )

  back <- marker_artefact_read(dir)
  expect_identical(back$manifest$contract$scenario, "J6")

  crate <- marker_export_rocrate(art, withr::local_tempdir("j6-crate-"))
  meta <- jsonlite::fromJSON(
    fs::path(crate, "ro-crate-metadata.json"),
    simplifyVector = FALSE
  )
  root <- Filter(function(e) identical(e[["@id"]], "./"), meta[["@graph"]])[[1L]]
  expect_identical(root$version, "tidymarkers/contract@1")
  expect_true("scenario:J6" %in% unlist(root$keywords))

  st26 <- fs::path(withr::local_tempdir("j6-st26-"), "sequence_listing.xml")
  suppressWarnings(marker_patent_export(solved, "st26", st26))
  expect_true(file.exists(st26))
  # the XSD is operator-supplied; without it the validator must name the gap
  # instead of claiming a pass (declared degradation, ADR 0011)
  xsd <- tryCatch(marker_patent_validate(st26), error = identity)
  expect_true(inherits(xsd, "error") || isTRUE(xsd))

  replay <- marker_replay(art)
  expect_s3_class(replay, "marker_replay")
  expect_true(replay$match)
  expect_true(all(replay$table$match))
})

test_that("J6: the text artefact is byte-identical and excludes the RDS", {
  skip_if_no_oligoarrayaux()
  genome <- read_genome_fasta(scenario_file("j1_genome.fa"))
  targets <- utils::read.csv(scenario_file("j1_target.csv"))
  solved <- marker_solve(
    marker_plan(genome) + marker_targets(targets) + marker_tech(schemes = "snp"),
    progress = FALSE
  )
  created <- "2026-01-01T00:00:00Z"
  a <- marker_artefact(solved, withr::local_tempdir("j6-a-"), created = created)
  b <- marker_artefact(solved, withr::local_tempdir("j6-b-"), created = created)

  text <- c("artefact.json", "decisions.csv", "checks.csv", "why.csv", "targets.csv")
  for (file in text) {
    expect_identical(
      readBin(fs::path(a$dir, file), "raw", file.info(fs::path(a$dir, file))$size),
      readBin(fs::path(b$dir, file), "raw", file.info(fs::path(b$dir, file))$size),
      info = file
    )
  }
  expect_identical(a$manifest$files$solved$name, "solved.rds")
  expect_identical(a$manifest$id, b$manifest$id)
})
