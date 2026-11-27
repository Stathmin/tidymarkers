# =============================================================================
#  ZY23: the artefact as an RO-Crate
# =============================================================================
#  A crate is the artefact plus a metadata descriptor a WorkflowHub or Galaxy
#  reader can follow.  These tests build one real artefact, export it, and
#  check three things: the crate carries exactly the files the manifest names
#  (and nothing else from the artefact directory), the JSON-LD graph says what
#  the crate claims, and an import re-checks the hashes before returning the
#  artefact.
#
#  There is no official RO-Crate validator in this tree - they are Python and
#  Node tools - so these are structural checks, not a validation.  The graph is
#  read with `simplifyVector = FALSE` and asserted field by field; a crate that
#  a real validator would reject on vocabulary grounds is still possible, but a
#  crate that disagrees with its own files is not.

rocrate_fixture <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      skip_if_no_oligoarrayaux()
      genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
      loci <- utils::read.csv(
        system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
        colClasses = "character"
      )[, c("chrom", "pos", "ref", "alt")]
      plan <- marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm")
      solved <- suppressWarnings(marker_solve(plan, progress = FALSE))
      artefact <- marker_artefact(solved, dir = tempfile("rocrate-ref"))
      value <<- list(artefact = artefact, solved = solved)
    }
    value
  }
})

rocrate_manifest_files <- function(artefact) {
  # an excluded file is a record (name + reason), not a sidecar the crate copies
  declared <- artefact$manifest$files
  declared <- declared[vapply(declared, is.character, logical(1))]
  files <- unlist(declared, use.names = FALSE)
  files <- as.character(files)
  c("artefact.json", files[!is.na(files) & nzchar(files)])
}

test_that("export copies exactly the files the manifest names", {
  art <- rocrate_fixture()$artefact
  # A stray file in the artefact directory is not part of the record.
  stray <- fs::path(art$dir, "stray.txt")
  writeLines("not part of the artefact", stray)

  crate <- marker_export_rocrate(art, dir = tempfile("crate"))
  expect_type(crate, "character")
  expect_length(crate, 1L)

  expected <- sort(unique(c(rocrate_manifest_files(art), "ro-crate-metadata.json")))
  actual <- sort(basename(fs::dir_ls(crate, type = "file")))
  expect_setequal(actual, expected)
  expect_false("stray.txt" %in% actual)
})

test_that("a manifest path cannot escape the artefact directory", {
  art <- rocrate_fixture()$artefact

  # the manifest is input, not trusted state: a path that leaves the directory
  # is refused before any file is copied
  tampered <- art
  tampered$manifest$files$decisions <- "../escape.csv"
  expect_error(
    marker_export_rocrate(tampered, dir = tempfile("crate")),
    "unsafe file path"
  )

  # the reader refuses the same name in a persisted artefact
  dst <- withr::local_tempdir()
  fs::file_copy(fs::dir_ls(art$dir, type = "file"), dst)
  manifest <- jsonlite::fromJSON(
    fs::path(art$dir, "artefact.json"),
    simplifyVector = FALSE
  )
  manifest$files$decisions <- "../escape.csv"
  writeLines(
    jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
    fs::path(dst, "artefact.json")
  )
  expect_error(marker_artefact_read(dst), "unsafe file path")
})

test_that("the artefact writer refuses to write through a link", {
  fixture <- rocrate_fixture()
  dir <- withr::local_tempdir()
  outside <- fs::path(withr::local_tempdir(), "decisions.csv")
  writeLines("keep", outside)
  skip_if_not(file.symlink(outside, fs::path(dir, "decisions.csv")), "this filesystem cannot make a symlink")

  expect_error(marker_artefact(fixture$solved, dir = dir), "which is a link", class = "rlang_error")
  expect_identical(readLines(outside), "keep")
})

test_that("an artefact sidecar cannot link out of the record", {
  dir <- withr::local_tempdir()
  outside <- fs::path(withr::local_tempdir(), "outside.csv")
  writeLines("label", outside)
  skip_if_not(file.symlink(outside, fs::path(dir, "decisions.csv")), "this filesystem cannot make a symlink")
  writeLines(
    '{"schema": "tidymarkers/artefact@2", "id": "x", "files": {"decisions": "decisions.csv"}}',
    fs::path(dir, "artefact.json")
  )

  # The name is one component but the bytes are not the record's: both the
  # reader and the crate writer refuse the link instead of following it.
  expect_error(marker_artefact_read(dir), "resolves to", class = "rlang_error")
  art <- structure(
    list(id = "x", dir = dir, manifest = list(files = list(decisions = "decisions.csv"))),
    class = "marker_artefact"
  )
  expect_error(marker_export_rocrate(art, dir = tempfile("crate")), "resolves to", class = "rlang_error")
})

test_that("a crate descriptor cannot link out of the crate", {
  dir <- withr::local_tempdir()
  outside <- fs::path(withr::local_tempdir(), "ro-crate-metadata.json")
  writeLines('{"@context":"https://w3id.org/ro/crate/1.1/context","@graph":[]}', outside)
  skip_if_not(file.symlink(outside, fs::path(dir, "ro-crate-metadata.json")), "this filesystem cannot make a symlink")

  # The descriptor is the crate's own file; a link would make the import read
  # and vouch for metadata the crate does not carry.
  expect_error(marker_import_rocrate(dir), "resolves to", class = "rlang_error")
})

test_that("the metadata graph describes one run over the manifest inputs", {
  art <- rocrate_fixture()$artefact
  crate <- marker_export_rocrate(art, dir = tempfile("crate"))
  meta <- rocrate_meta(crate)

  expect_identical(meta[["@context"]], "https://w3id.org/ro/crate/1.1/context")
  ids <- vapply(meta[["@graph"]], function(entry) as.character(entry[["@id"]]), character(1))
  expect_true(all(c("ro-crate-metadata.json", "./", "#run", "#agent", "#howto") %in% ids))

  descriptor <- rocrate_entry(meta, "ro-crate-metadata.json")
  expect_identical(as.character(descriptor[["conformsTo"]][["@id"]]), "https://w3id.org/ro/crate/1.1")
  expect_identical(as.character(descriptor[["about"]][["@id"]]), "./")

  root <- rocrate_entry(meta, "./")
  expect_identical(as.character(root[["identifier"]]), art$id)
  expect_setequal(
    vapply(root[["hasPart"]], function(part) as.character(part[["@id"]]), character(1)),
    rocrate_manifest_files(art)
  )

  run <- rocrate_entry(meta, "#run")
  expect_true(all(c("CreateAction", "http://www.w3.org/ns/prov#Activity") %in% as.character(run[["@type"]])))
  expect_identical(as.character(run[["http://www.w3.org/ns/prov#wasAssociatedWith"]][["@id"]]), "#agent")
  expect_identical(as.character(run[["instrument"]][["@id"]]), "#howto")
  expect_setequal(
    vapply(run[["result"]], function(part) as.character(part[["@id"]]), character(1)),
    rocrate_manifest_files(art)
  )
  used <- vapply(run[["http://www.w3.org/ns/prov#used"]], function(part) as.character(part[["@id"]]), character(1))
  expect_identical(used, paste0("#input-", seq_along(art$manifest$inputs)))

  layers <- as.character(art$manifest$plan$layers)
  howto <- rocrate_entry(meta, "#howto")
  steps <- vapply(howto[["step"]], function(step) as.character(step[["name"]]), character(1))
  expect_identical(steps, layers)
})

test_that("every File entity hashes the file the crate carries", {
  art <- rocrate_fixture()$artefact
  crate <- marker_export_rocrate(art, dir = tempfile("crate"))
  meta <- rocrate_meta(crate)

  for (file in rocrate_manifest_files(art)) {
    entry <- rocrate_entry(meta, file)
    expect_false(is.null(entry), info = file)
    expect_identical(as.character(entry[["sha256"]]), hash_file(fs::path(crate, file)))
    expect_identical(as.character(entry[["http://www.w3.org/ns/prov#wasGeneratedBy"]][["@id"]]), "#run")
  }
})

test_that("a crate round-trips back to the same artefact", {
  art <- rocrate_fixture()$artefact
  crate <- marker_export_rocrate(art, dir = tempfile("crate"))

  # The reference is the artefact as read back from its own directory: the
  # in-memory manifest has richer types than the JSON round-trip preserves.
  ref <- marker_artefact_read(art$dir)
  back <- marker_import_rocrate(crate)
  expect_s3_class(back, "marker_artefact")
  expect_identical(back$id, art$id)
  expect_equal(back$decisions, ref$decisions, ignore_attr = TRUE)
  expect_identical(back$manifest$plan$layers, ref$manifest$plan$layers)
  expect_identical(back$manifest$counts, ref$manifest$counts)
  expect_identical(back$manifest$inputs, ref$manifest$inputs)

  # the declared verdicts travel too, so an imported crate replays without its RDS
  expect_identical(fs::path_file(back$manifest$files$replay), "replay.csv")
  expect_true(fs::file_exists(fs::path(back$dir, "replay.csv")))

  # The descriptor inside a crate is also a valid entry point.
  from_meta <- marker_import_rocrate(fs::path(crate, "ro-crate-metadata.json"))
  expect_identical(from_meta$id, art$id)
})

test_that("a crate cannot declare a file outside its own directory", {
  art <- rocrate_fixture()$artefact
  crate <- marker_export_rocrate(art, dir = tempfile("crate"))

  # A crate is the record it carries: an id that leaves the directory would let
  # an import hash - and thus vouch for - a file the crate does not contain.
  outside <- fs::path(fs::path_dir(crate), "outside-the-crate.txt")
  writeLines("not part of the record", outside)
  meta <- rocrate_meta(crate)
  meta[["@graph"]][[length(meta[["@graph"]]) + 1L]] <- list(
    "@id" = paste0("../", fs::path_file(outside)),
    "@type" = "File",
    "sha256" = tidymarkers:::hash_file(outside),
    "encodingFormat" = "text/plain"
  )
  writeLines(
    jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA),
    fs::path(crate, "ro-crate-metadata.json")
  )
  expect_error(marker_import_rocrate(crate), "outside", class = "rlang_error")
})

test_that("a crate cannot carry a link out of itself", {
  art <- rocrate_fixture()$artefact
  crate <- marker_export_rocrate(art, dir = tempfile("crate"))

  # The id is inside the crate but the bytes are not: a link would let the
  # import vouch for whatever the link happens to point at.
  outside <- fs::path(fs::path_dir(crate), "linked-outside.txt")
  writeLines("not part of the record", outside)
  member <- fs::path(crate, "decisions.csv")
  fs::file_delete(member)
  skip_if_not(file.symlink(outside, member), "this filesystem cannot make a symlink")
  expect_error(marker_import_rocrate(crate), "resolves to", class = "rlang_error")
})

test_that("a tampered sidecar is named rather than read", {
  art <- rocrate_fixture()$artefact
  crate <- marker_export_rocrate(art, dir = tempfile("crate"))

  path <- fs::path(crate, "decisions.csv")
  bytes <- readBin(path, "raw", n = file.size(path))
  bytes[[1L]] <- as.raw(bitwXor(as.integer(bytes[[1L]]), 1L))
  writeBin(bytes, path)

  expect_error(marker_import_rocrate(crate), "decisions\\.csv", class = "rlang_error")
})

test_that("a crate whose root identifier disagrees with artefact.json is refused", {
  art <- rocrate_fixture()$artefact
  crate <- marker_export_rocrate(art, dir = tempfile("crate"))
  meta <- rocrate_meta(crate)

  for (i in seq_along(meta[["@graph"]])) {
    if (identical(as.character(meta[["@graph"]][[i]][["@id"]]), "./")) {
      meta[["@graph"]][[i]][["identifier"]] <- "not-this-artefact"
    }
  }
  writeLines(
    jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA),
    fs::path(crate, "ro-crate-metadata.json")
  )

  expect_error(marker_import_rocrate(crate), "identifier", class = "rlang_error")
})

test_that("export refuses a non-artefact and the artefact's own directory", {
  art <- rocrate_fixture()$artefact
  expect_error(marker_export_rocrate(list()), class = "rlang_error")
  expect_error(marker_export_rocrate(art, dir = art$dir), class = "rlang_error")
})
