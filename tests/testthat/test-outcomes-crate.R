# =============================================================================
#  ZY27: outcomes handed to a trustee without a network
# =============================================================================
#  `marker_export_outcomes()` takes a corpus and writes one RO-Crate per
#  artefact, each a real artefact crate with the recorded rows attached, plus an
#  identity map beside the export.  These tests check the three promises the
#  function makes: the crate round-trips through `marker_import_rocrate()` (so
#  the hashes still agree with the files), no absolute path and no artefact
#  identity travels inside a crate, and the pseudonym is reproducible from the
#  salt but not from anything the crate carries.
#
#  There is no official RO-Crate validator in this tree, so the JSON-LD is read
#  structurally and the round-trip is a hash check.  That is a structural check,
#  not a validation.

outcomes_corpus_fixture <- local({
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
      root <- tempfile("corpus")
      a <- marker_artefact(solved, dir = file.path(root, "batch-a"))
      b <- marker_artefact(solved, dir = file.path(root, "batch-b"))
      marker_record_outcome(a, "works")
      marker_record_outcome(b, "fails", notes = "no amplification")
      dir <- paste0(tempfile("handover"), "/out")
      map <- marker_export_outcomes(root, dir, trustee = "lab-a")
      value <<- list(root = root, dir = dir, map = map, corpus = marker_outcome_corpus(root))
    }
    value
  }
})

test_that("one crate per artefact, and each one imports to its own id", {
  fx <- outcomes_corpus_fixture()
  crates <- fs::dir_ls(fs::path(fx$dir, "crates"))
  expect_length(crates, 2L)
  expect_true(all(startsWith(fs::path_file(crates), "site-")))
  ids <- unname(vapply(crates, function(crate) marker_import_rocrate(crate)$id, character(1)))
  expect_setequal(ids, fx$map$artefact_id)
})

test_that("a crate carries the rows and no absolute path", {
  fx <- outcomes_corpus_fixture()
  crate <- fs::path(fx$dir, "crates", fx$map$pseudonym[[1L]])
  raw <- readr::read_file(fs::path(crate, "ro-crate-metadata.json"))
  expect_false(grepl(as.character(fs::path_abs(fx$root)), raw, fixed = TRUE))

  meta <- rocrate_meta(crate)
  urls <- vapply(
    meta[["@graph"]],
    function(entry) as.character(entry[["contentUrl"]] %||% NA_character_),
    character(1)
  )
  expect_false(any(fs::is_absolute_path(urls[!is.na(urls)])))

  rows <- fx$corpus[fx$corpus$artefact == fx$map$site[[1L]], , drop = FALSE]
  outcome_ids <- vapply(meta[["@graph"]], function(entry) {
    id <- as.character(entry[["@id"]])
    startsWith(id, "#outcome-") && !identical(id, "#outcome-index")
  }, logical(1))
  expect_length(which(outcome_ids), nrow(rows))

  first <- rocrate_entry(meta, "#outcome-1")
  expect_identical(as.character(first[["schema"]]), "tidymarkers/outcomes@1")
  expect_identical(as.character(first[["result"]]), rows$result[[1L]])
  expect_identical(as.character(first[[paste0(PROV_NS, "wasDerivedFrom")]][["@id"]]), "./")
  expect_null(first[["notes"]])
  expect_null(first[["target"]])
})

test_that("the copied outcomes.csv is hashed in the crate", {
  fx <- outcomes_corpus_fixture()
  crate <- fs::path(fx$dir, "crates", fx$map$pseudonym[[1L]])
  entity <- rocrate_entry(rocrate_meta(crate), "outcomes.csv")
  expect_identical(
    as.character(entity[["sha256"]]),
    hash_file(fs::path(crate, "outcomes.csv"))
  )
  root <- rocrate_entry(rocrate_meta(crate), "./")
  parts <- vapply(root[["hasPart"]], function(part) as.character(part[["@id"]]), character(1))
  expect_true("outcomes.csv" %in% parts)
})

test_that("the identity map stays outside the crate and maps back", {
  fx <- outcomes_corpus_fixture()
  path <- paste0(fx$dir, "-identity.csv")
  expect_true(fs::file_exists(path))
  map <- utils::read.csv(path, colClasses = "character")
  expect_named(map, c("site", "pseudonym", "artefact_id"))
  expect_setequal(map$site, c("batch-a", "batch-b"))
  inside <- fs::dir_ls(fs::path(fx$dir, "crates"), recurse = TRUE, type = "file")
  expect_false(any(grepl("identity", fs::path_file(inside), fixed = TRUE)))
  for (i in seq_len(nrow(map))) {
    rel <- map$site[[i]]
    expect_false(identical(map$pseudonym[[i]], fs::path_file(rel)))
    expect_true(fs::dir_exists(fs::path(fx$dir, "crates", map$pseudonym[[i]])))
    expect_identical(
      artefact_manifest(fs::path(fx$root, rel))$id,
      map$artefact_id[[i]]
    )
  }
})

test_that("the handover refuses to write its identity map through a link", {
  fx <- outcomes_corpus_fixture()
  dir <- paste0(withr::local_tempdir(), "/out")
  outside <- fs::path(withr::local_tempdir(), "identity.csv")
  writeLines("keep", outside)
  skip_if_not(
    file.symlink(outside, paste0(dir, "-identity.csv")),
    "this filesystem cannot make a symlink"
  )

  expect_error(
    marker_export_outcomes(fx$corpus, dir, trustee = "lab-a"),
    "which is a link"
  )
  expect_identical(readLines(outside), "keep")
  expect_false(fs::dir_exists(fs::path(dir, "crates")))
})

test_that("the pseudonym is reproducible from the salt alone", {
  fx <- outcomes_corpus_fixture()
  expected <- paste0("site-", substr(hash_object("lab-abatch-a"), 1L, 8L))
  expect_identical(fx$map$pseudonym[fx$map$site == "batch-a"], expected)

  other <- marker_export_outcomes(fx$corpus, paste0(tempfile("handover2"), "/out"), trustee = "lab-b")
  expect_false(any(other$pseudonym %in% fx$map$pseudonym))
  again <- marker_export_outcomes(fx$corpus, paste0(tempfile("handover3"), "/out"), trustee = "lab-a")
  expect_identical(again$pseudonym, fx$map$pseudonym)
})

test_that("a corpus row cannot name a directory outside the scanned root", {
  fx <- outcomes_corpus_fixture()
  tampered <- fx$corpus
  tampered$artefact[[1L]] <- "../outside-the-corpus"
  expect_error(
    marker_export_outcomes(tampered, paste0(tempfile("handover"), "/out")),
    "outside its scanned root"
  )
})

test_that("a corpus without a scanned root or an export inside it is refused", {
  fx <- outcomes_corpus_fixture()
  hand_built <- tibble::tibble(artefact = "x", artefact_id = "y")
  expect_error(
    marker_export_outcomes(hand_built, tempfile("handover")),
    "no scanned root"
  )
  expect_error(
    marker_export_outcomes(fx$root, fs::path(fx$root, "out")),
    "outside the scanned corpus root"
  )
})
