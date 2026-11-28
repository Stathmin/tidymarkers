# =============================================================================
#  ZY5: the trace schema is versioned, and a migration never drops a field
# =============================================================================
#  A trace is the flattened, publication-shaped projection of an artefact: the
#  record a pooled corpus of designs would be made of.  So it carries a version
#  of its own, and a record written by an older reader must still arrive with
#  every field it had.  These tests build a real artefact from the smallest plan
#  that needs neither OligoArrayAux nor the network, then assert the shape, the
#  round trip, the legacy migration and the abort paths.

# A synthetic poly-A reference: every candidate is doomed by the GC filter, which
# is exactly the row a decision log must never leave blank.  Memoised like the
# artefact fixture in test-plan.R, because a solve is the expensive part.
trace_fixture_solve <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("A", 600L)))
      plan <- marker_plan(genome) +
        marker_targets(chrom = "Chr1", pos = 300L, ref = "A", alt = "G") +
        marker_tech(schemes = "snp")
      value <<- suppressWarnings(
        marker_solve(plan, min_disc = 0.4, progress = FALSE)
      )
    }
    value
  }
})

# The 16 fields a trace carries, in order, spelled out rather than read back from
# the package constant: a field silently dropped from the constant must fail here.
TRACE_NAMES <- c(
  "schema", "contract", "created", "package", "platform", "engine_stamp", "inputs", "plan",
  "parameters", "overrides", "limitations", "counts", "files", "id",
  "decisions", "signals"
)

# The 13 fields an `artefact@1` manifest carried.  Enumerated so the migration
# test names every field it must keep, instead of trusting the converter.
MANIFEST_NAMES <- c(
  "schema", "created", "package", "platform", "engine_stamp", "inputs", "plan",
  "parameters", "overrides", "limitations", "counts", "files", "id"
)

# The decision table's columns, from `artefact_decisions()`; a reader of a trace
# depends on all 19, even when the record carries no row.
DECISION_COLUMNS <- c(
  "target", "chrom", "pos", "ref", "alt", "scheme", "analytical", "cost",
  "chosen", "outcome", "verdict", "risk", "n_failed", "failed", "relaxed",
  "violated", "stage", "reason", "discrimination"
)

empty_decisions <- function() {
  tidymarkers:::artefact_decisions(list(attempts = tibble::tibble()))
}

# `withr::local_tempdir()` would be cleaned when a helper returns, so the artefact
# goes to `marker_artefact()`'s own temp dir, which lasts for the session.
new_artefact <- function() {
  marker_artefact(trace_fixture_solve())
}

# --- the record --------------------------------------------------------------

test_that("a trace from an artefact carries every field and its decisions", {
  art <- new_artefact()
  tr <- marker_trace_from_artefact(art)

  expect_s3_class(tr, "marker_trace")
  expect_identical(class(tr), c("marker_trace", "list"))
  expect_named(tr, TRACE_NAMES)
  expect_identical(tr$schema, marker_trace_schema())
  expect_identical(marker_trace_schema(), "tidymarkers/trace@1")
  # the pointer an artefact manifest now carries is not a trace field: on a trace
  # `schema` already names the version, so the pointer would only disagree
  expect_false("trace_schema" %in% names(tr))

  expect_identical(tr$decisions, art$decisions)
  expect_identical(tr$id, art$manifest$id)
  expect_identical(tr$counts, art$manifest$counts)
  expect_identical(tr$plan, art$manifest$plan)
  expect_identical(tr$overrides, art$manifest$overrides)
  expect_identical(tr$limitations, art$manifest$limitations)
  expect_identical(tr$files, art$manifest$files)
})

test_that("an artefact manifest names the trace version it maps to", {
  art <- new_artefact()

  expect_identical(art$manifest$schema, "tidymarkers/artefact@2")
  expect_identical(art$manifest$trace_schema, marker_trace_schema())
  json <- jsonlite::fromJSON(
    readr::read_file(fs::path(art$dir, "artefact.json")),
    simplifyVector = FALSE
  )
  expect_identical(json$trace_schema, marker_trace_schema())
  # the id hashes a named core the pointer is not part of, so adding it must not
  # renumber an artefact
  without <- art$manifest
  without$trace_schema <- NULL
  expect_identical(tidymarkers:::artefact_id(without), art$manifest$id)
})

# The id pins the serialisation, not just its self-consistency: a change to
# `auto_unbox`/`null`/`digits` in `artefact_id()` renumbers every artefact.
test_that("the artefact id is a frozen hash of the decision core", {
  # a hand-built core: `artefact_id()` keeps these eight fields and drops the
  # rest, so this needs no solve, no genome and no temp dir
  core <- list(
    schema = "tidymarkers/artefact@2",
    contract = list(version = "tidymarkers/contract@1", scenario = "unspecified"),
    package = "tidymarkers",
    engine_stamp = list(tidymarkers = "0.0.0.9000", primer3 = NA),
    inputs = list(list(role = "reference", name = "ref.fa", sha256 = "abc", bytes = 600L)),
    plan = list(targets = "one SNP"),
    parameters = list(min_disc = 0.4),
    overrides = list()
  )

  # a frozen hash, not a value to update casually: it is only renumbered when the
  # decision core or its exact JSON serialisation changes on purpose
  expect_identical(
    tidymarkers:::artefact_id(core),
    "275eb1097b756ca3730bcb77093b6389f15d82cea300c9bee3a09486acd266cd"
  )
  # an unnamed extra field is not part of the core, so it cannot renumber the id
  expect_identical(tidymarkers:::artefact_id(c(list(noise = TRUE), core)), tidymarkers:::artefact_id(core))
})

test_that("a solved plan gives the same trace without writing an artefact", {
  solved <- trace_fixture_solve()
  art <- new_artefact()
  tr <- marker_trace_from_artefact(solved)

  expect_s3_class(tr, "marker_trace")
  expect_named(tr, TRACE_NAMES)
  expect_identical(tr$counts, art$manifest$counts)
  expect_identical(tr$decisions, art$decisions)
  expect_identical(tr$plan, art$manifest$plan)
  expect_identical(tr$overrides, art$manifest$overrides)
  # nothing was written, so there are no sidecars to name, and this plan carries
  # no signal table
  expect_null(tr$files)
  expect_null(tr$signals)
  expect_match(tr$id, "^[0-9a-f]{64}$")
})

# --- round trip --------------------------------------------------------------

test_that("a trace survives a write and a read unchanged", {
  tr <- marker_trace_from_artefact(new_artefact())
  path <- withr::local_tempfile(fileext = ".rds")

  expect_identical(marker_trace_write(tr, path), path)
  back <- marker_trace_read(path)

  expect_identical(back, tr)
  expect_named(back, TRACE_NAMES)
  expect_identical(back$schema, marker_trace_schema())
  expect_identical(back$decisions, tr$decisions)
  expect_identical(back$counts, tr$counts)
  expect_identical(back$id, tr$id)
})

# --- migration ---------------------------------------------------------------

test_that("a legacy artefact@1 record keeps every field it had", {
  art <- new_artefact()
  legacy <- art$manifest
  legacy$schema <- "tidymarkers/artefact@1"
  legacy$trace_schema <- NULL
  legacy$inputs[[1L]]$path <- "/home/someone/ref.fa"
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(legacy, path)

  tr <- marker_trace_read(path)

  expect_s3_class(tr, "marker_trace")
  expect_named(tr, TRACE_NAMES)
  expect_true(all(MANIFEST_NAMES %in% names(tr)))
  expect_false("trace_schema" %in% names(tr))
  expect_identical(tr$schema, marker_trace_schema())

  # every field the legacy record had, compared one by one: `schema` is the
  # version that changed and `inputs` is the path that had to be reduced
  for (nm in setdiff(MANIFEST_NAMES, c("schema", "inputs"))) {
    expect_identical(tr[[nm]], legacy[[nm]], info = nm)
  }
  expect_identical(tr$inputs[[1L]]$path, "ref.fa")
  expect_identical(
    tr$inputs[[1L]][c("role", "source", "name", "sha256", "bytes")],
    legacy$inputs[[1L]][c("role", "source", "name", "sha256", "bytes")]
  )

  # `artefact@1` kept its rows in decisions.csv, so a bare manifest cannot carry
  # them back: the field is present and empty, never invented and never absent
  expect_identical(tr$decisions, empty_decisions())
  expect_identical(names(tr$decisions), DECISION_COLUMNS)
  expect_identical(nrow(tr$decisions), 0L)
  expect_null(tr$signals)
})

test_that("a legacy record is written back at the current schema", {
  art <- new_artefact()
  legacy <- art$manifest
  legacy$schema <- "tidymarkers/artefact@1"
  legacy$trace_schema <- NULL
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(legacy, path)

  migrated <- marker_trace_read(path)
  out <- withr::local_tempfile(fileext = ".rds")
  marker_trace_write(migrated, out)
  back <- marker_trace_read(out)

  expect_identical(back, migrated)
  expect_identical(back$schema, marker_trace_schema())
})

# An artefact directory is the package's own on-disk record: `artefact@1` kept its
# decision rows in `decisions.csv` beside the manifest, so reading the directory
# back is what lets a migrated trace keep them instead of arriving empty.
legacy_artefact_dir <- function(art, dir) {
  manifest <- art$manifest
  manifest$schema <- "tidymarkers/artefact@1"
  manifest$trace_schema <- NULL
  manifest$inputs[[1L]]$path <- "/home/someone/ref.fa"
  writeLines(
    jsonlite::toJSON(
      manifest,
      auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA
    ),
    fs::path(dir, "artefact.json")
  )
  readr::write_csv(art$decisions, fs::path(dir, manifest$files$decisions))
  dir
}

test_that("an artefact directory migrates to a trace without losing its decisions", {
  art <- new_artefact()
  dir <- legacy_artefact_dir(art, withr::local_tempdir("legacy-artefact-"))

  back <- marker_artefact_read(dir)
  expect_s3_class(back, "marker_artefact")
  expect_named(back, c("id", "dir", "manifest", "decisions"))
  expect_identical(back$id, art$manifest$id)
  expect_identical(names(back$decisions), DECISION_COLUMNS)
  expect_identical(nrow(back$decisions), nrow(art$decisions))
  expect_true(nrow(back$decisions) > 0L)

  tr <- marker_trace_from_artefact(dir)
  expect_s3_class(tr, "marker_trace")
  expect_named(tr, TRACE_NAMES)
  expect_identical(tr$schema, marker_trace_schema())
  expect_identical(tr$id, art$manifest$id)
  # a bare manifest would give the empty table; the sidecar gives the rows back
  expect_identical(tr$decisions, back$decisions)
  expect_identical(tr$decisions$outcome, art$decisions$outcome)
  expect_identical(tr$inputs[[1L]]$path, "ref.fa")
  expect_length(tr$limitations, length(art$manifest$limitations))

  # the manifest file itself is accepted where the directory is
  expect_identical(marker_trace_from_artefact(fs::path(dir, "artefact.json")), tr)
})

test_that("an artefact reader refuses a schema it does not know and names it", {
  dir <- withr::local_tempdir("future-artefact-")
  writeLines('{"schema": "tidymarkers/artefact@3"}', fs::path(dir, "artefact.json"))

  err <- tryCatch(marker_artefact_read(dir), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "tidymarkers/artefact@3", fixed = TRUE)
  expect_match(conditionMessage(err), "Update tidymarkers", fixed = TRUE)

  err <- tryCatch(marker_trace_from_artefact(dir), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "tidymarkers/artefact@3", fixed = TRUE)

  empty <- withr::local_tempdir("empty-artefact-")
  err <- tryCatch(marker_artefact_read(empty), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "No artefact.json", fixed = TRUE)
})

test_that("migrating a trace@1 is identity and idempotent", {
  tr <- marker_trace_from_artefact(new_artefact())

  expect_identical(marker_trace_migrate(tr), tr)
  expect_identical(marker_trace_migrate(marker_trace_migrate(tr)), tr)
  # a record stored without the class keeps its shape and gains the class
  plain <- unclass(tr)
  expect_false(inherits(plain, "marker_trace"))
  expect_identical(marker_trace_migrate(plain), tr)
  expect_identical(marker_trace_migrate(marker_trace_migrate(plain)), tr)
})

test_that("an unknown or newer schema aborts and names the version", {
  err <- tryCatch(marker_trace_migrate(list(schema = "nope@1")), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "nope@1", fixed = TRUE)

  err <- tryCatch(marker_trace_migrate(list()), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "<none>", fixed = TRUE)

  err <- tryCatch(marker_trace_migrate(list(schema = "tidymarkers/trace@2")), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "tidymarkers/trace@2", fixed = TRUE)

  err <- tryCatch(marker_trace_migrate(list(schema = "tidymarkers/artefact@3")), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "tidymarkers/artefact@3", fixed = TRUE)
})

# --- the aggregation contract ------------------------------------------------

test_that("a trace carries no absolute path and no hostname", {
  tr <- marker_trace_from_artefact(new_artefact())
  flat <- as.character(unlist(tr, recursive = TRUE, use.names = FALSE))
  flat <- flat[!is.na(flat)]

  expect_false(any(grepl("^(/|[A-Za-z]:[/\\\\])", flat)))
  node <- Sys.info()[["nodename"]]
  if (!is.na(node) && nzchar(node)) {
    expect_false(any(grepl(node, flat, fixed = TRUE)))
  }
})

# --- refusals ----------------------------------------------------------------

test_that("reading something that is not a trace aborts and names the file", {
  path <- withr::local_tempfile(pattern = "not-a-trace", fileext = ".rds")
  saveRDS(list(schema = "something/else@1"), path)

  err <- tryCatch(marker_trace_read(path), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), fs::path_file(path), fixed = TRUE)
  expect_match(conditionMessage(err), "Not a trace", fixed = TRUE)

  missing <- withr::local_tempfile(pattern = "absent", fileext = ".rds")
  expect_false(fs::file_exists(missing))
  err <- tryCatch(marker_trace_read(missing), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "No trace at", fixed = TRUE)
})

test_that("writing refuses a missing directory and an incomplete trace", {
  tr <- marker_trace_from_artefact(new_artefact())
  absent <- fs::path(withr::local_tempdir(), "absent-subdirectory")

  err <- tryCatch(marker_trace_write(tr, fs::path(absent, "trace.rds")), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "No directory", fixed = TRUE)

  broken <- tr
  broken$counts <- NULL
  err <- tryCatch(
    marker_trace_write(broken, withr::local_tempfile(fileext = ".rds")),
    error = identity
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "counts", fixed = TRUE)
})

test_that("converting something that is neither an artefact nor a solve aborts", {
  err <- tryCatch(marker_trace_from_artefact(list()), error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "must be a marker_artefact", fixed = TRUE)
})
