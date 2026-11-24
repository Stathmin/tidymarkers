# =============================================================================
#  ZY22: a second writer reproduces an artefact byte for byte
# =============================================================================
#  The conformance set is the four files that make a design replayable:
#  `artefact.json` plus `decisions.csv`, `checks.csv` and `why.csv`.
#  `marker_trace_payload()` hands those to `inst/python/marker_trace.py` with the
#  in-memory manifest (so an `AsIs` vector and a scalar stay distinguishable) and
#  the sidecar rows read with the schema's types.  These tests build one real
#  artefact, re-render it in Python, and assert the bytes, the replay and the
#  schema.  The schema test compares the schema against the tables the package
#  builds - not against a file - so a column added in R fails here.

# A solve is the expensive part; the hrm scheme needs OligoArrayAux, so the whole
# file skips without it (the synthetic poly-A fixture writes no checks/why).
conformance_fixture <- local({
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
      artefact <- marker_artefact(solved, dir = tempfile("conformance-ref"))
      value <<- list(solved = solved, artefact = artefact)
    }
    value
  }
})

# The schema's own type names, derived from the in-memory column.
trace_type_of <- function(x) {
  if (is.factor(x)) {
    return("string")
  }
  if (is.integer(x)) {
    return("integer")
  }
  if (is.double(x)) {
    return("number")
  }
  if (is.logical(x)) {
    return("logical")
  }
  if (is.character(x)) {
    return("string")
  }
  NA_character_
}

trace_column_types <- function(tab) {
  stats::setNames(vapply(tab, trace_type_of, character(1)), names(tab))
}

test_that("the manifest names the values jsonlite renders inline", {
  art <- conformance_fixture()$artefact
  expect_identical(trace_inline_arrays(art$manifest), c("plan.layers", "limitations"))
})

test_that("the payload carries the manifest, its inline paths and typed rows", {
  art <- conformance_fixture()$artefact
  payload <- marker_trace_payload(art)
  expect_identical(
    names(payload),
    c("schema", "trace_schema", "manifest", "inline_arrays", "tables")
  )
  expect_identical(payload$schema, "tidymarkers/trace-payload@1")
  expect_identical(payload$trace_schema, TRACE_SCHEMA)
  expect_identical(payload$manifest, art$manifest)
  expect_identical(unclass(payload$inline_arrays), c("plan.layers", "limitations"))
  expect_identical(names(payload$tables), c("decisions", "checks", "why"))

  # A cell is a one-element list, so `auto_unbox` keeps every row an array and a
  # missing value stays `null` rather than vanishing from the row.
  first_row <- payload$tables$decisions$rows[[1]]
  expect_true(all(vapply(first_row, function(cell) is.list(cell) && length(cell) == 1L, logical(1))))
  expect_identical(payload$tables$decisions$columns, trace_conformance_schema()$tables[[1]]$columns)

  # A directory has already lost the atomic-vs-list distinction the payload needs.
  expect_error(marker_trace_payload(art$dir), class = "rlang_error")
})

test_that("the conformance schema names the columns and types the package builds", {
  fixture <- conformance_fixture()
  solved <- fixture$solved
  schema <- trace_conformance_schema()
  checks <- artefact_checks(solved)
  tables <- list(
    decisions = list(columns = artefact_decisions(solved), schema = schema$tables[[1]]),
    checks = list(columns = checks, schema = schema$tables[[2]]),
    why = list(columns = marker_explain(checks), schema = schema$tables[[3]])
  )
  for (name in names(tables)) {
    want_names <- vapply(tables[[name]]$schema$columns, function(col) col$name, character(1))
    want_types <- vapply(tables[[name]]$schema$columns, function(col) col$type, character(1))
    expect_identical(names(tables[[name]]$columns), want_names, info = name)
    expect_identical(trace_column_types(tables[[name]]$columns), stats::setNames(want_types, want_names), info = name)
  }
})

test_that("the Python writer reproduces the four conformance files byte for byte", {
  skip_if(Sys.which("python3") == "", "python3 is not available")
  fixture <- conformance_fixture()
  art <- fixture$artefact
  script <- system.file("python", "marker_trace.py", package = "tidymarkers")
  expect_true(nzchar(script))

  payload <- marker_trace_payload(art)
  payload_path <- tempfile("payload", fileext = ".json")
  writeLines(
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", na = "null", digits = 17),
    payload_path
  )
  outdir <- tempfile("conformance-py")
  status <- system2(
    "python3",
    c(script, payload_path, outdir, "--solved", fs::path(art$dir, "solved.rds")),
    stdout = TRUE, stderr = TRUE
  )
  expect_identical(attr(status, "status"), NULL, info = paste(status, collapse = "\n"))

  bytes <- function(path) readBin(path, "raw", n = file.size(path))
  for (file in c("artefact.json", "decisions.csv", "checks.csv", "why.csv")) {
    expect_identical(
      bytes(fs::path(outdir, file)),
      bytes(fs::path(art$dir, file)),
      info = file
    )
  }

  # The copy is not the part being trusted: replay re-solves from `solved.rds`.
  replayed <- marker_replay(outdir, progress = FALSE)
  reference <- marker_replay(art$dir, progress = FALSE)
  expect_true(replayed$match)
  expect_identical(replayed$table, reference$table)
})

test_that("the Python writer refuses a payload that diverges from the schema", {
  skip_if(Sys.which("python3") == "", "python3 is not available")
  art <- conformance_fixture()$artefact
  script <- system.file("python", "marker_trace.py", package = "tidymarkers")

  payload <- marker_trace_payload(art)
  payload$tables$decisions$columns[[1]]$type <- "date"
  payload_path <- tempfile("payload", fileext = ".json")
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", na = "null", digits = 17), payload_path)
  outdir <- tempfile("conformance-bad")

  divergent <- suppressWarnings(
    system2("python3", c(script, payload_path, outdir), stdout = TRUE, stderr = TRUE)
  )
  expect_identical(attr(divergent, "status"), 1L)
  expect_false(fs::dir_exists(outdir))

  missing_schema <- suppressWarnings(system2(
    "python3",
    c(script, payload_path, outdir, "--schema", tempfile("absent")),
    stdout = TRUE, stderr = TRUE
  ))
  expect_identical(attr(missing_schema, "status"), 2L)
})
