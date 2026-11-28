# =============================================================================
#  ZY4  --  an external corpus is queued, never merged
# =============================================================================
#  The registry only grows through marker_ingest_approve().  These tests pin the
#  two halves of that promise: ingestion sorts every source row into candidates
#  or rejected (no row is silently dropped), and a read never writes a file.
#  The thresholds exercised here are local policy, not a literature claim.

zy4_corpus <- function() {
  testthat::test_path("fixtures", "ingest-candidates.csv")
}

test_that("ingest sorts every source row into candidates or rejected", {
  q <- marker_ingest_literature(zy4_corpus(), min_confidence = 0.5)

  expect_s3_class(q, "marker_ingest_queue")
  expect_equal(nrow(q$candidates) + nrow(q$rejected), 6L)
  expect_equal(nrow(q$candidates), 2L)
  expect_equal(nrow(q$rejected), 4L)

  # the clean row and the first copy of the duplicate pair survive, in file order
  expect_equal(q$candidates$record_id, c("lit-001", "lit-005"))
  expect_equal(q$candidates$value, c(0.85, 0.55))
  expect_equal(q$candidates$doi, c("10.1000/aaa001", "10.1000/aaa005"))
  expect_identical(q$candidates$evidence, c(3L, 2L))
  expect_equal(q$candidates$confidence, c(0.9, 0.8))

  # the other four are kept with the first rule each one broke
  expect_equal(q$rejected$record_id, c("lit-002", "lit-003", "lit-004", "lit-005"))
  expect_match(q$rejected$reason[q$rejected$record_id == "lit-002"], "provenance")
  expect_match(q$rejected$reason[q$rejected$record_id == "lit-003"], "min_evidence")
  expect_match(q$rejected$reason[q$rejected$record_id == "lit-004"], "min_confidence")
  expect_match(q$rejected$reason[q$rejected$record_id == "lit-005"], "duplicate record_id")
  expect_identical(names(q$rejected)[[length(names(q$rejected))]], "reason")

  # the queue records what it applied and where it read from
  expect_identical(q$criteria$min_evidence, 1L)
  expect_identical(q$criteria$min_confidence, 0.5)
  expect_true(q$criteria$provenance_required)
  expect_identical(q$source, as.character(fs::path_norm(fs::path_abs(zy4_corpus()))))
  expect_equal(nrow(q$approved), 0L)
})

test_that("reading a registry that does not exist yields zero rows and writes nothing", {
  registry <- fs::path(withr::local_tempdir(), "registry.csv")

  reg <- marker_ingest_registry(registry)

  expect_equal(nrow(reg), 0L)
  expect_identical(
    names(reg),
    c(
      "record_id", "doi", "metric", "value", "evidence", "confidence", "notes",
      "approved_by", "approved_at"
    )
  )
  expect_false(fs::file_exists(registry))
})

test_that("approval refuses to write the registry through a link", {
  dir <- withr::local_tempdir()
  outside <- fs::path(dir, "outside.csv")
  readr::write_lines("keep", outside)
  registry <- fs::path(dir, "registry.csv")
  skip_if(!file.symlink(outside, registry), "symlinks are not supported")
  q <- marker_ingest_literature(zy4_corpus(), min_confidence = 0.5)

  expect_error(
    marker_ingest_approve(q, "lit-001", registry = registry),
    "which is a link"
  )
  expect_identical(readLines(outside), "keep")
})

test_that("approval is the only write path, it appends one stamped record", {
  registry <- fs::path(withr::local_tempdir(), "deep", "nested", "registry.csv")
  q <- marker_ingest_literature(zy4_corpus(), min_confidence = 0.5)

  approved <- withVisible(marker_ingest_approve(q, "lit-001", registry = registry))
  expect_false(approved$visible)
  q2 <- approved$value

  reg <- marker_ingest_registry(registry)
  expect_equal(nrow(reg), 1L)
  expect_equal(reg$record_id, "lit-001")
  expect_equal(reg$doi, "10.1000/aaa001")
  expect_equal(reg$metric, "LoD")
  expect_equal(reg$value, 0.85)
  expect_identical(reg$evidence, 3L)
  expect_equal(reg$confidence, 0.9)
  expect_equal(reg$notes, "clean row")
  expect_equal(reg$approved_by, "maintainer")
  expect_match(reg$approved_at, "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$")

  # the returned queue moved the row; the input queue is untouched
  expect_equal(q2$candidates$record_id, "lit-005")
  expect_equal(q2$approved$record_id, "lit-001")
  expect_equal(q$candidates$record_id, c("lit-001", "lit-005"))
  expect_equal(nrow(q$approved), 0L)

  # a second approval appends instead of replacing
  q3 <- marker_ingest_approve(q2, "lit-005", registry = registry)
  expect_equal(marker_ingest_registry(registry)$record_id, c("lit-001", "lit-005"))
  expect_equal(q3$approved$record_id, c("lit-001", "lit-005"))
  expect_equal(nrow(q3$candidates), 0L)

  # approving the same id again is an error, and nothing is written
  expect_error(
    marker_ingest_approve(q, "lit-001", registry = registry),
    "already in the registry"
  )
  expect_equal(nrow(marker_ingest_registry(registry)), 2L)

  expect_error(
    marker_ingest_approve(q, "lit-999", registry = registry),
    "not among this queue's candidates"
  )
  expect_error(marker_ingest_approve(list(), "lit-001", registry = registry), "marker_ingest_literature")
  expect_error(marker_ingest_approve(q, "", registry = registry), "record_id")
})

test_that("the JSON path sorts rows exactly like the CSV path", {
  tmp <- withr::local_tempdir()
  json <- fs::path(tmp, "corpus.json")
  dat <- readr::read_csv(
    zy4_corpus(),
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA"), progress = FALSE
  )
  jsonlite::write_json(dat, json, auto_unbox = TRUE, na = "null", pretty = TRUE)

  csv_queue <- marker_ingest_literature(zy4_corpus(), min_confidence = 0.5)
  json_queue <- marker_ingest_literature(json, min_confidence = 0.5)

  expect_equal(nrow(json_queue$candidates) + nrow(json_queue$rejected), 6L)
  expect_equal(json_queue$candidates, csv_queue$candidates)
  expect_equal(json_queue$rejected, csv_queue$rejected)

  # a top-level object is not a corpus
  object <- fs::path(tmp, "object.json")
  readr::write_lines('{"record_id": "r1"}', object)
  expect_error(marker_ingest_literature(object), "top-level array")
})

test_that("provenance_required = FALSE admits a doi-less row, leaving its doi NA", {
  q <- marker_ingest_literature(zy4_corpus(), min_confidence = 0.5, provenance_required = FALSE)

  expect_equal(q$candidates$record_id, c("lit-001", "lit-002", "lit-005"))
  doi <- q$candidates$doi[q$candidates$record_id == "lit-002"]
  expect_length(doi, 1L)
  expect_true(is.na(doi))
  expect_false("lit-002" %in% q$rejected$record_id)
  expect_false(q$criteria$provenance_required)
})

test_that("a malformed cell is rejected with the rule that failed", {
  tmp <- withr::local_tempdir()
  corpus <- fs::path(tmp, "corpus.csv")
  readr::write_lines(
    c(
      "record_id,doi,metric,value,evidence,confidence",
      "a,10.1000/aaa,LoD,high,2,0.9",
      "b,10.1000/bbb,LoD,0.9,two,0.9",
      "c,10.1000/ccc,LoD,0.9,2,1.2",
      "d,10.1000/ddd,LoD,0.9,,",
      "e,10.1000/eee,LoD,0.9,2,0.9"
    ),
    corpus
  )
  q <- marker_ingest_literature(corpus)

  expect_equal(q$rejected$record_id, c("a", "b", "c"))
  expect_match(q$rejected$reason[[1L]], "malformed value")
  expect_match(q$rejected$reason[[2L]], "malformed evidence")
  expect_match(q$rejected$reason[[3L]], "malformed confidence")

  # a blank evidence takes 1; a blank confidence is not filtered on
  expect_equal(q$candidates$record_id, c("d", "e"))
  expect_identical(q$candidates$evidence, c(1L, 2L))
  expect_true(is.na(q$candidates$confidence[[1L]]))
})

test_that("a clean source rejects nothing but still carries the reason column", {
  tmp <- withr::local_tempdir()
  corpus <- fs::path(tmp, "corpus.csv")
  readr::write_lines(c("record_id,doi,metric,value", "a,10.1000/aaa,LoD,0.9"), corpus)

  q <- marker_ingest_literature(corpus)

  expect_equal(nrow(q$rejected), 0L)
  expect_equal(
    names(q$rejected),
    c("record_id", "doi", "metric", "value", "evidence", "confidence", "notes", "reason")
  )
  expect_equal(nrow(q$candidates), 1L)
  expect_identical(q$candidates$evidence, 1L)
})

test_that("a bad source path or threshold is an error, not an empty queue", {
  tmp <- withr::local_tempdir()
  corpus <- fs::path(tmp, "corpus.csv")
  readr::write_lines(c("record_id,doi,metric", "a,10.1000/aaa,LoD"), corpus)

  expect_error(
    marker_ingest_literature(corpus),
    "`literature source` is missing column(s): value",
    fixed = TRUE
  )
  expect_error(marker_ingest_literature(fs::path(tmp, "absent.csv")), "is not one")
  txt <- fs::path(tmp, "corpus.txt")
  readr::write_lines("x", txt)
  expect_error(marker_ingest_literature(txt), "must end in .csv or .json")
  expect_error(marker_ingest_literature(zy4_corpus(), min_confidence = 2), "min_confidence")
  expect_error(marker_ingest_literature(zy4_corpus(), min_evidence = -1), "min_evidence")
  # a count of Inf is refused by the criterion, not accepted as "no bar"
  expect_error(marker_ingest_literature(zy4_corpus(), min_evidence = Inf), "min_evidence")
  expect_error(marker_ingest_literature(zy4_corpus(), provenance_required = NA), "provenance_required")
})
