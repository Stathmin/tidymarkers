# =============================================================================
#  ZY25  --  what was rejected, across batches
# =============================================================================
#  The registry is a view over sidecars an artefact already writes, so the
#  fixture is hand-written CSVs: no solve, no genome, nothing to install.
#  What matters is that a rejection survives the trip -- a run where every
#  chosen marker scored clean still had candidates that did not make it, and
#  `why.csv` cannot see them.

fixture_batch <- function(root, name, why = NULL, checks = NULL, decisions = NULL,
                          manifest = NULL) {
  dir <- file.path(root, name)
  dir.create(dir, recursive = TRUE)
  if (!is.null(checks)) {
    utils::write.csv(checks, file.path(dir, "checks.csv"), row.names = FALSE)
  }
  if (!is.null(why)) {
    utils::write.csv(why, file.path(dir, "why.csv"), row.names = FALSE)
  }
  if (!is.null(decisions)) {
    utils::write.csv(decisions, file.path(dir, "decisions.csv"), row.names = FALSE)
  }
  if (!is.null(manifest)) {
    writeLines(
      as.character(jsonlite::toJSON(manifest, auto_unbox = TRUE)),
      file.path(dir, "artefact.json")
    )
  }
  dir
}

why_row <- function(target, criterion, status, reason) {
  data.frame(
    target = target, criterion = criterion, status = status, reason = reason,
    stringsAsFactors = FALSE
  )
}

checks_row <- function(target, criterion, scheme, pass = FALSE, status = "fail") {
  data.frame(
    target = target, scheme = scheme, stage = "design", criterion = criterion,
    observed = "1", limit = "0", pass = pass, status = status,
    detail = "", level = "design", stringsAsFactors = FALSE
  )
}

dec_row <- function(target, scheme, outcome, failed = NA_character_,
                    reason = NA_character_) {
  data.frame(
    target = target, scheme = scheme, outcome = outcome, failed = failed,
    reason = reason, stringsAsFactors = FALSE
  )
}

test_that("ZY25: a run with clean chosen markers still reports its rejections", {
  root <- withr::local_tempdir()
  fixture_batch(root, "batch-a",
    why = why_row("t1", "gc", "pass", NA),
    checks = checks_row("t1", "gc", "hrm", pass = TRUE, status = "pass"),
    decisions = rbind(
      dec_row("t2", "hrm", "rejected", "tm_gap; gc", "no signal window"),
      dec_row("t3", "hrm", "chosen")
    )
  )
  reg <- marker_failure_registry(root)
  # why.csv has zero non-pass rows; the two rejected rules are the whole table
  expect_equal(nrow(reg), 2L)
  expect_setequal(reg$rule_id, c("tm_gap", "gc"))
  expect_true(all(reg$candidate_id == "batch-a"))
  expect_false(any(reg$target == "t3"))
})

test_that("ZY25: two artefacts of one scheme union without duplicates", {
  root <- withr::local_tempdir()
  for (b in c("batch-a", "batch-b")) {
    fixture_batch(root, b,
      decisions = dec_row("t1", "hrm", "rejected", "gc", "gc too high")
    )
  }
  reg <- marker_failure_registry(root)
  expect_equal(nrow(reg), 2L)
  expect_setequal(reg$candidate_id, c("batch-a", "batch-b"))
  expect_equal(unique(reg$scheme), "hrm")
})

test_that("ZY25: a missing reason and a missing rule are named, not blank", {
  root <- withr::local_tempdir()
  fixture_batch(root, "batch-a",
    decisions = dec_row("t9", "dpcr", "rejected", NA_character_, NA_character_)
  )
  reg <- marker_failure_registry(root)
  expect_equal(nrow(reg), 1L)
  expect_equal(reg$rule_id, "unrecorded")
  expect_equal(reg$reason, "unrecorded")
  expect_false(any(is.na(reg$rule_id) | is.na(reg$reason)))
})

test_that("ZY25: row count is why-failures plus rejected rules", {
  root <- withr::local_tempdir()
  fixture_batch(root, "batch-a",
    why = rbind(
      why_row("t1", "gc", "fail", "gc: 71 is outside 40-65"),
      why_row("t1", "tm_gap", "unknown", "Tm gap not measurable"),
      why_row("t1", "homopolymer", "pass", NA)
    ),
    checks = rbind(
      checks_row("t1", "gc", "lamp"),
      checks_row("t1", "tm_gap", "lamp", pass = TRUE, status = "unknown"),
      checks_row("t1", "homopolymer", "lamp", pass = TRUE, status = "pass")
    ),
    decisions = rbind(
      dec_row("t2", "lamp", "rejected", "gc; duplex", "both fired"),
      dec_row("t3", "lamp", "designable"),
      dec_row("t4", "lamp", "chosen")
    )
  )
  reg <- marker_failure_registry(root)
  # 2 non-pass why rows (the pass row is not a rejection) + 2 rejected rules
  expect_equal(nrow(reg), 4L)
  expect_equal(sum(reg$source_file == "batch-a/why.csv"), 2L)
  expect_equal(sum(reg$source_file == "batch-a/decisions.csv"), 2L)
})

test_that("ZY25: source_file is relative and source_line lands on the row", {
  root <- withr::local_tempdir()
  batch <- fixture_batch(root, "batch-a",
    decisions = rbind(
      dec_row("t1", "hrm", "chosen"),
      dec_row("t2", "hrm", "rejected", "gc", "high")
    )
  )
  reg <- marker_failure_registry(root)
  expect_equal(reg$source_file, "batch-a/decisions.csv")
  expect_false(any(grepl(root, reg$source_file, fixed = TRUE)))
  lines <- readLines(file.path(batch, "decisions.csv"))
  # header is line 1, so the second data row is line 3 and it names target t2
  expect_equal(reg$source_line, 3L)
  expect_match(lines[[3L]], "t2")
})

test_that("ZY25: the manifest id names the batch when one exists", {
  root <- withr::local_tempdir()
  fixture_batch(root, "batch-a",
    decisions = dec_row("t1", "hrm", "rejected", "gc", "high"),
    manifest = list(id = "cand-42")
  )
  expect_equal(marker_failure_registry(root)$candidate_id, "cand-42")
})

test_that("ZY25: the summary counts rules and batches, not rows twice", {
  root <- withr::local_tempdir()
  fixture_batch(root, "batch-a",
    decisions = rbind(
      dec_row("t1", "hrm", "rejected", "gc", "high"),
      dec_row("t2", "hrm", "rejected", "gc; tm_gap", "wide")
    )
  )
  fixture_batch(root, "batch-b",
    decisions = dec_row("t1", "hrm", "rejected", "gc", "high again")
  )
  summ <- marker_failure_summary(marker_failure_registry(root))
  expect_equal(summ$rule_id[1], "gc")
  expect_equal(summ$n[summ$rule_id == "gc"], 3L)
  expect_equal(summ$candidates[summ$rule_id == "gc"], 2L)
  expect_equal(summ$n[summ$rule_id == "tm_gap"], 1L)
  # the frequent rule sorts first
  expect_true(which(summ$rule_id == "gc") < min(which(summ$rule_id == "tm_gap")))
})

test_that("ZY25: an empty tree yields the declared columns, not an error", {
  root <- withr::local_tempdir()
  reg <- marker_failure_registry(root)
  expect_equal(nrow(reg), 0L)
  expect_equal(names(reg), failure_columns())
  expect_equal(nrow(marker_failure_summary(reg)), 0L)
  expect_error(marker_failure_registry(file.path(root, "absent")), "No directory")
})
