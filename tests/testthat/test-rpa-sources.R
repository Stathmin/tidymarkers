# =============================================================================
#  RPA mismatch calibration -- installed data must match the fetched papers
# =============================================================================
#  The three CSVs this file guards are written by
#    data-raw/rpa_mismatch_class.R, data-raw/asrpa_mdr.R, data-raw/asrpa_sickle.R
#  Each generator re-reads .tmp/lit/<key>/ and stop()s on disagreement; here we
#  re-assert the same arithmetic on the installed files so a stale or
#  hand-edited CSV fails the suite without needing the papers on disk.

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("rpa_mismatch_terminal: 12 classes, probabilities and ranking", {
  d <- ext("rpa_mismatch_terminal.csv")
  expect_equal(nrow(d), 12L)
  expect_setequal(names(d)[1:5], c(
    "pattern", "terminal_mismatch", "RS_logit",
    "RS_ci_95", "RS_p_star"
  ))
  expect_true(all(grepl("^[?ACGT]{4}-[?ACGT]{4}$", d$pattern)))
  # every terminal class is the unnamed base of both 4-letter blocks
  expect_equal(substr(d$pattern, 4L, 4L), substr(d$terminal_mismatch, 1L, 1L))
  # no Watson-Crick pair is a "mismatch"
  expect_false(any(d$terminal_mismatch %in% c("A-T", "T-A", "C-G", "G-C")))
  # the paper's two verbatim probabilities
  expect_equal(round(d$RS_probability[d$terminal_mismatch == "C-C"], 2), 0.59)
  expect_equal(round(d$RS_probability[d$terminal_mismatch == "G-A"], 2), 0.62)
  # and the ranking the table supports: C-C worst, G-A next
  expect_equal(d$terminal_mismatch[which.min(d$RS_logit)], "C-C")
  expect_equal(d$terminal_mismatch[order(d$RS_logit)][2L], "G-A")
  # only A-C and T-G are unstarred on reaction success
  expect_setequal(d$terminal_mismatch[!d$RS_starred], c("A-C", "T-G"))
  # every point estimate inside its own printed CI
  num <- regmatches(d$RS_ci_95, gregexpr("-?[0-9]+\\.[0-9]+", d$RS_ci_95))
  lo <- vapply(num, function(x) as.numeric(x[1L]), numeric(1))
  hi <- vapply(num, function(x) as.numeric(x[2L]), numeric(1))
  expect_true(all(d$RS_logit >= lo & d$RS_logit <= hi))
  expect_true(all(
    c("outcome_definition", "evidence_tier", "source_file") %in% names(d)
  ))
})

test_that("rpa_mismatch_position: 279 rows from the 36 printed position classes", {
  d <- ext("rpa_mismatch_position.csv")
  expect_equal(nrow(d), 279L) # every mismatch cell of the 36 printed classes
  expect_true(all(d$terminal_offset == 0L))
  expect_true(all(d$nonterminal_offset %in% 1:3))
  expect_true(all(grepl("^[?ACGT]{4}-[?ACGT]{4}$", d$pattern)))
  # the named offset really is where the deeper concrete base sits
  blk <- strsplit(d$pattern, "-", fixed = TRUE)
  first_concrete <- vapply(blk, function(b) {
    which(strsplit(b[1L], "")[[1L]] != "?")[1L]
  }, integer(1))
  expect_equal(4L - first_concrete, d$nonterminal_offset)
})

test_that("rpa_mismatch_claims: disagreements with the table are recorded", {
  d <- ext("rpa_mismatch_claims.csv")
  expect_true(nrow(d) >= 6L)
  expect_true(all(c(
    "claim", "paper_says", "table_shows", "agrees",
    "evidence_tier"
  ) %in% names(d)))
  # the abstract's C-T claim is contradicted by the table and must say so
  top <- d$claim == "terminal mismatch most detrimental to reaction success"
  expect_match(d$note[top], "C-T")
  # 159 of 315 does reproduce from the printed stars
  k <- which(grepl("P < 0.05 on reaction success", d$claim))
  expect_true(d$agrees[k])
})

test_that("rpa_asrpa_mdr: call rates, and they are labelled as call rates", {
  d <- ext("rpa_asrpa_mdr_calls.csv")
  expect_equal(nrow(d), 6L)
  expect_true(all(d$N == 141L))
  expect_true(all(d$TP + d$FP + d$FN + d$TN == 141L))
  expect_equal(round(100 * d$TN / (d$TN + d$FP), 2), d$specificity_pct)
  combined <- grepl("combined", d$assay)
  expect_equal(sum(d$FN[combined]), 8L)
  expect_equal(sum(d$TP[combined] + d$FN[combined]), 54L)
  expect_true(all(grepl(
    "CALL AGREEMENT, not a mismatch effect",
    d$outcome_definition,
    fixed = TRUE
  )))
})

test_that("rpa_asrpa_sickle: calls against genotype and the 91-sample denominator", {
  d <- ext("rpa_asrpa_sickle_calls.csv")
  expect_equal(nrow(d), 2L)
  expect_true(all(d$n_samples == 91L))
  expect_equal(
    round(100 * d$sensitivity_n / d$sensitivity_d, 1),
    d$sensitivity_pct
  )
  expect_equal(
    round(100 * d$specificity_n / d$specificity_d, 1),
    d$specificity_pct
  )
  expect_equal(round(100 * d$accuracy_n / d$accuracy_d, 1), d$accuracy_pct)
  expect_true(all(d$n_samples_excluded == 18L))
  expect_true(all(grepl("CALL AGREEMENT", d$outcome_definition)))
})

test_that("rpa generators reproduce their CSVs when the sources are present", {
  # Run each generator in a CHILD process with `root` bound to a scratch
  # library and the fetched source directory injected, then require the file it
  # writes to be byte-identical to the one we ship. A stale or hand-edited CSV
  # fails here even though the assertions inside the generator pass.
  repo <- test_path("..", "..")
  keys <- c(
    rpa_mismatch_class = "rpa_mismatch_terminal.csv",
    asrpa_mdr = "rpa_asrpa_mdr_calls.csv",
    asrpa_sickle = "rpa_asrpa_sickle_calls.csv"
  )
  for (key in names(keys)) {
    src <- file.path(repo, ".tmp", "lit", key)
    skip_if(!dir.exists(src), paste("source not fetched:", key))
    tmp <- withr::local_tempdir()
    out_dir <- file.path(tmp, "inst", "extdata")
    dir.create(out_dir, recursive = TRUE)

    body <- readLines(file.path(repo, "data-raw", paste0(key, ".R")),
      warn = FALSE
    )
    # Rebuild the child's environment: `root` points at a scratch library, the
    # script's `--file=` names the real generator (so it resolves its own
    # helper correctly), and the missing `--dir` makes it read .tmp/lit/<key>
    # from the repo root we setwd() into.
    runner <- file.path(tmp, "run.R")
    body <- readLines(file.path(repo, "data-raw", paste0(key, ".R")),
      warn = FALSE
    )
    start <- grep(
      "^(script <- grep|suppressPackageStartupMessages)",
      body
    )[1L]
    writeLines(c(
      sprintf("root <- %s", deparse(normalizePath(tmp))),
      sprintf("file_arg <- %s", deparse(normalizePath(file.path(
        repo, "data-raw", paste0(key, ".R")
      )))),
      sprintf("dir_arg <- %s", deparse(normalizePath(src))),
      paste0(
        "commandArgs <- function(...) c(paste0(\"--file=\", file_arg), ",
        "paste0(\"--dir=\", dir_arg))"
      ),
      sprintf("setwd(%s)", deparse(normalizePath(repo))),
      body[start:length(body)]
    ), runner)

    st <- system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(runner)),
      stdout = FALSE, stderr = FALSE
    )
    expect_equal(st, 0L, label = paste("generator exit status:", key))
    fresh <- readLines(file.path(out_dir, keys[[key]]), warn = FALSE)
    shipped <- readLines(test_path("../../inst/extdata", keys[[key]]),
      warn = FALSE
    )
    expect_identical(fresh, shipped)
  }
})
