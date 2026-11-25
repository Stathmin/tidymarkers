# =============================================================================
#  dev/success-model.R  --  which checks predict an accepted construction?
# =============================================================================
#  The perspectives draft closes with "a success model would need ...".  This
#  harness makes the first half of that sentence runnable: on the outcome corpus
#  built by `marker_outcome_corpus()` it asks whether any single design check
#  separates constructions a laboratory accepted from ones it rejected, and
#  whether a simple sum of the checks does better than each one alone.  It is the
#  multi-axis counterpart of `dev/single-score-baseline.R`: that scalar did not
#  separate the published panel, and this run asks which axis, if any, carries
#  the separation.
#
#  Unit of observation: one recorded outcome, i.e. a row of an `outcomes.csv`,
#  joined to the artefact's `checks.csv` by `target` -- or by the whole artefact
#  when the outcome names no target.  A check counts as passed for that row only
#  when every row written for that criterion passed; a fail anywhere is a fail;
#  anything else (`unknown`, `not_applicable`) leaves that check unused for the
#  row.  A row whose `result` is `"unknown"` settles nothing and is dropped from
#  the corpus before anything is counted.
#
#  Metric: AUC of a check's pass/fail against `result == "works"`.  0.5 is no
#  separation, 1.0 is a perfect one, and a tie counts one half.  A leave-one-out
#  mean and a min-max range sit beside it: one row can carry an AUC, so the
#  spread is what says whether the number is worth quoting.
#
#  Sum predictor: the count of passed checks, restricted to the checks present in
#  every surviving row.  Without that restriction a row with more evaluated
#  checks could outscore a better design by arithmetic alone.
#
#  Threshold: below SUCCESS_MIN_ROWS recorded outcomes, or below
#  SUCCESS_MIN_PER_CLASS of either class, a predictor is `not_run` with its
#  reason rather than a number.  An empty corpus is one `not_run` row, never a
#  `0` that reads as "measured and bad".
#
#  Blocker, stated in the output and not hidden: the tree holds no
#  `marker_record_outcome()` rows, so the default corpus is empty and every
#  predictor is `not_run`.  The blocker bars training, not the harness: the
#  arithmetic is exercised by a synthetic corpus in a temporary directory.
#
#  Run:  Rscript dev/success-model.R [--corpus=DIR] [--check]
#        --corpus=DIR  directory searched for outcomes.csv (default .tmp/corpus)
#        --check       recompute, and fail unless the table reproduces and the
#                      block in VALIDATION.md matches it; never writes
#
#  Writes the marked block of `## Success model` in VALIDATION.md (between
#  `<!-- success-model:begin -->` and `<!-- success-model:end -->`) and the same
#  table to .tmp/success-model.csv (git-ignored).
#
#  Verdicts: `read` (the corpus row), `measured` (an AUC was computed), and
#  `not_run` (a reason instead of a number).
# =============================================================================

# ---- plumbing ----------------------------------------------------------------

success_root <- function() {
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))))
  }
  wd <- getwd()
  repeat {
    if (file.exists(file.path(wd, "DESCRIPTION")) && dir.exists(file.path(wd, "R"))) {
      return(wd)
    }
    parent <- dirname(wd)
    if (identical(parent, wd)) {
      return(getwd())
    }
    wd <- parent
  }
}

success_path <- function(...) file.path(success_root(), ...)

# A path as the committed document should read it: relative to the repository
# when it lives there, and only otherwise as given.
success_display_path <- function(path) {
  rel <- try(fs::path_rel(path, success_root()), silent = TRUE)
  if (inherits(rel, "try-error") || grepl("^(\\.\\.|/)", as.character(rel))) {
    return(as.character(path))
  }
  as.character(rel)
}

SUCCESS_MIN_ROWS <- 12L
SUCCESS_MIN_PER_CLASS <- 4L
SUCCESS_DEFAULT_CORPUS <- ".tmp/corpus"
SUCCESS_DOC <- "VALIDATION.md"
SUCCESS_BEGIN <- "<!-- success-model:begin -->"
SUCCESS_END <- "<!-- success-model:end -->"

success_args <- function(args) {
  corpus <- grep("^--corpus=", args, value = TRUE)
  corpus <- if (length(corpus)) sub("^--corpus=", "", corpus[[1L]]) else SUCCESS_DEFAULT_CORPUS
  if (!grepl("^(/|~|[A-Za-z]:)", corpus)) {
    corpus <- success_path(corpus)
  }
  list(corpus = path.expand(corpus), check = "--check" %in% args)
}

success_corpus <- function(dir) {
  if (!fs::dir_exists(dir)) {
    return(structure(
      tibble::add_column(
        outcome_empty(),
        artefact_id = character(0), artefact = character(0),
        .before = 1
      ),
      root = as.character(fs::path_abs(dir))
    ))
  }
  marker_outcome_corpus(dir)
}

# ---- the corpus as rows ------------------------------------------------------

# One row per recorded outcome with a settled result; `unknown` and a missing
# result are dropped together, because neither says whether the design worked.
success_labels <- function(corpus) {
  result <- as.character(corpus$result)
  keep <- result %in% c("works", "fails")
  tibble::tibble(
    record = which(keep),
    label = as.integer(result[keep] == "works"),
    target = as.character(corpus$target[keep]),
    artefact = as.character(corpus$artefact[keep])
  )
}

success_checks_file <- function(path) {
  if (!fs::file_exists(path)) {
    return(NULL)
  }
  tbl <- utils::read.csv(path, colClasses = "character", stringsAsFactors = FALSE)
  if (!all(c("criterion", "status") %in% names(tbl))) {
    return(NULL)
  }
  tbl
}

# The value of every check for one recorded outcome: 1 passed, 0 failed, NA when
# the check was not evaluated for that row.
success_record_values <- function(root, rec) {
  tbl <- success_checks_file(file.path(root, rec$artefact, "checks.csv"))
  if (is.null(tbl)) {
    return(NULL)
  }
  if ("target" %in% names(tbl) && !is.na(rec$target)) {
    tbl <- tbl[tbl$target == rec$target, , drop = FALSE]
  }
  if (!nrow(tbl)) {
    return(NULL)
  }
  crit <- sort(unique(tbl$criterion))
  value <- vapply(crit, function(one) {
    status <- tbl$status[tbl$criterion == one]
    if (any(status == "fail")) {
      0L
    } else if (all(status == "pass")) {
      1L
    } else {
      NA_integer_
    }
  }, integer(1))
  tibble::tibble(record = rec$record, check = crit, value = unname(value))
}

# Records x checks, with NA where the check was not evaluated.
success_matrix <- function(root, labels) {
  parts <- lapply(seq_len(nrow(labels)), function(i) success_record_values(root, labels[i, ]))
  long <- dplyr::bind_rows(parts[!vapply(parts, is.null, logical(1))])
  if (!nrow(long)) {
    return(NULL)
  }
  checks <- sort(unique(long$check))
  m <- matrix(
    NA_integer_,
    nrow = nrow(labels), ncol = length(checks),
    dimnames = list(NULL, checks)
  )
  m[cbind(match(long$record, labels$record), match(long$check, checks))] <- long$value
  m
}

# ---- the metric --------------------------------------------------------------

success_auc <- function(value, label) {
  ok <- !is.na(value) & !is.na(label)
  value <- value[ok]
  label <- label[ok]
  n_works <- sum(label == 1L)
  n_fails <- sum(label == 0L)
  if (!n_works || !n_fails) {
    return(NA_real_)
  }
  concordant <- 0
  tied <- 0
  for (v in value[label == 1L]) {
    concordant <- concordant + sum(v > value[label == 0L])
    tied <- tied + sum(v == value[label == 0L])
  }
  (concordant + 0.5 * tied) / (n_works * n_fails)
}

success_loo <- function(value, label) {
  none <- c(mean = NA_real_, min = NA_real_, max = NA_real_)
  if (length(value) < 3L) {
    return(none)
  }
  one <- vapply(seq_along(value), function(i) success_auc(value[-i], label[-i]), numeric(1))
  one <- one[!is.na(one)]
  if (!length(one)) {
    return(none)
  }
  c(mean = mean(one), min = min(one), max = max(one))
}

success_predictor <- function(name, value, label, note = "") {
  value <- as.numeric(value)
  label <- as.integer(label)
  used <- !is.na(value)
  n <- sum(used)
  works <- sum(label[used] == 1L)
  fails <- sum(label[used] == 0L)
  reason <- if (n < SUCCESS_MIN_ROWS) {
    sprintf("%d recorded outcome(s) used, below the %d-row threshold", n, SUCCESS_MIN_ROWS)
  } else if (works < SUCCESS_MIN_PER_CLASS || fails < SUCCESS_MIN_PER_CLASS) {
    sprintf(
      "%d works and %d fails used, below %d per class",
      works, fails, SUCCESS_MIN_PER_CLASS
    )
  } else {
    ""
  }
  if (nzchar(reason)) {
    return(tibble::tibble(
      predictor = name, n = n, works = works, fails = fails,
      auc = NA_real_, loo_mean = NA_real_, loo_min = NA_real_, loo_max = NA_real_,
      verdict = "not_run", note = trimws(paste(reason, note))
    ))
  }
  auc <- success_auc(value, label)
  loo <- success_loo(value, label)
  tibble::tibble(
    predictor = name, n = n, works = works, fails = fails,
    auc = auc,
    loo_mean = unname(loo[["mean"]]),
    loo_min = unname(loo[["min"]]),
    loo_max = unname(loo[["max"]]),
    verdict = if (is.na(auc)) "not_run" else "measured",
    note = if (is.na(auc)) trimws(paste("only one class in the used rows", note)) else note
  )
}

# ---- the table ---------------------------------------------------------------

success_rows <- function(corpus) {
  root <- attr(corpus, "root")
  labels <- success_labels(corpus)
  dropped <- nrow(corpus) - nrow(labels)
  corpus_row <- tibble::tibble(
    predictor = "success-model:corpus",
    n = nrow(labels),
    works = sum(labels$label == 1L),
    fails = sum(labels$label == 0L),
    auc = NA_real_, loo_mean = NA_real_, loo_min = NA_real_, loo_max = NA_real_,
    verdict = if (nrow(labels) == 0L) "not_run" else "read",
    note = if (dropped > 0L) sprintf("%d dropped as unknown or unrecorded", dropped) else ""
  )
  no_corpus <- tibble::tibble(
    predictor = "success-model:corpus", n = 0L, works = 0L, fails = 0L,
    auc = NA_real_, loo_mean = NA_real_, loo_min = NA_real_, loo_max = NA_real_,
    verdict = "not_run",
    note = sprintf(
      "no outcomes.csv under %s; record some with marker_record_outcome()",
      success_display_path(root)
    )
  )
  if (nrow(corpus) == 0L) {
    return(no_corpus)
  }
  if (nrow(labels) == 0L) {
    return(corpus_row)
  }
  m <- success_matrix(root, labels)
  if (is.null(m)) {
    return(dplyr::bind_rows(corpus_row, tibble::tibble(
      predictor = "success-model:checks", n = 0L, works = 0L, fails = 0L,
      auc = NA_real_, loo_mean = NA_real_, loo_min = NA_real_, loo_max = NA_real_,
      verdict = "not_run",
      note = "no checks.csv beside the recorded outcomes, so no check can be scored"
    )))
  }
  per <- lapply(colnames(m), function(one) success_predictor(one, m[, one], labels$label))
  complete <- colnames(m)[colSums(is.na(m)) == 0L]
  sum_row <- if (length(complete)) {
    success_predictor(
      "success-model:sum", rowSums(m[, complete, drop = FALSE]), labels$label,
      note = sprintf("%d check(s) present in every recorded row", length(complete))
    )
  } else {
    tibble::tibble(
      predictor = "success-model:sum", n = 0L, works = 0L, fails = 0L,
      auc = NA_real_, loo_mean = NA_real_, loo_min = NA_real_, loo_max = NA_real_,
      verdict = "not_run",
      note = "no check is present in every recorded row, so a sum is not comparable"
    )
  }
  dplyr::bind_rows(corpus_row, per, sum_row)
}

# ---- rendering ---------------------------------------------------------------

success_fmt <- function(x) {
  ifelse(is.na(x), "\u2014", formatC(x, format = "f", digits = 3L))
}

success_range <- function(low, high) {
  bad <- is.na(low) | is.na(high)
  ifelse(bad, "\u2014", paste0(success_fmt(low), "\u2013", success_fmt(high)))
}

success_block <- function(x) {
  measured <- sum(x$verdict == "measured")
  lead <- if (x$n[[1L]] == 0L) {
    paste0(
      "_The corpus holds no recorded outcome, so no predictor can be scored. ",
      "This run reports `not_run` with its reason rather than a number._"
    )
  } else {
    paste0(
      sprintf(
        "_The corpus holds %d recorded outcome(s): %d works, %d fails. ",
        x$n[[1L]], x$works[[1L]], x$fails[[1L]]
      ),
      sprintf(
        "%d predictor(s) carry both classes in enough rows and are reported as AUC; ",
        measured
      ),
      "the rest are `not_run` with their reason._"
    )
  }
  head <- c("Predictor", "n", "works", "fails", "AUC", "LOO mean", "LOO min-max", "Verdict", "Reason")
  body <- cbind(
    paste0("`", x$predictor, "`"),
    as.character(x$n), as.character(x$works), as.character(x$fails),
    success_fmt(x$auc), success_fmt(x$loo_mean), success_range(x$loo_min, x$loo_max),
    paste0("`", x$verdict, "`"),
    ifelse(nzchar(x$note), x$note, "\u2014")
  )
  # Pad every column to its widest cell (at least three characters) so the
  # generated block matches what a markdown formatter leaves behind; otherwise
  # `--check` fails the first time an editor saves the file.
  cells <- rbind(head, body)
  widths <- pmax(3L, apply(nchar(cells, type = "width"), 2L, max))
  row_md <- function(v) {
    pad <- strrep(" ", widths - nchar(v, type = "width"))
    paste0("| ", paste0(v, pad, collapse = " | "), " |")
  }
  rule <- paste0("| ", paste0(strrep("-", widths), collapse = " | "), " |")
  c(
    "", lead, "", row_md(head), rule, apply(body, 1L, row_md), "",
    paste0(
      "The harness is `dev/success-model.R`; its header names the join rule, ",
      "the tie convention and the thresholds above."
    )
  )
}

success_current_block <- function(lines) {
  begin <- which(lines == SUCCESS_BEGIN)
  end <- which(lines == SUCCESS_END)
  if (length(begin) != 1L || length(end) != 1L || end - begin <= 1L) {
    return(character(0))
  }
  lines[(begin + 1L):(end - 1L)]
}

success_replace <- function(lines, block) {
  begin <- which(lines == SUCCESS_BEGIN)
  end <- which(lines == SUCCESS_END)
  if (length(begin) != 1L || length(end) != 1L || end < begin) {
    stop(sprintf(
      "%s needs exactly one %s / %s pair to hold the generated table",
      SUCCESS_DOC, SUCCESS_BEGIN, SUCCESS_END
    ), call. = FALSE)
  }
  c(lines[seq_len(begin)], block, lines[end:length(lines)])
}

# ---- driver ------------------------------------------------------------------

main <- function() {
  opt <- success_args(commandArgs(trailingOnly = TRUE))
  if (!"tidymarkers" %in% loadedNamespaces()) {
    suppressMessages(pkgload::load_all(success_root(), quiet = TRUE))
  }
  x <- success_rows(success_corpus(opt$corpus))
  block <- success_block(x)
  doc <- success_path(SUCCESS_DOC)
  lines <- readLines(doc, warn = FALSE, encoding = "UTF-8")
  message(sprintf(
    "success model: %d row(s); %d measured, %d not_run",
    nrow(x), sum(x$verdict == "measured"), sum(x$verdict == "not_run")
  ))
  for (i in seq_len(nrow(x))) {
    message(sprintf("  %s: %s", x$predictor[[i]], x$verdict[[i]]))
  }
  if (opt$check) {
    again <- success_block(success_rows(success_corpus(opt$corpus)))
    if (!identical(again, block)) {
      stop("--check: the AUC table moved between two runs", call. = FALSE)
    }
    if (!identical(success_current_block(lines), block)) {
      stop(
        "--check: ", SUCCESS_DOC, " does not match a fresh run; rerun without --check",
        call. = FALSE
      )
    }
    message("--check: reproduced")
    return(invisible(x))
  }
  csv <- success_path(".tmp", "success-model.csv")
  dir.create(dirname(csv), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(x, csv, row.names = FALSE, na = "")
  writeLines(success_replace(lines, block), doc, useBytes = TRUE)
  message("wrote ", csv, " and the marked block of ", SUCCESS_DOC)
  invisible(x)
}

if (sys.nframe() == 0L) {
  main()
}
