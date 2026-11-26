# =============================================================================
#  dev/duckdb-aggregation-bench.R  --  ZX5.2: dplyr vs duckplyr on one aggregate
# =============================================================================
#  ADR 0003 defers DuckDB/duckplyr until "the same aggregation run through
#  duckplyr and dplyr, with wall time and peak memory recorded in results/, and
#  the two answers compared on the same input".  This script is that
#  measurement.  It does not change the decision: duckdb stays outside
#  DESCRIPTION (dev-only if ever adopted).
#
#  There is no population-scale table in inst/extdata (the largest, the Steele
#  rice KASP calls, is one row per marker: 3366 rows), so the aggregate runs on
#  a scale ladder built by replicating that real schema to 1e5 / 1e6 / 1e7 rows,
#  with the real table as the anchor row.  The ladder is synthetic; the schema
#  and the query are not.
#
#  Aggregate, identical in both engines:
#    group_by(assay_id) -> sum(ref_calls, alt_calls, het_calls, missing_calls,
#    n_called), mean(call_rate), arrange(assay_id)
#
#  One child Rscript per (engine x scale) so each cell gets its own VmHWM.  The
#  parent times the child end to end and compares the answers across engines
#  with all.equal.  `elapsed_s` is the aggregation call alone (the data frame
#  build is common to both engines and is excluded); `peak_rss_mb` is the whole
#  child, so it does include the frame each engine must hold.
#
#  Note: both engines receive an in-memory data frame, so this measures query
#  cost, not an out-of-core CSV scan.  duckdb may download extensions on its
#  first run in a session (network); `duckdb(shared_home = TRUE)` keeps the
#  extension cache under ~/.duckdb instead of a fresh temp home.
#
#  Writes results/duckdb-aggregation.csv (results/ is git-ignored; the numbers
#  that must survive in the tree go into ROADMAP.md and ADR 0003).
#
#  Run:  Rscript dev/duckdb-aggregation-bench.R [--check]
#    --check  exit non-zero if the CSV is missing or a measured row is malformed
# =============================================================================

# ---- plumbing ----------------------------------------------------------------

bench_root <- function() {
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))))
  }
  getwd()
}

root <- bench_root()
setwd(root)
source(file.path(root, "dev", "proc-mem.R"))

argv <- commandArgs(trailingOnly = TRUE)
out_csv <- file.path(root, "results", "duckdb-aggregation.csv")
real_csv <- file.path(root, "inst", "extdata", "steele2025_rice_kasp_calls.csv")
scales <- c(3366L, 1e5L, 1e6L, 1e7L)
engines <- c("dplyr", "duckplyr")
columns <- c("scale_rows", "engine", "status", "elapsed_s", "peak_rss_mb", "answer_md5", "cause")

# ---- the aggregate -----------------------------------------------------------

read_real <- function() {
  d <- utils::read.csv(real_csv, stringsAsFactors = FALSE)
  d$assay_id <- as.character(d$assay_id)
  d
}

replicate_rows <- function(d, rows) {
  idx <- rep(seq_len(nrow(d)), length.out = as.integer(rows))
  out <- d[idx, , drop = FALSE]
  rownames(out) <- NULL
  out
}

duck_aggregate <- function(d, engine) {
  if (identical(engine, "dplyr")) {
    grouped <- dplyr::group_by(d, .data$assay_id)
    agg <- dplyr::summarise(
      grouped,
      ref_calls = sum(.data$ref_calls),
      alt_calls = sum(.data$alt_calls),
      het_calls = sum(.data$het_calls),
      missing_calls = sum(.data$missing_calls),
      n_called = sum(.data$n_called),
      call_rate = mean(.data$call_rate),
      .groups = "drop"
    )
    out <- dplyr::arrange(agg, .data$assay_id)
    as.data.frame(out, stringsAsFactors = FALSE)
  } else {
    tbl <- duckplyr::duckdb_tibble(d)
    grouped <- dplyr::group_by(tbl, .data$assay_id)
    agg <- dplyr::summarise(
      grouped,
      ref_calls = sum(.data$ref_calls),
      alt_calls = sum(.data$alt_calls),
      het_calls = sum(.data$het_calls),
      missing_calls = sum(.data$missing_calls),
      n_called = sum(.data$n_called),
      call_rate = mean(.data$call_rate)
    )
    out <- dplyr::arrange(agg, .data$assay_id)
    as.data.frame(out, stringsAsFactors = FALSE)
  }
}

# ---- child mode: one (engine x scale) cell -----------------------------------

if (length(argv) && identical(argv[[1L]], "--child")) {
  engine <- argv[[2L]]
  rows <- as.integer(argv[[3L]])
  outdir <- argv[[4L]]
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  suppressPackageStartupMessages(loadNamespace("dplyr"))
  if (identical(engine, "duckplyr")) {
    suppressPackageStartupMessages(loadNamespace("duckplyr"))
  }

  d <- replicate_rows(read_real(), rows)
  gc()
  t0 <- proc.time()[["elapsed"]]
  answer <- duck_aggregate(d, engine)
  elapsed <- proc.time()[["elapsed"]] - t0
  peak_kb <- proc_peak_rss_kb()

  utils::write.csv(answer, file.path(outdir, "answer.csv"), row.names = FALSE)
  utils::write.csv(
    data.frame(elapsed_s = elapsed, peak_rss_kb = peak_kb, scale_rows = nrow(d)),
    file.path(outdir, "stats.csv"),
    row.names = FALSE
  )
  quit(save = "no", status = 0L)
}

# ---- check mode --------------------------------------------------------------

if ("--check" %in% argv) {
  if (!file.exists(out_csv)) stop("missing ", out_csv, "; run without --check first", call. = FALSE)
  tab <- utils::read.csv(out_csv, stringsAsFactors = FALSE)
  missing_cols <- setdiff(columns, names(tab))
  if (length(missing_cols)) {
    stop("results CSV is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  measured <- tab$status == "measured"
  if (any(!is.finite(tab$elapsed_s[measured])) || any(!is.finite(tab$peak_rss_mb[measured]))) {
    stop("a measured row has a non-finite elapsed_s or peak_rss_mb", call. = FALSE)
  }
  if (any(!nzchar(tab$answer_md5[measured]))) {
    stop("a measured row has an empty answer_md5", call. = FALSE)
  }
  bad <- tab[!measured & is.na(tab$cause), ]
  if (nrow(bad)) stop(nrow(bad), " non-measured row(s) without a cause", call. = FALSE)
  cat("results/duckdb-aggregation.csv: ", nrow(tab), " rows, ", sum(measured), " measured\n", sep = "")
  quit(save = "no", status = 0L)
}

# ---- parent: run every cell --------------------------------------------------

if (!requireNamespace("dplyr", quietly = TRUE)) stop("dev/duckdb-aggregation-bench.R needs dplyr")
have_duck <- requireNamespace("duckplyr", quietly = TRUE) &&
  requireNamespace("duckdb", quietly = TRUE)
if (!have_duck) {
  cause <- "duckplyr/duckdb not installed (requireNamespace FALSE)"
  grid <- expand.grid(scale_rows = scales, engine = engines, stringsAsFactors = FALSE)
  tab <- data.frame(
    scale_rows = grid$scale_rows, engine = grid$engine, status = "not_run",
    elapsed_s = NA_real_, peak_rss_mb = NA_real_, answer_md5 = NA_character_, cause = cause
  )
  dir.create(file.path(root, "results"), showWarnings = FALSE)
  utils::write.csv(tab, out_csv, row.names = FALSE)
  cat("duckdb/duckplyr unavailable: wrote ", out_csv, " with status=not_run\n", sep = "")
  cat("install decision is an operator action (ADR 0003), not a pipeline step\n")
  quit(save = "no", status = 0L)
}

rscript <- file.path(R.home("bin"), "Rscript")
self <- normalizePath(file.path(root, "dev", "duckdb-aggregation-bench.R"))
work <- file.path(tempdir(), "duckdb-bench")
dir.create(work, showWarnings = FALSE, recursive = TRUE)

rows <- list()
answers <- list()
for (scale in scales) {
  for (engine in engines) {
    cell <- file.path(work, paste(engine, scale, sep = "-"))
    dir.create(cell, showWarnings = FALSE, recursive = TRUE)
    t0 <- proc.time()[["elapsed"]]
    status <- system2(
      rscript,
      c("--vanilla", shQuote(self), "--child", engine, as.character(scale), shQuote(cell)),
      stdout = file.path(cell, "child.log"), stderr = file.path(cell, "child.err")
    )
    total_s <- proc.time()[["elapsed"]] - t0
    stats_path <- file.path(cell, "stats.csv")
    answer_path <- file.path(cell, "answer.csv")
    ok <- identical(status, 0L) && file.exists(stats_path) && file.exists(answer_path)
    if (ok) {
      st <- utils::read.csv(stats_path, stringsAsFactors = FALSE)
      answers[[paste(engine, scale, sep = "-")]] <- utils::read.csv(answer_path, stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(
        scale_rows = st$scale_rows, engine = engine, status = "measured",
        elapsed_s = round(st$elapsed_s, 3), peak_rss_mb = round(st$peak_rss_kb / 1024, 1),
        answer_md5 = unname(tools::md5sum(answer_path)), cause = NA_character_
      )
      cat(sprintf(
        "  %-8s %9d rows  agg %7.2f s  child %7.2f s  peak %s\n",
        engine, st$scale_rows, st$elapsed_s, total_s, fmt_peak_rss(st$peak_rss_kb)
      ))
    } else {
      err <- paste(utils::head(readLines(file.path(cell, "child.err"), warn = FALSE), 3L), collapse = " | ")
      rows[[length(rows) + 1L]] <- data.frame(
        scale_rows = scale, engine = engine, status = "failed",
        elapsed_s = NA_real_, peak_rss_mb = NA_real_, answer_md5 = NA_character_,
        cause = paste0("child exit ", status, if (nzchar(err)) paste0(": ", err) else "")
      )
      cat(sprintf("  %-8s %9d rows  FAILED (exit %s)\n", engine, scale, status))
    }
  }
}

tab <- do.call(rbind, rows)
tab <- tab[order(tab$scale_rows, tab$engine), columns]

# Same input, same answer: compare the two engines per scale before reporting.
agree <- logical()
for (scale in unique(tab$scale_rows)) {
  key_a <- paste("dplyr", scale, sep = "-")
  key_b <- paste("duckplyr", scale, sep = "-")
  if (!is.null(answers[[key_a]]) && !is.null(answers[[key_b]])) {
    cmp <- isTRUE(all.equal(answers[[key_a]], answers[[key_b]], tolerance = 1e-8))
    agree <- c(agree, cmp)
    if (!cmp) {
      tab$status[tab$scale_rows == scale] <- "failed"
      tab$cause[tab$scale_rows == scale] <- "engines disagree on the answer (all.equal FALSE)"
    }
  }
}

dir.create(file.path(root, "results"), showWarnings = FALSE)
utils::write.csv(tab, out_csv, row.names = FALSE)

cat("\ndplyr vs duckplyr, one aggregate on the Steele rice KASP call schema\n")
print(tab[, c("scale_rows", "engine", "status", "elapsed_s", "peak_rss_mb")], row.names = FALSE)
cat("\nanswers agree across engines on every compared scale: ", all(agree), "\n", sep = "")
cat("wrote ", out_csv, " (results/ is git-ignored; durable numbers live in ROADMAP/ADR)\n", sep = "")
