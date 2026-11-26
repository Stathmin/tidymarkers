# =============================================================================
#  marker_failure_registry.R  --  what was rejected, across batches (ZY25)
# =============================================================================
#  An artefact records its own rejections: `decisions.csv` names every
#  (target, chemistry) attempt that did not make it, and `why.csv` names the
#  checks a chosen marker did not pass.  Both are per-run, so "what did we
#  reject, and why" has no answer that survives two batches.  This file is that
#  answer, and it is deliberately a view: it reads the sidecars an artefact
#  already writes and invents no field the run did not produce.
#
#  Two sources, because one is not enough.  `why.csv` only covers chosen
#  markers (`artefact_checks()` keeps the ones that passed the gate), so a run
#  where every chosen marker scored clean has no failing row there, while
#  `decisions.csv` still holds every rejected candidate.  A registry built from
#  `why.csv` alone would report an empty corpus on exactly the runs that
#  rejected the most.  `checks.csv` is read for one reason: it is the only
#  sidecar that carries `scheme` for a `why.csv` row.
#
#  A row is a (candidate, target, scheme, rule) rejection.  `candidate_id` names
#  the batch - the artefact's own `id`, or the directory basename when no
#  manifest was written - because a rejected attempt has no id of its own.
#  `source_file` is relative to the scanned root and `source_line` is the data
#  row's line number (the header is line 1), so a reader can open the sidecar
#  and land on the row.
# =============================================================================

#' The columns of the failure registry
#' @keywords internal
#' @noRd
failure_columns <- function() {
  c(
    "candidate_id", "target", "scheme", "rule_id", "reason",
    "source_file", "source_line"
  )
}

#' @keywords internal
#' @noRd
failure_empty <- function() {
  tibble::tibble(
    candidate_id = character(0), target = character(0), scheme = character(0),
    rule_id = character(0), reason = character(0), source_file = character(0),
    source_line = integer(0)
  )
}

#' The batch name of a sidecar's directory
#' @keywords internal
#' @noRd
failure_candidate_id <- function(dir) {
  man <- tryCatch(artefact_manifest(dir), error = function(e) NULL)
  id <- man$id %||% NULL
  if (is.character(id) && length(id) == 1L && !is.na(id) && nzchar(id)) {
    return(id)
  }
  basename(dir)
}

#' A sidecar as a character table, or `NULL` when it is absent
#' @keywords internal
#' @noRd
failure_sidecar <- function(path) {
  if (!fs::file_exists(path)) {
    return(NULL)
  }
  tryCatch(
    readr::read_csv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      na = c("", "NA"), show_col_types = FALSE, progress = FALSE
    ),
    error = function(e) NULL
  )
}

#' One sidecar cell, `NA` when the column or the row is absent
#' @keywords internal
#' @noRd
failure_cell <- function(tbl, name, i) {
  if (is.null(tbl) || !name %in% names(tbl)) {
    return(NA_character_)
  }
  as.character(tbl[[name]][[i]])
}

#' Split a `; `-joined rule list into its rules
#' @keywords internal
#' @noRd
failure_rules <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }
  rules <- trimws(strsplit(x, ";", fixed = TRUE)[[1L]])
  rules <- rules[nzchar(rules)]
  if (!length(rules)) NA_character_ else rules
}

#' A missing reason or rule is named, never left blank
#' @keywords internal
#' @noRd
failure_named <- function(x) {
  blank <- is.na(x) | !nzchar(x)
  x[blank] <- "unrecorded"
  x
}

#' Non-`pass` rows of one `why.csv`
#' @keywords internal
#' @noRd
failure_from_why <- function(why, checks, candidate_id, rel) {
  if (is.null(why) || !nrow(why)) {
    return(failure_empty())
  }
  col <- function(name) {
    vapply(seq_len(nrow(why)), function(i) failure_cell(why, name, i), character(1))
  }
  status <- col("status")
  keep <- !is.na(status) & status != "pass" & nzchar(status)
  if (!any(keep)) {
    return(failure_empty())
  }
  # `why.csv` has no `scheme`; `checks.csv` is the sidecar that does, and its
  # (target, criterion) pair is unique per chosen marker
  scheme_key <- NULL
  if (!is.null(checks) && nrow(checks) && all(c("target", "criterion", "scheme") %in% names(checks))) {
    scheme_key <- stats::setNames(
      as.character(checks$scheme),
      paste(checks$target, checks$criterion, sep = "\r")
    )
  }
  idx <- which(keep)
  target <- col("target")
  criterion <- col("criterion")
  reason <- col("reason")
  tibble::tibble(
    candidate_id = candidate_id,
    target = target[idx],
    scheme = if (is.null(scheme_key)) {
      rep(NA_character_, length(idx))
    } else {
      unname(scheme_key[paste(target[idx], criterion[idx], sep = "\r")])
    },
    rule_id = criterion[idx],
    reason = reason[idx],
    source_file = rel,
    source_line = as.integer(idx + 1L)
  )
}

#' Rejected rows of one `decisions.csv`, one row per failed rule
#' @keywords internal
#' @noRd
failure_from_decisions <- function(dec, candidate_id, rel) {
  if (is.null(dec) || !nrow(dec) || !"outcome" %in% names(dec)) {
    return(failure_empty())
  }
  outcome <- as.character(dec$outcome)
  idx <- which(!is.na(outcome) & outcome == "rejected")
  if (!length(idx)) {
    return(failure_empty())
  }
  purrr::map(idx, function(i) {
    rules <- failure_rules(failure_cell(dec, "failed", i))
    tibble::tibble(
      candidate_id = candidate_id,
      target = failure_cell(dec, "target", i),
      scheme = failure_cell(dec, "scheme", i),
      rule_id = rules,
      reason = failure_cell(dec, "reason", i),
      source_file = rel,
      source_line = as.integer(i + 1L)
    )
  }) |>
    purrr::list_rbind()
}

#' Every rejection the artefacts under a directory recorded
#'
#' Reads the sidecars an artefact already writes and binds them into one table,
#' so two batches are comparable without a manual diff.  Both sources are used:
#' non-`pass` rows of `why.csv` (a chosen marker whose check did not pass) and
#' `outcome = "rejected"` rows of `decisions.csv` (a candidate that never made
#' it).  A rejected candidate that failed several rules contributes one row per
#' rule, because a count per rule is the question the summary answers.
#'
#' A row is never dropped for being incomplete: an empty `reason` and a missing
#' `rule_id` are written as `"unrecorded"` rather than as a blank cell, so a
#' missing explanation is visible instead of invisible.
#' @param dir the directory to scan (default: the working directory).
#' @param recursive scan subdirectories?  A corpus of batched artefacts usually
#'   needs `TRUE` (the default).
#' @return a tibble with `candidate_id` (the artefact's `id`, else its directory
#'   name), `target`, `scheme`, `rule_id`, `reason`, `source_file` (relative to
#'   `dir`) and `source_line`.
#' @examples
#' root <- tempfile("batches")
#' batch <- file.path(root, "batch-a")
#' dir.create(batch, recursive = TRUE)
#' utils::write.csv(
#'   data.frame(
#'     target = "t1", scheme = "hrm", criterion = "gc",
#'     observed = 71, limit = "40-65", pass = FALSE, status = "fail"
#'   ),
#'   file.path(batch, "checks.csv"),
#'   row.names = FALSE
#' )
#' utils::write.csv(
#'   data.frame(
#'     target = "t1", criterion = "gc", status = "fail",
#'     reason = "gc: 71 is outside 40-65"
#'   ),
#'   file.path(batch, "why.csv"),
#'   row.names = FALSE
#' )
#' utils::write.csv(
#'   data.frame(
#'     target = c("t2", "t3"), scheme = "hrm",
#'     outcome = c("rejected", "chosen"),
#'     failed = c("tm_gap; gc", NA), reason = c("no signal window", NA)
#'   ),
#'   file.path(batch, "decisions.csv"),
#'   row.names = FALSE
#' )
#' reg <- marker_failure_registry(root)
#' reg
#' marker_failure_summary(reg)
#' @seealso [marker_artefact()], [marker_explain()]
#' @family plan
#' @export
marker_failure_registry <- function(dir = ".", recursive = TRUE) {
  root <- fs::path_abs(dir)
  if (!fs::dir_exists(root)) {
    cli::cli_abort(gettextf("No directory at %s.", shQuote(as.character(dir))))
  }
  files <- fs::dir_ls(root, recurse = recursive, type = "file", fail = FALSE)
  files <- as.character(files)
  files <- files[grepl("(^|/)(why|decisions)\\.csv$", files)]
  files <- sort(files)
  parts <- purrr::map(files, function(path) {
    here <- dirname(path)
    rel <- as.character(fs::path_rel(path, start = root))
    id <- failure_candidate_id(here)
    if (basename(path) == "why.csv") {
      failure_from_why(
        failure_sidecar(path),
        failure_sidecar(fs::path(here, "checks.csv")),
        id, rel
      )
    } else {
      failure_from_decisions(failure_sidecar(path), id, rel)
    }
  })
  out <- if (length(parts)) purrr::list_rbind(parts) else failure_empty()
  if (!nrow(out)) {
    return(out)
  }
  out$rule_id <- failure_named(out$rule_id)
  out$reason <- failure_named(out$reason)
  out <- dplyr::distinct(out)
  out[order(out$candidate_id, out$target, out$scheme, out$rule_id, out$source_line), , drop = FALSE]
}

#' How often each rule rejected a candidate
#'
#' One row per rule, most frequent first.  `candidates` counts the batches that
#' rejected on it, so a rule that fires in ten batches is not confused with ten
#' rejections of one batch.
#' @param registry a table from [marker_failure_registry()].
#' @return a tibble with `rule_id`, `n`, `candidates` and `schemes`.
#' @examples
#' root <- tempfile("batches")
#' batch <- file.path(root, "batch-a")
#' dir.create(batch, recursive = TRUE)
#' utils::write.csv(
#'   data.frame(
#'     target = "t1", scheme = "hrm", outcome = "rejected",
#'     failed = "gc; tm_gap", reason = "no signal window"
#'   ),
#'   file.path(batch, "decisions.csv"),
#'   row.names = FALSE
#' )
#' marker_failure_summary(marker_failure_registry(root))
#' @seealso [marker_failure_registry()]
#' @family plan
#' @export
marker_failure_summary <- function(registry) {
  if (!is.data.frame(registry)) {
    abort_msg("`registry` must be a table from marker_failure_registry().")
  }
  if (!nrow(registry)) {
    return(tibble::tibble(
      rule_id = character(0), n = integer(0),
      candidates = integer(0), schemes = character(0)
    ))
  }
  rules <- sort(unique(as.character(registry$rule_id)))
  purrr::map(rules, function(rule) {
    hit <- registry$rule_id == rule
    schemes <- sort(unique(as.character(registry$scheme[hit])))
    schemes <- schemes[!is.na(schemes) & nzchar(schemes)]
    tibble::tibble(
      rule_id = rule,
      n = as.integer(sum(hit)),
      candidates = as.integer(length(unique(registry$candidate_id[hit]))),
      schemes = if (length(schemes)) paste(schemes, collapse = "; ") else NA_character_
    )
  }) |>
    purrr::list_rbind() |>
    dplyr::arrange(dplyr::desc(.data$n), .data$rule_id)
}
