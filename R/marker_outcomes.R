# =============================================================================
#  marker_outcomes.R  --  assay outcomes fed back to a design (YP5)
# =============================================================================
#  The corpus holds denominators per study (Z8), but no table carries the
#  primer triple and what the assay did in a lab.  This file is the return leg:
#  a call names a stable artefact id, a result and, for a failure, the reason.
#
#  Three rules, and they are the same rules Z8 keeps:
#    * `result` is ordinal ("works", "fails", "unknown"), never a probability;
#    * every row is `calibrated = FALSE` - these are observations, not a
#      correction, and nothing in this package turns them into a success rate;
#    * a row may not be pooled with another study's row, which is why each row
#      carries its own artefact id and nothing aggregates across ids here.
#
#  Retesting appends a row rather than overwriting one: the same assay run
#  twice is two outcomes, and collapsing them would hide the disagreement.
# =============================================================================

#' The columns of the outcome table
#' @keywords internal
#' @noRd
outcome_columns <- function() {
  c("id", "recorded", "target", "result", "notes", "calibrated")
}

#' @keywords internal
#' @noRd
outcome_types <- function() {
  readr::cols(
    id = readr::col_character(),
    recorded = readr::col_character(),
    target = readr::col_character(),
    result = readr::col_character(),
    notes = readr::col_character(),
    calibrated = readr::col_logical()
  )
}

#' @keywords internal
#' @noRd
outcome_empty <- function() {
  tibble::tibble(
    id = character(0), recorded = character(0), target = character(0),
    result = character(0), notes = character(0), calibrated = logical(0)
  )
}

#' The targets an artefact recorded
#'
#' `targets.csv` is written by [marker_artefact()]; `decisions.csv` is the
#' fallback, because a solve whose targets all failed still has attempts.
#' @keywords internal
#' @noRd
artefact_targets <- function(dir) {
  for (name in c("targets.csv", "decisions.csv")) {
    path <- fs::path(dir, name)
    if (!fs::file_exists(path)) next
    tbl <- tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(tbl) && "target" %in% names(tbl)) {
      return(unique(as.character(tbl$target)))
    }
  }
  NULL
}

#' Record what a lab got out of a designed assay
#'
#' Appends one row to `outcomes.csv` beside the artefact's manifest, bound to
#' the stable `id` [marker_artefact()] computes.  A re-test appends another
#' row: the table keeps every observation instead of the last one.
#'
#' `result` is ordinal.  `"works"` and `"fails"` are what the assay did;
#' `"unknown"` says the run did not settle it.  A `"fails"` row **must** carry
#' `notes`, because a failure without a reason cannot be told apart from a
#' failed design, a failed DNA prep or a failed read.  Every row is written
#' with `calibrated = FALSE`: the table is evidence, not a correction, and no
#' function here turns it into a probability of success.
#' @param artefact a `marker_artefact` (see [marker_artefact()]) or the
#'   directory it was written to.
#' @param result one of `"works"`, `"fails"`, `"unknown"`.
#' @param target the target the assay tested, as it appears in the artefact's
#'   `targets.csv`.  `NULL` means the whole artefact.
#' @param notes free text; required when `result = "fails"`.
#' @return the outcome table after the append, invisibly (see
#'   [marker_outcomes()]).
#' @examples
#' \dontrun{
#' art <- marker_export(marker_solve(plan), dir = "results/panel")
#' marker_record_outcome(art, "works", target = "ENSVATH01332998")
#' marker_record_outcome(art, "fails",
#'   target = "ENSVATH02162468",
#'   notes = "no amplification at 60 C"
#' )
#' marker_outcomes(art)
#' }
#' @family KASP plan
#' @export
marker_record_outcome <- function(artefact,
                                  result = c("works", "fails", "unknown"),
                                  target = NULL,
                                  notes = NULL) {
  result <- match.arg(result)
  dir <- artefact_dir(artefact)
  manifest <- artefact_manifest(dir)

  if (!is.null(notes) && !(is.character(notes) && length(notes) == 1L)) {
    cli::cli_abort(gettext("`notes` must be a single string."))
  }
  notes <- trimws(notes %||% NA_character_)
  if (is.na(notes) || !nzchar(notes)) {
    notes <- NA_character_
  }
  if (identical(result, "fails") && is.na(notes)) {
    cli::cli_abort(gettext(
      "A failing outcome must carry `notes`: without a reason it cannot be told \\
       apart from a design failure or a failed run."
    ))
  }

  if (!is.null(target)) {
    if (!is.character(target) || length(target) != 1L || is.na(target)) {
      cli::cli_abort(gettext("`target` must be one label or NULL."))
    }
    known <- artefact_targets(dir)
    if (!is.null(known) && !target %in% known) {
      cli::cli_abort(gettextf(
        "Target %s is not in this artefact. Known targets: %s.",
        shQuote(target), paste(utils::head(known, 10L), collapse = ", ")
      ))
    }
  }

  row <- tibble::tibble(
    id = as.character(manifest$id %||% NA_character_),
    recorded = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    target = if (is.null(target)) NA_character_ else as.character(target),
    result = result,
    notes = notes,
    calibrated = FALSE
  )

  path <- fs::path(dir, "outcomes.csv")
  if (fs::file_exists(path)) {
    header <- names(utils::read.csv(path, nrows = 1L))
    if (!identical(header, outcome_columns())) {
      cli::cli_abort(gettextf(
        "The existing outcomes.csv does not have the expected columns (%s).",
        paste(outcome_columns(), collapse = ", ")
      ))
    }
    readr::write_csv(row, path, append = TRUE, col_names = FALSE)
  } else {
    readr::write_csv(row, path, col_names = TRUE)
  }

  invisible(marker_outcomes(dir))
}

#' The outcomes recorded for an artefact
#'
#' Reads `outcomes.csv` back as a tibble.  An artefact with no recorded
#' outcome returns the empty table rather than `NULL`, so a caller can append
#' to it without a special case.
#'
#' The rows are deliberately left unaggregated: a success rate computed here
#' would pool runs that a lab has not agreed are comparable, which is the
#' mistake [the Z8 tests](../../tests/testthat/test-z8-outcomes.R) exist to
#' prevent.
#' @param artefact a `marker_artefact` or its directory.
#' @return a tibble with `id`, `recorded`, `target`, `result`, `notes`, and
#'   `calibrated`.
#' @examples
#' \dontrun{
#' marker_outcomes("results/panel")
#' }
#' @family KASP plan
#' @export
marker_outcomes <- function(artefact) {
  dir <- artefact_dir(artefact)
  path <- fs::path(dir, "outcomes.csv")
  if (!fs::file_exists(path)) {
    return(outcome_empty())
  }
  readr::read_csv(path, col_types = outcome_types(), na = c("", "NA"))
}
