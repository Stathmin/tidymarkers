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
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' solved <- marker_solve(marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm"))
#' art <- marker_artefact(solved, dir = tempfile("artefact"))
#' marker_record_outcome(art, "works")
#' marker_record_outcome(art, "fails", notes = "no amplification at 60 C")
#' marker_outcomes(art)
#' @family plan
#' @export
marker_record_outcome <- function(artefact,
                                  result = c("works", "fails", "unknown"),
                                  target = NULL,
                                  notes = NULL) {
  call <- sys.call()
  result <- match.arg(result)
  dir <- artefact_dir(artefact)
  manifest <- artefact_manifest(dir)
  if (contracts_enabled()) {
    # a real strengthening: `id` was silently `NA` before, so the outcome could
    # not be bound back to the artefact it claims to describe
    marker_contract(
      is.character(manifest$id) && length(manifest$id) == 1L && !is.na(manifest$id) && nzchar(manifest$id),
      gettext("Contract violated in marker_record_outcome: the artefact id resolves to no existing artefact"),
      call = call
    )
  }

  if (!is.null(notes) && !(is.character(notes) && length(notes) == 1L)) {
    abort_msg("`notes` must be a single string.")
  }
  notes <- trimws(notes %||% NA_character_)
  if (is.na(notes) || !nzchar(notes)) {
    notes <- NA_character_
  }
  if (identical(result, "fails") && is.na(notes)) {
    abort_msg(
      "A failing outcome must carry `notes`: without a reason it cannot be told \\
       apart from a design failure or a failed run."
    )
  }

  if (!is.null(target)) {
    if (!is.character(target) || length(target) != 1L || is.na(target)) {
      abort_msg("`target` must be one label or NULL.")
    }
    known <- artefact_targets(dir)
    if (!is.null(known) && !target %in% known) {
      abort_msg(
        "Target {shQuote(target)} is not in this artefact. Known targets: {paste(utils::head(known, 10L), collapse = ', ')}."
      )
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
      abort_msg(
        "The existing outcomes.csv does not have the expected columns ({paste(outcome_columns(), collapse = ', ')})."
      )
    }
    readr::write_csv(row, path, append = TRUE, col_names = FALSE)
  } else {
    readr::write_csv(row, path, col_names = TRUE)
  }

  out <- marker_outcomes(dir)
  if (contracts_enabled()) {
    # `NA` target means "the whole artefact", so it must match as `NA`, not as a
    # string, or a whole-artefact call would look unreadable
    target_same <- (is.na(out$target) & is.na(row$target)) |
      (!is.na(out$target) & !is.na(row$target) & out$target == row$target)
    hit <- out$id == row$id & out$result == row$result & target_same & !out$calibrated
    marker_contract(
      isTRUE(any(hit)),
      gettext("Contract violated in marker_record_outcome: the outcome is not readable back by marker_outcomes()"),
      call = call
    )
  }
  invisible(out)
}

#' The outcomes recorded for an artefact
#'
#' Reads `outcomes.csv` back as a tibble.  An artefact with no recorded
#' outcome returns the empty table rather than `NULL`, so a caller can append
#' to it without a special case.
#'
#' The rows are deliberately left unaggregated: a success rate computed here
#' would pool runs that a lab has not agreed are comparable, which is the
#' mistake [the Z8 tests](https://github.com/Stathmin/tidymarkers/blob/main/tests/testthat/test-z8-outcomes.R)
#' exist to prevent.
#' @param artefact a `marker_artefact` or its directory.
#' @return a tibble with `id`, `recorded`, `target`, `result`, `notes`, and
#'   `calibrated`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' solved <- marker_solve(marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm"))
#' art <- marker_artefact(solved, dir = tempfile("artefact"))
#' marker_record_outcome(art, "works")
#' marker_outcomes(art)
#' @family plan
#' @export
marker_outcomes <- function(artefact) {
  dir <- artefact_dir(artefact)
  path <- fs::path(dir, "outcomes.csv")
  if (!fs::file_exists(path)) {
    return(outcome_empty())
  }
  readr::read_csv(path, col_types = outcome_types(), na = c("", "NA"))
}

#' The outcomes recorded across every artefact under a directory
#'
#' [marker_outcomes()] reads one artefact.  A study grows by batches, so the
#' question "what did the lab actually get" spans directories.  This walker
#' binds every `outcomes.csv` it finds under `dir` and adds the two columns the
#' single-artefact reader cannot know: `artefact_id`, the id of the containing
#' artefact's manifest, and `artefact`, its directory relative to `dir`.
#'
#' A row's own `id` is left exactly as it was recorded, even when it disagrees
#' with the manifest of the directory holding it: relabelling a laboratory's
#' observation to match a rewritten artefact is worse than showing both, so a
#' disagreement stays visible as two different columns.
#'
#' Nothing is aggregated here, for the reason [marker_outcomes()] gives: rows
#' from different designs are not comparable until a laboratory says they are.
#' The corpus is the raw material of that decision, not the decision.
#' @param dir a directory searched recursively for `outcomes.csv`.
#' @return a tibble with `artefact_id`, `artefact`, and the outcome columns of
#'   [marker_outcomes()].  A directory with no recorded outcome returns the
#'   empty table.  The absolute scanned root travels with the table as its
#'   `root` attribute, because the relative `artefact` column alone cannot
#'   locate a directory again; [marker_export_outcomes()] needs it.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' solved <- marker_solve(marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm"))
#' art <- marker_artefact(solved, dir = file.path(tempdir(), "run"))
#' marker_record_outcome(art, "works")
#' marker_outcome_corpus(tempdir())
#' @family plan
#' @export
marker_outcome_corpus <- function(dir) {
  if (!fs::dir_exists(dir)) {
    abort_msg("`dir` must be a directory searched for outcomes.csv.")
  }
  root <- fs::path_abs(dir)
  files <- fs::dir_ls(root, recurse = TRUE, type = "file")
  files <- files[fs::path_file(files) == "outcomes.csv"]

  parts <- lapply(files, function(path) {
    tbl <- readr::read_csv(
      path,
      col_types = outcome_types(), na = c("", "NA"), show_col_types = FALSE
    )
    adir <- fs::path_dir(path)
    manifest <- tryCatch(artefact_manifest(adir), error = function(e) NULL)
    mid <- if (is.null(manifest)) NA_character_ else manifest$id
    id <- if (length(mid) == 1L && !is.na(mid) && nzchar(mid)) {
      mid
    } else {
      fs::path_file(adir)
    }
    tibble::add_column(
      tbl,
      artefact_id = as.character(id),
      artefact = as.character(fs::path_rel(adir, root)),
      .before = 1
    )
  })

  if (length(parts) == 0L) {
    return(structure(
      tibble::add_column(
        outcome_empty(),
        artefact_id = character(0), artefact = character(0),
        .before = 1
      ),
      root = as.character(root)
    ))
  }
  # bind_rows drops attributes, so the root is attached after the bind.
  structure(dplyr::bind_rows(parts), root = as.character(root))
}
