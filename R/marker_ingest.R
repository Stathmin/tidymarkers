# =============================================================================
#  marker_ingest.R  --  literature corpus -> review queue -> approved registry
# =============================================================================
#  The chemistry performance registry grows from a corpus a maintainer drops in
#  `marker_cache_dir()`, never from a file shipped with the package (external
#  data does not travel with this package).  Ingestion is deliberately passive:
#  it reads that file, sorts every row into `candidates` or `rejected`, and
#  stops.  Nothing reaches the registry until a maintainer approves one record
#  by id, so a corpus row - or a row this module mis-parsed - can never edit a
#  stored performance figure on its own.
#
#  A rejected row is kept, with its `reason`, in input order.  A silently
#  dropped row is a decision nobody can audit; a rejected row is one a
#  maintainer can fix upstream and re-ingest.
# =============================================================================

#' Where the approved literature records are kept
#' @return a path under [marker_cache_dir()].  The file appears on the first
#'   [marker_ingest_approve()] and is never created by a read.
#' @keywords internal
#' @noRd
ingest_registry_path <- function() {
  fs::path(marker_cache_dir(), "chemistry_performance.csv")
}

#' The columns of one candidate row
#' @keywords internal
#' @noRd
ingest_columns <- function() {
  c("record_id", "doi", "metric", "value", "evidence", "confidence", "notes")
}

#' The columns of an approved registry row
#' @keywords internal
#' @noRd
ingest_registry_columns <- function() {
  c(ingest_columns(), "approved_by", "approved_at")
}

#' A zero-row registry, with the registry columns
#' @keywords internal
#' @noRd
ingest_registry_empty <- function() {
  tibble::tibble(
    record_id = character(0), doi = character(0), metric = character(0),
    value = numeric(0), evidence = integer(0), confidence = numeric(0),
    notes = character(0), approved_by = character(0), approved_at = character(0)
  )
}

#' Build one of the module's tibbles from column vectors
#'
#' Shared by `candidates`, `rejected` and a registry row, so the three can never
#' disagree about a column's type.
#' @keywords internal
#' @noRd
ingest_rows <- function(record_id, doi, metric, value, evidence, confidence, notes) {
  tibble::tibble(
    record_id = as.character(record_id),
    doi = as.character(doi),
    metric = as.character(metric),
    value = as.numeric(value),
    evidence = as.integer(evidence),
    confidence = as.numeric(confidence),
    notes = as.character(notes)
  )
}

#' Validate the three ingest thresholds
#' @keywords internal
#' @noRd
ingest_check_criteria <- function(min_evidence, min_confidence, provenance_required) {
  if (
    !is.numeric(min_evidence) || length(min_evidence) != 1L || is.na(min_evidence) ||
      min_evidence < 0 || min_evidence != round(min_evidence)
  ) {
    cli::cli_abort(gettext("`min_evidence` must be one non-negative whole number."))
  }
  if (
    !is.numeric(min_confidence) || length(min_confidence) != 1L || is.na(min_confidence) ||
      min_confidence < 0 || min_confidence > 1
  ) {
    cli::cli_abort(gettext("`min_confidence` must be one number in [0, 1]."))
  }
  if (!is.logical(provenance_required) || length(provenance_required) != 1L || is.na(provenance_required)) {
    cli::cli_abort(gettext("`provenance_required` must be TRUE or FALSE."))
  }
  invisible(NULL)
}

#' Normalise and check the corpus path
#' @keywords internal
#' @noRd
ingest_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    cli::cli_abort(gettext("`path` must be one existing .csv or .json file."))
  }
  path <- fs::path_abs(path)
  if (!fs::is_file(path)) {
    cli::cli_abort(gettextf(
      "`path` must be an existing .csv or .json file; %s is not one.",
      encodeString(as.character(path), quote = "`")
    ))
  }
  ext <- tolower(fs::path_ext(path))
  if (!ext %in% c("csv", "json")) {
    cli::cli_abort(gettextf(
      "`path` must end in .csv or .json, not .%s (%s).", ext,
      encodeString(as.character(path), quote = "`")
    ))
  }
  fs::path_norm(path)
}

#' Read a corpus as one character column per field
#'
#' Everything is read as text and parsed by this module: type guessing would
#' turn a column whose first row is blank into a logical column of `NA`s, and a
#' blank `value` in row 1 would then be indistinguishable from a bad one.
#' @keywords internal
#' @noRd
ingest_read <- function(path) {
  if (tolower(fs::path_ext(path)) == "json") {
    return(ingest_read_json(path))
  }
  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA"), progress = FALSE
  )
}

#' Read a JSON corpus: a top-level array of objects with the same fields
#' @keywords internal
#' @noRd
ingest_read_json <- function(path) {
  items <- jsonlite::fromJSON(readr::read_file(path), simplifyVector = FALSE)
  if (!is.list(items) || !is.null(names(items)) || !all(vapply(items, is.list, logical(1)))) {
    cli::cli_abort(gettext("A JSON literature source must be a top-level array of objects."))
  }
  if (!length(items)) {
    return(tibble::tibble(
      record_id = character(0), doi = character(0), metric = character(0),
      value = character(0), evidence = character(0), confidence = character(0),
      notes = character(0)
    ))
  }
  fields <- unique(unlist(lapply(items, names), use.names = FALSE))
  cols <- lapply(stats::setNames(fields, fields), function(f) {
    vapply(items, function(item) {
      v <- item[[f]]
      if (is.null(v) || length(v) != 1L) NA_character_ else as.character(v)
    }, character(1))
  })
  tibble::as_tibble(cols)
}

#' Check the required columns and fill the optional ones with `NA`
#' @keywords internal
#' @noRd
ingest_fields <- function(dat) {
  required <- c("record_id", "doi", "metric", "value")
  absent <- setdiff(required, names(dat))
  if (length(absent)) {
    cli::cli_abort(gettextf(
      "The literature source is missing required column(s): %s.",
      paste(absent, collapse = ", ")
    ))
  }
  n <- nrow(dat)
  fill <- function(name) {
    if (name %in% names(dat)) as.character(dat[[name]]) else rep(NA_character_, n)
  }
  tibble::tibble(
    record_id = fill("record_id"), doi = fill("doi"), metric = fill("metric"),
    value = fill("value"), evidence = fill("evidence"),
    confidence = fill("confidence"), notes = fill("notes")
  )
}

#' Sort the source rows into candidates and rejected rows
#'
#' The precedence is fixed, so the same file always yields the same queue: a row
#' is judged by the first rule it breaks.  A `record_id` is claimed by its first
#' occurrence in file order, so a later copy is a duplicate even when the first
#' copy is itself rejected; the queue never has to decide which of two rows "is"
#' the record.  Every input row comes out in exactly one of the two tibbles.
#' @keywords internal
#' @noRd
ingest_classify <- function(dat, min_evidence, min_confidence, provenance_required) {
  n <- nrow(dat)
  id <- dat$record_id
  doi <- dat$doi
  metric <- dat$metric
  raw_evidence <- dat$evidence
  raw_confidence <- dat$confidence
  value <- suppressWarnings(as.numeric(dat$value))
  evidence <- suppressWarnings(as.numeric(raw_evidence))
  confidence <- suppressWarnings(as.numeric(raw_confidence))
  reason <- rep(NA_character_, n)
  claimed <- character(n)
  n_claimed <- 0L

  for (i in seq_len(n)) {
    rid <- id[i]
    if (is.na(rid) || !nzchar(rid)) {
      reason[i] <- "missing record_id"
      next
    }
    if (rid %in% claimed[seq_len(n_claimed)]) {
      reason[i] <- paste0("duplicate record_id ", rid)
      next
    }
    n_claimed <- n_claimed + 1L
    claimed[n_claimed] <- rid

    if (is.na(doi[i]) || !nzchar(doi[i])) {
      if (provenance_required) {
        reason[i] <- "no provenance: a doi is required"
        next
      }
    } else if (!grepl("^10\\.", doi[i])) {
      reason[i] <- paste0("doi ", doi[i], " does not look like a DOI (no 10. prefix)")
      next
    }
    if (is.na(metric[i]) || !nzchar(metric[i])) {
      reason[i] <- "missing metric"
      next
    }
    if (is.na(value[i]) || !is.finite(value[i])) {
      reason[i] <- paste0("malformed value ", dat$value[i])
      next
    }
    if (!is.na(raw_evidence[i])) {
      if (
        is.na(evidence[i]) || !is.finite(evidence[i]) || evidence[i] != round(evidence[i]) ||
          evidence[i] < 0 || evidence[i] > .Machine$integer.max
      ) {
        reason[i] <- paste0("malformed evidence ", raw_evidence[i], " (a non-negative whole number)")
        next
      }
      if (evidence[i] < min_evidence) {
        reason[i] <- paste0("evidence ", evidence[i], " is below min_evidence ", min_evidence)
        next
      }
    }
    if (!is.na(raw_confidence[i])) {
      if (
        is.na(confidence[i]) || !is.finite(confidence[i]) ||
          confidence[i] < 0 || confidence[i] > 1
      ) {
        reason[i] <- paste0("malformed confidence ", raw_confidence[i], " (must be in [0, 1])")
        next
      }
      if (confidence[i] < min_confidence) {
        reason[i] <- paste0("confidence ", confidence[i], " is below min_confidence ", min_confidence)
        next
      }
    }
  }

  keep <- is.na(reason)
  # a row that omits `evidence` takes its documented default of 1; a rejected
  # row keeps the `NA` its source produced, so the queue never invents a field
  evidence[keep & is.na(evidence)] <- 1

  rejected <- ingest_rows(
    record_id = id[!keep], doi = doi[!keep], metric = metric[!keep],
    value = value[!keep], evidence = evidence[!keep], confidence = confidence[!keep],
    notes = dat$notes[!keep]
  )
  rejected$reason <- reason[!keep]

  list(
    candidates = ingest_rows(
      record_id = id[keep], doi = doi[keep], metric = metric[keep],
      value = value[keep], evidence = evidence[keep], confidence = confidence[keep],
      notes = dat$notes[keep]
    ),
    rejected = rejected
  )
}

#' The timestamp an approval is stamped with
#' @keywords internal
#' @noRd
ingest_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

#' Add the approval stamp to candidate rows
#' @keywords internal
#' @noRd
ingest_stamp <- function(rows, reviewer) {
  tibble::as_tibble(c(
    as.list(rows),
    list(
      approved_by = rep(as.character(reviewer), nrow(rows)),
      approved_at = rep(ingest_now(), nrow(rows))
    )
  ))
}

#' Readable types for the registry
#' @keywords internal
#' @noRd
ingest_registry_types <- function() {
  readr::cols(
    record_id = readr::col_character(),
    doi = readr::col_character(),
    metric = readr::col_character(),
    value = readr::col_double(),
    evidence = readr::col_integer(),
    confidence = readr::col_double(),
    notes = readr::col_character(),
    approved_by = readr::col_character(),
    approved_at = readr::col_character()
  )
}

#' Append approved rows to the registry file
#' @keywords internal
#' @noRd
ingest_write_registry <- function(rows, registry) {
  registry <- fs::path(registry)
  dir <- fs::path_dir(registry)
  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir, recurse = TRUE)
  }
  if (fs::file_exists(registry)) {
    header <- names(utils::read.csv(registry, nrows = 1L))
    if (!identical(header, names(rows))) {
      cli::cli_abort(gettextf(
        "The existing registry %s does not have the expected columns (%s).",
        encodeString(as.character(registry), quote = "`"),
        paste(ingest_registry_columns(), collapse = ", ")
      ))
    }
    readr::write_csv(rows, registry, append = TRUE, col_names = FALSE)
  } else {
    readr::write_csv(rows, registry, col_names = TRUE)
  }
  invisible(registry)
}

#' Ingest an external literature corpus into a review queue
#'
#' Reads a structured corpus (CSV or JSON) of chemistry performance records and
#' sorts every row into `candidates` or `rejected`.  Nothing is merged: a
#' candidate reaches the registry only through [marker_ingest_approve()], so a
#' downloaded file can never change a stored figure by itself.
#'
#' A row is rejected when it has no `record_id`, repeats an earlier row's
#' `record_id`, has no `doi` while `provenance_required` is `TRUE`, has a `doi`
#' that does not start with `10.`, has no `metric`, has a non-numeric `value`,
#' has a malformed `evidence` or `confidence`, or falls below `min_evidence` or
#' `min_confidence`.  The first rule a row breaks is the reason recorded, and a
#' rejected row keeps its parsed fields next to that reason.
#'
#' @param path path to an existing `.csv` (readr) or `.json` (jsonlite, a
#'   top-level array of objects with the same fields) file.  A corpus downloaded
#'   into [marker_cache_dir()] is the intended source; this function downloads
#'   nothing and writes nothing.
#' @param min_evidence smallest accepted `evidence`, a non-negative whole
#'   number.  A row that omits `evidence` takes `1`.
#' @param min_confidence smallest accepted `confidence`, a number in `[0, 1]`.
#'   A row that omits `confidence` is not filtered on it.
#' @param provenance_required when `TRUE` (the default) a row with a blank `doi`
#'   is rejected: no DOI, no provenance chain, no queue.
#' @return a `marker_ingest_queue`: a list with the `candidates` and `rejected`
#'   tibbles (`rejected` carries `reason`), the normalised `source` path, the
#'   `criteria` actually applied, and an empty `approved` tibble.
#' @examples
#' corpus <- fs::path(tempdir(), "ingest-corpus.csv")
#' readr::write_lines(
#'   c(
#'     "record_id,doi,metric,value,evidence,confidence",
#'     "r1,10.1000/example,LoD,0.9,2,0.8",
#'     "r2,,LoD,0.9,2,0.8"
#'   ),
#'   corpus
#' )
#' queue <- marker_ingest_literature(corpus)
#' queue$candidates
#' queue$rejected
#' @family input
#' @export
marker_ingest_literature <- function(path, min_evidence = 1L, min_confidence = 0,
                                     provenance_required = TRUE) {
  ingest_check_criteria(min_evidence, min_confidence, provenance_required)
  path <- ingest_path(path)
  parts <- ingest_classify(ingest_fields(ingest_read(path)), min_evidence, min_confidence, provenance_required)
  structure(
    list(
      candidates = parts$candidates,
      rejected = parts$rejected,
      source = as.character(path),
      criteria = list(
        min_evidence = as.integer(min_evidence),
        min_confidence = as.numeric(min_confidence),
        provenance_required = provenance_required
      ),
      approved = ingest_registry_empty()
    ),
    class = "marker_ingest_queue"
  )
}

#' Approve one queued record into the registry
#'
#' The only write path in this module, and the only way a corpus row reaches the
#' registry.  The row is written with `approved_by` and `approved_at` (UTC,
#' ISO-8601), and the queue it returns has it moved from `candidates` to
#' `approved`.  Approving an id that is already in the registry is an error
#' rather than an overwrite: two conflicting values for one record are a
#' maintainer's decision, not something an append should settle silently.
#'
#' @param queue a `marker_ingest_queue` from [marker_ingest_literature()].
#' @param record_id the `record_id` of the candidate to approve.
#' @param registry path to the registry CSV; created (with its parent
#'   directory) on the first approval, appended to afterwards.
#' @param reviewer the name stamped into `approved_by`.
#' @return the updated queue, invisibly.
#' @examples
#' corpus <- fs::path(tempdir(), "ingest-approve.csv")
#' readr::write_lines(c("record_id,doi,metric,value", "r1,10.1000/example,LoD,0.9"), corpus)
#' registry <- fs::path(tempdir(), "ingest-registry.csv")
#' unlink(registry)
#' queue <- marker_ingest_literature(corpus)
#' queue <- marker_ingest_approve(queue, "r1", registry = registry)
#' marker_ingest_registry(registry)
#' @family input
#' @export
marker_ingest_approve <- function(queue, record_id, registry = ingest_registry_path(),
                                  reviewer = "maintainer") {
  if (!inherits(queue, "marker_ingest_queue")) {
    cli::cli_abort(gettext("`queue` must come from marker_ingest_literature()."))
  }
  if (!is.character(record_id) || length(record_id) != 1L || is.na(record_id) || !nzchar(record_id)) {
    cli::cli_abort(gettext("`record_id` must be one non-empty string."))
  }
  if (!is.character(reviewer) || length(reviewer) != 1L || is.na(reviewer) || !nzchar(reviewer)) {
    cli::cli_abort(gettext("`reviewer` must be one non-empty string."))
  }
  rows <- queue$candidates[queue$candidates$record_id == record_id, , drop = FALSE]
  if (!nrow(rows)) {
    cli::cli_abort(gettextf(
      "Record %s is not among this queue's candidates.", encodeString(record_id, quote = "`")
    ))
  }
  if (record_id %in% marker_ingest_registry(registry)$record_id) {
    cli::cli_abort(gettextf(
      "Record %s is already in the registry %s.", encodeString(record_id, quote = "`"),
      encodeString(as.character(fs::path(registry)), quote = "`")
    ))
  }
  stamped <- ingest_stamp(rows, reviewer)
  ingest_write_registry(stamped, registry)
  queue$approved <- dplyr::bind_rows(queue$approved, stamped)
  queue$candidates <- queue$candidates[queue$candidates$record_id != record_id, , drop = FALSE]
  invisible(queue)
}

#' The approved literature records
#'
#' @param registry path to the registry CSV (see [marker_ingest_approve()]).
#' @return a tibble with [marker_ingest_approve()]'s columns; zero rows when
#'   nothing has been approved yet.  A missing file is not created here.
#' @examples
#' marker_ingest_registry(fs::path(tempdir(), "ingest-registry-absent.csv"))
#' @family input
#' @export
marker_ingest_registry <- function(registry = ingest_registry_path()) {
  registry <- fs::path(registry)
  if (!fs::file_exists(registry)) {
    return(ingest_registry_empty())
  }
  readr::read_csv(registry, col_types = ingest_registry_types(), na = c("", "NA"), progress = FALSE)
}
