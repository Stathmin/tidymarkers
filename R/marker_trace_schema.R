# =============================================================================
#  marker_trace_schema.R  --  the versioned record a pooled corpus is made of (ZY5)
# =============================================================================
#  An artefact is the record a single reader can replay.  A trace is the
#  flattened, publication-shaped projection of it: the same provenance, plan,
#  counts and decisions, with no replay copy and no absolute path.  That is what
#  a community collection would pool, so it needs a version of its own and a
#  migration path that never drops a field.
#
#  Two rules this file keeps deliberately:
#    * the schema version is a constant, not an option: a reader either knows the
#      version it was built for or says so and stops - an option would let two
#      readers of the same file disagree about what they are reading;
#    * a migration maps field by field and is idempotent.  A `modifyList()` guess
#      would let a renamed or vanished field pass silently, which is exactly the
#      loss this file exists to prevent.
#
#  ZY7 is the separate anonymisation promise; a trace only refuses to carry an
#  identifier in the first place.  Publishing a community collection (Zenodo,
#  ExperimentHub) is the blocked remainder of ZY5 and is not implemented here.
# =============================================================================

TRACE_SCHEMA <- "tidymarkers/trace@1"

# the newest version this reader understands, per record kind
TRACE_FORMATS <- c(trace = 1L, artefact = 1L)

# the shape of every trace, in order: an artefact manifest's fields plus the two
# things only a trace carries (`decisions`, `signals`), so a trace has one shape
# whether it came from an artefact or a solved object
TRACE_FIELDS <- c(
  "schema", "created", "package", "platform", "engine_stamp", "inputs", "plan",
  "parameters", "overrides", "limitations", "counts", "files", "id",
  "decisions", "signals"
)

# --- schema strings ----------------------------------------------------------

#' Split a schema identifier into its kind and version
#' @param schema a candidate schema string.
#' @return `list(kind, version)`, or `NULL` when `schema` is not one this
#'   package understands.
#' @keywords internal
#' @noRd
trace_schema_parts <- function(schema) {
  if (!is.character(schema) || length(schema) != 1L || is.na(schema)) {
    return(NULL)
  }
  m <- regmatches(schema, regexec("^tidymarkers/(trace|artefact)@([0-9]+)$", schema))[[1L]]
  if (!length(m)) {
    return(NULL)
  }
  list(kind = m[[2L]], version = as.numeric(m[[3L]]))
}

#' A schema string for a message, never `NULL`
#' @keywords internal
#' @noRd
trace_schema_label <- function(schema) {
  if (is.character(schema) && length(schema) == 1L && !is.na(schema)) {
    schema
  } else {
    "<none>"
  }
}

trace_abort_unknown <- function(schema) {
  cli::cli_abort(c(
    gettext("Not a trace: `schema` is {.val {trace_schema_label(schema)}}."),
    i = gettext(
      "Expected {.val {marker_trace_schema()}} or a legacy {.val {'tidymarkers/artefact@1'}} manifest."
    )
  ))
}

# --- the aggregation contract ------------------------------------------------

# A trace carries no absolute path, wherever the path sits: the input rows are
# the obvious place, but a knob or a designer's reason can quote one too.  The
# value becomes its file name, which still names the reference without naming
# the machine.  This is not anonymisation (ZY7); it only refuses the identifier.
trace_basename_paths <- function(x) {
  if (!length(x)) {
    return(x)
  }
  if (is.character(x)) {
    hit <- !is.na(x) & grepl("^(/|[A-Za-z]:[/\\\\])", x)
    x[hit] <- fs::path_file(x[hit])
    return(x)
  }
  if (is.data.frame(x)) {
    for (nm in names(x)) {
      x[[nm]] <- trace_basename_paths(x[[nm]])
    }
    return(x)
  }
  if (is.list(x)) {
    x[] <- lapply(x, trace_basename_paths)
  }
  x
}

# The decision table of a record that has none, with the columns a reader of a
# trace expects.  `artefact@1` kept these rows in `decisions.csv` beside the
# manifest, so a bare manifest cannot be made to carry them back.
trace_empty_decisions <- function() {
  artefact_decisions(list(attempts = tibble::tibble()))
}

# --- record construction -----------------------------------------------------

#' Flatten an artefact manifest into a trace
#'
#' Every field is mapped by name; nothing is copied wholesale and nothing is
#' guessed.  `trace_schema` is deliberately not carried: on a trace the schema is
#' the trace's own version, so the pointer would only disagree with `schema`.
#' @param manifest an `artefact@1` manifest (or its successor).
#' @param decisions the artefact's decision tibble, when the caller has it.
#' @param signals the artefact's signal table, or `NULL`.
#' @return a `marker_trace`.
#' @keywords internal
#' @noRd
trace_from_manifest <- function(manifest, decisions = NULL, signals = NULL) {
  out <- list(
    schema = marker_trace_schema(),
    created = manifest$created,
    package = manifest$package,
    platform = manifest$platform,
    engine_stamp = manifest$engine_stamp,
    inputs = manifest$inputs,
    plan = manifest$plan,
    parameters = manifest$parameters,
    overrides = manifest$overrides,
    limitations = manifest$limitations,
    counts = manifest$counts,
    files = manifest$files,
    id = manifest$id,
    decisions = if (is.data.frame(decisions)) decisions else trace_empty_decisions(),
    signals = signals
  )
  structure(trace_basename_paths(out), class = c("marker_trace", "list"))
}

#' The trace a solved plan would produce, without writing an artefact
#'
#' Mirrors the manifest `marker_artefact()` builds, minus its sidecars: no
#' `files` are written, so the field is `NULL`.  The id excludes `created`, so
#' two traces of the same plan agree.
#' @param solved a `marker_solved`.
#' @return a `marker_trace`.
#' @keywords internal
#' @noRd
trace_from_solved <- function(solved) {
  knobs <- solved$knobs %||% list()
  decisions <- artefact_decisions(solved)
  signals <- solved$signals %||% NULL
  tg <- solved$plan$targets
  out <- list(
    schema = marker_trace_schema(),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    package = list(
      name = "tidymarkers",
      version = as.character(utils::packageVersion("tidymarkers"))
    ),
    platform = list(
      r = as.character(getRversion()),
      os = R.version$platform
    ),
    engine_stamp = as.list(marker_engine_stamp()),
    inputs = artefact_inputs(solved$plan),
    plan = list(
      layers = solved$plan$layers,
      targets = nrow(tg),
      schemes_available = solved$schemes,
      costs = as.list(solved$plan$costs %||% NA)
    ),
    parameters = purrr::map(knobs, artefact_value),
    overrides = artefact_overrides(knobs, unique(solved$attempts$scheme)),
    limitations = I(as.character(artefact_limitations(solved))),
    counts = list(
      targets = nrow(solved$report),
      solved = sum(solved$report$ok),
      attempts = nrow(decisions),
      rejected = sum(decisions$outcome == "rejected"),
      signals = if (is.null(signals)) 0L else nrow(signals)
    ),
    files = NULL
  )
  out$id <- artefact_id(out)
  out$decisions <- decisions
  # `[<-` keeps a NULL signals in its slot; `$<- NULL` would drop the field
  out["signals"] <- list(signals)
  structure(trace_basename_paths(out), class = c("marker_trace", "list"))
}

#' Assert a record is a trace at the current schema
#' @param x a candidate trace.
#' @return `x`, invisibly.
#' @keywords internal
#' @noRd
trace_validate <- function(x) {
  if (!inherits(x, "marker_trace")) {
    cli::cli_abort(gettext("`x` is not a marker_trace; migrate it with marker_trace_migrate()."))
  }
  missing <- setdiff(TRACE_FIELDS, names(x))
  if (length(missing)) {
    cli::cli_abort(c(
      gettext("A trace must have every field."),
      x = gettext("Missing: {paste(missing, collapse = ', ')}.")
    ))
  }
  if (!identical(x$schema, marker_trace_schema())) {
    cli::cli_abort(c(
      gettext("A trace must carry the current schema."),
      x = gettext("Found {.val {trace_schema_label(x$schema)}}, expected {.val {marker_trace_schema()}}.")
    ))
  }
  if (!is.data.frame(x$decisions)) {
    cli::cli_abort(gettext("`decisions` in a trace must be a data frame."))
  }
  invisible(x)
}

# --- the public API ----------------------------------------------------------

#' The current trace schema identifier
#'
#' A trace is the flattened, publication-shaped projection of a design artefact:
#' what a pooled corpus of designs would be made of.  This function names the
#' version of that record, so a reader can say which records it understands
#' before it reads one.
#'
#' The version is a constant, not an option.  A reader either knows the version
#' it was built for or says so and stops; an option would let two readers of the
#' same file disagree about what they are reading.  Publishing a community
#' collection (Zenodo, ExperimentHub) is the blocked remainder of ZY5 and is not
#' implemented here - no archive is claimed to exist.
#' @return a length-one character vector, `"tidymarkers/trace@1"`.
#' @examples
#' marker_trace_schema()
#' @family plan
#' @export
marker_trace_schema <- function() {
  TRACE_SCHEMA
}

#' Build a trace from an artefact or a solved plan
#'
#' A trace carries the provenance, plan, counts, overrides, limitations and
#' decision table of a design, and nothing that identifies a person, a machine
#' hostname or an absolute path.  That is the aggregation contract, not
#' anonymisation: `ZY7` is the anonymisation promise.  Absolute paths are
#' reduced to their file name, and only the fields listed below travel.
#'
#' The fields are `schema`, `created`, `package` (name, version), `platform`
#' (r, os), `engine_stamp`, `inputs`, `plan`, `parameters`, `overrides`,
#' `limitations`, `counts`, `files`, `id`, `decisions` and `signals`.
#'
#' An artefact directory (or its `artefact.json`) is read with
#' [marker_artefact_read()], so the `decisions` and `signals` rows the sidecars
#' carry travel with the trace; a bare manifest in memory arrived without them.
#' @param x a `marker_artefact` (see [marker_artefact()]), the directory it was
#'   written to (or the `artefact.json` inside it), or a `marker_solved` (see
#'   [marker_solve()]).
#' @return a `marker_trace`: a plain list of class `c("marker_trace", "list")`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' plan <- marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm")
#' solved <- marker_solve(plan)
#' trace <- marker_trace_from_artefact(marker_artefact(solved, dir = tempfile("artefact")))
#' names(trace)
#' marker_trace_from_artefact(solved)$id
#' @family plan
#' @export
marker_trace_from_artefact <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    x <- marker_artefact_read(x)
  }
  if (inherits(x, "marker_artefact")) {
    files <- x$manifest$files %||% list()
    signals <- artefact_sidecar(x$dir, files$signals, function() NULL)
    return(trace_from_manifest(x$manifest, decisions = x$decisions, signals = signals))
  }
  if (inherits(x, "marker_solved")) {
    return(trace_from_solved(x))
  }
  cli::cli_abort(gettext(
    "`x` must be a marker_artefact, an artefact directory or a marker_solved."
  ))
}

#' Migrate a record to the current trace schema
#'
#' The migration path, in order:
#'
#' 1. a `trace@1` record is returned unchanged (or gains the class when it was
#'    stored without it), so the migration is idempotent;
#' 2. a legacy `artefact@1` manifest is mapped field by field: its schema becomes
#'    the trace schema and the fields a manifest cannot carry (`decisions`,
#'    `signals`) are added empty rather than invented - `artefact@1` kept its
#'    decision rows in `decisions.csv` beside the manifest;
#' 3. an absent or unknown `schema` aborts; a migration never guesses;
#' 4. a version newer than this reader understands aborts and names the version:
#'    the reader is too old, not the record wrong.
#'
#' `trace_schema`, when a record has it, is dropped: on a trace the schema field
#' already names the trace version, and keeping the old pointer would let the two
#' disagree.
#' @param x a record to migrate.
#' @return a `marker_trace` at [marker_trace_schema()].
#' @examples
#' legacy <- list(schema = "tidymarkers/artefact@1")
#' marker_trace_migrate(legacy)$schema
#' @family plan
#' @export
marker_trace_migrate <- function(x) {
  if (!is.list(x)) {
    trace_abort_unknown(NULL)
  }
  parts <- trace_schema_parts(x$schema)
  if (is.null(parts)) {
    trace_abort_unknown(x$schema)
  }
  current <- TRACE_FORMATS[[parts$kind]]
  if (parts$version > current) {
    cli::cli_abort(c(
      gettext("This record is {.val {x$schema}}: newer than this reader's {parts$kind} schema."),
      i = gettext("Update tidymarkers to read it.")
    ))
  }
  if (parts$version < current) {
    # no such version was ever written; treating it as current would guess
    trace_abort_unknown(x$schema)
  }
  if (parts$kind == "trace") {
    if (inherits(x, "marker_trace")) {
      return(x)
    }
    return(structure(x, class = c("marker_trace", "list")))
  }
  trace_from_manifest(x)
}

#' Read a trace from an RDS file
#'
#' Migration happens on read, so a file written by an older tidymarkers comes
#' back at [marker_trace_schema()] and is validated before it is returned.
#' @param path path to an RDS file written by [marker_trace_write()] or an older
#'   reader.
#' @return a `marker_trace`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' solved <- marker_solve(marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm"))
#' path <- tempfile(fileext = ".rds")
#' marker_trace_write(marker_trace_from_artefact(solved), path)
#' marker_trace_read(path)$id
#' @family plan
#' @export
marker_trace_read <- function(path) {
  if (!fs::file_exists(path)) {
    cli::cli_abort(gettext("No trace at {.path {path}}."))
  }
  out <- tryCatch(
    trace_validate(marker_trace_migrate(readRDS(path))),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    cli::cli_abort(c(
      gettext("{.path {path}} is not a readable trace."),
      x = gettext("The error was: {conditionMessage(out)}")
    ))
  }
  out
}

#' Write a trace, always at the current schema
#'
#' Migration and validation happen before the write, so a file this function
#' returns can be read by any reader of [marker_trace_schema()].
#' @param trace a `marker_trace` or a record [marker_trace_migrate()] can read.
#' @param path output file path; its directory must exist.
#' @return `path`, invisibly.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' solved <- marker_solve(marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm"))
#' path <- tempfile(fileext = ".rds")
#' marker_trace_write(marker_trace_from_artefact(solved), path)
#' file.exists(path)
#' @family plan
#' @export
marker_trace_write <- function(trace, path) {
  trace <- trace_validate(marker_trace_migrate(trace))
  dir <- fs::path_dir(path)
  if (!fs::dir_exists(dir)) {
    cli::cli_abort(gettext("No directory at {.path {dir}}; create it before writing a trace."))
  }
  saveRDS(trace, path)
  invisible(path)
}

# =============================================================================
#  Cross-language conformance (ZY22)
# =============================================================================
#  A trace is only worth pooling if a second implementation can write the same
#  record.  This section is the contract between the R writer and
#  `inst/python/marker_trace.py`: one schema file names the four conformance
#  files, their columns and their types, and both writers read it rather than
#  keeping a private copy of the table layout.
#
#  The payload is the transfer format, not a fifth artefact file.  It carries
#  the in-memory manifest (where an `AsIs` length-1 vector and a scalar still
#  differ, which JSON alone would lose) and the sidecar rows already pinned to
#  the schema types.  A second writer re-renders those; it does not re-derive
#  them, so a schema drift breaks the read instead of passing silently.

#' The cross-language conformance schema
#'
#' The schema names the four files a record must reproduce byte for byte
#' (`artefact.json`, `decisions.csv`, `checks.csv`, `why.csv`), each table's
#' columns in order, their types, and the CSV spellings both writers must use.
#' `targets.csv` and `signals.csv` are informational and outside the set.
#' @param path a schema file; the package's own when `NULL`.
#' @return the parsed schema as a nested list (`simplifyVector = FALSE`), so
#'   `columns` stays an ordered list of `list(name, type)`.
#' @keywords internal
#' @noRd
trace_conformance_schema <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", "trace-schema.json", package = "tidymarkers")
  }
  if (!nzchar(path) || !fs::file_exists(path)) {
    cli::cli_abort(gettext("No conformance schema at {.path {path}}."))
  }
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

# readr collectors for one schema table, keyed by the schema's own type names.
# `status` is a factor in memory and a string on disk, so it is pinned to `c`.
trace_pinned_cols <- function(table) {
  map <- c(string = "c", integer = "i", number = "d", logical = "l")
  types <- vapply(table$columns, function(col) col$type, character(1))
  unknown <- unique(types[!types %in% names(map)])
  if (length(unknown)) {
    cli::cli_abort(gettext(
      "Unknown column type {.val {unknown}} in the conformance schema for {.file {table$file}}."
    ))
  }
  do.call(
    readr::cols,
    stats::setNames(as.list(map[types]), vapply(table$columns, function(col) col$name, character(1)))
  )
}

#' Read one conformance table with the schema's column types
#'
#' The types come from the schema, never from readr's guess; the header must
#' name the same columns in the same order or the read aborts.
#' @keywords internal
#' @noRd
trace_pinned_read <- function(dir, table) {
  file <- fs::path(dir, table$file)
  out <- readr::read_csv(
    file,
    col_types = trace_pinned_cols(table),
    progress = FALSE,
    show_col_types = FALSE,
    na = c("NA", "")
  )
  want <- vapply(table$columns, function(col) col$name, character(1))
  if (!identical(names(out), want)) {
    cli::cli_abort(c(
      gettext("{.file {table$file}} does not match the conformance schema."),
      x = gettext("Schema columns: {.val {want}}"),
      i = gettext("File columns: {.val {names(out)}}")
    ))
  }
  out
}

#' Dot paths of the manifest values jsonlite renders inline
#'
#' `jsonlite::toJSON()` writes an atomic vector inline (`[1, 2]`) and a list
#' multiline, and an `I()` length-1 vector inline as a one-element array.  That
#' distinction lives in the R object, not in its JSON, so a second writer needs
#' it named explicitly.
#' @keywords internal
#' @noRd
trace_inline_arrays <- function(x, prefix = NULL) {
  if (!is.list(x) || is.null(names(x))) {
    return(character(0))
  }
  out <- character(0)
  for (nm in names(x)) {
    value <- x[[nm]]
    path <- if (is.null(prefix)) nm else paste0(prefix, ".", nm)
    inline <- (is.atomic(value) && (length(value) != 1L || inherits(value, "AsIs"))) ||
      (is.list(value) && inherits(value, "AsIs"))
    if (inline) {
      out <- c(out, path)
    } else if (is.list(value) && !is.data.frame(value) && length(value)) {
      out <- c(out, trace_inline_arrays(value, path))
    }
  }
  out
}

# One table cell as a JSON array element.  A one-element list keeps the cell an
# array under `auto_unbox`; a missing value becomes `null`, which is `NA` in the
# CSV and `null` in the JSON.
trace_payload_cell <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (length(x) != 1L || is.na(x)) {
    return(list(NULL))
  }
  list(unname(x))
}

#' Payload a second writer needs to reproduce an artefact byte for byte
#'
#' The payload is the transfer format between [marker_artefact()] and
#' `inst/python/marker_trace.py`; it is not written into an artefact.  It holds
#' the in-memory manifest (so an `AsIs` length-1 vector and a scalar stay
#' distinguishable), the dot paths jsonlite renders inline, and the sidecar rows
#' already read with the conformance schema's types.
#'
#' Write it with `jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null",
#' na = "null", digits = 17)`; `digits = 17` matters, because `digits = NA`
#' rounds to 15 significant digits and the CSV would then differ.
#' @param artefact a `marker_artefact` from [marker_artefact()].  Its sidecar
#'   files are read from `artefact$dir`.
#' @return a list with `schema`, `trace_schema`, `manifest`, `inline_arrays` and
#'   `tables`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' solved <- marker_solve(marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm"))
#' payload <- marker_trace_payload(marker_artefact(solved, dir = tempfile("artefact")))
#' names(payload)
#' names(payload$tables)
#' @family plan
#' @export
marker_trace_payload <- function(artefact) {
  if (!inherits(artefact, "marker_artefact")) {
    cli::cli_abort(gettext(
      "`artefact` must be a marker_artefact; a directory has already lost the in-memory manifest."
    ))
  }
  schema <- trace_conformance_schema()
  tables <- list()
  for (table in schema$tables) {
    if (!fs::file_exists(fs::path(artefact$dir, table$file))) {
      next
    }
    tab <- trace_pinned_read(artefact$dir, table)
    rows <- lapply(seq_len(nrow(tab)), function(i) {
      lapply(names(tab), function(nm) trace_payload_cell(tab[[nm]][[i]]))
    })
    tables[[table$name]] <- list(
      file = table$file,
      columns = table$columns,
      rows = rows
    )
  }
  list(
    schema = "tidymarkers/trace-payload@1",
    trace_schema = TRACE_SCHEMA,
    manifest = artefact$manifest,
    inline_arrays = I(trace_inline_arrays(artefact$manifest)),
    tables = tables
  )
}
