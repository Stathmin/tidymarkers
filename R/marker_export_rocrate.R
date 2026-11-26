# =============================================================================
#  marker_export_rocrate.R  --  the artefact as an RO-Crate (ZY23)
# =============================================================================
#  RO-Crate and PROV-O are the two vocabularies the methodology draft promises,
#  and the artefact is the record they would describe.  An artefact alone is a
#  directory of sidecars; a crate is that same record with a metadata descriptor
#  a WorkflowHub or Galaxy reader can follow, so the design can be published
#  rather than only replayed on the machine that made it.
#
#  Three rules this file keeps deliberately:
#    * the crate is self-contained: every file the manifest names is copied, and
#      a file the manifest names but the artefact lacks is an error rather than a
#      dangling reference;
#    * nothing else in the artefact directory is swept in - a stray file there is
#      not part of the record, and promoting it silently would make the crate
#      claim more than the artefact does;
#    * the import verifies the hashes the crate declares before it returns a
#      `marker_artefact`, so a tampered sidecar is named rather than read.
#
#  PROV terms are written as full IRIs, not through a `prov:` prefix.  The
#  RO-Crate 1.1 context does not map that prefix, and a JSON-LD consumer that
#  does not fetch the context has only the IRI to go on; shorthand there would
#  mean nothing.  RO-Crate's own terms (File, Dataset, CreateAction) stay in
#  their shorthand form, which is what the context and every crate in the wild
#  use.
#
#  "Decision graph" here means the run's decisions - the manifest plus
#  `decisions.csv` - not `marker_decision_graph(scheme)`, which is a criteria DAG
#  derived from code and carries no artefact.  `marker_import_rocrate()` returns
#  a `marker_artefact` for that reason.
#
#  No official RO-Crate validator ships in this tree (they are Python/Node), so
#  the tests check the JSON-LD structurally and the round-trip by hashing; that
#  is a structural check, not a validation, and the test says so.
# =============================================================================

PROV_NS <- "http://www.w3.org/ns/prov#"
ROCRATE_PROFILE <- "https://w3id.org/ro/crate/1.1"
ROCRATE_CONTEXT <- "https://w3id.org/ro/crate/1.1/context"

#' Files an artefact's crate carries
#'
#' The manifest names the sidecars, and a crate that only linked to them would
#' not travel; every non-`NA` entry is therefore copied, alongside the manifest
#' itself.
#' @param manifest an artefact manifest.
#' @return a character vector of file names, `artefact.json` first.
#' @keywords internal
#' @noRd
rocrate_files <- function(manifest) {
  files <- unlist(manifest$files %||% list(), use.names = FALSE)
  files <- as.character(files)
  files <- files[!is.na(files) & nzchar(files)]
  unique(c("artefact.json", files))
}

#' The media type of a sidecar
#' @keywords internal
#' @noRd
rocrate_encoding <- function(file) {
  switch(tools::file_ext(file),
    csv = "text/csv",
    json = "application/json",
    rds = "application/x-r-serialize",
    "application/octet-stream"
  )
}

#' One `File` entity, hashed from the crate's own copy
#'
#' The hash is taken from the copy, not from the artefact directory, so the crate
#' declares what it actually carries.
#' @keywords internal
#' @noRd
rocrate_file_entity <- function(dir, file, run = "#run") {
  path <- fs::path(dir, file)
  entity <- list(
    "@id" = file,
    "@type" = "File",
    "sha256" = hash_file(path),
    "contentSize" = as.integer(file.size(path)),
    "encodingFormat" = rocrate_encoding(file)
  )
  entity[[paste0(PROV_NS, "wasGeneratedBy")]] <- list("@id" = run)
  entity
}

#' One input as a PROV entity
#'
#' An input that came from a file keeps its path as `contentUrl`: the file itself
#' is not copied, because a reference genome is not part of the design record,
#' and the hash recorded beside it is what a reader re-checks.
#' @keywords internal
#' @noRd
rocrate_input_entity <- function(input, i) {
  prov_entity <- paste0(PROV_NS, "Entity")
  entity <- list(
    "@id" = paste0("#input-", i),
    "@type" = prov_entity,
    "name" = paste(input$role, " (", input$source, ")", sep = "")
  )
  path <- input$path %||% NA_character_
  if (is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)) {
    entity[["@type"]] <- c("File", prov_entity)
    entity[["contentUrl"]] <- path
  }
  sha <- input$sha256 %||% NA_character_
  if (is.character(sha) && length(sha) == 1L && !is.na(sha) && nzchar(sha)) {
    entity[["sha256"]] <- sha
  }
  bytes <- input$bytes %||% NA_integer_
  if (is.numeric(bytes) && length(bytes) == 1L && !is.na(bytes)) {
    entity[["contentSize"]] <- as.integer(bytes)
  }
  entity
}

#' The JSON-LD graph of one artefact
#' @param artefact a `marker_artefact`.
#' @param dir the crate directory the sidecars were copied into.
#' @return the graph as a list, `@context` included.
#' @keywords internal
#' @noRd
rocrate_graph <- function(artefact, dir) {
  m <- artefact$manifest
  files <- rocrate_files(m)
  inputs <- lapply(seq_along(m$inputs %||% list()), function(i) rocrate_input_entity(m$inputs[[i]], i))
  input_ids <- lapply(seq_along(inputs), function(i) list("@id" = paste0("#input-", i)))

  descriptor <- list(
    "@id" = "ro-crate-metadata.json",
    "@type" = "CreativeWork",
    "conformsTo" = list("@id" = ROCRATE_PROFILE),
    "about" = list("@id" = "./")
  )

  root <- list(
    "@id" = "./",
    "@type" = "Dataset",
    "name" = "tidymarkers design artefact",
    "identifier" = m$id,
    "datePublished" = m$created,
    "conformsTo" = list("@id" = ROCRATE_PROFILE),
    "hasPart" = lapply(files, function(file) list("@id" = file))
  )

  run <- list(
    "@id" = "#run",
    "@type" = c("CreateAction", paste0(PROV_NS, "Activity")),
    "name" = "marker_solve design run",
    "endTime" = m$created,
    "instrument" = list("@id" = "#howto"),
    "result" = lapply(files, function(file) list("@id" = file))
  )
  run[[paste0(PROV_NS, "startedAtTime")]] <- m$created
  run[[paste0(PROV_NS, "wasAssociatedWith")]] <- list("@id" = "#agent")
  if (length(input_ids)) {
    run[[paste0(PROV_NS, "used")]] <- input_ids
  }

  agent <- list(
    "@id" = "#agent",
    "@type" = c("SoftwareApplication", paste0(PROV_NS, "SoftwareAgent")),
    "name" = m$package$name %||% "tidymarkers",
    "softwareVersion" = m$package$version
  )
  git_sha <- m$package$git_sha %||% NA_character_
  if (is.character(git_sha) && length(git_sha) == 1L && !is.na(git_sha) && nzchar(git_sha)) {
    agent[["identifier"]] <- paste0("git:", git_sha)
  }

  layers <- as.character(m$plan$layers %||% character())
  howto <- list(
    "@id" = "#howto",
    "@type" = "HowTo",
    "name" = "marker plan layers"
  )
  if (length(layers)) {
    howto[["step"]] <- lapply(seq_along(layers), function(i) {
      list("@type" = "HowToStep", "position" = i, "name" = layers[[i]])
    })
  }

  graph <- c(
    list(descriptor, root, run, agent, howto),
    lapply(files, function(file) rocrate_file_entity(dir, file)),
    inputs
  )
  list("@context" = ROCRATE_CONTEXT, "@graph" = graph)
}

#' Write an artefact as an RO-Crate
#'
#' The crate is a self-contained directory: the files the manifest names are
#' copied into `dir` beside `ro-crate-metadata.json`, and nothing else from the
#' artefact directory is taken.  The metadata describes one run as a
#' `CreateAction` / PROV-O `Activity` and one PROV-O `Entity` per hashed input, so
#' a published crate carries the same claim the artefact does - which reference,
#' which variant table, which parameters, which verdicts.
#'
#' @param artefact a `marker_artefact` (see [marker_artefact()]).
#' @param dir the crate directory to write; it is created, and must not be the
#'   artefact's own directory.
#' @return the crate directory, as a length-one character vector.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' plan <- marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm")
#' art <- marker_artefact(marker_solve(plan, progress = FALSE), dir = tempfile("artefact"))
#' crate <- marker_export_rocrate(art, dir = tempfile("crate"))
#' list.files(crate)
#' @family plan
#' @export
marker_export_rocrate <- function(artefact, dir = tempfile("roc-crate")) {
  if (!inherits(artefact, "marker_artefact")) {
    abort_msg("`artefact` must be a marker_artefact (see marker_artefact()).")
  }
  src <- artefact_dir(artefact)
  if (identical(fs::path_abs(dir), fs::path_abs(src))) {
    abort_msg("`dir` must be a new directory, not the artefact's own directory.")
  }
  fs::dir_create(dir)
  files <- rocrate_files(artefact$manifest)
  for (file in files) {
    from <- fs::path(src, file)
    if (!fs::file_exists(from)) {
      abort_msg(
        "The manifest names {.file {file}}, but the artefact does not contain it.",
        i = "A crate that claims a file it does not carry is not self-contained."
      )
    }
    fs::file_copy(from, fs::path(dir, file), overwrite = TRUE)
  }
  graph <- rocrate_graph(artefact, dir)
  writeLines(
    jsonlite::toJSON(
      graph,
      auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA
    ),
    fs::path(dir, "ro-crate-metadata.json")
  )
  as.character(dir)
}

#' Read an artefact back from an RO-Crate
#'
#' Every `sha256` the crate declares is verified against the file on disk, and
#' the root `identifier` must equal the `id` in `artefact.json`; a mismatch aborts
#' naming the file rather than returning a record that no longer matches its
#' crate.  The descriptor must declare the RO-Crate 1.1 context this reader
#' understands.
#'
#' @param path the crate directory, or the `ro-crate-metadata.json` inside it.
#' @return a `marker_artefact`; see [marker_import_rocrate()]'s counterpart
#'   [marker_export_rocrate()].
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' plan <- marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm")
#' art <- marker_artefact(marker_solve(plan, progress = FALSE), dir = tempfile("artefact"))
#' back <- marker_import_rocrate(marker_export_rocrate(art))
#' identical(back$id, art$id)
#' @family plan
#' @export
marker_import_rocrate <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    abort_msg("`path` must be one crate directory or its ro-crate-metadata.json.")
  }
  dir <- if (fs::is_dir(path)) path else fs::path_dir(path)
  meta_path <- fs::path(dir, "ro-crate-metadata.json")
  if (!fs::file_exists(meta_path)) {
    abort_msg("No ro-crate-metadata.json at {.path {meta_path}}.")
  }
  meta <- jsonlite::fromJSON(readr::read_file(meta_path), simplifyVector = FALSE)
  context <- unlist(meta[["@context"]] %||% "", use.names = FALSE)
  if (!any(grepl("w3id.org/ro/crate/1\\.1", context))) {
    abort_msg(
      "Not an RO-Crate 1.1 metadata document.",
      x = "@context is {.val {paste(context, collapse = ', ')}}."
    )
  }
  graph <- meta[["@graph"]] %||% list()
  ids <- vapply(graph, function(entry) as.character(entry[["@id"]] %||% NA_character_), character(1))

  root_i <- match("./", ids)
  if (is.na(root_i)) {
    abort_msg("The crate has no root data entity (`\"@id\": \"./\"`).")
  }

  for (i in seq_along(graph)) {
    id <- ids[[i]]
    declared <- graph[[i]][["sha256"]]
    if (is.null(declared) || is.na(id) || startsWith(id, "#") || grepl("^[A-Za-z][A-Za-z0-9+.-]*:", id)) {
      next
    }
    file <- fs::path(dir, id)
    if (!fs::file_exists(file)) {
      abort_msg(
        "The crate declares {.file {id}}, but the file is not in {.path {dir}}.",
        i = "A crate is a directory, not a manifest of what once was there."
      )
    }
    actual <- hash_file(file)
    if (!identical(as.character(declared), actual)) {
      abort_msg(
        "{.file {id}} does not hash to what the crate declares.",
        x = "Declared sha256: {.val {declared}}",
        i = "File sha256: {.val {actual}}"
      )
    }
  }

  artefact <- marker_artefact_read(dir)
  declared_id <- graph[[root_i]][["identifier"]]
  if (is.null(declared_id) || !identical(as.character(declared_id), as.character(artefact$id))) {
    abort_msg(
      "The crate's root identifier does not match artefact.json.",
      x = "RO-Crate identifier: {.val {declared_id}}",
      i = "artefact.json id: {.val {artefact$id}}"
    )
  }
  artefact
}
