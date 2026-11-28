# =============================================================================
#  marker_export_outcomes.R  --  outcomes handed over without a network (ZY27)
# =============================================================================
#  An aggregate is not an anonymisation, and a crate that carries the run but
#  not the rows is not a handover.  This file turns a corpus (ZY26) into one
#  RO-Crate per participating artefact, so a laboratory can give a data trustee
#  the observations without giving it the identity behind them.
#
#  Three rules this file keeps deliberately:
#    * identity never travels with the data: the artefact directory name is
#      replaced by a salted pseudonym, and the mapping back to `artefact_id`
#      is written beside the export, outside every crate;
#    * provenance survives the pseudonym: each outcome row is a PROV-O `Entity`
#      derived from the crate's root dataset, and the inputs keep their hashes,
#      so the trustee can still re-check what the row was measured against;
#    * no network, no absolute paths: nothing here opens a connection, and the
#      sender's disk layout is rewritten to file names before the metadata is
#      written, because a crate that carries `/home/...` leaks the sender;
#      the hash beside an input still identifies it.
#
#  The crates round-trip through `marker_import_rocrate()` because each one is a
#  real artefact crate first and an outcomes crate second.  No official RO-Crate
#  validator ships in this tree, so the test checks the JSON-LD structurally and
#  the round-trip by hashing; that is a structural check, not a validation, and
#  the test says so.
# =============================================================================

#' Attach outcome rows to an artefact crate
#'
#' The crate is an artefact crate; this adds the pseudonymised `outcomes.csv`
#' and one PROV-O `Entity` per recorded row.  Absolute `contentUrl` values are
#' rewritten to file names first, so the crate cannot carry the sender's disk
#' layout.
#' @param crate the crate directory, already written by [marker_export_rocrate()].
#' @param rows the corpus rows belonging to this artefact.
#' @return the crate directory, invisibly.
#' @keywords internal
#' @noRd
outcomes_attach <- function(crate, rows) {
  meta_path <- fs::path(crate, "ro-crate-metadata.json")
  meta <- jsonlite::fromJSON(readr::read_file(meta_path), simplifyVector = FALSE)
  graph <- meta[["@graph"]]
  ids <- vapply(graph, function(entry) as.character(entry[["@id"]] %||% NA_character_), character(1))

  for (i in seq_along(graph)) {
    url <- graph[[i]][["contentUrl"]]
    if (is.character(url) && length(url) == 1L && !is.na(url) && fs::is_absolute_path(url)) {
      graph[[i]][["contentUrl"]] <- as.character(fs::path_file(url))
    }
  }

  if (!"outcomes.csv" %in% ids) {
    prov_entity <- paste0(PROV_NS, "Entity")
    row_ids <- lapply(seq_len(nrow(rows)), function(i) list("@id" = paste0("#outcome-", i)))
    entities <- lapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i, , drop = FALSE]
      entity <- list(
        "@id" = paste0("#outcome-", i),
        "@type" = prov_entity,
        "schema" = "tidymarkers/outcomes@1",
        "target" = as.character(row$target[[1L]]),
        "result" = as.character(row$result[[1L]]),
        "recorded" = as.character(row$recorded[[1L]]),
        "notes" = as.character(row$notes[[1L]]),
        "calibrated" = as.logical(row$calibrated[[1L]])
      )
      entity[[paste0(PROV_NS, "wasDerivedFrom")]] <- list("@id" = "./")
      entity[[paste0(PROV_NS, "wasGeneratedBy")]] <- list("@id" = "#run")
      entity
    })
    index <- list(
      "@id" = "#outcome-index",
      "@type" = "Dataset",
      "name" = "recorded outcomes",
      "schema" = "tidymarkers/outcomes@1",
      "hasPart" = row_ids
    )
    graph <- c(graph, list(rocrate_file_entity(crate, "outcomes.csv")), list(index), entities)
    root_i <- match("./", ids)
    graph[[root_i]][["hasPart"]] <- c(
      graph[[root_i]][["hasPart"]] %||% list(),
      list(list("@id" = "outcomes.csv"))
    )
    run_i <- match("#run", ids)
    graph[[run_i]][["result"]] <- c(
      graph[[run_i]][["result"]] %||% list(),
      list(list("@id" = "outcomes.csv"))
    )
  }

  meta[["@graph"]] <- graph
  writeLines(
    jsonlite::toJSON(
      meta,
      auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA
    ),
    meta_path
  )
  invisible(crate)
}

#' Export outcomes for a data trustee without a network
#'
#' [marker_outcome_corpus()] answers "what did the lab get"; this writes that
#' answer as crates a trustee can read.  One crate per artefact directory, each
#' a real artefact crate with the recorded rows attached, and a mapping from
#' pseudonym to `artefact_id` written beside the export so the owner keeps the
#' identity the trustee does not get.
#'
#' The pseudonym is `site-` plus eight hex characters of
#' `sha256(salt + artefact directory name)`.  `trustee` is the salt: the same
#' trustee over the same corpus gives the same pseudonyms, and a different salt
#' gives different ones.  Nothing is sent anywhere; the export is files on disk.
#'
#' @param corpus a corpus from [marker_outcome_corpus()], or a directory to scan
#'   with it.  A hand-built table has no scanned root and is refused.
#' @param dir the directory to write `<dir>/crates/<pseudonym>/` and
#'   `<dir>-identity.csv` into.
#' @param trustee an optional length-one salt for the pseudonyms.  When `NULL`
#'   a random 16-hex salt is used, so two exports are not linkable.
#' @return the identity map, invisibly: `site`, `pseudonym`, `artefact_id`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' solved <- marker_solve(marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm"))
#' root <- tempfile("corpus")
#' art <- marker_artefact(solved, dir = file.path(root, "batch-a"))
#' marker_record_outcome(art, "works")
#' handover <- tempfile("handover")
#' map <- marker_export_outcomes(root, handover, trustee = "lab-a")
#' map$pseudonym
#' @family plan
#' @export
marker_export_outcomes <- function(corpus, dir, trustee = NULL) {
  if (is.character(corpus)) {
    if (length(corpus) != 1L || is.na(corpus)) {
      abort_msg("A character `corpus` must be one directory path.")
    }
    corpus <- marker_outcome_corpus(corpus)
  }
  if (!is.data.frame(corpus)) {
    abort_msg("`corpus` must be a corpus table or a directory.")
  }
  root <- attr(corpus, "root")
  if (!is.character(root) || length(root) != 1L || is.na(root) || !nzchar(root)) {
    abort_msg(
      "`corpus` carries no scanned root.",
      i = "Build it with marker_outcome_corpus(), or pass the directory instead."
    )
  }
  root <- fs::path_abs(root)
  if (!is.character(dir) || length(dir) != 1L || is.na(dir) || !nzchar(dir)) {
    abort_msg("`dir` must be one directory path.")
  }
  if (identical(fs::path_abs(dir), root) || fs::path_has_parent(fs::path_abs(dir), root)) {
    abort_msg(
      "`dir` must sit outside the scanned corpus root.",
      i = "An export inside it would be read back as more corpus on the next scan."
    )
  }

  salt <- if (is.null(trustee)) {
    paste(sprintf("%02x", sample(0:255, 8, replace = TRUE)), collapse = "")
  } else {
    if (!is.character(trustee) || length(trustee) != 1L || is.na(trustee) || !nzchar(trustee)) {
      abort_msg("`trustee` must be NULL or one non-empty string.")
    }
    trustee
  }

  sites <- unique(as.character(corpus$artefact))
  # shortcut: the pseudonym is salted by the directory name alone; two artefact
  # directories with the same basename in different subtrees would collide, and
  # the guard below refuses rather than overwriting a crate silently.
  pseudonyms <- vapply(sites, function(rel) {
    paste0("site-", substr(hash_object(paste0(salt, fs::path_file(rel))), 1L, 8L))
  }, character(1), USE.NAMES = FALSE)
  if (anyDuplicated(pseudonyms)) {
    dup <- unique(pseudonyms[duplicated(pseudonyms)])
    cli::cli_abort(c(
      gettextf("Two artefact directories share a pseudonym: %s.", paste(dup, collapse = ", ")),
      i = gettext("Give them distinct names, or export them in separate runs.")
    ))
  }

  artefact <- as.character(corpus$artefact)
  artefact_id <- as.character(corpus$artefact_id)
  map <- tibble::tibble(
    site = as.character(sites),
    pseudonym = pseudonyms,
    artefact_id = vapply(sites, function(rel) {
      artefact_id[match(rel, artefact)]
    }, character(1), USE.NAMES = FALSE)
  )

  identity_path <- paste0(dir, "-identity.csv")
  refuse_links(fs::path_dir(identity_path), fs::path_file(identity_path))

  fs::dir_create(dir)
  for (i in seq_along(sites)) {
    rel <- sites[[i]]
    adir <- fs::path(root, rel)
    if (!identical(fs::path_abs(adir), root) && !fs::path_has_parent(fs::path_abs(adir), root)) {
      abort_msg(
        "The corpus names an artefact outside its scanned root: {.val {rel}}.",
        i = "A corpus row describes a directory under {.path {root}}."
      )
    }
    adir <- artefact_member_path(root, rel)
    crate <- marker_export_rocrate(
      marker_artefact_read(adir),
      fs::path(dir, "crates", pseudonyms[[i]])
    )
    refuse_links(crate, "outcomes.csv")
    fs::file_copy(
      artefact_member_path(adir, "outcomes.csv"),
      fs::path(crate, "outcomes.csv"),
      overwrite = TRUE
    )
    outcomes_attach(crate, corpus[corpus$artefact == rel, , drop = FALSE])
  }

  readr::write_csv(map, identity_path)
  invisible(map)
}
