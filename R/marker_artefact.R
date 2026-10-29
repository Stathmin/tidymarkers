# =============================================================================
#  marker_artefact.R  --  the self-documenting design artefact (YP3)
# =============================================================================
#  "We log everything" is a claim; an artefact is the thing that makes it
#  checkable.  `marker_artefact()` writes, beside the solved plan, the three
#  things a later reader needs in order to disagree with the result:
#
#    * inputs     -- which reference and which variant table, hashed;
#    * provenance -- package version, git SHA, engine versions, every parameter
#                    and, separately, the ones that differ from the default;
#    * decisions  -- one row per (target, chemistry) tried: the verdict, the
#                    risk, the relaxed gates and, for a rejection, the reason
#                    the designer gave.
#
#  `marker_replay()` re-runs the design from the RDS and compares, so a recorded
#  verdict is a claim that can be re-tested rather than trusted.  The RDS stays
#  the faithful copy; `artefact.json` and the CSVs are the sidecars a non-R
#  consumer can read.
#
#  Two rules this file keeps deliberately:
#    * an unavailable engine or an unknown default is written as NA, never
#      dropped - a manifest that silently omits BLAST looks exactly like one
#      taken on a machine that never had it (the stamp contract, ZW2);
#    * an override is a field of its own, not a value buried in `parameters`,
#      so a "pass" cannot hide a loosened threshold.
# =============================================================================

# --- hashing and provenance --------------------------------------------------

#' SHA-256 of a file
#' @param path a file path.
#' @return a length-one character vector, `NA` when the file is missing.
#' @keywords internal
#' @noRd
artefact_sha256_file <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NA_character_)
  }
  tryCatch(
    as.character(digest::digest(file = path, algo = "sha256")),
    error = function(e) NA_character_
  )
}

#' SHA-256 of an R object
#'
#' Kept apart from `artefact_sha256_file()` because the same character value is
#' a path in one call and data in the other, and guessing wrong hashes nothing.
#' @param x any R object to serialise and hash.
#' @return a length-one character vector, `NA` when hashing fails.
#' @keywords internal
#' @noRd
artefact_sha256 <- function(x) {
  tryCatch(
    as.character(digest::digest(x, algo = "sha256")),
    error = function(e) NA_character_
  )
}

#' The package source root, when there is one
#'
#' An installed package ships neither `.git` nor `renv.lock`, and `pkgload`
#' points `system.file()` at `<root>/inst`, so the root is the parent there.
#' @return a path, or `NA`.
#' @keywords internal
#' @noRd
artefact_pkg_root <- function() {
  dir <- suppressWarnings(system.file(package = "tidymarkers"))
  if (!nzchar(dir)) {
    return(NA_character_)
  }
  if (identical(basename(dir), "inst")) dirname(dir) else dir
}

#' The git commit of the source tree, when there is one
#'
#' An installed package has no `.git`, and inventing a SHA would be worse than
#' none: the field is `NA` there and the manifest says so.
#' @return a length-one character vector, or `NA`.
#' @keywords internal
#' @noRd
artefact_git_sha <- function() {
  sha <- utils::packageDescription("tidymarkers")[["RemoteSha"]]
  if (!is.null(sha) && nzchar(sha)) {
    return(as.character(sha))
  }
  root <- artefact_pkg_root()
  if (is.na(root) || !dir.exists(file.path(root, ".git"))) {
    return(NA_character_)
  }
  git <- Sys.which("git")
  if (!nzchar(git)) {
    return(NA_character_)
  }
  out <- tryCatch(
    suppressWarnings(system2(git, c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (!length(out) || !nzchar(out[[1L]])) NA_character_ else trimws(out[[1L]])
}

#' A value the manifest can carry
#'
#' JSON has no `Inf`, no functions and no environments, so anything the design
#' knobs could hold is narrowed to a scalar or a short array here.  A value that
#' cannot travel is written as its deparsed text rather than dropped - the
#' manifest must name what was passed, if not in its native type.
#' @keywords internal
#' @noRd
artefact_value <- function(x) {
  if (is.null(x)) {
    return(NA)
  }
  if (is.function(x) || is.environment(x) || is.language(x)) {
    return(paste(deparse(x), collapse = " "))
  }
  if (is.numeric(x) && any(!is.finite(x))) {
    return(paste(format(x), collapse = ", "))
  }
  if (is.atomic(x) && length(x) <= 20L) {
    return(x)
  }
  paste(utils::capture.output(utils::str(x)), collapse = " ")
}

#' The reference a plan designs against, as an hashed manifest row
#' @keywords internal
#' @noRd
artefact_genome_input <- function(genome) {
  if (inherits(genome, "marker_genome")) {
    if (!is.null(genome$bsgenome)) {
      return(list(
        role = "genome", source = "BSgenome", name = class(genome$bsgenome)[1L],
        path = NA_character_, sha256 = NA_character_, bytes = NA_integer_
      ))
    }
    path <- genome$path
    return(list(
      role = "genome", source = "fasta", name = fs::path_file(path),
      path = as.character(fs::path_abs(path)), sha256 = artefact_sha256_file(path),
      bytes = as.integer(file.size(path))
    ))
  }
  if (inherits(genome, "DNAStringSet")) {
    seqs <- as.character(genome)
    return(list(
      role = "genome", source = "DNAStringSet",
      name = paste0(length(seqs), " sequence(s)"),
      path = NA_character_,
      sha256 = artefact_sha256(paste(seqs, collapse = "")),
      bytes = as.integer(sum(nchar(seqs)))
    ))
  }
  list(
    role = "genome", source = class(genome)[1L], name = NA_character_,
    path = NA_character_, sha256 = NA_character_, bytes = NA_integer_
  )
}

#' The hashed inputs a plan read
#' @keywords internal
#' @noRd
artefact_inputs <- function(plan) {
  out <- list(artefact_genome_input(plan$genome))
  if (!is.null(plan$variants)) {
    out[[length(out) + 1L]] <- list(
      role = "variants", source = "table",
      name = paste0(nrow(plan$variants), " variant(s)"),
      path = NA_character_, sha256 = artefact_sha256(plan$variants),
      bytes = NA_integer_
    )
  }
  if (!is.null(plan$haplotypes)) {
    out[[length(out) + 1L]] <- list(
      role = "haplotypes", source = "table",
      name = paste0(nrow(plan$haplotypes), " haplotype row(s)"),
      path = NA_character_, sha256 = artefact_sha256(plan$haplotypes),
      bytes = NA_integer_
    )
  }
  out
}

# --- parameters and overrides ------------------------------------------------

#' The designer a scheme dispatches to, for reading its defaults
#'
#' Only used to answer "was this knob left at its default?"; a scheme this table
#' does not know simply reports `default = NA`, which reads as "unknown" rather
#' than as "unchanged".
#' @keywords internal
#' @noRd
artefact_designer <- function(scheme) {
  switch(scheme,
    snp = ,
    indel = ,
    multi = assemble_kasp_assay,
    caps = design_caps_marker,
    arms = design_arms_marker,
    hrm = design_hrm_marker,
    taqman = design_taqman_marker,
    dpcr = design_dpcr_marker,
    rpa = design_rpa_marker,
    extra_scheme_designer(scheme)
  )
}

#' The evaluated default of one designer argument
#' @keywords internal
#' @noRd
artefact_default <- function(fn, arg) {
  if (!is.function(fn) || !arg %in% names(formals(fn))) {
    return(NA)
  }
  d <- formals(fn)[[arg]]
  if (identical(d, quote(expr = ))) {
    return(NA)
  }
  tryCatch(artefact_value(eval(d, envir = environment(fn))), error = function(e) NA)
}

#' Knobs passed to the solver that are not the designers' defaults
#'
#' A knob is compared only against the schemes that declare it as a formal;
#' every other scheme in the run knows nothing about it, and treating that as an
#' override would flood the field with false flags.  A knob that **no** scheme
#' that ran declares is still listed, with `scheme` and `default` `NA`: passing
#' an argument nothing reads is exactly the sort of thing a manifest should
#' surface.
#' @param knobs the solver's `...` (named list).
#' @param schemes the schemes that were tried.
#' @return a tibble with `scheme`, `argument`, `value`, `default`, `overridden`.
#' @keywords internal
#' @noRd
artefact_overrides <- function(knobs, schemes) {
  empty <- tibble::tibble(
    scheme = character(0), argument = character(0), value = list(),
    default = list(), overridden = logical(0)
  )
  if (!length(knobs)) {
    return(empty)
  }
  declare <- function(sc, nm) {
    fn <- artefact_designer(sc)
    is.function(fn) && nm %in% names(formals(fn))
  }
  rows <- purrr::imap(knobs, function(v, nm) {
    declared <- purrr::keep(schemes, function(sc) declare(sc, nm))
    if (!length(declared)) {
      return(tibble::tibble(
        scheme = NA_character_, argument = nm, value = list(artefact_value(v)),
        default = list(NA), overridden = NA
      ))
    }
    out <- purrr::map(declared, function(sc) {
      d <- artefact_default(artefact_designer(sc), nm)
      tibble::tibble(
        scheme = sc, argument = nm, value = list(artefact_value(v)),
        default = list(d),
        overridden = if (length(d) == 1L && is.na(d)) NA else !identical(artefact_value(v), d)
      )
    }) |>
      purrr::list_rbind()
    out[is.na(out$overridden) | out$overridden, , drop = FALSE]
  }) |>
    purrr::list_rbind()
  if (is.null(rows)) empty else rows
}

# --- the decision log --------------------------------------------------------

#' `x` joined for a CSV cell, `NA` when there is nothing to say
#' @keywords internal
#' @noRd
artefact_join <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) NA_character_ else paste(x, collapse = "; ")
}

#' What one (target, chemistry) attempt decided
#'
#' The point of the log is that a rejection is never bare: a `marker_doomed`
#' carries the stage and the reason the designer gave, and a designable marker
#' carries the verdict, the risk and the gates the selection ladder relaxed.
#' @param m the designer's return for this attempt (may be `NULL` when it erred).
#' @return a one-row list.
#' @keywords internal
#' @noRd
artefact_candidate <- function(m) {
  out <- list(
    outcome = "rejected", verdict = NA_character_, risk = NA_character_,
    n_failed = NA_integer_, failed = NA_character_, relaxed = NA_character_,
    violated = NA_character_, stage = NA_character_, reason = NA_character_,
    discrimination = NA_real_
  )
  if (is.null(m)) {
    out$reason <- gettext("the designer returned no marker object")
    out$stage <- "design"
    return(out)
  }
  if (inherits(m, "marker_doomed")) {
    out$reason <- as.character(m$reason)
    out$stage <- as.character(m$stage)
    out$verdict <- "doomed"
    return(out)
  }
  out$verdict <- tryCatch(evaluation_verdict(m), error = function(e) NA_character_)
  rk <- tryCatch(marker_risk(m), error = function(e) NULL)
  if (!is.null(rk) && nrow(rk)) {
    out$risk <- as.character(rk$risk[[1L]])
    out$n_failed <- as.integer(rk$n_failed[[1L]])
    out$failed <- artefact_join(rk$failed[[1L]])
    out$discrimination <- as.numeric(rk$discrimination[[1L]])
  }
  sel <- m$selection
  if (is.list(sel)) {
    out$relaxed <- artefact_join(sel$relaxed)
    out$violated <- artefact_join(sel$violated)
  }
  out$outcome <- if (marker_ok(m)) "designable" else "rejected"
  if (identical(out$outcome, "rejected") && is.na(out$reason)) {
    out$reason <- if (is.na(out$failed)) {
      gettext("a gate rejected the candidate")
    } else {
      out$failed
    }
  }
  out
}

#' The cached design behind one attempt
#' @keywords internal
#' @noRd
artefact_lookup <- function(solved, label, scheme) {
  cache <- solved$cache
  if (!is.environment(cache)) {
    return(NULL)
  }
  tg <- solved$plan$targets
  i <- match(label, tg$label)
  if (is.na(i)) {
    return(NULL)
  }
  key <- paste(
    target_key(tg[i, , drop = FALSE]), scheme,
    knobs_key(solved$knobs %||% list()),
    sep = "||"
  )
  get0(key, envir = cache, inherits = FALSE)
}

#' One row per (target, chemistry) attempt, with its decision
#' @keywords internal
#' @noRd
artefact_decisions <- function(solved) {
  a <- solved$attempts
  if (!nrow(a)) {
    return(tibble::tibble(
      target = character(0), chrom = character(0), pos = integer(0),
      ref = character(0), alt = character(0), scheme = character(0),
      analytical = character(0), cost = numeric(0), chosen = logical(0),
      outcome = character(0),
      verdict = character(0), risk = character(0), n_failed = integer(0),
      failed = character(0), relaxed = character(0), violated = character(0),
      stage = character(0), reason = character(0), discrimination = numeric(0)
    ))
  }
  tg <- solved$plan$targets
  idx <- match(a$label, tg$label)
  # one report row per target: compare each attempt against its target's winner
  win <- match(a$label, solved$report$label)
  won <- !is.na(win) & !is.na(solved$report$scheme[win]) & solved$report$ok[win]
  won[is.na(won)] <- FALSE
  chosen <- won & a$scheme == solved$report$scheme[win]
  chosen[is.na(chosen)] <- FALSE
  purrr::map(seq_len(nrow(a)), function(k) {
    i <- idx[[k]]
    is_chosen <- chosen[[k]]
    m <- artefact_lookup(solved, a$label[[k]], a$scheme[[k]])
    d <- artefact_candidate(m)
    # Z8.5: the decision names the analytical quantities that must be reported
    # separately for this chemistry (NA when the scheme declares none).
    aq <- names(scheme_analytical(a$scheme[[k]]))
    tibble::tibble(
      target = as.character(a$label[[k]]),
      chrom = if (is.na(i)) NA_character_ else as.character(tg$chrom[[i]]),
      pos = if (is.na(i)) NA_integer_ else as.integer(tg$pos[[i]]),
      ref = if (is.na(i)) NA_character_ else as.character(tg$ref[[i]]),
      alt = if (is.na(i)) NA_character_ else as.character(tg$alt[[i]]),
      scheme = a$scheme[[k]],
      analytical = if (length(aq)) paste(aq, collapse = "; ") else NA_character_,
      cost = a$cost[[k]], chosen = is_chosen,
      outcome = if (is_chosen) "chosen" else d$outcome,
      verdict = d$verdict, risk = d$risk, n_failed = d$n_failed,
      failed = d$failed, relaxed = d$relaxed, violated = d$violated,
      stage = d$stage, reason = d$reason, discrimination = d$discrimination
    )
  }) |>
    purrr::list_rbind()
}

#' The gate table of every chosen marker
#' @keywords internal
#' @noRd
artefact_checks <- function(solved) {
  r <- solved$report
  keep <- which(r$ok & !vapply(r$marker, is.null, logical(1)))
  if (!length(keep)) {
    return(NULL)
  }
  purrr::map(keep, function(i) {
    ev <- r$marker[[i]]$evaluation
    if (is.null(ev)) {
      ev <- tryCatch(
        marker_design_checks(r$marker[[i]]),
        error = function(e) NULL
      )
    }
    if (is.null(ev) || !nrow(ev)) {
      return(NULL)
    }
    ev <- ev[, !vapply(ev, is.list, logical(1)), drop = FALSE]
    ev$target <- as.character(r$label[[i]])
    ev$scheme <- as.character(r$scheme[[i]])
    ev[, c("target", "scheme", setdiff(names(ev), c("target", "scheme"))), drop = FALSE]
  }) |>
    purrr::list_rbind() |>
    tibble::as_tibble()
}

# --- the artefact ------------------------------------------------------------

#' A stable id for the design, independent of when it was written
#' @keywords internal
#' @noRd
artefact_id <- function(manifest) {
  core <- manifest[c(
    "schema", "package", "engine_stamp", "inputs", "plan", "parameters", "overrides"
  )]
  artefact_sha256(jsonlite::toJSON(core, auto_unbox = TRUE, null = "null", digits = NA))
}

#' Write a self-documenting artefact for a solved plan
#'
#' Writes, into `dir`: `artefact.json` (the manifest), `decisions.csv` (one row
#' per target x chemistry tried, with the reason for every rejection),
#' `checks.csv` (the gate table of each chosen marker), `signals.csv` (for a
#' regional plan, each haplotype's readout and whether it holds inside the
#' phased block), `targets.csv` and, when it is not already there, `solved.rds`
#' (the faithful copy [marker_replay()] re-runs).
#'
#' The manifest records what the numbers cannot say for themselves: the hashed
#' reference and variant table, the package version and git SHA, every engine
#' version (an absent engine is `NA`, never omitted), the design knobs, and -
#' as its own field - the knobs that differ from the designer defaults.  An
#' overridden threshold is visible without diffing `parameters` against a
#' manual, because a "pass" obtained under a loosened limit must not look like
#' one obtained under the default.
#'
#' `id` is a SHA-256 of the parts that decide the design (not the timestamp), so
#' two runs of the same plan share it and a changed input or parameter does not.
#' @param solved a `marker_solved` (see [marker_solve()]).  Its `cache` is where
#'   the rejected candidates' reasons are read from; a solve that has been
#'   round-tripped through `saveRDS()` keeps it.
#' @param dir output directory (default: a temp dir).
#' @return a `marker_artefact`: `id`, `dir`, `manifest`, `decisions`.
#' @examples
#' \dontrun{
#' plan <- marker_plan(ath_genome()) + marker_targets(marker_example())
#' solved <- marker_solve(plan)
#' art <- marker_artefact(solved, dir = "results/artefact")
#' art
#' # re-run the design from the artefact and check the verdicts agree
#' marker_replay(art)
#' }
#' @family KASP plan
#' @export
marker_artefact <- function(solved, dir = tempdir()) {
  if (!inherits(solved, "marker_solved")) {
    cli::cli_abort(gettext("`solved` must be a marker_solved (see marker_solve())."))
  }
  fs::dir_create(dir)
  pkg_dir <- artefact_pkg_root()
  knobs <- solved$knobs %||% list()
  decisions <- artefact_decisions(solved)
  checks <- artefact_checks(solved)
  signals <- solved$signals %||% NULL

  readr::write_csv(decisions, fs::path(dir, "decisions.csv"))
  if (!is.null(checks)) {
    readr::write_csv(checks, fs::path(dir, "checks.csv"))
  }
  if (!is.null(signals) && nrow(signals)) {
    readr::write_csv(signals, fs::path(dir, "signals.csv"))
  }
  tg <- solved$plan$targets
  readr::write_csv(
    tg[, intersect(c("label", "chrom", "pos", "ref", "alt", "type", "id"), names(tg)), drop = FALSE],
    fs::path(dir, "targets.csv")
  )

  manifest <- list(
    schema = "tidymarkers/artefact@1",
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    package = list(
      name = "tidymarkers",
      version = as.character(utils::packageVersion("tidymarkers")),
      git_sha = artefact_git_sha(),
      renv_lock_sha256 = if (is.na(pkg_dir)) {
        NA_character_
      } else {
        artefact_sha256_file(fs::path(pkg_dir, "renv.lock"))
      }
    ),
    platform = list(
      r = as.character(getRversion()),
      os = R.version$platform,
      locale = Sys.getlocale("LC_COLLATE")
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
    counts = list(
      targets = nrow(solved$report),
      solved = sum(solved$report$ok),
      attempts = nrow(decisions),
      rejected = sum(decisions$outcome == "rejected"),
      signals = if (is.null(signals)) 0L else nrow(signals)
    ),
    files = list(
      decisions = "decisions.csv",
      checks = if (is.null(checks)) NA_character_ else "checks.csv",
      signals = if (is.null(signals) || !nrow(signals)) NA_character_ else "signals.csv",
      targets = "targets.csv",
      solved = "solved.rds"
    )
  )
  manifest$id <- artefact_id(manifest)
  writeLines(
    jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA),
    fs::path(dir, "artefact.json")
  )
  if (!fs::file_exists(fs::path(dir, "solved.rds"))) {
    saveRDS(solved, fs::path(dir, "solved.rds"))
  }
  structure(
    list(
      id = manifest$id, dir = as.character(dir), manifest = manifest,
      decisions = decisions
    ),
    class = "marker_artefact"
  )
}

#' @export
print.marker_artefact <- function(x, ...) {
  m <- x$manifest
  cat(glue::glue("<marker_artefact> {substr(m$id, 1L, 16L)}\n"))
  cat(glue::glue(
    "  {m$counts$targets} target(s), {m$counts$solved} solved, ",
    "{m$counts$attempts} attempt(s), {m$counts$rejected} rejected\n"
  ))
  cat(glue::glue("  package {m$package$version} (git {m$package$git_sha %||% 'NA'})\n"))
  if (length(m$overrides)) {
    cat(glue::glue("  overrides: {paste(m$overrides$argument, collapse = ', ')}\n"))
  }
  cat(glue::glue("  {x$dir}\n"))
  invisible(x)
}

# --- replay ------------------------------------------------------------------

#' The directory an artefact lives in
#' @keywords internal
#' @noRd
artefact_dir <- function(artefact) {
  if (inherits(artefact, "marker_artefact")) {
    return(artefact$dir)
  }
  if (is.character(artefact) && length(artefact) == 1L) {
    return(artefact)
  }
  cli::cli_abort(gettext("`artefact` must be a marker_artefact or a directory."))
}

#' Read an artefact's manifest back
#' @keywords internal
#' @noRd
artefact_manifest <- function(dir) {
  path <- fs::path(dir, "artefact.json")
  if (!fs::file_exists(path)) {
    cli::cli_abort(gettext("No artefact.json in {.path {dir}}; write one with marker_artefact()."))
  }
  jsonlite::fromJSON(readr::read_file(path), simplifyVector = FALSE)
}

#' Do the inputs still hash to what the artefact recorded?
#'
#' Kept separate from the replay itself: a design that no longer matches because
#' the reference was replaced is a different failure from one that no longer
#' matches because an engine drifted, and the caller should hear which.
#' @keywords internal
#' @noRd
artefact_inputs_match <- function(recorded, plan) {
  now <- artefact_inputs(plan)
  if (length(recorded) != length(now)) {
    return(FALSE)
  }
  all(vapply(seq_along(recorded), function(i) {
    a <- recorded[[i]]
    b <- now[[i]]
    identical(a$role, b$role) &&
      (is.null(a$sha256) || is.na(a$sha256) || identical(a$sha256, b$sha256))
  }, logical(1)))
}

#' Does the engine stamp still match what the artefact recorded?
#'
#' A replay that differs cannot say *why* it differs, and a version bump in an
#' engine changes the numbers a design rests on just as a model change does.  The
#' recorded stamp is therefore compared before the verdicts, so a reader is told
#' which engines moved instead of being left to blame the model (ZW2.3).  An
#' engine absent from either stamp counts as a change -- a stamp that dropped an
#' engine is not the same stamp.  A manifest that recorded no stamp at all is
#' `NA` (unknown), never a match.
#' @keywords internal
#' @noRd
artefact_engine_match <- function(recorded, now = marker_engine_stamp()) {
  rec <- unlist(recorded, use.names = TRUE)
  rec <- stats::setNames(as.character(rec), names(rec))
  if (!length(rec)) {
    return(list(match = NA, changed = character()))
  }
  engines <- union(names(rec), names(now))
  changed <- engines[!vapply(engines, function(e) {
    a <- unname(rec[e])
    b <- unname(now[e])
    (is.na(a) && is.na(b)) || identical(a, b)
  }, logical(1))]
  list(match = !length(changed), changed = changed)
}

#' Re-run the design an artefact recorded and compare
#'
#' Reads `solved.rds`, re-solves its plan with the recorded knobs, and compares
#' the chosen chemistry and verdict per target.  This is what turns "we logged
#' everything" into a claim that can be checked: a mismatch names the target and
#' both verdicts instead of silently returning a different design.
#'
#' The input hashes and the engine stamp are re-checked first, so a reference that
#' was replaced is reported as such (`inputs_match = FALSE`) and a library upgrade
#' as an engine change (`engine_match = FALSE`, `engines_changed`), rather than
#' either being read as a model change.
#' @param artefact a `marker_artefact` or the directory it was written to.
#' @param cache an environment to memoise the replay in; pass the original
#'   solve's `cache` to make replay free, `NULL` (the default) to force the
#'   design to run again - which is the only version that tests anything.
#' @param progress show a progress bar (default: not, the replay is a check).
#' @param tolerance accepted; reserved for numeric comparison of the markers.
#' @return a `marker_replay`: `id`, `match` (one logical for the whole run),
#'   `inputs_match`, `engine_match` (`NA` when the artefact recorded no stamp),
#'   `engines_changed` (names whose version differs), and `table` (one row per
#'   target: `target`, `recorded_scheme`, `replayed_scheme`,
#'   `recorded_verdict`, `replayed_verdict`, `match`).
#' @examples
#' \dontrun{
#' art <- marker_artefact(marker_solve(plan), dir = "results/artefact")
#' marker_replay(art) # $match
#' }
#' @family KASP plan
#' @export
marker_replay <- function(artefact, cache = NULL, progress = FALSE, tolerance = NULL) {
  dir <- artefact_dir(artefact)
  manifest <- artefact_manifest(dir)
  path <- fs::path(dir, "solved.rds")
  if (!fs::file_exists(path)) {
    cli::cli_abort(gettext("No solved.rds in {.path {dir}}; the artefact cannot be replayed."))
  }
  solved <- readRDS(path)
  inputs_match <- artefact_inputs_match(manifest$inputs, solved$plan)
  if (isFALSE(inputs_match)) {
    cli::cli_warn(gettext(
      "The reference or variant table no longer matches the artefact's hash; \\
       a replay difference below is explained by the input, not the engine."
    ))
  }
  engines <- artefact_engine_match(manifest$engine_stamp)
  if (isFALSE(engines$match)) {
    cli::cli_warn(c(
      gettext("Engine changed since the artefact was written: {.val {engines$changed}}."),
      "i" = gettext("A replay difference below may be engine drift, not a model change.")
    ))
  }
  knobs <- solved$knobs %||% list()
  again <- rlang::exec(
    marker_solve, solved$plan, !!!knobs,
    cache = cache, progress = isTRUE(progress)
  )
  rec <- solved$report
  rep_ <- again$report
  j <- match(rec$label, rep_$label)
  recorded_verdict <- purrr::map2_chr(rec$ok, rec$marker, artefact_verdict)
  replayed_verdict <- purrr::map2_chr(rep_$ok[j], rep_$marker[j], artefact_verdict)
  replayed_scheme <- rep_$scheme[j]
  tbl <- tibble::tibble(
    target = rec$label,
    recorded_scheme = rec$scheme, replayed_scheme = replayed_scheme,
    recorded_verdict = recorded_verdict, replayed_verdict = replayed_verdict
  )
  tbl$match <- identical_na(tbl$recorded_scheme, tbl$replayed_scheme) &
    identical_na(tbl$recorded_verdict, tbl$replayed_verdict)
  structure(
    list(
      id = manifest$id, dir = dir, match = all(tbl$match),
      inputs_match = isTRUE(inputs_match), engine_match = engines$match,
      engines_changed = engines$changed, table = tbl
    ),
    class = "marker_replay"
  )
}

#' The verdict of one target: a target nothing designed is doomed
#' @keywords internal
#' @noRd
artefact_verdict <- function(ok, marker) {
  if (!isTRUE(ok)) {
    return("doomed")
  }
  artefact_candidate(marker)$verdict %||% NA_character_
}

#' Element-wise equality that treats NA as equal to NA
#' @keywords internal
#' @noRd
identical_na <- function(a, b) {
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
}

#' @export
print.marker_replay <- function(x, ...) {
  cat(glue::glue("<marker_replay> {substr(x$id, 1L, 16L)}: {if (x$match) 'identical' else 'DIFFERENT'}\n"))
  cat(glue::glue("  inputs: {if (x$inputs_match) 'match' else 'changed since the artefact was written'}\n"))
  engines <- if (isTRUE(x$engine_match)) {
    "match"
  } else if (isFALSE(x$engine_match)) {
    paste("changed:", paste(x$engines_changed, collapse = ", "))
  } else {
    "not recorded"
  }
  cat(glue::glue("  engines: {engines}\n"))
  print(x$table)
  invisible(x)
}
