# =============================================================================
#  marker_artefact.R  --  the self-documenting design artefact (YP3)
# =============================================================================
#  "We log everything" is a claim; an artefact is the thing that makes it
#  checkable.  `marker_artefact()` writes, beside the solved plan, the three
#  things a later reader needs to disagree with the result:
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
    run_external(git, c("-C", root, "rev-parse", "HEAD"), timeout = 10),
    error = function(e) NULL
  )
  if (is.null(out) || out$status != 0L) {
    return(NA_character_)
  }
  sha <- trimws(strsplit(out$stdout, "\n", fixed = TRUE)[[1L]])
  if (!length(sha) || !nzchar(sha[[1L]])) NA_character_ else sha[[1L]]
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
      path = as.character(fs::path_abs(path)), sha256 = hash_file(path),
      bytes = as.integer(file.size(path))
    ))
  }
  if (inherits(genome, "DNAStringSet")) {
    seqs <- as.character(genome)
    return(list(
      role = "genome", source = "DNAStringSet",
      name = paste0(length(seqs), " sequence(s)"),
      path = NA_character_,
      sha256 = hash_object(paste(seqs, collapse = "")),
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
      path = NA_character_, sha256 = hash_object(plan$variants),
      bytes = NA_integer_
    )
  }
  targets <- plan$targets
  if (!is.null(targets) && nrow(targets)) {
    out[[length(out) + 1L]] <- list(
      role = "targets", source = "table",
      name = paste0(nrow(targets), " target(s)"),
      path = NA_character_, sha256 = hash_object(targets),
      bytes = NA_integer_
    )
  }
  if (!is.null(plan$haplotypes)) {
    ref_calls <- attr(plan$haplotypes, "reference_calls")
    ref_note <- if (is.null(ref_calls)) "" else paste0(", ", ref_calls, " reference call(s)")
    out[[length(out) + 1L]] <- list(
      role = "haplotypes", source = "table",
      name = paste0(nrow(plan$haplotypes), " haplotype row(s)", ref_note),
      path = NA_character_, sha256 = hash_object(plan$haplotypes),
      bytes = NA_integer_
    )
  }
  out
}

# --- parameters and overrides ------------------------------------------------

#' The designer a scheme dispatches to, for reading its defaults
#'
#' Only used to answer "was this knob left at its default?"; a scheme this table
#' does not know reports `default = NA`, which reads as "unknown" rather
#' than as "unchanged".
#' @keywords internal
#' @noRd
artefact_designer <- function(scheme) {
  # the one resolver in R/marker_scheme_registry.R, shared with the design seam,
  # so a new built-in scheme cannot land in one switch and miss the other
  scheme_designer(scheme)
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
artefact_candidate <- function(m, cause = NA_character_) {
  out <- list(
    outcome = "rejected", verdict = NA_character_, risk = NA_character_,
    n_failed = NA_integer_, failed = NA_character_, relaxed = NA_character_,
    violated = NA_character_, stage = NA_character_, reason = NA_character_,
    discrimination = NA_real_
  )
  if (is.null(m)) {
    # ZY2: the solve already knows why the designer gave nothing back; that
    # message names the threshold, `artefact_candidate()` cannot recover it
    if (!is.na(cause)) {
      # one line: a CSV row is the contract, an error message may not be
      out$reason <- gsub("[\r\n]+", " ", cause)
      out$stage <- "design"
      return(out)
    }
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
    cause <- if ("cause" %in% names(a)) a$cause[[k]] else NA_character_
    d <- artefact_candidate(m, cause)
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

# Z7.6: a limit of the readout model, not a design bug.  When a deliberate
# mismatch is used, product starts incorporating the mismatched primer into the
# amplicon and discrimination weakens over cycles (observed in tetra-primer
# ARMS).  The dG model describes hybridisation at one temperature, not that
# accumulation, so a late cycle reads as if it were an early one.
ARTEFACT_READOUT_LIMITS <- c(
  late_cycle_discrimination = paste(
    "late-cycle discrimination: once product accumulates, the primer carrying a",
    "deliberate mismatch is incorporated into the amplicon and discrimination",
    "weakens over cycles. The dG model describes hybridisation, not this",
    "accumulation, so a late cycle is read as if it were an early one."
  )
)

# TRUE when any chosen marker carries a non-NA `stab_at` anywhere in its tree
# (the ARMS inner primers - the only designer that places a deliberate mismatch).
artefact_uses_mismatch <- function(x) {
  # a marker carries its decision tree under `$decision`, which is a tbl_graph:
  # `x[["stab_at"]]` on a graph is a vertex lookup and aborts, so the walk stops
  # at any graph rather than descending into one
  if (!is.list(x) || inherits(x, "igraph")) {
    return(FALSE)
  }
  s <- x[["stab_at"]]
  if (!is.null(s) && length(s) == 1L && !is.na(s)) {
    return(TRUE)
  }
  any(vapply(x, artefact_uses_mismatch, logical(1)))
}

# The readout limitations this artefact must declare.  Empty when no chosen
# design uses a deliberate mismatch, so the manifest never claims a mechanism
# the design does not use.
artefact_limitations <- function(solved) {
  used <- any(vapply(
    solved$report$marker, artefact_uses_mismatch,
    logical(1)
  ))
  if (!used) {
    return(character(0))
  }
  unname(ARTEFACT_READOUT_LIMITS)
}

#' The recorded verdicts in a declared format
#'
#' One row per target, all three columns character, so a reader can check a
#' record without opening its `solved.rds` (ZV13.3).
#' @keywords internal
#' @noRd
artefact_replay_table <- function(solved) {
  r <- solved$report
  # A named `scheme` would carry its names into the table, so a replay from a
  # replay.csv (unnamed) and one from solved.rds (named) would not be identical.
  tibble::tibble(
    label = unname(r$label),
    scheme = unname(r$scheme),
    verdict = unname(purrr::map2_chr(r$ok, r$marker, artefact_verdict))
  )
}

#' The `plan` block shared by the artefact and the trace
#' @keywords internal
#' @noRd
artefact_plan <- function(solved, targets) {
  list(
    layers = solved$plan$layers,
    targets = nrow(targets),
    schemes_available = solved$schemes,
    costs = as.list(solved$plan$costs %||% NA),
    ld = solved$plan$ld %||% "not_applied",
    # written by marker_design_panel(), absent for a plain solved plan
    panel = solved$plan$panel
  )
}

#' A stable id for the design, independent of when it was written
#' @keywords internal
#' @noRd
artefact_id <- function(manifest) {
  core <- manifest[c(
    "schema", "contract", "package", "engine_stamp", "inputs", "plan",
    "parameters", "overrides"
  )]
  # Hash the JSON `artefact.json` actually stores, not an in-memory variant: a
  # character NA (an absent engine's version) would otherwise hash as "NA" while
  # the file says null, and the id could never be re-derived from the record.
  # The explicit `na` is what makes the id checkable; the default is not null.
  hash_object(jsonlite::toJSON(
    core,
    auto_unbox = TRUE, null = "null", na = "null", digits = NA
  ))
}

#' The manifest timestamp: explicit value, then SOURCE_DATE_EPOCH, then now
#' @keywords internal
#' @noRd
artefact_created <- function(created = NULL) {
  if (is.null(created)) {
    created <- Sys.getenv("SOURCE_DATE_EPOCH", unset = "")
  }
  if (is.null(created) || !length(created) || (is.character(created) && !nzchar(created[[1L]]))) {
    return(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  }
  if (!length(created) || length(created) != 1L) {
    abort_msg("`created` must be NULL, one POSIXct, one epoch, or one ISO timestamp.")
  }
  if (inherits(created, "POSIXt")) {
    if (is.na(created)) {
      abort_msg("`created` must not be NA.")
    }
    return(format(created, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  }
  epoch <- if (is.numeric(created)) {
    created
  } else if (is.character(created) && grepl("^[0-9]+$", created)) {
    as.numeric(created)
  } else {
    NULL
  }
  if (!is.null(epoch)) {
    if (!is.finite(epoch) || epoch < 0) {
      abort_msg("`created` epoch must be one finite, non-negative number.")
    }
    return(format(
      as.POSIXct(epoch, origin = "1970-01-01", tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    ))
  }
  if (!is.character(created) || is.na(created)) {
    abort_msg("`created` must be NULL, a POSIXct, an epoch, or an ISO timestamp.")
  }
  at <- as.POSIXct(
    created,
    format = "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )
  if (is.na(at)) {
    abort_msg("`created` must be an ISO timestamp such as 2026-01-01T00:00:00Z.")
  }
  format(at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
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
#' reference, variant and target table, the package version and git SHA, every engine
#' version (an absent engine is `NA`, never omitted), the design knobs, and -
#' as its own field - the knobs that differ from the designer defaults.  An
#' overridden threshold is visible without diffing `parameters` against a
#' manual, because a "pass" obtained under a loosened limit must not look like
#' one obtained under the default.  `limitations` names a known limit of the
#' readout model when the chosen design uses it (today: late-cycle discrimination
#' with a deliberate mismatch), so the limit is recorded rather than discovered.
#'
#' `id` is a SHA-256 of the parts that decide the design (not the timestamp), so
#' two runs of the same plan share it and a changed input or parameter does not.
#' @param solved a `marker_solved` (see [marker_solve()]).  Its `cache` is where
#'   the rejected candidates' reasons are read from; a solve that has been
#'   round-tripped through `saveRDS()` keeps it.
#' @param dir output directory (default: a temp dir).
#' @param created artefact timestamp: `NULL` reads `SOURCE_DATE_EPOCH`, then the
#'   current time; a `POSIXct`, an epoch string/number, or an ISO timestamp is
#'   used as given.  A fixed value makes the text files byte-identical on rerun.
#' @param scenario the executable scenario this record belongs to (`"J1"`–`"J7"`,
#'   or `"unspecified"`).  It is written into `contract` and hashed into `id`;
#'   a wrong attribution is worse than an explicit gap, so callers state it
#'   rather than guess.
#' @return a `marker_artefact`: `id`, `dir`, `manifest`, `decisions`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' plan <- marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm")
#' solved <- marker_solve(plan)
#' art <- marker_artefact(solved, dir = tempfile("artefact"))
#' art
#' # re-run the design from the artefact and check the verdicts agree
#' marker_replay(art)$match
#' @family plan
#' @export
marker_artefact <- function(solved, dir = tempdir(), created = NULL,
                            scenario = "unspecified") {
  if (inherits(solved, "marker_panel")) {
    solved <- as_marker_solved(solved)
  }
  if (!inherits(solved, "marker_solved")) {
    abort_msg("`solved` must be a marker_solved (see marker_solve()).")
  }
  if (
    !is.character(scenario) || length(scenario) != 1L || is.na(scenario) ||
      !scenario %in% c(paste0("J", 1:7), "unspecified")
  ) {
    abort_msg("`scenario` must be one of J1-J7 or \"unspecified\".")
  }
  fs::dir_create(dir)
  written <- c(
    "decisions.csv", "checks.csv", "why.csv", "signals.csv", "targets.csv",
    "replay.csv", "artefact.json", "solved.rds"
  )
  refuse_links(dir, written)
  pkg_dir <- artefact_pkg_root()
  knobs <- solved$knobs %||% list()
  decisions <- artefact_decisions(solved)
  checks <- artefact_checks(solved)
  why <- if (is.null(checks)) NULL else marker_explain(checks)
  signals <- solved$signals %||% NULL

  readr::write_csv(decisions, fs::path(dir, "decisions.csv"))
  if (!is.null(checks)) {
    readr::write_csv(checks, fs::path(dir, "checks.csv"))
    readr::write_csv(tibble::as_tibble(why), fs::path(dir, "why.csv"))
  }
  if (!is.null(signals) && nrow(signals)) {
    readr::write_csv(signals, fs::path(dir, "signals.csv"))
  }
  tg <- solved$plan$targets
  readr::write_csv(
    tg[, intersect(c("label", "chrom", "pos", "ref", "alt", "type", "id"), names(tg)), drop = FALSE],
    fs::path(dir, "targets.csv")
  )
  readr::write_csv(artefact_replay_table(solved), fs::path(dir, "replay.csv"))

  manifest <- list(
    schema = "tidymarkers/artefact@2",
    trace_schema = marker_trace_schema(),
    contract = list(
      version = "tidymarkers/contract@1",
      scenario = scenario
    ),
    created = artefact_created(created),
    package = list(
      name = "tidymarkers",
      version = as.character(utils::packageVersion("tidymarkers")),
      git_sha = artefact_git_sha(),
      renv_lock_sha256 = if (is.na(pkg_dir)) {
        NA_character_
      } else {
        hash_file(fs::path(pkg_dir, "renv.lock"))
      }
    ),
    platform = list(
      r = as.character(getRversion()),
      os = R.version$platform,
      locale = Sys.getlocale("LC_COLLATE")
    ),
    engine_stamp = as.list(marker_engine_stamp()),
    inputs = artefact_inputs(solved$plan),
    plan = artefact_plan(solved, tg),
    parameters = purrr::map(knobs, artefact_value),
    overrides = artefact_overrides(knobs, unique(solved$attempts$scheme)),
    # Z7.6: a limit of the readout model that is not a design bug.  Always a
    # list in the JSON (`I()` keeps a one-element vector an array), empty when
    # no chosen design uses a deliberate mismatch.
    limitations = I(as.character(artefact_limitations(solved))),
    counts = artefact_counts(solved, decisions, signals),
    files = list(
      decisions = "decisions.csv",
      checks = if (is.null(checks)) NA_character_ else "checks.csv",
      why = if (is.null(why)) NA_character_ else "why.csv",
      signals = if (is.null(signals) || !nrow(signals)) NA_character_ else "signals.csv",
      targets = "targets.csv",
      # the recorded verdicts in a declared format: a record from a stranger can
      # be checked by passing a trusted solve, without opening its RDS (ZV13.3)
      replay = "replay.csv",
      # excluded by name and reason: an igraph pointer and a random id make
      # the RDS byte-unstable, so it is not part of the reproducible set
      solved = list(name = "solved.rds", excluded = "igraph pointer and random id")
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

# A sidecar is named by the manifest and may be absent (`NA` in `files`); a
# missing file is an empty table, never an error -- a solve that produced no row
# still wrote a record.
artefact_safe_name <- function(file) {
  # the manifest is a file a stranger can hand us: never resolve a path out of
  # the artefact directory
  ok <- is.character(file) && length(file) == 1L && !is.na(file) && nzchar(file) &&
    !grepl("[/\\\\]", file) && !file %in% c(".", "..")
  if (!ok) {
    abort_msg("The manifest names an unsafe file path: {.val {file}}.")
  }
  file
}

# A sidecar name inside the record can still be a link that leaves it: the
# lexical check above is not enough for a record a stranger wrote.
artefact_member_path <- function(dir, file) {
  path <- fs::path(dir, file)
  if (fs::file_exists(path)) {
    target <- fs::path_real(path)
    root <- fs::path_real(dir)
    if (!identical(target, root) && !fs::path_has_parent(target, root)) {
      abort_msg(
        "The record names {.file {file}}, which resolves to {.path {target}} outside {.path {dir}}.",
        i = "A record carries files, not links out of itself."
      )
    }
  }
  path
}

#' Refuse to write a record member through an existing link
#'
#' The writer creates each member.  A link at the name would send the write to
#' its target, and a dangling link would create the target outside the record.
#' @keywords internal
#' @noRd
refuse_links <- function(dir, names) {
  linked <- names[fs::is_link(fs::path(dir, names))]
  if (length(linked)) {
    abort_msg(
      c(
        "The record writes {.file {linked[[1L]]}}, which is a link.",
        i = "A record creates its own files; remove the link first."
      )
    )
  }
  invisible(NULL)
}

artefact_sidecar <- function(dir, file, empty) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    return(empty())
  }
  path <- artefact_member_path(dir, artefact_safe_name(file))
  if (!fs::file_exists(path)) {
    return(empty())
  }
  readr::read_csv(path, show_col_types = FALSE)
}

# The recorded verdicts, when the manifest names them.  A record written before
# the sidecar existed returns NULL and the caller falls back to the RDS.
artefact_recorded_replay <- function(dir, manifest) {
  file <- manifest$files$replay
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    return(NULL)
  }
  path <- artefact_member_path(dir, artefact_safe_name(file))
  if (!fs::file_exists(path)) {
    return(NULL)
  }
  tab <- readr::read_csv(
    path,
    col_types = readr::cols(label = "c", scheme = "c", verdict = "c"),
    show_col_types = FALSE
  )
  if (!all(c("label", "scheme", "verdict") %in% names(tab))) {
    abort_msg(
      "The artefact's replay table {.val {file}} must carry `label`, `scheme` and `verdict`."
    )
  }
  if (anyNA(tab$label) || !all(nzchar(tab$label))) {
    abort_msg("The artefact's replay table {.val {file}} names no target on a row.")
  }
  dup <- unique(tab$label[duplicated(tab$label)])
  if (length(dup)) {
    # two verdicts for one target: the record does not say which one is its claim
    abort_msg(
      "The artefact's replay table {.val {file}} names a target twice: {.val {dup}}."
    )
  }
  tab[, c("label", "scheme", "verdict")]
}

#' Read a persisted artefact back
#'
#' [marker_artefact()] writes its decision rows to `decisions.csv` beside the
#' manifest, so a bare `artefact.json` cannot carry them.  Reading the directory
#' back is what lets [marker_trace_from_artefact()] migrate a persisted
#' `artefact@1` record without losing its decisions and signals.
#'
#' A record from a newer tidymarkers is refused and names its version: the reader
#' is too old, not the record wrong.
#' @param path the directory [marker_artefact()] wrote to, or the `artefact.json`
#'   inside it.
#' @return a `marker_artefact`: `id`, `dir`, `manifest`, `decisions`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' plan <- marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm")
#' dir <- tempfile("artefact")
#' marker_artefact(marker_solve(plan), dir = dir)
#' str(marker_artefact_read(dir), max.level = 1)
#' @family plan
#' @export
marker_artefact_read <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    abort_msg("`path` must be one directory or artefact.json path.")
  }
  dir <- if (fs::is_dir(path)) path else fs::path_dir(path)
  manifest <- artefact_manifest(dir)
  parts <- trace_schema_parts(manifest$schema)
  if (is.null(parts) || parts$kind != "artefact") {
    abort_msg(
      "Not an artefact: `schema` is {.val {trace_schema_label(manifest$schema)}}.",
      i = "Expected a {.val {'tidymarkers/artefact'}} record."
    )
  }
  current <- TRACE_FORMATS[["artefact"]]
  if (parts$version > current) {
    abort_msg(
      "This artefact is {.val {manifest$schema}}: newer than this reader's artefact schema.",
      i = "Update tidymarkers to read it."
    )
  }
  if (parts$version < 1L) {
    trace_abort_unknown(manifest$schema)
  }
  files <- manifest$files %||% list()
  structure(
    list(
      id = manifest$id,
      dir = as.character(dir),
      manifest = manifest,
      decisions = artefact_sidecar(dir, files$decisions, trace_empty_decisions)
    ),
    class = "marker_artefact"
  )
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
  abort_msg("`artefact` must be a marker_artefact or a directory.")
}

# The trace and the manifest publish the same counts; one definition keeps
# them from drifting apart.
artefact_counts <- function(solved, decisions, signals) {
  list(
    targets = nrow(solved$report),
    solved = sum(solved$report$ok),
    attempts = nrow(decisions),
    rejected = sum(decisions$outcome == "rejected"),
    signals = if (is.null(signals)) 0L else nrow(signals)
  )
}

#' Read an artefact's manifest back
#' @keywords internal
#' @noRd
artefact_manifest <- function(dir) {
  path <- artefact_member_path(dir, "artefact.json")
  if (!fs::file_exists(path)) {
    abort_msg("No artefact.json in {.path {dir}}; write one with marker_artefact().")
  }
  manifest <- jsonlite::fromJSON(readr::read_file(path), simplifyVector = FALSE)
  artefact_verify_id(manifest, path)
  manifest
}

# The id is the hash of the manifest's own core, so a reader can check it. An
# artefact@1 record predates `contract`, and its core cannot be re-derived;
# every artefact@2 was written by this package, so the check costs nothing there.
# A newer schema is refused by the reader for that reason alone, not for its id.
artefact_verify_id <- function(manifest, path) {
  parts <- trace_schema_parts(manifest$schema)
  verifiable <- !is.null(parts) && parts$kind == "artefact" &&
    parts$version >= 2L && parts$version <= TRACE_FORMATS[["artefact"]]
  if (!verifiable) {
    return(invisible(NULL))
  }
  expected <- artefact_id(manifest)
  stored <- manifest$id
  if (!is.character(stored) || length(stored) != 1L || is.na(stored)) {
    abort_msg(
      "The record in {.path {path}} names no id to check.",
      i = "A manifest truncated after it was written hashes to {.val {expected}}."
    )
  }
  if (!identical(stored, expected)) {
    abort_msg(
      "The id in {.path {path}} does not match its own manifest: {.val {stored}}.",
      i = "The manifest hashes to {.val {expected}}; it was edited after it was written."
    )
  }
  invisible(NULL)
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
#' Re-solves the recorded plan with the recorded knobs and compares the chosen
#' chemistry and verdict per target.  This is what turns "we logged everything"
#' into a claim that can be checked: a mismatch names the target and both
#' verdicts instead of silently returning a different design.
#'
#' The input hashes and the engine stamp are re-checked first, so a reference that
#' was replaced is reported as such (`inputs_match = FALSE`) and a library upgrade
#' as an engine change (`engine_match = FALSE`, `engines_changed`), rather than
#' either being read as a model change.
#'
#' `solved.rds` is a serialized R object, and `readRDS()` executes what the file
#' tells it to. Replay is therefore for an artefact this machine wrote, or one
#' from a source whose code you would run anyway. A record written by
#' [marker_artefact()] also names its verdicts in `replay.csv`; pass a trusted
#' `solved` and the replay compares against that file without opening the RDS.
#' [marker_import_rocrate()] never opens the RDS, so a crate that arrived from a
#' stranger can be inspected without being trusted with execution.
#' @param artefact a `marker_artefact` or the directory it was written to.
#' @param cache an environment to memoise the replay in; pass the original
#'   solve's `cache` to make replay free, `NULL` (the default) to force the
#'   design to run again - which is the only version that tests anything.
#' @param progress show a progress bar (default: not, the replay is a check).
#' @param tolerance accepted; reserved for numeric comparison of the markers.
#' @param solved a `marker_solved` you trust (see [marker_solve()]).  When given,
#'   its plan is re-solved and compared with the recorded `replay.csv`, and
#'   `solved.rds` is never read.  This is the way to replay a record from an
#'   untrusted source.  A record that names no `replay.csv` has no verdicts of
#'   its own to compare, so it is refused instead of passing on your own solve.
#' @return a `marker_replay`: `id`, `match` (one logical for the whole run),
#'   `inputs_match`, `engine_match` (`NA` when the artefact recorded no stamp),
#'   `engines_changed` (names whose version differs), and `table` (one row per
#'   target: `target`, `recorded_scheme`, `replayed_scheme`,
#'   `recorded_verdict`, `replayed_verdict`, `match`).
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' plan <- marker_plan(genome) + marker_targets(loci) + marker_tech(schemes = "hrm")
#' art <- marker_artefact(marker_solve(plan), dir = tempfile("artefact"))
#' marker_replay(art)$match
#' @family plan
#' @export
marker_replay <- function(artefact, cache = NULL, progress = FALSE, tolerance = NULL,
                          solved = NULL) {
  dir <- artefact_dir(artefact)
  manifest <- artefact_manifest(dir)
  recorded <- artefact_recorded_replay(dir, manifest)
  trusted <- !is.null(solved)
  if (!trusted) {
    path <- artefact_member_path(dir, "solved.rds")
    if (!fs::file_exists(path)) {
      abort_msg(
        "No solved.rds in {.path {dir}}; the artefact cannot be replayed.",
        i = if (is.null(recorded)) {
          NULL
        } else {
          "Pass `solved =` a marker_solved you trust to replay from `replay.csv` without opening an RDS."
        }
      )
    }
    solved <- readRDS(path)
  } else if (!inherits(solved, "marker_solved")) {
    abort_msg("`solved` must be a marker_solved (see marker_solve()).")
  }
  if (is.null(recorded)) {
    if (trusted) {
      # Without a replay.csv the record's verdicts live only in its RDS, and a
      # trusted `solved` replayed against itself always matches (ADR 0024).
      abort_msg(
        "The artefact records no verdicts of its own, so a trusted `solved` would be compared with itself.",
        i = "Replay it without `solved` to read its own `solved.rds`, or use a record that names `replay.csv`."
      )
    }
    recorded <- artefact_replay_table(solved)
  }
  inputs_match <- artefact_inputs_match(manifest$inputs, solved$plan)
  if (isFALSE(inputs_match)) {
    warn_msg("The reference or variant table no longer matches the artefact's hash; \\
       a replay difference below is explained by the input, not the engine.")
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
  rep_ <- again$report
  # A row missing from either side is a difference, not a silence: a record that
  # dropped a verdict must not pass by never mentioning that target.
  labels <- unique(c(recorded$label, rep_$label))
  rec <- recorded[match(labels, recorded$label), ]
  rep2 <- rep_[match(labels, rep_$label), ]
  replayed_verdict <- unname(purrr::map2_chr(rep2$ok, rep2$marker, artefact_verdict))
  tbl <- tibble::tibble(
    target = labels,
    recorded_scheme = rec$scheme, replayed_scheme = rep2$scheme,
    recorded_verdict = rec$verdict, replayed_verdict = replayed_verdict
  )
  tbl$match <- !is.na(match(labels, recorded$label)) &
    !is.na(match(labels, rep_$label)) &
    identical_na(tbl$recorded_scheme, tbl$replayed_scheme) &
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
