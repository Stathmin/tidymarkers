# =============================================================================
#  marker_solve.R  --  the solver: run a plan, cache, report
# =============================================================================
#  `marker_solve()` is the compute verb of the plan grammar (like
#  `ggplot2::ggplot_build()`).  It walks each target's allowed chemistries
#  cheapest-first (the tech layer's cost schema), stops at the first that designs,
#  keeps a memo cache so re-solving is free, and returns a `marker_solved` - a tidy
#  report plus the per-scheme attempt matrix - ready to plot, check for tube
#  compatibility, and export.
# =============================================================================

#' The variants table a target needs (a one-row tibble)
#' @keywords internal
#' @noRd
target_snp <- function(row) {
  out <- tibble::tibble(
    chrom = as.character(row$chrom[[1L]]), pos = as.integer(row$pos[[1L]]),
    ref = as.character(row$ref[[1L]]), alt = as.character(row$alt[[1L]])
  )
  if (!is.null(row$alt2) && !is.na(row$alt2[[1L]])) {
    out$alt2 <- as.character(row$alt2[[1L]])
  }
  if ("strand" %in% names(row)) {
    out$strand <- as.character(row$strand[[1L]])
  }
  out
}

#' A design-cache key for one target
#' @keywords internal
#' @noRd
target_key <- function(tg) paste(tg$chrom, tg$pos, tg$ref, tg$alt, sep = "|")

#' A cache key that also varies with the design knobs
#' @keywords internal
#' @noRd
knobs_key <- function(dots) {
  if (!length(dots)) {
    return("")
  }
  paste(names(dots), vapply(dots, function(x) paste(x, collapse = ","), character(1)), sep = "=", collapse = ";")
}

#' Is a designer's return a usable marker?
#'
#' Class-based (`is_target_designable()`) plus a scheme's own `verdict`: a
#' `marker_lax`/`marker_duplex` object can carry a `"fail"` verdict (its inner assay
#' was doomed), so the class alone is not enough.
#' @keywords internal
#' @noRd
marker_ok <- function(m) {
  if (!is_target_designable(m)) {
    return(FALSE)
  }
  if (!is.null(m$verdict)) {
    return(!identical(m$verdict, "fail"))
  }
  TRUE
}

#' Why a scheme attempt did not design its target
#'
#' `solve_one()` keeps only the message, not the condition: the attempt table is
#' written to disk, and a condition object carries an environment that must not
#' travel with the report.
#' @keywords internal
#' @noRd
attempt_cause <- function(m, ok) {
  if (ok) {
    return(NA_character_)
  }
  if (is.null(m)) {
    return(gettext("the designer returned no marker object (it may have errored)"))
  }
  if (inherits(m, "marker_doomed")) {
    return(paste(m$stage, m$reason, sep = ": "))
  }
  if (!is.null(m$verdict) && identical(m$verdict, "fail")) {
    return(gettext("the designed marker carries a fail verdict"))
  }
  gettext("no designable marker was returned")
}

#' The cost of a (possibly unchosen) scheme
#' @keywords internal
#' @noRd
scheme_cost <- function(scheme, costs) {
  if (length(scheme) != 1L || is.na(scheme)) {
    return(NA_real_)
  }
  v <- unname(costs[scheme])
  if (is.na(v)) Inf else v
}

#' Solve one target cheapest-first, recording every attempt
#' @keywords internal
#' @noRd
solve_one <- function(plan, row, schemes, dots, cache, all = FALSE) {
  snp <- target_snp(row)
  kk <- paste(target_key(row), sep = "")
  costs <- plan$costs %||% KASP_TECH_COSTS
  attempts <- list()
  chosen <- NA_character_
  marker <- NULL
  for (sc in schemes) {
    ck <- paste(kk, sc, knobs_key(dots), sep = "||")
    # ZY2: a design that errored is cached as NULL, so its message is kept beside
    # it - otherwise a re-solve of the same target loses the one thing that says
    # which threshold rejected the scheme
    killed <- paste0(ck, "||cause")
    known <- exists(ck, envir = cache, inherits = FALSE)
    cause <- if (exists(killed, envir = cache, inherits = FALSE)) {
      get(killed, envir = cache, inherits = FALSE)
    } else {
      NULL
    }
    m <- if (known) get(ck, envir = cache, inherits = FALSE) else NULL
    if (!known) {
      args <- c(list(plan$genome, snp, sc), dots)
      if (!is.null(plan$variants)) args$variants <- plan$variants
      if (!is.null(plan$haplotypes) && scheme_base(sc) %in% haplotype_schemes()) args$haplotypes <- plan$haplotypes
      # the error becomes a value, so `attempts` records it and the walk still
      # reaches the next scheme; a warning is muffled, not caught, so it cannot
      # change control flow
      m <- withCallingHandlers(
        tryCatch(rlang::exec(marker_design_at, !!!args), error = function(e) e),
        warning = function(w) invokeRestart("muffleWarning")
      )
      if (inherits(m, "error")) {
        cause <- conditionMessage(m)
        m <- NULL
      }
      assign(ck, m, envir = cache)
      if (!is.null(cause)) assign(killed, cause, envir = cache)
    }
    ok <- marker_ok(m)
    if (is.null(cause)) cause <- attempt_cause(m, ok)
    attempts[[length(attempts) + 1L]] <- tibble::tibble(
      scheme = sc, cost = scheme_cost(sc, costs), ok = ok, cause = cause
    )
    if (ok && is.na(chosen)) {
      chosen <- sc
      marker <- m
      if (!all) break
    }
  }
  at <- purrr::list_rbind(attempts)
  list(
    scheme = chosen, cost = scheme_cost(chosen, costs), ok = !is.na(chosen),
    marker = marker, tried = paste(at$scheme, collapse = ", "), attempts = at
  )
}

#' Solve a marker-design plan
#'
#' Runs the plan's targets against the plan's allowed chemistries, cheapest-first
#' (the tech cost schema), stopping at the first that designs each target.  The
#' expensive per-(target, chemistry) design is memoised in `cache`, so a repeated
#' target - or a re-solve with the same cache - costs nothing.  The result is a
#' tidy report (`report`), the per-scheme attempt matrix (`attempts`), and the
#' designed marker objects (`report$marker`), which [marker_tubes()], [marker_plot()]
#' and [marker_export()] consume.
#' @param plan a `marker_plan` (see [marker_plan()]).
#' @param ... passed to the per-scheme designer (e.g. `flank`, `cores`,
#'   `run_blast`).
#' @param cache an environment to memoise designs in; pass the same one across
#'   solves to reuse work.  `NULL` uses a fresh cache.
#' @param all evaluate every allowed scheme rather than stopping at the first
#'   that designs (a full attempt matrix; slower).
#' @param progress show a progress bar (default: when interactive).
#' @return a `marker_solved`: `report` (one row per target: `scheme`, `cost`, `ok`,
#'   `tried`, `marker`), `attempts` (one row per target x tried scheme, with
#'   `cause` - why that scheme did not design the target, `NA` when it did), the
#'   `plan`/`schemes`/`cache`, `knobs` - the `...` the design ran with, kept
#'   so [marker_artefact()] can record them and [marker_replay()] can re-run the
#'   same design instead of a similar one - and, for a plan with a `region`, a
#'   `signals` table naming each haplotype's readout (`within_phased_block`).
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' # `+` binds tighter than `|>`, so solve the whole plan in one call
#' solved <- marker_solve(
#'   marker_plan(genome) + marker_targets(loci) + marker_tech(preset = "lab")
#' )
#' solved
#' @family plan
#' @export
marker_solve <- function(plan, ..., cache = NULL, all = FALSE, progress = interactive()) {
  if (!inherits(plan, "marker_plan")) {
    abort_msg("`plan` must be a marker_plan (see marker_plan()).")
  }
  tg <- plan$targets
  if (!is.null(plan$region)) {
    tg <- dplyr::bind_rows(tg, region_targets(plan))
  }
  tg <- filter_targets_by_depth(tg, plan$depth)
  if (!nrow(tg)) {
    abort_msg(
      "The plan has no targets; add them with marker_targets() or name a `region`."
    )
  }
  schemes <- plan_schemes(plan)
  if (is.null(cache)) cache <- new.env(parent = emptyenv())
  dots <- rlang::list2(...)

  key <- target_key(tg)
  idx <- which(!duplicated(key))
  nb <- if (isTRUE(progress)) cli::cli_progress_bar("Designing", total = length(idx))
  sols <- purrr::map(idx, function(i) {
    if (isTRUE(progress)) cli::cli_progress_update(id = nb)
    solve_one(plan, tg[i, , drop = FALSE], schemes, dots, cache, all = all)
  })
  if (isTRUE(progress)) cli::cli_progress_done(id = nb)
  names(sols) <- key[idx]
  recs <- sols[key]

  report <- tibble::tibble(
    label = tg$label, chrom = tg$chrom, pos = tg$pos, ref = tg$ref, alt = tg$alt, type = tg$type,
    scheme = purrr::map_chr(recs, "scheme"),
    cost = purrr::map_dbl(recs, "cost"),
    ok = purrr::map_lgl(recs, "ok"),
    tried = purrr::map_chr(recs, "tried"),
    marker = purrr::map(recs, "marker")
  )
  # a region's answer is only a statement about haplotypes inside its phased
  # block; the table says so on every row (roadmap Y2.1)
  signals <- if (is.null(plan$region)) NULL else region_signals(plan, report)
  attempts <- purrr::map2(
    tg$label, recs,
    function(lb, r) {
      a <- r$attempts
      a$label <- lb
      a[, c("label", "scheme", "cost", "ok", "cause")]
    }
  ) |>
    purrr::list_rbind()

  structure(
    list(
      report = report, attempts = attempts, plan = plan, schemes = schemes,
      cache = cache, knobs = dots, signals = signals
    ),
    class = "marker_solved"
  )
}

#' @export
print.marker_solved <- function(x, ...) {
  r <- x$report
  cat(glue::glue("<marker_solved> {nrow(r)} target(s), {sum(r$ok)} solved\n"))
  print(r[setdiff(names(r), "marker")])
  invisible(x)
}

#' Plot a solved plan: which chemistry takes each target
#'
#' The attempt matrix as tiles: a row per target, a column per tried scheme,
#' shaded by cost, labelled with the outcome - the cheapest chemistry that
#' designs each target is the leftmost "solved" tile.
#' @param solved a `marker_solved` (see [marker_solve()]).
#' @param ... reserved for future options.
#' @return a `ggplot`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' solved <- marker_solve(
#'   marker_plan(genome) + marker_targets(loci) + marker_tech(preset = "lab")
#' )
#' marker_plot(solved)
#' @family plan
#' @export
marker_plot <- function(solved, ...) {
  if (!inherits(solved, "marker_solved")) {
    abort_msg("`solved` must be a marker_solved (see marker_solve()).")
  }
  rlang::check_installed("ggplot2", reason = gettext("to plot a solved plan"))
  a <- solved$attempts
  a$outcome <- dplyr::if_else(a$ok, "solved", "tried")
  ggplot2::ggplot(a, ggplot2::aes(x = .data$scheme, y = .data$label)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$cost), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$outcome), size = 2.6, colour = "grey20") +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey93") +
    ggplot2::labs(
      x = "chemistry (cheapest first)", y = NULL, fill = "cost",
      title = "Which chemistry takes each target"
    ) +
    ggplot2::theme_minimal()
}

#' Write a solved plan to disk
#'
#' Writes the target report, the attempt matrix, every solved marker's primers,
#' the whole object, and - through [marker_artefact()] - the manifest
#' (`artefact.json`), the per-candidate decision log (`decisions.csv`), the gate
#' table of the chosen markers (`checks.csv`) and the target table.
#'
#' The artefact is the part a later reader can disagree with: it names the
#' hashed inputs, the engine versions and the overridden knobs, and
#' [marker_replay()] re-runs the design from `solved.rds` and compares.  Pass
#' `artefact = FALSE` for the bare tables.
#' @param solved a `marker_solved` (see [marker_solve()]).
#' @param dir output directory (default: a temp dir).
#' @param artefact also write the self-documenting manifest and decision log
#'   (default `TRUE`).
#' @return `dir`, invisibly.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' solved <- marker_solve(
#'   marker_plan(genome) + marker_targets(loci) + marker_tech(preset = "lab")
#' )
#' marker_export(solved)
#' @family plan
#' @export
marker_export <- function(solved, dir = tempdir(), artefact = TRUE) {
  if (inherits(solved, "marker_panel")) {
    solved <- as_marker_solved(solved)
  }
  if (!inherits(solved, "marker_solved")) {
    abort_msg("`solved` must be a marker_solved (see marker_solve()).")
  }
  fs::dir_create(dir)
  # the same rule as marker_artefact() (ADR 0025 p.5): create the files, never
  # write through a link that would put them outside the directory
  refuse_links(dir, c("markers.csv", "attempts.csv", "primers.csv", "solved.rds"))
  readr::write_csv(
    solved$report[setdiff(names(solved$report), "marker")], fs::path(dir, "markers.csv")
  )
  readr::write_csv(solved$attempts, fs::path(dir, "attempts.csv"))
  readr::write_csv(marker_primers_table(solved), fs::path(dir, "primers.csv"))
  saveRDS(solved, fs::path(dir, "solved.rds"))
  if (isTRUE(artefact)) {
    marker_artefact(solved, dir = dir)
  }
  cli::cli_alert_success(gettext(
    "Wrote {.path {dir}} ({paste(sort(list.files(dir)), collapse = ', ')})."
  ))
  invisible(dir)
}

#' The primers of every solved marker, as a tidy table
#' @keywords internal
#' @noRd
marker_primers_table <- function(solved) {
  r <- solved$report
  purrr::map(which(r$ok), function(i) {
    pm <- marker_primers(r$marker[[i]])
    if (!length(pm)) {
      return(NULL)
    }
    tibble::tibble(
      target = r$label[i], scheme = r$scheme[i],
      role = names(pm), sequence = unname(pm)
    )
  }) |>
    purrr::list_rbind()
}
