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
    if (exists(ck, envir = cache, inherits = FALSE)) {
      m <- get(ck, envir = cache, inherits = FALSE)
    } else {
      args <- c(list(plan$genome, snp, sc), dots)
      if (!is.null(plan$variants)) args$variants <- plan$variants
      if (!is.null(plan$haplotypes) && sc %in% haplotype_schemes()) args$haplotypes <- plan$haplotypes
      m <- tryCatch(
        suppressWarnings(rlang::exec(marker_design_at, !!!args)),
        error = function(e) NULL
      )
      assign(ck, m, envir = cache)
    }
    ok <- marker_ok(m)
    attempts[[length(attempts) + 1L]] <- tibble::tibble(
      scheme = sc, cost = scheme_cost(sc, costs), ok = ok
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
#'   `tried`, `marker`), `attempts` (one row per target x tried scheme), the
#'   `plan`/`schemes`/`cache`, `knobs` - the `...` the design ran with, kept
#'   so [marker_artefact()] can record them and [marker_replay()] can re-run the
#'   same design instead of a similar one - and, for a plan with a `region`, a
#'   `signals` table naming each haplotype's readout (`within_phased_block`).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' solved <- marker_plan(genome) +
#'   marker_targets(marker_example("snp")) +
#'   marker_tech(preset = "lab") |>
#'   marker_solve()
#' solved
#' }
#' @family KASP plan
#' @export
marker_solve <- function(plan, ..., cache = NULL, all = FALSE, progress = interactive()) {
  if (!inherits(plan, "marker_plan")) {
    cli::cli_abort(gettext("`plan` must be a marker_plan (see marker_plan())."))
  }
  tg <- plan$targets
  if (!is.null(plan$region)) {
    tg <- dplyr::bind_rows(tg, region_targets(plan))
  }
  tg <- filter_targets_by_depth(tg, plan$depth)
  if (!nrow(tg)) {
    cli::cli_abort(gettext(
      "The plan has no targets; add them with marker_targets() or name a `region`."
    ))
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
      a[, c("label", "scheme", "cost", "ok")]
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
#' \dontrun{
#' marker_plot(marker_solve(plan))
#' }
#' @family KASP plan
#' @export
marker_plot <- function(solved, ...) {
  if (!inherits(solved, "marker_solved")) {
    cli::cli_abort(gettext("`solved` must be a marker_solved (see marker_solve())."))
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
#' \dontrun{
#' marker_export(marker_solve(plan), dir = "results")
#' }
#' @family KASP plan
#' @export
marker_export <- function(solved, dir = tempdir(), artefact = TRUE) {
  if (!inherits(solved, "marker_solved")) {
    cli::cli_abort(gettext("`solved` must be a marker_solved (see marker_solve())."))
  }
  fs::dir_create(dir)
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
