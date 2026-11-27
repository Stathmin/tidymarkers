# =============================================================================
#  marker_batch.R  --  many loci x many traits in one call
# =============================================================================
#  `design_markers()` already walks a table of loci, but it returns a summary:
#  no assay objects, no timing, and `suppressWarnings()` + `tryCatch(NULL)` collapse
#  the warnings and the error message, so a doomed locus cannot say why.  This
#  runner returns the design itself per locus, with status, timing and message.

#' Run a per-locus design (or any single-locus call) over many loci, in parallel
#'
#' `design_markers()` is the summary of a design table; this is the batch that
#' keeps the result.  For a data frame of loci x traits it runs one design per
#' locus - [design_auto_marker()] on the rescue ladder, or one scheme, or an
#' arbitrary function of the locus row - and returns one row per locus with the
#' status (`"ok"` / `"fail"`), the message, the elapsed time and the object
#' itself.  A locus that errors is captured (`status = "fail"`, the error message
#' in `message`); it never aborts the batch.
#'
#' The results are identical to looping the same call by hand, and the loop is
#' literally the serial path: `BPPARAM = NULL` runs in this session.  Passing a
#' `MulticoreParam(k)` or `SnowParam(k)` runs the loci on `k` workers; the
#' inner call is forced to one core so the machine is not oversubscribed.
#'
#' The per-locus work runs with `cores = 1`: an outer pool that also forked inside
#' each worker would oversubscribe the machine (`k^2` processes), the same rule
#' [design_kasp_markers()] follows.
#'
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param loci data frame/tibble of loci; needs `chrom`, `pos`, `ref`, `alt` when
#'   `design` is a scheme or `NULL`.
#' @param design a scheme id (e.g. `"snp"`), `NULL` for the rescue ladder
#'   ([design_auto_marker()]), or a function of one locus row returning the object
#'   to keep.
#' @param BPPARAM optional `BiocParallel` backend; `NULL` (the default) runs
#'   serially in this session.
#' @param ... passed to the per-locus call (e.g. `flank`).
#' @return a tibble, one row per input locus: the `loci` columns, `status`,
#'   `scheme`, `message`, `elapsed` (s) and `result` (list-column of the
#'   returned object, `NULL` when the locus failed).
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' batch <- design_marker_batch(genome, loci)
#' batch[, c("chrom", "pos", "status", "scheme")]
#' @seealso [design_kasp_markers()], [design_markers()], [design_auto_marker()]
#' @family design
#' @export
design_marker_batch <- function(genome, loci, design = NULL, BPPARAM = NULL, ...) {
  if (is.null(design) || is.character(design)) {
    need_cols(loci, c("chrom", "pos", "ref", "alt"), arg = "loci")
  }
  dots <- rlang::list2(...)
  one <- function(i) {
    row <- loci[i, , drop = FALSE]
    t0 <- proc.time()[[3]]
    warned <- character(0)
    # The outcome is *returned*, never assigned into an enclosing list: under
    # `pkgload::load_all()` an assignment made inside a `tryCatch` handler does
    # not reach the enclosing frame, so a box-based version silently lost the
    # error.  Warnings are collected here and turned into the row message.
    outcome <- withCallingHandlers(
      tryCatch(
        {
          res <- marker_batch_design_one(genome, row, design, dots)
          scheme <- if (is.character(design) && length(design) == 1L) {
            design
          } else if (is.null(design)) {
            res$scheme
          } else {
            NA_character_
          }
          # the ladder treats a doomed locus as a valid (negative) answer; a named
          # scheme is only ok when it actually produced a designable marker
          ok <- if (is.null(design) || is.function(design)) {
            TRUE
          } else {
            is_target_designable(res)
          }
          list(
            status = if (ok) "ok" else "fail", scheme = as.character(scheme),
            message = if (ok) NA_character_ else "not designable for this scheme",
            result = res
          )
        },
        error = function(e) {
          list(
            status = "fail",
            scheme = if (is.character(design) && length(design) == 1L) design else NA_character_,
            message = conditionMessage(e), result = NULL
          )
        }
      ),
      warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    if (length(warned)) {
      seen <- if (is.na(outcome$message)) character(0) else outcome$message
      outcome$message <- paste(c(seen, unique(warned)), collapse = "; ")
    }
    tibble::tibble(
      status = outcome$status, scheme = outcome$scheme, message = outcome$message,
      elapsed = proc.time()[[3]] - t0, result = list(outcome$result)
    )
  }

  # `serial_backend()` says `MulticoreParam` is serial (the workers do not get a
  # BiocParallel pool); that is right for the *inner* call and wrong here, where
  # the backend is the outer pool.  Only a backend with > 1 worker parallelises.
  if (is.null(BPPARAM) || marker_ncores(BPPARAM) == 1L) {
    rows <- purrr::map(seq_len(nrow(loci)), one)
  } else {
    # one parent build, then every worker only opens an existing shard
    index_prebuild(genome)
    rows <- bp_lapply(seq_len(nrow(loci)), function(i) {
      worker_biostrings()
      one(i)
    }, BPPARAM)
  }
  # `cbind()` loses the tibble class (tibble has no `cbind.tbl_df`) and cannot
  # hold a list-column at all, so the design columns are assigned in place
  out <- tibble::as_tibble(loci)
  out$status <- purrr::map_chr(rows, "status")
  out$scheme <- purrr::map_chr(rows, "scheme")
  out$message <- purrr::map_chr(rows, "message")
  out$elapsed <- vapply(rows, function(r) r$elapsed, numeric(1))
  out$result <- lapply(rows, function(r) r$result[[1L]])
  out
}

#' One locus through the batch runner's dispatch, with the per-locus cores forced
#' @param genome,loc,design,dots see `design_marker_batch()`.
#' @return the designed object (or the `design` function's return value).
#' @keywords internal
#' @noRd
marker_batch_design_one <- function(genome, loc, design, dots) {
  if (is.function(design)) {
    return(rlang::exec(design, loc))
  }
  # the built-in designers get the same canonical row `design_markers()` builds
  # (`alt2`/`strand` filled in); a caller's own `design` function above still gets
  # its row untouched, because its input shape is the caller's contract
  snp <- locus_row(loc)
  dots$cores <- 1L
  if (is.null(design)) {
    return(rlang::exec(design_auto_marker, genome, snp, !!!dots))
  }
  rlang::exec(marker_design_at, genome, snp, design, !!!dots)
}

#' Continue a batch report as a solved object
#'
#' J2 ends with a table of per-locus verdicts.  The scenarios continue that run
#' through the typed pipeline, so the report becomes a `marker_solved`: the same
#' rows as `report`, one attempt per locus in `attempts`, and the assays the
#' batch built in `report$marker`.  `plan` is `NULL`, because a batch design was
#' not driven by a plan; `marker_artefact()` therefore needs a plan and the bare
#' tables are exported with `marker_export(artefact = FALSE)`.
#'
#' @param x a `marker_batch` report from [design_kasp_markers()], a
#'   `marker_solved` (returned unchanged), or a `marker_panel` (its `solved`).
#' @return a `marker_solved`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' solved <- as_marker_solved(design_kasp_markers(genome, loci[1, ], flank = 150L))
#' solved
#' @family design
#' @export
as_marker_solved <- function(x) {
  UseMethod("as_marker_solved")
}

#' @export
as_marker_solved.marker_solved <- function(x) x

#' @export
as_marker_solved.marker_panel <- function(x) {
  out <- x$solved
  attr(out, "panel") <- x$artefact$manifest$plan$panel
  out
}

#' @export
as_marker_solved.marker_batch <- function(x) {
  markers <- attr(x, "markers")
  if (is.null(markers) || length(markers) != nrow(x)) {
    abort_msg("This batch report carries no assays; rebuild it with design_kasp_markers().")
  }
  ok <- x$verdict %in% c("pass", "warn") & !vapply(markers, is.null, logical(1))
  report <- tibble::tibble(
    label = paste0(x$chrom, ":", x$pos, ":", x$ref, ">", x$alt),
    chrom = x$chrom, pos = x$pos, ref = x$ref, alt = x$alt, type = x$type,
    scheme = ifelse(ok, "snp", NA_character_),
    cost = ifelse(ok, unname(KASP_TECH_COSTS[["snp"]]), NA_real_),
    ok = ok,
    tried = "snp",
    marker = markers
  )
  attempts <- tibble::tibble(
    label = report$label, scheme = report$scheme, cost = report$cost,
    ok = report$ok,
    cause = ifelse(report$ok, NA_character_, x$note)
  )
  structure(
    list(
      report = report, attempts = attempts, plan = NULL, schemes = "snp",
      cache = NULL, knobs = list(), signals = NULL
    ),
    class = "marker_solved"
  )
}

#' @export
as_marker_solved.default <- function(x) {
  abort_msg("`x` must be a marker_batch report, a marker_solved or a marker_panel.")
}
