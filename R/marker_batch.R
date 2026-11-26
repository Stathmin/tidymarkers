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
    missing <- setdiff(c("chrom", "pos", "ref", "alt"), names(loci))
    if (length(missing)) {
      abort_msg("`loci` is missing column(s): {paste(missing, collapse = ', ')}")
    }
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
