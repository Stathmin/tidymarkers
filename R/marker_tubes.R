# =============================================================================
#  marker_tubes.R  --  which designed markers can share one tube
# =============================================================================
#  Two perfectly good markers can be useless together: their primers cross-dimerise
#  or hairpin against each other.  `marker_primers()` is the one accessor that
#  pulls the physical oligos out of any designed marker, and `marker_tubes()` runs
#  the joint compatibility screen ([primer_set_compatibility()]) over every pair of
#  solved markers - the multiplex layer on top of the per-target design.
# =============================================================================

#' The oligos a designed marker would put in a tube
#'
#' One accessor over every chemistry's object, so the tube layer never has to know
#' each object's shape.  The default handles the shared flanking-pair shape
#' (`forward`/`reverse`); the methods cover the ones that differ (the KASP triple,
#' probes, the LAMP set, ...).  Returns a named character vector (empty when the
#' marker exposes no oligos).
#' @param marker a designed marker (a `marker_*` object).
#' @param ... unused.
#' @return a named character vector of oligos (may be empty).
#' @examples
#' \dontrun{
#' marker_primers(assemble_kasp_assay(genome, snp))
#' }
#' @family KASP plan
#' @export
marker_primers <- function(marker, ...) UseMethod("marker_primers")

#' @export
marker_primers.default <- function(marker, ...) {
  p <- marker$primers
  if (is.list(p) && all(c("forward", "reverse") %in% names(p))) {
    return(c(forward = p$forward, reverse = p$reverse))
  }
  character(0)
}

#' @export
marker_primers.marker_assay <- function(marker, ...) {
  a <- marker$assay
  c(AS_ref = a$as_ref[1L], AS_alt = a$as_alt[1L], COM = a$com[1L])
}

#' @export
marker_primers.marker_lax <- function(marker, ...) {
  if (!is.null(marker$pool)) {
    return(marker$pool)
  }
  if (!is.null(marker$base)) {
    return(marker_primers(marker$base))
  }
  # the no-haplotype form: one plain KASP assay per (single) ALT
  purrr::imap(marker$assays, function(a, i) {
    p <- marker_primers(a)
    stats::setNames(p, paste0(names(p), "_", i))
  }) |>
    unlist(use.names = TRUE)
}

#' @export
marker_primers.marker_duplex <- function(marker, ...) {
  purrr::imap(marker$assays, function(a, i) {
    p <- marker_primers(a)
    stats::setNames(p, paste0(names(p), "_", i))
  }) |>
    unlist(use.names = TRUE)
}

#' @export
marker_primers.marker_taqman <- function(marker, ...) {
  c(
    forward = marker$primers$forward, reverse = marker$primers$reverse,
    stats::setNames(marker$probes$probe, paste0("probe_", marker$probes$allele))
  )
}

#' @export
marker_primers.marker_dpcr <- function(marker, ...) {
  c(
    forward = marker$primers$forward, reverse = marker$primers$reverse,
    stats::setNames(marker$probes$probe, paste0("probe_", marker$probes$allele))
  )
}

#' @export
marker_primers.marker_castpcr <- function(marker, ...) {
  c(
    AS = marker$primers$as_primer, forward = marker$primers$forward,
    reverse = marker$primers$reverse, probe = marker$probe$sequence
  )
}

#' @export
marker_primers.marker_rpa <- function(marker, ...) {
  out <- c(forward = marker$primers$forward, reverse = marker$primers$reverse)
  if (!is.null(marker$probe)) out[["probe"]] <- marker$probe$sequence
  out
}

#' @export
marker_primers.marker_snapshot <- function(marker, ...) {
  c(
    forward = marker$amplification$forward, reverse = marker$amplification$reverse,
    extension = marker$extension$primer
  )
}

#' @export
marker_primers.marker_lamp <- function(marker, ...) unlist(marker$primers)

#' @export
marker_primers.marker_arms <- function(marker, ...) {
  c(
    outer_forward = marker$outer$forward, outer_reverse = marker$outer$reverse,
    AS_ref = marker$inner$ref$primer, AS_alt = marker$inner$alt$primer
  )
}

#' @export
marker_primers.marker_crispr <- function(marker, ...) {
  c(forward = marker$amplicon$forward, reverse = marker$amplicon$reverse, guide = marker$guide)
}

#' Which solved markers can share one tube without interfering
#'
#' The multiplex layer: for every pair of *solved* markers it pools their oligos
#' and runs the joint compatibility screen ([primer_set_compatibility()] over all
#' of them at once), so a pair whose primers cross-dimerise or hairpin is flagged
#' before it reaches the bench.  Markers exposing no oligos are skipped.
#' @param solved a `marker_solved` (see [marker_solve()]).
#' @param min_cross_run,max_selfcomp the base-pair rules (see
#'   [primer_set_compatibility()]).
#' @return a tibble, one row per marker pair: `a`, `b`, `n_pairs` (inter-marker
#'   oligo pairs checked), `worst_run` (the longest inter-primer complementary
#'   stretch) and `pass`.
#' @examples
#' \dontrun{
#' marker_tubes(marker_solve(plan))
#' }
#' @family KASP plan
#' @export
marker_tubes <- function(solved, min_cross_run = 6L, max_selfcomp = 5L) {
  if (!inherits(solved, "marker_solved")) {
    cli::cli_abort(gettext("`solved` must be a marker_solved (see marker_solve())."))
  }
  r <- solved$report
  keep <- which(r$ok)
  empty <- tibble::tibble(
    a = character(), b = character(), n_pairs = integer(), worst_run = integer(), pass = logical()
  )
  if (length(keep) < 2L) {
    return(empty)
  }
  labels <- r$label[keep]
  pools <- purrr::map(r$marker[keep], marker_primers)
  combos <- utils::combn(length(keep), 2L)
  out <- purrr::map(seq_len(ncol(combos)), function(j) {
    i1 <- combos[1L, j]
    i2 <- combos[2L, j]
    p1 <- pools[[i1]]
    p2 <- pools[[i2]]
    if (!length(p1) || !length(p2)) {
      return(NULL)
    }
    pooled <- c(
      stats::setNames(p1, paste0(labels[i1], ":", names(p1))),
      stats::setNames(p2, paste0(labels[i2], ":", names(p2)))
    )
    # the tube verdict is the *inter-marker* compatibility; each primer's own
    # self-complementarity is a single-primer property already screened at design
    # time, so it does not belong in a pairwise "can these share a tube" answer
    cc <- primer_set_compatibility(pooled, min_cross_run = min_cross_run, max_selfcomp = max_selfcomp)
    cross <- cc[cc$kind == "cross", , drop = FALSE]
    tibble::tibble(
      a = labels[i1], b = labels[i2], n_pairs = nrow(cross),
      worst_run = if (nrow(cross)) max(cross$run) else 0L, pass = all(cross$pass)
    )
  }) |>
    purrr::list_rbind()
  if (nrow(out)) out else empty
}
