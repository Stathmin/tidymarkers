# =============================================================================
#  marker_explain.R  --  "why this marker", the readable layer (Y7)
# =============================================================================
#  marker_design_checks() answers "what the marker scored"; this renders each
#  row as a sentence and orders failures first, so a reader sees why a marker
#  was accepted or rejected without reading a numeric table.  It is a pure view
#  over the check table, so it cannot disagree with the verdict.
# =============================================================================

#' Why this marker: the design checks as ordered, readable reasons
#'
#' A human-readable layer over [marker_design_checks()].  Failure and unknown
#' rows come first, so the first line a reader sees is the thing that stands
#' between the marker and acceptance; passing rows follow.  The `headline`
#' attribute summarises the counts in one line.
#'
#' @param x a designed marker, or a check table already returned by
#'   [marker_design_checks()].
#' @param scheme optional scheme id, used only when `x` is a marker and the
#'   marker does not carry its own.
#' @return a tibble with `target` (when the input has one), `criterion`,
#'   `status` and `reason`, classed `marker_explain`, with a `headline`
#'   attribute.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(system.file("extdata", "example_region_targets.csv", package = "tidymarkers"), colClasses = "character")[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' assay <- design_hrm_marker(genome, loci[1, ])
#' head(marker_explain(assay))
#' @seealso [marker_design_checks()], [marker_artefact()]
#' @family checks
#' @export
marker_explain <- function(x, scheme = NULL) {
  chk <- if (is.data.frame(x)) x else marker_design_checks(x, scheme = scheme)
  if (!nrow(chk)) {
    out <- tibble::tibble(
      criterion = character(), status = character(), reason = character()
    )
    return(structure(out, class = c("marker_explain", class(out)), headline = "no checks"))
  }
  status <- as.character(chk$status)
  rank <- match(status, c("fail", "unknown", "pass", "not_applicable"))
  keep <- !is.na(rank) & status != "not_applicable"
  chk <- chk[keep, , drop = FALSE]
  rank <- rank[keep]
  has_target <- "target" %in% names(chk)
  ord <- if (has_target) order(chk$target, rank) else order(rank)
  chk <- chk[ord, , drop = FALSE]
  status <- as.character(chk$status)

  reason <- vapply(
    seq_len(nrow(chk)),
    function(i) explain_reason(chk[i, , drop = FALSE], status[[i]]),
    character(1)
  )
  out <- chk[, intersect(c("target", "criterion", "status"), names(chk)), drop = FALSE]
  out$reason <- reason
  out <- tibble::as_tibble(out)

  tab <- table(factor(status, levels = c("pass", "fail", "unknown")))
  headline <- sprintf(
    "%d passed, %d failed, %d unknown", tab[["pass"]], tab[["fail"]], tab[["unknown"]]
  )
  structure(out, class = c("marker_explain", class(out)), headline = headline)
}

#' One check row as a sentence (internal)
#' @keywords internal
#' @noRd
explain_reason <- function(row, status) {
  crit <- as.character(row$criterion)
  observed <- suppressWarnings(as.numeric(row$observed))
  limit <- as.character(row$limit %||% NA_character_)
  detail <- as.character(row$detail %||% NA_character_)
  shown <- if (is.finite(observed)) {
    format(observed, trim = TRUE)
  } else if (!is.na(detail) && nzchar(detail)) {
    detail
  } else {
    NA_character_
  }
  inside <- !is.na(shown) && !is.na(limit) && nzchar(limit)
  switch(status,
    pass = if (inside) {
      sprintf("%s: %s is inside %s", crit, shown, limit)
    } else {
      sprintf("%s: satisfied", crit)
    },
    fail = if (inside) {
      sprintf("%s: %s is outside %s", crit, shown, limit)
    } else {
      sprintf("%s: failed", crit)
    },
    unknown = sprintf(
      "%s: not evaluated%s", crit,
      if (!is.na(detail) && nzchar(detail)) sprintf(" (%s)", detail) else ""
    ),
    sprintf("%s: %s", crit, status)
  )
}

#' @export
print.marker_explain <- function(x, ...) {
  cat(attr(x, "headline"), "\n", sep = "")
  NextMethod()
}
