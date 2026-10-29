# =============================================================================
#  marker_performance.R  --  analytical performance, one quantity per field
# =============================================================================
#  "Detectable" is not one number (Z8.5).  A dPCR assay reports at least a
#  limit of blank, a limit of detection and a limit of quantification, and a
#  method paper may add a minor-allele fraction, a false-positive droplet rate
#  and a partition count.  They answer different questions against different
#  denominators, so folding them into one "sensitivity" produces a number with
#  no meaning.  The decision report names each quantity on its own line
#  (`marker_decision_text()`) and the artefact records the names in
#  `decisions.csv`, so a later reader knows which limit a claim refers to.
#
#  The measured values stay in the sourced tables (jak2v617f_ddpcr_limits.csv
#  and the other `dpcr_*_limits.csv` files), where each row already carries its
#  own definition and evidence tier.  This file holds only the vocabulary, so
#  the report cannot invent a merged column.
# =============================================================================

#' Analytical quantities, one definition each
#'
#' Kept as a named vector so a scheme declares *which* quantities it reports and
#' the report prints the same definition everywhere.  A new quantity must be
#' added here first; [scheme_analytical()] then rejects a scheme that names an
#' unknown one.
#' @keywords internal
#' @noRd
ANALYTICAL_QUANTITIES <- c(
  LoB = "limit of blank: the highest signal a blank sample reaches",
  LoD = "limit of detection: the lowest level still called positive",
  LoQ = "limit of quantification: the lowest level quantified within a stated CV",
  MAF = "minor-allele fraction: the lowest variant fraction the assay resolves",
  false_positive_rate = "false-positive rate: the fraction of negatives called positive",
  partition_count = "partition count: accepted droplets per reaction"
)

#' The analytical quantities a scheme's report must keep separate
#'
#' Reads the scheme's `analytical` field from its `scheme_def_*()` spec, because
#' [marker_schemes()] deliberately exposes only the fixed comparison columns.
#' @param scheme a scheme id (see [marker_schemes()]).
#' @return a named character vector (quantity -> definition), or `NULL` when the
#'   scheme declares none.
#' @keywords internal
#' @noRd
scheme_analytical <- function(scheme) {
  defs <- extra_scheme_defs()
  d <- if (scheme %in% names(defs)) defs[[scheme]] else NULL
  ids <- d$analytical
  if (!length(ids)) {
    return(NULL)
  }
  unknown <- setdiff(ids, names(ANALYTICAL_QUANTITIES))
  if (length(unknown)) {
    stop("scheme '", scheme, "' declares unknown analytical quantities: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  ANALYTICAL_QUANTITIES[ids]
}
