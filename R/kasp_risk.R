# =============================================================================
#  kasp_risk.R  --  one experimenter-facing risk verdict over the checks
# =============================================================================
#  kasp_design_checks() is a long flat table; an experimenter wants a single
#  call.  Not every failed check is equal: no allele discrimination, an
#  off-target product, a variant in the AS 3' window or a saturated footprint
#  mean the assay can simply fail, whereas a shorter homopolymer or a small Tm
#  drift is a nuisance.  So the verdict is tiered, not a mean of everything
#  (matching the grain of the selection policy, which also gates the grave and
#  ranks the minor).  The reasons are always returned - never a bare label.
# =============================================================================

# A failed check here means the assay can fail outright, not just underperform.
KASP_RISK_HIGH <- c(
  "discrimination", "in-silico PCR (REF 1 / ALT 0)", "off-target exact primer hits (n)",
  "AS 3' window variants (n)", "AS 3' window dropout",
  "AS footprint dropout", "COM footprint dropout",
  "target site class", "ALT alleles (n)"
)

#' An experimenter-facing risk verdict for an assay
#'
#' Collapses [kasp_design_checks()] into one `low` / `medium` / `high` verdict
#' with the reasons.  A failure of a *grave* check (no discrimination, an
#' off-target product, a variant in the AS 3' window, a saturated footprint,
#' an unusable site) is `high`; any other failure (composition, structure, a
#' milder footprint or PCR deviation) is `medium`; all passing is `low`.  This
#' is a roll-up for a decision, not a probability of genotyping success.
#' @param assay result of [assemble_kasp_assay()] (or a `kasp_doomed` object,
#'   whose partial `evaluation` is used when present).
#' @param limits thresholds (see [kasp_design_checks()]).
#' @return a one-row tibble: `risk` (ordered factor `low` < `medium` < `high`),
#'   `n_failed` and `failed` (list-column of the failed criteria).
#' @examples
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' kasp_risk(assay)
#' }
#' @family KASP checks
#' @export
kasp_risk <- function(assay, limits = KASP_LIMITS) {
  ev <- assay$evaluation
  if (is.null(ev) || !nrow(ev)) {
    ev <- tryCatch(suppressWarnings(kasp_design_checks(assay, limits)), error = function(e) NULL)
  }
  if (is.null(ev) || !nrow(ev)) {
    return(tibble::tibble(
      risk = factor(NA_character_, levels = c("low", "medium", "high")),
      n_failed = NA_integer_, failed = list(character(0))
    ))
  }
  failed <- ev$criterion[!ev$pass]
  risk <- if (any(failed %in% KASP_RISK_HIGH)) {
    "high"
  } else if (length(failed)) {
    "medium"
  } else {
    "low"
  }
  tibble::tibble(
    risk = factor(risk, levels = c("low", "medium", "high")),
    n_failed = length(failed), failed = list(failed)
  )
}
