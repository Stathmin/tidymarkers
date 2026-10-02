# =============================================================================
#  kasp_checks.R  --  acceptance checks (driven by the criterion registry)
# =============================================================================

#' Longest homopolymer run in a sequence
#' @keywords internal
#' @noRd
longest_run <- function(seq) {
  max(nchar(stringr::str_extract_all(seq, "([ACGT])\\1*")[[1L]]))
}

#' Check an assembled assay against the KASP acceptance criteria
#'
#' Turns the prose design manual into a tidy pass/fail table so "good" and
#' "bad" assays are graded the same way every time (thresholds: `KASP_LIMITS`).
#' The criterion list, levels and limits come from [kasp_criteria()] (its
#' `role == "check"` rows), so the checks and the decision tree share one
#' vocabulary.  Reads only the assay object produced by [assemble_kasp_assay()],
#' so it is pure and cheap to test.
#'
#' @param assay result of [assemble_kasp_assay()].
#' @param limits list of thresholds; defaults to `KASP_LIMITS`.
#' @return tibble with `stage`, `criterion`, `observed`, `limit`, `pass`,
#'   `detail` and `level`.  Criteria cover the target *site*, each *feature*
#'   oligo, the *alignment*, and *off-target* specificity.
#' @examples
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' kasp_design_checks(assay)
#' }
#' @family KASP checks
#' @export
kasp_design_checks <- function(assay, limits = KASP_LIMITS) {
  aa <- assay$assay[1L, ]
  th <- assay$thermo
  d <- assay$complementarity
  disc <- assay$efficiency$discrimination[1L]
  tol <- limits$footprint_dropout %||% 0.05

  as_tm <- mean(th$tm[th$name %in% c("AS_ref", "AS_alt")])
  com_tm <- aa$com_tm
  hair <- min(th$hairpin_dg)
  cross <- min(d$dg[d$a != d$b]) # worst hetero-dimer
  homo <- min(d$dg[d$a == d$b]) # worst self-dimer
  run <- max(longest_run(aa$as_ref), longest_run(aa$com))

  # in-silico PCR (only if computed): REF primer -> 1 product, ALT primer -> 0
  amp <- assay$amplicons
  amp_obs <- if (is.null(amp)) NA_real_ else max(amp$n_products)
  amp_pass <- if (is.null(amp)) {
    FALSE
  } else {
    ref <- amp$n_products[amp$primer == "AS_ref"]
    alt <- amp$n_products[amp$primer == "AS_alt"]
    length(ref) == 1L && length(alt) == 1L && ref == 1L && alt == 0L
  }

  # off-target: the AS primers' 3'-seed must be unique (one locus), when screened
  sp <- assay$specificity
  ot_obs <- if (is.null(sp) || !"hits" %in% names(sp)) NA_real_ else suppressWarnings(max(sp$hits, na.rm = TRUE))
  ot_pass <- is.na(ot_obs) || ot_obs <= 1

  # site-level facts (absent on an older assay)
  site <- assay$site
  site_class <- if (is.null(site)) "biallelic" else as.character(site$class[1L])
  site_disc <- if (is.null(site)) TRUE else isTRUE(site$disc[1L])
  n_alt <- if (is.null(site)) 1L else as.integer(site$n_alt[1L])
  linked_skipped <- as.integer(assay$selection$linked_skipped %||% 0L)

  obs <- c(
    NA_real_, as.numeric(n_alt), as.numeric(linked_skipped),
    aa$amplicon, as_tm, com_tm, com_tm - as_tm, aa$as_gc, aa$com_gc,
    run, sum(lengths(find_repeats(c(aa$as_ref, aa$com)))),
    sum(purrr::map_int(find_palindromes(c(aa$as_ref, aa$as_alt, aa$com)), nrow)),
    hair, cross, homo, aa$as_n_3p, aa$n_var, aa$as_dropout, aa$com_dropout,
    disc, amp_obs, ot_obs
  )
  pass <- c(
    site_disc, n_alt >= 1L, linked_skipped == 0L,
    dplyr::between(aa$amplicon, limits$amplicon[1], limits$amplicon[2]),
    dplyr::between(as_tm, limits$as_tm[1], limits$as_tm[2]),
    dplyr::between(com_tm, limits$com_tm[1], limits$com_tm[2]),
    com_tm >= as_tm,
    dplyr::between(aa$as_gc, limits$gc[1], limits$gc[2]),
    dplyr::between(aa$com_gc, limits$gc[1], limits$gc[2]),
    run < limits$max_run,
    !has_repeat(aa$as_ref) && !has_repeat(aa$com),
    !any(has_palindrome(c(aa$as_ref, aa$as_alt, aa$com))),
    hair >= limits$hairpin_dg,
    cross >= limits$duplex_dg,
    homo >= limits$homodimer_dg,
    !is.na(aa$as_dropout) | aa$as_n_3p == 0L,
    !is.na(aa$com_dropout) | aa$n_var == 0L,
    is.na(aa$as_dropout) | aa$as_dropout <= tol,
    is.na(aa$com_dropout) | aa$com_dropout <= tol,
    disc >= limits$discrimination,
    amp_pass, ot_pass
  )
  detail <- c(site_class, rep(NA_character_, length(obs) - 1L))

  crit <- kasp_criteria(limits)
  crit <- crit[crit$role == "check", ]
  tibble::tibble(
    stage = crit$stage, criterion = crit$label, observed = obs,
    limit = crit$limit, pass = pass, detail = detail, level = crit$level
  )
}
