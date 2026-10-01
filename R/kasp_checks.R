# =============================================================================
#  kasp_checks.R  --  acceptance checks
# =============================================================================
# =============================================================================
#  9. Acceptance checks - "is this a good KASP?"
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
#' Reads only the assay object produced by [assemble_kasp_assay()], so it is
#' pure and cheap to test.
#'
#' @param assay result of [assemble_kasp_assay()].
#' @param limits list of thresholds; defaults to `KASP_LIMITS`.
#' @return tibble with `criterion`, `observed`, `limit` and `pass`.
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

  tibble::tibble(
    criterion = c(
      "amplicon length (bp)", "AS Tm (C)", "COM Tm (C)",
      "COM Tm >= AS Tm", "AS GC", "COM GC",
      "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
      "inverted repeat (arms>=4, loop<=6)", "worst hairpin dG37",
      "worst cross-dimer dG37", "worst homodimer dG37",
      "AS 3' window variants (n)", "COM footprint variants (n)",
      "AS 3' window dropout", "COM footprint dropout",
      "discrimination", "in-silico PCR (REF 1 / ALT 0)"
    ),
    observed = c(
      aa$amplicon, as_tm, com_tm, com_tm - as_tm, aa$as_gc, aa$com_gc,
      run, sum(lengths(find_repeats(c(aa$as_ref, aa$com)))),
      sum(purrr::map_int(find_palindromes(c(aa$as_ref, aa$as_alt, aa$com)), nrow)),
      hair, cross, homo, aa$as_n_3p, aa$n_var, aa$as_dropout, aa$com_dropout, disc, amp_obs
    ),
    limit = c(
      glue::glue("[{limits$amplicon[1]}, {limits$amplicon[2]}]"),
      glue::glue("[{limits$as_tm[1]}, {limits$as_tm[2]}]"),
      glue::glue("[{limits$com_tm[1]}, {limits$com_tm[2]}]"),
      ">= 0",
      glue::glue("[{limits$gc[1]}, {limits$gc[2]}]"),
      glue::glue("[{limits$gc[1]}, {limits$gc[2]}]"),
      glue::glue("< {limits$max_run}"),
      "0",
      "0",
      glue::glue(">= {limits$hairpin_dg}"),
      glue::glue(">= {limits$duplex_dg}"),
      glue::glue(">= {limits$homodimer_dg}"),
      "0",
      "0",
      glue::glue("<= {tol}"),
      glue::glue("<= {tol}"),
      glue::glue(">= {limits$discrimination}"),
      "REF=1, ALT=0"
    ),
    pass = c(
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
      amp_pass
    )
  )
}
