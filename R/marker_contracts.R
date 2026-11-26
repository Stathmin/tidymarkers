# =============================================================================
#  ZY12: design by contract on the critical paths
# =============================================================================
#  Pre/post conditions that must fail at the call site, not somewhere downstream:
#  a malformed target row should name the function that accepted it, not surface
#  later as an `NA` inside a primer search.  Base R only - this package adds a
#  dependency only for a real call, and `checkmate`/`assertthat` are not one.
#
#  The checks are gated: a session that loaded the package from source
#  (`pkgload::load_all()`) checks, an installed build does not.  The option
#  overrides both ways, which is also the escape hatch for a dev session (see
#  the roxygen of `contracts_enabled()`).
#
#  Every helper returns its input unchanged when the contracts are off, so a
#  call site reads `x <- contract_x(x, call)` whether or not the gate is open.

#' Is the design-by-contract checking on?
#'
#' On by default in a session that loaded the package from source (dev mode:
#' `"pkgload" %in% loadedNamespaces()`), off in an installed build.
#'
#' The escape hatch is the `tidymarkers.contracts` option, which overrides both
#' ways: set it to `FALSE` to switch the checks off while developing a caller
#' against the released behaviour, or to `TRUE` to keep them on in a deployed
#' build.  Leave it unset to get the dev-mode default.
#' @return a single logical.
#' @keywords internal
#' @noRd
contracts_enabled <- function() {
  isTRUE(getOption("tidymarkers.contracts", in_dev_mode()))
}

#' Is the package loaded from source (dev mode)?
#' @keywords internal
#' @noRd
in_dev_mode <- function() {
  "pkgload" %in% loadedNamespaces()
}

#' Abort on a broken contract, with a classed error naming the caller
#'
#' `isTRUE()` collapses `NA` and length > 1 to a violation.  The message is a
#' static `gettext()` string: `cli_abort()` interpolates `{...}` in it.
#' @param condition a single logical; `TRUE` means the invariant holds.
#' @param message the translated, static failure message.
#' @param ... named cli bullet hints; appended to `message` so `cli_abort()`
#'   renders them as bullets rather than dropping them into condition fields.
#' @param class condition class placed first, so a caller can catch it.
#' @param call the call shown to the user; the exported function's own call.
#' @keywords internal
#' @noRd
marker_contract <- function(condition, message, ...,
                            class = "tidymarkers_contract_error",
                            call = sys.call(-1L)) {
  if (isTRUE(condition)) {
    return(invisible(TRUE))
  }
  bullets <- list(...)
  if (length(bullets)) {
    nm <- names(bullets)
    if (is.null(nm)) nm <- rep("", length(bullets))
    vals <- vapply(bullets, function(x) as.character(gettext(x)), character(1))
    message <- c(message, stats::setNames(vals, nm))
  }
  cli::cli_abort(message, class = class, call = call)
}

#' Abort when a table or named list lacks required columns
#'
#' The one required-column guard: it names `arg` and the missing columns in the
#' order the caller listed them, and carries both condition classes so a caller
#' can catch either.
#' @param x a data frame, tibble or named list.
#' @param cols the required column names.
#' @param arg the name the message uses for `x`.
#' @param call the call shown to the user; defaults to the caller's frame.
#' @param ... optional cli bullet hints. Each hint is translated once by
#'   `marker_contract()` and rendered as a bullet; a `{...}` interpolation in a
#'   hint resolves in `marker_contract()`'s frame, not the caller's, so a hint
#'   that reads a caller value must be pre-formatted there. A hint cannot be
#'   named `x`, `cols`, `arg`, `call`, `class`, `message` or `condition`: those
#'   names match the formals of this helper or of `marker_contract()`.
#' @keywords internal
#' @noRd
need_cols <- function(x, cols, arg = "x", call = rlang::caller_env(), ...) {
  missing <- setdiff(cols, names(x))
  marker_contract(
    !length(missing),
    gettextf("`%s` is missing column(s): %s", arg, paste(missing, collapse = ", ")),
    ...,
    class = c("tidymarkers_contract_error", "tidymarkers_missing_columns"),
    call = call
  )
  invisible(TRUE)
}

# Post-condition of `pick_target_snp()`: the returned row names one locus whose
# AS window is designable at all.  `marker_site_class()` is the one definition of
# "designable" (no second one is invented here); the attached `"marker"`, when
# present, must itself be designable.  Every input is in the row, so no genome
# scan runs.
contract_target_row <- function(row, call) {
  if (!contracts_enabled()) {
    return(row)
  }
  marker_contract(
    is.data.frame(row) && nrow(row) == 1L && all(c("pos", "ref", "alt") %in% names(row)),
    gettext("Contract violated in pick_target_snp: the returned target is not a one-row pos/ref/alt locus"),
    call = call
  )
  pos <- row$pos[[1L]]
  ref <- as.character(row$ref)[1L]
  alt <- as.character(row$alt)[1L]
  marker_contract(
    !is.na(pos) && !is.na(ref) && !is.na(alt) && nzchar(ref) && nzchar(alt),
    gettext("Contract violated in pick_target_snp: the returned target has an unset pos/ref/alt"),
    call = call
  )
  marker_contract(
    isTRUE(marker_site_class(ref, marker_snp_alts(row)$alts_in)$disc[[1L]]),
    gettext("Contract violated in pick_target_snp: the returned SNP has no designable AS window"),
    call = call
  )
  marker <- attr(row, "marker")
  marker_contract(
    is.null(marker) || is_target_designable(marker),
    gettext("Contract violated in pick_target_snp: the attached marker is not designable"),
    call = call
  )
  row
}

# Pre-condition of `assemble_kasp_assay()`: a genome and a single locus whose
# REF/ALT can be designed at all.  `ref` is checked before `marker_site_class()`,
# which reads an `NA` REF as an indel and so would pass it.
contract_assay_input <- function(genome, snp, call) {
  if (!contracts_enabled()) {
    return(invisible(TRUE))
  }
  marker_contract(
    !is.null(genome),
    gettext("Contract violated in assemble_kasp_assay: no genome supplied"),
    call = call
  )
  marker_contract(
    is.data.frame(snp) && nrow(snp) == 1L && all(c("chrom", "pos", "ref", "alt") %in% names(snp)),
    gettext("Contract violated in assemble_kasp_assay: `snp` is not a one-row chrom/pos/ref/alt locus"),
    call = call
  )
  pos <- snp$pos[[1L]]
  marker_contract(
    is.numeric(pos) && length(pos) == 1L && !is.na(pos) && pos > 0,
    gettext("Contract violated in assemble_kasp_assay: the locus position is not a positive number"),
    call = call
  )
  ref <- as.character(snp$ref)[1L]
  marker_contract(
    !is.na(ref) && nzchar(ref),
    gettext("Contract violated in assemble_kasp_assay: the locus REF allele is unset"),
    call = call
  )
  alts_in <- marker_snp_alts(snp)$alts_in
  marker_contract(
    length(alts_in) > 0L && isTRUE(marker_site_class(ref, alts_in)$disc[[1L]]),
    gettext("Contract violated in assemble_kasp_assay: the locus ALT alleles are not designable"),
    call = call
  )
  invisible(TRUE)
}

# Post-condition of `assemble_kasp_assay()`: the assay carries the fields the
# function promises - a scheme, labelled AS/COM primers, an amplicon length.  A
# `marker_doomed` return is a legitimate answer, not a violation.
contract_assay_output <- function(out, call) {
  if (!contracts_enabled() || inherits(out, "marker_doomed")) {
    return(out)
  }
  marker_contract(
    inherits(out, "marker_assay"),
    gettext("Contract violated in assemble_kasp_assay: the result is not a marker_assay object"),
    call = call
  )
  marker_contract(
    is.character(out$scheme) && length(out$scheme) == 1L && !is.na(out$scheme) && nzchar(out$scheme),
    gettext("Contract violated in assemble_kasp_assay: the assay carries no scheme"),
    call = call
  )
  assay <- out$assay
  filled <- function(col) {
    v <- assay[[col]]
    is.character(v) && length(v) >= 1L && !is.na(v[[1L]]) && nzchar(v[[1L]])
  }
  marker_contract(
    is.data.frame(assay) && nrow(assay) >= 1L && filled("as_ref") && filled("as_alt") && filled("com"),
    gettext("Contract violated in assemble_kasp_assay: the assay has an empty AS/COM primer"),
    call = call
  )
  marker_contract(
    is.numeric(assay$amplicon) && length(assay$amplicon) >= 1L && is.finite(assay$amplicon[[1L]]),
    gettext("Contract violated in assemble_kasp_assay: the assay has no finite amplicon length"),
    call = call
  )
  tails <- out$tails
  labelled <- c("AS_ref", "AS_alt", "COM")
  marker_contract(
    is.character(tails) && !is.null(names(tails)) && all(labelled %in% names(tails)) &&
      all(!is.na(tails[labelled]) & nzchar(tails[labelled])),
    gettext("Contract violated in assemble_kasp_assay: the assay has no labelled AS/COM primer tails"),
    call = call
  )
  out
}
