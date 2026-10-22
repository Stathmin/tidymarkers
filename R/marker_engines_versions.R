# =============================================================================
#  marker_engines_versions.R  --  what produced a number (ZW2)
# =============================================================================
#  An evaluation is reproducible only if the engines that made it are named and
#  versioned.  DECIPHER's efficiency model, OligoArrayAux's dG, BLAST's hit
#  counts and TmCalculator's nearest-neighbour Tm all drift between releases --
#  and a drifted model can silently move a verdict.  Recording the versions
#  turns "the results changed" into "DECIPHER 3.8.1 -> 3.10.0 changed them".
#
#  The table is deliberately forgiving: an engine that is not installed is
#  reported with `version = NA` and `available = FALSE` rather than omitted, so
#  a missing engine is visible in the manifest instead of disappearing.
#
#  `marker_engine_stamp()` keeps that contract too.  An earlier version dropped
#  unavailable engines, which made a stamp that had silently lost BLAST
#  indistinguishable from one taken on a machine that never had it -- exactly
#  the failure this file exists to prevent.  It now always carries all seven
#  entries, with NA for a missing engine, and `marker_engine_stamp_field()`
#  flattens it to ONE string for a CSV or JSON column.
#
#  Roles are fixed here on purpose: the stamp format must not depend on which
#  engines happen to be installed, or two manifests would not line up.
# =============================================================================

#' Versions of every engine an Evaluation depends on
#'
#' @return a tibble with `engine`, `component`, `version`, `available` and
#'   `role` (what the engine decides in an evaluation).
#' @examples
#' marker_engine_versions()
#' @family KASP utilities
#' @export
marker_engine_versions <- function() {
  pkg_version <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      return(NA_character_)
    }
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
  }
  exe_version <- function(cmd, args = "--version") {
    path <- Sys.which(cmd)
    if (!nzchar(path)) {
      return(NA_character_)
    }
    out <- tryCatch(suppressWarnings(system2(path, args, stdout = TRUE, stderr = TRUE)),
      error = function(e) character()
    )
    if (!length(out)) {
      return(NA_character_)
    }
    # first line that carries a number ("blastn: 2.16.0+", "hybrid-min 3.8.1")
    hit <- grep("[0-9]+\\.[0-9]+", out, value = TRUE)
    if (!length(hit)) {
      return(NA_character_)
    }
    sub("^.*?([0-9]+\\.[0-9][0-9.]*[^ ]*).*$", "\\1", hit[[1L]])
  }
  oa <- tryCatch(oligoarrayaux_paths(), error = function(e) NULL)
  oa_bin <- if (is.null(oa) || !length(oa)) NA_character_ else oa[[1L]]

  tibble::tribble(
    ~engine, ~version, ~component, ~role,
    "DECIPHER", pkg_version("DECIPHER"),
    "CalculateEfficiencyPCR", "3' mismatch efficiency (discrimination, the Z7 model)",
    "OligoArrayAux", if (is.na(oa_bin)) NA_character_ else exe_version(oa_bin, "--version"),
    "hybrid-min", "duplex / hairpin / dimer free energy",
    "TmCalculator", pkg_version("TmCalculator"),
    "tm_calculate (tm_nn)", "nearest-neighbour Tm (salt- and DMSO-corrected)",
    "BLAST+", exe_version("blastn"),
    "blastn / makeblastdb", "off-target hit counts, in-silico PCR products",
    "Biostrings", pkg_version("Biostrings"),
    "vmatchPattern / PDict", "sequence I/O and pattern matching",
    "Rsamtools", pkg_version("Rsamtools"),
    "FaFile / faidx", "indexed FASTA access for a marker_genome",
    "R", as.character(getRversion()),
    "base + dplyr", "the language (a release can change dplyr tie-breaking)"
  ) |>
    dplyr::mutate(available = !is.na(.data$version), .after = "version")
}

# The engines a manifest must always name, in a stable order.  Adding one here
# changes the stamp format, so it is a deliberate act, not a side effect.
marker_engine_names <- c(
  "DECIPHER", "OligoArrayAux", "TmCalculator", "BLAST+",
  "Biostrings", "Rsamtools", "R"
)

#' A version stamp for an Evaluation or an artifact manifest
#'
#' `stamp` is the tight, serialisable form: a named character vector of
#' `engine = version`, plus the R version.  Every engine is present even when
#' it is not installed, in which case its value is `NA` -- a stamp that silently
#' omitted an engine could not be told apart from one taken without it.
#' Writing it next to a result makes a later disagreement attributable.
#'
#' @return a named character vector, one entry per engine (never empty).
#' @examples
#' marker_engine_stamp()
#' @family KASP utilities
#' @export
marker_engine_stamp <- function() {
  v <- marker_engine_versions()
  # index by name, not by row order, so a missing engine yields NA rather than
  # shifting its neighbour's version into its slot
  rlang::set_names(
    v$version[match(marker_engine_names, v$engine)],
    marker_engine_names
  )
}
#' Engines that produced an evaluation
#'
#' The stamp is stored as an attribute of the `$evaluation` table, so an
#' evaluation travels with the versions that made it.  `marker_evaluation_stamp()`
#' reads it back and returns `NULL` when the table carries none (an evaluation
#' built by hand, or by an older version); it never invents a stamp, because a
#' guessed stamp would claim provenance the numbers do not have.
#'
#' A `dplyr` verb or `bind_rows()` drops attributes, so use
#' [marker_with_stamp()] to put the stamp back after combining tables.
#'
#' @param x a checks tibble (see [marker_design_checks()]).
#' @return a named character vector (see [marker_engine_stamp()]), or `NULL`.
#' @examples
#' g <- ath_genome()
#' snp <- pick_target_snp(ensembl_variants(ensembl_gene("AT4G00650")$chrom, 1e6, 2e6), g, ensembl_gene("AT4G00650"))
#' marker_evaluation_stamp(marker_design_checks(snp))
#' @family KASP utilities
#' @export
marker_evaluation_stamp <- function(x) {
  attr(x, "engine_stamp", exact = TRUE)
}
#' Attach an engine stamp to a checks table
#'
#' @param x a checks tibble.
#' @param stamp a stamp from [marker_engine_stamp()]; computed at the top of the
#'   evaluation.  When extending a table that already carries one, pass it
#'   explicitly: a newer stamp on an older table would claim an engine ran on
#'   numbers it did not produce.
#' @param from a table being extended; its stamp is inherited when `stamp` is
#'   `NULL`.  The constructor stamps first and extends afterwards, so a table
#'   built row-by-row from `marker_evaluation_table()` can never carry none.
#' @return `x` with the stamp attribute set.
#' @keywords internal
#' @noRd
marker_with_stamp <- function(x, stamp = NULL, from = NULL) {
  if (is.null(x)) {
    return(x)
  }
  stamp <- stamp %||% marker_evaluation_stamp(from)
  if (is.null(stamp)) {
    # a table built from scratch (from = NULL) is stamped by whoever built it;
    # extending a table (from = ) that carries none stays unstamped rather than
    # acquiring a stamp for numbers it does not own
    if (is.null(from)) stamp <- marker_engine_stamp() else return(x)
  }
  attr(x, "engine_stamp") <- stamp
  x
}

#' A one-column form of the engine stamp
#'
#' `fmt = "csv"` gives `DECIPHER=3.8.1;OligoArrayAux=3.8.1;...`, safe to embed as
#' a single CSV cell.  `fmt = "json"` gives a JSON object, e.g. to nest under a
#' manifest's `"engines"` key.  An unavailable engine is written as the literal
#' `NA` rather than dropped, for the reason above.
#'
#' @param fmt `"csv"` (default) or `"json"`.
#' @param stamp a stamp from [marker_engine_stamp()]; recomputed if `NULL`.
#' @return a length-one character vector (JSON may be `NA` if it cannot encode).
#' @examples
#' marker_engine_stamp_field()
#' marker_engine_stamp_field("json")
#' @family KASP utilities
#' @export
marker_engine_stamp_field <- function(fmt = c("csv", "json"), stamp = NULL) {
  fmt <- match.arg(fmt)
  if (is.null(stamp)) stamp <- marker_engine_stamp()
  # NA -> "NA" so a missing engine stays visible in the serialised form
  vals <- stamp
  vals[is.na(vals)] <- "NA"
  if (identical(fmt, "csv")) {
    return(paste0(names(vals), "=", vals, collapse = ";"))
  }
  tryCatch(
    as.character(jsonlite::toJSON(as.list(vals), auto_unbox = TRUE)),
    error = function(e) NA_character_
  )
}
