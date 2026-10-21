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
    if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
  }
  exe_version <- function(cmd, args = "--version") {
    path <- Sys.which(cmd)
    if (!nzchar(path)) return(NA_character_)
    out <- tryCatch(suppressWarnings(system2(path, args, stdout = TRUE, stderr = TRUE)),
                    error = function(e) character())
    if (!length(out)) return(NA_character_)
    # first line that carries a number ("blastn: 2.16.0+", "hybrid-min 3.8.1")
    hit <- grep("[0-9]+\\.[0-9]+", out, value = TRUE)
    if (!length(hit)) return(NA_character_)
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

#' A version stamp for an Evaluation or an artifact manifest
#'
#' `stamp` is the tight, serialisable form: a named character vector of
#' `engine = version` (unavailable engines are omitted), plus the R version.
#' Writing it next to a result makes a later disagreement attributable.
#' @return a named character vector.
#' @examples
#' marker_engine_stamp()
#' @family KASP utilities
#' @export
marker_engine_stamp <- function() {
  v <- marker_engine_versions()
  v <- v[v$available, , drop = FALSE]
  rlang::set_names(v$version, v$engine)
}
