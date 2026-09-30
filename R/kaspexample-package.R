#' kaspexample: KASP marker design for Arabidopsis thaliana
#'
#' An opinionated, reproducible pipeline for designing and quality-controlling
#' Kompetitive Allele-Specific PCR (KASP) markers.  Everything is done in R:
#' thermodynamics come from `TmCalculator`, `DECIPHER` and OligoArrayAux, sequence
#' work from `Biostrings`/`Rsamtools`, and the visuals from dedicated plotting
#' libraries (`ggmsa`, `gggenes`, `ggrepel`, `patchwork`).  See `vignette("kasp-intro")`.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

# Column names used inside dplyr/ggplot2 non-standard evaluation (covered by
# rlang's `.data` pronoun); listed so `R CMD check` does not flag them.
utils::globalVariables(c(".data"))
