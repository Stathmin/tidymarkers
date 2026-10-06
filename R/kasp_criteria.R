# =============================================================================
#  kasp_criteria.R  --  the one registry of everything we evaluate on a candidate
# =============================================================================
#  Selection (gates + ranking tiers) and acceptance (the "is it a good KASP?"
#  checks) used to be two disjoint vocabularies: the decision tree knew only the
#  gates, while kasp_design_checks() hardcoded its own list and several checks
#  (SSR, run, palindromes, off-targets) were recomputed for the winner only.
#
#  This registry is the single source: `kasp_criteria()` declares every criterion
#  with its stage, role, level, hardness and limit.  The selector, the acceptance
#  checks, the decision graph and the text export all read it, so a new check
#  appears in the tree and in the results at once.
# =============================================================================

# Human-readable limits for the selection gates (for the registry / text export).
KASP_GATE_LIMIT <- c(
  gate_discrimination = "discrimination >= min_disc",
  gate_structure = "hairpin/homodimer dG37 >= limit",
  gate_tm = "COM Tm >= AS Tm",
  gate_amplicon = "amplicon within [lo, hi]",
  gate_com_footprint = "COM dropout <= tol",
  gate_as_window = "no unlinked variant in AS 3' window",
  gate_wobble = "no weak G:T / A:C 3' mismatch",
  gate_as_footprint = "AS footprint dropout <= tol",
  gate_duplex_targets = "homologous (dKASP) or non-homologous (mxKASP)",
  gate_duplex_tm = "all primers Tm in [59, 65] C",
  gate_duplex_complement = "no >= 6 bp inter-pair complementarity",
  gate_duplex_selfcomp = "<= 5 self-complementary bases",
  gate_duplex_amplicon = "each amplicon <= 150 bp",
  gate_asq_cassette = "cassette per target (<= 4-plex)",
  gate_asq_readout = "allele resolved in real time",
  gate_caps_site = "an enzyme separates the REF and ALT amplicons",
  gate_caps_resolvable = "closest diagnostic bands >= 50 bp apart",
  gate_caps_amplicon = "amplicon within [120, 1200]",
  gate_caps_tm = "both primers Tm within [lo, hi]"
)

#' The registry of every candidate evaluation
#'
#' One row per criterion the pipeline applies: the selection gates and ranking
#' tiers of each stage, and the acceptance checks (composition, structure,
#' heteroduplexes, repeats, off-target specificity, in-silico PCR, site facts).
#' This is the single vocabulary the selector, [kasp_design_checks()], the
#' decision graph and [kasp_decision_text()] all read.  Each row carries the
#' marker `schemes` it serves ([kasp_schemes()]), so the same registry yields one
#' decision schema per chemistry.
#' @param limits thresholds (see [kasp_design_checks()]).
#' @param scheme optional scheme id (see [kasp_schemes()]); when given, only the
#'   criteria that scheme serves are returned.
#' @return tibble with `stage`, `role` (`gate`/`tier`/`check`), `id`, `label`,
#'   `level` (`site`/`feature`/`alignment`/`specificity`), `hard`, `limit`,
#'   `order`, `schemes` (list-column of scheme ids).
#' @examples
#' kasp_criteria()
#' kasp_criteria(scheme = "duplex")
#' @family KASP checks
#' @export
kasp_criteria <- function(limits = KASP_LIMITS, scheme = NULL) {
  tol <- limits$footprint_dropout %||% 0.05
  stage_rows <- function(stage, gates, relax, tiers, level, schemes) {
    n_gate <- length(gates)
    n_tier <- length(tiers)
    dplyr::bind_rows(
      tibble::tibble(
        stage = stage, role = "gate", id = paste0(stage, ":gate:", gates),
        label = unname(KASP_GATE_SHORT[gates]), level = level,
        hard = !(gates %in% relax), limit = unname(KASP_GATE_LIMIT[gates]),
        order = seq_along(gates), schemes = rep(list(schemes), n_gate)
      ),
      tibble::tibble(
        stage = stage, role = "tier", id = paste0(stage, ":tier:", names(tiers)),
        label = unname(tiers), level = level, hard = FALSE,
        limit = "maximise (per-tier bottleneck min)", order = seq_along(tiers),
        schemes = rep(list(schemes), n_tier)
      )
    )
  }
  checks <- tibble::tibble(
    stage = c(rep("site", 3L), rep("assay", 21L)),
    role = "check",
    label = c(
      "target site class", "ALT alleles (n)", "linked variants skipped (n)",
      "amplicon length (bp)", "AS Tm (C)", "COM Tm (C)", "COM Tm >= AS Tm",
      "AS GC", "COM GC", "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
      "inverted repeat (arms>=4, loop<=6)", "worst hairpin dG37",
      "worst cross-dimer dG37", "worst homodimer dG37",
      "AS 3' window variants (n)", "AS footprint variants (n)", "COM footprint variants (n)",
      "AS 3' window dropout", "AS footprint dropout", "COM footprint dropout", "discrimination",
      "in-silico PCR (REF 1 / ALT 0)", "off-target exact primer hits (n)"
    ),
    level = c(
      "site", "site", "site",
      check_level(c(
        "amplicon length (bp)", "AS Tm (C)", "COM Tm (C)", "COM Tm >= AS Tm",
        "AS GC", "COM GC", "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
        "inverted repeat (arms>=4, loop<=6)", "worst hairpin dG37",
        "worst cross-dimer dG37", "worst homodimer dG37",
        "AS 3' window variants (n)", "AS footprint variants (n)", "COM footprint variants (n)",
        "AS 3' window dropout", "AS footprint dropout", "COM footprint dropout", "discrimination",
        "in-silico PCR (REF 1 / ALT 0)", "off-target exact primer hits (n)"
      ))
    ),
    hard = c(
      TRUE, FALSE, FALSE,
      rep(FALSE, 21L)
    ),
    limit = c(
      "discriminable", ">= 1", "0",
      glue::glue("[{limits$amplicon[1]}, {limits$amplicon[2]}]"),
      glue::glue("[{limits$as_tm[1]}, {limits$as_tm[2]}]"),
      glue::glue("[{limits$com_tm[1]}, {limits$com_tm[2]}]"),
      ">= 0",
      glue::glue("[{limits$gc[1]}, {limits$gc[2]}]"),
      glue::glue("[{limits$gc[1]}, {limits$gc[2]}]"),
      glue::glue("< {limits$max_run}"),
      "0", "0",
      glue::glue(">= {limits$hairpin_dg}"),
      glue::glue(">= {limits$duplex_dg}"),
      glue::glue(">= {limits$homodimer_dg}"),
      "0", "0", "0",
      glue::glue("<= {tol}"), glue::glue("<= {tol}"), glue::glue("<= {tol}"),
      glue::glue(">= {limits$discrimination}"),
      "REF=1, ALT=0", "1 locus"
    )
  )
  checks$id <- paste0("check:", gsub("[^a-z0-9]+", "_", tolower(checks$label)))
  checks$order <- seq_len(nrow(checks))
  checks$schemes <- check_schemes(checks$label)
  out <- dplyr::bind_rows(
    stage_rows("as", KASP_GATES_AS, KASP_RELAX_AS, KASP_TIERS_AS, "feature", KASP_SCHEME_CORE),
    stage_rows("pair", KASP_GATES_PAIR, KASP_RELAX_PAIR, KASP_TIERS_PAIR, "alignment", KASP_SCHEME_CORE),
    stage_rows("duplex", KASP_GATES_DUPLEX, KASP_RELAX_DUPLEX, KASP_TIERS_DUPLEX, "alignment", "duplex"),
    stage_rows("asq", KASP_GATES_ASQ, KASP_RELAX_ASQ, KASP_TIERS_ASQ, "alignment", "asq"),
    stage_rows("caps", KASP_GATES_CAPS, KASP_RELAX_CAPS, KASP_TIERS_CAPS, "alignment", "caps"),
    checks
  )
  if (!is.null(scheme)) {
    keep <- vapply(out$schemes, function(s) scheme %in% s, logical(1))
    out <- out[keep, , drop = FALSE]
  }
  out
}

#' Evaluate every registered acceptance check on an assay, once
#'
#' The single entry point behind [kasp_design_checks()]: runs the registry's
#' `role == "check"` rows and returns one tidy row per criterion (`stage`,
#' `criterion`, `observed`, `limit`, `pass`, `detail`, `level`, `hard`).  Called
#' once at design time and stored on the assay (`assay$evaluation`), so plots and
#' reports read the results instead of recomputing them.
#' @keywords internal
#' @noRd
evaluate_assay <- function(assay, limits = KASP_LIMITS) {
  kasp_design_checks(assay, limits)
}
