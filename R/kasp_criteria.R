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

# The human-readable gate limits (`KASP_GATE_LIMIT`) live in R/kasp_select.R with
# the other gate lookups, so each scheme's entries sit in its own scheme file.

# How each not-gated *check* is named inside a gate (its id / label / limit),
# normalised the same way the DAG helper's C10 check does.  Only a check whose
# stem appears in one of the scheme's gates is really gated.
KASP_UNGATED_STEM <- c(
  "inverted repeat (arms>=4, loop<=6)" = "invertedrepeat",
  "in-silico PCR (REF 1 / ALT 0)" = "insilico",
  "off-target exact primer hits (n)" = "offtarget"
)

# Acceptance checks the *selection* does not gate (they need the genome or an
# in-silico PCR pass, computed only after the pair is chosen).  A failed one is
# reported as a warning by assemble_kasp_assay(), not returned in silence.
#
# The set is per-scheme (C10): a scheme whose stage *does* gate one of them
# (indel_size's uniqueness gate covers the PCR check, lamp's structure gate the
# inverted repeat, ssr's the off-target) must not still call it ungated - that
# would claim the selection ignores a criterion it enforces.  With no scheme the
# full un-filtered set is returned.
KASP_UNGATED_CHECKS <- function(scheme = NULL) {
  all <- names(KASP_UNGATED_STEM)
  if (is.null(scheme)) {
    return(all)
  }
  gate <- kasp_criteria(scheme = scheme)
  gate <- gate[gate$role == "gate", ]
  evidence <- gsub(
    "[^a-z0-9]", "",
    tolower(paste(c(gate$id, gate$label, gate$limit), collapse = " "))
  )
  all[!vapply(all, function(l) grepl(KASP_UNGATED_STEM[[l]], evidence, fixed = TRUE), logical(1))]
}

# The amplicon acceptance window a family is graded against.
#
# The `amplicon length (bp)` check is the marker-agnostic primer one, but a single
# 50-100 bp EasyKASP optimum is wrong for every chemistry that reads the product
# by *size*: a CAPS/ARMS band, a Sanger read, an SSR / InDel-size fragment, an
# isothermal LAMP/RPA/CRISPR product, or a short dPCR / TaqMan amplicon.  Each
# family's window is read from its own constant (or its designer's declared
# default), so the selection gate and the acceptance check cannot drift (C2/C5).
#
# A function, not a top-level constant: R sources kasp_criteria.R *before* the
# scheme files, so reading `KASP_LIMITS$amplicon` here at load time fails (the
# bug that already bit an agent); resolving the lookup lazily defers it until the
# namespace is complete.  `kasp` / `multiplex` / `asq` / `primer` - and `hrm` /
# `castpcr`, whose gate window *is* the KASP optimum - keep the marker-agnostic
# window, so the SNP/indel/multi/duplex/asq/primer behaviour is unchanged.
kasp_amplicon_window <- function(family = "kasp", limits = KASP_LIMITS) {
  product_window <- function(f) {
    fm <- formals(f)
    c(eval(fm$min_product), eval(fm$max_product))
  }
  switch(family,
    caps = caps_limits()$amplicon,
    arms = as.numeric(KASP_ARMS_WINDOW$control),
    crispr = product_window(design_crispr_marker),
    dpcr = product_window(design_dpcr_marker),
    indel_size = KASP_INDEL_SIZE_LIMITS$amplicon,
    lamp = product_window(design_lamp_marker),
    rpa = RPA_AMPLICON,
    sanger = KASP_SANGER_LIMITS$amplicon,
    snapshot = SNAPSHOT_AMPLICON,
    ssr = KASP_SSR_LIMITS$amplicon,
    taqman = product_window(design_taqman_marker),
    limits$amplicon
  )
}

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
  # `stage` / `role` / `level` are recycled over the gate (and tier) count so a
  # stage with `gates = character(0)` (a documentation-only stage, e.g. ASQ)
  # yields a 0-row gate block instead of a length-mismatch error.
  stage_rows <- function(stage, gates, relax, tiers, level, schemes) {
    n_gate <- length(gates)
    n_tier <- length(tiers)
    # `id` is built with an explicit `if` so an empty gate/tier vector gives a
    # length-0 id, never the dangling `"<stage>:gate:"` a bare paste0 would make
    dplyr::bind_rows(
      tibble::tibble(
        stage = rep(stage, n_gate), role = rep("gate", n_gate),
        id = if (n_gate) paste0(stage, ":gate:", gates) else character(0),
        label = unname(KASP_GATE_SHORT[gates]), level = rep(level, n_gate),
        hard = !(gates %in% relax), limit = unname(KASP_GATE_LIMIT[gates]),
        order = seq_along(gates), schemes = rep(list(schemes), n_gate)
      ),
      tibble::tibble(
        stage = rep(stage, n_tier), role = rep("tier", n_tier),
        id = if (n_tier) paste0(stage, ":tier:", names(tiers)) else character(0),
        label = unname(tiers), level = rep(level, n_tier), hard = FALSE,
        limit = "maximise (per-tier bottleneck min)", order = seq_along(tiers),
        schemes = rep(list(schemes), n_tier)
      )
    )
  }
  checks <- tibble::tibble(
    stage = c(rep("site", 5L), rep("assay", 21L)),
    role = "check",
    label = c(
      "target site class", "ALT alleles (n)", "linked variants skipped (n)",
      "ALT primer covers all ALTs", "design relies on linkage",
      "amplicon length (bp)", "AS Tm (C)", "COM Tm (C)", "COM Tm >= AS Tm",
      "AS GC", "COM GC", "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
      "inverted repeat (arms>=4, loop<=6)", "worst hairpin dG37",
      "worst cross-dimer dG37", "worst homodimer dG37",
      "AS 3' window variants (n)", "AS footprint variants (n)", "COM footprint variants (n)",
      "AS 3' window dropout", "AS footprint dropout", "COM footprint dropout", "discrimination",
      "in-silico PCR (REF 1 / ALT 0)", "off-target exact primer hits (n)"
    ),
    level = c(
      "site", "site", "site", "site", "site",
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
      TRUE, FALSE, FALSE, FALSE, FALSE,
      rep(FALSE, 21L)
    ),
    limit = c(
      "discriminable", ">= 1", "0", "0", "0",
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
  kasp_core <- c(KASP_SCHEME_CORE, extra_family_schemes("kasp"))
  out <- dplyr::bind_rows(
    stage_rows("as", KASP_GATES_AS, KASP_RELAX_AS, KASP_TIERS_AS, "feature", kasp_core),
    stage_rows("pair", KASP_GATES_PAIR, KASP_RELAX_PAIR, KASP_TIERS_PAIR, "alignment", kasp_core),
    stage_rows("duplex", KASP_GATES_DUPLEX, KASP_RELAX_DUPLEX, KASP_TIERS_DUPLEX, "alignment", "duplex"),
    stage_rows("asq", KASP_GATES_ASQ, KASP_RELAX_ASQ, KASP_TIERS_ASQ, "alignment", "asq"),
    stage_rows("caps", KASP_GATES_CAPS, KASP_RELAX_CAPS, KASP_TIERS_CAPS, "alignment", "caps"),
    stage_rows("arms", KASP_GATES_ARMS, KASP_RELAX_ARMS, KASP_TIERS_ARMS, "alignment", "arms"),
    stage_rows("hrm", KASP_GATES_HRM, KASP_RELAX_HRM, KASP_TIERS_HRM, "alignment", "hrm"),
    checks,
    extra_stage_rows()
  )
  if (!is.null(scheme)) {
    keep <- vapply(out$schemes, function(s) scheme %in% s, logical(1))
    out <- out[keep, , drop = FALSE]
    # the amplicon acceptance window is per-family (C2/C5): a product read by size
    # is graded against the scheme's own window, never the KASP 50-100 bp optimum
    win <- kasp_amplicon_window(kasp_scheme(scheme)$family, limits)
    i <- match("amplicon length (bp)", out$label)
    if (!is.na(i)) out$limit[i] <- glue::glue("[{win[1L]}, {win[2L]}]")
  }
  out
}
