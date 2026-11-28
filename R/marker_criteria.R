# =============================================================================
#  marker_criteria.R  --  the one registry of everything we evaluate on a candidate
# =============================================================================
#  Selection (gates + ranking tiers) and acceptance (the "is it a good KASP?"
#  checks) used to be two disjoint vocabularies: the decision tree knew only the
#  gates, while marker_design_checks() hardcoded its own list and several checks
#  (SSR, run, palindromes, off-targets) were recomputed for the winner only.
#
#  This registry is the single source: `marker_criteria()` declares every criterion
#  with its stage, role, level, hardness and limit.  The selector, the acceptance
#  checks, the decision graph and the text export all read it, so a new check
#  appears in the tree and in the results at once.
#
#  Z2: every criterion is also a *data object* (`criterion_spec()`), so the
#  threshold itself - not only its prose - lives in one place.  `limit` is the
#  rendering of that object at the thresholds in force; `bounds` is the threshold.
#  The gate label / check-limit maps below are views over the object, never a
#  second copy of the numbers.
# =============================================================================

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
  gate <- marker_criteria(scheme = scheme)
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
# A function, not a top-level constant: R sources marker_criteria.R *before* the
# scheme files, so reading `KASP_LIMITS$amplicon` here at load time fails (the
# bug that already bit an agent); resolving the lookup lazily defers it until the
# namespace is complete.  `kasp` / `multiplex` / `asq` / `primer` - and `hrm` /
# `castpcr`, whose gate window *is* the KASP optimum - keep the marker-agnostic
# window, so the SNP/indel/multi/duplex/asq/primer behaviour is unchanged.
marker_amplicon_window <- function(family = "kasp", limits = KASP_LIMITS) {
  product_window <- function(f) {
    fm <- formals(f)
    c(eval(fm$min_product), eval(fm$max_product))
  }
  switch(family,
    caps = caps_limits()$amplicon,
    arms = as.numeric(KASP_ARMS_WINDOW$control),
    asq = ASQ_AMPLICON,
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

# ---------------------------------------------------------------------------
#  the criterion data object (Z2)
# ---------------------------------------------------------------------------
# Every criterion is one row: an `id` (stable, snake_case), where it applies
# (`scope`: a scheme id, `"kasp"` for the shared KASP family stages, or
# `"family"` for the per-family acceptance checks), the `stage`/`role` that place
# it in the pipeline, the `level` it reports at, a `kind` (`numeric` /
# `categorical` / `boolean`), the `bounds` it is judged against, the `units` they
# are in, where the number comes from (`source`), the strength of the evidence
# behind it (`evidence`, see `CRITERION_EVIDENCE`) and the `evaluator` that
# produces the observed value.  `calibrated` is a deprecated boolean view of the
# `evidence` rung (ZW4).
#
# `bounds` is the *threshold*, a list-column resolved lazily from the limit
# object so a caller's `limits=` is honoured exactly as before:
#   numeric     c(min, max) or c(bound, NA) in `units` (a `>=` / `<` / `<=` bar);
#               c(bound, "exclusive") is the strict `<` the homopolymer bar uses;
#               "hard" marks a gate with no floor (it can never be relaxed);
#   categorical the allowed levels (their text is the `limit` string when there
#               is no tighter constraint, e.g. "TTTV");
#   boolean     the tested condition as text (the `limit` IS the predicate).
# The selector, the checker, the plotter and the decision graph read this object:
# `marker_criteria()`'s `limit` column is only its rendering, so the same
# criterion cannot carry two different bars.

# The one textual limit a ranking tier carries (the tiers are maximised, so their
# bar is the ranking itself, not a threshold on the observed value).
CRITERION_TIER_LIMIT <- "maximise (per-tier bottleneck min)"

# ZW4.6: a `*:tier:*` row is two quantities, and only one of them is citable.
#   * the TARGET - the reference point the tier ranks candidates toward (a
#     window, a dG floor, a design rule).  A paper or a vendor manual does print
#     this, so it carries the row's `evidence` rung and is stored in `target`.
#   * the AGGREGATION POLICY - bottleneck `min` within a tier, tier membership,
#     lexicographic order across tiers.  No source prints it, so it stays
#     `rule_of_thumb` by construction.  It is declared once, here, and
#     `order_by_margins()` / `select_with_ladder()` implement it.
# Before the split the two shared one `limit` string, which pinned every tier at
# `rule_of_thumb` even when its target was fully sourced.
CRITERION_TIER_POLICY_EVIDENCE <- "rule_of_thumb"

# Criterion kinds: a numeric bar, a set of allowed levels, or a tested predicate.
CRITERION_KINDS <- c("numeric", "categorical", "boolean")

# The evidence ladder, ascending: where a threshold comes from.  The order is the
# comparison (`>=` is meaningful), so the rungs are held in exactly this order.
CRITERION_EVIDENCE <- c(
  "rule_of_thumb", "vendor_guideline", "primary_literature", "fitted_here"
)

# The evidence tiers a provenance table may record (the `evidence_tier` column of
# inst/extdata/snapshot_multiplex_design_rules.csv), and the rung each maps to.
# The mapping only ever goes DOWN the ladder (a tier never upgrades a rung):
# `none` states a rule with no quantitative outcome, so it is a rule of thumb.
CRITERION_TIER_EVIDENCE <- c(
  moderate = "primary_literature",
  weak = "vendor_guideline",
  none = "rule_of_thumb"
)

#' The evidence rung of a criterion
#'
#' A threshold is only as good as where it came from, so every criterion row
#' carries its evidence on an ordered ladder:
#' `"rule_of_thumb"` < `"vendor_guideline"` < `"primary_literature"` <
#' `"fitted_here"`.  This is the single definition of that vocabulary.
#'
#' The result is an ordered factor with `levels = CRITERION_EVIDENCE`, so `>`
#' compares rungs and printing names the one in force.  `NULL` and `NA` are `NA`
#' (unset); anything else must be one of the four machine values - an unknown
#' string aborts rather than silently degrading to the weakest rung.
#'
#' **The migration from `calibrated`:** `calibrated = TRUE` means
#' `"fitted_here"`, and `calibrated = FALSE` means `"rule_of_thumb"` (nothing
#' weaker was ever recorded).  `calibrated = NA` means `NA`.  The boolean survives
#' only as a deprecated read of this ladder.
#' @param x a character vector of evidence rungs, or an existing factor.
#' @return an ordered factor with levels `CRITERION_EVIDENCE`.
#' @keywords internal
#' @noRd
criterion_evidence <- function(x = NULL) {
  if (is.null(x)) {
    x <- NA_character_
  }
  x <- as.character(x)
  bad <- unique(x[!is.na(x) & !x %in% CRITERION_EVIDENCE])
  if (length(bad)) {
    abort_msg(
      "Unknown criterion evidence: {.val {bad}}.",
      i = "Allowed values, weakest first: {.val {CRITERION_EVIDENCE}}."
    )
  }
  factor(x, levels = CRITERION_EVIDENCE, ordered = TRUE)
}

# The one migration rule `calibrated` -> `evidence` (see `criterion_evidence()`).
criterion_evidence_from_flag <- function(calibrated) {
  calibrated <- as.logical(calibrated)
  criterion_evidence(ifelse(calibrated, "fitted_here", "rule_of_thumb"))
}


# A whole number without the trailing `.0` a threshold never needs.
criterion_num <- function(x) {
  if (length(x) != 1L || is.na(x)) {
    return(NA_character_)
  }
  if (is.numeric(x) && x == round(x)) {
    return(as.character(as.integer(round(x))))
  }
  as.character(x)
}

# The bar a numeric bound states: an interval, or a one-sided bound (min, MAX) /
# (MIN, max) / strictly below max.  This is the one place a `bounds` vector
# becomes prose, so `limit` and every label derived from it agree (Z2: the
# check totals below are pinned by tests to the old prose exactly).
criterion_bound_text <- function(bounds) {
  if (length(bounds) < 2L) {
    return(paste0("= ", criterion_num(bounds[[1L]])))
  }
  lo <- bounds[[1L]]
  hi <- bounds[[2L]]
  if (isTRUE(hi == "exclusive") && length(bounds) == 3L) {
    return(glue::glue("< {criterion_num(bounds[[3L]])}"))
  }
  # An open end is `NA` in a parsed bound and `Inf` in a live `bounds` vector;
  # `criterion_num(Inf)` coerced to `NA_integer_` and rendered `[20, NA]`.
  if (is.na(lo) || is.infinite(lo)) {
    return(glue::glue("<= {criterion_num(hi)}"))
  }
  if (is.na(hi) || is.infinite(hi)) {
    return(glue::glue(">= {criterion_num(lo)}"))
  }
  if (hi == lo) {
    return(criterion_num(lo))
  }
  glue::glue("[{criterion_num(lo)}, {criterion_num(hi)}]")
}

# The text a criterion row shows: its `limit` in the pipeline's vocabulary.  A
# categorical bound that is its own constraint (e.g. "TTTV", "REF=1, ALT=0")
# shows as written; a boolean bound is the predicate the evaluator tests.
criterion_limit_text <- function(kind, bounds) {
  if (is.null(bounds) || !length(bounds) || all(is.na(bounds))) {
    return(NA_character_)
  }
  switch(kind,
    numeric = criterion_bound_text(bounds),
    boolean = as.character(bounds[[1L]]),
    categorical = {
      b <- as.character(bounds)
      if (length(b) == 1L) b else paste(b, collapse = " / ")
    }
  )
}

# Registry rows for a set of criteria that share their placement (one stage, one
# role, one level): `ids` are the stable ids, `values` the constant's entries
# (names + values, or the declared prose), `bounds` the thresholds.  A `bounds`
# element may be NULL, in which case the declared prose is the limit.
#
# `calibrated` is deprecated (ZW4): pass `evidence` instead.  It is still
# accepted and still emitted (so existing readers keep working), derived from
# the rung by the migration rule; `lifecycle` is not a dependency of this
# package, so the old argument works silently rather than warning.
criterion_tbl <- function(id, scope, stage, role, label, level, kind, bounds,
                          units, source, calibrated = NULL, evaluator,
                          values = NULL, evidence = NULL) {
  n <- length(id)
  stopifnot(all(kind %in% CRITERION_KINDS))
  if (is.null(evidence)) {
    rung <- if (is.null(calibrated)) {
      criterion_evidence(rep(NA_character_, n))
    } else {
      criterion_evidence_from_flag(calibrated)
    }
  } else {
    rung <- criterion_evidence(rep(evidence, length.out = n))
  }
  tibble::tibble(
    id = id,
    scope = rep(scope, length.out = n),
    stage = rep(stage, length.out = n),
    role = rep(role, length.out = n),
    label = label,
    level = rep(level, length.out = n),
    kind = kind,
    bounds = bounds,
    units = units,
    # the reference point a tier ranks toward; NA on every other row (ZW4.6)
    target = rep(NA_character_, n),
    source = rep(source, length.out = n),
    calibrated = rep(as.logical(calibrated %||% NA), length.out = n),
    evidence = rung,
    evaluator = evaluator,
    hard = FALSE,
    limit = mapply(criterion_limit_text, kind, bounds, USE.NAMES = FALSE) |>
      unname()
  ) |>
    dplyr::mutate(limit = ifelse(
      is.na(.data$limit),
      as.character(values %||% rep(NA_character_, n)),
      .data$limit
    ))
}

# The KASP-family acceptance checks, in registry order, with the check class each
# reports under (the family that owns it selects it through KASP_FAMILY_CHECKS).
#
# The class is DERIVED from the three `KASP_CHECK_*` label vectors, which are the
# single declaration of the partition (Z2.1/ZA6 step 1).  Before this, the same
# class was written twice - as ids here and as labels there - and a criterion
# added to only one of the two silently got `scope = NA` or "no scheme".  Now a
# row in no class, and a class label with no row, are both errors.
marker_check_classes <- function() {
  ways <- marker_check_ways()
  klass <- rep(NA_character_, nrow(ways))
  klass[ways$label %in% KASP_CHECK_PRIMER] <- "primer"
  klass[ways$label %in% KASP_CHECK_SITE] <- "site"
  klass[ways$label %in% KASP_CHECK_AS] <- "as"
  stray <- setdiff(c(KASP_CHECK_PRIMER, KASP_CHECK_SITE, KASP_CHECK_AS), ways$label)
  if (length(stray)) {
    stop("KASP_CHECK_* labels with no criterion row: ",
      paste(stray, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(klass)) {
    stop("criterion row(s) in no KASP_CHECK_* class: ",
      paste(ways$id[is.na(klass)], collapse = ", "),
      call. = FALSE
    )
  }
  tibble::tibble(id = ways$id, class = klass)
}

# The KASP-family acceptance checks: id, label, reporting level, kind, the
# literal bound when it is not a numeric threshold, its units, where the
# threshold comes from, and the evaluator.  The numeric bars are filled from the
# limit object in force by criterion_spec(), so a caller's `limits=` reaches the
# registry (and every view of it) unchanged.
#
# `affects` names, for a bar the check does not decide on its own, the other
# check it defers to: the footprint-variant counts are only judged when the
# footprint's dropout is unknown (a non-NA dropout is already the finer measure),
# so their bound restates that dropout bar and a test keeps the two equal.
marker_check_ways <- function() {
  tibble::tribble(
    ~id, ~label, ~level, ~kind, ~bounds, ~units, ~source, ~evaluator, ~affects,
    "site_class", "target site class", "site", "boolean", "discriminable", NA,
    "package assumption", "site_disc", NA,
    "site_alt_alleles", "ALT alleles (n)", "site", "numeric", NULL, "count",
    "package assumption", "site_n_alt", NA,
    "site_linked_skipped", "linked variants skipped (n)", "site", "numeric", NULL, "count",
    "package assumption", "site_linked_skipped", NA,
    "site_alt_covered", "ALT primer covers all ALTs", "site", "numeric", NULL, "count",
    "package assumption", "site_alt_covered", NA,
    "site_linked_dependent", "design relies on linkage", "site", "boolean", "= 0", NA,
    "package assumption", "site_linked_dependent", NA,
    "amplicon", "amplicon length (bp)", "alignment", "numeric", NULL, "bp",
    "package assumption", "assay_amplicon", NA,
    "as_tm", "AS Tm (C)", "feature", "numeric", NULL, "C",
    "package assumption", "assay_as_tm", NA,
    "com_tm", "COM Tm (C)", "feature", "numeric", NULL, "C",
    "package assumption", "assay_com_tm", NA,
    "com_tm_ge_as_tm", "COM Tm >= AS Tm", "alignment", "numeric", NULL, "C",
    "package assumption", "assay_com_tm_ge_as_tm", NA,
    "as_gc", "AS GC", "feature", "numeric", NULL, "fraction",
    "package assumption", "assay_as_gc", NA,
    "com_gc", "COM GC", "feature", "numeric", NULL, "fraction",
    "package assumption", "assay_com_gc", NA,
    "max_run", "longest homopolymer (bp)", "feature", "numeric", NULL, "bp",
    "package assumption", "longest_run", NA,
    "ssr", "SSR (2-4 nt, >=3 copies)", "feature", "numeric", NULL, "count",
    "package assumption", "has_repeat", NA,
    "inverted_repeat", "inverted repeat (arms>=4, loop<=6)", "feature", "numeric", NULL, "count",
    "package assumption", "has_palindrome", NA,
    "hairpin_dg", "worst hairpin dG37", "feature", "numeric", NULL, "dG37 (kcal/mol)",
    "package assumption", "min_hairpin_dg", NA,
    "cross_dimer_dg", "worst cross-dimer dG37", "feature", "numeric", NULL, "dG37 (kcal/mol)",
    "package assumption", "min_cross_dimer_dg", NA,
    "homodimer_dg", "worst homodimer dG37", "feature", "numeric", NULL, "dG37 (kcal/mol)",
    "package assumption", "min_homodimer_dg", NA,
    "as_n_3p", "AS 3' window variants (n)", "feature", "numeric", NULL, "count",
    "package assumption", "assay_as_n_3p", "as_dropout",
    "as_n_var", "AS footprint variants (n)", "feature", "numeric", NULL, "count",
    "package assumption", "assay_as_n_var", "as_foot_dropout",
    "com_n_var", "COM footprint variants (n)", "feature", "numeric", NULL, "count",
    "package assumption", "assay_com_n_var", "com_dropout",
    "as_dropout", "AS 3' window dropout", "feature", "numeric", NULL, "fraction",
    "package assumption", "assay_as_dropout", NA,
    "as_foot_dropout", "AS footprint dropout", "feature", "numeric", NULL, "fraction",
    "package assumption", "assay_as_foot_dropout", NA,
    "com_dropout", "COM footprint dropout", "feature", "numeric", NULL, "fraction",
    "package assumption", "assay_com_dropout", NA,
    "discrimination", "discrimination", "alignment", "numeric", NULL, "fraction",
    "package assumption", "assay_discrimination", NA,
    "as_mismatch_type", "AS 3' mismatch base pair", "feature", "categorical",
    c("strong", "medium", "weak"), NA, "primary literature", "as_mismatch_type", NA,
    "in_silico_pcr", "in-silico PCR (REF 1 / ALT 0)", "alignment", "boolean",
    "REF=1, ALT=0", NA, "package assumption", "in_silico_pcr", NA,
    "offtarget_hits", "off-target exact primer hits (n)", "specificity", "numeric", NULL,
    "count", "package assumption", "ideogram_hits", NA,
    "in_silico_product_count", "in-silico PCR products (n)", "specificity", "numeric", NULL,
    "count", "package assumption", "assay_amplicons", NA
  )
}

# Where each criterion's bar comes from, for the rows whose threshold a
# published source actually STATES.  A criterion is absent when no source states
# its number, and keeps the `source` its constructor gave it
# ("package assumption") - naming a vendor default that merely resembles the
# bar would turn a rule of thumb into a false citation (ZW4.6).
#
# Keys are criterion ids.  This is the citable field of `criterion_spec()`:
# it names the key of the row in inst/extdata/scientific_sources.csv a reader
# turns to, never the threshold itself (bounds stay in `bounds`).
criterion_sources <- function() {
  c(
    # The mutation-genotyping designs pick a 3' base by Kwok's mismatch rank: a
    # weak terminal mismatch is what makes the allele-specific primer specific.
    "as:gate:gate_wobble" = "kwok1990",
    "as:gate:gate_discrimination" = "kwok1990",
    "arms:gate:gate_arms_disc" = "kwok1990",
    "castpcr:gate:gate_castpcr_as3" = "kwok1990",
    "castpcr:gate:gate_castpcr_mismatch" = "kwok1990",
    "castpcr:gate:gate_castpcr_blocker" = "kwok1990",
    # `as_mismatch_type` names its source at its own row in marker_check_ways();
    # repeating it here keeps the map a complete view of the sourced rows.
    as_mismatch_type = "kwok1990",
    # The dG *values* use SantaLucia 1998's nearest-neighbour parameters, but
    # the *threshold* a design has to clear is printed in the papers below:
    # DeGenPrime filters a combined hairpin/dimer dG below -3 kcal/mol (the
    # same number as each *_structure bar), and an allele-specific qPCR protocol
    # keeps primers/probes with a secondary-structure dG >= -9 kcal/mol (the
    # cross_dimer_dg / homodimer_dg bars).
    hairpin_dg = "degenprime2024",
    cross_dimer_dg = "vavilov2024",
    homodimer_dg = "vavilov2024",
    "pair:gate:gate_structure" = "degenprime2024",
    "indel_size:gate:gate_indel_size_structure" = "degenprime2024",
    "dpcr:gate:gate_dpcr_probe_structure" = "degenprime2024",
    "taqman:gate:gate_taqman_probe_structure" = "degenprime2024",
    "snapshot:gate:gate_snapshot_structure" = "degenprime2024",
    "rpa:gate:gate_rpa_structure" = "degenprime2024",
    # The Barley ASQ method paper states its own amplicon (50-100 bp) and the
    # forward primers' Tm window (55-56 C, p. 7) for these two rows.
    amplicon = "fret_asqpcr2021",
    # ZW4.6 A/B: the 63 C ASQ optimum (range 62-64 C) is printed here, not in
    # the FRET method paper, which states the amplicon window only.
    "asq:tier:tier1" = "kalendar2022_competitive",
    # Lobato & O'Sullivan 2018 (a review) states the RPA primer-length,
    # amplicon and nfo-probe geometry bars the RPA scheme inherited.  The
    # primer-length sentence is the historical belief the review then softens
    # (normal PCR primers also work), so the row stays a rule of thumb.
    "rpa:gate:gate_rpa_primer_len" = "lobato2018",
    family_rpa_primer_length = "lobato2018",
    "rpa:gate:gate_rpa_amplicon" = "lobato2018",
    "rpa:gate:gate_rpa_probe" = "lobato2018",
    # ZW4.6 round 4.  The FAASA protocol (a peer-reviewed Bio-protocol) prints
    # the 40-60 % GC window its allele-specific primers are designed to.
    as_gc = "faasa2026",
    com_gc = "faasa2026",
    # EasyKASP states the KASP amplicon range as a typical design property
    # (50-100 bp); see its `access_detail` for the "typically" caveat.
    "pair:gate:gate_amplicon" = "easykasp2025",
    # Ichida 2018 prints the 50-100 bp product window its HRM primers were
    # designed to produce.
    "hrm:gate:gate_hrm_amplicon" = "jsbbs2018",
    # De Keyzer 2026 applies the Piepenburg/TwistDx manual's homopolymer rule
    # (no more than three identical nucleotides = run < 4) to its RPA primers.
    max_run = "dekeyzer2026",
    # A Cas12a editing protocol prints the canonical TTTV PAM (the family row's
    # whole bar).  The combined per-scheme row (Cas12a + Cas13a + Cas13b) is
    # NOT sourced: no single text prints all three tokens.
    family_crispr_pam = "cas12a_beech2026",
    # ZW4.6 round 4, relaxed same-quantity rule: a source that prints a design
    # rule for the same quantity qualifies even when its range is wider or
    # offset; the printed number and the offset live in the source's
    # `access_detail`, and the package `bounds` never move.
    #
    # FAASA (57-60 C for each primer / TSP-CGSP within 1-2 C).
    as_tm = "faasa2026",
    com_tm = "faasa2026",
    # EasyKASP common-primer Tm +2 C above the specific primers.
    com_tm_ge_as_tm = "easykasp2025",
    "pair:gate:gate_tm" = "easykasp2025",
    # mxKASP targets two non-homologous markers in one reaction (the dKASP
    # homologous branch is stated separately in the literature).
    "duplex:gate:gate_duplex_targets" = "mxkasp2025",
    # LAMP review: F2/B2 arm 60-65 C; target <=300 nt, typically 130-260;
    # no secondary structure.
    "lamp:gate:gate_lamp_inner_tm" = "nemeth2025",
    "lamp:gate:gate_lamp_amplicon" = "nemeth2025",
    "lamp:gate:gate_lamp_structure" = "nemeth2025",
    # SNP-LAMP puts the allele-specific base at the 3' terminus of an inner
    # primer (B2/BIP) - the same principle as the package's FIP rule.
    "lamp:gate:gate_lamp_variant" = "snp_lamp_lfd2025",
    # Cas12a crRNA design: 19-24 nt spacer and 30-70 % GC.
    "crispr:gate:gate_crispr_spacer" = "cas12a_beech2026",
    "crispr:gate:gate_crispr_gc" = "cas12a_beech2026",
    # ddPCR: amplicon 70-150 bp, primer Tm 58-62, probe Tm +5-10, probe
    # GC 40-60 %.
    "dpcr:gate:gate_dpcr_probe_comp" = "ppr_ddpcr2026",
    "dpcr:gate:gate_dpcr_primer_tm" = "ppr_ddpcr2026",
    "dpcr:gate:gate_dpcr_amplicon" = "ppr_ddpcr2026",
    "dpcr:gate:gate_dpcr_probe_tm" = "ppr_ddpcr2026",
    # TaqMan allele-specific qPCR: probe GC 30-70, no 5' G, no 4+ G run;
    # primer Tm 58-60; amplicon 80-120; probe +6-8 C.
    "taqman:gate:gate_taqman_probe_comp" = "vavilov2024",
    "taqman:gate:gate_taqman_primer_tm" = "vavilov2024",
    "taqman:gate:gate_taqman_amplicon" = "vavilov2024",
    # castPCR is TaqMan-probe allele-specific chemistry, so the same qPCR
    # criteria apply to its ASP, probe and amplicon.
    "castpcr:gate:gate_castpcr_as_tm" = "vavilov2024",
    "castpcr:gate:gate_castpcr_probe_tm" = "vavilov2024",
    "castpcr:gate:gate_castpcr_amplicon" = "vavilov2024",
    family_probe_tm_delta = "vavilov2024",
    # Sanger read margin: >=50 nt from the sequencing primer to the variant.
    "sanger:gate:gate_sanger_dist" = "vavilov2024",
    family_sanger_read_distance = "vavilov2024",
    # RPA: a mismatch at the primer 3' end suppresses amplification.
    "rpa:gate:gate_rpa_variant" = "lobato2018",
    # Avoid inverted repeats in a primer.
    inverted_repeat = "kalendar2025",
    # Gel band resolution >=10 bp; product-size delta >=10 bp; SSR amplicon
    # 50-150 bp.
    "caps:gate:gate_caps_resolvable" = "ssr_gel2026",
    "arms:gate:gate_arms_bands" = "ssr_gel2026",
    family_caps_band_gap = "ssr_gel2026",
    family_size_delta = "ssr_gel2026",
    "ssr:gate:gate_ssr_amplicon" = "ssr_gel2026",
    # SSR repeat-copy floor (motif-dependent; trinucleotide floor is 5).
    "ssr:gate:gate_ssr_polymorphic" = "moso_bamboo2026",
    # InDel amplicon window printed as 200-320 bp.
    "indel_size:gate:gate_indel_size_amplicon" = "indel_melon2025",
    "indel_size:gate:gate_indel_size_gap" = "ssr_gel2026",
    # General PCR-primer Tm window (58-60 C) printed under "Design of primers"
    # in an allele-specific qPCR protocol; the same quantity as each scheme's
    # primer-Tm bar.
    "hrm:gate:gate_hrm_tm" = "vavilov2024",
    "indel_size:gate:gate_indel_size_tm" = "vavilov2024",
    "ssr:gate:gate_ssr_tm" = "vavilov2024",
    "caps:gate:gate_caps_tm" = "vavilov2024",
    "sanger:gate:gate_sanger_tm" = "vavilov2024",
    "duplex:gate:gate_duplex_tm" = "vavilov2024",
    # Repeat-copy floor (motif-dependent; trinucleotide floor is 5).
    family_size_repeat_copies = "moso_bamboo2026",
    # Standard TaqMan guideline (preferred probe Tm 68-70 C) restated by a
    # paper whose own probe fell outside it.
    "taqman:gate:gate_taqman_probe_tm" = "cov_a4_qpcr2026",
    # ZW4.6 round 4 continued: multiplex KASP prints the exact primer-dimer
    # and self-complementarity bars; the TaqMan protocol prints the probe
    # variant-position rule; two more amplicon/probe windows are printed
    # verbatim.
    "duplex:gate:gate_duplex_complement" = "wang2025_multikasp",
    "duplex:gate:gate_duplex_selfcomp" = "wang2025_multikasp",
    "caps:gate:gate_caps_amplicon" = "liu2016_caps",
    "taqman:gate:gate_taqman_centrality" = "vavilov2024",
    # The SNaPshot extension primer inherits the same composition rule as a
    # probe (GC floor, homopolymer ceiling, no long SSR): a forensic InDel
    # paper prints it for its probes.
    "snapshot:gate:gate_snapshot_comp" = "dekeyzer2026",
    # The RPA amplicon window the CRISPR scheme inherits (80-250 bp) is the
    # RPA guide's 100-200 bp window, stated for the same quantity.
    "crispr:gate:gate_crispr_amplicon" = "lobato2018",
    # SNaPshot/SBE design: a protocol places the extension primer one base
    # upstream of the SNP; a 14-plex cSNP SNaPshot assay prints the extension
    # primer Tm window. An InDel-marker paper prints the single-locus rule.
    "snapshot:gate:gate_snapshot_placement" = "tso2022_snapshot",
    "snapshot:gate:gate_snapshot_tm" = "zhao2026_snapshot",
    "indel_size:gate:gate_indel_size_unique" = "li2026_indel",
    # A probe must bind only one region; the same protocol prints the 3'-end
    # homopolymer / secondary-structure primer rule.
    "taqman:gate:gate_taqman_probe_specific" = "vavilov2024",
    "sanger:gate:gate_sanger_primer" = "faasa2026",
    # CAPS: the polymorphic enzyme cut sits in the middle of the amplicon;
    # SSR: one in-silico product and primers from unique flanking sequence;
    # SNP-LAMP: the variant sits at the FIP 3' end (or BIP 5' end).
    "caps:gate:gate_caps_site" = "liu2016_caps",
    "ssr:gate:gate_ssr_offtarget" = "zheng2025_ssr",
    "ssr:gate:gate_ssr_flanks" = "zheng2025_ssr",
    # CAPS: the digest fragments of an allele sum to its undigested product
    # (printed in the Fig. 2 legend of the marker-development paper).
    family_caps_digest_spans = "liu2016_caps",
    # CAPS: VCF2CAPS keeps only markers with a single restriction site - the
    # target, polymorphic one - inside the amplified sequence.
    family_caps_site_conserved = "vcf2caps2021",
    family_lamp_fip_alt = "nemeth2025",
    family_lamp_fip_variant = "nemeth2025",
    # CAPS: the polymorphic cut site is defined within the amplicon and
    # separates the alleles (distinct bands after digestion).  LAMP: FIP/BIP
    # 3' ends lie inside the F3/B3 footprint.  TaqMan: the probe spans the
    # variant (site in the middle third of an 18-30 nt probe).
    family_caps_site_defined = "liu2016_caps",
    family_caps_site_separates = "liu2016_caps",
    family_lamp_inner_inside = "nemeth2025",
    family_probe_spans = "vavilov2024",
    # generic primer QC: avoid repeats inside the primer, and require no
    # non-specific genome binding (the InDel paper drops such primers).
    ssr = "faasa2026",
    offtarget_hits = "li2026_indel",
    # SSR repeat-unit length: the sugar-beet marker paper prints 1-6 bp
    # (the package's 2-4 nt is the subset it designs for).
    family_size_repeat_unit = "ssr_gel2026",
    # ZW4.6 phases A/B: the TARGET half of each ranking tier.  A tier's `limit`
    # is its aggregation policy (a package decision, `rule_of_thumb` by
    # construction); these keys name the source that prints the reference point
    # the tier ranks toward.  Printed number and any offset are in `target` and
    # the source's `access_detail`.
    "as:tier:tier2" = "vavilov2024",
    "pair:tier:tier1" = "degenprime2024",
    "pair:tier:tier2" = "easykasp2025",
    "pair:tier:tier3" = "easykasp2025",
    "duplex:tier:tier1" = "wang2025_multikasp",
    "duplex:tier:tier2" = "wang2025_multikasp",
    "caps:tier:tier1" = "ssr_gel2026",
    "caps:tier:tier2" = "vavilov2024",
    "caps:tier:tier3" = "liu2016_caps",
    "arms:tier:tier2" = "ssr_gel2026",
    "arms:tier:tier3" = "vavilov2024",
    "hrm:tier:tier1" = "jsbbs2018",
    "hrm:tier:tier2" = "vavilov2024",
    "castpcr:tier:tier2" = "vavilov2024",
    "castpcr:tier:tier3" = "vavilov2024",
    "crispr:tier:guide" = "cas12a_beech2026",
    "dpcr:tier:tier1" = "ppr_ddpcr2026",
    "dpcr:tier:tier2" = "degenprime2024",
    "dpcr:tier:tier3" = "ppr_ddpcr2026",
    "indel_size:tier:tier1" = "vavilov2024",
    "indel_size:tier:tier3" = "vavilov2024",
    "lamp:tier:tier2" = "nemeth2025",
    "lamp:tier:tier3" = "nemeth2025",
    "lamp:tier:tier4" = "nemeth2025",
    "rpa:tier:tier2" = "lobato2018",
    "rpa:tier:tier3" = "lobato2018",
    "rpa:tier:tier4" = "lobato2018",
    "sanger:tier:tier2" = "vavilov2024",
    "snapshot:tier:tier1" = "zhao2026_snapshot",
    "snapshot:tier:tier2" = "degenprime2024",
    "ssr:tier:tier1" = "vavilov2024",
    "ssr:tier:tier2" = "moso_bamboo2026",
    "ssr:tier:tier3" = "vavilov2024",
    "taqman:tier:tier2" = "cov_a4_qpcr2026",
    "taqman:tier:tier3" = "ppr_ddpcr2026"
  )
}

# The reference point each ranking tier ranks candidates toward (ZW4.6).  This
# is the citable half of the row; the aggregation policy - the bottleneck `min`,
# tier membership and order - is the package's, is declared once by
# `CRITERION_TIER_LIMIT` / `CRITERION_TIER_POLICY_EVIDENCE`, and is never
# disguised as a threshold.  A target is present for every tier, sourced or not;
# when a source states it, `criterion_sources()` names the key and
# `criterion_evidence_sources()` sets the rung.
#
# Keys are criterion ids.  Printed numbers are given in the source's own units;
# where the package bar is offset or narrower, the offset is repeated in the
# source's `access_detail`.
