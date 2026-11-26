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
criterion_tier_targets <- function() {
  c(
    "as:tier:tier1" = "mutant/WT discrimination margin (computed in-silico)",
    "as:tier:tier2" = "AS primer Tm 58-60 C (package bar 58-63 C)",
    "as:tier:tier3" = "3'-terminal mismatch weakness (Kwok rank)",
    "as:tier:tier4" = "footprint load (simulated dropout)",
    "pair:tier:tier1" = "hairpin + homodimer dG >= -3 kcal/mol",
    "pair:tier:tier2" = "common-primer Tm >= AS Tm, ideally +2 C",
    "pair:tier:tier3" = "amplicon 50-100 bp",
    "pair:tier:tier4" = "footprint load (simulated dropout)",
    "pair:tier:tier5" = "weakest COM 3'-terminal 5-mer dG (package rule)",
    "duplex:tier:tier1" = "<= 5 self-complementary bases",
    "duplex:tier:tier2" = "shared amplicon <= 150 bp",
    "asq:tier:tier1" = "AS primer Tm 63 C (range 62-64 C)",
    "caps:tier:tier1" = ">= 10 bp band gap (gel resolvability)",
    "caps:tier:tier2" = "primer pair dTm <= 2 C, Tm 58-60 C",
    "caps:tier:tier3" = "CAPS product 200-1000 bp",
    "arms:tier:tier1" = "inner-primer discrimination margin (Kwok rank)",
    "arms:tier:tier2" = ">= 10 bp band separation",
    "arms:tier:tier3" = "primer Tm spread <= 2 C across the tube",
    "hrm:tier:tier1" = "variant at the centre of the amplified sequence",
    "hrm:tier:tier2" = "primer pair dTm <= 2 C",
    "castpcr:tier:tier1" = "ASP mutant/WT discrimination margin (Kwok rank)",
    "castpcr:tier:tier2" = "AS/probe Tm 58-60 C, dTm <= 2 C",
    "castpcr:tier:tier3" = "amplicon 80-120 nt",
    "crispr:tier:guide" = "Cas12a TTTV PAM; spacer 19-24 nt; GC 30-70 %",
    "crispr:tier:seed" = "variant inside the PAM-proximal seed window",
    "crispr:tier:specificity" = "guide off-target count (computed)",
    "crispr:tier:amplicon" = "pre-amplification product window (scheme default)",
    "dpcr:tier:tier1" = "probe Tm 5-10 C above the primer Tm",
    "dpcr:tier:tier2" = "probe hairpin/homodimer dG >= -3 kcal/mol",
    "dpcr:tier:tier3" = "amplicon 70-150 bp",
    "indel_size:tier:tier1" = "primer pair dTm <= 2 C, Tm 58-60 C",
    "indel_size:tier:tier2" = "amplicon wide enough to resolve the InDel (scheme default)",
    "indel_size:tier:tier3" = "<= 2 G/C in the last 5 nt; no T at the 3' end",
    "lamp:tier:tier1" = "variant at the FIP 3' end / BIP 5' end (mismatch position)",
    "lamp:tier:tier2" = "F2/B2 arm Tm 60-65 C",
    "lamp:tier:tier3" = "no secondary structure in the primers",
    "lamp:tier:tier4" = "target 130-260 nt, at most 300 nt",
    "lax:tier:tier1" = "resolvable ALT haplotype frequency (computed)",
    "rpa:tier:tier1" = "3'-anchor mismatch class (classification)",
    "rpa:tier:tier2" = "primer length 30-35 nt",
    "rpa:tier:tier3" = "THF probe: >= 30 nt from the 5' end, >= 15 nt from the 3' end",
    "rpa:tier:tier4" = "amplicon 100-200 bp",
    "sanger:tier:tier1" = "variant central in the amplicon",
    "sanger:tier:tier2" = "primer pair dTm <= 2 C",
    "sanger:tier:tier3" = "amplicon length (scheme default)",
    "snapshot:tier:tier1" = "extension-primer Tm 53-60 C",
    "snapshot:tier:tier2" = "extension primer hairpin/homodimer dG >= -3 kcal/mol",
    "snapshot:tier:tier3" = "amplicon near the ideal (scheme default)",
    "ssr:tier:tier1" = "primer pair dTm <= 2 C, Tm 58-60 C",
    "ssr:tier:tier2" = "amplicon 100-300 bp",
    "ssr:tier:tier3" = "<= 2 G/C in the last 5 nt; no T at the 3' end",
    "taqman:tier:tier1" = "probe-variant mismatch dG margin (package rule)",
    "taqman:tier:tier2" = "probe Tm 68-70 C, 8-10 C above the primers",
    "taqman:tier:tier3" = "amplicon 70-150 bp"
  )
}

# The rung behind each raised row.  `criterion_sources()` says WHERE the number
# is printed; this says HOW STRONGLY it is fixed, and only ever raises a row that
# has a source (a rung with no citation would be provenance without a page).
# Rows absent here keep the constructor's `rule_of_thumb`.
criterion_evidence_sources <- function() {
  c(
    as_gc = "primary_literature",
    com_gc = "primary_literature",
    "pair:gate:gate_amplicon" = "primary_literature",
    "hrm:gate:gate_hrm_amplicon" = "primary_literature",
    # the number is a vendor (TwistDx/Piepenburg) default reproduced in a paper
    max_run = "vendor_guideline",
    family_crispr_pam = "primary_literature",
    # ZW4.6 round 4 expansion: a peer-reviewed paper or protocol prints the
    # same quantity as the bar (see each source's `access_detail` for the
    # printed number and its offset from the package `bounds`).
    as_tm = "primary_literature",
    com_tm = "primary_literature",
    com_tm_ge_as_tm = "primary_literature",
    "pair:gate:gate_tm" = "primary_literature",
    "duplex:gate:gate_duplex_targets" = "primary_literature",
    "lamp:gate:gate_lamp_inner_tm" = "primary_literature",
    "lamp:gate:gate_lamp_amplicon" = "primary_literature",
    "lamp:gate:gate_lamp_structure" = "primary_literature",
    "lamp:gate:gate_lamp_variant" = "primary_literature",
    "crispr:gate:gate_crispr_spacer" = "primary_literature",
    "crispr:gate:gate_crispr_gc" = "primary_literature",
    "dpcr:gate:gate_dpcr_probe_comp" = "primary_literature",
    "dpcr:gate:gate_dpcr_primer_tm" = "primary_literature",
    "dpcr:gate:gate_dpcr_amplicon" = "primary_literature",
    "dpcr:gate:gate_dpcr_probe_tm" = "primary_literature",
    "taqman:gate:gate_taqman_probe_comp" = "primary_literature",
    "taqman:gate:gate_taqman_primer_tm" = "primary_literature",
    "taqman:gate:gate_taqman_amplicon" = "primary_literature",
    "castpcr:gate:gate_castpcr_as_tm" = "primary_literature",
    "castpcr:gate:gate_castpcr_probe_tm" = "primary_literature",
    "castpcr:gate:gate_castpcr_amplicon" = "primary_literature",
    family_probe_tm_delta = "primary_literature",
    "sanger:gate:gate_sanger_dist" = "primary_literature",
    family_sanger_read_distance = "primary_literature",
    "rpa:gate:gate_rpa_variant" = "primary_literature",
    inverted_repeat = "primary_literature",
    "caps:gate:gate_caps_resolvable" = "primary_literature",
    "arms:gate:gate_arms_bands" = "primary_literature",
    family_caps_band_gap = "primary_literature",
    family_size_delta = "primary_literature",
    "ssr:gate:gate_ssr_amplicon" = "primary_literature",
    "ssr:gate:gate_ssr_polymorphic" = "primary_literature",
    "indel_size:gate:gate_indel_size_amplicon" = "primary_literature",
    "indel_size:gate:gate_indel_size_gap" = "primary_literature",
    "hrm:gate:gate_hrm_tm" = "primary_literature",
    "indel_size:gate:gate_indel_size_tm" = "primary_literature",
    "ssr:gate:gate_ssr_tm" = "primary_literature",
    "caps:gate:gate_caps_tm" = "primary_literature",
    "sanger:gate:gate_sanger_tm" = "primary_literature",
    "duplex:gate:gate_duplex_tm" = "primary_literature",
    family_size_repeat_copies = "primary_literature",
    # a vendor default restated inside a paper
    "taqman:gate:gate_taqman_probe_tm" = "vendor_guideline",
    # ZW4.6 round 4 continued: the secondary-structure thresholds and the
    # primer-dimer / self-complementarity bars are all printed (see each
    # source's `access_detail` for the number and its offset).
    hairpin_dg = "primary_literature",
    cross_dimer_dg = "primary_literature",
    homodimer_dg = "primary_literature",
    "pair:gate:gate_structure" = "primary_literature",
    "indel_size:gate:gate_indel_size_structure" = "primary_literature",
    "dpcr:gate:gate_dpcr_probe_structure" = "primary_literature",
    "taqman:gate:gate_taqman_probe_structure" = "primary_literature",
    "snapshot:gate:gate_snapshot_structure" = "primary_literature",
    "rpa:gate:gate_rpa_structure" = "primary_literature",
    "duplex:gate:gate_duplex_complement" = "primary_literature",
    "duplex:gate:gate_duplex_selfcomp" = "primary_literature",
    "caps:gate:gate_caps_amplicon" = "primary_literature",
    "taqman:gate:gate_taqman_centrality" = "primary_literature",
    "rpa:gate:gate_rpa_amplicon" = "primary_literature",
    "rpa:gate:gate_rpa_probe" = "primary_literature",
    amplicon = "primary_literature",
    "snapshot:gate:gate_snapshot_comp" = "primary_literature",
    "crispr:gate:gate_crispr_amplicon" = "primary_literature",
    "snapshot:gate:gate_snapshot_placement" = "primary_literature",
    "snapshot:gate:gate_snapshot_tm" = "primary_literature",
    "indel_size:gate:gate_indel_size_unique" = "primary_literature",
    "taqman:gate:gate_taqman_probe_specific" = "primary_literature",
    "sanger:gate:gate_sanger_primer" = "primary_literature",
    "caps:gate:gate_caps_site" = "primary_literature",
    "ssr:gate:gate_ssr_offtarget" = "primary_literature",
    "ssr:gate:gate_ssr_flanks" = "primary_literature",
    family_lamp_fip_alt = "primary_literature",
    family_lamp_fip_variant = "primary_literature",
    family_caps_site_defined = "primary_literature",
    family_caps_site_separates = "primary_literature",
    family_caps_digest_spans = "primary_literature",
    family_caps_site_conserved = "primary_literature",
    family_lamp_inner_inside = "primary_literature",
    family_probe_spans = "primary_literature",
    # the review prints the 30-35 nt rule before noting PCR primers also work
    "rpa:gate:gate_rpa_primer_len" = "primary_literature",
    family_rpa_primer_length = "primary_literature",
    ssr = "primary_literature",
    offtarget_hits = "primary_literature",
    family_size_repeat_unit = "primary_literature",
    # ZW4.6 phases A/B: the target's rung.  Every tier above has a printed
    # reference point; the aggregation policy stays rule_of_thumb by
    # construction (CRITERION_TIER_POLICY_EVIDENCE), never promoted here.
    "as:tier:tier2" = "primary_literature",
    "pair:tier:tier1" = "primary_literature",
    "pair:tier:tier2" = "primary_literature",
    "pair:tier:tier3" = "primary_literature",
    "duplex:tier:tier1" = "primary_literature",
    "duplex:tier:tier2" = "primary_literature",
    "asq:tier:tier1" = "primary_literature",
    "caps:tier:tier1" = "primary_literature",
    "caps:tier:tier2" = "primary_literature",
    "caps:tier:tier3" = "primary_literature",
    "arms:tier:tier2" = "primary_literature",
    "arms:tier:tier3" = "primary_literature",
    "hrm:tier:tier1" = "primary_literature",
    "hrm:tier:tier2" = "primary_literature",
    "castpcr:tier:tier2" = "primary_literature",
    "castpcr:tier:tier3" = "primary_literature",
    "crispr:tier:guide" = "primary_literature",
    "dpcr:tier:tier1" = "primary_literature",
    "dpcr:tier:tier2" = "primary_literature",
    "dpcr:tier:tier3" = "primary_literature",
    "indel_size:tier:tier1" = "primary_literature",
    "indel_size:tier:tier3" = "primary_literature",
    "lamp:tier:tier2" = "primary_literature",
    "lamp:tier:tier3" = "primary_literature",
    "lamp:tier:tier4" = "primary_literature",
    "rpa:tier:tier2" = "primary_literature",
    "rpa:tier:tier3" = "primary_literature",
    "rpa:tier:tier4" = "primary_literature",
    "sanger:tier:tier2" = "primary_literature",
    "snapshot:tier:tier1" = "primary_literature",
    "snapshot:tier:tier2" = "primary_literature",
    "ssr:tier:tier1" = "primary_literature",
    "ssr:tier:tier2" = "primary_literature",
    "ssr:tier:tier3" = "primary_literature",
    "taqman:tier:tier2" = "vendor_guideline",
    "taqman:tier:tier3" = "primary_literature"
  )
}

# Where a metamorphic invariant deliberately does NOT hold, and why.  A
# criterion is absent when the rule holds for it; a row that carries a reason
# is one whose two swapped inputs are not interchangeable, so asserting
# symmetry over it would be asserting a falsehood (ZW0.1).
#
# Keys are criterion ids, values one sentence a reviewer can check against the
# design code - the reason names the structural fact that breaks the symmetry,
# never a bare "exception".
criterion_asymmetry <- function() {
  c(
    # KASP: `AS_ref` is built on the FAM tail and `AS_alt` on the HEX tail
    # (R/marker_assemble.R), while both sit at the same AS concentration
    # (R/marker_mix.R).  The gates grade the chosen allele-specific primer, not
    # the dye it lands on, so swapping REF/ALT changes which allele a fluorophore
    # reports even though neither primer changed.
    "as:gate:gate_discrimination" = paste0(
      "KASP fluorophore asymmetry: AS_ref carries the FAM tail and AS_alt the HEX ",
      "tail (marker_assemble.R), so swapping REF/ALT moves an allele to another ",
      "fluorophore though the AS recipe is identical."
    ),
    "as:gate:gate_wobble" = paste0(
      "KASP fluorophore asymmetry: the wobble gate grades the chosen AS primer ",
      "once, while its allele is pinned to FAM or HEX by the tail, so the pair is ",
      "not mirrorable across REF/ALT."
    ),
    "as:gate:gate_as_footprint" = paste0(
      "KASP fluorophore asymmetry: the footprint gate is judged on AS_ref/AS_alt ",
      "as placed on the FAM/HEX tails, so a REF/ALT swap re-assigns the graded ",
      "primer to the other dye slot."
    ),
    "as:gate:gate_as_window" = paste0(
      "KASP fluorophore asymmetry: the AS window gate reports the primer in its ",
      "fixed FAM or HEX slot, so the REF/ALT image of the probe set is a different ",
      "probe set, not the same one relabelled."
    ),
    # LAMP: FIP is built with the ALT allele at the 3' arm and its 3' end anchors
    # BIP (R/marker_lamp.R); `lamp_fuse()` also requires each inner primer's 5'
    # loop arm to be hotter than its 3' displacement arm.  FIP and BIP are
    # therefore neither interchangeable with each other nor mirrorable.
    "family_lamp_fip_alt" = paste0(
      "LAMP FIP/BIP asymmetry: only FIP reads the variant, with the ALT base at ",
      "the 3' terminus of its displacement arm, so BIP cannot host the same check."
    ),
    "family_lamp_inner_tm_balance" = paste0(
      "LAMP FIP/BIP asymmetry: each inner primer must have a 5' loop arm hotter ",
      "than its 3' displacement arm (lamp_fuse), so the FIP and BIP Tm patterns are ",
      "constrained separately and cannot be exchanged."
    ),
    "family_lamp_fip_variant" = paste0(
      "LAMP FIP/BIP asymmetry: the variant is a fixed feature of FIP's 3' arm, ",
      "and BIP is placed downstream of FIP (lamp_bip from fip$arm3), so the two ",
      "inner primers are positionally and functionally distinct."
    ),
    "lamp:gate:gate_lamp_variant" = paste0(
      "LAMP FIP/BIP asymmetry: discrimination is scored at the FIP 3' terminus ",
      "(or a loop primer), never at BIP, which lies downstream of the variant."
    ),
    "lamp:gate:gate_lamp_inner_tm" = paste0(
      "LAMP FIP/BIP asymmetry: FIP's 3' base is the ALT allele so it cannot always ",
      "be a G/C clamp, making the inner-primer Tm gate deliberately weaker for FIP ",
      "than for BIP."
    )
  )
}

# The numeric KASP-family bars, read from the limit object in force.
marker_check_bounds <- function(limits) {
  tol <- limits$footprint_dropout %||% 0.05
  list(
    site_alt_alleles = c(1, NA),
    site_linked_skipped = c(NA, 0),
    site_alt_covered = c(NA, 0),
    amplicon = limits$amplicon,
    as_tm = limits$as_tm,
    com_tm = limits$com_tm,
    com_tm_ge_as_tm = c(0, NA),
    as_gc = limits$gc,
    com_gc = limits$gc,
    max_run = c(NA, "exclusive", limits$max_run),
    ssr = c(NA, 0),
    inverted_repeat = c(NA, 0),
    hairpin_dg = c(limits$hairpin_dg, NA),
    cross_dimer_dg = c(limits$duplex_dg, NA),
    homodimer_dg = c(limits$homodimer_dg, NA),
    as_n_3p = c(NA, 0),
    as_n_var = c(NA, 0),
    com_n_var = c(NA, 0),
    as_dropout = c(NA, tol),
    as_foot_dropout = c(NA, tol),
    com_dropout = c(NA, tol),
    discrimination = c(limits$discrimination, NA),
    # a categorical check: the allowed levels ARE the bound (Z0.4)
    as_mismatch_type = c("strong", "medium", "weak"),
    offtarget_hits = c(NA, 1),
    in_silico_product_count = c(NA, 1)
  )
}

# The family-specific checks: a label may serve several families (the probe and
# size checks), so the vocabulary lives here once rather than in a flat vector
# per family.  The `MARKER_CHECK_*` vectors are views of this table.
#
# `limit` is the published prose; `bounds` re-states the threshold the prose
# carries (`c(lo, hi)` with `Inf` for an open side, `NA` where the bar depends on
# the scheme's product window, or the literal level), and a test pins the two
# together, so a family check cannot show one bar and test another.
marker_family_check_ways <- function() {
  tibble::tribble(
    ~family, ~label, ~level, ~kind, ~bounds, ~units, ~limit, ~source, ~evaluator,
    "caps", "restriction site defined", "alignment", "categorical", "present", NA,
    "present", "package assumption", "caps_site_defined",
    "caps", "restriction site separates alleles", "alignment", "categorical",
    "REF != ALT", NA, "REF != ALT", "package assumption", "caps_site_separates",
    "caps", "band gap (bp)", "alignment", "numeric", "c(20, Inf)", "bp",
    ">= 20 bp", "package assumption", "caps_band_gap",
    "caps", "digest sums to amplicon", "feature", "categorical", "REF/ALT = product", NA,
    "REF/ALT = product", "package assumption", "caps_digest_spans",
    "caps", "restriction site conserved (population)", "specificity", "categorical",
    "yes", NA, "yes", "package assumption", "caps_site_conserved",
    "lamp", "FIP 3' ALT", "feature", "categorical", "ALT at FIP 3'", NA,
    "ALT at FIP 3'", "package assumption", "lamp_fip_alt",
    "lamp", "inner arm Tm balance (C)", "feature", "numeric", "c(NA, 5)", "C",
    "<= 5 C", "package assumption", "lamp_inner_tm_balance",
    "lamp", "inner arms inside amplicon", "alignment", "categorical", "inside", NA,
    "inside", "package assumption", "lamp_inner_inside",
    "lamp", "variant at FIP 3' end", "feature", "categorical", "yes", NA,
    "yes", "package assumption", "lamp_fip_variant",
    "rpa", "primer length (nt)", "feature", "numeric", "c(30, 35)", "count",
    "30-35 nt", "package assumption", "rpa_primer_length",
    "rpa", "THF off variant", "feature", "categorical", "yes", NA,
    "yes", "package assumption", "rpa_thf_off",
    "rpa", "RPA probe spans variant", "alignment", "categorical", "yes", NA,
    "yes", "package assumption", "rpa_probe_spans",
    "rpa", "variant under primer/probe", "feature", "categorical", "yes", NA,
    "yes", "package assumption", "rpa_variant_under",
    "crispr", "PAM valid", "feature", "categorical", "TTTV", NA,
    "TTTV", "package assumption", "crispr_pam",
    "crispr", "variant in seed", "feature", "categorical", "yes", NA,
    "yes", "package assumption", "crispr_seed",
    "crispr", "guide spans variant", "alignment", "categorical", "yes", NA,
    "yes", "package assumption", "crispr_guide_spans",
    "crispr", "guide/PAM conserved (population)", "specificity", "categorical", "yes", NA,
    "yes", "package assumption", "crispr_conserved",
    "crispr", "off-target protospacers (n)", "specificity", "numeric", "c(NA, 0)", "count",
    "0", "package assumption", "crispr_offtarget",
    "probe", "probe spans variant", "alignment", "categorical", "yes", NA,
    "yes", "package assumption", "probe_spans",
    "probe", "probe Tm - primer Tm (C)", "feature", "numeric", "c(5, Inf)", "C",
    ">= 5 C", "package assumption", "probe_tm_delta",
    "probe", "ASP 3' discrimination", "feature", "categorical", "> 0", NA,
    "> 0", "package assumption", "probe_asp_discrimination",
    "probe", "extension 3' at variant", "feature", "categorical", "adjacent", NA,
    "adjacent", "package assumption", "probe_extension",
    "sanger", "variant-read distance (bp)", "alignment", "numeric", "c(50, Inf)", "bp",
    ">= 50 bp", "package assumption", "sanger_read_distance",
    "sanger", "read primers inside amplicon", "alignment", "categorical", "inside", NA,
    "inside", "package assumption", "sanger_primers_inside",
    "sanger", "read primer pair defined", "feature", "categorical", "pair", NA,
    "pair", "package assumption", "sanger_pair_defined",
    "size", "repeat unit (nt)", "feature", "numeric", "c(2, 4)", "count",
    "2-4 nt", "package assumption", "size_repeat_unit",
    "size", "repeat copies", "feature", "numeric", "c(5, Inf)", "count",
    ">= 5", "package assumption", "size_repeat_copies",
    "size", "REF/ALT size difference (bp)", "alignment", "numeric", "c(3, Inf)", "bp",
    ">= 3 bp", "package assumption", "size_delta",
    "size", "fragment size (bp)", "alignment", "numeric", NA, "bp",
    "scheme window", "package assumption", "size_fragment",
    # §W: the scheme's design function refuses a multi-base target (see the
    # `cli_inform` in each `design_*_marker()`), so the scheme's inability is
    # recorded here as a criterion instead of only as console text.  A marker
    # that reaches grading with a multi-base target is `unknown` (the design was
    # never evaluated for it), never a silent pass.  One shared row: the families
    # that declare it in MARKER_FAMILY_CHECKS do the scoping.
    "single_base", "target is a single-base substitution", "site", "boolean", "= 1", NA,
    "REF/ALT 1 nt each", "package assumption", "single_base"
  )
}

# A family check's threshold as a bound object: the numeric `c()` expression its
# text spells out, or the literal level (a `=` predicate for a boolean-style one).
criterion_parse_bound <- function(text, kind) {
  if (length(text) != 1L || is.na(text)) {
    return(NULL)
  }
  if (kind == "numeric" && grepl("^c\\(", text)) {
    out <- tryCatch(eval(parse(text = text)), error = function(e) NULL)
    if (is.numeric(out)) {
      return(out)
    }
  }
  text
}

# The family checks as registry rows, one per (family, check), scoped to the
# family that declares it (`scope = "family"`, the `family` column naming it).
# The fragment-size window is the scheme's own product window, so it is filled
# lazily (`NA` bounds, prose only) exactly as it is published.
marker_family_check_spec <- function() {
  ways <- marker_family_check_ways()
  tibble::tibble(
    family = ways$family,
    id = paste0("family_", ways$evaluator),
    scope = "family",
    stage = "assay",
    role = "check",
    label = ways$label,
    level = ways$level,
    kind = ways$kind,
    bounds = unname(mapply(
      criterion_parse_bound, ways$bounds, ways$kind,
      SIMPLIFY = FALSE
    )),
    units = ways$units,
    source = ways$source,
    calibrated = FALSE,
    evidence = "rule_of_thumb",
    evaluator = ways$evaluator,
    hard = FALSE,
    limit = ways$limit
  )
}

# `marker_family_check_rows()`' legacy shape (stage/role/label/level/hard/limit),
# built from the registry so the family vocabulary cannot fork.
marker_family_check_rows <- function() {
  fam <- marker_family_check_ways()
  tibble::tibble(
    stage = "assay",
    role = "check",
    label = fam$label,
    level = fam$level,
    hard = FALSE,
    limit = fam$limit
  )
}

# Every rule in the SNaPshot provenance table, as a registry row.  One row per
# CSV rule id, `moderate` / `weak` / `none` mapped DOWN the ladder by
# CRITERION_TIER_EVIDENCE (never upgraded).  `source` names the registry key of
# the paper in inst/extdata/scientific_sources.csv; the DOI / PMCID live on that
# row, not in the key.
#
# `label` is the rule id: a paper's prose is a sentence, not a criterion name, and
# the CSV holds no short label per rule.  Role `rule`, not `check`: no designer
# evaluates these, so they must not enter the decision graph (which draws every
# `stage == "assay"` row) nor a verdict (which reads every `role == "check"` row).
# `stage = "snapshot"` is the scheme's own stage id, the same shape a discovered
# scheme's stage rows take.
# shortcut: the table is read on every `criterion_spec()` call, install it as
# package data if that call ever shows up in a profile.
snapshot_rule_spec <- function() {
  # full provenance strings are held once, so each rule row stays readable and
  # under the line-length limit
  src_atanasovska <- "atanasovska2012"
  src_namdar <- "namdar2023"
  src_coutinho <- "coutinho2014"
  src_zhao <- "zhao2026_snapshot"
  src_dicristofaro <- "dicristofaro2010"
  src_fanis <- "fanis2014"
  src_mrazkova <- "mrazkova2021"
  rule <- function(id, level, tier, source) {
    criterion_tbl(
      id = paste0("family_", id), scope = "family", stage = "snapshot",
      role = "rule", label = id, level = level, kind = "categorical",
      bounds = lapply(id, function(.) "stated"),
      units = NA_character_, source = source, calibrated = FALSE,
      evidence = unname(CRITERION_TIER_EVIDENCE[[tier]]), evaluator = id
    ) |>
      dplyr::mutate(family = "snapshot", .after = "scope")
  }
  rows <- dplyr::bind_rows(
    rule("sbe_3prime_adjacent", "feature", "moderate", src_atanasovska),
    rule("sbe_3prime_no_flexibility", "feature", "moderate", src_atanasovska),
    rule("sbe_core_length", "feature", "moderate", src_namdar),
    rule("sbe_spacing_min", "alignment", "moderate", src_atanasovska),
    rule("sbe_spacing_used_4bp", "alignment", "weak", src_coutinho),
    rule("sbe_spacing_used_ct5bp", "alignment", "weak", src_zhao),
    rule("sbe_tail_5prime", "feature", "moderate", src_namdar),
    rule("sbe_tail_keeps_tm", "feature", "moderate", src_atanasovska),
    rule("sbe_amp_primer_tm", "feature", "weak", src_namdar),
    rule("sbe_sbe_tm_gc_window", "feature", "weak", src_zhao),
    rule("sbe_inter_primer_tm_delta", "feature", "weak", src_zhao),
    rule("sbe_primer_ratio_reported", "feature", "none", src_dicristofaro),
    rule("sbe_peak_height_floor", "alignment", "weak", src_dicristofaro),
    rule("sbe_het_balance_floor", "alignment", "weak", src_zhao),
    rule("sbe_dye_yield_asymmetry", "alignment", "weak", src_dicristofaro),
    rule("sbe_structure_check", "feature", "moderate", src_atanasovska),
    rule("sbe_nonhomologous_window", "alignment", "weak", src_fanis),
    rule("sbe_no_template_control", "alignment", "moderate", src_mrazkova),
    rule("sbe_apply_format", "feature", "weak", src_dicristofaro)
  )
  # the rung is read back from the CSV, so the code's copy of the tier can never
  # silently disagree with the recorded provenance
  csv <- utils::read.csv(
    system.file("extdata", "snapshot_multiplex_design_rules.csv", package = "tidymarkers"),
    stringsAsFactors = FALSE
  )
  tier <- csv$evidence_tier[match(sub("^family_", "", rows$id), csv$rule_id)]
  rows$evidence <- unname(CRITERION_TIER_EVIDENCE[tier])
  rows
}

# One acceptance check as a registry row: its placement (the scheme stage that
# selects it) comes from the class it reports under, its limit from the bound.
marker_check_spec <- function(ways, classes, limits) {
  stage <- ifelse(ways$level == "site", "site", "assay")
  # a literal bound where the check declares one, else NULL for the numeric bar
  # filled from the limit object in force
  b <- lapply(seq_len(nrow(ways)), function(i) {
    lit <- ways$bounds[[i]]
    if (is.null(lit) || (length(lit) == 1L && is.na(lit))) NULL else lit
  })
  num <- marker_check_bounds(limits)
  fill <- ways$id %in% names(num)
  b[fill] <- unname(num[ways$id[fill]])
  scope <- unname(c(primer = "family", site = "family", as = "kasp")[
    classes$class[match(ways$id, classes$id)]
  ])
  criterion_tbl(
    id = ways$id, scope = scope, stage = stage, role = "check",
    label = ways$label, level = ways$level, kind = ways$kind, bounds = b,
    units = ways$units, source = ways$source,
    calibrated = FALSE,
    evaluator = ways$evaluator
  )
}

#' Construct a criterion table
#'
#' `criterion_spec()` returns a tibble subclass, not a bare tibble, so the
#' columns every consumer relies on are checked once at construction instead of
#' failing as `NULL` at a distant call site.  The class carries no behaviour of
#' its own: it is the typed table ZA6 step 1 asks for, and `vctrs` keeps it
#' through subsetting and row-binding.
#' @param x a list or data frame of criterion rows.
#' @return a `marker_criteria` tibble.
#' @keywords internal
#' @noRd
new_marker_criteria <- function(x) {
  # `family` must be here: the family check rows are the only ones that carry it,
  # and dropping it in the constructor would hide the column the spec documents
  # and every `scope == "family"` reader expects (the bind to the stage rows
  # leaves it `NA` there, never absent).
  fields <- c(
    "id", "scope", "family", "stage", "role", "label", "level", "kind",
    "bounds", "units", "source", "calibrated", "evidence", "evaluator", "hard",
    "limit", "order", "schemes", "asymmetry_reason", "target"
  )
  need_cols(x, fields, arg = "criterion table")
  vctrs::new_data_frame(
    x[, fields, drop = FALSE],
    class = c("marker_criteria", "tbl_df", "tbl", "data.frame")
  )
}

#' @exportS3Method vctrs::vec_restore
vec_restore.marker_criteria <- function(x, to, ...) new_marker_criteria(x)

# The double-dispatch method names below are fixed by vctrs (`generic.class.class`)
# and cannot be shortened without breaking dispatch, so they are exempt from the
# 30-character name linter.
# nolint start: object_length_linter.
#' @exportS3Method vctrs::vec_ptype_abbr
vec_ptype_abbr.marker_criteria <- function(x, ...) "crtr"

#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.marker_criteria.marker_criteria <- function(x, y, ...) {
  new_marker_criteria(tibble::as_tibble(x)[0, , drop = FALSE])
}

#' @exportS3Method vctrs::vec_cast
vec_cast.marker_criteria.marker_criteria <- function(x, to, ...) new_marker_criteria(x)

# A criterion table stays a criterion table when row-bound through vctrs to a
# plain tibble, and a plain tibble is downgraded rather than silently upgraded,
# so the common type of the two is a `tbl_df`, never a half-typed result.
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.marker_criteria.tbl_df <- function(x, y, ...) tibble::as_tibble(x)[0, , drop = FALSE]

#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.tbl_df.marker_criteria <- function(x, y, ...) tibble::as_tibble(y)[0, , drop = FALSE]

#' @exportS3Method vctrs::vec_cast
vec_cast.tbl_df.marker_criteria <- function(x, to, ...) tibble::as_tibble(x)

#' @exportS3Method vctrs::vec_cast
vec_cast.marker_criteria.tbl_df <- function(x, to, ...) new_marker_criteria(x)
# nolint end

#' The criterion registry: one data object per threshold
#'
#' Every gate, ranking tier and acceptance check the pipeline applies, as one
#' row each: its stable `id`, the `scope` it applies to (a scheme id, `"kasp"`
#' for the shared KASP-family stages, or `"family"` for the per-family
#' acceptance checks), `stage` and `role` (`gate`/`tier`/`check`), the `label`
#' used in the decision tree and the results, the `level` it reports at, its
#' `kind` (`numeric`, `categorical` or `boolean`), the `bounds` it is judged
#' against (a list-column: `c(min, max)` with `NA` for an open side,
#' `c(min, "exclusive", x)` for a strict `<`, the allowed levels, or the literal
#' predicate), the `units`, where the threshold comes from (`source`;
#' `"package assumption"` where no literature value was claimed), the
#' `evidence` rung behind it (an ordered factor, see
#' `criterion_evidence()`; every row on this base is `"rule_of_thumb"` because no
#' source was ever claimed for these numbers), the `evaluator` that produces
#' the observed value, and the `asymmetry_reason` that names the structural
#' fact breaking a metamorphic symmetry for the rows where one deliberately
#' does not hold (`NA` everywhere else, see `criterion_asymmetry()`; ZW0.1).
#'
#' This is internal because `bounds` is only meaningful together with the
#' `limits` in force; [marker_criteria()] is the public view, and its `limit`
#' column is this object rendered.  The selector, the acceptance checks and the
#' decision graph all read it, so the same criterion cannot carry two different
#' bars (Z2).
#' @param limits thresholds (see [marker_design_checks()]).
#' @return tibble, one row per registered criterion (a `family` column names the
#'   family a `scope = "family"` row belongs to).  The `calibrated` column is
#'   **deprecated** (ZW4): it is the boolean view of `evidence`, still emitted so
#'   existing readers keep working; read `evidence` instead.
#' @keywords internal
#' @noRd
criterion_spec <- function(limits = KASP_LIMITS) {
  # The registry is a pure function of `limits` *and* the discovered schemes
  # (their gate/tier rows are bound in), and it is rebuilt on every
  # `marker_family_check_labels()` call - which `marker_generic_checks()` used to
  # call once per check row, so a design run rebuilt it thousands of times.
  # Measured on test-falsifiable.R: 68.5 s of a 95 s file was under this one
  # call.  Memoised under the same key shape `marker_criteria()` uses.
  session_memo(
    paste0(
      "criterion_spec|", scheme_registry_token(), "|",
      paste(deparse(limits), collapse = "")
    ),
    criterion_spec_build(limits)
  )
}

#' The criterion registry, built (the uncached body of `criterion_spec()`)
#' @param limits thresholds.
#' @keywords internal
#' @noRd
criterion_spec_build <- function(limits = KASP_LIMITS) {
  # a gate's scope is the scheme that declares its constant; the `as`/`pair`
  # stages are shared by the KASP family, the rest belong to their own scheme.
  # One spec for all of them: marker_builtin_stages() (R/marker_stages.R).
  builtins <- marker_builtin_stages()
  gate_labels <- lapply(builtins, marker_stage_constants, what = "gates")
  relax_labels <- lapply(builtins, marker_stage_constants, what = "relax")
  tier_labels <- lapply(builtins, marker_stage_constants, what = "tiers")
  gate_scope <- vapply(builtins, `[[`, character(1), "scope")
  level_of <- vapply(builtins, `[[`, character(1), "level")
  # A gate's bar is the scheme's published limit text, and its relaxation margin
  # (the floor) is a ladder property, not a threshold on the observed value, so a
  # gate carries `bounds = NULL` and its prose as `values`.  The rows whose
  # `bounds` *is* the threshold are the acceptance checks below.
  # A gate vector is unnamed: its *values* are the gate ids (the tier vectors are
  # the named ones, where the name is the tier key).
  gate_rows <- function(stage, g) {
    out <- criterion_tbl(
      id = paste0(stage, ":gate:", g), scope = gate_scope[[stage]],
      stage = stage, role = "gate", label = scheme_label("short", g),
      level = level_of[[stage]] %||% "alignment",
      kind = rep("boolean", length(g)),
      bounds = vector("list", length(g)),
      units = rep(NA_character_, length(g)),
      source = "package assumption", calibrated = FALSE,
      evidence = "rule_of_thumb",
      evaluator = rep("selector", length(g)),
      values = scheme_label("limit", g)
    )
    out$hard <- !(g %in% relax_labels[[stage]])
    out$order <- seq_along(g)
    out
  }
  tier_rows <- function(stage, t) {
    out <- criterion_tbl(
      id = paste0(stage, ":tier:", names(t)), scope = gate_scope[[stage]],
      stage = stage, role = "tier", label = unname(t),
      level = level_of[[stage]] %||% "alignment",
      kind = rep("categorical", length(t)),
      bounds = vector("list", length(t)),
      units = rep(NA_character_, length(t)),
      source = "package assumption", calibrated = FALSE,
      evidence = "rule_of_thumb",
      evaluator = rep("tier", length(t)),
      values = rep(CRITERION_TIER_LIMIT, length(t))
    )
    out$order <- seq_along(t)
    out
  }
  # A stage is run by a *family* (see marker_builtin_stages()), and its gate/tier
  # rows are tagged with every scheme of that family - which is what lets a
  # composed technology inherit its base's stage instead of losing it.
  stage_scheme <- lapply(builtins, function(st) stage_schemes(st$schemes))
  stages <- lapply(names(builtins), function(nm) {
    g <- gate_labels[[nm]]
    t <- tier_labels[[nm]]
    d <- dplyr::bind_rows(
      if (length(g)) gate_rows(nm, g) else NULL,
      if (length(t)) tier_rows(nm, t) else NULL
    )
    if (!nrow(d)) {
      return(NULL)
    }
    member <- stage_scheme[[nm]]
    d$schemes <- rep(list(member), nrow(d))
    d
  })
  stages <- purrr::compact(stages)

  ways <- marker_check_ways()
  checks <- marker_check_spec(ways, marker_check_classes(), limits)
  # SNaPshot literature rows (ROADMAP §U) are sourced claims, not evaluations:
  # they keep `role = "rule"` for their whole life, so no reader that selects
  # `role == "check"` (or an `assay`-stage row) can show one as an evaluation
  checks <- dplyr::bind_rows(checks, marker_family_check_spec())
  # `family` must exist on every row before the subset: tibble subsetting treats
  # a missing first column as all-`NA`, which `[` then drops
  checks$family <- checks$family %||% rep(NA_character_, nrow(checks))
  checks <- checks[checks$role %in% c("gate", "check"), , drop = FALSE]
  checks$order <- seq_len(nrow(checks))
  # `schemes` is decided by the check's *class* (marker-agnostic / site / allele
  # specific / family), which the class rules recognise by their registered
  # labels; family membership itself is by id (see `check_schemes()`).
  checks$schemes <- check_schemes(checks$label)

  # the SNaPshot provenance rows (inst/extdata/snapshot_multiplex_design_rules.csv)
  # keep the rung the CSV records; they are bound last so the per-scheme defaults
  # of `extra` above cannot overwrite them
  extra <- extra_stage_rows()
  if (!is.null(extra)) {
    n <- nrow(extra)
    extra$kind <- "boolean"
    extra$bounds <- vector("list", n)
    extra$units <- NA_character_
    extra$source <- "package assumption"
    extra$calibrated <- FALSE
    extra$evidence <- "rule_of_thumb"
    extra$evaluator <- extra$role
  }
  rules <- snapshot_rule_spec()
  # set here, after `check_schemes()` ran on the checks: the rules all serve one
  # scheme and are not on any class the membership map knows
  rules$schemes <- rep(list("snapshot"), nrow(rules))
  rules$order <- nrow(checks) + seq_len(nrow(rules))
  out <- dplyr::bind_rows(stages, checks, extra, rules)
  out <- out[out$role %in% c("gate", "tier", "check", "rule"), , drop = FALSE]
  # ZW4.6: overlay the rows whose bar a published source states; every other
  # row keeps the constructor's "package assumption".
  src <- criterion_sources()
  hit <- match(out$id, names(src))
  out$source[!is.na(hit)] <- src[hit[!is.na(hit)]]
  # the rung behind a raised row: `evidence` may be an ordered factor at this
  # point (snapshot rules), so coerce before replacing and re-impose below
  out$evidence <- as.character(out$evidence)
  ev <- criterion_evidence_sources()
  hit <- match(out$id, names(ev))
  out$evidence[!is.na(hit)] <- unname(ev[hit[!is.na(hit)]])
  # ZW4.6: the reference point a tier ranks toward (the citable half of the row;
  # the aggregation policy stays `CRITERION_TIER_POLICY_EVIDENCE` by construction)
  out$target <- rep(NA_character_, nrow(out))
  tg <- criterion_tier_targets()
  hit <- match(out$id, names(tg))
  out$target[!is.na(hit)] <- unname(tg[hit[!is.na(hit)]])
  # ZW0.1: declare where a metamorphic invariant deliberately does not hold;
  # every other row is symmetric and stays `NA`.
  asym <- criterion_asymmetry()
  hit <- match(out$id, names(asym))
  out$asymmetry_reason <- rep(NA_character_, nrow(out))
  out$asymmetry_reason[!is.na(hit)] <- unname(asym[hit[!is.na(hit)]])
  # `evidence` is appended right after the deprecated `calibrated` so that no
  # existing column moves (downstream code indexes the spec positionally).
  # `asymmetry_reason` is appended last for the same reason (ZW0.1).
  cols <- c(
    "id", "scope", "family", "stage", "role", "label", "level", "kind",
    "bounds", "units", "source", "calibrated", "evidence", "evaluator", "hard",
    "limit", "order", "schemes", "asymmetry_reason", "target"
  )
  # `bind_rows()` drops an ordered factor's class when the other side of the bind
  # is a plain character column, so the ladder is re-imposed once, on the whole
  # spec: every constructor above names a rung, this makes it the ordered factor.
  out$evidence <- criterion_evidence(out$evidence)
  new_marker_criteria(out[, intersect(cols, names(out)), drop = FALSE])
}

#' The registry of every candidate evaluation
#'
#' One row per criterion the pipeline applies: the selection gates and ranking
#' tiers of each stage, and the acceptance checks (composition, structure,
#' heteroduplexes, repeats, off-target specificity, in-silico PCR, site facts).
#' This is the single vocabulary the selector, [marker_design_checks()], the
#' decision graph and [marker_decision_text()] all read.  Each row carries the
#' marker `schemes` it serves ([marker_schemes()]), so the same registry yields one
#' decision schema per chemistry.
#'
#' The thresholds themselves live in the `bounds` column of `criterion_spec()`;
#' this view renders them into `limit`, so the same criterion cannot carry two
#' different bars.
#' @param limits thresholds (see [marker_design_checks()]).
#' @param scheme optional scheme id (see [marker_schemes()]); when given, only the
#'   criteria that scheme serves are returned.
#' @return tibble with `stage`, `role` (`gate`/`tier`/`check`), `id`, `label`,
#'   `level` (`site`/`feature`/`alignment`/`specificity`), `hard`, `limit`,
#'   `target` (a tier's printed reference point, `NA` elsewhere), `order`,
#'   `schemes` (list-column of scheme ids).
#' @examples
#' marker_criteria()
#' marker_criteria(scheme = "duplex")
#' @family checks
#' @export
marker_criteria <- function(limits = KASP_LIMITS, scheme = NULL) {
  # `asq` is graded against its own windows (100-200 bp, 62-64 C); a bare
  # `marker_criteria(scheme = "asq")` shows them too.  An assay's own `$limits`
  # (passed by `marker_design_checks()`) always wins over this default.
  if (!is.null(scheme) && identical(scheme, "asq") && identical(limits, KASP_LIMITS)) {
    limits <- asq_limits()
  }
  # a product read by size is graded against the scheme's own window, never the
  # KASP 50-100 bp optimum (C2/C5).  Overriding `amplicon` here keeps the gate,
  # the check, its `bounds` and the rendered `limit` in step.
  if (!is.null(scheme)) {
    limits$amplicon <- marker_amplicon_window(marker_scheme(scheme)$family, limits)
  }
  # ~0.2 s per call, and the decision schema asks for it once per marker; the
  # answer is a pure function of `limits` (a flat numeric list) and of the
  # discovered schemes (`extra_stage_rows()`), so both are in the key.
  key <- paste0("criteria|", scheme_registry_token(), "|", paste(deparse(limits), collapse = ""))
  out <- session_memo(key, {
    spec <- criterion_spec(limits)
    tibble::tibble(
      stage = spec$stage,
      role = spec$role,
      id = spec$id,
      label = spec$label,
      level = spec$level,
      hard = spec$hard,
      limit = spec$limit,
      target = spec$target,
      order = spec$order,
      schemes = spec$schemes
    )
  })
  if (is.null(scheme)) {
    return(out)
  }
  out[vapply(out$schemes, function(s) scheme %in% s, logical(1)), , drop = FALSE]
}
