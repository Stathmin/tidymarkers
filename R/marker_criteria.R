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
# are in, where the number comes from (`source`) and whether it was `calibrated`
# against held-out data (almost none are), and the `evaluator` that produces the
# observed value.
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

# Criterion kinds: a numeric bar, a set of allowed levels, or a tested predicate.
CRITERION_KINDS <- c("numeric", "categorical", "boolean")



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
  if (is.na(lo)) {
    return(glue::glue("<= {criterion_num(hi)}"))
  }
  if (is.na(hi)) {
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
criterion_tbl <- function(id, scope, stage, role, label, level, kind, bounds,
                          units, source, calibrated, evaluator, values = NULL) {
  n <- length(id)
  stopifnot(all(kind %in% CRITERION_KINDS))
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
    source = rep(source, length.out = n),
    calibrated = rep(calibrated, length.out = n),
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
# This is the one place the class vectors and the registry's ids are tied.
marker_check_classes <- function() {
  tibble::tribble(
    ~id, ~class,
    "site_class", "site",
    "site_alt_alleles", "site",
    "site_linked_skipped", "site",
    "site_alt_covered", "site",
    "site_linked_dependent", "site",
    "amplicon", "primer",
    "as_tm", "as",
    "com_tm", "as",
    "com_tm_ge_as_tm", "as",
    "as_gc", "as",
    "com_gc", "as",
    "max_run", "primer",
    "ssr", "primer",
    "inverted_repeat", "primer",
    "hairpin_dg", "primer",
    "cross_dimer_dg", "primer",
    "homodimer_dg", "primer",
    "as_n_3p", "as",
    "as_n_var", "as",
    "com_n_var", "as",
    "as_dropout", "as",
    "as_foot_dropout", "as",
    "com_dropout", "as",
    "discrimination", "as",
    "in_silico_pcr", "primer",
    "offtarget_hits", "primer"
  )
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
    "in_silico_pcr", "in-silico PCR (REF 1 / ALT 0)", "alignment", "boolean",
    "REF=1, ALT=0", NA, "package assumption", "in_silico_pcr", NA,
    "offtarget_hits", "off-target exact primer hits (n)", "specificity", "numeric", NULL,
    "count", "package assumption", "ideogram_hits", NA
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
    offtarget_hits = c(NA, 1)
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
    "scheme window", "package assumption", "size_fragment"
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
    units = ways$units, source = ways$source, calibrated = FALSE,
    evaluator = ways$evaluator
  )
}

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
#' `"package assumption"` where no literature value was claimed), whether it was
#' `calibrated` on held-out data, and the `evaluator` that produces the observed
#' value.
#'
#' This is internal because `bounds` is only meaningful together with the
#' `limits` in force; [marker_criteria()] is the public view, and its `limit`
#' column is this object rendered.  The selector, the acceptance checks and the
#' decision graph all read it, so the same criterion cannot carry two different
#' bars (Z2).
#' @param limits thresholds (see [marker_design_checks()]).
#' @return tibble, one row per registered criterion (a `family` column names the
#'   family a `scope = "family"` row belongs to).
#' @keywords internal
#' @noRd
criterion_spec <- function(limits = KASP_LIMITS) {
  # a gate's scope is the scheme that declares its constant; the `as`/`pair`
  # stages are shared by the KASP family, the rest belong to their own scheme.
  # Resolved lazily (`mget`) because marker_criteria.R is sourced before them.
  gate_scope <- c(
    as = "kasp", pair = "kasp", duplex = "duplex", asq = "asq",
    caps = "caps", arms = "arms", hrm = "hrm"
  )
  # Only the built-in stages are listed here.  A discovered scheme's stage
  # (lax, ssr, sanger, ...) is appended by extra_stage_rows(), so listing it here
  # too would emit its gates twice.
  stage_gates <- c(
    as = "KASP_GATES_AS", pair = "KASP_GATES_PAIR", duplex = "KASP_GATES_DUPLEX",
    asq = "KASP_GATES_ASQ", caps = "KASP_GATES_CAPS", arms = "KASP_GATES_ARMS",
    hrm = "KASP_GATES_HRM"
  )
  stage_relax <- c(
    as = "KASP_RELAX_AS", pair = "KASP_RELAX_PAIR", duplex = "KASP_RELAX_DUPLEX",
    asq = "KASP_RELAX_ASQ", caps = "KASP_RELAX_CAPS", arms = "KASP_RELAX_ARMS",
    hrm = "KASP_RELAX_HRM"
  )
  stage_tiers <- c(
    as = "KASP_TIERS_AS", pair = "KASP_TIERS_PAIR", duplex = "KASP_TIERS_DUPLEX",
    asq = "KASP_TIERS_ASQ", caps = "KASP_TIERS_CAPS", arms = "KASP_TIERS_ARMS",
    hrm = "KASP_TIERS_HRM"
  )
  # `mget()` names its result after the constants, so re-key it by stage id
  gate_labels <- stats::setNames(
    mget(unname(stage_gates), envir = asNamespace("tidymarkers")),
    names(stage_gates)
  )
  relax_labels <- stats::setNames(
    mget(unname(stage_relax), envir = asNamespace("tidymarkers")),
    names(stage_relax)
  )
  tier_labels <- stats::setNames(
    mget(unname(stage_tiers), envir = asNamespace("tidymarkers")),
    names(stage_tiers)
  )
  # the AS stage's criteria report at feature level; every other stage's at
  # alignment level (as the old stage_rows() calls declared one by one)
  level_of <- c(
    as = "feature", pair = "alignment", duplex = "alignment",
    asq = "alignment", caps = "alignment", arms = "alignment",
    hrm = "alignment"
  )
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
      evaluator = rep("tier", length(t)),
      values = rep(CRITERION_TIER_LIMIT, length(t))
    )
    out$order <- seq_along(t)
    out
  }
  # The `as` / `pair` stages are shared by the whole KASP family; every other
  # built-in stage belongs to one scheme.
  marker_core <- c(KASP_SCHEME_CORE, extra_family_schemes("kasp"))
  stage_scheme <- c(
    as = "kasp", pair = "kasp", duplex = "duplex", asq = "asq",
    caps = "caps", arms = "arms", hrm = "hrm"
  )
  stages <- purrr::map(names(stage_gates), function(stage) {
    g <- gate_labels[[stage]]
    t <- tier_labels[[stage]]
    d <- dplyr::bind_rows(
      if (length(g)) gate_rows(stage, g) else NULL,
      if (length(t)) tier_rows(stage, t) else NULL
    )
    if (!nrow(d)) {
      return(NULL)
    }
    member <- if (stage %in% c("as", "pair")) marker_core else stage_scheme[[stage]]
    d$schemes <- rep(list(member), nrow(d))
    d
  })
  stages <- purrr::compact(stages)

  ways <- marker_check_ways()
  checks <- marker_check_spec(ways, marker_check_classes(), limits)
  checks <- dplyr::bind_rows(checks, marker_family_check_spec())
  checks$order <- seq_len(nrow(checks))
  # `schemes` is decided by the check's *class* (marker-agnostic / site / allele
  # specific / family), which the class rules recognise by their registered
  # labels; family membership itself is by id (see `check_schemes()`).
  checks$schemes <- check_schemes(checks$label)

  extra <- extra_stage_rows()
  if (!is.null(extra)) {
    n <- nrow(extra)
    extra$kind <- "boolean"
    extra$bounds <- vector("list", n)
    extra$units <- NA_character_
    extra$source <- "package assumption"
    extra$calibrated <- FALSE
    extra$evaluator <- extra$role
  }

  out <- dplyr::bind_rows(stages, checks, extra)
  cols <- c(
    "id", "scope", "family", "stage", "role", "label", "level", "kind",
    "bounds", "units", "source", "calibrated", "evaluator", "hard", "limit",
    "order", "schemes"
  )
  out[, intersect(cols, names(out)), drop = FALSE]
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
#'   `order`, `schemes` (list-column of scheme ids).
#' @examples
#' marker_criteria()
#' marker_criteria(scheme = "duplex")
#' @family KASP checks
#' @export
marker_criteria <- function(limits = KASP_LIMITS, scheme = NULL) {
  # a product read by size is graded against the scheme's own window, never the
  # KASP 50-100 bp optimum (C2/C5).  Overriding `amplicon` here keeps the gate,
  # the check, its `bounds` and the rendered `limit` in step.
  if (!is.null(scheme)) {
    limits$amplicon <- marker_amplicon_window(marker_scheme(scheme)$family, limits)
  }
  spec <- criterion_spec(limits)
  out <- tibble::tibble(
    stage = spec$stage,
    role = spec$role,
    id = spec$id,
    label = spec$label,
    level = spec$level,
    hard = spec$hard,
    limit = spec$limit,
    order = spec$order,
    schemes = spec$schemes
  )
  if (is.null(scheme)) {
    return(out)
  }
  out[vapply(out$schemes, function(s) scheme %in% s, logical(1)), , drop = FALSE]
}
