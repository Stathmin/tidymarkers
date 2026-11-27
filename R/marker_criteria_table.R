# =============================================================================
#  marker_criteria_table.R  --  bounds, family rows, vctrs class, spec build
# =============================================================================
#  Split from marker_criteria.R (roadmap ZT T7).
# =============================================================================
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
