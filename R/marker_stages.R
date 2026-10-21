# =============================================================================
#  marker_stages.R  --  the built-in selection stages, declared once
# =============================================================================
#  A selection stage is one row here: its gates (in the order the method
#  applies them), its relaxation ladder, its tiers, the level its criteria
#  report at, the scope its rows take, and the schemes that own it.  Before
#  this list existed the same seven stages were spelled four times - three
#  named constant maps plus a `level_of` map in marker_criteria.R, and six
#  near-identical `if (family == ...)` blocks in marker_decision.R - so a new
#  built-in stage meant four edits that could silently disagree.
#
#  A function, not a top-level constant: marker_criteria.R is sourced before the
#  scheme files, so reading `KASP_GATES_*` at load time fails (the same reason
#  marker_amplicon_window() is a function).  The constants are the single source
#  of the values; this list only says which stage owns which.
# =============================================================================

#' The built-in selection stages, in pipeline order
#'
#' One entry per stage, each a list of
#' - `label` - the stage's display name;
#' - `gates`, `relax`, `tiers` - the per-scheme constants (see `KASP_GATES_*`),
#'   named as strings so the lookup stays lazy;
#' - `level` - the level the stage's criteria rows report at (`"feature"` for the
#'   AS stage, `"alignment"` otherwise);
#' - `schemes` - the scheme ids that run the stage.  The `as`/`pair` stages are
#'   named `"kasp"` in `scope` (their criterion rows belong to the KASP family)
#'   but `duplex` and `asq` run them too, so membership is listed, never
#'   inferred from the family.
#'
#' Only built-in stages are listed.  A discovered scheme's stage
#' (lax, ssr, sanger, ...) is appended by `extra_stage_rows()` /
#' `extra_decision_stages()`, so listing it here would emit it twice.
#' @return a named list of stage specs.
#' @keywords internal
#' @noRd
marker_builtin_stages <- function() {
  stage <- function(label, gates, relax, tiers, level, scope, schemes) {
    list(
      label = label, gates = gates, relax = relax, tiers = tiers,
      level = level, scope = scope, schemes = schemes
    )
  }
  list(
    as = stage(
      "AS primer", "KASP_GATES_AS", "KASP_RELAX_AS", "KASP_TIERS_AS",
      "feature", "kasp", c("kasp", "multiplex", "asq")
    ),
    pair = stage(
      "AS x COM pair", "KASP_GATES_PAIR", "KASP_RELAX_PAIR",
      "KASP_TIERS_PAIR", "alignment", "kasp", c("kasp", "multiplex", "asq")
    ),
    duplex = stage(
      "duplex compatibility", "KASP_GATES_DUPLEX", "KASP_RELAX_DUPLEX",
      "KASP_TIERS_DUPLEX", "alignment", "duplex", "multiplex"
    ),
    asq = stage(
      "ASQ primer geometry", "KASP_GATES_ASQ", "KASP_RELAX_ASQ",
      "KASP_TIERS_ASQ", "alignment", "asq", "asq"
    ),
    caps = stage(
      "CAPS digest", "KASP_GATES_CAPS", "KASP_RELAX_CAPS",
      "KASP_TIERS_CAPS", "alignment", "caps", "caps"
    ),
    arms = stage(
      "ARMS bands", "KASP_GATES_ARMS", "KASP_RELAX_ARMS",
      "KASP_TIERS_ARMS", "alignment", "arms", "arms"
    ),
    hrm = stage(
      "HRM amplicon", "KASP_GATES_HRM", "KASP_RELAX_HRM",
      "KASP_TIERS_HRM", "alignment", "hrm", "hrm"
    )
  )
}

#' Resolve a stage spec's named constants (lazily, by name)
#' @keywords internal
#' @noRd
marker_stage_constants <- function(stage, what = c("gates", "relax", "tiers")) {
  what <- match.arg(what)
  nm <- stage[[what]]
  if (is.null(nm)) {
    return(NULL)
  }
  get(nm, envir = asNamespace("tidymarkers"))
}
