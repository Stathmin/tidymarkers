# =============================================================================
#  marker_plan.R  --  the design plan: a layered grammar of marker design
# =============================================================================
#  Designing markers is a *plan* built up in layers, then solved: start from the
#  genome with marker_plan(), add the target loci (marker_targets()), the variant /
#  haplotype context (marker_context()) and the available chemistries with their
#  costs (marker_tech()), then run marker_solve(plan).
#  `marker_plan()` is the base (like `ggplot()`), `marker_targets()`/`marker_context()`/
#  `marker_tech()` are layers added with `+` (like `geom_*()`), and `marker_solve()` is
#  the compute verb (like `ggplot_build()`).  The plan holds only *what* to do, so
#  it is cheap to build and inspect; the cost of designing lives in the solve.
# =============================================================================

# Relative cost of each chemistry (bigger = harder / more expensive to run).  The
# solver walks the allowed schemes cheapest-first and stops at the first that
# designs, so KASP is tried long before fragment analysis, and fragment analysis
# is never run when a KASP already answers the target.
KASP_TECH_COSTS <- c(
  snp = 1, lax = 1.5, multi = 1.5, duplex = 2, hrm = 3, arms = 4, caps = 5,
  castpcr = 6, taqman = 6, indel_size = 10, ssr = 12, snapshot = 15, dpcr = 20,
  lamp = 25, rpa = 25, crispr = 30, sanger = 40
)

# Bundles of chemistries a workflow typically has on hand.  `full` is not a
# preset id; leaving `schemes` unset already means "every designable scheme".
KASP_TECH_PRESETS <- list(
  kasp = c("snp", "lax", "multi", "duplex"),
  gel = c("snp", "lax", "multi", "caps", "arms", "indel_size", "ssr"),
  lab = c("snp", "lax", "multi", "duplex", "hrm", "caps", "arms", "taqman", "castpcr", "dpcr"),
  field = c("snp", "lax", "multi", "lamp", "rpa", "crispr")
)

# schemes whose designer understands a haplotype table (the KASP assembly family)
haplotype_schemes <- function() c("snp", "indel", "multi", "lax")

# The columns every target table carries (an `alt2` keeps a multi-allelic site's
# second ALT, `type` its class; `label` falls back to the panel's `id`)
target_cols <- function() c("label", "chrom", "pos", "ref", "alt", "alt2", "type")

#' An empty target table
#' @keywords internal
#' @noRd
empty_targets <- function() {
  tibble::tibble(
    label = character(), chrom = character(), pos = integer(),
    ref = character(), alt = character(), alt2 = character(), type = character()
  )
}

#' Normalise a target table (derive a `label` if absent)
#' @keywords internal
#' @noRd
normalise_targets <- function(x) {
  x <- tibble::as_tibble(x)
  need <- c("chrom", "pos", "ref", "alt")
  missing <- setdiff(need, names(x))
  if (length(missing)) {
    cli::cli_abort(glue::glue(
      gettext("`targets` is missing column(s): {paste(missing, collapse = ', ')}.")
    ))
  }
  x$chrom <- as.character(x$chrom)
  x$pos <- as.integer(x$pos)
  x$ref <- toupper(as.character(x$ref))
  x$alt <- toupper(as.character(x$alt))
  if (!"label" %in% names(x)) {
    x$label <- if ("id" %in% names(x)) {
      as.character(x$id)
    } else {
      glue::glue("{x$chrom}:{x$pos} {x$ref}>{x$alt}")
    }
  }
  if (!"alt2" %in% names(x)) x$alt2 <- NA_character_
  x$alt2 <- as.character(x$alt2)
  if (!"type" %in% names(x)) x$type <- NA_character_
  x[, target_cols()]
}

#' One row (or a table) of target loci for [marker_targets()]
#' @keywords internal
#' @noRd
target_rows <- function(chrom, pos, ref, alt, data, type, label) {
  if (!is.null(data)) {
    return(normalise_targets(data))
  }
  if (is.null(chrom) || is.null(pos) || is.null(ref) || is.null(alt)) {
    cli::cli_abort(gettext("Give `chrom`, `pos`, `ref`, `alt`, or a `data` table."))
  }
  out <- tibble::tibble(
    chrom = as.character(chrom), pos = as.integer(pos),
    ref = toupper(as.character(ref)), alt = toupper(as.character(alt))
  )
  out$label <- if (is.null(label)) {
    glue::glue("{out$chrom}:{out$pos} {out$ref}>{out$alt}")
  } else {
    as.character(label)
  }
  out$alt2 <- NA_character_
  out$type <- if (is.null(type)) NA_character_ else as.character(type)
  out[, target_cols()]
}

#' A plan layer (internal; layers are added with `+`)
#' @keywords internal
#' @noRd
marker_layer <- function(apply, label) {
  structure(list(apply = apply, label = label), class = "marker_layer")
}

#' Start a marker-design plan
#'
#' The base of the grammar, like `ggplot()`.  Add layers with `+`:
#' [marker_targets()] (which loci), [marker_context()] (variants / haplotypes the
#' design should know about) and [marker_tech()] (which chemistries are available,
#' and how expensive each is).  Then [marker_solve()] runs the design.
#'
#' The plan is a plain list: it is cheap to build, print and inspect, and it
#' carries no results until you solve it.
#' @param genome a `marker_genome` or `DNAStringSet` (see [marker_genome()]).
#' @param targets optional target table (`chrom`, `pos`, `ref`, `alt`, and an
#'   optional `label`/`type`) to seed the plan; otherwise add them with
#'   [marker_targets()].
#' @param region optional single region (`GRanges`, or a data frame with
#'   `chrom`, `start`, `end`).  At solve time it expands to every variant the
#'   context places inside it, so the plan sweeps sites x chemistries without a
#'   target list (roadmap Y2.1); a variant outside the block is never a target.
#' @return a `marker_plan`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' plan <- marker_plan(genome) +
#'   marker_targets(marker_example("snp")) +
#'   marker_tech(preset = "lab")
#' plan
#' }
#' @family KASP plan
#' @export
marker_plan <- function(genome, targets = NULL, region = NULL) {
  if (is.null(genome)) {
    cli::cli_abort(gettext("`genome` is required (see marker_genome())."))
  }
  structure(
    list(
      genome = genome,
      targets = if (is.null(targets)) empty_targets() else normalise_targets(targets),
      region = if (is.null(region)) NULL else normalise_region(region),
      variants = NULL, haplotypes = NULL,
      schemes = NULL, costs = NULL, exclude = NULL,
      layers = character(0)
    ),
    class = "marker_plan"
  )
}

#' Add target loci to a plan
#'
#' A target layer: either a table (`data=`, e.g. [marker_example()] or your own),
#' or explicit columns.  Added twice, the rows accumulate (so several
#' [marker_targets()] layers compose).
#' @param data a target table (`chrom`, `pos`, `ref`, `alt`; optional `label`,
#'   `type`, `alt2`, `id`) used instead of the explicit columns.
#' @param chrom,pos,ref,alt one or more loci (vectors, recycled).
#' @param type optional free-text class per target (e.g. `"snp"`, `"indel"`).
#' @param label optional display label; derived from the coordinates otherwise.
#' @return a plan layer (add it with `+`).
#' @examples
#' \dontrun{
#' marker_plan(genome) + marker_targets(marker_example())
#' marker_plan(genome) + marker_targets(chrom = "Chr4", pos = 269116, ref = "G", alt = "A")
#' }
#' @family KASP plan
#' @export
marker_targets <- function(data = NULL, chrom = NULL, pos = NULL, ref = NULL,
                           alt = NULL, type = NULL, label = NULL) {
  marker_layer(
    function(plan) {
      plan$targets <- dplyr::bind_rows(plan$targets, target_rows(chrom, pos, ref, alt, data, type, label))
      plan
    },
    label = "targets"
  )
}

#' Add variant / haplotype context to a plan
#'
#' The design should know what else is polymorphic where it works: `variants`
#' lets it count and weight other variants under a primer, `haplotypes` (a tidy
#' haplotype table, see [assemble_kasp_assay()]) lets it bake co-segregating
#' sites into an ALT primer and mark a design as linkage-dependent.
#' @param variants a tidy variant table (`chrom`, `pos`, and either `alt` or
#'   `ref_genome`; optional `pop_freq`), e.g. from [ensembl_variants()].
#' @param haplotypes a tidy haplotype table (`haplotype`, `freq`, `chrom`,
#'   `pos`, `ref`, `alt`).
#' @return a plan layer (add it with `+`).
#' @examples
#' \dontrun{
#' marker_plan(genome) + marker_context(variants = ensembl_variants("Chr4", 268901, 271503))
#' }
#' @family KASP plan
#' @export
marker_context <- function(variants = NULL, haplotypes = NULL) {
  marker_layer(
    function(plan) {
      if (!is.null(variants)) plan$variants <- variants
      if (!is.null(haplotypes)) plan$haplotypes <- haplotypes
      plan
    },
    label = "context"
  )
}

#' Declare which chemistries a plan may use, and their costs
#'
#' The tech layer.  `schemes` is an allow-list (leaving it unset means every
#' designable scheme); `preset` names a bundle (`KASP_TECH_PRESETS`); `exclude`
#' drops schemes from the allowed set.  `costs` is a cost schema - a relative
#' effort per chemistry, bigger meaning harder - and the solver walks the allowed
#' schemes cheapest-first, stopping at the first that designs.
#'
#' A chemistry is "available" if its designer is installed (see [marker_schemes()]);
#' a custom chemistry you wrote as `design_<name>_marker()` + `scheme_def_<name>()`
#' is discovered automatically and can be named here like any other.
#' @param schemes scheme ids to allow (see [marker_schemes()]); `NULL` = all.
#' @param costs a named cost per scheme, e.g. [marker_costs()]; `NULL` = defaults.
#' @param preset a bundle name in `KASP_TECH_PRESETS` (`"kasp"`, `"gel"`,
#'   `"lab"`, `"field"`).
#' @param exclude scheme ids to drop from the allowed set.
#' @return a plan layer (add it with `+`).
#' @examples
#' \dontrun{
#' marker_plan(genome) + marker_tech(preset = "lab", costs = marker_costs(hrm = 1.2))
#' }
#' @family KASP plan
#' @export
marker_tech <- function(schemes = NULL, costs = NULL, preset = NULL, exclude = NULL) {
  if (!is.null(preset)) {
    preset <- match.arg(preset, names(KASP_TECH_PRESETS))
    schemes <- KASP_TECH_PRESETS[[preset]]
  }
  marker_layer(
    function(plan) {
      if (!is.null(schemes)) plan$schemes <- schemes
      if (!is.null(costs)) plan$costs <- costs
      if (!is.null(exclude)) plan$exclude <- exclude
      plan
    },
    label = "tech"
  )
}

#' A cost schema for the chemistries
#'
#' Returns the default relative-cost vector, optionally overridden.  Bigger means
#' harder: the solver will not run a more expensive chemistry when a cheaper one
#' already designs the target.  Overriding `snp = 1, caps = 10` says "KASP is ten
#' times easier than a CAPS assay" - so a KASP win pre-empts the CAPS search.
#' @param ... named cost overrides, e.g. `marker_costs(hrm = 1.2, caps = 8)`.
#' @return a named numeric vector (one cost per scheme).
#' @examples
#' marker_costs()
#' marker_costs(hrm = 1.2)
#' @family KASP plan
#' @export
marker_costs <- function(...) {
  ov <- c(...)
  if (length(ov) && (is.null(names(ov)) || any(!nzchar(names(ov))))) {
    cli::cli_abort(gettext("`...` must be named costs, e.g. marker_costs(hrm = 1.2)."))
  }
  base <- KASP_TECH_COSTS
  base[names(ov)] <- as.numeric(ov)
  base
}

#' @export
"+.marker_plan" <- function(e1, e2) {
  if (!inherits(e2, "marker_layer")) {
    cli::cli_abort(gettext(
      "Only marker_*() layers (marker_targets(), marker_context(), marker_tech()) can be added to a marker_plan."
    ))
  }
  out <- e2$apply(e1)
  out$layers <- c(out$layers, e2$label)
  out
}

#' @export
print.marker_plan <- function(x, ...) {
  n <- nrow(x$targets)
  cat(glue::glue("<marker_plan> {n} target(s); layers: {paste(x$layers, collapse = ' + ')}\n"))
  if (!is.null(x$region)) {
    cat(glue::glue(
      "  region: {x$region$chrom}:{x$region$start}-{x$region$end}\n"
    ))
  }
  if (n) {
    print(utils::head(x$targets[, c("label", "chrom", "pos", "ref", "alt")], 8L))
  }
  tech <- if (is.null(x$schemes)) "all designable" else paste(x$schemes, collapse = ", ")
  cat(glue::glue("  context: {if (is.null(x$variants)) 'none' else 'variants'}", "\n"))
  cat(glue::glue("  tech: {tech}", "\n"))
  invisible(x)
}

#' The schemes a plan may try, cheapest first
#'
#' Resolves the tech layer to an ordered character vector: the allowed schemes
#' (explicit, preset, or every designable scheme), minus any `exclude`, ordered by
#' cost ascending (ties keep registry order).
#' @keywords internal
#' @noRd
plan_schemes <- function(plan) {
  allowed <- plan$schemes
  if (is.null(allowed)) allowed <- per_locus_schemes()
  allowed <- setdiff(allowed, plan$exclude)
  unknown <- setdiff(allowed, marker_schemes()$scheme)
  if (length(unknown)) {
    cli::cli_abort(glue::glue(
      gettext("Unknown scheme(s) in the plan: {paste(unknown, collapse = ', ')}; see marker_schemes().")
    ))
  }
  costs <- plan$costs %||% KASP_TECH_COSTS
  cost <- unname(costs[allowed])
  cost[is.na(cost)] <- Inf
  reg <- marker_schemes()$scheme
  allowed[order(cost, match(allowed, reg))]
}
