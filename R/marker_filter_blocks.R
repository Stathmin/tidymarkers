# =============================================================================
#  marker_filter_blocks.R  --  the blocks a technology is composed from
# =============================================================================
#  Roadmap Y1.1: composition is not omnivorous.  A block names one constraint a
#  technology may add (avoid long homopolymer runs, anchor the 3' end, ...), and
#  each constraint is meaningful for some chemistries only.  "Avoid repeats" is
#  the whole point of an SSR and nonsense as a filter over its target window; a
#  guide-only CRISPR assay has no primer pair whose composition or 3' end could
#  be constrained.  So every block declares `applies_to` (the chemistry families
#  it can act on), and `marker_technology()` refuses an incompatible mix with the
#  list of offending pairs instead of building a technology that cannot mean what
#  it says.
#
#  A block is *declarative*: it records a target-region / design constraint on
#  the technology.  The base designer's generic primer QC runs regardless; a
#  filter block states which of those concerns the technology claims, and the
#  region prefilter (roadmap Y2) reads the same set.  The registry is the one
#  place the "does this combination make chemical sense?" question is answered.
# =============================================================================

# The chemistry families a block may act on: every family in `marker_schemes()`
# except the guide-only CRISPR, which has no primer pair to compose.
MARKER_PRIMER_FAMILIES <- c(
  "primer", "kasp", "multiplex", "asq", "caps", "arms", "hrm", "castpcr",
  "dpcr", "indel_size", "lamp", "rpa", "sanger", "snapshot", "ssr", "taqman"
)

# A size marker's *target* is a repeat tract (SSR) or a repeat-like size
# difference (indel_size): a block that rejects repeats is self-defeating there.
MARKER_REPEAT_TARGET_FAMILIES <- c("ssr", "indel_size")

# Chemistries whose discrimination is anchored at a primer's 3' terminal base
# (the KASP / ARMS / ASP families; a CAPS cuts at a site, a guide has a seed).
MARKER_ANCHOR_FAMILIES <- c(
  "kasp", "multiplex", "asq", "arms", "castpcr", "dpcr", "snapshot", "taqman"
)

#' The blocks a technology is composed from
#'
#' A catalogue of the composable axes `marker_technology()` honours.  Each row is
#' one block: its `kind`, the `values` it takes (or, for a filter, the constraint
#' it states), the chemistry families it `applies_to`, and a note on why.
#'
#' - `chemistry`, `policy` and `designer` rows are axes the constructor already
#'   exposes as its `base`, `specificity` and `params` arguments.
#' - `filter` rows are the extra constraints a caller requests through the
#'   constructor's `blocks` argument.  A filter block is refused on a family it
#'   cannot act on (roadmap Y1.1); `applies_to` is `"all"` only for an axis that
#'   makes sense for every chemistry.
#'
#' `applies_to` names **families**, not scheme ids, so a composed technology
#' inherits its base's answer (and a discovered chemistry inherits its family's).
#' @return a tibble with `block`, `kind`, `values`, `applies_to` (a `" | "`
#'   joined family list, or `"all"`) and `note`.
#' @examples
#' marker_filter_blocks()
#' @family KASP technologies
#' @export
marker_filter_blocks <- function() {
  primer <- MARKER_PRIMER_FAMILIES
  not_repeat_target <- setdiff(primer, MARKER_REPEAT_TARGET_FAMILIES)
  join <- function(x) if (length(x)) paste(x, collapse = " | ") else "all"
  tibble::tribble(
    ~block, ~kind, ~values, ~applies_to, ~note,
    "base", "chemistry", "any scheme id (see marker_schemes())", "all",
    "the designer + DAG stage the technology extends",
    "specificity", "policy", paste(KASP_SPECIFICITY_MODES, collapse = " | "), "all",
    "how a genome-level off-target is treated (Y3)",
    "params", "designer", "any named argument of the base designer", "all",
    "e.g. as_lengths for longer primers, tm_range, gc_range",
    "readout", "policy", "the chemistry's own readout model", "all",
    "a technology may not invent a readout the chemistry does not have",
    "composition", "filter", "primer length, GC window, GC clamp", join(primer),
    "composition of the oligo pair (a guide-only CRISPR has none)",
    "hairpins", "filter", "worst hairpin dG37", join(primer),
    "secondary structure of a primer",
    "dimers", "filter", "worst homodimer / cross-dimer dG37", join(primer),
    "primer-primer structure",
    "homopolymers", "filter", "longest homopolymer run (bp)",
    join(not_repeat_target),
    "a repeat target is excluded: the run is what the marker reads",
    "repeats", "filter", "SSR / tandem repeat in the target window",
    join(not_repeat_target),
    "an SSR's target *is* a repeat, so the filter cannot apply to it",
    "three_prime_anchor", "filter", "3' terminal base / anchor distance",
    join(MARKER_ANCHOR_FAMILIES),
    "only a chemistry that discriminates at a primer's 3' end",
    "amplicon", "filter", "amplicon length window (bp)", join(primer),
    "a guide-only CRISPR assay has no amplicon of its own"
  )
}

#' Does a block act on a chemistry family?
#' @param applies_to a `marker_filter_blocks()` `applies_to` value.
#' @param family a scheme family.
#' @return `TRUE`/`FALSE`.
#' @keywords internal
#' @noRd
filter_block_applies <- function(applies_to, family) {
  identical(applies_to, "all") ||
    family %in% strsplit(applies_to, " | ", fixed = TRUE)[[1L]]
}

#' Validate the requested blocks against a chemistry family (roadmap Y1.1)
#'
#' Refuses a block the registry does not offer, and a block that cannot act on
#' `family` - with every offending pair listed, not just the first.  This is the
#' constructor's refusal point: a chemically meaningless mix is an error, not a
#' warning.
#' @param blocks the requested block ids.
#' @param family the base scheme's family.
#' @return the de-duplicated block ids, on success.
#' @keywords internal
#' @noRd
marker_technology_blocks <- function(blocks, family) {
  if (!is.character(blocks) || anyNA(blocks)) {
    cli::cli_abort(gettext("`blocks` must be a character vector of block ids."))
  }
  blocks <- unique(blocks)
  reg <- marker_filter_blocks()
  requestable <- reg$block[reg$kind == "filter"]
  unknown <- setdiff(blocks, requestable)
  if (length(unknown)) {
    cli::cli_abort(c(
      gettext("Unknown technology block{?s}: {.val {unknown}}."),
      i = gettext("Only the {.field filter} blocks of {.fn marker_filter_blocks} can be requested.")
    ))
  }
  applies <- reg$applies_to[match(blocks, reg$block)]
  bad <- blocks[!vapply(applies, filter_block_applies, logical(1), family = family)]
  if (length(bad)) {
    pairs <- stats::setNames(
      vapply(
        bad,
        function(b) gettextf("block `%s` does not apply to the `%s` family", b, family),
        character(1)
      ),
      rep("x", length(bad))
    )
    cli::cli_abort(c(
      gettext("This block/chemistry mix cannot hold."),
      pairs,
      i = gettext("See {.fn marker_filter_blocks} for each block's {.field applies_to}.")
    ))
  }
  blocks
}
