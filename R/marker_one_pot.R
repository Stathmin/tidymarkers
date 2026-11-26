# =============================================================================
#  marker_one_pot.R  --  one-pot cross-chemistry multiplexing
# =============================================================================
#  A "pot" is one reaction mix several assays share, whatever their chemistry.
#  What may share a pot is declared as *data*: `one_pot_matrix()` returns one row
#  per rule (field, comparison, threshold) and the grouper evaluates those rows.
#  The thresholds therefore live in the matrix - edit a `value` and the grouping
#  changes - and no threshold is hard-coded in the placement loop.
#
#  What the matrix does not cover is the wet-lab remainder: reagent compatibility
#  (master-mix chemistry, Mg/ionic strength, enzyme) and detection-channel overlap
#  beyond the declared `channel` column are not modelled, which is why a result
#  carries `experimental = TRUE`.
#
#  shortcut: the first-fit greedy below is not an ILP optimum; it is deterministic
#  and order-stable, which a wet-lab mix needs more than the minimum bin count.
#  Replace it with a solver when the optimal packing (or a proven group count) is
#  required.
# =============================================================================

#' One-pot compatibility matrix
#'
#' The rules a set of assays must satisfy to share a single reaction mix, as
#' data: one row per rule, naming the `assays` field the rule reads, how that
#' field is compared, and the threshold it is compared against.
#' [one_pot_compatible()] and [marker_one_pot()] evaluate these rows and know no
#' threshold of their own, so editing a `value` (or removing a row) changes the
#' grouping.
#'
#' @return A tibble with one row per rule and the columns:
#' * `rule` - a stable id for the rule;
#' * `field` - the `assays` column the rule reads (`tm`, `amplicon`, `channel`);
#' * `comparison` - how the field is compared within a candidate set:
#'   `"range_max"` (max - min across the set), `"min_gap"` (smallest pairwise
#'   difference) or `"unique"` (all non-missing values distinct);
#' * `value` - the threshold, `NA` for a comparison that needs none;
#' * `source` - what the threshold rests on, or `"local policy"`.
#'
#' @details The shipped thresholds are local policy: a 3 °C ceiling on the `tm`
#'   span, a 10 bp minimum amplicon gap for assays that share a detection channel
#'   (a missing `channel` is unassigned and is treated as sharing), and one assay
#'   per channel.  No row of `inst/extdata/scientific_sources.csv` states a
#'   one-pot threshold of this kind - the multiplex KASP record
#'   (`wang2025_multikasp`) anchors primer-dimer screening, not an assay-level
#'   `tm` span - so the `source` column says `"local policy"` rather than
#'   borrowing a citation.  Which *pairs* the `min_gap` rule applies to belongs to
#'   the comparison's definition; only the threshold is data.
#'
#' @examples
#' one_pot_matrix()
#' @family multiplex
#' @export
one_pot_matrix <- function() {
  tibble::tibble(
    rule = c("tm_span", "amplicon_gap", "channel_unique"),
    field = c("tm", "amplicon", "channel"),
    comparison = c("range_max", "min_gap", "unique"),
    value = c(3, 10, NA_real_),
    source = rep("local policy", 3L)
  )
}

#' Can these assays share one reaction mix?
#'
#' Evaluates the compatibility matrix against a candidate set of assays.  Every
#' rule row is applied to the set as a whole; the matrix is the only place a
#' threshold exists, so a set cannot pass on a criterion the matrix does not
#' declare.
#'
#' @param assays A data frame or tibble with at least the columns `label`,
#'   `scheme`, `tm`, `amplicon` and `channel`.  `channel` may be `NA`
#'   (unassigned); `tm` and `amplicon` must be numeric, and a missing value
#'   counts as "not proved", never as a pass.
#' @param matrix A compatibility matrix from [one_pot_matrix()], or a tibble with
#'   the same columns.
#' @param require `"all"` (default) to return a single logical: `TRUE` only when
#'   every matrix row passes.  `"any"` returns one named logical per matrix row,
#'   i.e. that row's verdict, for diagnosing which rule rejects a set.
#' @return With `require = "all"`, one logical.  With `require = "any"`, a
#'   logical vector named by the matrix's `rule` column.
#' @examples
#' assays <- data.frame(
#'   label = c("a", "b"), scheme = c("kasp", "hrm"),
#'   tm = c(60, 60.5), amplicon = c(100, 130), channel = c("FAM", "HEX")
#' )
#' one_pot_compatible(assays)
#' one_pot_compatible(assays, require = "any")
#' @family multiplex
#' @export
one_pot_compatible <- function(assays, matrix = one_pot_matrix(), require = "all") {
  if (!is.character(require) || length(require) != 1L || !require %in% c("all", "any")) {
    abort_msg('`require` must be "all" or "any", not {.val {require}}.')
  }
  one_pot_check_assays(assays)
  one_pot_check_matrix(matrix)
  verdict <- vapply(
    seq_len(nrow(matrix)),
    function(k) one_pot_rule_pass(assays, matrix[k, ]),
    logical(1)
  )
  names(verdict) <- matrix$rule
  if (identical(require, "any")) {
    return(verdict)
  }
  all(verdict)
}

#' Group assays into one-pot multiplexes
#'
#' Packs assays into the fewest reaction mixes this greedy pass can find, under
#' the compatibility matrix.  Assays are visited in input order; each one joins
#' the first existing group that stays compatible ([one_pot_compatible()]) and
#' has fewer than `max_markers` members, otherwise it opens a new group.  The
#' pass is deterministic and order-stable: the same input gives the same groups.
#'
#' `max_markers` is the trade-off knob between "as many markers as possible per
#' multiplex" and "as few multiplex groups as possible": a smaller value yields
#' more groups (a larger value packs more per pot, and is likelier to leave an
#' assay out where the matrix forbids a pair).
#'
#' @param assays A data frame or tibble with at least the columns `label`,
#'   `scheme`, `tm`, `amplicon` and `channel` (see [one_pot_compatible()]).
#' @param matrix The compatibility matrix to enforce, default
#'   [one_pot_matrix()].  It is carried into the result, so the grouping can be
#'   re-checked against the matrix it was built with.
#' @param max_markers Maximum assays per group (a positive whole number).
#' @return An S3 `marker_one_pot` list:
#' * `groups` - tibble `group`, `label`, `scheme`, `tm`, `amplicon`, `channel`,
#'   one row per placed assay, `group` numbering the pots in creation order;
#' * `unassigned` - tibble `label`, `reason` for assays no pot can hold, with a
#'   non-empty reason naming the matrix rule(s) that failed;
#' * `matrix` - the matrix the grouping was built with;
#' * `experimental` - `TRUE`: reagent compatibility and detection-channel overlap
#'   beyond the `channel` column are not modelled, so a pot is a proposal for the
#'   bench, not a validated mix;
#' * `n_groups` - the number of groups.
#'
#' @details Every group, including a group of one, must pass the matrix on its
#'   own; an assay that cannot stand alone under the matrix (for example a `tm`
#'   that is `NA`, so the `tm` span cannot be proved) is reported in `unassigned`
#'   rather than dropped.  An assay that merely clashes with its neighbours gets
#'   its own group: that is a separate pot, not a failure.  The grouping is a
#'   greedy first fit, not an integer-programming optimum.
#'
#'   `experimental = TRUE` is the wet-lab caveat: reagent compatibility
#'   (master-mix chemistry, ionic strength, enzyme) and detection-channel overlap
#'   beyond the declared `channel` column are not modelled by this matrix.
#' @examples
#' assays <- data.frame(
#'   label = c("a", "b", "c", "d"),
#'   scheme = c("kasp", "kasp", "hrm", "caps"),
#'   tm = c(60, 60.5, 61, 59.5),
#'   amplicon = c(100, 130, 160, 190),
#'   channel = c("FAM", "HEX", "FAM", "HEX")
#' )
#' marker_one_pot(assays)
#' marker_one_pot(assays, max_markers = 2L)$n_groups
#' @family multiplex
#' @export
marker_one_pot <- function(assays, matrix = one_pot_matrix(), max_markers = 4L) {
  one_pot_check_assays(assays)
  one_pot_check_matrix(matrix)
  if (
    !is.numeric(max_markers) || length(max_markers) != 1L || is.na(max_markers) ||
      max_markers < 1 || max_markers %% 1 != 0
  ) {
    abort_msg("`max_markers` must be one positive whole number, not {.val {max_markers}}.")
  }
  # a group cannot hold more assays than exist, and the cap keeps the integer cast safe
  max_markers <- as.integer(min(max_markers, nrow(assays)))

  groups <- list()
  unassigned_label <- character(0)
  unassigned_reason <- character(0)
  for (i in seq_len(nrow(assays))) {
    verdict <- one_pot_compatible(assays[i, , drop = FALSE], matrix, require = "any")
    if (!all(verdict)) {
      unassigned_label <- c(unassigned_label, assays$label[[i]])
      unassigned_reason <- c(unassigned_reason, glue::glue(
        "cannot share any pot: {paste(names(verdict)[!verdict], collapse = ', ')} fails for the assay alone"
      ))
      next
    }
    placed <- FALSE
    for (g in seq_along(groups)) {
      if (length(groups[[g]]) >= max_markers) next
      candidate <- c(groups[[g]], i)
      if (isTRUE(one_pot_compatible(assays[candidate, , drop = FALSE], matrix))) {
        groups[[g]] <- candidate
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      groups[[length(groups) + 1L]] <- i
    }
  }

  index <- unlist(groups, use.names = FALSE)
  if (is.null(index)) index <- integer(0)
  rows <- assays[index, , drop = FALSE]
  structure(
    list(
      groups = tibble::tibble(
        group = rep(seq_along(groups), lengths(groups)),
        label = rows$label,
        scheme = rows$scheme,
        tm = rows$tm,
        amplicon = rows$amplicon,
        channel = rows$channel
      ),
      unassigned = tibble::tibble(label = unassigned_label, reason = unassigned_reason),
      matrix = matrix,
      experimental = TRUE,
      n_groups = length(groups)
    ),
    class = "marker_one_pot"
  )
}

#' @export
print.marker_one_pot <- function(x, ...) {
  cat(glue::glue(
    "<marker_one_pot> {x$n_groups} group(s); {nrow(x$groups)} assay(s); ",
    "{nrow(x$unassigned)} unassigned\n"
  ))
  if (nrow(x$groups)) {
    print(utils::head(x$groups, 12L))
  }
  if (nrow(x$unassigned)) {
    cat(glue::glue("  unassigned: {paste(x$unassigned$label, collapse = ', ')}\n"))
  }
  if (isTRUE(x$experimental)) {
    cat("  experimental: reagent compatibility and channel overlap are not modelled\n")
  }
  invisible(x)
}

# ---- internals ---------------------------------------------------------------

# The columns every caller of the matrix needs; `channel` is the only one that
# may be missing-valued, so it is not type-checked.
ONE_POT_COLUMNS <- c("label", "scheme", "tm", "amplicon", "channel")

one_pot_check_assays <- function(assays) {
  if (!is.data.frame(assays)) {
    abort_msg("`assays` must be a data frame or tibble, not {.cls {class(assays)}}.")
  }
  need_cols(
    assays, ONE_POT_COLUMNS,
    arg = "assays",
    i = gettext("One-pot grouping needs {.field {ONE_POT_COLUMNS}}.")
  )
  if (!nrow(assays)) {
    abort_msg("`assays` has no rows, so there is nothing to place in a pot.")
  }
  not_numeric <- c("tm", "amplicon")[!vapply(assays[c("tm", "amplicon")], is.numeric, logical(1))]
  if (length(not_numeric)) {
    abort_msg(
      "{.field {not_numeric}} must be numeric in `assays`, not {.cls {class(assays[[not_numeric[1L]]])}}."
    )
  }
  invisible(assays)
}

one_pot_check_matrix <- function(matrix) {
  needed <- c("rule", "field", "comparison", "value", "source")
  if (!is.data.frame(matrix) || length(setdiff(needed, names(matrix)))) {
    abort_msg(
      "`matrix` must be a one-pot compatibility matrix with {.field {needed}}.",
      i = "Build one with `one_pot_matrix()` and edit its rows."
    )
  }
  if (!nrow(matrix)) {
    abort_msg("`matrix` has no rules, so every set would pass vacuously.")
  }
  invisible(matrix)
}

one_pot_rule_pass <- function(assays, rule) {
  field <- assays[[rule$field]]
  comparison <- rule$comparison
  if (identical(comparison, "range_max")) {
    one_pot_need_numeric_field(field, rule)
    return(isTRUE(max(field) - min(field) <= one_pot_threshold(rule)))
  }
  if (identical(comparison, "min_gap")) {
    one_pot_need_numeric_field(field, rule)
    return(one_pot_pairs_pass(field, assays$channel, one_pot_threshold(rule)))
  }
  if (identical(comparison, "unique")) {
    return(anyDuplicated(field[!is.na(field)]) == 0L)
  }
  abort_msg(
    "Unknown `comparison` {.val {comparison}} on one-pot rule {.val {rule$rule}}."
  )
}

one_pot_need_numeric_field <- function(field, rule) {
  if (!is.numeric(field)) {
    abort_msg(
      "One-pot rule {.val {rule$rule}} compares {.field {rule$field}}, which is not numeric in `assays`."
    )
  }
  invisible(field)
}

one_pot_threshold <- function(rule) {
  value <- rule$value
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    abort_msg(
      "One-pot rule {.val {rule$rule}} needs one non-missing numeric `value`, not {.val {value}}."
    )
  }
  as.numeric(value)
}

# A `min_gap` rule separates pairs that a detector cannot tell apart: two assays
# on the same channel, or where either channel is `NA` (unassigned).  Pairs on
# different known channels need no size gap.  An unprovable gap (a missing
# amplicon) fails: a pot is only accepted when it is proved compatible.
one_pot_pairs_pass <- function(values, channels, gap) {
  n <- length(values)
  if (n < 2L) {
    return(TRUE)
  }
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      distinguishable <- !is.na(channels[i]) && !is.na(channels[j]) && !isTRUE(channels[i] == channels[j])
      if (distinguishable) next
      if (!isTRUE(abs(values[i] - values[j]) >= gap)) {
        return(FALSE)
      }
    }
  }
  TRUE
}
