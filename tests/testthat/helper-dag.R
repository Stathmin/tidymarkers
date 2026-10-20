# =============================================================================
#  helper-dag.R  --  the shared contract check for one scheme's decision DAG
# =============================================================================
#  One line asserts the statically checkable invariants C1-C11 of
#  dev/dag-rebuild-plan.md section 3 for a scheme, so a scheme agent tests its
#  own DAG by calling the helper inside a test whose single expression is that
#  call, with the scheme id as the argument.
#
#  It reads only the generated schema (marker_decision_schema()), the registry
#  (marker_criteria()) and the scheme's declared stage spec, so it needs neither a
#  genome nor an assembled assay, and it is cheap enough for the inner loop.
#
#  Scope: the invariants that are *statically* visible from the declared data are
#  enforced here.  A few are only partly visible without running the designer
#  (C5 floors vs KASP_LIMITS, C7 "actually evaluated", C8 hardness justified);
#  those are checked where determinable and left to the scheme agent's explicit
#  assertions, exactly as dev/dag-rebuild-plan.md section 6d requires.
# =============================================================================

# A scheme's declared stages, in pipeline order (built-in + discovered).
dag_stages <- function(scheme) marker_decision_stages(scheme)

# Normalise a criterion label for the "same criterion shown twice" test: lower
# case, keep only alphanumerics.  So the gate "COM Tm >= AS Tm" and the check of
# the same name collapse to one key, and "amplicon length" is a prefix of
# "amplicon length (bp)".
dag_norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# The ungated set for one scheme: KASP_UNGATED_CHECKS may be a plain vector, a
# per-scheme list, or a function of the scheme once C10 makes it scheme-aware.
# Guarded so a future rename of the constant skips C10 rather than erroring.
dag_ungated <- function(scheme) {
  if (!exists("KASP_UNGATED_CHECKS", inherits = TRUE)) {
    return(character(0))
  }
  u <- KASP_UNGATED_CHECKS
  if (is.function(u)) {
    u(scheme)
  } else if (is.list(u)) {
    u[[scheme]] %||% character(0)
  } else {
    u
  }
}

#' Assert the statically checkable DAG invariants for one scheme
#'
#' Checks C1 (one criterion, one node), C2 (family-correct check membership),
#' C4 (a failed hard gate dooms), C6 (both orders declared), C7 (declared gates
#' are emitted), C9 (no `[lo, hi]` placeholder), C10 (the ungated list is honest)
#' and C11 (acyclic, strictly downward, every endpoint known, a verdict).
#' @param scheme a scheme id (see [marker_schemes()]).
#' @return the scheme id, invisibly.
#' @keywords internal
#' @noRd
expect_dag_invariants <- function(scheme) {
  info <- marker_scheme(scheme)
  crit <- marker_criteria(scheme = scheme)
  s <- marker_decision_schema(scheme)
  n <- s$nodes
  e <- s$edges
  stages <- dag_stages(scheme)
  at <- function(msg) paste0(scheme, ": ", msg)

  # ---- C11: acyclic, strictly downward, endpoints known, a verdict ----------
  expect_true(all(e$from %in% n$name), info = at("edge from an unknown node"))
  expect_true(all(e$to %in% n$name), info = at("edge to an unknown node"))
  expect_equal(anyDuplicated(n$name), 0L, info = at("duplicate node name"))
  ord <- stats::setNames(n$order, n$name)
  expect_true(all(ord[e$to] > ord[e$from]), info = at("an edge is not strictly downward"))
  expect_true("assay:verdict" %in% n$name, info = at("no assay:verdict node"))

  # ---- C1: one criterion, one node ------------------------------------------
  # a criterion must not be shown both as a gate and as a check (the graph shows
  # it once); matched after normalisation so "amplicon length" vs "amplicon
  # length (bp)" is caught too
  crit_nodes <- n[n$kind %in% c("gate", "check"), ]
  keys <- dag_norm(crit_nodes$label)
  dup <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  expect_false(any(dup), info = at(paste0(
    "a criterion is shown more than once: ",
    paste(unique(crit_nodes$label[dup]), collapse = "; ")
  )))

  # ---- C2: family-correct check membership -----------------------------------
  classes <- KASP_FAMILY_CHECKS[[info$family]]
  if (is.null(classes)) classes <- character(0)
  allowed <- c(
    KASP_CHECK_PRIMER,
    if ("site" %in% classes) KASP_CHECK_SITE,
    if ("as" %in% classes) KASP_CHECK_AS,
    MARKER_FAMILY_CHECKS[[info$family]] %||% character(0)
  )
  chk <- crit$label[crit$role == "check"]
  rogue <- setdiff(chk, allowed)
  expect_length(rogue, 0L) # a scheme carries only the checks its family obeys
  # the DAG displays exactly the registry's checks for this scheme
  expect_setequal(n$label[n$kind == "check"], chk)

  # ---- C4: a failed hard gate dooms ------------------------------------------
  # a stage with a hard gate must route that failure to DOOMED straight - from the
  # stage's decision (or a hard gate itself), never only through the relax branch
  if (any(n$kind == "gate" & n$hard)) {
    expect_true("doomed" %in% n$name, info = at("hard gate but no DOOMED"))
    for (st in stages) {
      hard <- st$gates[!(st$gates %in% st$relax)]
      if (!length(hard)) next
      from <- c(paste0(st$id, ":gate:", hard), paste0(st$id, ":decision"))
      expect_true(
        any(e$kind == "doomed" & e$from %in% from),
        info = at(paste0("stage ", st$id, " hard gate not routed to DOOMED"))
      )
    }
  }

  # ---- C6: both orders are declared and labelled ------------------------------
  for (st in stages) {
    # a gate-less stage (documentation-only, e.g. ASQ) names no gate node; a bare
    # paste0 would otherwise yield the dangling `"<stage>:gate:"`
    ids <- if (length(st$gates)) paste0(st$id, ":gate:", st$gates) else character(0)
    # gate order = the spec's vector: the schema's gate orders are increasing
    expect_true(
      all(diff(n$order[match(ids, n$name)]) > 0L),
      info = at(paste0("stage ", st$id, " gate order not declared order"))
    )
    # relaxation order is a declared subset of the gates
    expect_true(
      all(st$relax %in% st$gates),
      info = at(paste0("stage ", st$id, " relax not a subset of the gates"))
    )
  }
  if (any(lengths(lapply(stages, `[[`, "relax")) > 0L)) {
    expect_true(
      any(grepl("relaxation ladder \\(least grave first\\)", marker_decision_text(scheme))),
      info = at("the relaxation ladder is not labelled")
    )
  }

  # ---- C7: declared gates are emitted (and no decorative ones) ----------------
  declared <- unname(unlist(lapply(stages, function(st) {
    if (length(st$gates)) paste0(st$id, ":gate:", st$gates) else character(0)
  })))
  if (is.null(declared)) declared <- character(0)
  expect_setequal(n$name[n$kind == "gate"], declared)

  # ---- C9: no `[lo, hi]` placeholders ----------------------------------------
  placeholder <- grepl("\\[lo, ?hi\\]", c(crit$limit, n$detail))
  expect_false(any(placeholder), info = at("a limit is a `[lo, hi]` placeholder"))

  # ---- C10: the ungated list is honest ---------------------------------------
  ungated <- intersect(dag_ungated(scheme), chk)
  hard_evidence <- dag_norm(paste(
    stats::na.omit(crit$id[crit$role == "gate" & crit$hard]),
    stats::na.omit(crit$label[crit$role == "gate" & crit$hard]),
    stats::na.omit(crit$limit[crit$role == "gate" & crit$hard])
  ))
  for (u in ungated) {
    stem <- KASP_UNGATED_STEM[[u]]
    if (is.null(stem)) next
    expect_false(
      any(grepl(stem, hard_evidence, fixed = TRUE)),
      info = at(paste0("calls '", u, "' ungated but a hard gate enforces it"))
    )
  }

  invisible(scheme)
}

# The shared KASP core DAG block: every SNP-like scheme (snp, indel, multi, ...)
# declares the same AS + pair gate set, the same two hard AS gates, and exactly
# one DOOMED terminal.  Call it after `expect_dag_invariants(scheme)`.
#' Assert a scheme uses the shared KASP core stages
#'
#' @param scheme a scheme id (see [marker_schemes()]).
#' @return the scheme id, invisibly.
#' @keywords internal
#' @noRd
expect_kasp_core_dag <- function(scheme) {
  s <- marker_decision_schema(scheme)
  expect_identical(
    s$nodes$name[s$nodes$kind == "gate"],
    c(paste0("as:gate:", KASP_GATES_AS), paste0("pair:gate:", KASP_GATES_PAIR))
  )
  hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
  expect_setequal(hard, c("as:gate:gate_discrimination", "as:gate:gate_as_window"))
  expect_equal(sum(s$nodes$kind == "doomed"), 1L)
  invisible(scheme)
}
