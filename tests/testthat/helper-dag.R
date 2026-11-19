# =============================================================================
#  helper-dag.R  --  the shared contract check for one scheme's decision DAG
# =============================================================================
#  One line asserts the statically checkable invariants C1-C11, C14 of
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
#  (C5 floors vs KASP_LIMITS, C7 "actually evaluated"); those are checked where
#  determinable and left to the scheme agent's explicit assertions, exactly as
#  dev/dag-rebuild-plan.md section 6d requires.  C13 (every per-locus scheme is
#  designable) is asserted in test-decision.R, where the ladder is exercised.
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

# The soft gates whose declared floor is legitimately the placeholder 0: the
# footprint margin depends on the limits/backend, so grade_assay_pairs() fills it
# in at run time (R/marker_select.R) and the static constant must stay 0.
DAG_FLOOR_ZERO_OK <- "gate_com_footprint"

# The numeric `KASP_GATE_FLOOR_*` constant each scheme owns, spelled out rather
# than derived from the scheme id: the two disagree often enough (snp -> _PAIR,
# asq -> _AS, indel_size -> _INDEL_SIZE) that a derived name matches nothing and
# a wholly mis-keyed constant would be filtered out instead of surfacing as
# `stray`.  The `_LABEL*` constants are character and are not floors.  The C8
# block checks that every value exists in the namespace and every key is a
# registered scheme.
DAG_SCHEME_FLOOR_CONSTANTS <- c(
  snp = "KASP_GATE_FLOOR_PAIR",
  asq = "KASP_GATE_FLOOR_AS",
  indel_size = "KASP_GATE_FLOOR_INDEL_SIZE",
  duplex = "KASP_GATE_FLOOR_DUPLEX",
  lax = "KASP_GATE_FLOOR_LAX",
  arms = "KASP_GATE_FLOOR_ARMS",
  caps = "KASP_GATE_FLOOR_CAPS"
)

# The relaxed floors that apply to one scheme, keyed by gate id.  The scheme's
# own sources are always read - the constant DAG_SCHEME_FLOOR_CONSTANTS names for
# it (plus a derived `KASP_GATE_FLOOR_<SCHEME>` fallback that matches nothing
# today) and the `gate_floor` map in its discovered spec - and so is every other
# `KASP_GATE_FLOOR_*` numeric vector that names at least one of the scheme's
# gates (a shared constant applies to whoever owns that gate; a per-scheme
# constant must not be read by a scheme that does not own it).  Reading the own
# sources even when every key is wrong is what makes a mistyped key visible as
# `stray` instead of silently filtered out.
# `stray` names declared keys that are not this scheme's gates, `conflict` names
# keys two sources disagree about, and `counts` says how many sources declare
# each key.
dag_scheme_floors <- function(scheme, gate_ids) {
  ns <- asNamespace("tidymarkers")
  nms <- grep("^KASP_GATE_FLOOR_", ls(ns, all.names = TRUE), value = TRUE)
  raw <- stats::setNames(lapply(nms, get, envir = ns), nms)
  is_keyed_num <- function(x) {
    is.numeric(x) && length(x) > 0L && length(names(x)) == length(x)
  }
  srcs <- c(raw, lapply(extra_scheme_defs(), `[[`, "gate_floor"))
  srcs <- srcs[vapply(srcs, is_keyed_num, logical(1))]
  own <- c(
    unname(DAG_SCHEME_FLOOR_CONSTANTS[scheme]),
    paste0("KASP_GATE_FLOOR_", gsub("[^A-Za-z0-9]", "", toupper(scheme))),
    scheme
  )
  own <- own[!is.na(own)]
  shared <- vapply(srcs, function(x) any(names(x) %in% gate_ids), logical(1))
  srcs <- srcs[shared | names(srcs) %in% own]
  empty <- list(
    values = stats::setNames(numeric(0), character(0)),
    stray = character(0), conflict = character(0),
    counts = stats::setNames(integer(0), character(0))
  )
  if (!length(srcs)) {
    return(empty)
  }
  flat <- unlist(unname(srcs), use.names = TRUE)
  stray <- setdiff(unique(names(flat)), gate_ids)
  flat <- flat[names(flat) %in% gate_ids]
  # `counts` is the number of *sources* declaring a key (a repeat is a copy),
  # while `conflict` is the number of distinct values (a disagreement).
  counts <- table(names(flat))
  spread <- tapply(unname(flat), names(flat), function(v) length(unique(v)))
  list(
    values = flat[!duplicated(names(flat))],
    stray = stray,
    conflict = names(spread)[spread > 1L],
    counts = counts
  )
}

#' Assert the statically checkable DAG invariants for one scheme
#'
#' Checks C1 (one criterion, one node), C2 (family-correct check membership),
#' C4 (a failed hard gate dooms), C6 (both orders declared), C7 (declared gates
#' are emitted), C8 (hardness justified, soft gates floored), C9 (no `[lo, hi]`
#' placeholder), C10 (the ungated list is honest), C11 (acyclic, strictly
#' downward, every endpoint known, a verdict) and C14 (a tier is not a hard gate).
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
    # the family constants hold criterion ids; the labels are their names
    unname(names(MARKER_FAMILY_CHECKS[[info$family]])) %||% character(0)
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
    expect_true(any(grepl(":doomed$", n$name)), info = at("hard gate but no DOOMED"))
    for (st in stages) {
      hard <- st$gates[!(st$gates %in% st$relax)]
      if (!length(hard)) next
      from <- c(paste0(st$id, ":gate:", hard), paste0(st$id, ":decision"))
      expect_true(
        any(e$kind == "doomed" & e$from %in% from & e$to == paste0(st$id, ":doomed")),
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

  # ---- C8: hardness is justified ---------------------------------------------
  # a hard gate must say why it cannot be relaxed; a soft gate must declare a
  # floor, and a numeric floor is finite and strictly negative (a floor of 0
  # would mean "never relaxed", i.e. hard) unless it is the run-time placeholder
  # named in DAG_FLOOR_ZERO_OK.  A scheme with no soft gate passes trivially.
  gate_nodes <- n[n$kind == "gate", ]
  gate_ids <- sub("^[^:]+:gate:", "", gate_nodes$name)
  hard_ids <- gate_ids[gate_nodes$hard]
  soft_ids <- gate_ids[!gate_nodes$hard]
  hard_lbl <- scheme_label("hard", hard_ids)
  bad <- is.na(hard_lbl) | !nzchar(hard_lbl) # is.na first: nzchar(NA) is TRUE
  expect_false(any(bad), info = at(paste0(
    "hard gate without a reason: ", paste(hard_ids[bad], collapse = "; ")
  )))
  # the floor map must name real constants and real schemes: a typo would turn
  # the `own` safety net into a name that reads nothing, so a wholly mis-keyed
  # constant would be dropped instead of surfacing as `stray`
  ns_names <- ls(asNamespace("tidymarkers"), all.names = TRUE)
  missing_const <- DAG_SCHEME_FLOOR_CONSTANTS[
    !(DAG_SCHEME_FLOOR_CONSTANTS %in% ns_names)
  ]
  expect_true(!length(missing_const), info = at(paste0(
    "DAG_SCHEME_FLOOR_CONSTANTS names a missing constant: ",
    paste(unname(missing_const), collapse = "; ")
  )))
  expect_true(
    all(names(DAG_SCHEME_FLOOR_CONSTANTS) %in% marker_schemes()$scheme),
    info = at("DAG_SCHEME_FLOOR_CONSTANTS names an unknown scheme")
  )
  fl <- dag_scheme_floors(scheme, gate_ids)
  # a floor keyed to a gate the scheme does not own is a copy from another
  # scheme: it would relax a gate this DAG never evaluates
  expect_true(!length(fl$stray), info = at(paste0(
    "floors for gates the scheme does not own: ", paste(fl$stray, collapse = "; ")
  )))
  # two applicable maps disagreeing about one gate is a bug: the ladder would
  # relax to whichever constant the design path happened to read
  expect_true(!length(fl$conflict), info = at(paste0(
    "conflicting floors for: ", paste(fl$conflict, collapse = "; ")
  )))
  # one declaration per soft gate (a hard gate repeats the 0 convention in each
  # stage's map, so only soft ids are tested here)
  soft_multi <- intersect(names(fl$counts)[fl$counts > 1L], soft_ids)
  expect_true(!length(soft_multi), info = at(paste0(
    "soft gate declared a floor twice: ", paste(soft_multi, collapse = "; ")
  )))
  soft_lbl <- scheme_label("floor", soft_ids)
  for (i in seq_along(soft_ids)) {
    gid <- soft_ids[i]
    hit <- match(gid, names(fl$values))
    has_lbl <- !is.na(soft_lbl[i]) && nzchar(soft_lbl[i])
    expect_true(
      !is.na(hit) || has_lbl,
      info = at(paste0("soft gate '", gid, "' declares no floor"))
    )
    if (!is.na(hit)) {
      margin <- unname(fl$values[hit])
      expect_true(
        is.finite(margin) && (margin < 0 || gid %in% DAG_FLOOR_ZERO_OK),
        info = at(paste0(
          "soft gate '", gid, "' floor must be finite and strictly negative"
        ))
      )
    }
  }

  # ---- C14: a tier label never collides with a hard gate ----------------------
  # a tier is a rung that may be traded off; restating a hard gate as a tier
  # advertises a choice the DAG does not have
  tier_lbl <- n$label[n$kind == "tier"]
  tier_lbl <- tier_lbl[!is.na(tier_lbl)]
  hard_display <- c(hard_ids, unname(hard_lbl))
  hard_display <- hard_display[!is.na(hard_display)]
  keys <- dag_norm(hard_display)
  hit <- dag_norm(tier_lbl) %in% keys
  expect_false(any(hit), info = at(paste0(
    "tier restates a hard gate: ",
    paste0(tier_lbl[hit], " ~ ", hard_display[match(dag_norm(tier_lbl[hit]), keys)],
      collapse = "; "
    )
  )))

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
# declares the same AS + pair gate set, the same two hard AS gates, and one DOOMED
# sink per killing stage (AS and pair), plus the assay's failed-check exit.
# Call it after `expect_dag_invariants`.
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
  expect_equal(sum(s$nodes$kind == "doomed"), 3L)
  expect_setequal(
    s$nodes$name[s$nodes$kind == "doomed"],
    c("as:doomed", "pair:doomed", "assay:doomed")
  )
  invisible(scheme)
}
