# The decision schema as data: every scheme's graph must mirror its declared
# stage spec, so a new gate / tier / rung shows up automatically.  These tests
# iterate over kasp_schemes() and read each scheme's spec (built-in or
# discovered), so no scheme agent edits them.

# a scheme's declared stages, in pipeline order (the generator's view)
stages_of <- function(scheme) kasp_decision_stages(scheme)

# the node names a stage's gate / tier ids map to, guarded so a gate-less or
# tier-less stage (a documentation-only ASQ, or a scheme under rework) contributes
# nothing instead of the dangling `"<stage>:gate:"` / `"<stage>:tier:"` a bare
# paste0 would yield for an empty id vector
stage_names <- function(stages, ids_of, kind) {
  unname(unlist(lapply(stages, function(st) {
    ids <- ids_of(st)
    if (length(ids)) paste0(st$id, ":", kind, ":", ids) else character(0)
  })))
}

test_that("each scheme's schema mirrors its declared stages and is a DAG", {
  for (sc in kasp_schemes()$scheme) {
    s <- kasp_decision_schema(sc)
    stages <- stages_of(sc)
    gate <- s$nodes$name[s$nodes$kind == "gate"]
    tier <- s$nodes$name[s$nodes$kind == "tier"]
    want_gate <- stage_names(stages, function(st) st$gates, "gate")
    want_tier <- stage_names(stages, function(st) names(st$tiers), "tier")
    expect_equal(sort(gate), sort(want_gate %||% character(0)), info = sc)
    expect_equal(sort(tier), sort(want_tier %||% character(0)), info = sc)
    # each stage shows exactly one refinement node iff it has a ladder, and a
    # survivors node always
    for (st in stages) {
      n_relax <- sum(s$nodes$kind == "relax" & s$nodes$stage == st$id)
      expect_equal(n_relax, if (length(st$relax)) 1L else 0L, info = sc)
      expect_true(paste0(st$id, ":survivors") %in% s$nodes$name, info = sc)
    }
    # hard gates are exactly the ones absent from the ladder, and are flagged
    hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
    soft_ids <- stage_names(stages, function(st) st$relax, "gate")
    expect_equal(sort(hard), sort(setdiff(gate, soft_ids)), info = sc)
    # one shared DOOMED terminal exactly when the scheme has stages
    expect_equal(sum(s$nodes$kind == "doomed"), if (length(stages)) 1L else 0L, info = sc)
    # consecutive stages are chained best -> next gate (or, for a gate-less
    # documentation-only stage such as ASQ, best -> its survivors node)
    if (length(stages) > 1L) {
      for (k in seq_len(length(stages) - 1L)) {
        nxt <- stages[[k + 1L]]
        from <- paste0(stages[[k]]$id, ":best")
        to <- if (length(nxt$gates)) paste0(nxt$id, ":gate:", nxt$gates[1L]) else paste0(nxt$id, ":survivors")
        expect_true(paste(from, to, sep = " -> ") %in% paste(s$edges$from, s$edges$to, sep = " -> "), info = sc)
      }
    }
  }
})

test_that("each scheme's text export names its gates, ladder and the doomed exit", {
  for (sc in kasp_schemes()$scheme) {
    txt <- kasp_decision_text(sc)
    expect_type(txt, "character")
    stages <- stages_of(sc)
    if (length(stages)) {
      expect_true(any(grepl("gates \\(hard = never relaxed\\)", txt)), info = sc)
      expect_true(any(grepl("DOOMED", txt)), info = sc)
    }
    if (any(lengths(lapply(stages, `[[`, "relax")) > 0L)) {
      expect_true(any(grepl("relaxation ladder", txt)), info = sc)
      expect_true(any(grepl("floor:", txt)), info = sc)
    }
    expect_true(any(grepl("assay acceptance checks", txt)), info = sc)
  }
})

test_that("every scheme yields a downward-flowing DAG", {
  for (sc in kasp_schemes()$scheme) {
    s <- kasp_decision_schema(sc)
    expect_true(all(s$edges$from %in% s$nodes$name), info = sc)
    expect_true(all(s$edges$to %in% s$nodes$name), info = sc)
    ord <- stats::setNames(s$nodes$order, s$nodes$name)
    expect_true(all(ord[s$edges$to] > ord[s$edges$from]), info = sc)
    expect_true("assay:verdict" %in% s$nodes$name, info = sc)
  }
})

test_that("the marker-agnostic primer DAG obeys the full DAG contract", {
  # the one scheme with no gates and no DOOMED, so it exercises the shared
  # helper (helper-dag.R) inside the gate; scheme agents add their own
  # expect_dag_invariants("<scheme>") call after fixing their DAG
  expect_dag_invariants("primer")
})

test_that("a scheme with its own stage carries it, not the KASP core", {
  dup <- kasp_decision_schema("duplex")
  expect_true(all(paste0("duplex:gate:", stages_of("duplex")$duplex$gates) %in% dup$nodes$name))
  expect_true("pair:best" %in% dup$nodes$name)
  # the ASQ stage is gate-less (documentation-only, no engine): it declares no
  # gate, emits no gate / decision node and never dooms - only its ranking tail
  # hangs off the previous stage (C7/C13)
  asq <- kasp_decision_schema("asq")
  expect_length(stages_of("asq")$asq$gates, 0L)
  expect_false(any(grepl("^asq:gate:", asq$nodes$name)))
  expect_false("asq:decision" %in% asq$nodes$name)
  expect_false(any(grepl("^asq:", asq$edges$from[asq$edges$kind == "doomed"])))
  expect_true("asq:survivors" %in% asq$nodes$name)
  expect_true(all(paste0("asq:tier:", names(KASP_TIERS_ASQ)) %in% asq$nodes$name))
  # the bare primer scheme has no allele-specific stages and no DOOMED
  pri <- kasp_decision_schema("primer")
  expect_false(any(c("as", "pair") %in% pri$nodes$stage))
  expect_false("doomed" %in% pri$nodes$name)
  expect_setequal(pri$nodes$label[pri$nodes$kind == "check"], KASP_CHECK_PRIMER)
})

test_that("an unknown scheme is rejected by name", {
  expect_error(kasp_scheme("nope"), "Unknown marker scheme")
})

test_that("the graph is a tidygraph and can be drawn", {
  g <- kasp_decision_graph()
  expect_s3_class(g, "tbl_graph")
  skip_if_not_installed("ggraph")
  expect_s3_class(plot_decision_tree(), "ggraph")
  expect_s3_class(plot_decision_tree(scheme = "duplex"), "ggraph")
})

test_that("an assay's outcome is stamped on the gates", {
  fake <- list(selection = list(relaxed = "gate_tm", violated = c("gate_tm", "gate_com_footprint"), as = "gate_wobble"))
  g <- kasp_decision_graph(fake)
  st <- function(n) tibble::as_tibble(g)[["status"]][tibble::as_tibble(g)[["name"]] == n]
  expect_equal(st("pair:gate:gate_tm"), "violated")
  expect_equal(st("pair:gate:gate_com_footprint"), "violated")
  expect_equal(st("pair:gate:gate_amplicon"), "ok")
  expect_equal(st("as:gate:gate_wobble"), "violated")
  expect_equal(st("as:gate:gate_discrimination"), "ok")
})

test_that("a discovered scheme's own gate nodes are stamped too (C13)", {
  # CAPS records its relaxation / violation by its *own* gate ids; the node name
  # is stage-prefixed, so the stamping must be per stage, not only as / pair
  fake <- list(scheme = "caps", selection = list(relaxed = "gate_caps_amplicon", violated = "gate_caps_resolvable"))
  g <- tibble::as_tibble(kasp_decision_graph(fake))
  st <- function(n) g[["status"]][g[["name"]] == n]
  expect_equal(st("caps:gate:gate_caps_amplicon"), "relaxed")
  expect_equal(st("caps:gate:gate_caps_resolvable"), "violated")
  expect_equal(st("caps:gate:gate_caps_site"), "ok")
})

test_that("site-level checks read the assay's site and linked skips", {
  p <- c(
    AS_ref = "CATGGACCTAGCATTGACGTA",
    AS_alt = "CATGGACCTAGCATTGACGTT",
    COM = "GCTTAGCCATGGTCAATGCTA"
  )
  good <- kasp_design_checks(list(
    assay = tibble::tibble(
      amplicon = 90L, as_ref = p[[1]], as_alt = p[[2]], com = p[[3]],
      com_tm = 60, as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L,
      as_n_var = 0L, as_dropout = 0, as_foot_dropout = 0, com_dropout = 0
    ),
    thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), hairpin_dg = c(-1, -1, -1)),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3), b = rep(names(p), times = 3),
      dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
    ),
    amplicons = tibble::tibble(primer = c("AS_ref", "AS_alt"), n_products = c(1L, 0L)),
    efficiency = tibble::tibble(discrimination = 0.9),
    site = kasp_site_class("A", c("T", "G")), # degenerate
    selection = list(linked_skipped = 2L)
  ))
  expect_equal(good$detail[good$criterion == "target site class"], "degenerate")
  expect_true(good$pass[good$criterion == "target site class"])
  expect_false(good$pass[good$criterion == "linked variants skipped (n)"])
  expect_true("site" %in% check_level(good$criterion))
})

# ---- the shared DAG contract, for every scheme ------------------------------
# The invariants of dev/dag-rebuild-plan.md section 3 that are visible from the
# declared data.  Each scheme agent still adds its own explicit assertions; this
# is the foundation's regression guard over all 20 at once.

test_that("every scheme's DAG obeys the static contract (C1-C11)", {
  for (sc in kasp_schemes()$scheme) {
    expect_dag_invariants(sc)
  }
})

test_that("a hard gate dooms even in a mixed stage, beside the soft ladder (C4)", {
  # the AS stage mixes a hard gate (discrimination / AS 3' window) with soft ones
  # (wobble / footprint), so with only `decision --fail--> refine` drawn a hard
  # failure looked relaxable.  Every hard gate now routes straight to DOOMED.
  for (sc in kasp_schemes()$scheme) {
    s <- kasp_decision_schema(sc)
    sources <- s$edges$from[s$edges$kind == "doomed"]
    for (st in stages_of(sc)) {
      hard <- setdiff(st$gates, st$relax)
      if (!length(hard)) next
      from <- paste0(st$id, ":gate:", hard)
      expect_true(
        any(from %in% sources) || paste0(st$id, ":decision") %in% sources,
        info = sc
      )
    }
  }
})

test_that("both declared orders are labelled per stage (C6)", {
  for (sc in kasp_schemes()$scheme) {
    s <- kasp_decision_schema(sc)
    ord <- stats::setNames(s$nodes$order, s$nodes$name)
    txt <- kasp_decision_text(sc)
    for (st in stages_of(sc)) {
      # the gate chain follows the declared order (grave / cheap first)
      idx <- match(paste0(st$id, ":gate:", st$gates), s$nodes$name)
      expect_true(all(diff(ord[idx]) > 0L), info = sc)
      # the ladder is a declared subset, tried least grave first
      expect_true(all(st$relax %in% st$gates), info = sc)
      expect_true(
        any(grepl("in application order (grave / cheap first)", txt, fixed = TRUE)),
        info = sc
      )
      if (length(st$relax)) {
        expect_true(
          any(grepl("relaxation ladder (least grave first)", txt, fixed = TRUE)),
          info = sc
        )
      }
    }
  }
})

test_that("every gate keeps its short/limit label and its floor or hard label", {
  # the label lookups are assembled from each scheme's own file (see
  # kasp_select.R), so this guards that extraction: no gate is left unnamed in
  # the tree, a soft gate names its floor, and a hard gate says why it is hard.
  # (The long `label` used by the relaxation warning is optional - it falls back
  # to the gate id - so it is not asserted here.)
  for (sc in kasp_schemes()$scheme) {
    for (st in stages_of(sc)) {
      ids <- st$gates
      expect_false(any(is.na(scheme_label("short", ids))), info = sc)
      expect_false(any(is.na(scheme_label("limit", ids))), info = sc)
      soft <- st$relax
      hard <- setdiff(ids, soft)
      if (length(soft)) expect_false(any(is.na(scheme_label("floor", soft))), info = sc)
      if (length(hard)) expect_false(any(is.na(scheme_label("hard", hard))), info = sc)
    }
  }
})

test_that("the runtime graph is the assay's own scheme (C12)", {
  # an assay carries the chemistry it was assembled as, so it is never drawn (or
  # graded) as another scheme
  dup <- tibble::as_tibble(kasp_decision_graph(list(scheme = "duplex", selection = NULL)))
  expect_true("duplex" %in% dup$stage)
  snp <- tibble::as_tibble(kasp_decision_graph(list(scheme = "snp", selection = NULL)))
  expect_false("duplex" %in% snp$stage)
  # a scheme-less object keeps the historical default
  bare <- tibble::as_tibble(kasp_decision_graph(list(selection = NULL)))
  expect_false("duplex" %in% bare$stage)
  # the assembled assay stores its scheme, so a plot/report can match it
  a <- suppressWarnings(structure(list(scheme = "indel"), class = "kasp_assay"))
  expect_identical(a$scheme, "indel")
})

test_that("the ungated-check list is honest per scheme (C10)", {
  # a criterion the selection gates is not also reported as ungated, or the
  # warning would claim the selection ignores something it enforces
  all <- names(KASP_UNGATED_STEM)
  expect_setequal(KASP_UNGATED_CHECKS(), all) # with no scheme, the full set
  for (sc in kasp_schemes()$scheme) {
    u <- KASP_UNGATED_CHECKS(sc)
    expect_true(all(u %in% all), info = sc)
    gate <- kasp_criteria(scheme = sc)
    gate <- gate[gate$role == "gate", ]
    evidence <- gsub(
      "[^a-z0-9]", "",
      tolower(paste(c(gate$id, gate$label, gate$limit), collapse = " "))
    )
    for (l in u) {
      expect_false(grepl(KASP_UNGATED_STEM[[l]], evidence, fixed = TRUE), info = paste(sc, l))
    }
  }
  # the two conflicts the audit found, and a scheme that really does not gate them
  expect_false("off-target exact primer hits (n)" %in% KASP_UNGATED_CHECKS("indel_size"))
  expect_false("inverted repeat (arms>=4, loop<=6)" %in% KASP_UNGATED_CHECKS("lamp"))
  expect_true("off-target exact primer hits (n)" %in% KASP_UNGATED_CHECKS("snp"))
})

test_that("a relaxed winner is graded at the relaxed bar, not the ideal (C5)", {
  p <- c(
    AS_ref = "CATGGACCTAGCATTGACGTA",
    AS_alt = "CATGGACCTAGCATTGACGTT",
    COM = "GCTTAGCCATGGTCAATGCTA"
  )
  mock <- function(sel) {
    list(
      assay = tibble::tibble(
        amplicon = 115L, as_ref = p[[1]], as_alt = p[[2]], com = p[[3]],
        com_tm = 60, as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L,
        as_n_var = 0L, as_dropout = 0, as_foot_dropout = 0, com_dropout = 0
      ),
      thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), hairpin_dg = c(-1, -1, -1)),
      complementarity = tibble::tibble(
        a = rep(names(p), each = 3), b = rep(names(p), times = 3),
        dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
      ),
      amplicons = tibble::tibble(primer = c("AS_ref", "AS_alt"), n_products = c(1L, 0L)),
      efficiency = tibble::tibble(discrimination = 0.9),
      selection = sel
    )
  }
  # 115 bp is outside the ideal window, so an un-relaxed winner fails that check
  plain <- kasp_design_checks(mock(NULL))
  expect_false(plain$pass[plain$criterion == "amplicon length (bp)"])
  expect_false(grepl("relaxed", plain$limit[plain$criterion == "amplicon length (bp)"]))
  # ... but the ladder may widen it by 20 bp, so a relaxed (and warned) winner
  # passes, at the bar the selection actually applied (no contradiction, C5)
  relaxed <- kasp_design_checks(mock(list(relaxed = "gate_amplicon", violated = "gate_amplicon")))
  expect_true(relaxed$pass[relaxed$criterion == "amplicon length (bp)"])
  expect_true(grepl("relaxed", relaxed$limit[relaxed$criterion == "amplicon length (bp)"]))
})
