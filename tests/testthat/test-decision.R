# The decision schema as data: every scheme's graph must mirror its declared
# stage spec, so a new gate / tier / rung shows up automatically.  These tests
# iterate over marker_schemes() and read each scheme's spec (built-in or
# discovered), so no scheme agent edits them.

# a scheme's declared stages, in pipeline order (the generator's view)
stages_of <- function(scheme) marker_decision_stages(scheme)

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
  for (sc in marker_schemes()$scheme) {
    s <- marker_decision_schema(sc)
    stages <- stages_of(sc)
    gate <- s$nodes$name[s$nodes$kind == "gate"]
    tier <- s$nodes$name[s$nodes$kind == "tier"]
    want_gate <- stage_names(stages, function(st) st$gates, "gate")
    want_tier <- stage_names(stages, function(st) names(st$tiers), "tier")
    expect_equal(sort(gate), sort(want_gate %||% character(0)), info = sc)
    expect_equal(sort(tier), sort(want_tier %||% character(0)), info = sc)
    # each stage shows one rung per ladder step (in ladder order) iff it has a
    # ladder, and a survivors node always
    for (st in stages) {
      rungs <- s$nodes$name[s$nodes$kind == "relax" & s$nodes$stage == st$id]
      want_rung <- stage_names(list(st), function(x) x$relax, "relax")
      expect_identical(rungs, want_rung, info = sc)
      expect_true(paste0(st$id, ":survivors") %in% s$nodes$name, info = sc)
    }
    # hard gates are exactly the ones absent from the ladder, and are flagged
    hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
    soft_ids <- stage_names(stages, function(st) st$relax, "gate")
    expect_equal(sort(hard), sort(setdiff(gate, soft_ids)), info = sc)
    # one DOOMED sink per killing stage - a hard gate, or a spent ladder - so the
    # structure names the stage that killed the candidate, plus the assay's own
    # failed-check exit (`evaluation_verdict()` dooms a marker whose check fails)
    can_kill <- vapply(stages, function(st) length(st$gates) > 0L, logical(1))
    want_doomed <- if (any(can_kill)) {
      paste0(vapply(stages[can_kill], `[[`, character(1), "id"), ":doomed")
    } else {
      character(0)
    }
    if (any(marker_criteria(scheme = sc)$stage == "assay")) {
      want_doomed <- c(want_doomed, "assay:doomed")
    }
    expect_setequal(s$nodes$name[s$nodes$kind == "doomed"], want_doomed)
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
  for (sc in marker_schemes()$scheme) {
    txt <- marker_decision_text(sc)
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
  for (sc in marker_schemes()$scheme) {
    s <- marker_decision_schema(sc)
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
  dup <- marker_decision_schema("duplex")
  expect_true(all(paste0("duplex:gate:", stages_of("duplex")$duplex$gates) %in% dup$nodes$name))
  expect_true("pair:best" %in% dup$nodes$name)
  # the ASQ stage is gate-less (documentation-only, no engine): it declares no
  # gate, emits no gate / decision node and never dooms - only its ranking tail
  # hangs off the previous stage (C7/C13)
  asq <- marker_decision_schema("asq")
  expect_length(stages_of("asq")$asq$gates, 0L)
  expect_false(any(grepl("^asq:gate:", asq$nodes$name)))
  expect_false("asq:decision" %in% asq$nodes$name)
  expect_false(any(grepl("^asq:", asq$edges$from[asq$edges$kind == "doomed"])))
  expect_true("asq:survivors" %in% asq$nodes$name)
  expect_true(all(paste0("asq:tier:", names(KASP_TIERS_ASQ)) %in% asq$nodes$name))
  # the bare primer scheme has no allele-specific stages and no DOOMED
  pri <- marker_decision_schema("primer")
  expect_false(any(c("as", "pair") %in% pri$nodes$stage))
  expect_false("doomed" %in% pri$nodes$name)
  expect_setequal(pri$nodes$label[pri$nodes$kind == "check"], KASP_CHECK_PRIMER)
})

test_that("an unknown scheme is rejected by name", {
  expect_error(marker_scheme("nope"), "Unknown marker scheme")
})

test_that("the graph is a tidygraph and can be drawn", {
  g <- marker_decision_graph()
  expect_s3_class(g, "tbl_graph")
  # the picture is SVG markup: the solver places, the renderer writes ink, and
  # no graphics device is involved (so this draws headless, in tests)
  svg <- plot_decision_tree()
  expect_type(svg, "character")
  expect_true(any(grepl("^<svg ", svg)))
  expect_true(any(grepl("decision graph</title>", svg)))
})

test_that("an assay's outcome is stamped on the gates", {
  fake <- list(selection = list(relaxed = "gate_tm", violated = c("gate_tm", "gate_com_footprint"), as = "gate_wobble"))
  g <- marker_decision_graph(fake)
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
  g <- tibble::as_tibble(marker_decision_graph(fake))
  st <- function(n) g[["status"]][g[["name"]] == n]
  expect_equal(st("caps:gate:gate_caps_amplicon"), "relaxed")
  expect_equal(st("caps:gate:gate_caps_resolvable"), "violated")
  expect_equal(st("caps:gate:gate_caps_site"), "ok")
})

test_that("an unknown check is stamped unknown, never ok", {
  fake <- list(
    scheme = "snp",
    selection = list(),
    evaluation = tibble::tibble(
      criterion = c("amplicon length (bp)", "AS Tm (C)", "COM GC"),
      pass = c(FALSE, FALSE, TRUE),
      status = marker_check_status(c("unknown", "fail", "pass"))
    )
  )
  nd <- tibble::as_tibble(marker_decision_graph(fake))
  st <- function(lab) nd$status[match(lab, nd$label)]
  expect_equal(st("amplicon length (bp)"), "unknown")
  expect_equal(st("AS Tm (C)"), "violated")
  expect_equal(st("COM GC"), "ok")
  # a check with no evaluation row is neutral, not silently green
  expect_equal(st("AS GC"), "neutral")
})

test_that("a not_applicable check is neutral, not ok", {
  fake <- list(
    scheme = "snp",
    selection = list(),
    evaluation = tibble::tibble(
      criterion = "amplicon length (bp)",
      pass = FALSE,
      status = marker_check_status("not_applicable")
    )
  )
  nd <- tibble::as_tibble(marker_decision_graph(fake))
  expect_equal(nd$status[match("amplicon length (bp)", nd$label)], "neutral")
})

test_that("a legacy evaluation without a status still maps NA to unknown", {
  fake <- list(
    scheme = "snp",
    selection = list(),
    evaluation = tibble::tibble(
      criterion = c("amplicon length (bp)", "AS Tm (C)"),
      pass = c(NA, FALSE)
    )
  )
  nd <- tibble::as_tibble(marker_decision_graph(fake))
  st <- function(lab) nd$status[match(lab, nd$label)]
  expect_equal(st("amplicon length (bp)"), "unknown")
  expect_equal(st("AS Tm (C)"), "violated")
})

test_that("site-level checks read the assay's site and linked skips", {
  p <- kasp_primers()
  good <- marker_design_checks(list(
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
    site = marker_site_class("A", c("T", "G")), # degenerate
    selection = list(linked_skipped = 2L)
  ))
  expect_equal(good$detail[good$criterion == "target site class"], "degenerate")
  expect_true(good$pass[good$criterion == "target site class"])
  expect_false(good$pass[good$criterion == "linked variants skipped (n)"])
  expect_true("site" %in% marker_criteria()$level)
})

# ---- the shared DAG contract, for every scheme ------------------------------
# The invariants of dev/dag-rebuild-plan.md section 3 that are visible from the
# declared data.  Each scheme agent still adds its own explicit assertions; this
# is the foundation's regression guard over all 20 at once.

test_that("every scheme's DAG obeys the static contract (C1-C11)", {
  for (sc in marker_schemes()$scheme) {
    expect_dag_invariants(sc)
  }
})

test_that("the solver places every node once and routes every edge", {
  skip_if_not(decision_solver_available(), "the layout engine is not installed")
  for (sc in c("snp", "duplex", "primer")) {
    g <- marker_decision_graph(scheme = sc)
    o <- decision_solve(g, decision_style())
    expect_false(is.null(o), info = sc)
    # the node set survives the solver: nothing dropped, nothing duplicated
    expect_setequal(o$nodes$id, tidygraph::as_tibble(g)$name)
    expect_false(anyDuplicated(o$nodes$id) > 0L)
    expect_true(all(is.finite(o$nodes$x) & is.finite(o$nodes$y)), info = sc)

    # boxes never overlap
    for (i in seq_len(nrow(o$nodes) - 1L)) {
      for (j in (i + 1L):nrow(o$nodes)) {
        a <- o$nodes[i, ]
        b <- o$nodes[j, ]
        ix <- min(a$x + a$w, b$x + b$w) - max(a$x, b$x)
        iy <- min(a$y + a$h, b$y + b$h) - max(a$y, b$y)
        expect_false(ix > 0.01 && iy > 0.01, info = paste(sc, a$id, b$id))
      }
    }

    # every edge is routed, and no segment crosses a box it does not touch
    e <- decision_edges(g)
    expect_equal(nrow(o$edges), nrow(e))
    for (k in seq_len(nrow(o$edges))) {
      pts <- o$edges$points[[k]]
      expect_gte(nrow(pts), 2L)
      for (m in seq_len(nrow(pts) - 1L)) {
        for (j in seq_len(nrow(o$nodes))) {
          nb <- o$nodes[j, ]
          if (nb$id %in% c(o$edges$from[k], o$edges$to[k])) next
          x0 <- min(pts[m, 1], pts[m + 1L, 1])
          x1 <- max(pts[m, 1], pts[m + 1L, 1])
          y0 <- min(pts[m, 2], pts[m + 1L, 2])
          y1 <- max(pts[m, 2], pts[m + 1L, 2])
          crosses <- x1 > nb$x + 0.5 && x0 < nb$x + nb$w - 0.5 &&
            y1 > nb$y + 0.5 && y0 < nb$y + nb$h - 0.5
          expect_false(
            crosses,
            info = paste(sc, o$edges$from[k], "->", o$edges$to[k], "through", nb$id)
          )
        }
      }
    }
  }
})

test_that("the drawing is near square, not a strip", {
  # the old failure was 1:4.5 (and 1:10 with dot): one column per stage, so a
  # scheme with many checks grew sideways *and* downwards
  skip_if_not(decision_solver_available(), "the layout engine is not installed")
  for (sc in c("snp", "duplex")) {
    o <- decision_solve(marker_decision_graph(scheme = sc), decision_style())
    ar <- o$width / o$height
    expect_gt(ar, 0.3)
    expect_lt(ar, 3)
  }
})

test_that("a render never errors when the layout engine is missing", {
  withr::with_envvar(c(TIDYKASP_NO_ELK = "1"), {
    svg <- NULL
    expect_message(svg <- plot_decision_tree(scheme = "snp"), "without the layout engine")
    expect_true(any(grepl("^<svg ", svg)))
    # the fallback still names every node (wrapped inside its box, so the first
    # word is what survives intact) and draws one line per edge
    g <- marker_decision_graph(scheme = "snp")
    n <- tidygraph::as_tibble(g)
    first <- sub("[[:space:]].*$", "", n$label)
    for (k in seq_along(first)) {
      expect_true(any(grepl(first[k], svg, fixed = TRUE)), info = n$label[k])
    }
    # one <path> per edge (the arrow heads are paths too, inside <defs>)
    routed <- grep('^<path d="M', svg, value = TRUE)
    expect_equal(length(routed), nrow(decision_edges(g)))
  })
})

test_that("a hard gate dooms even in a mixed stage, beside the soft ladder (C4)", {
  # the AS stage mixes a hard gate (discrimination / AS 3' window) with soft ones
  # (wobble / footprint), so with only `decision --fail--> refine` drawn a hard
  # failure looked relaxable.  Every hard gate now routes straight to DOOMED.
  for (sc in marker_schemes()$scheme) {
    s <- marker_decision_schema(sc)
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
  for (sc in marker_schemes()$scheme) {
    s <- marker_decision_schema(sc)
    ord <- stats::setNames(s$nodes$order, s$nodes$name)
    txt <- marker_decision_text(sc)
    for (st in stages_of(sc)) {
      # the gate nodes are listed in the declared order (presentation only: the
      # selector ANDs them, it does not apply them one after another)
      idx <- match(paste0(st$id, ":gate:", st$gates), s$nodes$name)
      expect_true(all(diff(ord[idx]) > 0L), info = sc)
      # the ladder is a declared subset, tried least grave first
      expect_true(all(st$relax %in% st$gates), info = sc)
      expect_true(
        any(grepl("all applied to the same candidates (AND)", txt, fixed = TRUE)),
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
  # marker_select.R), so this guards that extraction: no gate is left unnamed in
  # the tree, a soft gate names its floor, and a hard gate says why it is hard.
  # (The long `label` used by the relaxation warning is optional - it falls back
  # to the gate id - so it is not asserted here.)
  for (sc in marker_schemes()$scheme) {
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
  dup <- tibble::as_tibble(marker_decision_graph(list(scheme = "duplex", selection = NULL)))
  expect_true("duplex" %in% dup$stage)
  snp <- tibble::as_tibble(marker_decision_graph(list(scheme = "snp", selection = NULL)))
  expect_false("duplex" %in% snp$stage)
  # a scheme-less object keeps the historical default
  bare <- tibble::as_tibble(marker_decision_graph(list(selection = NULL)))
  expect_false("duplex" %in% bare$stage)
  # the assembled assay stores its scheme, so a plot/report can match it
  a <- suppressWarnings(structure(list(scheme = "indel"), class = "marker_assay"))
  expect_identical(a$scheme, "indel")
})

test_that("the ungated-check list is honest per scheme (C10)", {
  # a criterion the selection gates is not also reported as ungated, or the
  # warning would claim the selection ignores something it enforces
  all <- names(KASP_UNGATED_STEM)
  expect_setequal(KASP_UNGATED_CHECKS(), all) # with no scheme, the full set
  for (sc in marker_schemes()$scheme) {
    u <- KASP_UNGATED_CHECKS(sc)
    expect_true(all(u %in% all), info = sc)
    gate <- marker_criteria(scheme = sc)
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
  p <- kasp_primers()
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
  plain <- marker_design_checks(mock(NULL))
  expect_false(plain$pass[plain$criterion == "amplicon length (bp)"])
  expect_false(grepl("relaxed", plain$limit[plain$criterion == "amplicon length (bp)"]))
  # ... but the ladder may widen it by 20 bp, so a relaxed (and warned) winner
  # passes, at the bar the selection actually applied (no contradiction, C5)
  relaxed <- marker_design_checks(mock(list(relaxed = "gate_amplicon", violated = "gate_amplicon")))
  expect_true(relaxed$pass[relaxed$criterion == "amplicon length (bp)"])
  expect_true(grepl("relaxed", relaxed$limit[relaxed$criterion == "amplicon length (bp)"]))
})

# ---- Z2: family membership is held as criterion ids ---------------------------
# `MARKER_CHECK_*` / `MARKER_FAMILY_CHECKS` name registry ids; the label printed
# anywhere is resolved through the registry, so the two cannot drift.
# each constant against one scheme family that carries it (`taqman`/`ssr` stand
# in for the probe / size checks, which those families share)
family_constants <- function() {
  list(
    caps = MARKER_CHECK_CAPS, lamp = MARKER_CHECK_LAMP, rpa = MARKER_CHECK_RPA,
    crispr = MARKER_CHECK_CRISPR, taqman = MARKER_CHECK_PROBE,
    sanger = MARKER_CHECK_SANGER, ssr = MARKER_CHECK_SIZE
  )
}

test_that("family checks are registry ids, resolved to labels through it", {
  spec <- criterion_spec()
  consts <- family_constants()
  for (fam in names(consts)) {
    ids <- as.vector(consts[[fam]])
    # an id, not a human label: it is a registry id and not a registry label
    expect_true(all(ids %in% spec$id), info = fam)
    expect_length(intersect(ids, spec$label), 0L)
    expect_true(all(grepl("^family_[a-z0-9_]+$", ids)), info = fam)
    # the label is the registry's own row, in the same order
    expect_equal(
      marker_family_check_labels(fam),
      spec$label[match(ids, spec$id)],
      info = fam
    )
    # every family in MARKER_FAMILY_CHECKS resolves to a non-empty label vector
    for (f in names(MARKER_FAMILY_CHECKS)) {
      expect_true(length(marker_family_check_labels(f)) > 0L, info = f)
    }
  }
  # a family that shares a vector resolves through the same table
  expect_identical(marker_family_check_labels("taqman"), marker_family_check_labels("dpcr"))
  expect_identical(marker_family_check_labels("ssr"), marker_family_check_labels("indel_size"))
  # an unknown family is empty, never an error
  expect_length(marker_family_check_labels("nope"), 0L)
})

test_that("the printed label of every family id is unchanged (no figure drift)", {
  # the committed figures and `dev/check.R` pin this text: a criterion's display
  # label must still be the one the family list used to carry, byte for byte
  labels <- c(
    "restriction site defined", "restriction site separates alleles",
    "band gap (bp)", "digest sums to amplicon",
    "restriction site conserved (population)",
    "FIP 3' ALT", "inner arm Tm balance (C)", "inner arms inside amplicon",
    "variant at FIP 3' end",
    "primer length (nt)", "THF off variant", "RPA probe spans variant",
    "variant under primer/probe",
    "PAM valid", "variant in seed", "guide spans variant",
    "guide/PAM conserved (population)", "off-target protospacers (n)",
    "probe spans variant", "probe Tm - primer Tm (C)",
    "ASP 3' discrimination", "extension 3' at variant",
    "variant-read distance (bp)", "read primers inside amplicon",
    "read primer pair defined",
    "repeat unit (nt)", "repeat copies",
    "REF/ALT size difference (bp)", "fragment size (bp)"
  )
  printed <- unique(unlist(lapply(names(MARKER_FAMILY_CHECKS), marker_family_check_labels)))
  expect_setequal(printed, labels)

  # the same label reaches the DAG and the text export for every scheme that
  # owns the check, so neither can print an id
  spec <- criterion_spec()
  for (sc in marker_schemes()$scheme) {
    fam <- marker_scheme(sc)$family
    ids <- as.vector(MARKER_FAMILY_CHECKS[[fam]] %||% character(0))
    if (!length(ids)) {
      next
    }
    labels_of <- stats::setNames(marker_family_check_labels(fam), ids)
    nodes <- marker_decision_schema(sc)$nodes
    txt <- marker_decision_text(sc)
    for (id in ids) {
      expect_true(id %in% nodes$name, info = paste(sc, id))
      expect_identical(nodes$label[match(id, nodes$name)], unname(labels_of[[id]]), info = paste(sc, id))
      # the text export prints the label, and never the id
      expect_true(any(grepl(unname(labels_of[[id]]), txt, fixed = TRUE)), info = paste(sc, id))
    }
  }
})
