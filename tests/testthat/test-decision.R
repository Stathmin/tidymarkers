# The selection policy as data: the graph must mirror the engine's constants,
# so a new gate / tier / rung shows up automatically.

test_that("the decision schema mirrors the selection constants and is a DAG", {
  s <- kasp_decision_schema()
  gate <- function(st) s$nodes$name[s$nodes$kind == "gate" & s$nodes$stage == st]
  relax <- function(st) s$nodes$name[s$nodes$kind == "relax" & s$nodes$stage == st]
  tier <- function(st) s$nodes$name[s$nodes$kind == "tier" & s$nodes$stage == st]

  expect_equal(gate("as"), paste0("as:gate:", KASP_GATES_AS))
  expect_equal(gate("pair"), paste0("pair:gate:", KASP_GATES_PAIR))
  # one refinement node per stage (the soft-gate ladder), not one node per rung
  expect_equal(relax("as"), "as:refine")
  expect_equal(relax("pair"), "pair:refine")
  expect_equal(tier("as"), paste0("as:tier:", names(KASP_TIERS_AS)))
  expect_equal(tier("pair"), paste0("pair:tier:", names(KASP_TIERS_PAIR)))

  # hard gates are exactly the ones absent from the ladder, and are flagged
  hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
  soft_ids <- c(paste0("as:gate:", KASP_RELAX_AS), paste0("pair:gate:", KASP_RELAX_PAIR))
  expect_setequal(hard, setdiff(c(gate("as"), gate("pair")), soft_ids))
  expect_true(all(c("doomed", "as:survivors", "pair:survivors") %in% s$nodes$name))
  # one shared DOOMED terminal, reached from both refinement nodes
  expect_equal(sum(s$nodes$kind == "doomed"), 1L)

  # endpoints are known, the stages are chained, and it only flows downward
  expect_true(all(s$edges$from %in% s$nodes$name))
  expect_true(all(s$edges$to %in% s$nodes$name))
  expect_true("as:best -> pair:gate:gate_structure" %in% paste(s$edges$from, s$edges$to, sep = " -> "))
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(ord[s$edges$to] > ord[s$edges$from]))
  # the refinement branch rejoins the ranking, not the outcome, and has a doomed exit
  expect_true("pair:survivors" %in% s$edges$to[s$edges$kind == "relax"])
  expect_true(any(s$edges$kind == "doomed"))
})

test_that("the decision text export names hard gates, floors and the doomed exit", {
  txt <- kasp_decision_text()
  expect_type(txt, "character")
  expect_true(any(grepl("hard = never relaxed", txt)))
  expect_true(any(grepl("floor:", txt)))
  expect_true(any(grepl("DOOMED", txt)))
  expect_true(any(grepl("relaxation ladder", txt)))
})

test_that("the graph is a tidygraph and can be drawn", {
  g <- kasp_decision_graph()
  expect_s3_class(g, "tbl_graph")
  skip_if_not_installed("ggraph")
  expect_s3_class(plot_decision_tree(), "ggraph")
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
