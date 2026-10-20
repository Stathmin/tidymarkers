# ASQ / PACE: documentation-only in this build.
#
# ASQ (Kalendar 2022) reads the same competitive allele-specific extension as KASP
# but places the variant at the penultimate 3' base of the allele-specific primer
# and reports through universal probes against a single common quencher.  The
# reporter/plexity are panel/format properties, so they are not selection gates;
# the one sequence-level rule is not expressible by the shared 3'-terminal
# AS-primer engine, so `asq` carries no designer - and therefore *no gate at all*.
# These tests pin that decision and the honest, gate-less DAG that follows from it
# (see the header of R/marker_asq.R).

test_that("asq is documentation-only: no engine, so it never takes a locus", {
  expect_dag_invariants("asq")
  # no designer is registered, so marker_design_at() falls through to NULL ...
  expect_null(marker_design_at(NULL, NULL, "asq"))
  expect_false(exists("design_asq_marker", inherits = TRUE))
  # ... and the per-locus rescue ladder can therefore never select asq (C13)
  expect_false("asq" %in% per_locus_schemes())
})

test_that("the asq stage declares no gate, because none is computable (C7)", {
  # a gate no code path computes is decorative.  The two panel/format pseudo-gates
  # the audit first flagged are gone, and so is the uncomputed placement "gate"
  # (the ASQ-defining rule): neither is evaluated by any engine.
  expect_length(KASP_GATES_ASQ, 0L)
  s <- marker_decision_schema("asq")
  expect_false(any(grepl("^asq:gate:", s$nodes$name)))
  expect_false(any(grepl("gate_asq_cassette|gate_asq_readout|gate_asq_asp_geometry", s$nodes$name)))
  # the method's real rule lives in the file header (Kalendar 2022); the tree
  # names the stage by what it describes and shows the scheme's own readout, but
  # no gate - a declaration-only gate would be decorative (C7)
  txt <- marker_decision_text("asq")
  expect_true(any(grepl("ASQ primer geometry", txt, fixed = TRUE)))
  expect_true(any(grepl(marker_scheme("asq")$discrimination, txt, fixed = TRUE)))
})

test_that("the gate-less asq stage emits only its ranking tail, never a DOOMED edge (C13)", {
  s <- marker_decision_schema("asq")
  # no gate node, no decision node and no doomed source from this stage
  expect_false("asq:decision" %in% s$nodes$name)
  expect_false(any(grepl("^asq:", s$edges$from[s$edges$kind == "doomed"])))
  expect_true("asq:survivors" %in% s$nodes$name)
  # the tier is kept as documentation, chained off the previous stage
  expect_true(all(paste0("asq:tier:", names(KASP_TIERS_ASQ)) %in% s$nodes$name))
  ed <- paste(s$edges$from, s$edges$to, sep = " -> ")
  expect_true("pair:best -> asq:survivors" %in% ed)
  # a placement cannot be half-relaxed, so there is no rung and no refine node
  expect_length(KASP_RELAX_ASQ, 0L)
  expect_false(any(s$nodes$kind == "relax" & s$nodes$stage == "asq"))
})

test_that("asq keeps the KASP family's checks, because its core is KASP's (C2)", {
  chk <- marker_criteria(scheme = "asq")
  chk <- chk$label[chk$role == "check"]
  expect_true(all(KASP_CHECK_PRIMER %in% chk)) # the marker-agnostic set
  expect_true(all(KASP_CHECK_AS %in% chk)) # the allele-specific core is shared
  # the AS/COM Tm and GC windows moved into the `as` class with ASQ, not away
  expect_true(all(c("AS Tm (C)", "COM Tm (C)", "AS GC", "COM GC") %in% chk))
  # the site checks follow the family too
  expect_true(all(KASP_CHECK_SITE %in% chk))
})
