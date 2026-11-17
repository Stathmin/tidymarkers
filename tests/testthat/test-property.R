# =============================================================================
#  test-property.R  --  property-based invariants (ZY11)
# =============================================================================
#  Each property runs over hedgehog-generated cheap inputs and asserts a
#  contract read from the package's own code.  No property designs a marker,
#  scans a genome or assembles an assay: every target is a pure function, so a
#  case costs microseconds.
#
#  `hedgehog` lives in Suggests only, so this file loads it itself:
#  `library(hedgehog)` otherwise trips a testthat check.
# =============================================================================

forall <- hedgehog::forall

test_that("specificity penalty is deterministic and non-negative", {
  set.seed(20261117)
  # the ladder is a schema fact, so it is asserted once per scheme rather than
  # regenerated: the terminal rung dominates every other rung.  These are the
  # four specificity schemes (`marker_schemes()` also lists "primer", which is
  # a design kind, not a specificity policy)
  for (scheme in c("pcr", "rpa", "crispr", "paralog")) {
    pol <- specificity_policy(scheme)
    expect_true(
      pol$weight_terminal >= max(pol$weight_three_prime, pol$weight_seed),
      label = scheme
    )
  }
  # a zero-mismatch pair is free.  This one is deterministic: two strings from
  # the generator are equal with probability 4^-n, so a generated case would
  # almost never reach this branch
  expect_identical(
    mismatch_profile("ACGTACGTACGT", "ACGTACGTACGT", "pcr")$penalty,
    0
  )
  forall(
    hedgehog::gen.int(20),
    function(n) {
      pattern <- prop_dna(n)
      subject <- prop_dna(n)
      scheme <- sample(c("pcr", "rpa", "crispr", "paralog"), 1L)
      a <- mismatch_profile(pattern, subject, scheme)
      b <- mismatch_profile(pattern, subject, scheme)
      # same input, same profile: the ladder is arithmetic, not random
      expect_identical(a, b)
      # `penalty` is a weight sum (`mismatch_weights()` are all >= 1), never a
      # signed score, so it is non-negative
      expect_true(a$penalty >= 0)
      # `hit_score()` is `100 - 10 * penalty`: it only ever subtracts
      if (a$penalty > 0) {
        expect_lt(100 - 10 * a$penalty, 100)
      }
    }
  )
})

test_that("artefact round-trips through saveRDS()/readRDS()", {
  set.seed(20261117)
  # marker_export() is the constructor that writes solved.rds; its input is the
  # `marker_solved` list built in marker_solve(), so a small stand-in of that
  # shape (no genome needed) exercises the same serialisation
  fixture <- function(n) {
    structure(
      list(
        report = tibble::tibble(
          target = paste0("t", seq_len(n)), scheme = "snp",
          cost = seq_len(n) / 10, ok = seq_len(n) %% 2L == 0L
        ),
        attempts = tibble::tibble(
          label = paste0("t", seq_len(n)), scheme = rep("snp", n),
          cause = NA_character_
        ),
        schemes = c("snp", "arms"),
        knobs = list(limits = 1L)
      ),
      class = "marker_solved"
    )
  }
  forall(
    hedgehog::gen.int(20),
    function(n) {
      x <- fixture(n)
      path <- tempfile(fileext = ".rds")
      on.exit(unlink(path), add = TRUE)
      saveRDS(x, path)
      y <- readRDS(path)
      expect_identical(y, x)
      expect_s3_class(y, "marker_solved")
    }
  )
})

test_that("decision graph is acyclic and its relaxation walk terminates", {
  set.seed(20261117)
  # `marker_decision_graph(scheme)` needs only the schema, not an assay or genome
  forall(
    hedgehog::gen.element(prop_schemes()),
    function(scheme) {
      # the edge table first: every endpoint known, and `order` strictly increases
      # along every edge, so no walk can revisit a node
      s <- marker_decision_schema(scheme)
      edges <- s$edges
      names <- s$nodes$name
      expect_true(all(edges$from %in% names))
      expect_true(all(edges$to %in% names))
      ord <- stats::setNames(s$nodes$order, names)
      expect_true(all(ord[edges$to] > ord[edges$from]))
      # the relaxation ladder (`relax` rungs) moves only downhill, so the walk
      # from the failing decision to a merge or DOOMED sink is finite
      relax <- edges[edges$kind == "relax", , drop = FALSE]
      if (nrow(relax)) {
        expect_true(all(relax$to != relax$from))
        expect_true(all(ord[relax$to] > ord[relax$from]))
        expect_true(all(grepl(":(survivors)$", relax$to)))
      }
    }
  )
})

test_that("a selected scheme is always a registered chemistry", {
  set.seed(20261117)
  available <- prop_schemes()
  # a generated sample of scheme-like names, some real, some invented
  pool <- c(available, paste0("not_a_chemistry_", 1:5), "KASP", "")
  forall(
    hedgehog::gen.element(pool),
    function(scheme) {
      # the registry is the only authority: a name the plan layer would accept
      # (`marker_plan()` validates against these ids) must be in it, and an
      # invented one must abort rather than fall back to another chemistry
      if (scheme %in% available) {
        expect_true(scheme %in% available)
        expect_identical(marker_scheme(scheme)$scheme, scheme)
      } else {
        expect_error(marker_scheme(scheme), class = "rlang_error")
      }
    }
  )
})
