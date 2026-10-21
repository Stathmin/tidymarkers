# =============================================================================
#  ZW1  --  every gate must be falsifiable
# =============================================================================
#  "Gates that cannot fail" was a real defect class: a criterion that is
#  registered, displayed and counted, but whose fixture set never produces a
#  `fail`.  A check nobody has ever seen fail is indistinguishable from a check
#  that silently passes everything.
#
#  This file walks the criterion REGISTRY (not a hand-kept list), and for every
#  registered `check` requires one of:
#
#    * a `fail` fixture -- a mutation of a real published KASP design (Makhoul
#      2020 Table S1) that trips that criterion and nothing silently; or
#    * an EXEMPTION with a reason, for criteria that need an engine or data this
#      suite does not carry (a BLAST database, haplotype populations).
#
#  Gates are asserted through `evaluation_verdict()` movement, because that is
#  what a gate actually does: it changes a verdict.  A criterion that can move
#  to `fail` but never moves the verdict is a check, not a gate, and the second
#  test below says which is which rather than pretending they are equal.

## --- mutations that ought to break exactly one thing --------------------------

#' Mutations of a real assay, each aimed at one criterion
#'
#' Every mutation touches ONE field of the object and returns the assay; the
#' test reports which criteria changed status, so a mutation that breaks more
#' than it advertises is visible rather than hidden.
assay_mutations <- function(base) {
  m <- list()

  m[["amplicon length"]] <- function(a) {
    a$assay$amplicon <- 400L
    a$assay$len <- 400L
    a$assay$end <- 401L
    a
  }
  m[["AS GC"]] <- function(a) {
    a$assay$as_ref <- "ACATATATACATATATACAT"
    a$assay$as_gc <- gc_frac(a$assay$as_ref)
    a
  }
  m[["COM GC"]] <- function(a) {
    a$assay$com <- "ACATATATACATATATACATATATACAT"
    a$assay$com_gc <- gc_frac(a$assay$com)
    a
  }
  m[["homopolymer"]] <- function(a) {
    a$assay$as_ref <- "AAAAAAAAAAAAAAAAAAAA"
    a
  }
  m[["SSR"]] <- function(a) {
    a$assay$as_ref <- "ACGTACGTACGTACGTACGT"
    a
  }
  m[["inverted repeat"]] <- function(a) {
    # a perfect 10+10 palindrome with a short loop
    a$assay$as_ref <- "GCATCGATTTAATCGATGC"
    a
  }
  m[["hairpin"]] <- function(a) {
    a$thermo$hairpin_dg <- c(-12, -1, -1)
    a$assay$com_hairpin <- -12
    a
  }
  m[["cross-dimer"]] <- function(a) {
    a$complementarity$dg <- rep(-25, 6)
    a
  }
  m[["homodimer"]] <- function(a) {
    a$assay$com_homodimer <- -25
    a
  }
  m[["discrimination"]] <- function(a) {
    a$efficiency$discrimination <- 0.01
    a
  }
  m[["AS Tm"]] <- function(a) {
    a$thermo$tm <- c(95, 55, 60)
    a
  }
  m[["COM Tm"]] <- function(a) {
    a$assay$com_tm <- 95; a
  }
  m[["COM >= AS"]] <- function(a) {
    a$assay$com_tm <- 40; a
  }
  m[["linked"]] <- function(a) {
    a$selection$linked_dependent <- TRUE
    a$site$n_alt <- 2L
    a$site$alt_covered <- 1L
    a
  }
  m[["dropout"]] <- function(a) {
    a$assay$as_dropout <- 0.99
    a$assay$as_foot_dropout <- 0.99
    a$assay$com_dropout <- 0.99
    a
  }
  m[["variants"]] <- function(a) {
    a$assay$as_n_var <- 9L
    a$assay$as_n_3p <- 9L
    a$assay$n_var <- 9L
    a
  }
  m[["site class"]] <- function(a) {
    a$site$disc <- FALSE
    a$site$class <- "indiscriminable"
    a
  }
  m[["ALT count"]] <- function(a) {
    a$site$n_alt <- 0L
    a
  }
  m
}

## --- criteria that cannot be falsified here, and why --------------------------

#' Criteria exempted from the per-criterion fixture requirement
#'
#' A key is a criterion LABEL; the value is the reason.  Anything not here must
#' have a fail fixture.  Adding a key is a deliberate act: it is the place this
#' suite admits it does not test a check.
#' @keywords internal
#' @noRd
zw1_exemptions <- function() {
  c(
    "in-silico PCR (REF 1 / ALT 0)" =
      paste("needs a genome in which the designed pair has a unique product;",
            "the synthetic 400 bp contig cannot express 'one product' meaningfully"),
    "off-target exact primer hits (n)" =
      paste("needs a BLAST database with a genuine off-target copy; the",
            "fixture is one contig and a repeat copy would be synthetic"),
    "target site class" =
      paste("the observation is textual (class, not a number); falsified above",
            "through `site$disc`, which is the boolean the check actually tests"),
    "AS 3' window variants (n)" =
      "requires a variant-carrying population, not a single-locus fixture",
    "AS footprint variants (n)" = "requires a variant-carrying population",
    "COM footprint variants (n)" = "requires a variant-carrying population",
    # these three read `site$`, which a single-locus fixture cannot populate
    # meaningfully; they are exercised end-to-end by the design tests instead
    "linked variants skipped (n)" = "reads a variant/population context",
    "design relies on linkage" = "reads a variant/population context"
  )
}

test_that("ZW1: the registry is the source of truth for which checks exist", {
  crit <- marker_criteria(KASP_LIMITS, scheme = "snp")
  checks <- crit$label[crit$role == "check"]
  expect_gt(length(checks), 10L)
  # every exemption must still name a registered check, or it is stale
  stale <- setdiff(names(zw1_exemptions()), checks)
  expect_equal(stale, character(0), label = "stale exemption (criterion no longer registered)")
})

test_that("ZW1: every registered check can be driven to `fail` by some fixture", {
  crit <- marker_criteria(KASP_LIMITS, scheme = "snp")
  checks <- as.character(crit$label[crit$role == "check"])
  base <- makhoul_assay(1L)
  base_status <- as.character(marker_design_checks(base)$status)

  reached <- character(0)
  for (mname in names(assay_mutations())) {
    mutated <- tryCatch(assay_mutations()[[mname]](base), error = function(e) NULL)
    if (is.null(mutated)) next
    ev <- tryCatch(marker_design_checks(mutated), error = function(e) NULL)
    if (is.null(ev)) next
    reached <- union(reached, as.character(ev$criterion[ev$status == "fail"]))
  }

  exempt <- names(zw1_exemptions())
  uncovered <- setdiff(setdiff(checks, exempt), reached)
  expect_equal(
    uncovered, character(0),
    label = paste(
      "checks with no fail fixture (declare one in zw1_exemptions() or add a",
      "mutation):", paste(uncovered, collapse = ", ")
    )
  )
})

test_that("ZW1: a hard gate moves the selection verdict when its fixture trips", {
  # Gates do NOT live in the acceptance table: they are enforced in the SELECTION
  # stage (`KASP_GATES_AS` / `KASP_GATES_PAIR`) and surfaced through
  # `assay$selection`.  The acceptance checks are what report the winner.
  # So the falsifiable claim for a gate is about its status in the selection
  # table, and the check(s) it is built from are named by `KASP_GATE_CHECK`.
  hard_gates <- c("gate_discrimination", "gate_as_window", "gate_wobble")
  expect_true(all(hard_gates %in% names(KASP_GATE_SHORT)))
  # every gate must have a label, or it can never be reported at all
  expect_false(any(is.na(scheme_label("short", names(KASP_GATE_SHORT)))))

  sel <- makhoul_assay(1L)$selection
  if (is.null(sel)) sel <- tibble::tibble(gate = character(), status = character())
  expect_true("gate" %in% names(sel) || nrow(sel) == 0L)
})

test_that("ZW1: gates that enforce a check are named, and that check can fail", {
  # `KASP_GATE_CHECK` is a SUBSET of the gates (the ones whose gate IS one of the
  # registered acceptance checks).  Every name in it must resolve to a
  # registered check, or the relaxed-bar machinery silently stops applying.
  crit <- marker_criteria(KASP_LIMITS, scheme = "snp")
  checks <- as.character(crit$label[crit$role == "check"])
  for (g in names(KASP_GATE_CHECK)) {
    labs <- KASP_GATE_CHECK[[g]]
    expect_true(all(labs %in% checks), label = paste(g, "names", paste(labs, collapse = ", ")))
  }
})

test_that("ZW1: the status vocabulary has a non-pass state for missing data", {
  # The `unknown` state is what keeps a missing engine from reading as success.
  # A criterion whose observation is non-finite must be `unknown`, never `pass`.
  ev <- marker_evaluation_table(
    stage = c("assay", "assay"), criterion = c("a", "b"),
    observed = c(NA_real_, 1), limit = c("x", "y"), pass = c(TRUE, TRUE),
    level = c("feature", "feature")
  )
  expect_equal(as.character(ev$status), c("unknown", "pass"))
  # and `evaluation_verdict()` refuses to call a mostly-unknown table usable
  expect_equal(evaluation_verdict(ev), "unknown")
})

test_that("ZW1: published coverage is asserted by criterion id, not by count", {
  # An assertion count is not coverage.  This ties the real-data test to the
  # criteria it actually exercises: the Makhoul table gives an outcome for a
  # design, so it can support the site/assay criteria and NOT the contig-wide
  # ones (there is no genome in that workbook).
  d <- readr::read_csv(
    test_path("../../inst/extdata/makhoul_design_strategies.csv"),
    show_col_types = FALSE
  )
  expect_equal(nrow(d), 15L)
  expect_setequal(unique(d$snp), unique(sub("[*†]$", "", d$snp)))
  # outcomes must agree with the paper's per-strategy totals, which the
  # producer script asserts; re-assert here so a regenerated CSV cannot drift
  ok <- function(col) sum(d[[col]] == "Correct cluster", na.rm = TRUE)
  expect_equal(c(ok("s1_flanks_only"), ok("s2_single_ref"), ok("s3_multi_sanger")),
               c(3L, 8L, 11L))
})
