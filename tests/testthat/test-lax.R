# Lax KASP: one ALT primer per distinct ALT haplotype, on a shared backbone.
# Exercised on a real local chunk of the SAT2 transcripts (no synthetic fixture).
#
# The scheme's *own* stage (the lax "haplotype set" stage) decides the tube, so
# these tests assert its gate set, hardness, both orders, its floors - and, per
# dev/dag-rebuild-plan.md section 6, that the designer actually evaluates it (C7)
# instead of declaring it.

lax_hp <- local({
  f <- system.file("extdata", "sat2_transcripts.fa", package = "tidymarkers")
  if (!nzchar(f)) {
    return(NULL)
  }
  haplotypes_from_alignment(read_transcripts(f), chrom = "SAT2")
})

# a 301 bp chunk around SAT2:431 (window-local position 151)
lax_chunk <- function() {
  if (is.null(lax_hp)) {
    return(NULL)
  }
  marker_chunk(lax_hp, from = 281L, to = 581L)
}

# The design is deterministic and takes tens of seconds: build each path once and
# let the assertions read the same object instead of re-running the search.
lax_design <- local({
  cache <- new.env(parent = emptyenv())
  function(with_haplotypes = TRUE) {
    key <- if (with_haplotypes) "haplotypes" else "plain"
    if (!is.null(cache[[key]])) {
      return(cache[[key]])
    }
    ch <- lax_chunk()
    snp <- tibble::tibble(chrom = "SAT2", pos = 151L, ref = "A", alt = "G")
    args <- list(chunk_genome(ch), snp, flank = 150L, cores = 1L)
    if (with_haplotypes) {
      args$haplotypes <- chunk_haplotypes(ch)
    }
    cache[[key]] <- suppressWarnings(do.call(design_lax_marker, args))
    cache[[key]]
  }
})

test_that("design_lax_marker groups matching ALT haplotypes into one primer each", {
  skip_if(is.null(lax_hp), "bundled SAT2 FASTA not found")
  skip_if_no_oligoarrayaux()
  x <- lax_design()
  expect_s3_class(x, "marker_lax")
  # every ALT haplotype is covered, but haplotypes with the same linked variants
  # collapse to a single ALT primer group
  expect_gte(length(x$alt), 5L)
  n_groups <- length(unique(x$haplotype_groups$group))
  expect_lt(n_groups, length(x$alt))
  expect_equal(nrow(x$haplotype_groups), length(x$alt))
  # the tube: one shared REF primer and one shared COM, plus the ALT primers
  expect_equal(sum(names(x$pool) == "AS_ref"), 1L)
  expect_equal(sum(names(x$pool) == "COM"), 1L)
  # *every* ALT haplotype designs: the coverage gate is hard, so none may drop out
  expect_false(anyNA(x$alt_primers))
  expect_equal(sum(!is.na(x$alt_primers)), length(x$alt))
  # the *joint* pool is screened, and it is clean here
  expect_true(is.data.frame(x$compatibility))
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(x$compatibility)))
  expect_true(all(x$compatibility$pass))
  expect_output(print(x), "marker_lax")
})

test_that("design_lax_marker leaves a plain ALT primer without haplotypes", {
  skip_if(is.null(lax_hp), "bundled SAT2 FASTA not found")
  skip_if_no_oligoarrayaux()
  x <- lax_design(with_haplotypes = FALSE)
  expect_s3_class(x, "marker_lax")
  expect_equal(x$alt, "ALT")
  expect_length(x$assays, 1L)
  expect_output(print(x), "plain ALT primer")
  # even this degenerate (single ALT primer) path runs the lax stage, so the
  # declared gates are never decoration (C7) - and nothing is at risk
  expect_equal(x$selection$coverage, 0L)
  expect_equal(x$selection$cross_priming, 0L)
  expect_length(x$selection$violated, 0L)
})

test_that("the lax scheme is registered and has its own decision stage", {
  expect_true("lax" %in% marker_schemes()$scheme)
  s <- marker_decision_schema("lax")
  expect_true("lax" %in% s$nodes$stage)
  # the scheme keeps the whole KASP family machinery (it *is* a KASP: one AS
  # primer per haplotype on a shared backbone) plus its own set-level stage
  expect_true(all(c("as", "pair") %in% s$nodes$stage))
})

test_that("lax declares its gate set, hardness, both orders and real floors", {
  # the scheme carries the KASP family stages too, so read only its own stage
  crit <- marker_criteria(scheme = "lax")
  gates <- crit[crit$role == "gate" & crit$stage == "lax", ]
  expect_setequal(gates$id, paste0("lax:gate:", KASP_GATES_LAX))
  # application order (grave / cheap first) is the declared vector
  expect_equal(gates$id[order(gates$order)], paste0("lax:gate:", KASP_GATES_LAX))
  # only coverage is hard (C8): a haplotype no ALT primer binds cannot be called
  expect_equal(gates$id[gates$hard], "lax:gate:gate_lax_coverage")
  expect_setequal(gates$id[!gates$hard], paste0("lax:gate:", KASP_RELAX_LAX))
  # both declared orders exist and differ (C6): the graver wild-type-dropout risk
  # is applied first but relaxed last
  expect_identical(
    KASP_GATES_LAX,
    c("gate_lax_coverage", "gate_lax_alt_specific", "gate_lax_cross_priming")
  )
  expect_identical(KASP_RELAX_LAX, c("gate_lax_cross_priming", "gate_lax_alt_specific"))
  expect_false(identical(KASP_RELAX_LAX, intersect(KASP_GATES_LAX, KASP_RELAX_LAX)))
  # the stage no longer restates the AS stage's discrimination gate (C14): the
  # per-primer REF/ALT discrimination is already the hard as-stage gate
  expect_false(any(grepl("discriminat", KASP_GATES_LAX)))
  # the soft floors are real and *bounded* (C8): below their threshold, never
  # unbounded, and the coverage floor is the hard 0
  expect_equal(KASP_GATE_FLOOR_LAX[["gate_lax_coverage"]], 0)
  expect_lt(KASP_GATE_FLOOR_LAX[["gate_lax_alt_specific"]], 0)
  expect_lt(KASP_GATE_FLOOR_LAX[["gate_lax_cross_priming"]], 0)
  # one tube is one candidate, so the single tier reports the tube's *quality*
  # (the resolvable haplotype frequency); it must not restate the hard gate (C14)
  expect_setequal(names(scheme_def_lax()$tiers), "tier1")
  expect_equal(unname(KASP_TIERS_LAX), "resolvable ALT haplotype frequency")
})

test_that("lax's DAG shows the hard coverage gate, a two-rung ladder and its tier", {
  s <- marker_decision_schema("lax")
  n <- s$nodes
  lax_gate <- n$kind == "gate" & n$stage == "lax"
  expect_equal(n$name[lax_gate & n$hard], "lax:gate:gate_lax_coverage")
  # a failed hard gate dooms directly, beside the soft ladder (C4)
  ed <- paste(s$edges$from, s$edges$kind, s$edges$to, sep = " ")
  expect_true("lax:gate:gate_lax_coverage doomed doomed" %in% ed)
  # one refinement node for the two soft gates, a survivors/best tail, one tier
  expect_equal(sum(n$kind == "relax" & n$stage == "lax"), 1L)
  expect_equal(sum(n$kind == "tier" & n$stage == "lax"), 1L)
  expect_true(all(c("lax:survivors", "lax:best") %in% n$name))
  # the ladder is labelled least grave first (C6)
  txt <- marker_decision_text("lax")
  expect_true(any(grepl("haplotype cross-priming", txt, fixed = TRUE)))
  expect_true(any(grepl("relaxation ladder \\(least grave first\\)", txt)))
})

test_that("design_lax_marker evaluates the lax gates and records them (C7)", {
  skip_if(is.null(lax_hp), "bundled SAT2 FASTA not found")
  skip_if_no_oligoarrayaux()
  x <- lax_design()
  expect_setequal(
    names(x$selection),
    c(
      "relaxed", "violated", "doomed", "coverage", "alt_specific", "cross_priming",
      "resolvable", "as", "pair"
    )
  )
  # the recorded lax gate ids are the scheme's own (never another stage's)
  expect_true(all(c(x$selection$relaxed, x$selection$violated) %in% KASP_GATES_LAX))
  # the defining promise holds here: every ALT haplotype is covered (hard gate)
  # and both baked linked alleles are ALT-specific
  expect_equal(x$selection$coverage, 0L)
  expect_equal(x$selection$alt_specific, 0)
  expect_false(x$selection$doomed)
  # this target's two groups differ only 6 nt from the 3' end - *outside* the
  # 5-nt 3'-proximal window - so their ALT primers cross-prime: the soft gate is
  # relaxed to its bounded floor and reported, not hidden
  expect_equal(x$selection$cross_priming, 1L)
  expect_equal(x$selection$violated, "gate_lax_cross_priming")
  # ... so the tube resolves none of the ALT haplotypes by group (the tier margin)
  expect_equal(x$selection$resolvable, 0)
  expect_equal(x$verdict, "warn")
  expect_output(print(x), "lax stage")
})

test_that("lax_ref_leakage catches a linked allele the wild type also carries", {
  snp <- tibble::tibble(chrom = "Chr1", pos = 10L, ref = "A", alt = "G")
  # ref1 and alt1 both carry T at 20, so that "linked" allele is not ALT-specific
  h <- tibble::tibble(
    haplotype = c("ref1", "ref2", "alt1", "alt1", "ref1"),
    freq = c(0.5, 0.3, 0.2, 0.2, 0.5),
    chrom = "Chr1",
    pos = c(30L, 40L, 10L, 20L, 20L),
    ref = c("A", "G", "A", "C", "C"),
    alt = c("T", "A", "G", "T", "T")
  )
  # nothing baked, and a position no ALT carrier holds: nothing to leak
  expect_equal(lax_ref_leakage(h, snp, integer(0)), 0)
  expect_equal(lax_ref_leakage(h, snp, 40L), 0)
  # the whole wild type that carries the baked allele at 20 is dropped out
  expect_equal(lax_ref_leakage(h, snp, 20L), 0.5)
})

test_that("lax's set-level metrics count shared windows and weigh them by frequency", {
  snp <- tibble::tibble(chrom = "Chr1", pos = 100L, ref = "A", alt = "G")
  h <- tibble::tibble(
    haplotype = c("a", "a", "b", "c", "c"),
    freq = c(0.4, 0.4, 0.3, 0.3, 0.3),
    chrom = "Chr1",
    pos = c(100L, 120L, 100L, 100L, 103L),
    ref = c("A", "A", "A", "A", "C"),
    alt = c("G", "T", "G", "G", "T")
  )
  # a design stub: only `rep` and its AS 3'-proximal window are read
  stub <- function(rep) {
    list(
      rep = rep, fp = c(101L, 130L),
      design = list(assay = tibble::tibble(win3p_start = 101L, win3p_end = 105L))
    )
  }
  cover <- function(...) {
    list(
      designs = list(...),
      groups = tibble::tibble(haplotype = c("a", "b", "c"), group = c("", "", "103:C:T"))
    )
  }
  # `a`'s linked variant is at 120 and `b` has none, so their windows are
  # identical (one collapsible pair); `c`'s 103:T is inside the window, so it is
  # separable from both
  expect_equal(lax_cross_priming_pairs(h, snp, list(stub("a"), stub("b"), stub("c"))), 1L)
  expect_equal(lax_cross_priming_pairs(h, snp, list(stub("a"), stub("c"))), 0L)
  expect_equal(lax_cross_priming_pairs(h, snp, list(stub("a"))), 0L)
  # only `c`'s window is unique, so only `c`'s 30 % of the population is resolvable
  expect_equal(lax_resolvable_freq(h, snp, cover(stub("a"), stub("b"), stub("c"))), 0.3)
  # with `b` gone the other two separate, so the whole ALT population resolves
  expect_equal(lax_resolvable_freq(h, snp, cover(stub("a"), stub("c"))), 1)
  expect_equal(lax_resolvable_freq(h, snp, list(designs = list())), 0)
})

test_that("the lax stage dooms on its hard gate and on an exhausted floor", {
  clear <- suppressWarnings(lax_stage_pick())
  expect_false(is.null(clear$pick$best))
  expect_length(clear$pick$violated, 0L)
  # a haplotype no design covers cannot be called: the gate is hard (C8)
  expect_null(suppressWarnings(lax_stage_pick(n_uncovered = 1L))$pick$best)
  # a soft gate relaxes only as far as its bounded floor, then dooms
  soft <- suppressWarnings(lax_stage_pick(n_cross = KASP_LAX_CROSS_PRIMING_FLOOR))
  expect_false(is.null(soft$pick$best))
  expect_equal(soft$pick$violated, "gate_lax_cross_priming")
  expect_null(suppressWarnings(lax_stage_pick(n_cross = KASP_LAX_CROSS_PRIMING_FLOOR + 1L))$pick$best)
  leak <- suppressWarnings(lax_stage_pick(leakage = 2 * KASP_LAX_ALT_SPECIFIC_TOL))
  expect_false(is.null(leak$pick$best))
  expect_equal(leak$pick$violated, "gate_lax_alt_specific")
  expect_null(
    suppressWarnings(lax_stage_pick(leakage = 2 * KASP_LAX_ALT_SPECIFIC_TOL + 0.05))$pick$best
  )
})

test_that("lax DAG obeys the contract", expect_dag_invariants("lax"))
