# The duplex engine: two loci in one tube + the compatibility stage.
#
# The compatibility rules and their provenance are in the header of
# R/marker_duplex.R (Wang 2025 Multi-KASP screen; Jiang 2023 dKASP; Meilink 2025
# mxKASP).  The gate-level tests isolate one rung at a time by passing explicit
# limits, so they need no genome; `design_duplex()` is exercised end-to-end at the
# bottom (it does need OligoArrayAux, like the rest of the KASP design).

# a minimal assembled-assay stub: the duplex stage reads only these four fields
mk_assay <- function(as_ref, as_alt, com, amplicon = 90L) {
  list(assay = tibble::tibble(as_ref = as_ref, as_alt = as_alt, com = com, amplicon = amplicon))
}

# Isolate the self-complementarity rung: hold the other three gates at passing
# bars (a wide Tm window, no inter-pair complementarity possible) so the ladder
# can only be moved by the soft cap.
ladder <- function(a1, a2, max_selfcomp = 5L) {
  duplex_selection(list(a1, a2),
    readout = "presence-absence", homologous = FALSE,
    tm_window = c(0, 200), min_cross_run = 100L, max_selfcomp = max_selfcomp
  )
}

# 20-mers with known self-complementarity (checked in R): run 6 / run 8, neither
# folding back at the 3' end, and a 3'-anchored 4 bp fold-back at run 4.
P_RUN6 <- "GGTAAGGAGATCTGCCCTTG"
P_RUN8 <- "ATTCGGTGACGCGTCCAGGT"
P_FOLD <- "GATCGTACCAAGTTCCGAT"
P_CLEAN <- "TGCACCATGCATCGATCAG"

test_that("duplex_selfcomp applies both halves of the multiplex self-comp rule", {
  # Wang 2025 rule 1, first half: the longest self-complementary run is capped
  expect_equal(duplex_selfcomp(P_RUN6)$run, 6L)
  expect_false(duplex_selfcomp(P_RUN6)$fold)
  expect_equal(duplex_selfcomp(P_RUN8)$run, 8L)
  # ... and the second half is a *separate* test: a 4 bp 3'-anchored fold-back
  # fails even though 4 <= 5
  expect_equal(duplex_selfcomp(P_FOLD)$run, 4L)
  expect_true(duplex_selfcomp(P_FOLD)$fold)
})

test_that("the duplex gates are the four computed ones, with no amplicon pseudo-gate", {
  a <- mk_assay(P_RUN6, P_RUN6, P_CLEAN)
  cc <- duplex_compatibility(list(a, a),
    readout = "presence-absence", homologous = FALSE,
    tm_window = c(0, 200), min_cross_run = 100L
  )
  expect_setequal(cc$gate, KASP_GATES_DUPLEX)
  # the multiplex <= 150 bp budget is *not* a stage gate: the per-locus KASP
  # amplicon window [50, 100] is stricter, so the budget lives in a tier instead
  expect_false(any(grepl("amplicon", cc$gate)))
  expect_gte(DUPLEX_MAX_AMPLICON, KASP_LIMITS$amplicon[2L])
})

test_that("the 3' fold-back fails gate_duplex_selfcomp inside the five-base cap", {
  a <- mk_assay(P_FOLD, P_FOLD, P_CLEAN)
  cc <- duplex_compatibility(list(a, a),
    readout = "presence-absence", homologous = FALSE,
    tm_window = c(0, 200), min_cross_run = 100L, max_selfcomp = 5L
  )
  expect_false(cc$pass[cc$gate == "gate_duplex_selfcomp"])
  expect_match(cc$detail[cc$gate == "gate_duplex_selfcomp"], "fold")
})

test_that("a hard duplex gate dooms (C4); a soft one relaxes within a bounded floor", {
  clean <- mk_assay(P_RUN6, P_RUN6, P_CLEAN)
  # the soft self-complementarity rung: run 6 against a cap of 5 is margin -1,
  # inside the floor (-2), so the ladder admits it *and* reports it as violated
  soft <- ladder(clean, clean)
  expect_false(is.null(soft$pick$best))
  expect_identical(soft$pick$violated, "gate_duplex_selfcomp")
  expect_identical(soft$pick$relaxed, "gate_duplex_selfcomp")
  # ... but the floor bounds it: run 8 is margin -3, past -2, so the pairing dooms
  beyond <- mk_assay(P_RUN8, P_RUN8, P_CLEAN)
  expect_null(ladder(beyond, beyond)$pick$best)

  # a *hard* gate dooms immediately, even with the soft cap held passing:
  # a co-dominant (shared-cassette dKASP) readout needs homologous loci
  hard <- duplex_selection(list(clean, clean),
    readout = "co-dominant", homologous = FALSE,
    tm_window = c(0, 200), min_cross_run = 100L, max_selfcomp = 20L
  )
  expect_false(hard$metrics$tgt_pass)
  expect_null(hard$pick$best)

  # ... and so does inter-pair cross-talk (Wang 2025 rule 2): locus 2's common
  # primer is the reverse complement of locus 1's, a full-length hetero-dimer
  rcp <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(P_CLEAN)))
  x1 <- mk_assay(P_RUN6, P_RUN6, P_CLEAN)
  x2 <- mk_assay(P_RUN6, P_RUN6, rcp)
  cross <- duplex_selection(list(x1, x2),
    readout = "presence-absence", homologous = FALSE,
    tm_window = c(0, 200), max_selfcomp = 20L
  )
  expect_gt(cross$metrics$cross, DUPLEX_MIN_CROSS_RUN)
  expect_null(cross$pick$best)
})

test_that("duplex_compatibility flags inter-locus complementarity and the topology", {
  a1 <- mk_assay("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGA", "TTTTGCGCGCGCGCGCGCAA")
  # locus 2's `com` is the reverse complement of locus 1's -> a long complementary run
  a2 <- mk_assay("GGGGCCCCGGGGCCCCGGGG", "GGGGCCCCGGGGCCCCGGGA", "TTGCGCGCGCGCGCGCAAAA")
  cc <- duplex_compatibility(list(a1, a2), readout = "co-dominant", homologous = TRUE)
  expect_setequal(cc$gate, KASP_GATES_DUPLEX)
  expect_true(cc$pass[cc$gate == "gate_duplex_targets"]) # co-dominant needs homologous loci
  expect_false(cc$pass[cc$gate == "gate_duplex_complement"]) # the revcomp pair
  # a co-dominant call on non-homologous loci is a targets-gate failure
  bad <- duplex_compatibility(list(a1, a2), readout = "co-dominant", homologous = FALSE)
  expect_false(bad$pass[bad$gate == "gate_duplex_targets"])
  # presence/absence accepts non-homologous loci
  pa <- duplex_compatibility(list(a1, a2), readout = "presence-absence", homologous = FALSE)
  expect_true(pa$pass[pa$gate == "gate_duplex_targets"])
})

test_that("the duplex DAG obeys the contract and shows the real, computed gates", {
  expect_dag_invariants("duplex")
  s <- marker_decision_schema("duplex")
  expect_setequal(
    s$nodes$name[s$nodes$kind == "gate" & s$nodes$stage == "duplex"],
    paste0("duplex:gate:", KASP_GATES_DUPLEX)
  )
  # the three fatal rules are hard; only the self-complementarity cap is relaxable
  expect_identical(KASP_RELAX_DUPLEX, "gate_duplex_selfcomp")
  expect_equal(
    sum(s$nodes$hard[s$nodes$kind == "gate" & s$nodes$stage == "duplex"]),
    length(KASP_GATES_DUPLEX) - length(KASP_RELAX_DUPLEX)
  )
  # the ladder has a numeric, negative floor behind it (it used to declare none)
  expect_true(is.numeric(KASP_GATE_FLOOR_DUPLEX))
  expect_true(all(KASP_GATE_FLOOR_DUPLEX < 0))
  # gate order is the declared application order (grave / cheap first)
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(diff(ord[paste0("duplex:gate:", KASP_GATES_DUPLEX)]) > 0L))
  # every gate is named, the soft one names its floor, the hard ones say why
  hard <- setdiff(KASP_GATES_DUPLEX, KASP_RELAX_DUPLEX)
  expect_false(any(is.na(scheme_label("short", KASP_GATES_DUPLEX))))
  expect_false(any(is.na(scheme_label("limit", KASP_GATES_DUPLEX))))
  expect_false(any(is.na(scheme_label("hard", hard))))
  expect_false(any(is.na(scheme_label("floor", KASP_RELAX_DUPLEX))))
})

test_that("design_duplex assembles two loci and reports a verdict", {
  skip_if_no_oligoarrayaux()
  W <- bseq(800L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  bases <- c(substr(W, 250L, 250L), substr(W, 550L, 550L))
  loci <- tibble::tibble(
    chrom = "Chr1", pos = c(250L, 550L), ref = bases,
    alt = vapply(bases, function(b) setdiff(c("A", "C", "G", "T"), b)[1L], character(1))
  )
  d <- suppressWarnings(design_duplex(genome, loci, flank = 150L, cores = 1L))
  expect_s3_class(d, "marker_duplex")
  expect_true(d$verdict %in% c("pass", "warn", "fail"))
  expect_length(d$assays, 2L)
  expect_output(print(d), "marker_duplex")
})

test_that("design_duplex maps compatibility to pass / fail, never a hard failure to warn (C4)", {
  # a controlled per-locus design, so the verdict depends only on the duplex stage
  # (random genomes rarely yield *two* valid KASP pairs, so this is stubbed)
  stub <- function(as_ref, as_alt, com) {
    function(genome, snp, ...) {
      list(assay = tibble::tibble(as_ref = as_ref, as_alt = as_alt, com = com, amplicon = 90L))
    }
  }
  loci <- tibble::tibble(chrom = "Chr1", pos = c(1L, 2L), ref = "A", alt = "G")
  local_mocked_bindings(
    assemble_kasp_assay = stub(
      "CCGCGTAGGTGGGCGGATTTA", "CCGCGTAGGTGGGCGGATTTA", "AGTGTGTCGGGTGGGTGTTGT"
    ),
    .package = "tidymarkers"
  )
  # a clean pair (shared Tm window, no cross-talk, within the self-comp cap) passes
  expect_identical(
    design_duplex(NULL, loci, readout = "presence-absence", homologous = FALSE)$verdict,
    "pass"
  )
  # a *hard* topology failure dooms the pairing - it must read "fail", not "warn"
  # (the defect this replaces returned warn for a failed hard gate)
  hard <- design_duplex(NULL, loci, readout = "co-dominant", homologous = FALSE)
  expect_identical(hard$verdict, "fail")
  expect_null(hard$selection$pick$best)
  expect_false(hard$compatibility$pass[hard$compatibility$gate == "gate_duplex_targets"])
  # one doomed locus dooms the whole pairing (there is no half-duplex)
  local_mocked_bindings(
    assemble_kasp_assay = function(genome, snp, ...) marker_doomed("x", "as", snp, NULL),
    .package = "tidymarkers"
  )
  expect_identical(design_duplex(NULL, loci, homologous = FALSE)$verdict, "fail")
})

test_that("design_duplex rejects anything other than two loci", {
  expect_error(design_duplex(NULL, tibble::tibble(chrom = "1", pos = 1L, ref = "A", alt = "G")))
})
