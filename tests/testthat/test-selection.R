# The selection logic: gate the grave sins, then rank on tier margins.

# A minimal candidate table with the columns grade_assay_pairs()/select need.
pair_table <- function(discrimination, com_homodimer = -1, com_hairpin = -0.5,
                       com_tm = 60, amplicon = 80, n_var = 0L, as_n_3p = 0L,
                       run = 0L, ssr = 0L) {
  tibble::tibble(
    discrimination = discrimination,
    com_homodimer = com_homodimer, com_hairpin = com_hairpin,
    com_tm = com_tm, amplicon = amplicon, n_var = as.integer(n_var),
    as_n_3p = as.integer(as_n_3p), run = as.integer(run), ssr = as.integer(ssr),
    com_start = seq_along(discrimination)
  )
}

test_that("grade_assay_pairs builds gates and bottleneck tier margins", {
  d <- pair_table(c(0.9, 0.9), com_homodimer = c(-4, -8), com_hairpin = c(-1, -6))
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_true(all(c("tier1", "tier2", "m_disc") %in% names(g$data)))
  # tier1 is the bottleneck of the two structure margins
  expect_equal(
    g$data$tier1,
    pmin(c(-4, -8) - KASP_LIMITS$homodimer_dg, c(-1, -6) - KASP_LIMITS$hairpin_dg)
  )
  expect_true(all(g$data$gate_discrimination))
})

test_that("a clean pairing beats a prettier one with a COM footprint variant", {
  d <- pair_table(
    discrimination = c(0.99, 0.60),
    com_tm = c(59.95, 61.5),
    n_var = c(2L, 0L) # the first is otherwise perfect but sits on known SNPs
  )
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  pick <- select_with_ladder(g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties, g$gate_margin, g$gate_floor)
  expect_equal(pick$best$n_var, 0L)
  expect_length(pick$relaxed, 0L)
  expect_length(pick$violated, 0L)
})

test_that("the ladder relaxes only as far as needed, least grave first", {
  # every candidate sits on known SNPs: the COM-footprint gate must give way,
  # and only after the milder gates have been tried.  The two rows are otherwise
  # equal, so the least-laden pairing wins on the footprint tie-break.
  d <- pair_table(discrimination = c(0.9, 0.9), n_var = c(3L, 1L))
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  pick <- select_with_ladder(g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties, g$gate_margin, g$gate_floor)
  expect_true("gate_com_footprint" %in% pick$relaxed)
  expect_false("gate_as_window" %in% pick$relaxed) # the gravest is never reached
  expect_true("gate_com_footprint" %in% pick$violated)
  expect_equal(pick$best$n_var, 1L) # still the least-laden pairing
})

test_that("margins are compared within a tier, not summed", {
  # same tier1 bottleneck (min = 1) but very different totals -> still tied
  d <- pair_table(
    discrimination = c(0.9, 0.9),
    com_homodimer = c(-8, -8), # margin 1 both
    com_hairpin = c(-0.1, -2.0) # margin 2.9 vs 1.0; bottleneck is 1 for both
  )
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_equal(g$data$tier1[1], g$data$tier1[2])
})

test_that("the AS 3' window is an AS-stage concern, not re-gated per pair", {
  # a variant in the AS 3' window is fixed by the AS placement, so the pair
  # grader no longer gates it (it was constant across pairs there)
  expect_false("gate_as_window" %in% KASP_GATES_PAIR)
  expect_false("gate_discrimination" %in% KASP_GATES_PAIR)
  expect_true("gate_as_window" %in% KASP_GATES_AS)
  expect_true("gate_as_window" %in% setdiff(KASP_GATES_AS, KASP_RELAX_AS)) # hard
})

test_that("weighted load ranks one common variant below three rare ones", {
  d <- pair_table(discrimination = c(0.9, 0.9), n_var = c(1L, 3L))
  d$w_var <- c(0.40, 0.03) # one common SNP is the graver load here
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_equal(g$data$m_footprint, -c(0.40, 0.03))
  # the gate is still the plain count: both rows have variants under the footprint
  expect_false(any(g$data$gate_com_footprint))
})

test_that("primer_dropout is position- and frequency-aware", {
  skip_if_no_oligoarrayaux()
  p <- "CCGACGACAACGAACTCCTT"
  flip <- function(i) setdiff(c("A", "C", "G", "T"), substr(p, i, i))[1L]
  # a variant at the primer's 3' end costs far more than one at its 5' end
  d3 <- primer_dropout(p, p, TRUE, nchar(p), flip(nchar(p)), 1)
  d5 <- primer_dropout(p, p, TRUE, 1L, flip(1L), 1)
  expect_gt(d3, d5)
  # weights (population frequency) scale the expected dropout
  half <- primer_dropout(p, p, TRUE, nchar(p), flip(nchar(p)), 0.5)
  expect_equal(half, 0.5 * d3)
  # nothing to lose when there is no variant
  expect_equal(primer_dropout(p, p, TRUE, integer(0), character(0), numeric(0)), 0)
})

test_that("zone_dropouts pools candidates and matches the per-zone dropout", {
  skip_if_no_oligoarrayaux()
  tmpl <- list(seq = Biostrings::DNAString(paste(rep("ACGT", 40L), collapse = "")))
  off <- 1L
  vpos <- c(30L, 60L, 61L, 90L)
  valt <- c("T", "G", "T", "C")
  vweight <- c(1, 0.5, 1, 0.2)
  primers <- c("ACGTACGTACGTACGTACGT", "TTTTACGTACGTACGTACGT")
  starts <- c(20L, 50L)
  ends <- c(40L, 70L)
  forwards <- c(TRUE, FALSE)
  nvars <- c(2L, 2L)
  ref <- vapply(seq_along(primers), function(i) {
    zone_dropout(primers[i], starts[i], ends[i], forwards[i], nvars[i], tmpl, off, vpos, integer(0), valt, vweight)
  }, numeric(1))
  got <- zone_dropouts(primers, starts, ends, forwards, nvars, tmpl, off, vpos, integer(0), valt, vweight)
  expect_equal(got, ref)
})

test_that("the dropout gate tolerates benign footprint variants (dropout <= tol)", {
  d <- pair_table(discrimination = c(0.9, 0.9))
  d$n_var <- c(2L, 2L)
  d$com_dropout <- c(0.02, 0.60) # same count, very different cost
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_true(g$data$gate_com_footprint[1])
  expect_false(g$data$gate_com_footprint[2])
})

test_that("the pair screen swaps a top pairing with an off-target duplicate for a clean sibling", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  W <- f$W
  genome <- f$genome
  snp <- f$snp
  # the ranked winner with the screen off, so we know which COM to sabotage
  base <- suppressWarnings(
    assemble_kasp_assay(genome, snp, flank = 150L, cores = 1L, pair_screen = 0L)
  )
  com0 <- base$assay$com[1L]
  # duplicate the winner COM's exact site on a second chromosome: it now has an
  # off-target exact hit, so the pair screen must prefer a clean sibling
  genome2 <- Biostrings::DNAStringSet(c(
    Chr1 = W,
    Chr2 = paste0(strrep("A", 60L), com0, strrep("A", 60L))
  ))
  # regression: the screen's grouped scan must search a primer's reverse
  # complement too.  Here the plus strand carries only the Chr2 duplicate and the
  # reverse complement the real Chr1 site, so a plus-only batch hides one of the
  # two and the pairing then reads clean (which is the bug this test guards).
  rc0 <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(com0)))
  hits <- pattern_hits_batch(c(com0, rc0), genome2, 0L)
  hit_chroms <- function(p) names(Filter(length, hits[[p]]))
  expect_identical(hit_chroms(com0), "Chr2")
  expect_identical(hit_chroms(rc0), "Chr1")
  msgs <- character()
  a <- withCallingHandlers(
    assemble_kasp_assay(genome2, snp, flank = 150L, cores = 1L, pair_screen = 3L),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(inherits(a, "marker_doomed"))
  expect_true(any(grepl("Top-ranked pairing", msgs, fixed = TRUE)))
  expect_false(identical(a$assay$com[1L], com0))
  chk <- marker_design_checks(a)
  expect_true(chk$pass[chk$criterion == "off-target exact primer hits (n)"])
})

test_that("a strict specificity policy dooms instead of returning an off-target marker (Y3)", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  W <- f$W
  genome <- f$genome
  snp <- f$snp
  # the AS oligo of the ranked winner (screen off), then duplicate its exact site
  # on a second chromosome: every pairing now has an off-target exact hit
  base <- suppressWarnings(
    assemble_kasp_assay(genome, snp, flank = 150L, cores = 1L, pair_screen = 0L)
  )
  as0 <- base$assay$as_ref[1L]
  genome2 <- Biostrings::DNAStringSet(c(
    Chr1 = W,
    Chr2 = paste0(strrep("A", 60L), as0, strrep("A", 60L))
  ))
  # an explicit strict policy -> no clean pairing -> DOOMED with the reason
  d <- suppressWarnings(
    assemble_kasp_assay(genome2, snp,
      flank = 150L, cores = 1L, pair_screen = 3L,
      on_offtarget = "doom"
    )
  )
  expect_true(inherits(d, "marker_doomed"))
  expect_match(d$reason, "off-target")
  # the default policy is tolerant ("don't be heavy-handed"): the warned leader
  # is kept, not doomed
  a <- suppressWarnings(
    assemble_kasp_assay(genome2, snp, flank = 150L, cores = 1L, pair_screen = 3L)
  )
  expect_false(inherits(a, "marker_doomed"))
})

test_that("marker_specificity reads the scheme policy, tolerant by default (Y3)", {
  expect_identical(marker_specificity("snp"), "multicopy")
  expect_true(marker_specificity("crispr") %in% KASP_SPECIFICITY_MODES)
  expect_setequal(KASP_SPECIFICITY_STRICT, c("unique", "homeolog-aware"))
  expect_false("multicopy" %in% KASP_SPECIFICITY_STRICT)
})

test_that("the two declared orders are separate and explicit (C6)", {
  # gate order = how the method applies the gates (grave / cheap first); relax
  # order = how the ladder lowers them (least grave first).  They need not agree.
  expect_identical(KASP_GATES_AS, c("gate_discrimination", "gate_wobble", "gate_as_footprint", "gate_as_window"))
  # the wobble (the manual's WORST mismatch) is the graver of the two soft AS
  # gates, so the ladder spends the milder footprint trade-off first
  expect_identical(KASP_RELAX_AS, c("gate_as_footprint", "gate_wobble"))
  expect_identical(KASP_GATES_PAIR, c("gate_structure", "gate_tm", "gate_amplicon", "gate_com_footprint"))
  expect_identical(KASP_RELAX_PAIR, c("gate_amplicon", "gate_structure", "gate_tm", "gate_com_footprint"))
  # both are declared subsets, and the pair ladder genuinely reorders the gates
  expect_true(all(KASP_RELAX_AS %in% KASP_GATES_AS))
  expect_true(all(KASP_RELAX_PAIR %in% KASP_GATES_PAIR))
  expect_false(identical(KASP_RELAX_PAIR, intersect(KASP_GATES_PAIR, KASP_RELAX_PAIR)))
  # the gate order and the relax order differ for the AS stage too
  expect_false(identical(KASP_RELAX_AS, intersect(KASP_GATES_AS, KASP_RELAX_AS)))
})

test_that("a soft gate's acceptance check relaxes to the gate's floor (C5)", {
  # the check must read the bar the ladder may reach, or a relaxed winner (which
  # the selector deliberately accepted, with a warning) would read as a failure
  amp <- marker_gate_relaxed_bar("amplicon length (bp)", KASP_LIMITS)
  expect_true(amp$test(KASP_LIMITS$amplicon[2L] + 20)) # at the floor
  expect_false(amp$test(KASP_LIMITS$amplicon[2L] + 21)) # past it
  expect_false(amp$test(KASP_LIMITS$amplicon[1L] - 21))
  tm <- marker_gate_relaxed_bar("COM Tm >= AS Tm", KASP_LIMITS)
  expect_true(tm$test(-2))
  expect_false(tm$test(-3))
  hair <- marker_gate_relaxed_bar("worst hairpin dG37", KASP_LIMITS)
  homo <- marker_gate_relaxed_bar("worst homodimer dG37", KASP_LIMITS)
  expect_true(hair$test(KASP_LIMITS$hairpin_dg - 3))
  expect_true(homo$test(KASP_LIMITS$homodimer_dg - 3))
  # every bar is labelled as a relaxed one, so the table never hides the trade-off
  for (lab in unlist(KASP_GATE_CHECK)) {
    expect_true(grepl("relaxed", marker_gate_relaxed_bar(lab, KASP_LIMITS)$limit), info = lab)
  }
})

test_that("only a soft gate maps to an acceptance check (C5)", {
  # a hard gate can never be relaxed, so it has no relaxed bar to read
  hard <- unname(unlist(lapply(marker_decision_stages("snp"), function(st) setdiff(st$gates, st$relax))))
  expect_length(intersect(hard, names(KASP_GATE_CHECK)), 0L)
  # ... and every soft gate the ladder lists has its check mapped
  soft <- unname(unlist(lapply(marker_decision_stages("snp"), function(st) st$relax)))
  expect_setequal(setdiff(soft, "gate_wobble"), intersect(soft, names(KASP_GATE_CHECK)))
})

# ---- the AS stage as a ranked oligo (grade_as_candidates) --------------------
# A minimal AS-candidate table with the columns grade_as_candidates() reads.
as_table <- function(discrimination, tm = 58, dg_3p = -1, wobble = FALSE,
                     n_var = 0L, w_var = 0, n_3p = 0L) {
  tibble::tibble(
    discrimination = discrimination, tm = tm, dg_3p = dg_3p,
    wobble = wobble, n_var = as.integer(n_var), w_var = as.numeric(w_var),
    n_3p = as.integer(n_3p)
  )
}

pick_as <- function(d) {
  g <- grade_as_candidates(d, KASP_LIMITS, min_disc = 0.5)
  select_with_ladder(g$data, g$gates, KASP_RELAX_AS, g$tiers, g$ties, g$gate_margin, g$gate_floor)
}

test_that("the AS 3' end is ranked for weakness, not stability (KASP rationale)", {
  # KASP holds because Taq fails to extend the primer on the wrong allele, so a
  # WEAK 3' anchor is wanted; a strong one lets the wrong allele extend anyway.
  expect_true("3' end weakness" %in% marker_decision_stages("snp")$as$tiers)
  # equal on every other tier: the weaker 3' end (dg_3p closer to 0) must win
  expect_equal(pick_as(as_table(c(0.9, 0.9), dg_3p = c(-6, -1)))$best$dg_3p, -1)
})

test_that("the AS tie-break is content-only, so it cannot rank by chromosome position", {
  # A tie on seed_hits + as_length used to fall through to as_start, an ABSOLUTE
  # coordinate, so the last word on which AS primer wins was where it happened to
  # sit.  The same locus handed over at a different position then flipped the
  # winner.  Only content keys may decide.
  expect_identical(KASP_TIES_AS, c("seed_hits", "as_length"))
  expect_identical(KASP_TIES_PAIR, c("run", "ssr"))

  # Two placements tying on every ranked margin, listed 5'-most first: same
  # length, opposite strands, seed_hits UNKNOWN (the seed screen never reached
  # them, which is where all 25 observed as_start ties lived).  Without as_start
  # the earlier row wins at every position, i.e. from content and input order; a
  # chromosome coordinate would have overruled both.
  placements <- function(pos) {
    tibble::tibble(
      discrimination = c(0.9, 0.9), tm = 58, dg_3p = -1, wobble = FALSE,
      n_var = 0L, w_var = 0, n_3p = 0L,
      seed_hits = c(NA_integer_, NA_integer_), as_length = c(20L, 20L),
      strand = c("-", "+"), as_id = c("leftmost", "rightmost"),
      # `-` starts at pos, `+` ends at pos, so row 1 is the 5'-most placement
      as_start = ifelse(c("-", "+") == "+", pos - 20L + 1L, pos)
    )
  }
  winners <- vapply(c(100L, 200L, 300L), function(pos) pick_as(placements(pos))$best$as_id, character(1))
  expect_identical(unname(winners), rep("leftmost", 3L))
  # and the winner's identity is independent of where the locus sits
  expect_length(unique(winners), 1L)
})

test_that("the AS seed screen is partial: an unscreened row is NA and cannot outrank a screened one", {
  # Only the top `specificity_candidates` placements are seed-screened, so
  # `seed_hits` is a PARTIAL tie-break: NA means "not screened", not "not
  # specific".  Two rows tying on every margin differ only in that column; the
  # (lower = more specific) must win, i.e. NA sorts last rather than dropping the
  # row or dropping the column.
  placements <- function(hits) {
    tibble::tibble(
      discrimination = c(0.9, 0.9), tm = 58, dg_3p = -1, wobble = FALSE,
      n_var = 0L, w_var = 0, n_3p = 0L,
      seed_hits = hits, as_length = c(20L, 20L), as_id = c("row1", "row2")
    )
  }
  expect_identical(pick_as(placements(c(50L, NA_integer_)))$best$as_id, "row1")
  expect_identical(pick_as(placements(c(NA_integer_, 50L)))$best$as_id, "row2")
  # 50 hits is the more specific row here, so it also beats a screened 80
  expect_identical(pick_as(placements(c(80L, 50L)))$best$as_id, "row2")
})

test_that("marker_pick_as screens exactly `specificity_candidates` placements", {
  skip_if_no_oligoarrayaux()
  # a random genome so the seeds are unique; three placements that tie on every
  # margin, so the reported NA count is the screen width, not the ranking
  withr::local_seed(11)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  genome <- Biostrings::DNAStringSet(c(chr1 = rand(400L)))
  as_cand <- tibble::tibble(
    ref_primer = vapply(1:3, function(i) rand(24L), character(1)),
    discrimination = 0.9, tm = 58, dg_3p = -1,
    wobble = FALSE, n_var = 0L, w_var = 0, n_3p = 0L
  )
  pick <- function(k) {
    marker_pick_as(as_cand, genome,
      min_disc = 0.5, specificity_candidates = k, BPPARAM = NULL,
      snp = list(chrom = "chr1", pos = 200L), tmpl = NULL
    )
  }
  short <- pick(2L)$as_cand
  expect_equal(sum(!is.na(short$seed_hits)), 2L) # two screened of three rows
  expect_equal(sum(is.na(pick(0L)$as_cand$seed_hits)), 3L) # 0 = no screen at all
})

test_that("the COM primer is also ranked for a weaker 3' end (README section 3)", {
  expect_true("COM 3' end weakness" %in% marker_decision_stages("snp")$pair$tiers)
  # the two pairings differ only in the COM primer's 3'-terminal dG37
  d <- pair_table(c(0.9, 0.9))
  d$com_dg_3p <- c(-8, -1)
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_true("tier5" %in% names(g$data))
  pick <- select_with_ladder(g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties, g$gate_margin, g$gate_floor)
  expect_equal(pick$best$com_dg_3p, -1) # the weaker COM 3' end wins
})

test_that("a footprint trade-off is spent before the wobble is touched", {
  # only the footprint blocks: the ladder must lower it and leave the wobble
  # (the graver sin) alone
  d <- as_table(c(0.9, 0.9), n_var = c(3L, 1L), w_var = c(0.3, 0.1))
  pick <- pick_as(d)
  expect_true("gate_as_footprint" %in% pick$relaxed)
  expect_false("gate_wobble" %in% pick$relaxed)
  expect_equal(pick$best$n_var, 1L) # still the least-laden placement
})

test_that("gate_wobble is really evaluated and relaxable for the AS stage (C7)", {
  # a wobble at the 3' base is a soft failure: it blocks, then the ladder may
  # lower it (so the gate is not decoration - it can fail and be relaxed)
  d <- as_table(c(0.9, 0.9), wobble = c(TRUE, TRUE))
  g <- grade_as_candidates(d, KASP_LIMITS, min_disc = 0.5)
  expect_false(any(g$data$gate_wobble))
  pick <- pick_as(d)
  expect_true("gate_wobble" %in% pick$relaxed)
  expect_true(is.logical(d$wobble))
})

test_that("a relaxed structure floor still passes the structure checks (C5)", {
  # the amplicon/structure floors widen the bar; the checks read the same bar, so
  # "relaxed" and "accepted" cannot contradict
  expect_equal(-KASP_GATE_FLOOR_PAIR[["gate_amplicon"]], 20L)
  expect_equal(-KASP_GATE_FLOOR_PAIR[["gate_structure"]], 3L)
  for (lab in c("worst hairpin dG37", "worst homodimer dG37")) {
    bar <- marker_gate_relaxed_bar(lab, KASP_LIMITS)
    base <- if (lab == "worst hairpin dG37") KASP_LIMITS$hairpin_dg else KASP_LIMITS$homodimer_dg
    expect_true(bar$test(base - 3))
    expect_false(bar$test(base - 4))
  }
})

# ---- the snp DAG, explicitly (gate set, order, hardness) --------------------
test_that("the snp DAG obeys the contract and declares its gates and hardness", {
  expect_dag_invariants("snp")
  # gate set in application order (C6/C7): the AS stage then the pair stage;
  # hardness (C8): only the AS discrimination and the AS 3' window are hard;
  # and one DOOMED terminal (C4)
  expect_kasp_core_dag("snp")
})

test_that("the marker-agnostic primer scheme is stage-less by design (C13)", {
  # `primer` is not a genotyping chemistry: it carries no allele, so it has no
  # selection stage and no locus-level designer (documentation-only).  Its
  # marker-agnostic checks are still evaluable on any assembled pair.
  expect_length(marker_decision_stages("primer"), 0L)
  expect_null(marker_design_at(NULL, NULL, "primer"))
  s <- marker_decision_schema("primer")
  expect_false(any(s$nodes$kind %in% c("gate", "relax")))
  # the one DOOMED is the shared assay-verdict terminal: any assembled pair can
  # fail an acceptance check, so it does not mean a selection stage came back
  expect_equal(s$nodes$name[s$nodes$kind == "doomed"], "assay:doomed")
  expect_setequal(s$nodes$label[s$nodes$kind == "check"], KASP_CHECK_PRIMER)
  expect_dag_invariants("primer")
})
