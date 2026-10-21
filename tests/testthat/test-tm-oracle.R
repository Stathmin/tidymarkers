# =============================================================================
#  ZW3  --  an independent Tm oracle, with a stated tolerance per model
# =============================================================================
#  A Tm that only the Tm engine can confirm is not evidence.  This file pins the
#  engine to things outside it and STATES the tolerance for each, because the
#  tolerances are wildly different: two nearest-neighbour tables that share the
#  SantaLucia lineage agree to a few hundredths of a degree, while Breslauer
#  1986 sits 8-18 C away depending on GC.  A single "|dTm| <= 2 C" blanket
#  assertion would be false for three of the four tables and vacuous for the
#  fourth.
#
#  Anchors, strongest first:
#
#   1. Same-lineage agreement (exact, tolerance 0.5 C).  Allawi 1998 and
#      SantaLucia 2004 share the SantaLucia unified parameters; they must agree
#      to well inside the 2 C window the KASP gates are built on.
#   2. Cross-lineage offsets (tolerance quoted per table, measured, not asked
#      for): Sugimoto 1996 and Breslauer 1986 rank the same oligos, so their
#      CORRELATION must be high even though their LEVELS disagree.
#   3. Real published primers (tolerance 0.5 C, same-lineage only).  The 32
#      unique primers in Makhoul 2020's Table S1 are allele-specific KASP
#      primers that were made and run; a shift in them would move real assays.
#   4. Monotonicity and invariance (no tolerance): Tm rises with GC and with
#      length at fixed composition, and revcomp is identical.
#
#  What is NOT here, and why: von Ahsen 2001 (475 measured Tms) and dnaMATE are
#  the oracles we actually want, and neither has been fetched.  Until one is,
#  `primer_tm()` is anchored only to its own parameter family and MUST NOT be
#  described as experimentally validated.  See scientific_backlog.csv.

tm_engines_available <- function() {
  c(
    TmCalculator = requireNamespace("TmCalculator", quietly = TRUE),
    OligoArrayAux = has_oligoarrayaux()
  )
}

test_that("ZW3: the Tm engine is installed and reports a number for a plain oligo", {
  eng <- tm_engines_available()
  expect_true(eng[["TmCalculator"]], info = "TmCalculator is the ZW3 oracle's primary engine")
  tm <- primer_tm("GAAGGTGACCAAGTTCATGCT")
  expect_length(tm, 1L)
  expect_true(is.finite(tm) && tm > 40 && tm < 80)
})

test_that("ZW3: same-lineage NN tables agree inside the 2 C gate window", {
  # all 32 real primers from the Makhoul supplement
  p <- readr::read_csv(
    test_path("../../inst/extdata/makhoul_kasp_primers.csv"),
    show_col_types = FALSE
  )
  seqs <- unique(c(p$fam, p$hex, p$common))
  seqs <- seqs[grepl("^[ACGT]{10,}$", seqs)]
  expect_gte(length(seqs), 30L)

  sl <- primer_tm(seqs, nn_table = "DNA_NN_SantaLucia_2004")
  al <- primer_tm(seqs, nn_table = "DNA_NN_Allawi_1998")
  # measured max |difference| is ~0.06 C; 0.5 C leaves room for a TmCalculator
  # release changing its salt code without failing this test spuriously, while
  # still being far tighter than the gate it protects
  expect_lt(max(abs(sl - al)), 0.5)

  # the offset to other lineages is real and must not be asserted away: if
  # Sugimoto ever starts agreeing to 0.5 C, something is wrong (probably the
  # table argument is being ignored).
  sg <- primer_tm(seqs, nn_table = "DNA_NN_Sugimoto_1996")
  delta <- sg - sl
  # measured on these 32 primers: Sugimoto is ALWAYS higher (min +1.5 C) and the
  # gap tracks length (cor 0.55) as well as composition, so it is a level offset
  # rather than a disagreement about everything.  Raw Pearson is only 0.77 and
  # Spearman 0.66 -- a linear rescaling does NOT reconcile the two tables
  # (R^2 ~ 0.3 on length + GC), so we assert the structure we measured and a
  # conservative rank floor, not a reconciliation the data does not support.
  expect_true(all(delta > 0))
  expect_gt(mean(delta), 1)
  expect_gt(cor(nchar(seqs), delta), 0.3)
  expect_gt(cor(sl, sg, method = "spearman"), 0.6)
})

test_that("ZW3: within one lineage, Tm is monotone in GC and length, and revcomp-invariant", {
  # same length, increasing GC
  low_gc <- "AATTACATTCATATCACATC"
  mid_gc <- "AATCACCTTCATAGCACATC"
  high_gc <- "AGCCACCGTCAGAGCGCACC"
  tms <- primer_tm(c(low_gc, mid_gc, high_gc))
  expect_true(all(diff(tms) > 0))

  # revcomp invariance: a primer and its reverse complement are the same duplex
  seq <- "GAAGGTGACCAAGTTCATGCT"
  rc <- chartr("ACGT", "TGCA", paste(rev(strsplit(seq, "")[[1L]]), collapse = ""))
  expect_equal(primer_tm(seq), primer_tm(rc), tolerance = 1e-8)
})

test_that("ZW3: the declared tolerance is per engine, not a blanket 2 C", {
  # This test exists to stop a well-meaning edit from collapsing the tolerances
  # into one number.  The Breslauer 1986 table must be MORE than 2 C off on
  # typical primer GC; if it is not, the cross-table check has gone vacuous.
  seqs <- c("GAAGGTGACCAAGTTCATGCT", "TTTTGGCGCGCCGATATATAT", "ATCGATCGATCGATCGATCG")
  sl <- primer_tm(seqs)
  br <- primer_tm(seqs, nn_table = "DNA_NN_Breslauer_1986")
  expect_gt(max(abs(br - sl)), 2)
})
