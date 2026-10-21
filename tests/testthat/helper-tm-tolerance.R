# The per-model-pair Tm tolerances, and the oligo set that probes the edges.
#
# This file is the single place both live: a test asks `tm_tolerance()` for a
# pair instead of carrying its own constant, and asks `tm_edge_oligos()` for the
# GC-edge / homopolymer probes instead of restating those sequences.  Nobody
# keeps a second copy of these numbers.
#
# Why a pair and not one constant: `primer_tm()` exposes four published DNA
# nearest-neighbour tables.  SantaLucia (2004) and Allawi (1998) share a
# lineage and agree to well under a degree; Sugimoto (1996) and Breslauer
# (1986) are older lineages whose gap to SantaLucia is large and
# composition-dependent.  A blanket "|dTm| <= 2 C" is false in both directions
# (it is too loose for the same lineage, far too tight across lineages), so the
# tolerance is keyed on the model pair (ROADMAP ZW3).
#
# Measured here (R 4.6.1, TmCalculator 1.1.2, default reaction_conditions():
# 50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTPs, 0 DMSO, 20-nt oligos) unless
# `tm_tolerances()` says otherwise.

# `tm_tolerances()` says otherwise.
#
# `dev/check.R` also loads this file into the package environment (to lint it),
# so keep this file free of helper names that could shadow the package.

NN_SANTA_LUCIA <- "DNA_NN_SantaLucia_2004"
NN_ALLAWI <- "DNA_NN_Allawi_1998"
NN_SUGIMOTO <- "DNA_NN_Sugimoto_1996"
NN_BRESLAUER <- "DNA_NN_Breslauer_1986"

# The four tables the oracle compares.  `primer_tm()`'s default must be one of
# them, which test-tm-oracle.R asserts.
tm_models <- function() {
  c(NN_SANTA_LUCIA, NN_ALLAWI, NN_SUGIMOTO, NN_BRESLAUER)
}

# One row per ordered-insensitive model pair.
#
# `basis` states where the number came from, so a reader can re-measure it:
#   * SantaLucia vs Allawi: measured max |dTm| 0.46 C over the 32-oligo set
#     (including the edge probes); the 0.5 C cap absorbs a TmCalculator
#     salt-code change without admitting a real model difference.
#   * Sugimoto vs SantaLucia: not a bound but a one-sided floor, asserted only
#     on the unremarkable 20-mers of the oracle set (min 2.2 C there).  Across
#     all DNA 20-mers the gap falls to about -1 C at 70-80 % GC (25 of the
#     3501 20-mers sampled are < 0), so it is NOT a global one-sided rule.
#   * Breslauer vs SantaLucia: asserted only on the GC-rich subset, where the
#     measured min is 9.5 C (11 same-composition 19-mers: 7.9 C).
tm_tolerances <- function() {
  tibble::tribble(
    ~model_a, ~model_b, ~max_delta_c, ~basis,
    NN_SANTA_LUCIA, NN_ALLAWI, 0.5,
    paste(
      "same lineage; max |dTm| 0.46 C measured on the 32-oligo oracle set",
      "(TmCalculator 1.1.2, 20-nt oligos, default reaction_conditions); the",
      "0.5 C cap absorbs a TmCalculator salt-code change"
    ),
    NN_SUGIMOTO, NN_SANTA_LUCIA, 2.0,
    paste(
      "one-sided floor on the unremarkable 20-mers only, NOT a bound:",
      "min 2.2 C / mean 4.1 C measured there; across all 20-mers the gap",
      "reaches about -1 C at 70-80 % GC, so it is not global"
    ),
    NN_BRESLAUER, NN_SANTA_LUCIA, 9.5,
    paste(
      "GC-rich subset (GC >= 65 %) only, where the gap is real: min 9.5 C",
      "measured there; on A/T-rich 20-mers it is smaller (5.7 C on a balanced",
      "A/T probe, up to 10.9 C on an A homopolymer), so no A/T-side bound is pinned"
    )
  )
}

# GC content in percent, per oligo.
gc_content <- function(seq) {
  100 * nchar(gsub("[^GC]", "", seq)) / nchar(seq)
}

# Look up one pair's tolerance; order-insensitive, errors on an unknown pair.
tm_tolerance <- function(model_a, model_b) {
  tab <- tm_tolerances()
  hit <- (tab$model_a == model_a & tab$model_b == model_b) |
    (tab$model_a == model_b & tab$model_b == model_a)
  if (!any(hit)) {
    rlang::abort(paste0("no Tm tolerance for ", model_a, " vs ", model_b))
  }
  tab[hit, , drop = FALSE]
}

# The GC-edge / homopolymer probes: the sequences the "explicit tolerance" rule
# is easiest to get wrong on.  Named, so a failure says which edge broke.
# `tm_oligos()` (in test-tm-oracle.R) is the main set; this is the subset that
# exercises the edges, and it is defined here once.
tm_edge_oligos <- function() {
  c(
    # 6/20 G/C: low but not degenerate, so it is the unremarkable A/T edge
    at_rich = "AATATTACCGTCGATACTAT", # GC = 30 % (< 35 %)
    gc_rich = "GCGCGCGCGGGCGCGGCGCG", # GC = 100 % (> 65 %)
    homo_a = "GATGATAAAAAAAAATGATC", # 9-nt A run, mid-oligo
    homo_g = "CATCATGGGGGGGGATCATG", # 8-nt G run, mid-oligo
    three_prime_a = "ATTAAACAACACACAAAAAA", # 6-nt A run AT the 3' end
    three_prime_g = "GGGGGGGGATGCATGCATCA", # 8-nt G run opening at the 5' end
    homo5 = "CAACAGCAAAAACAGACAAC", # exactly a 5-nt run, for the 5+ boundary
    poly_a = "AAAAAAAAAAAAAAAAAAAA" # degenerate: 0 % GC, a 20-mer homopolymer
  )
}
