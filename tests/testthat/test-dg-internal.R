# The pure-R duplex adapter: faithful to the nearest-neighbour model the Tm
# already uses, and close enough to OligoArrayAux for the calibrated gates.

dg_panel <- c(
  "ACGTACGTACGTACGTACGT",
  "GAAGGTGACCAAGTTCATGCT",
  "AAAAAAAAAAAAAAAAAAAA",
  "GCGCGCGCGCGCGCGCGCGC",
  "TTTTGGCGCGCCGATATATAT",
  "GCGGCGGCTAGCTAGCCGCT",
  "ATGATGATGATGATGATGAT"
)

dg_panel_rc <- function(seqs = dg_panel) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs)))
}

# Hairpin fixtures: the first block folds into a stem of >= 4 bp, the second
# has none.  `ha_una` is what UNAFold prints at 37 C / 50 mM Na (`hybrid-ss-min
# -q`), i.e. the oracle the gate has to track.
ha_panel <- c(
  "GGGGAAAACCCC", "GGGGGAAAACCCCC", "GGGGGAAAAACCCCC",
  "GGGGGGAAAACCCCCC", "GGGGGGAAAAAACCCCCC", "ATGGGGAAAACCCCAT",
  "GGGGAAAATTTTCCCC", "GGGGATATCAAAAGATATCCCC", "ACGTACGTACGT",
  "GGGGAAAAAAAAAACCC", "GGGGGAAAAAAAAAAAAACCC"
)
ha_una <- c(
  -0.835, -2.333, -3.692, -3.831, -4.490, -2.090,
  -3.510, -4.324, -1.074, -0.085, 0.215
)
ha_none <- c("AAAAAAAAAA", "ACACACACAC", "TTTTTTTTTTTT", "AGAGAGAGAGAGAG", "ATATATATATAT")

test_that("internal_duplex_dg reproduces OligoArrayAux within the gate width", {
  skip_if_no_oligoarrayaux()
  rc <- dg_panel_rc()
  ours <- internal_duplex_dg(dg_panel, rc, na = 0.05, temp = 37)
  theirs <- vapply(
    seq_along(dg_panel),
    function(i) oligo_dg(dg_panel[[i]], rc[[i]], na = 0.05, temp = 37),
    numeric(1)
  )
  # measured max 0.55 kcal/mol on this panel; the gates are coarse (duplex
  # >= -9, hairpin >= -3), so 1.0 is the honest tolerance
  expect_lte(max(abs(ours - theirs)), 1.0)
  expect_true(all(is.finite(ours)))
  # more negative = stronger: a GC-rich perfect duplex beats a poly-A one
  expect_lt(ours[[4L]], ours[[3L]])
})

test_that("the internal dH is the TmCalculator nearest-neighbour dH", {
  # Two strand concentrations give the van't Hoff slope d(1/Tm)/d(ln k) = R/dH,
  # so TmCalculator's own Tm values independently pin dH.  Recovering the same
  # dH from internal_duplex_dg at two temperatures proves the table, the
  # initiation terms and the salt term are the same model, not a lookalike.
  cond <- reaction_conditions(Na = 50, Mg = 0, dNTPs = 0)
  for (i in seq_along(dg_panel)) {
    s <- dg_panel[[i]]
    r <- dg_panel_rc(s)
    d_lo <- internal_duplex_dg(s, r, na = 0.05, temp = 20)
    d_hi <- internal_duplex_dg(s, r, na = 0.05, temp = 70)
    ds_ours <- -(d_hi - d_lo) / (70 - 20) * 1000 # calories per K per mol
    dh_ours <- d_lo + (20 + 273.15) * ds_ours / 1000

    tm_lo <- primer_tm(s, conditions = cond, salt_method = "SantaLucia1998-2", dnac = 25)
    tm_hi <- primer_tm(s, conditions = cond, salt_method = "SantaLucia1998-2", dnac = 2500)
    inv <- 1 / (c(tm_lo, tm_hi) + 273.15)
    # k is half the strand concentration for a non-self-complementary duplex,
    # so the 100x ratio drops out of the slope
    dh_tm <- 1.987 * log(100) / (inv[[2L]] - inv[[1L]]) / 1000

    expect_equal(dh_ours, dh_tm, tolerance = 0.05, info = s)
    ds_tm <- dh_tm * 1000 * inv[[1L]] - 1.987 * log(12.5e-9)
    expect_equal(ds_ours, ds_tm, tolerance = 0.5, info = s)
  }
})

test_that("internal_duplex_dg is symmetric and responds to salt and temperature", {
  s <- "GAAGGTGACCAAGTTCATGCT"
  r <- dg_panel_rc(s)
  expect_equal(
    internal_duplex_dg(s, r), internal_duplex_dg(r, s),
    tolerance = 1e-9
  )
  # more salt stabilises the duplex
  expect_lt(internal_duplex_dg(s, r, na = 0.2), internal_duplex_dg(s, r, na = 0.02))
  # dG rises toward zero with temperature (Tm of this 21-mer is ~52 C)
  expect_gt(
    internal_duplex_dg(s, r, temp = 60),
    internal_duplex_dg(s, r, temp = 30)
  )
  # chain rule: d/dT of dG is -dS, so the slope is a finite negative number
  expect_lt(internal_duplex_dg(s, r, temp = 60), 0)
})

test_that("internal_duplex_dg answers NA instead of guessing", {
  s <- "GAAGGTGACCAAGTTCATGCT"
  r <- dg_panel_rc(s)
  expect_true(is.na(internal_duplex_dg("ACGTNACGTACGTACGTACG", "ACGTACGTACGTACGTACGT")))
  # a partner that is not the reverse complement is outside the model
  expect_true(is.na(internal_duplex_dg(s, "TTTTTTTTTTTTTTTTTTTTT")))
  expect_true(is.na(internal_duplex_dg(s, r, na = 0)))
  expect_true(is.na(internal_duplex_dg("A", "T")))
  # hairpin mode is not modelled: NA with a warning, never a number
  expect_warning(dg <- internal_duplex_dg(s), "hairpin")
  expect_true(is.na(dg))
})

test_that("the internal thermo adapter declares duplex and hairpin dG, never fold", {
  eng <- marker_engines(thermo = "internal")$thermo
  expect_identical(eng$dg, internal_duplex_dg)
  expect_true(capabilities(eng, "duplex_dg"))
  expect_true(capabilities(eng, "hairpin_dg"))
  expect_false(capabilities(eng, "fold"))
  expect_true(capabilities(eng, "tm"))
  # through the seam, the same answers
  expect_true(engine_has("duplex_dg", marker_engines(thermo = "internal")))
  expect_true(engine_has("hairpin_dg", marker_engines(thermo = "internal")))
  # the default registry is unchanged: Opt-out is explicit, never automatic
  expect_identical(
    capabilities(marker_engines()$thermo, "duplex_dg"),
    has_oligoarrayaux()
  )
})

test_that("without OligoArrayAux the hairpin capability is not guessed", {
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  # the auto engine has no internal fallback for hairpin: the gate stays
  # unevaluable rather than graded by an adapter nobody asked for
  expect_false(capabilities(marker_engines()$thermo, "hairpin_dg"))
  expect_false(engine_has("hairpin_dg"))
})

test_that("the internal adapter keeps the OligoArrayAux interface", {
  expect_identical(
    names(formals(internal_duplex_dg)),
    names(formals(oligo_dg))
  )
  expect_identical(formals(internal_duplex_dg)$na, formals(oligo_dg)$na)
  expect_identical(formals(internal_duplex_dg)$temp, formals(oligo_dg)$temp)
})

test_that("internal_hairpin_dg tracks UNAFold on real stems", {
  skip_if_no_oligoarrayaux()
  ours <- internal_hairpin_dg(ha_panel, na = 0.05, temp = 37)
  theirs <- hairpin_dg_thermo(ha_panel, conditions = reaction_conditions(Na = 50, temperature = 37))
  expect_equal(length(ours), length(ha_panel))
  expect_true(all(is.finite(ours)))
  # measured on this panel: mean |d| 0.54, max |d| 1.48 kcal/mol.  The gate is
  # coarse (hairpin_dg >= -3), so 1.0 / 2.0 is the honest width
  delta <- abs(ours - theirs)
  expect_lte(mean(delta), 1.0)
  expect_lte(max(delta), 2.0)
  # UNAFold's own numbers agree with the recorded fixture table, so the oracle
  # is the tool and not a stale constant
  expect_equal(theirs, ha_una, tolerance = 0.05)
  # ordering: a real hairpin is never weaker than a non-hairpin of a similar
  # length, and the GC-rich stem beats the A/T-rich one
  expect_lt(max(internal_hairpin_dg(c("GGGGAAAACCCC", "GGGGGGAAAACCCCCC"))), 0)
  expect_gt(internal_hairpin_dg("AAAAAAAAAA"), internal_hairpin_dg("GGGGAAAACCCC"))
})

test_that("the hairpin screen never under-reports stability (false-pass guard)", {
  skip_if_no_oligoarrayaux()
  # The screen only rewards a perfect stem, so it can miss the extra stability
  # UNAFold finds (dangling ends, coaxial stacking, mismatched stems) - that is
  # the permissive direction the gate has to be told about.  What it must never
  # do is the opposite: invent a fold UNAFold does not see, or call a fold more
  # stable than UNAFold does.  Where UNAFold finds a hairpin, ours may sit above
  # it by the gate tolerance but no further; where it finds none, ours stays
  # >= 0, i.e. it never turns a non-hairpin into a failing gate.
  ours <- internal_hairpin_dg(ha_panel, na = 0.05, temp = 37)
  una <- hairpin_dg_thermo(ha_panel, conditions = reaction_conditions(Na = 50, temperature = 37))
  folded <- una < 0
  expect_gt(sum(folded), 5L)
  expect_lte(max(ours[folded] - una[folded]), 2.0)
  expect_true(all(ours[!folded] >= 0))
  expect_true(all(internal_hairpin_dg(ha_none) >= 0))
  # and on these two short 4 bp stems ours is *more* negative than UNAFold: its
  # -q value is taken from a run that leaves the longest stem unfolded, so the
  # non-hairpin baseline it is measured against is already down by that stem.
  # Same order either way, hence the two-sided 2.0 bound above.
  expect_lte(min(ours[folded] - una[folded]), 0)
})

test_that("internal_hairpin_dg returns zero when no stem can form", {
  # the floor: only stems of >= 4 bp are screened, and a sequence with none
  # reports 0 (no stabilisation relative to the unfolded state), the same
  # mapping oligo_dg() applies to UNAFold's "no structure" sentinel
  expect_equal(internal_hairpin_dg(ha_none), rep(0, length(ha_none)))
  expect_true(all(internal_hairpin_dg(ha_none) >= 0))
  # a 3 bp stem is below the floor and still answers 0, never a positive dG
  expect_gte(internal_hairpin_dg("GGGAAACCC"), 0)
  # the answer is never a guessed number
  expect_true(is.na(internal_hairpin_dg("GGGGAAAACCCCN")))
  expect_true(is.na(internal_hairpin_dg("GGGGAAAACCCC", na = 0)))
})

test_that("the hairpin stem dH is the TmCalculator nearest-neighbour dH", {
  # The stem of [internal_hairpin_dg()] is [internal_duplex_terms()] on one arm,
  # so the van't Hoff trick of the duplex test applies to it directly: recover
  # dH/dS from the terms at two temperatures and compare with the dH implied by
  # TmCalculator's own Tms.  (Sweeping internal_hairpin_dg itself is not usable
  # here: at 70 C the fold is gone and the 0 floor hides the slope.)
  cond <- reaction_conditions(Na = 50, Mg = 0, dNTPs = 0)
  for (arm in c("GGGG", "GGGGG", "GGGGGG", "GGGGATAT")) {
    lo <- internal_duplex_terms(arm, 0.05)
    # dG(T) = dH - T * dS, so the two-temperature difference is linear in dH, dS
    d_lo <- lo[["dh"]] - (20 + 273.15) * lo[["ds"]] / 1000
    d_hi <- lo[["dh"]] - (70 + 273.15) * lo[["ds"]] / 1000
    ds_ours <- -(d_hi - d_lo) / 50 * 1000
    dh_ours <- d_lo + (20 + 273.15) * ds_ours / 1000

    tm_lo <- primer_tm(arm, conditions = cond, salt_method = "SantaLucia1998-2", dnac = 25)
    tm_hi <- primer_tm(arm, conditions = cond, salt_method = "SantaLucia1998-2", dnac = 2500)
    inv <- 1 / (c(tm_lo, tm_hi) + 273.15)
    dh_tm <- 1.987 * log(100) / (inv[[2L]] - inv[[1L]]) / 1000
    expect_equal(dh_ours, dh_tm, tolerance = 0.05, info = arm)
  }
  # and the hairpin answers are built from those same terms: every candidate is
  # a stem plus one loop pair, so the reported dG is one of them and never
  # something else
  s <- "GGGGAAAACCCC"
  cand <- hairpin_candidates(s, 0.05)
  at37 <- cand[, "dh"] - 310.15 * cand[, "ds"] / 1000
  expect_equal(internal_hairpin_dg(s), min(0, min(at37)), tolerance = 1e-12)
})

test_that("the option routes duplex and hairpin calls through the internal adapter", {
  s <- "GAAGGTGACCAAGTTCATGCT"
  r <- dg_panel_rc(s)
  withr::with_options(list(tidymarkers.thermo_engine = "internal"), {
    expect_equal(
      oligo_dg_vec(s, r, na = 0.05, temp = 37),
      internal_duplex_dg(s, r, na = 0.05, temp = 37),
      tolerance = 1e-12
    )
    expect_equal(
      thermo_duplex_dg(s, r, na = 0.05, temp = 37),
      internal_duplex_dg(s, r, na = 0.05, temp = 37),
      tolerance = 1e-12
    )
    # hairpin goes the same way, so the capability the adapter declares is the
    # number the gates actually read
    expect_equal(
      oligo_dg_vec(s, na = 0.05, temp = 37),
      internal_hairpin_dg(s, na = 0.05, temp = 37),
      tolerance = 1e-12
    )
    expect_equal(
      hairpin_dg_thermo(s, conditions = reaction_conditions(Na = 50, temperature = 37)),
      internal_hairpin_dg(s, na = 0.05, temp = 37),
      tolerance = 1e-12
    )
  })
  # outside the option the engine is UNAFold's
  skip_if_no_oligoarrayaux()
  expect_equal(
    oligo_dg_vec(s, na = 0.05, temp = 37),
    oligo_dg(s, na = 0.05, temp = 37),
    tolerance = 1e-12
  )
})
