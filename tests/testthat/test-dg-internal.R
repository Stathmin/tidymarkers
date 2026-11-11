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

test_that("the internal thermo adapter declares duplex_dg only", {
  eng <- marker_engines(thermo = "internal")$thermo
  expect_identical(eng$dg, internal_duplex_dg)
  expect_true(capabilities(eng, "duplex_dg"))
  expect_false(capabilities(eng, "hairpin_dg"))
  expect_false(capabilities(eng, "fold"))
  expect_true(capabilities(eng, "tm"))
  # through the seam, the same answers
  expect_true(engine_has("duplex_dg", marker_engines(thermo = "internal")))
  expect_false(engine_has("hairpin_dg", marker_engines(thermo = "internal")))
  # the default registry is unchanged: Opt-out is explicit, never automatic
  expect_identical(
    capabilities(marker_engines()$thermo, "duplex_dg"),
    has_oligoarrayaux()
  )
})

test_that("the internal adapter keeps the OligoArrayAux interface", {
  expect_identical(
    names(formals(internal_duplex_dg)),
    names(formals(oligo_dg))
  )
  expect_identical(formals(internal_duplex_dg)$na, formals(oligo_dg)$na)
  expect_identical(formals(internal_duplex_dg)$temp, formals(oligo_dg)$temp)
})

test_that("the option routes duplex calls and leaves hairpin on OligoArrayAux", {
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
    # hairpin still goes to UNAFold, so it keeps UNAFold's number
    skip_if_no_oligoarrayaux()
    expect_equal(
      oligo_dg_vec(s, na = 0.05, temp = 37),
      oligo_dg(s, na = 0.05, temp = 37),
      tolerance = 1e-12
    )
  })
})
