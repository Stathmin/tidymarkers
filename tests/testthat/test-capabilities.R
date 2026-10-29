# The engine capability seam (ZA2.1): design code asks `engine_has()` instead of
# probing the tool, so a run context's own adapter decides what is available.

test_that("engine_has() reads the live adapter when the context carries none", {
  ctx <- marker_run_context()
  expect_identical(engine_has("duplex_dg", ctx), has_oligoarrayaux())
  expect_identical(engine_has("duplex_dg", NULL), has_oligoarrayaux())
})

test_that("engine_has() reports an undeclared capability as FALSE, never TRUE", {
  ctx <- marker_run_context(engines = marker_engines())
  expect_false(engine_has("not_a_capability", ctx))
  # a capability the live table genuinely declares still answers TRUE
  expect_true(engine_has("tm", ctx))
})

test_that("engine_has() honours the adapter supplied in the run context", {
  off <- marker_engines()
  off$thermo$capabilities$duplex_dg <- FALSE
  ctx <- marker_run_context(engines = off)
  expect_false(engine_has("duplex_dg", ctx))
  expect_true(engine_has("tm", ctx))
})

test_that("a context without duplex_dg makes the mismatch dG NA, not a silent pass", {
  f <- snp_fixture()
  off <- marker_engines()
  off$thermo$capabilities$duplex_dg <- FALSE
  ctx <- marker_run_context(genome = f$genome, engines = off)
  m <- suppressMessages(marker_design_at(f$genome, f$snp, "taqman", context = ctx))
  # the mismatch dG is a hard gate: with the capability off it is NA, so the
  # candidate is doomed rather than passing a gate nobody could evaluate
  expect_null(m)

  # and the dG routine itself answers NA rather than calling the tool
  expect_true(is.na(taqman_mismatch_dg("ACGTACGTACGTACGTACGT", 10L, "A", engines = off)))
})

test_that("the default context leaves every thermo verdict unchanged", {
  f <- snp_fixture()
  on <- design_taqman_marker(f$genome, f$snp, flank = 150L)
  skip_if(is.null(on), "no TaqMan marker designs at this fixture")
  via_ctx <- marker_design_at(
    f$genome, f$snp, "taqman",
    context = marker_run_context(genome = f$genome, engines = marker_engines())
  )
  expect_identical(via_ctx$probes, on$probes)
  expect_identical(via_ctx$selection, on$selection)
})

test_that("a machine without OligoArrayAux never claims a stable 3' end", {
  fip <- list(
    primer = "AGACCGCTAGAGTGCGCCAACGACAATATCCGCGAGTGAGATTC",
    arm5 = c(238L, 257L), arm3 = c(277L, 300L)
  )
  bip <- list(
    primer = "AGGTGCGGACATAGCTGTAGGATGCTCAACACAAATTCCCGCTGC",
    arm5 = c(324L, 346L), arm3 = c(363L, 384L)
  )
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  g <- lamp_geometry(fip, bip, max_product = 300L)
  # the clamp/fold checks need no engine and still pass
  expect_true(g$f2_b2_ok)
  expect_true(g$loop_ok)
  # the dG half is unknown, not a pass, and is named as unknown
  expect_true(is.na(g$fip_3p_stable))
  expect_true(is.na(g$bip_3p_stable))
  expect_setequal(g$unknown, c("fip_3p_stability", "bip_3p_stability"))
  expect_length(g$failed, 0L)
  expect_true(is.na(g$ok))
})

test_that("a machine without OligoArrayAux cannot grade SNaPshot structure", {
  f <- snp_fixture(n = 800L, seed = 1L, pos = 400L)
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  # the structure margin is NA for every candidate: the design reports the
  # unevaluable gate by name and returns no marker - never a pass
  expect_message(
    m <- design_snapshot_marker(f$genome, f$snp, flank = 250L),
    "gate_snapshot_structure"
  )
  expect_null(m)
})

test_that("generic checks report NA, never TRUE, without the structure engine", {
  f <- snp_fixture()
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  m <- design_sanger_marker(f$genome, f$snp, flank = 250L)
  skip_if(is.null(m), "no Sanger design at this fixture")
  e <- marker_design_checks(m, scheme = "sanger")
  rows <- e[e$criterion %in% c(
    "worst hairpin dG37", "worst cross-dimer dG37", "worst homodimer dG37"
  ), ]
  expect_gt(nrow(rows), 0L)
  # the table's `status` is "unknown" when the value could not be computed
  expect_true(all(is.na(rows$observed)))
  expect_true(all(rows$status == "unknown"))
  expect_false(any(rows$status == "pass"))
})
