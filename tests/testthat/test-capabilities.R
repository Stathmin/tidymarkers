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
