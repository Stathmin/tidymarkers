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

test_that("a per-call engines override decides the gate, not the session registry", {
  off <- marker_engines()
  off$thermo$capabilities$duplex_dg <- FALSE
  # the session registry claims the capability, so only the supplied adapter
  # can make these gates answer FALSE
  local_mocked_bindings(has_oligoarrayaux = function() TRUE, .package = "tidymarkers")
  on <- marker_engines()
  on$thermo$capabilities$duplex_dg <- TRUE
  # Sanger pair screen: with the capability off the pair is not graded at all
  # (honest TRUE, the readers are never reached); with it on, a failing dG fails -
  # both dG routines answer what OligoArrayAux would downgrade below -9
  local_mocked_bindings(
    oligo_dg_vec = function(...) rep(-20, 2L),
    thermo_duplex_dg = function(...) -20,
    .package = "tidymarkers"
  )
  pair <- c("GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC")
  expect_true(pair_structure_ok(pair[1L], pair[2L], engines = off))
  expect_false(pair_structure_ok(pair[1L], pair[2L], engines = on))
  expect_identical(pair_structure_ok(pair[1L], pair[2L]), pair_structure_ok(pair[1L], pair[2L], engines = on))

  f <- snp_fixture(n = 800L, seed = 1L, pos = 400L)
  # same override through the public designer: the structure gate becomes
  # unevaluable and is reported by name, never passed
  expect_message(
    m <- suppressWarnings(design_snapshot_marker(f$genome, f$snp, flank = 250L, engines = off)),
    "gate_snapshot_structure"
  )
  expect_null(m)
})

test_that("generic checks report NA, never TRUE, without the structure engine", {
  f <- snp_fixture()
  local_mocked_bindings(has_oligoarrayaux = function() FALSE, .package = "tidymarkers")
  expect_warning(
    m <- design_sanger_marker(f$genome, f$snp, flank = 250L),
    "allele dropout is not checked"
  )
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
