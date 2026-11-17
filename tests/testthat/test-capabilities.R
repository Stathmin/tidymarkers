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

test_that("the InDel size designer forwards a per-call engines override to the structure screen", {
  f <- snp_fixture(n = 800L, seed = 1L, pos = 400L)
  indel <- tibble::tibble(chrom = "Chr1", pos = f$snp$pos, ref = f$ref, alt = paste0(f$ref, "TTTTTT"))
  # a string is not an engine registry: reached by the structure gate it aborts
  # there; dropped from the forward it would fall through `...` into
  # in_silico_pcr() and be swallowed by its tryCatch instead
  expect_error(
    suppressMessages(design_indel_size_marker(f$genome, indel, flank = 250L, engines = "not an adapter")),
    "invalid for atomic vectors"
  )
  off <- marker_engines()
  off$thermo$capabilities$duplex_dg <- FALSE
  on <- marker_engines()
  on$thermo$capabilities$duplex_dg <- TRUE
  # both on and off must reach the structure gate as the override says
  expect_no_error(
    invisible(suppressMessages(design_indel_size_marker(f$genome, indel, flank = 250L, engines = on)))
  )
  expect_no_error(
    invisible(suppressMessages(design_indel_size_marker(f$genome, indel, flank = 250L, engines = off)))
  )
})

test_that("marker_primer_checks() forwards a per-call engines override to the thermo screen", {
  local_mocked_bindings(has_oligoarrayaux = function() TRUE, .package = "tidymarkers")
  x <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "G"),
      amplicon = tibble::tibble(start = 250L, end = 350L),
      primers = list(
        forward = "ACGTACGTACGTACGTACGT", reverse = "TGCATGCATGCATGCATGCA"
      )
    ),
    class = "marker_caps"
  )
  # a dG the screen rejects: the tool answers, the adapter decides whether the
  # screen that calls it runs at all
  local_mocked_bindings(
    oligo_dg_vec = function(...) rep(-20, 2L),
    thermo_duplex_dg = function(...) -20,
    .package = "tidymarkers"
  )
  off <- marker_engines()
  off$thermo$capabilities$hairpin_dg <- FALSE
  off$thermo$capabilities$duplex_dg <- FALSE
  # the dG rows are emitted either way; the forwarded override decides their
  # status: "unknown" when the thermo screen did not run, otherwise a verdict
  dg_status <- function(chk) {
    rows <- chk[chk$name %in% c("hairpin", "dimer", "homodimer"), ]
    expect_equal(nrow(rows), 3L)
    rows$status
  }
  expect_true(all(dg_status(marker_primer_checks(x, "caps", engines = off)) == "unknown"))
  expect_true(all(dg_status(marker_primer_checks(x, "caps", engines = marker_engines())) == "broken"))
  expect_true(all(dg_status(marker_primer_checks(x, "caps")) == "broken"))
})

test_that("a composed technology forwards a per-call engines override to its designer", {
  f <- snp_fixture()
  # the generic (non-KASP) evaluator is where the adapter decides the dG rows:
  # it grades by calling the engines, so the override must reach the designer
  aux <- has_oligoarrayaux()
  skip_if_not(aux, "the structure rows need OligoArrayAux")
  local_mocked_bindings(
    oligo_dg_vec = function(...) rep(-20, 2L),
    thermo_duplex_dg = function(...) -20,
    .package = "tidymarkers"
  )
  on <- marker_engines()
  on$thermo$capabilities$hairpin_dg <- TRUE
  on$thermo$capabilities$duplex_dg <- TRUE
  off <- marker_engines()
  off$thermo$capabilities$hairpin_dg <- FALSE
  off$thermo$capabilities$duplex_dg <- FALSE
  marker_technology("caps_probe_test", base = "caps", overwrite = TRUE)
  on.exit(rm("caps_probe_test", envir = .marker_tech_store), add = TRUE)
  # a TaqMan technology over the KASP core: its designer is assemble_kasp_assay
  marker_technology("snp_probe_test", base = "snp", overwrite = TRUE)
  on.exit(rm("snp_probe_test", envir = .marker_tech_store), add = TRUE)
  # a probe the designer cannot accept: the wrapper's explicit forward makes
  # `...` a duplicate argument; without it the probe would just be ignored
  for (s in c("caps_probe_test", "snp_probe_test")) {
    expect_error(
      scheme_designer(s)(f$genome, f$snp, flank = 150L, engines = "not an adapter", engines = NULL),
      "formal argument \"engines\" matched by multiple actual arguments",
      label = s
    )
  }
  # the else branch: a primer-base technology delegates to its base's designer,
  # and the index it returns is one OligoArrayAux takes; the adapter it is handed
  # shows in the checks the design path runs on the result
  d <- scheme_designer("caps_probe_test")
  expect_equal(names(formals(d)), c("genome", "snp", "...", "engines"))
  m <- suppressWarnings(d(f$genome, f$snp, flank = 150L, engines = on))
  skip_if(is.null(m), "no composed CAPS design at this fixture")
  dg_status <- function(eng) {
    rows <- marker_primer_checks(m, "caps", engines = eng)
    rows <- rows[rows$name %in% c("hairpin", "dimer", "homodimer"), ]
    expect_equal(nrow(rows), 3L)
    rows$status
  }
  expect_true(all(dg_status(off) == "unknown"))
  expect_true(all(dg_status(on) == "broken"))
})

test_that("the KASP designer's own evaluation honours a per-call engines override", {
  skip_if_not(has_oligoarrayaux(), "the structure row needs OligoArrayAux")
  f <- snp_fixture()
  off <- marker_engines()
  off$thermo$capabilities$hairpin_dg <- FALSE
  off$thermo$capabilities$duplex_dg <- FALSE
  expect_message(
    a <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp, flank = 150L, engines = off)),
    "gate_structure"
  )
  # an excluded structure engine dooms the pairing, and the doom names the gate
  expect_s3_class(a, "marker_doomed")
  expect_match(a$reason, "gate_structure")
})

test_that("an all-NA structure margin never wins the ladder and names the gate", {
  d <- tibble::tibble(
    discrimination = c(0.9, 0.9),
    com_homodimer = c(NA, -1), com_hairpin = c(NA, -0.5),
    com_tm = 60, amplicon = 80, n_var = 0L, as_n_3p = 0L,
    run = 0L, ssr = 0L, com_start = 1:2, com_dg_3p = -1
  )
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  # NA is unknown: the gate does not reject it
  expect_true(all(g$data$gate_structure))
  pick <- select_with_ladder(g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties, g$gate_margin, g$gate_floor)
  # ...but the ladder never ranks it: the evaluated row wins and no gate was relaxed
  expect_identical(pick$best$.row, 2L)
  expect_false("gate_structure" %in% pick$violated)
  # with no margin at all the pairing is doomed, never ranked
  d$com_homodimer[] <- NA_real_
  d$com_hairpin[] <- NA_real_
  g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
  expect_null(
    select_with_ladder(g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties, g$gate_margin, g$gate_floor)$best
  )
})
