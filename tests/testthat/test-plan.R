# The layered plan grammar: build with +, solve, tubes, plot, export.
# Solved against the bundled, real Arabidopsis example panel (marker_example()).

# The shipped panel needs the TAIR10 genome; use it when it is available
# (cached, so offline), otherwise skip - a synthetic genome would not exercise
# the real targets.
# The TAIR10 genome is ~35 MB compressed: read it once per file, not once per plan.
example_genome <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      value <<- tryCatch(ath_genome(), error = function(e) NULL)
    }
    value
  }
})

example_plan <- function(ids = c("ENSVATH01332998", "ENSVATH02162468"), preset = "kasp") {
  ex <- marker_example()
  ex <- ex[ex$id %in% ids, , drop = FALSE]
  genome <- example_genome()
  if (is.null(genome) || !nrow(ex)) {
    return(NULL)
  }
  marker_plan(genome) + marker_targets(data = ex) + marker_tech(preset = preset)
}

# Solving is deterministic, and every report/plot test needs the same solved
# panel.  Solve it once so the tests assert the same design without repeating
# the search (and keep its cache for the cache-reuse test).
solved_example <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      plan <- example_plan()
      if (is.null(plan)) {
        return(NULL)
      }
      value <<- suppressWarnings(
        marker_solve(plan, flank = 150L, cores = 1L, progress = FALSE)
      )
    }
    value
  }
})

# a tiny stand-in genome for the pure-logic layers (nothing is designed)
dummy_plan <- function() {
  marker_plan(Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGT")))
}

test_that("the real example panel targets carry id labels and their alt2", {
  ex <- marker_example()
  p <- dummy_plan() + marker_targets(data = ex)
  expect_equal(nrow(p$targets), nrow(ex))
  expect_identical(p$targets$label, ex$id)
  expect_true(all(c("label", "chrom", "pos", "ref", "alt", "alt2", "type") %in% names(p$targets)))
  # the multi-allelic target keeps its second ALT
  multi <- ex$id[ex$type == "SNP" & !is.na(ex$alt2)][1L]
  expect_false(is.na(p$targets$alt2[p$targets$label == multi]))
})

test_that("a plan composes in layers with + and prints without solving", {
  p <- dummy_plan() +
    marker_targets(chrom = "Chr1", pos = 4L, ref = "T", alt = "A") +
    marker_context(variants = tibble::tibble(chrom = "Chr1", pos = 6L, alt = "G"))
  expect_s3_class(p, "marker_plan")
  expect_equal(nrow(p$targets), 1L)
  expect_output(print(p), "marker_plan")
  expect_error(p + 1, "layer")
})

test_that("legacy marker_context() calls still build a plan layer", {
  v <- tibble::tibble(chrom = "Chr1", pos = 6L, alt = "G")
  h <- tibble::tibble(
    haplotype = "H1", freq = 1, chrom = "Chr1", pos = 6L, ref = "C", alt = "G"
  )
  layer <- marker_context(variants = v, haplotypes = h)
  expect_true(inherits(layer, "marker_layer"))
  p <- dummy_plan() + layer
  expect_identical(p$variants, v)
  expect_identical(p$haplotypes, h)
})

test_that("marker_tech orders schemes cheapest-first and honours exclude", {
  p <- dummy_plan() + marker_tech(preset = "kasp")
  s <- plan_schemes(p)
  expect_identical(s[1L], "snp")
  expect_false(is.unsorted(unname(KASP_TECH_COSTS[s])))
  p2 <- dummy_plan() + marker_tech(preset = "lab", exclude = "caps")
  expect_false("caps" %in% plan_schemes(p2))
  expect_true("hrm" %in% plan_schemes(p2))
})

test_that("marker_costs overrides merge over the defaults", {
  c1 <- marker_costs(hrm = 1.2)
  expect_equal(unname(c1[["hrm"]]), 1.2)
  expect_equal(unname(c1[["snp"]]), 1)
  expect_error(marker_costs(1), "named")
})

test_that("a region expands only to its in-block variants and sweeps them", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 40L)))
  variants <- tibble::tibble(
    chrom = "Chr1", pos = c(20L, 40L, 140L),
    ref = "T", alt = c("G", "C", "A")
  )
  region <- GenomicRanges::GRanges("Chr1", IRanges::IRanges(10L, 60L))
  p <- marker_plan(genome, region = region) +
    marker_context(variants = variants) +
    marker_tech(schemes = "snp")
  expect_equal(p$region$start, 10L)
  expect_output(print(p), "region: Chr1:10-60")
  tg <- region_targets(p)
  expect_setequal(tg$pos, c(20L, 40L))
  expect_true(all(tg$type == "region"))
  # seed the cache so the sweep is checked without running a designer
  cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(tg))) {
    assign(
      paste(target_key(tg[i, , drop = FALSE]), "snp", "", sep = "||"),
      marker_doomed("test", "as", tg[i, , drop = FALSE], NULL),
      envir = cache
    )
  }
  solved <- marker_solve(p, cache = cache, progress = FALSE)
  expect_setequal(solved$report$pos, c(20L, 40L))
  expect_false(140L %in% solved$report$pos)
  # a region with no variant in it is an error, never an empty sweep
  empty <- marker_plan(
    genome,
    region = GenomicRanges::GRanges("Chr1", IRanges::IRanges(100L, 120L))
  ) + marker_context(variants = variants)
  expect_error(region_targets(empty), "No variant")
  # and a region without context is refused
  bare <- marker_plan(genome, region = region)
  expect_error(region_targets(bare), "variant context")
})

test_that("a solved region carries a haplotype -> signal table inside the block", {
  W <- withr::with_seed(4L, paste(sample(c("A", "C", "G", "T"), 500L, TRUE), collapse = ""))
  ref <- function(p) substring(W, p, p)
  alt_of <- function(p) setdiff(c("A", "C", "G", "T"), ref(p))[1L]
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 250L, ref = ref(250L), alt = alt_of(250L)
  )
  haps <- tibble::tibble(
    haplotype = c("H1", "H2", "H2", "H3", "H3", "H2"),
    freq = c(0.5, 0.3, 0.3, 0.2, 0.2, 0.3),
    chrom = "Chr1", pos = c(240L, 250L, 260L, 250L, 260L, 400L),
    ref = c(ref(240L), ref(250L), ref(260L), ref(250L), ref(260L), ref(400L)),
    alt = c("A", alt_of(250L), "A", alt_of(250L), "A", alt_of(400L))
  )
  region <- GenomicRanges::GRanges("Chr1", IRanges::IRanges(200L, 300L))
  plan <- marker_plan(genome, region = region) +
    marker_context(haplotypes = haps) +
    marker_tech(schemes = "snp")
  tg <- region_targets(plan)
  # the out-of-block variant is never a target
  expect_setequal(tg$pos, c(240L, 250L, 260L))
  # seed a haplotype-aware winner for the target and a doomed marker elsewhere
  hv <- haplotype_variants(haps, snp)
  hits <- haplotype_hits(haps, snp, as_start = 235L, as_end = 255L, alt_hap = hv$alt)
  winner <- structure(list(haplotype_hits = hits), class = "marker_assay")
  cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(tg))) {
    m <- if (tg$pos[[i]] == 250L) {
      winner
    } else {
      marker_doomed("test", "as", tg[i, , drop = FALSE], NULL)
    }
    assign(paste(target_key(tg[i, , drop = FALSE]), "snp", "", sep = "||"), m, envir = cache)
  }
  solved <- marker_solve(plan, cache = cache, progress = FALSE)
  sig <- solved$signals
  expect_true(all(c(
    "target", "haplotype", "target_allele", "primer", "signal",
    "within_phased_block"
  ) %in% names(sig)))
  expect_true(nrow(sig) > 0L)
  expect_true(all(nzchar(sig$signal)))
  expect_true(all(sig$within_phased_block))
  expect_setequal(sig$haplotype, c("H1", "H2", "H3"))
  # the artefact carries the table, with a non-empty signal column
  art <- marker_artefact(solved, dir = withr::local_tempdir())
  path <- fs::path(art$dir, "signals.csv")
  expect_true(fs::file_exists(path))
  back <- readr::read_csv(path, show_col_types = FALSE)
  expect_true("signal" %in% names(back))
  expect_true(all(nzchar(back$signal)))
  expect_identical(art$manifest$counts$signals, nrow(back))
})

test_that("marker_solve designs the real example targets and reports", {
  skip_if_no_oligoarrayaux()
  solved <- solved_example()
  skip_if(is.null(solved), "TAIR10 genome not available")
  expect_s3_class(solved, "marker_solved")
  expect_equal(nrow(solved$report), 2L)
  expect_true(all(solved$report$ok))
  expect_true(all(c("label", "scheme", "cost", "ok", "tried", "marker") %in% names(solved$report)))
  expect_output(print(solved), "marker_solved")
  expect_gte(nrow(solved$attempts), 2L)
})

test_that("marker_solve reuses its cache", {
  skip_if_no_oligoarrayaux()
  s1 <- solved_example()
  skip_if(is.null(s1), "TAIR10 genome not available")
  # re-solving with the cache the first solve used must return the same answer
  s2 <- suppressWarnings(
    marker_solve(s1$plan, flank = 150L, cores = 1L, progress = FALSE, cache = s1$cache)
  )
  expect_identical(s1$report$scheme, s2$report$scheme)
  expect_identical(s1$report$marker[[1L]], s2$report$marker[[1L]])
})

test_that("marker_primers pulls the oligos out of a solved marker", {
  skip_if_no_oligoarrayaux()
  solved <- solved_example()
  skip_if(is.null(solved), "TAIR10 genome not available")
  pm <- marker_primers(solved$report$marker[[1L]])
  expect_true(all(c("AS_ref", "AS_alt", "COM") %in% names(pm)))
  expect_true(all(nzchar(pm)))
})

test_that("marker_tubes screens every pair of solved markers", {
  skip_if_no_oligoarrayaux()
  solved <- solved_example()
  skip_if(is.null(solved), "TAIR10 genome not available")
  tb <- marker_tubes(solved)
  expect_true(all(c("a", "b", "n_pairs", "worst_run", "pass") %in% names(tb)))
  expect_equal(nrow(tb), 1L) # one pair of two solved targets
})

test_that("marker_plot and marker_export work off a solved plan", {
  skip_if_no_oligoarrayaux()
  skip_if_not_installed("ggplot2")
  solved <- solved_example()
  skip_if(is.null(solved), "TAIR10 genome not available")
  expect_s3_class(marker_plot(solved), "ggplot")
  dir <- withr::local_tempdir()
  marker_export(solved, dir = dir)
  for (f in c(
    "markers.csv", "attempts.csv", "primers.csv", "solved.rds",
    # the artefact sidecars: the manifest, the decision log, the gate table
    "artefact.json", "decisions.csv", "checks.csv", "targets.csv"
  )) {
    expect_true(fs::file_exists(fs::path(dir, f)), info = f)
  }
  expect_true("solved.rds" %in% list.files(dir))
})

# --- the artefact (YP3): provenance, decisions, replay ----------------------
# A synthetic poly-A reference, so the rejection path runs without OligoArrayAux,
# the TAIR10 download or a full chemistry ladder: every candidate is doomed by
# the GC filter, which is exactly the row a decision log must never leave blank.
rejected_solve <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("A", 600L)))
      plan <- marker_plan(genome) +
        marker_targets(chrom = "Chr1", pos = 300L, ref = "A", alt = "G") +
        marker_tech(schemes = "snp")
      value <<- suppressWarnings(
        marker_solve(plan, min_disc = 0.4, progress = FALSE)
      )
    }
    value
  }
})

test_that("a knob that differs from the designer default is its own manifest field", {
  ov <- artefact_overrides(list(min_disc = 0.4, flank = 150L, cores = 1L), "snp")
  expect_equal(ov$argument, "min_disc")
  expect_equal(ov$default[[1L]], 0.5)
  expect_true(ov$overridden[[1L]])
  # a knob is compared only against the scheme that declares it: `probe_gc` is
  # a dPCR argument, so KASP is not reported as having left it alone
  one <- artefact_overrides(list(probe_gc = 0.5), c("snp", "dpcr"))
  expect_equal(unique(one$scheme), "dpcr")
  # a knob no designer declares is kept, with an unknown default, so a typo
  # shows up in the manifest instead of silently changing nothing
  odd <- artefact_overrides(list(not_a_knob = 1), "snp")
  expect_equal(odd$argument, "not_a_knob")
  expect_true(is.na(odd$default[[1L]]))
  expect_true(is.na(odd$scheme[[1L]]))
})

test_that("the artefact records its inputs, its id and every rejection's reason", {
  solved <- rejected_solve()
  art <- marker_artefact(solved, dir = withr::local_tempdir())
  expect_s3_class(art, "marker_artefact")
  expect_match(art$id, "^[0-9a-f]{64}$")
  m <- art$manifest
  expect_identical(m$schema, "tidymarkers/artefact@1")
  expect_equal(m$counts$targets, 1L)
  expect_equal(m$counts$rejected, 1L)
  # the reference is named and hashed, and an installed engine is never dropped
  expect_equal(m$inputs[[1L]]$role, "genome")
  expect_true(nzchar(m$inputs[[1L]]$sha256))
  expect_identical(names(m$engine_stamp), marker_engine_names)
  # every rejected candidate says why, and at which stage
  dec <- art$decisions
  expect_true(nrow(dec) >= 1L)
  rej <- dec[dec$outcome == "rejected", , drop = FALSE]
  expect_true(all(nzchar(rej$reason)))
  expect_true(all(nzchar(rej$stage)))
  # the override travels as a field of its own, not buried in `parameters`
  ov <- jsonlite::fromJSON(fs::path(art$dir, "artefact.json"), simplifyVector = FALSE)$overrides
  expect_equal(purrr::map_chr(ov, "argument"), "min_disc")
  expect_equal(purrr::map_dbl(ov, "value"), 0.4)
})

test_that("a decision row is chosen only for its own target's winner", {
  # the attempt matrix is longer than the report (one row per target), so the
  # `chosen` flag must be built per target: a recycled report vector would mark
  # a different target's row
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 40L)))
  plan <- marker_plan(genome) +
    marker_targets(chrom = "Chr1", pos = 20L, ref = "T", alt = "G") +
    marker_targets(chrom = "Chr1", pos = 40L, ref = "T", alt = "C") +
    marker_tech(schemes = c("snp", "caps"))
  tg <- plan$targets
  cache <- new.env(parent = emptyenv())
  winner <- structure(list(scheme = "snp"), class = "marker_assay")
  assign(paste(target_key(tg[1L, , drop = FALSE]), "snp", "", sep = "||"), winner, envir = cache)
  assign(
    paste(target_key(tg[2L, , drop = FALSE]), "snp", "", sep = "||"),
    marker_doomed("no AS primer passed", "as", tg[2L, , drop = FALSE], NULL),
    envir = cache
  )
  solved <- structure(
    list(
      report = tibble::tibble(
        label = tg$label, chrom = tg$chrom, pos = tg$pos, ref = tg$ref,
        alt = tg$alt, type = tg$type, scheme = c("snp", NA_character_),
        cost = c(1, NA), ok = c(TRUE, FALSE), tried = c("snp", "snp, caps"),
        marker = list(winner, NULL)
      ),
      attempts = tibble::tibble(
        label = c(tg$label[1L], tg$label[2L], tg$label[2L]),
        scheme = c("snp", "snp", "caps"), cost = c(1, 1, 10),
        ok = c(TRUE, FALSE, FALSE)
      ),
      plan = plan, schemes = c("snp", "caps"), cache = cache, knobs = list()
    ),
    class = "marker_solved"
  )
  dec <- artefact_decisions(solved)
  expect_equal(nrow(dec), 3L)
  expect_equal(sum(dec$chosen), 1L)
  expect_identical(dec$target[dec$chosen], as.character(tg$label[1L]))
  expect_true(all(nzchar(dec$reason[!dec$chosen])))
  expect_identical(dec$stage[!dec$chosen][[1L]], "as")
})

test_that("marker_replay re-runs the design and compares the verdicts", {
  solved <- rejected_solve()
  art <- marker_artefact(solved, dir = withr::local_tempdir())
  # a fresh cache: a replay that reads the original cache tests nothing
  rp <- suppressWarnings(marker_replay(art))
  expect_s3_class(rp, "marker_replay")
  expect_true(rp$inputs_match)
  expect_true(rp$match)
  expect_equal(nrow(rp$table), 1L)
  expect_identical(rp$table$recorded_verdict, rp$table$replayed_verdict)
  expect_output(print(rp), "marker_replay")
  # the id is a property of the design, not of when it was written
  again <- marker_artefact(solved, dir = withr::local_tempdir())
  expect_identical(art$id, again$id)
  # YP5: an outcome is bound to that id, not to a file path
  out <- marker_record_outcome(art, "fails", notes = "synthetic reference")
  expect_identical(out$id, art$id)
  expect_false(any(out$calibrated))
})
