# The scheme-aware target pick: "doomed" depends on the marker type.

test_that("is_target_designable covers every marker class", {
  expect_false(is_target_designable(NULL))
  expect_false(is_target_designable(structure(list(), class = c("marker_doomed", "marker"))))
  expect_true(is_target_designable(structure(list(), class = c("marker_assay", "marker"))))
  expect_true(is_target_designable(structure(list(), class = c("marker_caps", "marker"))))
  expect_true(is_target_designable(structure(list(), class = c("marker_arms", "marker"))))
})

test_that("marker_design_at dispatches by scheme and is NULL for unsupported ones", {
  expect_null(marker_design_at(NULL, NULL, "primer"))
})

test_that("pick_target_snp with a scheme warns and returns the best when none designs", {
  withr::local_envvar(LANGUAGE = "en")
  seq <- strrep("ACGT", 20L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = seq))
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L)
  variants <- tibble::tibble(
    id = "v", chrom = "Chr1", pos = 20L,
    ref_genome = list(c("T", "A")), consequence = "missense_variant"
  )
  # "primer" has no designer, so nothing designs: warn + best-ranked
  expect_warning(
    tgt <- pick_target_snp(variants, genome, gene, scheme = "primer"),
    "designable for scheme primer"
  )
  expect_equal(tgt$id, "v")
  # without a scheme the same locus is the plain pick
  expect_equal(pick_target_snp(variants, genome, gene)$id, "v")
})

test_that("the target pick names its whole tuning counts before the scan", {
  seq <- strrep("ACGT", 20L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = seq))
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L)
  variants <- tibble::tibble(
    id = "v", chrom = "Chr1", pos = 20L,
    ref_genome = list(c("T", "A")), consequence = "missense_variant"
  )
  # IRanges cast `clean_3p` and `as.integer()` cast `max_tries`, so 2.9 built
  # the 2-base guard and 0.5 still ran one attempt; a fractional
  # `specificity_candidates` reached dplyr's own "`n` must be a round number"
  # (ADR 0061)
  expect_error(
    pick_target_snp(variants, genome, gene, specificity_candidates = 0L, clean_3p = 2.9),
    "`clean_3p` must be one positive whole number of nucleotides"
  )
  expect_error(
    pick_target_snp(variants, genome, gene, specificity_candidates = 1.9),
    "`specificity_candidates` must be one non-negative whole number"
  )
  expect_error(
    pick_target_snp(variants, genome, gene, specificity_candidates = 0L, max_tries = 0.5),
    "`max_tries` must be one positive whole number"
  )
  expect_error(
    pick_target_snp(variants, genome, gene, specificity_candidates = 0L, max_tries = Inf),
    "`max_tries` must be one positive whole number"
  )
  expect_equal(pick_target_snp(variants, genome, gene)$id, "v")
})

test_that("design_auto_marker walks the chemistries and returns the first that works", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  genome <- f$genome
  snp <- f$snp
  d <- suppressWarnings(design_auto_marker(genome, snp, flank = 150L))
  expect_setequal(d$tried$scheme, c("snp", "caps", "arms", "hrm"))
  expect_true(d$scheme %in% c("snp", "caps", "arms", "hrm"))
  expect_true(is_target_designable(d$marker))
  # restricted to HRM, it reports HRM (or nothing) and says so
  h <- suppressWarnings(design_auto_marker(genome, snp, schemes = "hrm", flank = 150L))
  expect_equal(h$tried$scheme, "hrm")
})

test_that("a context does not change dispatch or the no-context result", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "T", alt = "A")
  ctx <- marker_run_context(genome)
  # "primer" has no designer: both calls report the same null verdict
  expect_null(marker_design_at(NULL, NULL, "primer"))
  expect_null(marker_design_at(NULL, NULL, "primer", context = ctx))
  a <- suppressWarnings(design_auto_marker(genome, snp, schemes = "primer"))
  b <- suppressWarnings(design_auto_marker(genome, snp, schemes = "primer", context = ctx))
  expect_identical(a, b)
  expect_false(a$tried$ok)
})

test_that("design_marker is a deprecated alias of design_auto_marker", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "T", alt = "A")
  expect_warning(d <- design_marker(genome, snp, schemes = "primer"), "deprecated")
  expect_identical(d, suppressWarnings(design_auto_marker(genome, snp, schemes = "primer")))
})

test_that("design_markers accepts a context and keeps its result", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  loci <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "T", alt = "A")
  ctx <- marker_run_context(genome)
  a <- design_markers(genome, loci, scheme = "primer")
  b <- design_markers(genome, loci, scheme = "primer", context = ctx)
  expect_identical(a, b)
  expect_false(b$ok)
})

test_that("a plan layer cannot masquerade as a design context", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "T", alt = "A")
  layer <- marker_context(variants = tibble::tibble(chrom = "Chr1", pos = 1L, alt = "A"))
  expect_error(
    marker_design_at(genome, snp, "asq", context = layer),
    "marker_run_context\\(genome"
  )
})

# Ensembl answers 403, not only 429, when a client exceeds its rate limit; a
# retry set that omits it turns a throttle into a hard failure. Offline: the
# predicate reads only the status code, so a bare response object suffices.
test_that("ensembl_transient retries the statuses Ensembl refuses a client with", {
  expect_true(ensembl_transient(httr2::response(403L)))
  expect_true(ensembl_transient(httr2::response(408L)))
  expect_true(ensembl_transient(httr2::response(429L)))
  expect_true(ensembl_transient(httr2::response(500L)))
  expect_true(ensembl_transient(httr2::response(503L)))
  expect_false(ensembl_transient(httr2::response(200L)))
  expect_false(ensembl_transient(httr2::response(404L)))
  expect_false(ensembl_transient(httr2::response(400L)))
})

# A cache file is written by the call, so a link at the cache name is refused
# rather than followed: before the fix the response replaced the outside file.
test_that("the Ensembl cache refuses to write through a link", {
  httr2::local_mocked_responses(function(req) {
    httr2::response(200L, body = charToRaw("{}"))
  })
  cache <- withr::local_tempdir()
  calls <- list(
    "arabidopsis_thaliana_AT4G00650.json" = function() {
      ensembl_gene("AT4G00650", cache = cache)
    },
    "arabidopsis_thaliana_var_4:268901-271503.json" = function() {
      ensembl_variants("Chr4", 268901, 271503, cache = cache)
    }
  )
  for (name in names(calls)) {
    outside <- fs::path(withr::local_tempdir(), name)
    skip_if(
      !file.symlink(outside, fs::path(cache, name)),
      "this filesystem cannot make a symlink"
    )
    expect_error(calls[[name]](), "which is a link", info = name)
    expect_false(fs::file_exists(outside), info = name)
  }
})

# The read side of the same rule: a linked name is not a cache hit, so a valid
# file behind the link is refused rather than parsed as the Ensembl response.
test_that("the Ensembl cache does not read a hit through a link", {
  cache <- withr::local_tempdir()
  hits <- list(
    "arabidopsis_thaliana_AT4G00650.json" = function() {
      ensembl_gene("AT4G00650", cache = cache)
    },
    "arabidopsis_thaliana_var_4:268901-271503.json" = function() {
      ensembl_variants("Chr4", 268901, 271503, cache = cache)
    }
  )
  for (name in names(hits)) {
    outside <- fs::path(withr::local_tempdir(), name)
    writeLines('{"id": "foreign"}', outside)
    skip_if(
      !file.symlink(outside, fs::path(cache, name)),
      "this filesystem cannot make a symlink"
    )
    expect_error(hits[[name]](), "which is a link", info = name)
    expect_equal(readLines(outside), '{"id": "foreign"}', info = name)
    fs::file_delete(fs::path(cache, name))
  }
})
