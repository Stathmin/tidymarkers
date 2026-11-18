# Batching the whole-genome passes of `assemble_kasp_assay()`: the counted
# primitive is `pattern_hits_batch()` (one call = one grouped genome pass), so
# the test pins the pass count and that the batched paths answer what the
# per-pair calls they replace answer.

# the tracer evaluates inside the traced function's frames, so the counter lives
# in the package's own (unlocked) state environment
with_pass_counter <- function(code) {
  ns <- asNamespace("tidymarkers")
  state <- get(".marker_state", envir = ns)
  counter <- new.env()
  counter$n <- 0L
  state$pass_counter <- counter
  on.exit(
    {
      untrace("pattern_hits_batch", where = ns)
      rm("pass_counter", envir = state)
    },
    add = TRUE
  )
  trace("pattern_hits_batch",
    tracer = quote(.marker_state$pass_counter$n <- .marker_state$pass_counter$n + 1L),
    print = FALSE, where = ns
  )
  force(code)
  counter$n
}

test_that("one KASP design issues a fixed, batched number of genome passes", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  # the first design warms the per-pattern hit cache, so the counted one is the
  # steady state (a cold run issues more passes, never fewer)
  m <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp))
  n <- with_pass_counter(suppressWarnings(assemble_kasp_assay(f$genome, f$snp)))
  # the pair screen (exact, then one near pool), the final assay PCR's exact
  # pool and the ideogram's one grouped scan; one pass per pairing would be far
  # more.
  expect_lte(n, 8L)
})

test_that("the pooled virtual PCR equals the per-pair calls it replaces", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  m <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp))
  pairs <- as_pcr_pairs(m)
  pooled <- in_silico_pcr_pool(pairs, f$genome)
  one <- lapply(pairs, function(p) {
    in_silico_pcr(p[["F"]], p[["R"]], f$genome, max_mismatch = 0L)
  })
  expect_identical(pooled, stats::setNames(one, names(pooled)))
})

test_that("the batched ideogram hits equal the per-primer scan", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  m <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp))
  primers <- c(AS_ref = m$assay$as_ref[1], AS_alt = m$assay$as_alt[1], COM = m$assay$com[1])
  per_primer <- primer_hits(primers, f$genome)
  expect_identical(m$ideogram$hits[c("primer", "chrom", "pos")], per_primer)
  # one grouped pass answers all of them (the two allele-specific primers are two
  # distinct patterns, so a pass is per distinct pattern, not per primer)
  n <- with_pass_counter(suppressWarnings(batch_primer_hits(primers, f$genome)))
  expect_equal(n, 1L)
})
