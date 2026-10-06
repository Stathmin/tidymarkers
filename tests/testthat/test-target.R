# The scheme-aware target pick: "doomed" depends on the marker type.

test_that("is_target_designable covers every marker class", {
  expect_false(is_target_designable(NULL))
  expect_false(is_target_designable(structure(list(), class = "kasp_doomed")))
  expect_true(is_target_designable(structure(list(), class = "kasp_assay")))
  expect_true(is_target_designable(structure(list(), class = "kasp_caps")))
  expect_true(is_target_designable(structure(list(), class = "kasp_arms")))
})

test_that("kasp_design_at dispatches by scheme and is NULL for unsupported ones", {
  expect_null(kasp_design_at(NULL, NULL, "asq"))
  expect_null(kasp_design_at(NULL, NULL, "primer"))
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
  # "asq" has no engine in this build, so nothing designs: warn + best-ranked
  expect_warning(
    tgt <- pick_target_snp(variants, genome, gene, scheme = "asq"),
    "designable for scheme asq"
  )
  expect_equal(tgt$id, "v")
  # without a scheme the same locus is the plain pick
  expect_equal(pick_target_snp(variants, genome, gene)$id, "v")
})

test_that("design_marker walks the chemistries and returns the first that works", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  d <- suppressWarnings(design_marker(genome, snp, flank = 150L))
  expect_setequal(d$tried$scheme, c("snp", "caps", "arms", "hrm"))
  expect_true(d$scheme %in% c("snp", "caps", "arms", "hrm"))
  expect_true(is_target_designable(d$marker))
  # restricted to HRM, it reports HRM (or nothing) and says so
  h <- suppressWarnings(design_marker(genome, snp, schemes = "hrm", flank = 150L))
  expect_equal(h$tried$scheme, "hrm")
})
