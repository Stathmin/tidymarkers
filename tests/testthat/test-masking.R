# Masking flanking variants: the helper, and the hard footprint constraint.

test_that("mask_positions blanks the requested bases and ignores out-of-range", {
  expect_equal(as.character(mask_positions("ACGTACGT", c(2L, 5L))), "ANGTNCGT")
  expect_equal(as.character(mask_positions("ACGT", integer(0))), "ACGT")
  expect_equal(as.character(mask_positions("ACGT", c(0L, 99L))), "ACGT")  # ignored
})

test_that("max_footprint_variants enforces N-masking as a hard constraint", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  set.seed(4)
  W <- Biostrings::DNAString(paste(sample(c("A", "C", "G", "T"), 400L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  clean <- seq(160L, 240L); clean <- clean[vapply(clean, function(p)
    !has_run(as.character(Biostrings::subseq(W, p - 35L, p + 35L))), logical(1))]
  pos <- clean[1L]
  ref <- as.character(Biostrings::subseq(W, pos, pos))
  snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = ref,
                        alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
  # decoys on BOTH flanks, so every AS-primer orientation covers one
  variants <- tibble::tibble(chrom = "Chr1", pos = c(pos, pos - 3L, pos + 3L),
                             ref = c(ref, "A", "A"), alt = c(snp$alt, "C", "T"))

  a <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = variants, flank = 150L))
  expect_true(all(a$as_candidates$n_var >= 1L))          # variants are seen under the footprint

  expect_error(
    suppressWarnings(assemble_kasp_assay(genome, snp, variants = variants, flank = 150L,
                                         max_footprint_variants = 0L)),
    "max_footprint_variants")
})
