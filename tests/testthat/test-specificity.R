# Specificity seed counting: a known seed must be counted on both strands.

test_that("specificity_seed counts seed occurrences on both strands", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  primers <- c(P = "TTTTTTAAACCCGGG") # 3'-terminal 9 nt = AAACCCGGG
  # forward "AAACCCGGG" x2 + reverse-complement "CCCGGGTTT" x2 = 4
  res <- specificity_seed(primers, genome, k = 9L, max_mismatch = 0L)
  expect_equal(res$hits, 4L)
  expect_equal(res$seed, "AAACCCGGG")
})

test_that("specificity_seed returns zero hits for an absent seed", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  primers <- c(P = "TTTTTTACGACGAC")
  expect_equal(specificity_seed(primers, genome, k = 9L, max_mismatch = 0L)$hits, 0L)
})

test_that("the mismatch scan is invariant to the parallel backend", {
  skip_if_not_installed("BiocParallel")
  genome <- Biostrings::DNAStringSet(c(
    chr1 = "AAACCCGGGTTTAAACCCGGGTTT",
    chr2 = "GGGTTTAAACCC"
  ))
  primers <- c(P = "TTTTTTAAACCCGGG", Q = "TTTTTTGGGAAACCC")
  serial <- specificity_seed(primers, genome, k = 9L, max_mismatch = 1L)
  par <- specificity_seed(primers, genome,
    k = 9L, max_mismatch = 1L,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_equal(par$hits, serial$hits)
  skip_on_os("windows")
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  mc <- specificity_seed(primers, genome,
    k = 9L, max_mismatch = 1L,
    BPPARAM = bp
  )
  expect_equal(mc$hits, serial$hits)
})

test_that("specificity counts the reference genome, not the variant/allele list", {
  # the seed is looked up in the reference genome only, so a position that merely
  # *varies* across accessions (a row in the variant table) is one reference locus,
  # not a phantom repeat; a genuine duplication in the genome is what adds a hit.
  seed <- "GTCAGTCC"
  primer <- paste0("AAAA", seed)
  one <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTGTCAGTCCACGTACGT"))
  two <- Biostrings::DNAStringSet(c(chr1 = "GTCAGTCCACGTACGTGTCAGTCC"))
  expect_equal(specificity_seed(c(P = primer), one, k = 8L, max_mismatch = 0L)$hits, 1L)
  expect_equal(specificity_seed(c(P = primer), two, k = 8L, max_mismatch = 0L)$hits, 2L)
  # the specificity API takes a genome, never a variant table - so it cannot
  # double-count the alternative alleles of a multi-variant feature
  expect_false("variants" %in% names(formals(specificity_seed)))
})

test_that("the exact PDict path matches a manual both-strand count", {
  genome <- Biostrings::DNAStringSet(c(
    chr1 = "AAACCCGGGTTTAAACCCGGGTTT",
    chr2 = "GGGAAACCC"
  ))
  primers <- c(P = "TTTTTTAAACCCGGG", Q = "TTTTTTGGGAAACCC")
  res <- specificity_seed(primers, genome, k = 9L, max_mismatch = 0L)
  manual <- vapply(unname(substr(primers, 8L, 16L)), function(s) {
    d <- Biostrings::DNAString(s)
    r <- Biostrings::reverseComplement(d)
    sum(Biostrings::vcountPattern(d, genome)) + sum(Biostrings::vcountPattern(r, genome))
  }, numeric(1))
  expect_equal(res$hits, as.integer(manual))
})

test_that("a BiocParallel BPPARAM gives the same mismatch counts as serial", {
  skip_if_not_installed("BiocParallel")
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  primers <- c(P = "TTTTTTAAACCCGGG")
  serial <- specificity_seed(primers, genome, k = 9L, max_mismatch = 1L)
  par <- specificity_seed(primers, genome,
    k = 9L, max_mismatch = 1L,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_equal(par$hits, serial$hits)
})
