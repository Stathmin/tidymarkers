# Batch design: one row per locus with a verdict; a bad locus must not abort.

bseq <- function(n, seed) {
  set.seed(seed)
  paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
}

test_that("variant_type labels SNV, MNP and indels", {
  expect_equal(variant_type("A", "G"), "SNP")
  expect_equal(variant_type("AC", "GT"), "MNP")
  expect_equal(variant_type("A", "AGG"), "insertion")
  expect_equal(variant_type("AGG", "A"), "deletion")
})

test_that("design_kasp_markers returns a tidy per-locus verdict without aborting", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  W <- Biostrings::DNAString(bseq(500, 4))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  cand <- seq(160L, 340L, by = 1L)
  clean <- cand[vapply(cand, function(p)
    !has_run(as.character(Biostrings::subseq(W, p - 35L, p + 35L))), logical(1))]
  p1 <- clean[1L]; p2 <- clean[which(abs(clean - p1) > 60L)[1L]]
  ref1 <- as.character(Biostrings::subseq(W, p1, p1))
  alt1 <- setdiff(c("A", "C", "G", "T"), ref1)[1L]
  ref2 <- as.character(Biostrings::subseq(W, p2, p2))
  loci <- tibble::tibble(
    chrom = c("Chr1", "Chr1", "Chr1"),
    pos   = c(p1, p2, 5L),                                   # last one is un-designable
    ref   = c(ref1, ref2, "A"),
    alt   = c(alt1, paste0(ref2, "GG"), "G"))
  res <- design_kasp_markers(genome, loci, flank = 150L)
  expect_equal(nrow(res), 3L)
  expect_true(all(c("chrom","pos","ref","alt","type","verdict","note") %in% names(res)))
  expect_equal(res$type, c("SNP", "insertion", "SNP"))
  # the un-designable edge locus is captured as a failure, not thrown ...
  expect_equal(res$verdict[3L], "fail")
  expect_true(is.na(res$strand[3L]))
  expect_true(nzchar(res$note[3L]))
  # ... and the verdict wiring is consistent: an assay exists iff the locus did
  # not fail, with finite coordinates and discrimination
  built <- res$verdict != "fail"
  expect_gte(sum(built), 1L)
  expect_false(any(is.na(res$strand[built])))
  expect_true(all(is.finite(res$discrimination[built])))
  expect_false(any(is.na(res$amplicon[built])))
})

test_that("a BiocParallel backend gives the same batch verdicts", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  skip_if_not_installed("BiocParallel")
  W <- Biostrings::DNAString(bseq(400, 4))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- as.character(Biostrings::subseq(W, 210L, 210L))
  loci <- tibble::tibble(chrom = "Chr1", pos = 210L, ref = ref,
                         alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
  serial <- suppressWarnings(design_kasp_markers(genome, loci, flank = 150L))
  par    <- suppressWarnings(design_kasp_markers(genome, loci, flank = 150L,
                                                 BPPARAM = BiocParallel::SerialParam()))
  expect_equal(par, serial)
})
