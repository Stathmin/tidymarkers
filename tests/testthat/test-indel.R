# Indel engine: allele templates, the junction scan, and end-to-end assembly.

rand_seq <- function(n, seed) {
  withr::with_seed(seed, paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = ""))
}

test_that("allele_templates keeps the left flank and shifts the right one", {
  W <- Biostrings::DNAString(rand_seq(120, 7))
  snp_i <- 60L
  ref <- as.character(Biostrings::subseq(W, snp_i, snp_i + 2L)) # e.g. deletion block
  alt <- substr(ref, 1L, 1L)
  at <- allele_templates(W, snp_i, ref, alt)
  expect_equal(length(at$alt), length(W) + at$delta)
  expect_equal(
    as.character(Biostrings::subseq(at$alt, 1L, snp_i - 1L)),
    as.character(Biostrings::subseq(W, 1L, snp_i - 1L))
  )
  # the ALT block now sits where the REF block was
  expect_equal(as.character(Biostrings::subseq(at$alt, snp_i, snp_i)), alt)
})

test_that("the junction scan proposes discriminating primers for an insertion", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(rand_seq(160, 11))
  snp_i <- 80L
  ref <- as.character(Biostrings::subseq(W, snp_i, snp_i))
  alt <- paste0(ref, "TTG") # 2-nt insertion, VCF-anchored
  d <- design_as_primers_indel(W, snp_i, ref, alt, lengths = 20:22)
  expect_setequal(unique(d$allele), c("ref", "alt"))
  expect_gt(max(d$discrimination), 0)
  expect_true(all(c("primer", "on_eff", "off_eff", "as_start", "as_end") %in% names(d)))
})

test_that("the junction scan proposes discriminating primers for a deletion", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(rand_seq(160, 13))
  snp_i <- 80L
  ref <- as.character(Biostrings::subseq(W, snp_i, snp_i + 2L))
  alt <- substr(ref, 1L, 1L)
  d <- design_as_primers_indel(W, snp_i, ref, alt, lengths = 20:22)
  expect_gt(max(d$discrimination), 0)
})

test_that("design_as_primers_indel refuses equal-length alleles", {
  withr::local_envvar(LANGUAGE = "en")
  W <- Biostrings::DNAString(rand_seq(160, 17))
  expect_error(design_as_primers_indel(W, 80L, "A", "G"), "unequal")
})

test_that("assemble_kasp_assay designs an indel end-to-end", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(rand_seq(400, 4))
  # pick a locus whose neighbourhood has no long homopolymer (a real designer
  # would too - a run-rich junction is legitimately rejected by the checks)
  cand <- seq(160L, 240L, by = 1L)
  ok <- vapply(
    cand, function(p) !has_run(as.character(Biostrings::subseq(W, p - 35L, p + 35L))),
    logical(1)
  )
  pos <- cand[which(ok)[1L]]
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- as.character(Biostrings::subseq(W, pos, pos))
  snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = ref, alt = paste0(ref, "GGG"))
  a <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = NULL, flank = 150L))
  expect_equal(a$snp$ref, ref)
  expect_false(nchar(a$snp$ref) == nchar(a$snp$alt)) # indel => alleles differ in length
  expect_false(a$assay$as_ref[1] == a$assay$as_alt[1]) # the two AS primers differ
  expect_true(all(c("as_start", "as_alt_start") %in% names(a$assay)))
  chk <- kasp_design_checks(a)
  expect_true(all(c("criterion", "observed", "limit", "pass") %in% names(chk)))
  expect_true(is.logical(chk$pass))
})
