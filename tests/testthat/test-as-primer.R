# Invariants for the KASP core: the SNP must land on the 3'-terminal base, and
# the DECIPHER "target" must be the reverse complement of the primer.  A silent
# swap here would invert the discrimination metric, so it is tested explicitly.

W     <- Biostrings::DNAString(paste(rep("ACGT", 19L), collapse = ""))  # 76 bp
SNP_I <- 30L

seg_with_allele <- function(L, strand, allele) {
  seg <- as.character(as_window(W, SNP_I, L, strand))
  idx <- if (identical(strand, "+")) L else 1L
  substr(seg, idx, idx) <- allele
  Biostrings::DNAString(seg)
}

test_that("make_as_primer puts the requested allele at the 3' end on both strands", {
  L <- 22L
  for (st in c("+", "-")) {
    for (allele in c("A", "C")) {
      p <- make_as_primer(W, SNP_I, L, allele, strand = st)
      # "+": the primer carries the plus-strand base itself;
      # "-": the primer binds the plus strand, so its 3' base is complementary.
      expect_equal(nchar(p), L)
      expect_equal(substr(p, L, L),
                   if (identical(st, "+")) allele else
                     as.character(Biostrings::complement(Biostrings::DNAString(allele))))
    }
  }
})

test_that("the DECIPHER target is the reverse complement of the primer (no stabiliser)", {
  L <- 21L
  for (st in c("+", "-")) {
    for (allele in c("A", "C", "G")) {
      p   <- make_as_primer(W, SNP_I, L, allele, strand = st)
      tgt <- as.character(as_target(seg_with_allele(L, st, allele), st))
      expect_equal(as.character(Biostrings::reverseComplement(Biostrings::DNAString(p))), tgt)
    }
  }
})

test_that("stabilising mismatch replaces the intended base but keeps the 3' allele", {
  L <- 24L
  p0  <- make_as_primer(W, SNP_I, L, "A", strand = "+")
  b0  <- substr(p0, L - 3L + 1L, L - 3L + 1L)
  alt <- setdiff(c("A", "C", "G", "T"), b0)[1L]
  p1  <- make_as_primer(W, SNP_I, L, "A", strand = "+", stab_at = 3L, stab_base = alt)
  expect_equal(nchar(p1), L)
  expect_equal(substr(p1, L - 3L + 1L, L - 3L + 1L), alt)
  expect_equal(substr(p1, L, L), "A")
})
