# Shared fixtures for tests that need a reproducible DNA sequence.

# A reproducible random ACGT string of length `n` from a fixed `seed`,
# GC-balanced on average - the workhorse template for design tests.
bseq <- function(n, seed) {
  withr::with_seed(seed, paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = ""))
}

# The shared one-SNP fixture: a reproducible `n`-bp template with a substitution
# at `pos`, plus the `DNAStringSet` and the VCF-style locus row that every
# single-locus designer test builds by hand.  `W`/`ref`/`alt` are returned for
# the tests that need the raw template or the alleles.
snp_fixture <- function(n = 600L, seed = 2L, pos = 300L) {
  W <- bseq(n, seed)
  ref <- substr(W, pos, pos)
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  snp <- tibble::tibble(chrom = "Chr1", pos = as.integer(pos), ref = ref, alt = alt)
  list(W = W, genome = genome, snp = snp, ref = ref, alt = alt)
}

# The shared KASP primer trio: the same three synthetic oligos every assay mock
# uses (AS_ref / AS_alt / COM).
kasp_primers <- function() {
  c(
    AS_ref = "CATGGACCTAGCATTGACGTA",
    AS_alt = "CATGGACCTAGCATTGACGTT",
    COM = "GCTTAGCCATGGTCAATGCTA"
  )
}
