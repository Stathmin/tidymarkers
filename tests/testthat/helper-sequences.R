# Shared fixtures for tests that need a reproducible DNA sequence.

# A reproducible random ACGT string of length `n` from a fixed `seed`,
# GC-balanced on average - the workhorse template for design tests.
bseq <- function(n, seed) {
  withr::with_seed(seed, paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = ""))
}
