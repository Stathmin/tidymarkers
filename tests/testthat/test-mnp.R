# Substitution blocks (MNP): the same 3'-anchored engine must work when REF/ALT
# are strings longer than one base, and the anchor must follow the block.

test_that("as_anchor follows the variant block per strand", {
  expect_equal(as_anchor(100L, 1L, "+"), 100L)   # SNP reduces to the old case
  expect_equal(as_anchor(100L, 1L, "-"), 100L)
  expect_equal(as_anchor(100L, 3L, "+"), 102L)   # plus strand: block's last base
  expect_equal(as_anchor(100L, 3L, "-"), 100L)   # minus strand: block's first base
})

test_that("block_in_seg marks the 3' (plus) or 5' (minus) end of the window", {
  expect_equal(c(IRanges::start(block_in_seg(3L, 20L, "+")), IRanges::end(block_in_seg(3L, 20L, "+"))), c(18L, 20L))
  expect_equal(c(IRanges::start(block_in_seg(3L, 20L, "-")), IRanges::end(block_in_seg(3L, 20L, "-"))), c(1L, 3L))
})

test_that("make_as_primer inserts the whole MNP allele at the 3' end", {
  W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))  # 80 bp
  ref <- "GTA"; alt <- "ACG"; L <- 22L
  p_plus <- make_as_primer(W, as_anchor(39L, 3L, "+"), L, ref, strand = "+")
  expect_equal(substr(p_plus, L - 2L, L), ref)             # block at 3' end, plus strand
  p_minus <- make_as_primer(W, as_anchor(39L, 3L, "-"), L, alt, strand = "-")
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(alt)))
  expect_equal(substr(p_minus, L - 2L, L), rc)             # 3' block is the complement
})

test_that("design_as_primers discriminates a 3-nt substitution block", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))  # 80 bp
  d <- design_as_primers(W, 39L, ref = "GTA", alt = "ACG",
                         lengths = 20:22, stab_ats = NA)
  expect_gt(nrow(d), 0L)
  expect_gt(max(d$discrimination), 0)                      # a real block mismatch discriminates
  expect_gt(d$ref_on[1L], d$ref_off[1L])
})

test_that("pick_target_snp accepts an MNP and keeps REF/ALT strings", {
  seq <- "TTTTTTTTTTTTTTTTTTGTACCCCCCCCCCCCCCCCCCCCCCCCCC"
  genome <- Biostrings::DNAStringSet(c(Chr1 = seq))
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L)
  variants <- tibble::tibble(
    id = "mnp1", chrom = "Chr1", pos = 19L,
    ref_genome = list(c("GTA", "ACG")), consequence = "missense")
  tgt <- pick_target_snp(variants, genome, gene)
  expect_equal(tgt$ref, "GTA")
  expect_equal(tgt$alt, "ACG")
  expect_equal(nchar(tgt$ref), nchar(tgt$alt))
})

test_that("pick_target_snp skips indels and ambiguous alleles", {
  seq <- "ACGTACGTACGTACGTACGTACGTACGT"
  genome <- Biostrings::DNAStringSet(c(Chr1 = seq))
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L)
  variants <- tibble::tibble(
    id = c("del", "amb"),
    chrom = "Chr1", pos = c(5L, 9L),
    ref_genome = list(c("AC", "A"), c("A", "N")),
    consequence = "missense")
  expect_error(pick_target_snp(variants, genome, gene), "No usable substitution")
})
