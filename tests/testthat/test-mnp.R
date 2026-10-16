# Substitution blocks (MNP): the same 3'-anchored engine must work when REF/ALT
# are strings longer than one base, and the anchor must follow the block.

test_that("as_anchor follows the variant block per strand", {
  expect_equal(as_anchor(100L, 1L, "+"), 100L) # SNP reduces to the old case
  expect_equal(as_anchor(100L, 1L, "-"), 100L)
  expect_equal(as_anchor(100L, 3L, "+"), 102L) # plus strand: block's last base
  expect_equal(as_anchor(100L, 3L, "-"), 100L) # minus strand: block's first base
})

test_that("block_in_seg marks the 3' (plus) or 5' (minus) end of the window", {
  expect_equal(c(IRanges::start(block_in_seg(3L, 20L, "+")), IRanges::end(block_in_seg(3L, 20L, "+"))), c(18L, 20L))
  expect_equal(c(IRanges::start(block_in_seg(3L, 20L, "-")), IRanges::end(block_in_seg(3L, 20L, "-"))), c(1L, 3L))
})

test_that("make_as_primer inserts the whole MNP allele at the 3' end", {
  W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = "")) # 80 bp
  ref <- "GTA"
  alt <- "ACG"
  L <- 22L
  p_plus <- make_as_primer(W, as_anchor(39L, 3L, "+"), L, ref, strand = "+")
  expect_equal(substr(p_plus, L - 2L, L), ref) # block at 3' end, plus strand
  p_minus <- make_as_primer(W, as_anchor(39L, 3L, "-"), L, alt, strand = "-")
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(alt)))
  expect_equal(substr(p_minus, L - 2L, L), rc) # 3' block is the complement
})

test_that("design_as_primers discriminates a 3-nt substitution block", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = "")) # 80 bp
  d <- design_as_primers(W, 39L,
    ref = "GTA", alt = "ACG",
    lengths = 20:22, stab_ats = NA
  )
  expect_gt(nrow(d), 0L)
  expect_gt(max(d$discrimination), 0) # a real block mismatch discriminates
  expect_gt(d$ref_on[1L], d$ref_off[1L])
})

test_that("pick_target_snp accepts an MNP and keeps REF/ALT strings", {
  seq <- "TTTTTTTTTTTTTTTTTTGTACCCCCCCCCCCCCCCCCCCCCCCCCC"
  genome <- Biostrings::DNAStringSet(c(Chr1 = seq))
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L)
  variants <- tibble::tibble(
    id = "mnp1", chrom = "Chr1", pos = 19L,
    ref_genome = list(c("GTA", "ACG")), consequence = "missense"
  )
  tgt <- pick_target_snp(variants, genome, gene)
  expect_equal(tgt$ref, "GTA")
  expect_equal(tgt$alt, "ACG")
  expect_equal(nchar(tgt$ref), nchar(tgt$alt))
})

test_that("pick_target_snp skips indels and ambiguous alleles", {
  withr::local_envvar(LANGUAGE = "en")
  seq <- "ACGTACGTACGTACGTACGTACGTACGT"
  genome <- Biostrings::DNAStringSet(c(Chr1 = seq))
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L)
  variants <- tibble::tibble(
    id = c("del", "amb"),
    chrom = "Chr1", pos = c(5L, 9L),
    ref_genome = list(c("AC", "A"), c("A", "N")),
    consequence = "missense"
  )
  expect_error(pick_target_snp(variants, genome, gene), "No usable substitution")
})

# ---- the multi (degenerate) DAG and its ALT-set consensus gate ----------------
test_that("the multi DAG obeys the contract and shares the KASP core stages", {
  expect_dag_invariants("multi")
  s <- kasp_decision_schema("multi")
  expect_identical(
    s$nodes$name[s$nodes$kind == "gate"],
    c(paste0("as:gate:", KASP_GATES_AS), paste0("pair:gate:", KASP_GATES_PAIR))
  )
  hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
  expect_setequal(hard, c("as:gate:gate_discrimination", "as:gate:gate_as_window"))
})

test_that("the wobble is the worst case over ALL ALTs, not just the first", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(bseq(600L, 2L))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- as.character(Biostrings::subseq(W, 200L, 200L))
  # the transition is the wobble ALT (A>G / C>T / G>A / T>C); the transversion is not
  wob <- c(A = "G", C = "T", G = "A", T = "C")[[ref]]
  non <- setdiff(c("A", "C", "G", "T"), c(ref, wob))[1L]
  snp <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = ref, alt = non, alt2 = wob)
  a <- suppressWarnings(assemble_kasp_assay(genome, snp, flank = 150L, cores = 1L))
  expect_equal(a$site$class[1], "degenerate")
  # the FIRST ALT is not a wobble, but the ALT set contains one, so the primer's
  # wobble flag must be TRUE (the old code looked at alt[1] and missed it)
  expect_false(wobble_mismatch(ref, non))
  expect_true(all(a$as_candidates$wobble))
})

test_that("the degenerate ALT gate verifies the 3' code is the consensus of ALL ALTs", {
  snp <- tibble::tibble(chrom = "Chr1", pos = 1L, ref = "T", alt = "A")
  # not the consensus: the primer's 3' base is a plain ALT base, not the code R
  concrete <- tibble::tibble(alt_degen = "ACGTACGTACGTACGTACGT")
  expect_true(inherits(kasp_degenerate_alt_gate(concrete, c("A", "G"), snp, list()), "kasp_doomed"))
  # the true consensus (3' base R over A/G, REF T excluded) passes
  consensus <- tibble::tibble(alt_degen = c("ACGTACGTACGTACGTACGTR", "acgt"))
  expect_false(inherits(kasp_degenerate_alt_gate(consensus, c("A", "G"), snp, list()), "kasp_doomed"))
  # a code whose members ALSO include REF could extend the reference template, so
  # it discriminates nothing and is doomed even though it is a valid consensus
  ref_amb <- tibble::tibble(alt_degen = "ACGTACGTACGTACGTACGTK")
  expect_true(inherits(kasp_degenerate_alt_gate(ref_amb, c("G", "T"), snp, list()), "kasp_doomed"))
  # a biallelic site has one ALT and nothing to cover
  expect_identical(kasp_degenerate_alt_gate(concrete, "A", snp, list()), concrete)
})
