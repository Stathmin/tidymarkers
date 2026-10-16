# InDel size-polymorphism genotyping: one flanking pair, a resolvable size gap.

test_that("an insertion gives a REF/ALT product pair differing by the indel length", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  # a VCF-anchored 10 bp insertion: ref = anchor, alt = anchor + 10 bases
  ins <- paste0(ref, "TTTTTTTTTT")
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = ref, alt = ins)
  m <- design_indel_size_marker(genome, snp, flank = 250L, min_gap = 5L)
  expect_s3_class(m, "marker_indel_size")
  expect_equal(m$gap, 10L)
  # the pair straddles the indel, and the two alleles differ by the indel length
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  expect_equal(m$amplicon$alt_size - m$amplicon$ref_size, 10L)
  # both products sit in the accepted window (applied to both alleles)
  expect_gte(m$amplicon$ref_size, 80L)
  expect_lte(m$amplicon$alt_size, 300L)
  expect_true(
    grepl("^[ACGT]+$", m$primers$forward) && grepl("^[ACGT]+$", m$primers$reverse)
  )
  # the pair must amplify exactly one locus (hard uniqueness gate)
  expect_equal(nrow(in_silico_pcr(m$primers$forward, m$primers$reverse, genome)), 1L)
  expect_output(print(m), "marker_indel_size")
})

test_that("a deletion gives an ALT product shorter by the deleted length", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 3L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 312L) # a 13 bp REF block
  alt <- substr(ref, 1L, 4L) # delete 9 bp (VCF-anchored)
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = ref, alt = alt)
  m <- design_indel_size_marker(genome, snp, flank = 250L)
  expect_s3_class(m, "marker_indel_size")
  expect_equal(m$gap, 9L)
  expect_equal(m$amplicon$ref_size - m$amplicon$alt_size, 9L)
  expect_lt(m$amplicon$alt_size, m$amplicon$ref_size)
})

test_that("an equal-length (substitution) target is declined with a message", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(400L, 4L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = "A", alt = "G")
  expect_message(design_indel_size_marker(genome, snp), "same length")
})

test_that("a sub-resolvable size difference is declined with a message", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(400L, 5L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = "A", alt = "AT") # 1 bp
  # below even the 3 bp capillary floor, so no ladder can reach it
  expect_message(design_indel_size_marker(genome, snp, min_gap = 5L), "min_gap")
})

test_that("the resolvable-gap gate is evaluated and relaxes to the capillary floor (C7)", {
  skip_if_no_oligoarrayaux()
  W <- bseq(300L, 8L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 150L, 150L)
  # a 3 bp insertion: below the 5 bp agarose bar but at the capillary floor, so
  # the ladder must lower gate_indel_size_gap and record it, not reject
  snp <- tibble::tibble(chrom = "Chr1", pos = 150L, ref = ref, alt = paste0(ref, "TTT"))
  m <- design_indel_size_marker(genome, snp, flank = 140L, min_gap = 5L)
  expect_s3_class(m, "marker_indel_size")
  expect_equal(m$gap, 3L)
  expect_true("gate_indel_size_gap" %in% m$selection$relaxed)
})

test_that("scheme_def_indel_size declares the InDel-size gates, tiers and labels", {
  d <- scheme_def_indel_size()
  expect_equal(d$scheme, "indel_size")
  expect_equal(d$family, "indel_size")
  expect_equal(d$loci, 1L)
  expect_setequal(
    d$gates,
    c(
      "gate_indel_size_gap", "gate_indel_size_amplicon", "gate_indel_size_tm",
      "gate_indel_size_structure", "gate_indel_size_unique"
    )
  )
  expect_setequal(names(d$gate_short), d$gates)
  expect_setequal(names(d$gate_limit), d$gates)
  expect_identical(names(d$tiers), c("tier1", "tier2", "tier3"))
  expect_true(all(d$relax %in% d$gates))
  # the gap gate is soft (the lattice may relax it to the capillary floor); the
  # structure and uniqueness gates are hard and never relaxed
  expect_true("gate_indel_size_gap" %in% d$relax)
  expect_setequal(
    names(d$gate_hard_label),
    c("gate_indel_size_structure", "gate_indel_size_unique")
  )
  expect_setequal(names(d$gate_floor_label), setdiff(d$gates, names(d$gate_hard_label)))
})

test_that("a duplicated locus dooms the design through the hard uniqueness gate (C7)", {
  skip_if_no_oligoarrayaux()
  # two tandem copies: every primer's 3' seed occurs twice, so no candidate can
  # pass the hard uniqueness gate and the design must give up
  W <- bseq(400L, 4L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = paste0(W, W)))
  ref <- substr(W, 200L, 200L)
  snp <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = ref, alt = paste0(ref, "TTTTTTTTTT"))
  expect_message(
    expect_null(design_indel_size_marker(genome, snp, flank = 150L)),
    "unique 3' ends"
  )
})

test_that("indel_size DAG obeys the contract", {
  expect_dag_invariants("indel_size")
})
