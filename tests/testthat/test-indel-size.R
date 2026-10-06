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
  expect_s3_class(m, "kasp_indel_size")
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
  expect_output(print(m), "kasp_indel_size")
})

test_that("a deletion gives an ALT product shorter by the deleted length", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 3L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 312L) # a 13 bp REF block
  alt <- substr(ref, 1L, 4L) # delete 9 bp (VCF-anchored)
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = ref, alt = alt)
  m <- design_indel_size_marker(genome, snp, flank = 250L)
  expect_s3_class(m, "kasp_indel_size")
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
  expect_message(design_indel_size_marker(genome, snp, min_gap = 5L), "min_gap")
})

test_that("scheme_def_indel_size declares the InDel-size gates, tiers and labels", {
  d <- scheme_def_indel_size()
  expect_equal(d$scheme, "indel_size")
  expect_equal(d$family, "indel_size")
  expect_equal(d$loci, 1L)
  expect_true(all(c(
    "gate_indel_size_gap", "gate_indel_size_amplicon", "gate_indel_size_tm",
    "gate_indel_size_structure", "gate_indel_size_unique"
  ) %in% d$gates))
  expect_setequal(names(d$gate_short), d$gates)
  expect_setequal(names(d$gate_limit), d$gates)
  expect_identical(names(d$tiers), c("tier1", "tier2", "tier3"))
  expect_true(all(d$relax %in% d$gates))
  # the hard-gate labels and the relaxation floors cover exactly the gates
  expect_setequal(names(d$gate_hard_label), setdiff(d$gates, d$relax))
  expect_true(all(names(d$gate_floor_label) %in% d$relax))
})
