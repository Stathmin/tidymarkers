# dPCR as a marker scheme: a hydrolysis probe over a short amplicon.

test_that("design_dpcr_marker returns a short amplicon with a hydrolysis probe over the variant", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 1L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_dpcr_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "kasp_dpcr")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  # a short amplicon with the probe between the primers and over the variant
  expect_gte(m$primers$product, 60L)
  expect_lte(m$primers$product, 150L)
  expect_gt(m$probe$start, m$amplicon$start)
  expect_lt(m$probe$end, m$amplicon$end)
  expect_lte(m$probe$start, 300L)
  expect_gte(m$probe$end, 300L)
  # the probe's composition rules (18-25 nt, GC 40-60%, no 5' G, no run > 4)
  plen <- nchar(m$probe$sequence)
  expect_gte(plen, 18L)
  expect_lte(plen, 25L)
  expect_match(m$probe$sequence, "^[ACGT]+$")
  expect_false(substr(m$probe$sequence, 1L, 1L) == "G")
  expect_true(gc_frac(m$probe$sequence) >= 0.40 && gc_frac(m$probe$sequence) <= 0.60)
  expect_false(has_run(m$probe$sequence, n = 5L))
  # the probe out-competes the primers (the ~8-10 C rule)
  expect_gte(m$probe$tm - max(m$primers$tm_f, m$primers$tm_r), 8)
  expect_match(m$primers$forward, "^[ACGT]+$")
  expect_match(m$primers$reverse, "^[ACGT]+$")
  expect_output(print.kasp_dpcr(m), "kasp_dpcr")
})

test_that("design_dpcr_marker informs and returns NULL when the flank is too small", {
  W <- bseq(400L, 5L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 200L, 200L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 200L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  expect_message(m <- design_dpcr_marker(genome, snp, flank = 20L), "flanking primer pair")
  expect_null(m)
})

test_that("the dPCR scheme definition carries its gates, relax and tiers", {
  d <- scheme_def_dpcr()
  expect_identical(d$scheme, "dpcr")
  expect_identical(d$family, "dpcr")
  expect_identical(d$loci, 1L)
  expect_identical(d$verdict, "usable dPCR assay")
  gates <- c("gate_dpcr_amplicon", "gate_dpcr_primer_tm", "gate_dpcr_probe_tm", "gate_dpcr_probe_comp")
  expect_setequal(d$gates, gates)
  expect_true(all(d$relax %in% d$gates))
  expect_identical(names(d$tiers), c("tier1", "tier2", "tier3"))
  expect_true(all(d$gates %in% names(d$gate_short)))
  expect_true(all(d$gates %in% names(d$gate_limit)))
  expect_true(all(d$relax %in% names(d$gate_floor_label)))
  expect_true(all(setdiff(d$gates, d$relax) %in% names(d$gate_hard_label)))
  expect_identical(d$stage_label, "dPCR (droplet digital PCR)")
})
