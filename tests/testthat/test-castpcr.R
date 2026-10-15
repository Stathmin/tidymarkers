# castPCR / allele-specific qPCR: one AS primer + a locus probe + a common primer.

test_that("design_castpcr_marker returns an AS primer, a probe and a short amplicon", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = ref, alt = alt)
  m <- design_castpcr_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "kasp_castpcr")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))

  # one AS primer whose 3' base sits on the variant base
  as <- m$primers$as_primer
  expect_true(grepl("^[ACGT]+$", as))
  last <- substr(as, nchar(as), nchar(as))
  complement_ref <- as.character(Biostrings::complement(Biostrings::DNAString(ref)))
  on_plus <- if (identical(m$primers$as_role, "forward")) ref else complement_ref
  expect_equal(last, on_plus)

  # a locus probe that covers the variant ...
  expect_lte(m$probe$start, 300L)
  expect_gte(m$probe$end, 300L)
  # ... and a short amplicon that contains both primers and the probe
  expect_gte(m$primers$product, 60L)
  expect_lte(m$primers$product, 150L)
  expect_true(m$probe$start >= m$amplicon$start && m$probe$end <= m$amplicon$end)
  expect_true(grepl("^[ACGT]+$", m$primers$forward) && grepl("^[ACGT]+$", m$primers$reverse))
  expect_output(print(m), "kasp_castpcr")
})

test_that("design_castpcr_marker declines non-substitution alleles", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "A", alt = "AA")
  expect_message(
    suppressWarnings(design_castpcr_marker(genome, snp, flank = 15L)),
    "single-base"
  )
})

test_that("scheme_def_castpcr describes the chemistry's gates, tiers and stage", {
  d <- scheme_def_castpcr()
  expect_equal(d$scheme, "castpcr")
  expect_equal(d$family, "castpcr")
  expect_setequal(d$gates, KASP_GATES_CASTPCR)
  expect_true(all(d$relax %in% d$gates))
  expect_setequal(names(d$gate_short), d$gates)
  expect_setequal(names(d$gate_limit), d$gates)
  expect_true(all(names(d$gate_floor_label) %in% d$relax))
  expect_true(all(names(d$gate_hard_label) %in% setdiff(d$gates, d$relax)))
  expect_setequal(names(d$tiers), c("tier1", "tier2", "tier3"))
  expect_true(nzchar(d$stage_label))
  expect_equal(d$loci, 1L)
})
