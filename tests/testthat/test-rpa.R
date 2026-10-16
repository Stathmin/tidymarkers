# RPA (isothermal, recombinase-driven): two long primers, a short amplicon and an
# optional exo probe.  Design rules are grounded in the file header of marker_rpa.R.

test_that("design_rpa_marker returns long primers and a short amplicon", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = ref, alt = alt)
  m <- design_rpa_marker(genome, snp, flank = 200L)
  expect_s3_class(m, "marker_rpa")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))

  # RPA needs a long recombinase filament: both primers are in the length window
  expect_true(grepl("^[ACGT]+$", m$primers$forward) && grepl("^[ACGT]+$", m$primers$reverse))
  expect_true(dplyr::between(nchar(m$primers$forward), 30L, 35L))
  expect_true(dplyr::between(nchar(m$primers$reverse), 30L, 35L))

  # the amplicon is short and the variant sits under one of the primers
  expect_gte(m$primers$product, 100L)
  expect_lte(m$primers$product, 200L)
  expect_true(m$primers$variant_under %in% c("forward", "reverse"))
  expect_true(m$amplicon$start < 300L && m$amplicon$end > 300L)
  expect_null(m$probe)
  expect_output(print(m), "marker_rpa")
})

test_that("design_rpa_marker can add an exo probe that covers the variant", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = ref, alt = alt)
  m <- design_rpa_marker(genome, snp, flank = 200L, design_probe = TRUE)
  expect_s3_class(m, "marker_rpa")
  # the primers now flank the variant; the probe carries the discrimination
  expect_identical(m$primers$variant_under, "probe")
  expect_true(m$amplicon$start < 300L && m$amplicon$end > 300L)
  # an exo probe: 46-52 nt, a THF ("X") inside the published 5'/3' arm windows,
  # and a dT-fluorophore / dT-quencher within the allowed gap
  expect_true(nzchar(m$probe$sequence) && grepl("^[ACGTX]+$", m$probe$sequence))
  expect_true(grepl("X", m$probe$sequence, fixed = TRUE))
  expect_true(m$probe$start <= 300L && m$probe$end >= 300L)
  expect_true(m$probe$start >= m$amplicon$start && m$probe$end <= m$amplicon$end)
  expect_true(dplyr::between(m$probe$len, 46L, 52L))
  expect_gte(m$probe$thf_arm5, 30L)
  expect_gte(m$probe$thf_arm3, 15L)
  expect_lte(m$probe$fq_gap5, 4L)
  expect_lte(m$probe$fq_gap3, 4L)
  expect_output(print(m), "exo probe")
})

test_that("design_rpa_marker declines non-substitution alleles", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "A", alt = "AA")
  expect_message(
    suppressWarnings(design_rpa_marker(genome, snp, flank = 15L)),
    "single-base"
  )
})

test_that("scheme_def_rpa describes the chemistry's gates, tiers and stage", {
  d <- scheme_def_rpa()
  expect_equal(d$scheme, "rpa")
  expect_equal(d$family, "rpa")
  expect_setequal(d$gates, KASP_GATES_RPA)
  expect_true(all(d$relax %in% d$gates))
  expect_setequal(names(d$gate_short), d$gates)
  expect_setequal(names(d$gate_limit), d$gates)
  expect_true(all(names(d$gate_floor_label) %in% d$relax))
  expect_true(all(names(d$gate_hard_label) %in% setdiff(d$gates, d$relax)))
  expect_setequal(names(d$tiers), c("tier1", "tier2", "tier3", "tier4"))
  expect_true(nzchar(d$stage_label))
  expect_equal(d$loci, 1L)
})

test_that("the rpa DAG obeys the contract", expect_dag_invariants("rpa"))
