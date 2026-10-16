# Sanger as a marker scheme: one flanking pair around a well-placed variant,
# plus the sequencing primers that read into it.

test_that("design_sanger_marker returns a unique amplicon with the variant inside", {
  W <- bseq(1200L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 600L, 600L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 600L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_sanger_marker(genome, snp, flank = 500L)
  expect_s3_class(m, "kasp_sanger")
  # the variant is inside the amplicon and clear of both read starts
  expect_gt(m$amplicon$start, 0L)
  expect_lt(m$amplicon$start, snp$pos)
  expect_gt(m$amplicon$end, snp$pos)
  expect_gte(snp$pos - m$amplicon$start, 50L)
  expect_gte(m$amplicon$end - snp$pos, 50L)
  # the amplicon sits in the Sanger window and the primers are plain ACGT
  expect_gte(m$primers$product, 400L)
  expect_lte(m$primers$product, 800L)
  expect_true(grepl("^[ACGT]+$", m$primers$forward) && grepl("^[ACGT]+$", m$primers$reverse))
  expect_equal(nrow(in_silico_pcr(m$primers$forward, m$primers$reverse, genome)), 1L)
  # a sequencing primer is designed on each strand, inside the amplicon
  expect_true(m$seq_primers$nested)
  expect_true(grepl(m$seq_primers$forward, W, fixed = TRUE))
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(m$seq_primers$reverse)))
  expect_true(grepl(rc, W, fixed = TRUE))
})

test_that("print.kasp_sanger summarises the marker and its sequencing primers", {
  m <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 10L, ref = "A", alt = "G"),
      primers = list(forward = "ACGTACGTACGTACGTAC", reverse = "TGCATGCATGCATGCATG", product = 600L),
      amplicon = tibble::tibble(start = 5L, end = 605L),
      seq_primers = list(forward = "AAAACCCCGGGGTTTT", reverse = "TTTTGGGGCCCCAAAA", nested = TRUE),
      selection = list(relaxed = character(0), violated = character(0))
    ),
    class = "kasp_sanger"
  )
  expect_output(print(m), "kasp_sanger")
  expect_output(print(m), "600 bp")
  expect_output(print(m), "sequencing primers")
})

test_that("design_sanger_marker declines when no pair keeps the variant off the read start", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(300L, 3L)))
  ref <- substr(as.character(genome$Chr1), 25L, 25L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 25L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  # the hard `dist` gate dooms every candidate: no pair leaves >= 50 bp on both
  # sides of the variant in a 300 bp template
  expect_message(design_sanger_marker(genome, snp, flank = 20L), "read starts")
})

test_that("the amplicon gate is really evaluated and relaxes to its floor (C7)", {
  withr::local_envvar(LANGUAGE = "en")
  W <- bseq(300L, 9L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 150L, 150L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 150L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  # a 300 bp template cannot hold the ideal 400-800 bp amplicon, so the ladder
  # must lower gate_sanger_amplicon to its floor and say so
  m <- design_sanger_marker(genome, snp, flank = 140L)
  expect_s3_class(m, "kasp_sanger")
  expect_true("gate_sanger_amplicon" %in% m$selection$relaxed)
  expect_gte(m$primers$product, 200L) # the relaxed floor: 400 - 200
})

test_that("scheme_def_sanger describes the chemistry and its gates", {
  d <- scheme_def_sanger()
  expect_identical(d$scheme, "sanger")
  expect_identical(d$family, "sanger")
  expect_identical(d$loci, 1L)
  expect_setequal(d$gates, names(d$gate_short))
  expect_setequal(d$gates, names(d$gate_limit))
  expect_true(all(d$relax %in% d$gates))
  expect_true(all(d$gates %in% c(names(d$gate_floor_label), names(d$gate_hard_label))))
  expect_length(d$tiers, 3L)
  # exactly the placement and uniqueness gates are hard
  expect_setequal(
    names(d$gate_hard_label),
    c("gate_sanger_dist", "gate_sanger_unique")
  )
  # gate order is the declared order (grave / cheap first)
  expect_identical(d$gates[1L], "gate_sanger_dist")
  expect_identical(d$gates[length(d$gates)], "gate_sanger_unique")
})

test_that("a duplicated locus dooms the design through the hard uniqueness gate (C7)", {
  # two tandem copies of the template: every primer's 3' seed occurs twice, so no
  # candidate can pass the hard `unique` gate and the design must give up
  W <- bseq(400L, 4L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = paste0(W, W)))
  ref <- substr(W, 200L, 200L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 200L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  expect_message(
    expect_null(design_sanger_marker(genome, snp, flank = 150L)),
    "unique 3' ends"
  )
})

test_that("sanger DAG obeys the contract", {
  expect_dag_invariants("sanger")
})
