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
  expect_s3_class(m, "marker_sanger")
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
  # the read primers are nested, not the PCR primers drawn at the amplicon ends
  expect_gt(m$seq_primers$forward_start, m$amplicon$start)
  expect_lt(m$seq_primers$forward_end, snp$pos)
  expect_gt(m$seq_primers$reverse_start, snp$pos)
  expect_lt(m$seq_primers$reverse_end, m$amplicon$end)
  # both reads reach the variant inside the usable window [min_dist, read_length]
  expect_true(all(m$seq_primers$usable))
  expect_gte(snp$pos - m$seq_primers$forward_start, KASP_SANGER_LIMITS$min_dist)
  expect_lte(snp$pos - m$seq_primers$forward_start, KASP_SANGER_LIMITS$read_length)
  expect_gte(m$seq_primers$reverse_end - snp$pos, KASP_SANGER_LIMITS$min_dist)
  expect_lte(m$seq_primers$reverse_end - snp$pos, KASP_SANGER_LIMITS$read_length)
  spec <- marker_map_spec(m, "sanger")
  expect_setequal(spec$feature, c("amplicon", "PCR F", "PCR R", "read F", "read R"))
  expect_equal(spec$start[spec$feature == "read F"], m$seq_primers$forward_start)
  expect_equal(spec$end[spec$feature == "read R"], m$seq_primers$reverse_end)
  oligos <- marker_primers(m)
  expect_true(all(c("PCR_F", "PCR_R", "read_F", "read_R") %in% names(oligos)))
})

test_that("print.marker_sanger summarises the marker and its sequencing primers", {
  m <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 10L, ref = "A", alt = "G"),
      primers = list(forward = "ACGTACGTACGTACGTAC", reverse = "TGCATGCATGCATGCATG", product = 600L),
      amplicon = tibble::tibble(start = 5L, end = 605L),
      seq_primers = list(forward = "AAAACCCCGGGGTTTT", reverse = "TTTTGGGGCCCCAAAA", nested = TRUE),
      selection = list(relaxed = character(0), violated = character(0))
    ),
    class = "marker_sanger"
  )
  expect_output(print(m), "marker_sanger")
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
  expect_s3_class(m, "marker_sanger")
  expect_true("gate_sanger_amplicon" %in% m$selection$relaxed)
  expect_gte(m$primers$product, 200L) # the relaxed floor: 400 - 200
})

test_that("a variant beyond the usable read is warned about, not silently paired", {
  withr::local_envvar(LANGUAGE = "en")
  W <- bseq(1200L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 600L, 600L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 600L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  # a 400+ bp amplicon cannot put the variant within a 100 bp read of both ends,
  # so the design must say which side cannot call it
  expect_warning(
    m <- design_sanger_marker(genome, snp, flank = 500L, read_length = 100L),
    "No usable Sanger read"
  )
  expect_s3_class(m, "marker_sanger")
  expect_true(any(!m$seq_primers$usable))
  d <- c(
    forward = snp$pos - m$seq_primers$forward_start,
    reverse = m$seq_primers$reverse_end - snp$pos
  )
  expect_true(all(d[!m$seq_primers$usable] > 100L))
  expect_output(print(m), "unusable read")
})

test_that("a read primer falls back to the flanking primer when no nested one qualifies", {
  W <- bseq(900L, 11L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  start_f <- 110L
  end_r <- 790L
  flank_f <- substr(W, start_f, start_f + 19L)
  flank_r <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(substr(W, end_r - 19L, end_r))
  ))
  # an impossible Tm window rejects every nested candidate, so the routine must
  # fall back to the PCR pair (inside the 700 bp reach) instead of returning NULL
  sp <- sanger_seq_primers(
    Biostrings::DNAString(W),
    si = 450L, start_f = start_f, end_r = end_r,
    flank_f = flank_f, flank_r = flank_r, min_dist = 50L, genome = genome,
    lengths = 18:25, tm_range = c(0, 1), gc_range = c(0, 1)
  )
  expect_false(sp$nested)
  expect_identical(sp$forward, flank_f)
  expect_identical(sp$reverse, flank_r)
  expect_equal(sp$forward_start, start_f)
  expect_equal(sp$reverse_end, end_r)
  expect_true(all(sp$usable))
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
