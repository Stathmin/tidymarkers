# Sanger as a marker scheme: one flanking pair around a well-placed variant,
# plus the sequencing primers that read into it.

test_that("design_sanger_marker returns a unique amplicon with the variant inside", {
  withr::local_envvar(LANGUAGE = "en")
  W <- bseq(1200L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 600L, 600L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 600L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  # no population variants: the primer-over-variant check did not run, which must
  # be said out loud rather than passed silently
  expect_warning(
    m <- design_sanger_marker(genome, snp, flank = 500L),
    "allele dropout is not checked"
  )
  expect_identical(m$allele_dropout_risk, NA)
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
  m <- suppressWarnings(design_sanger_marker(genome, snp, flank = 140L))
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
    m <- withCallingHandlers(
      design_sanger_marker(genome, snp, flank = 500L, read_length = 100L),
      warning = function(w) {
        if (grepl("allele dropout", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
      }
    ),
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

# ---------------------------------------------------------------------------
# The read-out half: deconvolving an observed 4-channel trace
# ---------------------------------------------------------------------------

# An unambiguous primary read (no repeats, so a shift cannot align by accident).
trace_primary <- "ACGATCGTAGCTAGCATCGATTGAC"
# A minor allele that differs from the primary at *every* base (the complement),
# so a heterozygote can never be mistaken for a homozygote.
het_minor <- chartr("ACGT", "TGCA", trace_primary)

test_that("the detection floor is one uncalibrated constant", {
  expect_identical(SANGER_DETECTION$floor, 0.15)
  expect_identical(SANGER_DETECTION$mixed, 0.20)
  expect_false(SANGER_DETECTION$calibrated)
  expect_gte(SANGER_DETECTION$mixed, SANGER_DETECTION$floor)
  d <- sanger_deconvolve(sanger_example_trace(trace_primary))
  expect_false(any(d$calibrated))
})

test_that("a 50/50 heterozygote is called mixed at every position", {
  secondary <- het_minor
  tr <- sanger_example_trace(trace_primary, secondary = secondary, fraction = 0.5)
  d <- sanger_deconvolve(tr)
  expect_true(all(d$mixed))
  expect_true(all(d$status == "mixed"))
  expect_equal(unique(d$alt_fraction), 0.5)
  # the primary and minor calls are the two bases actually in the trace
  expect_equal(d$ref, strsplit(trace_primary, "")[[1L]])
  expect_equal(d$alt, strsplit(secondary, "")[[1L]])
  expect_equal(unique(d$call_quality), 0)
  f <- sanger_allele_fraction(tr, 1L)
  expect_equal(f$alt_fraction, 0.5)
  expect_true(f$above_floor)
})

test_that("a 20/80 minor allele is mixed, an 18/82 one is only unknown", {
  secondary <- het_minor
  tr <- sanger_example_trace(trace_primary, secondary = secondary, fraction = 0.2)
  d <- sanger_deconvolve(tr)
  expect_true(all(d$status == "mixed"))
  expect_true(all(abs(d$alt_fraction - 0.2) < 1e-3))
  expect_true(sanger_allele_fraction(tr, 1L)$above_floor)
  # just under the heterozygote threshold: a real minor peak, but not a call
  tr2 <- sanger_example_trace(trace_primary, secondary = secondary, fraction = 0.18)
  d2 <- sanger_deconvolve(tr2)
  expect_false(any(d2$mixed))
  expect_true(all(d2$status == "unknown"))
  expect_true(all(abs(d2$alt_fraction - 0.18) < 1e-3))
  expect_false(sanger_allele_fraction(tr2, 1L)$above_floor)
})

test_that("a 5/95 trace is below the floor and is never called a confident homozygote", {
  secondary <- het_minor
  tr <- sanger_example_trace(trace_primary, secondary = secondary, fraction = 0.05)
  d <- sanger_deconvolve(tr)
  expect_false(any(d$mixed))
  # the minor peak is reported (its fraction is not hidden), but the call is
  # `noise`: below the floor a minor allele is not callable from the trace
  expect_true(all(abs(d$alt_fraction - 0.05) < 1e-3))
  expect_true(all(d$status == "noise"))
  f <- sanger_allele_fraction(tr, 1L)
  expect_false(f$above_floor)
  expect_equal(f$status, "noise")
  expect_equal(f$alt_fraction, 0.05)
})

test_that("a homozygote has no minor allele, a tie-free single peak", {
  tr <- sanger_example_trace(trace_primary)
  d <- sanger_deconvolve(tr)
  expect_true(all(d$status == "homozygous"))
  expect_false(any(d$mixed))
  expect_equal(d$alt, rep("", nrow(d)))
  expect_true(all(d$alt_fraction == 0))
  expect_equal(d$ref, strsplit(trace_primary, "")[[1L]])
  expect_true(all(d$call_quality == 1))
  # no mixed stretch either, so an indel cannot be claimed
  expect_message(i <- sanger_indel_estimate(tr), "no visible heterozygous indel")
  expect_true(is.na(i$indel_size))
  expect_true(is.na(i$pos))
})

test_that("an unknown position is reported, not guessed", {
  withr::local_envvar(LANGUAGE = "en")
  expect_message(expect_null(sanger_allele_fraction(sanger_example_trace(trace_primary), 999L)), "not in")
})

test_that("a matrix trace and a pos column are accepted; rows stay aligned", {
  tr <- sanger_example_trace(trace_primary, secondary = het_minor, fraction = 0.5)
  m <- as.matrix(tr[, c("pos", "A", "C", "G", "T")])
  d <- sanger_deconvolve(m)
  expect_s3_class(d, "tbl_df")
  expect_equal(d$pos, tr$pos)
  expect_true(all(d$mixed))
  d2 <- sanger_deconvolve(tr[, c("A", "C", "G", "T")])
  expect_equal(d2$pos, seq_len(nrow(tr)))
  expect_error(sanger_deconvolve(tr[, c("A", "C")]), "missing channel")
})

test_that("bleed-through is corrected by the supplied matrix", {
  # a clean 50/50 trace: A+T against C+G, with a constant instrument floor so
  # each channel's baseline is that floor and the minor peak survives it
  base <- data.frame(
    pos = 1:4,
    A = c(500, 10, 10, 500),
    C = c(10, 500, 500, 10),
    G = c(10, 500, 500, 10),
    T = c(500, 10, 10, 500),
    base = c("A", "C", "C", "A")
  )
  # 10% of every dye leaks into the next channel
  leak <- diag(4)
  leak[cbind(2:4, 1:3)] <- 0.1
  leaked <- base
  leaked[, c("A", "C", "G", "T")] <- as.matrix(base[, c("A", "C", "G", "T")]) %*% t(leak)
  uncorrected <- sanger_deconvolve(leaked, bleed = NULL)
  corrected <- sanger_deconvolve(leaked, bleed = solve(t(leak)))
  # the leak shifts the ratio away from 0.5; the correction brings it back
  expect_gt(max(abs(uncorrected$alt_fraction - 0.5)), 0.01)
  expect_lt(
    max(abs(corrected$alt_fraction - 0.5)),
    max(abs(uncorrected$alt_fraction - 0.5))
  )
  expect_true(all(corrected$status == "mixed"))
})

test_that("a heterozygous deletion is located and sized heuristically", {
  # the minor allele is the primary with bases 12-15 deleted: up to position 11
  # the two alleles agree, from 12 on the trace is mixed (the deletion's first
  # *visible* position is the base after the deleted stretch)
  minor <- paste0(substr(trace_primary, 1L, 11L), substr(trace_primary, 16L, nchar(trace_primary)))
  tr <- sanger_example_trace(trace_primary, secondary = minor, fraction = 0.5, shift_at = 12L)
  i <- sanger_indel_estimate(tr)
  expect_equal(i$pos, 12L)
  expect_equal(i$indel_size, -4L)
  expect_equal(i$confidence, "heuristic")
  expect_false(i$calibrated)
  expect_equal(i$heterozygous_fraction, 0.5)
  # a clean trace has no site to report
  expect_message(sanger_indel_estimate(sanger_example_trace(trace_primary)), "no visible")
})

test_that("a clean trace yields no indel estimate, and a below-floor mixture is invisible", {
  # 5% minor allele: the mix is under the floor, so the stretch is not mixed
  minor <- paste0(substr(trace_primary, 1L, 11L), substr(trace_primary, 16L, nchar(trace_primary)))
  tr <- sanger_example_trace(trace_primary, secondary = minor, fraction = 0.05, shift_at = 12L)
  expect_message(i <- sanger_indel_estimate(tr), "no visible heterozygous indel")
  expect_true(is.na(i$indel_size))
})

test_that("design_sanger_marker flags a population variant under a primer, not a silent pass", {
  withr::local_envvar(LANGUAGE = "en")
  W <- bseq(1200L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 600L, 600L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 600L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  base <- suppressWarnings(design_sanger_marker(genome, snp, flank = 500L))
  expect_identical(base$allele_dropout_risk, NA)
  # a known variant inside the forward PCR primer, on a haplotype that also
  # carries the target ALT (so the population is not monomorphic at the target)
  under <- tibble::tibble(
    haplotype = c("X", "X"), freq = c(0.2, 0.2), chrom = "Chr1",
    pos = c(base$amplicon$start + 2L, snp$pos), ref = c("A", ref), alt = c("G", snp$alt)
  )
  expect_warning(m <- design_sanger_marker(genome, snp, flank = 500L, population = under), "allele-dropout")
  expect_true(m$allele_dropout_risk)
  expect_equal(nrow(m$allele_dropout), 1L)
  expect_equal(m$allele_dropout$pos, base$amplicon$start + 2L)
  expect_equal(m$allele_dropout$primer, "PCR_F")
  expect_output(print(m), "allele-dropout risk")
  # a variant clear of every primer is not a risk
  clear <- tibble::tibble(
    haplotype = c("Y", "Y"), freq = c(0.3, 0.3), chrom = "Chr1",
    pos = c(base$amplicon$start + 100L, snp$pos), ref = c("A", ref), alt = c("G", snp$alt)
  )
  m2 <- design_sanger_marker(genome, snp, flank = 500L, population = clear)
  expect_false(m2$allele_dropout_risk)
  expect_equal(nrow(m2$allele_dropout), 0L)
  expect_output(print(m2), "marker_sanger")
})
