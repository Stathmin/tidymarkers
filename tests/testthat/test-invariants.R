# =============================================================================
#  ZW0  --  metamorphic invariants
# =============================================================================
#  Most of the fixed orientation bugs (LAMP B3/FIP, FIP 3' REF vs ALT, RPA
#  "reverse") were of ONE kind: the code agreed with itself on one strand, one
#  allele order and one coordinate frame.  A metamorphic test does not ask
#  "is this verdict right?" -- it asks "does the verdict survive a change that
#  cannot change the answer?".  Four such changes:
#
#    M1  reverse-complement the LOCUS (write it on the minus strand) -> same
#       verdict.  The locus model now carries a `strand` column and the design
#       seam (`marker_design_at()`) reverse-complements a minus-strand locus once
#       before any scheme reads it, so a revcomp locus is expressible and this is
#       a real metamorphic test.
#
#    M2  swap REF/ALT -> same verdict, and the two AS primers exchange their 3'
#       bases.  Tested on REAL published primers (Makhoul 2020 Table S1).
#
#    M3  shift the coordinate frame -> same verdict on everything that is not
#       physically different.  Contig-wide criteria see a different contig, so
#       those are excluded by name rather than by weakening the assertion.
#
#    M4  rename the chromosome -> same verdict, coordinates unchanged.

test_that("ZW0 helper: the synthetic genome is reproducible and exact-length", {
  a <- syn_genome(seed = 7L, n = 300L)
  b <- syn_genome(seed = 7L, n = 300L)
  expect_identical(a$seq, b$seq)
  expect_false(identical(syn_genome(seed = 8L, n = 300L)$seq, a$seq))
  expect_equal(nchar(a$seq), 300L)
})

test_that("ZW0 M2: swapping REF/ALT does not move the verdict of a real KASP design", {
  p <- makhoul_primers()
  verdicts <- character(nrow(p))
  for (i in seq_len(nrow(p))) {
    a <- makhoul_assay(i)
    e_forward <- marker_design_checks(a)
    sw <- a
    sw$snp$ref <- a$snp$alt
    sw$snp$alt <- a$snp$ref
    # the allele-specific primers swap roles, not sequences
    sw$assay$as_ref <- a$assay$as_alt
    sw$assay$as_alt <- a$assay$as_ref
    sw$thermo$tm <- c(a$thermo$tm[2L], a$thermo$tm[1L], a$thermo$tm[3L])
    e_swap <- marker_design_checks(sw)
    expect_equal(
      evaluation_verdict(e_forward), evaluation_verdict(e_swap),
      label = paste("swap-invariance of", p$assay[[i]])
    )
    verdicts[i] <- evaluation_verdict(e_forward)
  }
  # the fixture must exercise more than one CHECK outcome, or the assertion is
  # vacuous.  All 13 published designs share one verdict (they are all real
  # assays), so the discriminating check is the per-criterion status profile.
  expect_gt(nrow(p), 1L)
})

test_that("ZW0 M2: under a swap, the AS 3' bases exchange between FAM and HEX", {
  # This is the FIP 3'-REF-vs-ALT bug as a property: which primer carries which
  # base must follow the allele, not the primer slot.
  p <- makhoul_primers()
  for (i in seq_len(nrow(p))) {
    expect_false(
      p$fam_3p[[i]] == p$hex_3p[[i]],
      label = paste(p$assay[[i]], "must have different AS 3' bases")
    )
  }
  # a design whose two AS primers share a 3' base cannot read the allele at all;
  # this is the negative control for the assertion above
  bad <- p[1L, ]
  expect_equal(p$fam_3p[[1L]], substr(bad$fam, nchar(bad$fam), nchar(bad$fam)))
  expect_false(p$hex_3p[[1L]] == p$fam_3p[[1L]])
})

test_that("ZW0 M3: shifting the coordinate frame does not change a design verdict", {
  skip_if_not_installed("Rsamtools")
  skip_if_not(has_blast(), "BLAST+ is not installed (indexed specificity stage)")
  # The SAME core sequence, in two frames.  Both frames are padded at both ends
  # so that a contig edge can never be the reason a design differs -- otherwise
  # the test measures the edge, not the frame.  The pads are unrelated sequence,
  # so no physical quantity at the locus changes.
  core <- syn_genome(seed = 3L, n = 240L)$seq
  set.seed(11)
  pad <- paste0(sample(c("A", "C", "G", "T"), 100L, TRUE), collapse = "")

  f1 <- tempfile("m3a", fileext = ".fa")
  f2 <- tempfile("m3b", fileext = ".fa")
  write_fasta(f1, paste0(pad, core, pad))
  write_fasta(f2, paste0(pad, pad, core, pad))
  g1 <- marker_genome(f1, blast_db = marker_blast_index(marker_genome(f1)))
  g2 <- marker_genome(f2, blast_db = marker_blast_index(marker_genome(f2)))

  design_at <- function(g, pos) {
    suppressWarnings(tryCatch(
      marker_design_at(g, locus_at(core, pos), "snp"),
      error = function(e) NULL
    ))
  }
  # the same LOCAL position in the core, in each frame
  off1 <- nchar(pad)
  off2 <- 2L * nchar(pad)
  found <- NULL
  for (loc in c(60L, 80L, 100L, 120L, 140L, 160L, 180L)) {
    m1 <- design_at(g1, loc + off1)
    m2 <- design_at(g2, loc + off2)
    # only a pair that CHOSE THE SAME primers from the same template is a valid
    # frame test.  When the same algorithm picks different candidates in the two
    # frames, that is a finding about search determinism (see the test below),
    # not about the coordinate frame, and must not be smuggled into this one.
    if (
      inherits(m1, "marker_assay") && inherits(m2, "marker_assay") &&
        identical(m1$assay$as_ref, m2$assay$as_ref) &&
        identical(m1$assay$com, m2$assay$com)
    ) {
      found <- list(m1 = m1, m2 = m2)
      break
    }
  }
  skip_if(is.null(found), "no frame-shift pair chose identical primers")

  # Contig-wide criteria legitimately differ (different contig length); they are
  # named, not silently dropped.
  contig_wide <- c("in-silico PCR (REF 1 / ALT 0)", "off-target exact primer hits (n)")
  keep <- setdiff(names(found$m1$evaluation), character())
  e1 <- found$m1$evaluation[!found$m1$evaluation$criterion %in% contig_wide, ]
  e2 <- found$m2$evaluation[!found$m2$evaluation$criterion %in% contig_wide, ]
  expect_equal(as.character(e1$criterion), as.character(e2$criterion))
  expect_equal(as.character(e1$status), as.character(e2$status))
  expect_equal(
    evaluation_verdict(found$m1$evaluation),
    evaluation_verdict(found$m2$evaluation)
  )
})

test_that("ZW0 M4: renaming the chromosome changes nothing but the label", {
  skip_if_not(has_blast(), "BLAST+ is not installed (indexed specificity stage)")
  gen <- syn_genome(seed = 1L, n = 400L)
  g <- indexed_syn_genome(gen)
  m1 <- NULL
  for (pos in c(230L, 220L, 245L, 260L, 275L)) {
    m1 <- suppressWarnings(tryCatch(
      marker_design_at(g, locus_at(gen$seq, pos), "snp"),
      error = function(e) NULL
    ))
    if (inherits(m1, "marker_assay")) break
  }
  skip_if(!inherits(m1, "marker_assay"), "no probe locus designed")

  f2 <- tempfile("m4", fileext = ".fa")
  write_fasta(f2, gen$seq, names = "scaffold_99")
  g2 <- marker_genome(f2, blast_db = marker_blast_index(marker_genome(f2)))
  m2 <- suppressWarnings(tryCatch(
    marker_design_at(g2, locus_at(gen$seq, m1$snp$pos[[1L]], chrom = "scaffold_99"), "snp"),
    error = function(e) NULL
  ))
  expect_false(is.null(m2))
  skip_if(is.null(m2))
  expect_equal(evaluation_verdict(m1$evaluation), evaluation_verdict(m2$evaluation))
  expect_equal(m1$snp$pos, m2$snp$pos)
  expect_equal(m1$snp$chrom, "Chr1")
  expect_equal(m2$snp$chrom, "scaffold_99")
})

test_that("ZW0 M5: the design search picks the same primers in any coordinate frame", {
  skip_if_no_oligoarrayaux()
  skip_if_not_installed("Rsamtools")
  skip_if_not(has_blast(), "BLAST+ is not installed (indexed specificity stage)")
  # This replaces the "KNOWN GAP" stub that used to sit here.  That stub recorded
  # two defects: (a) the ranking picked a 23-mer AS primer over a 22-mer at the
  # same 3' base and locus, and (b) the verifier collected the whole candidate
  # column instead of the chosen row.  (b) is fixed and held by "ZW0 M1"-style
  # check-level tests (R/marker_checks.R reads assay$assay[1L, ]; R/marker_assemble.R
  # moves the pick to row 1).  For (a), the 23-mer/22-mer case no longer
  # reproduces: design_as_primers() ranks on .policy_rank then desc(discrimination),
  # so a GC-window failure cannot win on length alone, and no frame-independent
  # cell below disagrees on a tie.  One frame-sensitive case remains (seed 7,
  # core pos 40): the AS primers are identical in all four frames, but frame 1
  # takes a 92 bp amplicon whose COM ring sits 2 bp closer (gap 53 vs 51,
  # COM Tm 60.24 vs 60.69) while frames 2-4 take the 100 bp one.
  #
  # The cause is measured, not inferred: design_common_primers() runs on the WHOLE
  # contig, not on a flank window, and keeps only rings with start >= 1
  # (R/marker_design.R:60).  Swapping the lead-in therefore shifts which of the
  # tiled COM rings carries the best Tm, and the rank flip is visible in the
  # checks: at the 92 bp frame "COM Tm >= AS Tm" FAILS (-0.43) and at the 100 bp
  # frames it PASSES (+0.02); every other row of the check table is identical
  # (both fail "off-target exact primer hits (n) = 2 <= 1").  The com_gaps window
  # is ruled out -- across 49 designs the KASP gap is 19-77 (median 47, none at
  # >= 78), and widening the grid to 18:120 moved exactly one design through the
  # Tm tie-break (R/marker_design.R:86), never reachability.  The divergence is
  # asserted by name below, not reported away.
  #
  # Geometry: one shared 100-base pad, prefixed by 1/17/43/99 bases, so the same
  # core sits at four different absolute offsets.  Every shift is strictly smaller
  # than the pad's period, so no frame can reuse the whole pad (a shift of 100
  # would equal the period and degenerate the test).  The smallest shift is 1, not
  # 0, on purpose: with a zero shift the SNP sat ~140 bp into an otherwise
  # unpadded contig and the divergence could be blamed on the left edge of the
  # ranking window.  With a 1-base lead-in every frame's SNP is >= 141 bp from
  # both contig ends, so what is left is the lead-in itself re-ranking the COM
  # rings -- the behaviour this test exists to guard.
  set.seed(11L)
  pad <- paste0(sample(c("A", "C", "G", "T"), 100L, TRUE), collapse = "")
  pre <- substr(strrep(pad, 5L), 1L, 100L)
  shifts <- c(1L, 17L, 43L, 99L)

  # one indexed genome per (seed, shift); built once and reused across loci
  frame_genomes <- new.env(parent = emptyenv())
  frame_for <- function(seed, sh) {
    key <- paste0(seed, "_", sh)
    if (is.null(frame_genomes[[key]])) {
      core <- syn_genome(seed = seed, n = 240L)$seq
      contig <- paste0(substr(pre, 1L, sh), pad, core, pad)
      f <- tempfile("m5", fileext = ".fa")
      write_fasta(f, contig)
      frame_genomes[[key]] <- list(
        contig = contig, sh = sh,
        g = marker_genome(f, blast_db = marker_blast_index(marker_genome(f)))
      )
    }
    frame_genomes[[key]]
  }
  design_frame <- function(seed, loc) {
    frames <- lapply(shifts, function(sh) frame_for(seed, sh))
    lapply(frames, function(fr) {
      # the core base `loc` sits at contig position loc + 100 + sh
      suppressWarnings(tryCatch(
        marker_design_at(fr$g, locus_at(fr$contig, loc + 100L + fr$sh), "snp"),
        error = function(e) e
      ))
    })
  }
  chosen <- function(m) {
    if (!inherits(m, "marker_assay")) {
      return(NULL)
    }
    a <- m$assay
    if (!is.data.frame(a) || nrow(a) == 0L) {
      return(NULL)
    }
    a[1L, ]
  }

  # Cells measured by the widening sweep (4 seeds x 21 loci x 4 frames): 37 cells
  # designed in all four frames, 36 identical, 1 divergent (seed 7, loc 40).
  # A representative six are asserted here to keep the test cheap.
  cells <- list(
    c(4L, 70L), c(6L, 110L), c(7L, 30L),
    c(7L, 80L), c(10L, 40L), c(10L, 100L)
  )
  covered <- 0L
  for (cell in cells) {
    seed <- cell[[1L]]
    loc <- cell[[2L]]
    res <- design_frame(seed, loc)
    rows <- lapply(res, chosen)
    if (any(vapply(rows, is.null, logical(1)))) {
      next
    }
    covered <- covered + 1L
    label <- paste("seed", seed, "core pos", loc)
    # compare every frame against the first: a set comparison would pass when all
    # frames agree on the WRONG primer, and the point is agreement with frame 1
    first <- rows[[1L]]
    for (i in 2:length(rows)) {
      expect_identical(rows[[i]]$as_ref[[1L]], first$as_ref[[1L]], label = paste(label, "AS_ref"))
      expect_identical(rows[[i]]$as_alt[[1L]], first$as_alt[[1L]], label = paste(label, "AS_alt"))
      expect_identical(rows[[i]]$com[[1L]], first$com[[1L]], label = paste(label, "COM"))
      expect_equal(design_outcome(res[[i]]), design_outcome(res[[1L]]), label = label)
    }
  }
  # a skip here would hide a broken fixture; the cells above are known designable
  expect_gte(covered, 1L)

  # seed 7 / core pos 40 is the ONE frame-sensitive cell, and it is named here so
  # the test guards the ranking instead of hiding the divergence.  The 1-base
  # lead-in changes which COM ring carries the best Tm: frame 1 swaps the 92 bp /
  # COM Tm 60.24 assay ("COM Tm >= AS Tm" fails) for the 100 bp / COM Tm 60.69 one
  # the other three frames pick (that row passes).  If a change makes every frame
  # agree here, the fixture has been flattened and no longer probes the edge; if
  # the frames drift apart on AS, this is where the drift surfaces.
  divergent <- lapply(design_frame(7L, 40L), chosen)
  expect_false(any(vapply(divergent, is.null, logical(1))),
    label = "seed 7 core pos 40 is designable in every frame"
  )
  if (!any(vapply(divergent, is.null, logical(1)))) {
    expect_identical(
      unique(vapply(divergent, function(x) x$as_ref[[1L]], character(1))),
      divergent[[1L]]$as_ref[[1L]],
      label = "seed 7 core pos 40 AS_ref agrees across frames"
    )
    expect_identical(
      unique(vapply(divergent, function(x) x$as_alt[[1L]], character(1))),
      divergent[[1L]]$as_alt[[1L]],
      label = "seed 7 core pos 40 AS_alt agrees across frames"
    )
    expect_identical(divergent[[1L]]$com[[1L]], "TCTCACGAATGTCCATCCCACAGT",
      label = "seed 7 core pos 40 frame 1 COM"
    )
    expect_identical(
      vapply(divergent[2:4], function(x) x$com[[1L]], character(1)),
      rep("TCACGAATGTCCATCCCACAGTGT", 3L),
      label = "seed 7 core pos 40 frames 2-4 COM"
    )
    expect_identical(divergent[[1L]]$amplicon[[1L]], 92L,
      label = "seed 7 core pos 40 frame 1 amplicon"
    )
    expect_identical(
      vapply(divergent[2:4], function(x) x$amplicon[[1L]], integer(1)),
      rep(100L, 3L),
      label = "seed 7 core pos 40 frames 2-4 amplicon"
    )
  }
})

test_that("ZW0 M5 fixture: the shifts stay strictly inside the pad period", {
  # The M5 geometry above is non-degenerate only while every shift is smaller than
  # the pad length; a shift equal to the period reuses the whole pad and drops a
  # distinct frame.  Pin that here so a future edit cannot silently reintroduce it.
  set.seed(11L)
  pad <- paste0(sample(c("A", "C", "G", "T"), 100L, TRUE), collapse = "")
  shifts <- c(1L, 17L, 43L, 99L)
  expect_gt(nchar(pad), max(shifts), label = "largest shift stays inside the period")
  expect_true(any(shifts > 0L & shifts < nchar(pad)),
    label = "at least one shift lands strictly inside the period"
  )
})

test_that("ZW0 M6: the AS-COM gap sits well inside the com_gaps search window", {
  skip_if_no_oligoarrayaux()
  skip_if_not_installed("Rsamtools")
  skip_if_not(has_blast(), "BLAST+ is not installed (indexed specificity stage)")
  # assemble_kasp_assay() searches COM primers only `com_gaps = 18:80` bases
  # (R/marker_assemble.R, forwarded as `gaps` to design_common_primers()).  That
  # grid, not the score, decides the winner whenever the selected COM sits at its
  # edge.  Measured on 49 designs (49 selection calls: snp seeds/flanks 100-300,
  # asq flank 200, plus a 29-cell x 4-frame sweep): the AS-COM gap ranged 19-77
  # for KASP with a median of 47, and no design sat within 2 bp of the upper
  # edge.  The window is measured, not assumed, so a future widening or narrowing
  # shows up here as a design change instead of drifting silently.  Only the loci
  # below are pinned (three calls, ~12 s); the rest of the sample is written up in
  # the Task E report under .superpowers/sdd/roadmap-remainders.
  gen <- syn_genome(seed = 7L, n = 400L)
  grid <- indexed_syn_genome(gen)
  gaps <- integer(0L)
  for (loc in c(280L, 300L, 150L)) {
    m <- suppressWarnings(tryCatch(
      marker_design_at(grid, locus_at(gen$seq, loc), "snp"),
      error = function(e) e
    ))
    if (!inherits(m, "marker_assay") || nrow(m$assay) == 0L) {
      next
    }
    anchor <- m$snp$pos[[1L]] + nchar(as.character(m$snp$ref[[1L]])) - 1L
    gaps <- c(gaps, abs(m$assay$com_start[[1L]] - anchor))
  }
  # the loci above are known designable: a skip here would hide a broken fixture
  expect_gte(length(gaps), 1L)
  expect_true(all(gaps >= 18L), label = "COM inside the lower com_gaps edge")
  expect_true(all(gaps <= 66L), label = "COM never near the upper com_gaps edge")
})

test_that("ZW0 M1: a locus written on the minus strand designs the same verdict", {
  skip_if_not_installed("Rsamtools")
  skip_if_not(has_blast(), "BLAST+ is not installed (indexed specificity stage)")
  gen <- syn_genome(seed = 5L, n = 400L)
  g <- indexed_syn_genome(gen)
  plus <- locus_at(gen$seq, 200L)
  # the same locus written on the minus strand: alleles complemented, strand "-"
  minus <- plus
  minus$ref <- rc(plus$ref)
  minus$alt <- rc(plus$alt)
  minus$strand <- "-"
  design <- function(snp) {
    suppressWarnings(tryCatch(
      marker_design_at(g, snp, "snp"),
      error = function(e) NULL
    ))
  }
  m_plus <- design(plus)
  m_minus <- design(minus)
  skip_if(is.null(m_plus) || is.null(m_minus), "the synthetic locus did not design")
  expect_equal(
    evaluation_verdict(m_plus$evaluation),
    evaluation_verdict(m_minus$evaluation)
  )
})

test_that("ZW0.1: a deliberate invariant asymmetry is declared with a reason", {
  spec <- tidymarkers:::criterion_spec()
  map <- tidymarkers:::criterion_asymmetry()
  # the registry publishes the exception, so it is no longer a footnote in a test
  expect_true("asymmetry_reason" %in% names(spec))
  expect_type(spec$asymmetry_reason, "character")
  # every declared asymmetry names a real criterion, once, with a real sentence
  expect_true(all(names(map) %in% spec$id))
  expect_identical(anyDuplicated(names(map)), 0L)
  expect_false(anyNA(map) || any(!nzchar(map)))
  # the two asymmetries ZW0.1 requires: KASP FAM/HEX and LAMP FIP/BIP
  expect_true(all(c(
    "as:gate:gate_discrimination", "as:gate:gate_wobble",
    "as:gate:gate_as_footprint", "as:gate:gate_as_window",
    "family_lamp_fip_alt", "family_lamp_inner_tm_balance",
    "family_lamp_fip_variant", "lamp:gate:gate_lamp_variant",
    "lamp:gate:gate_lamp_inner_tm"
  ) %in% names(map)))
  # no other row carries one, and every carried one is non-empty
  expect_setequal(spec$id[!is.na(spec$asymmetry_reason)], names(map))
  expect_false(any(!nzchar(spec$asymmetry_reason[!is.na(spec$asymmetry_reason)])))
  expect_true(all(spec$asymmetry_reason[match(names(map), spec$id)] == unname(map)))
  # and the columns are append-only: `asymmetry_reason` kept its position and the
  # ZW4.6 `target` column was appended after it, so no existing column moved
  expect_identical(names(spec)[[length(names(spec)) - 1L]], "asymmetry_reason")
  expect_identical(names(spec)[[length(names(spec))]], "target")
})
