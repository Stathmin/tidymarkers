# TaqMan 5'-nuclease allelic discrimination: one flanking pair + two probes.

# The scheme's declared gates in application order (grave / cheap first) and the
# soft subset (the ladder, least grave first).
taqman_gates <- c(
  "gate_taqman_probe_comp", "gate_taqman_probe_specific", "gate_taqman_probe_structure",
  "gate_taqman_probe_tm", "gate_taqman_primer_tm", "gate_taqman_discrimination",
  "gate_taqman_centrality", "gate_taqman_amplicon"
)
taqman_soft <- c("gate_taqman_centrality", "gate_taqman_amplicon", "gate_taqman_probe_structure")

test_that("design_taqman_marker returns a short pair with two allele-specific probes", {
  f <- snp_fixture()
  genome <- f$genome
  ref <- f$ref
  alt <- f$alt
  snp <- f$snp

  m <- design_taqman_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "marker_taqman")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  # a short amplicon straddling the variant
  expect_lte(m$amplicon$product, 150L)
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  # exactly one probe per allele, different fluorophores
  expect_equal(nrow(m$probes), 2L)
  expect_setequal(m$probes$allele, c("REF", "ALT"))
  expect_equal(length(unique(m$probes$fluorophore)), 2L)
  # the two probes are the same window, differing only at the variant base
  expect_equal(m$probes$start[1L], m$probes$start[2L])
  expect_equal(m$probes$end[1L], m$probes$end[2L])
  # the REF probe carries the reference base at the variant position
  vp <- m$probes$variant_pos[1L]
  expect_equal(substr(m$probes$probe[m$probes$allele == "REF"], vp, vp), ref)
  expect_equal(substr(m$probes$probe[m$probes$allele == "ALT"], vp, vp), alt)
  # whole probes are ACGT, 5'-G free, and inside the amplicon
  expect_true(all(grepl("^[ACGT]+$", m$probes$probe)))
  expect_true(all(substr(m$probes$probe, 1L, 1L) != "G"))
  expect_true(all(m$probes$start >= m$amplicon$start & m$probes$end <= m$amplicon$end))
  expect_output(print(m), "marker_taqman")
  # the method itself (generic dispatch needs the package NAMESPACE, regenerated
  # by roxygen once the scheme is wired in)
  expect_output(print.marker_taqman(m), "<marker_taqman>")
  expect_output(print.marker_taqman(m), "FAM")
  expect_output(print.marker_taqman(m), "HEX")
})

test_that("probes lie strictly between the primer 3' ends, not merely in the amplicon", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  genome <- f$genome
  snp <- f$snp
  m <- design_taqman_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "marker_taqman")
  # the probe may not overlap either primer: it must sit inside the region
  # bounded by the two primer 3' ends
  inner_lo <- m$primers$start + nchar(m$primers$forward)
  inner_hi <- m$primers$end - nchar(m$primers$reverse)
  expect_true(all(m$probes$start >= inner_lo & m$probes$end <= inner_hi))
})

test_that("an ALT probe with a second perfect genomic site fails the specificity gate", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  W <- f$W
  genome <- f$genome
  snp <- f$snp
  m <- design_taqman_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "marker_taqman")
  ref_probe <- m$probes$probe[1L]
  alt_probe <- m$probes$probe[2L]
  expect_false(alt_probe == ref_probe)
  # the ALT probe's only perfect site is the ALT allele itself, so the reference
  # genome carries no off-target copy
  expect_equal(exact_oligo_hits(ref_probe, genome), 1L)
  expect_equal(exact_oligo_hits(alt_probe, genome), 0L)
  # append one off-target copy outside the design window (positions 151-450):
  # together with the ALT allele's on-target site that is a second genomic site,
  # so the hard specificity gate must refuse the design
  genome2 <- Biostrings::DNAStringSet(c(Chr1 = paste0(W, alt_probe)))
  expect_equal(exact_oligo_hits(alt_probe, genome2), 1L)
  expect_null(suppressMessages(design_taqman_marker(genome2, snp, flank = 150L)))
})

test_that("every declared TaqMan gate is evaluated by the designer (C7)", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(600L, 2L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "C")
  m <- design_taqman_marker(genome, snp, flank = 150L)
  expect_setequal(names(m$selection$gates), scheme_def_taqman()$gates)
  expect_true(all(m$selection$gates))
  expect_length(m$selection$relaxed, 0L)
  # the probe really is the one perfect genomic site (the specificity gate)
  expect_equal(exact_oligo_hits(m$probes$probe[1L], genome), 1L)
})

test_that("the mismatch dG gate is the real discrimination, not only centrality", {
  skip_if_no_oligoarrayaux()
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(600L, 2L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "C")
  m <- design_taqman_marker(genome, snp, flank = 150L)
  ddg <- taqman_mismatch_dg(
    m$probes$probe[1L], m$probes$variant_pos[1L], "C"
  )
  # the mismatched template must cost the probe duplex real stability
  expect_gt(ddg, 0)
  expect_true(m$selection$gates[["gate_taqman_discrimination"]])
})

test_that("the centrality soft gate relaxes to its interior floor, never past it (C8)", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(600L, 2L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "C")
  # demand a centrality the locus cannot quite reach: the winner sits inside the
  # floor and the ladder records the rung instead of failing silently
  m <- design_taqman_marker(genome, snp, flank = 150L, min_centrality = 0.99)
  expect_s3_class(m, "marker_taqman")
  expect_lt(m$probes$centrality[1L], 0.99)
  expect_true("gate_taqman_centrality" %in% m$selection$relaxed)
  # past the floor the design is refused, so the rung is bounded
  expect_null(suppressMessages(
    design_taqman_marker(genome, snp, flank = 150L, min_centrality = 1.5)
  ))
})

test_that("probes are centred over the variant (maximal destabilisation)", {
  f <- snp_fixture()
  genome <- f$genome
  snp <- f$snp
  m <- design_taqman_marker(genome, snp, flank = 150L)
  expect_gte(m$probes$centrality[1L], 0.5)
  # centrality 1 = exactly central; the variant must sit near, not at, a probe end
  center <- (nchar(m$probes$probe[1L]) + 1L) / 2
  expect_lte(abs(m$probes$variant_pos[1L] - center), 2)
})

test_that("design_taqman_marker refuses non-substitutions with a message", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(400L, 3L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = "AC", alt = "GT")
  expect_null(suppressMessages(design_taqman_marker(genome, snp)))
})

test_that("the TaqMan DAG obeys the contract", {
  expect_dag_invariants("taqman")
})

test_that("the taqman scheme definition has the gates, tiers and labels wired", {
  d <- scheme_def_taqman()
  expect_equal(d$scheme, "taqman")
  expect_equal(d$family, "taqman")
  expect_equal(d$loci, 1L)
  expect_true(all(c(
    "scheme", "family", "label", "discrimination", "readout", "zygosity",
    "amplicon", "loci", "verdict", "gates", "relax", "tiers",
    "gate_short", "gate_limit", "gate_floor", "gate_floor_label",
    "gate_hard_label", "stage_label"
  ) %in% names(d)))
  expect_setequal(d$gates, taqman_gates)
  expect_identical(d$relax, taqman_soft)
  expect_setequal(names(d$gate_short), d$gates)
  expect_setequal(names(d$gate_limit), d$gates)
  expect_true(all(names(d$gate_floor_label) %in% d$relax))
  expect_true(all(names(d$gate_hard_label) %in% setdiff(d$gates, d$relax)))
  expect_setequal(names(d$gate_floor), d$relax)
  expect_true(all(d$gate_floor < 0))
  expect_equal(unname(d$tiers[1L]), "probe discrimination (mismatch dG margin)")
  # the probe-Tm *ordering* is a ranking (tier2), never a gate: the header's
  # "8-10 C above the primers" rule cannot be hard on this Tm scale
  expect_false(grepl("hottest", d$gate_limit[["gate_taqman_probe_tm"]], fixed = TRUE))
})

test_that("the TaqMan DAG shows the gates in application order, ordering only as a tier", {
  s <- marker_decision_schema("taqman")
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(diff(ord[paste0("taqman:gate:", taqman_gates)]) > 0L))
  hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
  expect_setequal(hard, paste0("taqman:gate:", setdiff(taqman_gates, taqman_soft)))
  # the Tm ordering survives as tier2 (non-vacuous: a continuous margin, C14)
  expect_match(unname(scheme_def_taqman()$tiers[2L]), "Tm ordering")
  expect_true(any(grepl(
    "relaxation ladder \\(least grave first\\)",
    marker_decision_text("taqman")
  )))
})

test_that("the mismatch dG threshold is a package rule and says so", {
  # Q1: no open paper in this package prints a dG threshold for TaqMan probe
  # discrimination, nor a matched-vs-mismatched dCt.  The number stays (it is
  # what this package scores against) but it is never presented as measured.
  expect_equal(TAQMAN_MISMATCH_DG, 1.5)
  # the criterion registry is the one place every reader (gate, check, plot)
  # takes the limit from, so the disclaimer lives there
  cr <- marker_criteria(scheme = "taqman")
  lim <- cr$limit[cr$id == "taqman:gate:gate_taqman_discrimination"]
  expect_length(lim, 1L)
  expect_match(lim, "package rule, not measured", fixed = TRUE)
})
