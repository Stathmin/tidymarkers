# castPCR / competitive allele-specific TaqMan: the assay detects the
# low-abundance MUTANT (ALT), so the allele-specific primer reads the mutant base
# and a WT-specific MGB blocker suppresses the wild type (Roma 2013; Barbano
# 2015).  The hydrolysis probe is *locus-specific*, not allele-specific.

castpcr_fixture <- function() {
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  list(genome = genome, snp = snp)
}

# the base a plus-strand oligo carries at a position, on the ASP's own strand
on_strand <- function(base, role) {
  if (identical(role, "forward")) base else as.character(Biostrings::complement(Biostrings::DNAString(base)))
}

test_that("design_castpcr_marker detects the MUTANT allele, with a WT blocker and a locus probe", {
  skip_if_no_oligoarrayaux()
  f <- castpcr_fixture()
  m <- design_castpcr_marker(f$genome, f$snp, flank = 150L)
  expect_s3_class(m, "kasp_castpcr")
  expect_equal(m$scheme, "castpcr")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))

  ref <- f$snp$ref
  alt <- f$snp$alt
  role <- m$primers$as_role
  expect_true(role %in% c("forward", "reverse"))

  # THE mutant anchor: the ASP's 3' base is the ALT base (and so not the WT base)
  as <- m$primers$as_primer
  last <- substr(as, nchar(as), nchar(as))
  expect_equal(last, on_strand(alt, role))
  expect_false(identical(last, on_strand(ref, role)))

  # the WT-specific blocker reads the WT base, so its own 3' mismatch spares the
  # mutant; it is a different oligo from the ASP
  expect_equal(substr(m$blocker$sequence, nchar(m$blocker$sequence), nchar(m$blocker$sequence)), on_strand(ref, role))
  expect_false(identical(m$blocker$sequence, as))
  expect_equal(m$blocker$allele, "WT")

  # the probe is locus-specific: inside the amplicon, between the primers and
  # never overlapping the ASP - so it cannot cover the variant the ASP anchors
  fwd_end <- m$amplicon$start + nchar(m$primers$forward) - 1L
  rev_start <- m$amplicon$end - nchar(m$primers$reverse) + 1L
  expect_equal(if (identical(role, "forward")) fwd_end else rev_start, 300L)
  expect_gte(m$probe$start, m$amplicon$start)
  expect_lte(m$probe$end, m$amplicon$end)
  if (identical(role, "forward")) {
    expect_gt(m$probe$start, fwd_end)
  } else {
    expect_lt(m$probe$end, rev_start)
  }
  expect_true(m$probe$start > 300L || m$probe$end < 300L)
  # the probe melts above the ASP: the hard ordering gate, the same rule dPCR
  # hard-gates (gate_castpcr_probe_tm <-> gate_dpcr_probe_tm)
  expect_gte(m$probe$tm, m$primers$as_tm + CASTPCR_PROBE_MARGIN - 1e-9)

  # the amplicon is short: inside the strict window, or the widened ladder bar it
  # exchanged, and the record says which
  lo <- KASP_LIMITS$amplicon[1L] - CASTPCR_AMPLICON_FLOOR
  hi <- KASP_LIMITS$amplicon[2L] + CASTPCR_AMPLICON_FLOOR
  expect_true(dplyr::between(m$primers$product, lo, hi))
  expect_true(all(m$selection$relaxed %in% c("gate_castpcr_amplicon", "gate_castpcr_as_tm")))
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

test_that("the castpcr DAG gates the mutant anchor, the WT blocker and the probe ordering", {
  d <- scheme_def_castpcr()
  s <- kasp_decision_schema("castpcr")
  # a single AS primer cannot resolve REF/ALT zygosity: castPCR is mutant
  # detection (dominant), not co-dominant
  expect_equal(d$zygosity, "dominant (mutant detection)")
  expect_false(grepl("co-dominant", d$zygosity))
  expect_equal(d$loci, 1L)
  expect_match(d$label, "mutant detection")
  # the WT blocker is a real gate, and the probe ordering is HARD - the same
  # hardness the physically identical dPCR probe gate carries
  expect_true("gate_castpcr_blocker" %in% d$gates)
  expect_false("gate_castpcr_probe_tm" %in% d$relax)
  expect_false("gate_dpcr_probe_tm" %in% scheme_def_dpcr()$relax)
  hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
  expect_true("castpcr:gate:gate_castpcr_probe_tm" %in% hard)
  expect_true("castpcr:gate:gate_castpcr_blocker" %in% hard)
  # the mutant anchor and its 3' mismatch are hard: without them the assay is not
  # allele-specific, so they are never relaxed
  expect_setequal(setdiff(d$gates, d$relax), c(
    "gate_castpcr_as3", "gate_castpcr_mismatch", "gate_castpcr_blocker", "gate_castpcr_probe_tm"
  ))
  # the amplicon gate's STRICT window is the acceptance check's window (C5)
  crit <- kasp_criteria(scheme = "castpcr")
  expect_match(crit$limit[crit$id == "castpcr:gate:gate_castpcr_amplicon"], "\\[50, 100\\] bp")
  # the gate order is the declared application order (grave / cheap first)
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(diff(ord[paste0("castpcr:gate:", d$gates)]) > 0L))
  # tiers rank the minor margins, they do not restate a hard gate (C14)
  expect_false(any(grepl("blocker|anchor", d$tiers, ignore.case = TRUE)))
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

test_that("the castpcr DAG obeys the contract", {
  expect_dag_invariants("castpcr")
})
