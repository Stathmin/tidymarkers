# dPCR as a marker scheme: a hydrolysis probe over a short amplicon.

# The scheme's declared gates, in application order (grave / cheap first), and the
# subscript of the ladder (least grave first).  Kept here so the DAG tests below
# read one vocabulary.
dpcr_gates <- c(
  "gate_dpcr_probe_comp", "gate_dpcr_probe_specific", "gate_dpcr_probe_structure",
  "gate_dpcr_probe_tm", "gate_dpcr_primer_tm", "gate_dpcr_amplicon"
)
dpcr_soft <- c("gate_dpcr_amplicon", "gate_dpcr_probe_structure")

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
  # two allele-specific probes (one per allele, FAM/HEX), both spanning the variant
  expect_equal(nrow(m$probes), 2L)
  expect_setequal(m$probes$allele, c("REF", "ALT"))
  expect_setequal(m$probes$fluorophore, c("FAM", "HEX"))
  expect_gt(m$probes$start[1L], m$amplicon$start)
  expect_lt(m$probes$end[1L], m$amplicon$end)
  expect_lte(m$probes$start[1L], 300L)
  expect_gte(m$probes$end[1L], 300L)
  # the ALT probe differs from the REF probe at exactly the variant position
  expect_equal(m$probes$variant_pos[1L], m$probes$variant_pos[2L])
  differs <- which(strsplit(m$probes$probe[1L], "")[[1L]] != strsplit(m$probes$probe[2L], "")[[1L]])
  expect_equal(differs, m$probes$variant_pos[1L])
  # the probe's composition rules (18-25 nt, GC 40-60%, no 5' G, no run > 4)
  plen <- nchar(m$probes$probe[1L])
  expect_gte(plen, 18L)
  expect_lte(plen, 25L)
  expect_match(m$probes$probe[1L], "^[ACGT]+$")
  expect_false(substr(m$probes$probe[1L], 1L, 1L) == "G")
  expect_true(gc_frac(m$probes$probe[1L]) >= 0.40 && gc_frac(m$probes$probe[1L]) <= 0.60)
  expect_false(has_run(m$probes$probe[1L], n = 5L))
  # the probe out-competes the primers (the ~8-10 C rule)
  expect_gte(m$probes$tm[1L] - max(m$primers$tm_f, m$primers$tm_r), 8)
  expect_match(m$primers$forward, "^[ACGT]+$")
  expect_match(m$primers$reverse, "^[ACGT]+$")
  expect_output(print.kasp_dpcr(m), "kasp_dpcr")
})

test_that("every declared dPCR gate is evaluated by the designer (C7)", {
  skip_if_no_oligoarrayaux()
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(600L, 1L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "T", alt = "A")
  m <- design_dpcr_marker(genome, snp, flank = 150L)
  # the engine scores exactly the scheme's declared gates, and a clean design
  # clears every one of them without relaxing anything
  expect_setequal(names(m$selection$gates), scheme_def_dpcr()$gates)
  expect_true(all(m$selection$gates))
  expect_length(m$selection$relaxed, 0L)
  expect_length(m$selection$violated, 0L)
  # the probe really is a single perfect genomic site (the specificity gate)
  expect_equal(exact_oligo_hits(m$probes$probe[1L], genome), 1L)
})

test_that("the amplicon soft gate relaxes to its floor instead of silently failing (C6/C7)", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(600L, 1L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "T", alt = "A")
  # ask for a window outside the declared [60, 150] bp: the winner is inside the
  # requested window but outside the gate's, so the ladder must record the rung
  m <- design_dpcr_marker(genome, snp, flank = 150L, min_product = 155L, max_product = 160L)
  expect_s3_class(m, "kasp_dpcr")
  expect_true(m$primers$product > 150L)
  expect_true("gate_dpcr_amplicon" %in% m$selection$relaxed)
  expect_true("gate_dpcr_amplicon" %in% m$selection$violated)
  # and the gate's floor really bounds it: a window past the floor is refused
  expect_null(suppressMessages(
    design_dpcr_marker(genome, snp, flank = 250L, min_product = 200L, max_product = 220L)
  ))
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

test_that("the dPCR DAG obeys the contract", {
  expect_dag_invariants("dpcr")
})

test_that("the dPCR scheme definition carries its gates, relax and tiers", {
  d <- scheme_def_dpcr()
  expect_identical(d$scheme, "dpcr")
  expect_identical(d$family, "dpcr")
  expect_identical(d$loci, 1L)
  expect_identical(d$verdict, "usable dPCR assay")
  expect_setequal(d$gates, dpcr_gates)
  # only the amplicon window and the probe's self-structure may be lowered, and
  # the ladder is declared least grave first
  expect_identical(d$relax, dpcr_soft)
  expect_identical(names(d$tiers), c("tier1", "tier2", "tier3"))
  expect_true(all(d$gates %in% names(d$gate_short)))
  expect_true(all(d$gates %in% names(d$gate_limit)))
  expect_true(all(d$relax %in% names(d$gate_floor_label)))
  expect_true(all(setdiff(d$gates, d$relax) %in% names(d$gate_hard_label)))
  # the floors are bounded and reachable (never -Inf), keyed by the soft gates
  expect_setequal(names(d$gate_floor), d$relax)
  expect_true(all(d$gate_floor < 0))
  expect_identical(d$stage_label, "dPCR (droplet digital PCR)")
  # the readout does not claim a partitioning/Poisson model this package never
  # computes (see kasp_readout())
  expect_false(grepl("Poisson", d$readout, fixed = TRUE))
})

test_that("the dPCR DAG shows the gates in application order and the ladder least grave first", {
  s <- kasp_decision_schema("dpcr")
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  # gate order follows the declared (grave / cheap first) order
  expect_true(all(diff(ord[paste0("dpcr:gate:", dpcr_gates)]) > 0L))
  # hardness is exactly "not on the ladder"
  hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
  expect_setequal(hard, paste0("dpcr:gate:", setdiff(dpcr_gates, dpcr_soft)))
  soft <- s$nodes$name[s$nodes$kind == "gate" & !s$nodes$hard & grepl("^dpcr:gate:", s$nodes$name)]
  expect_setequal(soft, paste0("dpcr:gate:", dpcr_soft))
  # the refine node names the rungs in the declared relax order
  refine <- s$nodes$detail[s$nodes$kind == "relax"]
  expect_match(refine, "amplicon")
  expect_true(grepl("least grave first", paste(kasp_decision_text("dpcr"), collapse = "\n"), fixed = TRUE))
})
