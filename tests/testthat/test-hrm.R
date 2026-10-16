# HRM as a marker scheme: one flanking pair over a short amplicon.
#
# Design rules (see R/marker_hrm.R): the short, single-domain amplicon is HRM's
# *core* requirement (Liew et al., Clin Chem 2004, 50:1156-64, "small
# amplicons"; Reed & Wittwer, Clin Chem 2004, 50:1748-54), so its gate is HARD
# and its window is exactly the marker-agnostic amplicon check's window (C5).
# The flanking pair's Tm balance is the soft gate; the ladder may widen it.  The
# variant must sit in the melt domain and under NEITHER primer - that is enforced
# structurally by the flanking geometry, not by a gate.

# the fixture every test designs against: a reproducible 600 bp template with a
# SNP at the centre
hrm_fixture <- function() {
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  list(genome = genome, snp = snp)
}

test_that("design_hrm_marker returns a short flanking pair with the variant inside", {
  skip_if_no_oligoarrayaux()
  f <- hrm_fixture()
  m <- design_hrm_marker(f$genome, f$snp, flank = 150L)
  expect_s3_class(m, "marker_hrm")
  # the variant is inside the amplicon and under neither primer (structural):
  # the forward primer ends before it, the reverse primer starts after it
  # (`$primers$start/end` are template coordinates; the genomic footprint follows
  # from the amplicon bounds plus the primer lengths)
  fwd_end <- m$amplicon$start + nchar(m$primers$forward) - 1L
  rev_start <- m$amplicon$end - nchar(m$primers$reverse) + 1L
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  expect_lt(fwd_end, 300L)
  expect_gt(rev_start, 300L)
  # HRM wants a *short* melt domain: the product is inside the marker-agnostic
  # amplicon window (the hard gate's window, C5)
  expect_true(dplyr::between(m$primers$product, KASP_LIMITS$amplicon[1L], KASP_LIMITS$amplicon[2L]))
  expect_true(grepl("^[ACGT]+$", m$primers$forward) && grepl("^[ACGT]+$", m$primers$reverse))
  # the strict design relaxes nothing, and the record says so
  expect_equal(m$selection$relaxed, character(0))
  # the marker carries its scheme, so a graph/report rebuilds the right tree (C12)
  expect_equal(m$scheme, "hrm")
  expect_output(print(m), "marker_hrm")
})

test_that("the amplicon gate is hard and the Tm balance is the ladder's only rung", {
  # the swap the audit asked for: HRM's core requirement (the short amplicon) is
  # the HARD gate, the Tm balance the soft one
  expect_equal(KASP_GATES_HRM, c("gate_hrm_amplicon", "gate_hrm_tm"))
  expect_equal(KASP_RELAX_HRM, "gate_hrm_tm")
  s <- marker_decision_schema("hrm")
  hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
  expect_equal(hard, "hrm:gate:gate_hrm_amplicon")
  # a failed hard amplicon gate dooms directly (C4)
  expect_true(any(s$edges$kind == "doomed" & s$edges$from == "hrm:gate:gate_hrm_amplicon"))
  # the amplicon gate's window is the acceptance check's window, so a marker the
  # gate admits is never rejected by the check (C5)
  crit <- marker_criteria(scheme = "hrm")
  lim <- crit$limit[crit$id == "hrm:gate:gate_hrm_amplicon"]
  expect_match(lim, "\\[50, 100\\] bp")
  # the gate order is the declared application order (grave / cheap first): the
  # amplicon (the grave one) before the Tm window
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_lt(ord[["hrm:gate:gate_hrm_amplicon"]], ord[["hrm:gate:gate_hrm_tm"]])
  # tiers rank what is left after the hard gate, they do not restate it (C14)
  expect_false(any(grepl("amplicon", KASP_TIERS_HRM, ignore.case = TRUE)))
})

test_that("the HRM ladder widens the Tm window and reports it", {
  skip_if_no_oligoarrayaux()
  f <- hrm_fixture()
  # a one-degree Tm window admits no pair, so the ladder must widen it by
  # HRM_TM_FLOOR and record the relaxation (the gate is evaluated, not decorative)
  m <- design_hrm_marker(f$genome, f$snp, flank = 150L, tm_range = c(63, 63))
  expect_s3_class(m, "marker_hrm")
  expect_equal(m$selection$relaxed, "gate_hrm_tm")
  # the hard amplicon gate is never relaxed: an impossible window dooms the design
  expect_message(
    expect_null(design_hrm_marker(f$genome, f$snp, flank = 150L, min_product = 300L)),
    "amplicon window"
  )
})

test_that("the hrm scheme is registered with its own decision stage and is designable", {
  expect_true("hrm" %in% marker_schemes()$scheme)
  s <- marker_decision_schema("hrm")
  expect_true("hrm" %in% s$nodes$stage)
  expect_true(is_target_designable(structure(list(), class = "marker_hrm")))
})

test_that("the hrm DAG obeys the contract", {
  expect_dag_invariants("hrm")
})

test_that("plot_hrm_amplicon draws the marker and rejects other objects", {
  skip_if_no_oligoarrayaux()
  f <- hrm_fixture()
  m <- design_hrm_marker(f$genome, f$snp, flank = 150L)
  expect_s3_class(plot_hrm_amplicon(m), "ggplot")
  expect_error(plot_hrm_amplicon(list()))
})
