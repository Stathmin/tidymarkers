# SNaPshot / single-base extension as a marker scheme: one flanking pair plus an
# extension primer whose 3' end abuts the variant.

# The scheme's declared gates in application order (grave / cheap first) and the
# soft subset (the ladder, least grave first).
snapshot_gates <- c(
  "gate_snapshot_placement", "gate_snapshot_comp", "gate_snapshot_unique",
  "gate_snapshot_structure", "gate_snapshot_tm", "gate_snapshot_amplicon"
)
snapshot_soft <- c("gate_snapshot_amplicon", "gate_snapshot_tm", "gate_snapshot_structure")

# A reproducible template with a mid-region substitution (seed 1 gives a valid
# design at pos 400 with flank 250 across the extension and amplification stages).
snap_W <- bseq(800L, 1L)
snap_genome <- Biostrings::DNAStringSet(c(Chr1 = snap_W))
snap_pos <- 400L
snap_flank <- 250L
snap_ref <- substr(snap_W, snap_pos, snap_pos)
snap_snp <- tibble::tibble(
  chrom = "Chr1", pos = snap_pos, ref = snap_ref,
  alt = setdiff(c("A", "C", "G", "T"), snap_ref)[1L]
)

test_that("the extension primer's 3' end sits one base 5' of the variant", {
  skip_if_no_oligoarrayaux()
  m <- design_snapshot_marker(snap_genome, snap_snp, flank = snap_flank)
  expect_s3_class(m, "kasp_snapshot")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  expect_equal(m$extension$placement, 1L)
  expect_match(m$extension$primer, "^[ACGT]+$")
  expect_gte(nchar(m$extension$primer), 18L)
  expect_lte(nchar(m$extension$primer), 25L)
  # whichever strand the extension primer is on, its 3' base abuts the variant
  if (m$extension$strand == "+") {
    expect_equal(m$extension$end, snap_snp$pos - 1L)
    expect_equal(m$extension$primer, substr(snap_W, m$extension$start, m$extension$end))
  } else {
    expect_equal(m$extension$start, snap_snp$pos + 1L)
    expect_equal(
      m$extension$primer,
      as.character(Biostrings::reverseComplement(
        Biostrings::DNAString(substr(snap_W, m$extension$start, m$extension$end))
      ))
    )
  }
  # the extension Tm lands in the declared window (or its bounded floor)
  expect_true(m$extension$tm >= 45 && m$extension$tm <= 60)
  expect_output(print(m), "kasp_snapshot")
})

test_that("every declared SNaPshot gate is evaluated by the designer (C7)", {
  skip_if_no_oligoarrayaux()
  m <- design_snapshot_marker(snap_genome, snap_snp, flank = snap_flank)
  expect_setequal(names(m$selection$gates), scheme_def_snapshot()$gates)
  expect_true(all(m$selection$gates))
  expect_length(m$selection$relaxed, 0L)
  # the extension primer really is a single perfect genomic site (uniqueness gate)
  expect_equal(exact_oligo_hits(m$extension$primer, snap_genome), 1L)
})

test_that("the extension-Tm soft gate relaxes to its floor, never past it (C8)", {
  # ask for a Tm window below the declared one: the winner's Tm is inside the
  # requested window but outside the gate's, so the ladder records the rung
  m <- design_snapshot_marker(snap_genome, snap_snp, flank = snap_flank, ext_tm = c(45, 47))
  expect_s3_class(m, "kasp_snapshot")
  expect_true("gate_snapshot_tm" %in% m$selection$relaxed)
  # past the floor the design is refused, so the rung is bounded
  expect_null(suppressMessages(
    design_snapshot_marker(snap_genome, snap_snp, flank = snap_flank, ext_tm = c(30, 32))
  ))
})

test_that("the flanking pair gives a short amplicon enclosing the variant", {
  skip_if_no_oligoarrayaux()
  m <- design_snapshot_marker(snap_genome, snap_snp, flank = snap_flank)
  expect_gte(m$amplification$product, 100L)
  expect_lte(m$amplification$product, 300L)
  expect_lt(m$amplicon$start, snap_snp$pos)
  expect_gt(m$amplicon$end, snap_snp$pos)
  # the flanking primers lie outside the variant and sequence back to the genome
  expect_equal(
    m$amplification$forward,
    substr(snap_W, m$amplification$start, m$amplification$start + nchar(m$amplification$forward) - 1L)
  )
  expect_equal(
    m$amplification$reverse,
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(
      substr(snap_W, m$amplification$end - nchar(m$amplification$reverse) + 1L, m$amplification$end)
    )))
  )
  expect_lt(m$amplification$start + nchar(m$amplification$forward) - 1L, snap_snp$pos)
  expect_gt(m$amplification$end - nchar(m$amplification$reverse) + 1L, snap_snp$pos)
})

test_that("a non-substitution target is refused with a message", {
  indel <- tibble::tibble(chrom = "Chr1", pos = 400L, ref = "AC", alt = "A")
  expect_message(m <- design_snapshot_marker(snap_genome, indel, flank = snap_flank))
  expect_null(m)
})

test_that("the SNaPshot DAG obeys the contract", {
  expect_dag_invariants("snapshot")
})

test_that("scheme_def_snapshot describes the SNaPshot scheme and its gates", {
  d <- scheme_def_snapshot()
  expect_equal(d$scheme, "snapshot")
  expect_equal(d$family, "snapshot")
  expect_equal(d$loci, 1L)
  expect_true(all(grepl("^gate_snapshot_", d$gates)))
  expect_setequal(d$gates, snapshot_gates)
  expect_setequal(d$gates, names(d$gate_short))
  expect_setequal(d$gates, names(d$gate_limit))
  # the ladder is declared least grave first and bounded by real floors
  expect_identical(d$relax, snapshot_soft)
  expect_setequal(names(d$gate_floor), d$relax)
  expect_true(all(d$gate_floor < 0))
  expect_setequal(names(d$tiers), c("tier1", "tier2", "tier3"))
  expect_setequal(names(d$gate_hard_label), setdiff(d$gates, d$relax))
  expect_true("gate_snapshot_placement" %in% names(d$gate_hard_label))
  expect_true("gate_snapshot_unique" %in% names(d$gate_hard_label))
  # the declared Tm window is the narrow one; the widening lives in the floor, so
  # the limit and the designer can no longer disagree
  expect_match(d$gate_limit[["gate_snapshot_tm"]], "50, 55")
})

test_that("the SNaPshot DAG shows no tier restating a hard gate, and its orders are declared", {
  s <- kasp_decision_schema("snapshot")
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  d <- scheme_def_snapshot()
  expect_true(all(diff(ord[paste0("snapshot:gate:", snapshot_gates)]) > 0L))
  hard <- s$nodes$name[s$nodes$kind == "gate" & s$nodes$hard]
  expect_setequal(hard, paste0("snapshot:gate:", setdiff(snapshot_gates, snapshot_soft)))
  # C14: tier1 no longer restates the hard placement gate
  expect_false(grepl("placement", unname(d$tiers[1L]), fixed = TRUE))
  txt <- kasp_decision_text("snapshot")
  expect_true(any(grepl("relaxation ladder \\(least grave first\\)", txt)))
  expect_true(any(grepl("in application order \\(grave / cheap first\\)", txt)))
})
