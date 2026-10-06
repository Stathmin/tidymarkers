# SNaPshot / single-base extension as a marker scheme: one flanking pair plus an
# extension primer whose 3' end abuts the variant.

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
  # the extension Tm lands in the window (or its single documented relaxation)
  expect_true(m$extension$tm >= 48 && m$extension$tm <= 60)
  expect_output(print(m), "kasp_snapshot")
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

test_that("scheme_def_snapshot describes the SNaPshot scheme and its gates", {
  d <- scheme_def_snapshot()
  expect_equal(d$scheme, "snapshot")
  expect_equal(d$family, "snapshot")
  expect_equal(d$loci, 1L)
  expect_true(all(grepl("^gate_snapshot_", d$gates)))
  expect_setequal(d$gates, names(d$gate_short))
  expect_setequal(d$gates, names(d$gate_limit))
  expect_true(all(d$relax %in% d$gates))
  expect_setequal(names(d$tiers), c("tier1", "tier2", "tier3"))
  expect_true("gate_snapshot_placement" %in% names(d$gate_hard_label))
})
