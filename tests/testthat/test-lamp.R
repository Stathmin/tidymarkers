# LAMP as a marker scheme: F3/B3 outside, FIP/BIP inside, the variant at FIP's 3' end.

test_that("design_lamp_marker returns a 4-primer LAMP set with the variant at FIP's 3' end", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  genome <- f$genome
  snp <- f$snp
  m <- design_lamp_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "marker_lamp")
  expect_true(all(c("F3", "B3", "FIP", "BIP") %in% names(m$primers)))
  # inner (two-arm fusion) primers are long, the outer F3/B3 short (LAMP rule)
  expect_gte(nchar(m$primers$FIP), 40L)
  expect_lte(nchar(m$primers$FIP), 50L)
  expect_gte(nchar(m$primers$BIP), 40L)
  expect_lte(nchar(m$primers$BIP), 50L)
  expect_gte(nchar(m$primers$F3), 21L)
  expect_lte(nchar(m$primers$F3), 24L)
  # the four primers are screened jointly, not each one alone
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  expect_true(all(grepl("^[ACGT]+$", unlist(m$primers))))
  # amplicon within the LAMP window and straddling the variant
  product <- m$amplicon$end - m$amplicon$start + 1L
  expect_gte(product, 120L)
  expect_lte(product, 300L)
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  # the variant is FIP's 3'-terminal base; the F2-B2 loop lies downstream of it
  expect_equal(m$variant_primer, "FIP")
  expect_gte(m$loop$start, 301L)
  # the inner primers are two-arm fusions, longer than the outer arms
  expect_gt(nchar(m$primers$FIP), nchar(m$primers$F3))
  expect_gt(nchar(m$primers$BIP), nchar(m$primers$B3))
  expect_output(print.marker_lamp(m), "marker_lamp")
  expect_output(print.marker_lamp(m), "FIP:")
})

test_that("lamp_fuse reverse-complements the 5' loop arm in both orientations", {
  # non-palindromic arms: the 5' arm is A/C-only, the 3' arm A/C-only, so each
  # fused primer stays free of inverted repeats (the G/T complement cannot pair).
  five_seq <- "AAAACCCC"
  three_seq <- "ACACCAACAC"
  five <- tibble::tibble(start = 1L, end = 8L, primer = five_seq, tm = 62)
  three <- tibble::tibble(start = 31L, end = 40L, primer = three_seq, tm = 58)

  fwd <- lamp_fuse(five, three, target_tm = 60, min_gap = 16L, reverse = FALSE)
  rev <- lamp_fuse(five, three, target_tm = 60, min_gap = 16L, reverse = TRUE)

  # FIP = revcomp(F1 site) + F2 site: its 3' displacement arm stays plus-strand
  expect_equal(fwd$primer, paste0(lamp_revcomp(five_seq), three_seq))
  expect_true(startsWith(fwd$primer, lamp_revcomp(five_seq)))
  expect_true(endsWith(fwd$primer, three_seq))
  # BIP = revcomp(B1 site) + revcomp(B2 site): both arms are minus-strand
  expect_equal(rev$primer, paste0(lamp_revcomp(five_seq), lamp_revcomp(three_seq)))
  expect_true(startsWith(rev$primer, lamp_revcomp(five_seq)))
  expect_true(endsWith(rev$primer, lamp_revcomp(three_seq)))
  # the old plus+plus FIP (five + three) must no longer be produced
  expect_false(identical(fwd$primer, paste0(five_seq, three_seq)))
  expect_false(identical(fwd$primer, rev$primer))
})

test_that("B3 is a minus-strand oligo and FIP carries the ALT base at its 3' end", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  W <- f$W
  genome <- f$genome
  ref <- f$ref
  snp <- f$snp
  m <- design_lamp_marker(genome, snp, flank = 150L)

  # (a) B3 anneals to the plus strand: it is the reverse complement of its
  # plus-strand window, whose end is the amplicon end.
  b3_end <- m$amplicon$end
  b3_start <- b3_end - nchar(m$primers$B3) + 1L
  expect_equal(
    m$primers$B3,
    as.character(Biostrings::reverseComplement(
      Biostrings::DNAString(substr(W, b3_start, b3_end))
    ))
  )

  # (b) FIP's 3' arm (F2) stays plus-strand and ends at the variant, carrying
  # the ALT base - the allele-specific block.
  fip_3p <- substr(m$primers$FIP, nchar(m$primers$FIP), nchar(m$primers$FIP))
  expect_equal(fip_3p, m$snp$alt)
  expect_false(identical(fip_3p, m$snp$ref))

  # (c) the fusion did not introduce non-ACGT characters or change the lengths
  expect_match(m$primers$FIP, "^[ACGT]+$")
  expect_match(m$primers$BIP, "^[ACGT]+$")
  expect_gte(nchar(m$primers$FIP), 40L)
  expect_lte(nchar(m$primers$FIP), 50L)
  expect_gte(nchar(m$primers$BIP), 40L)
  expect_lte(nchar(m$primers$BIP), 50L)
})

test_that("print.marker_lamp names the four primers and the loop", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  genome <- f$genome
  snp <- f$snp
  m <- design_lamp_marker(genome, snp, flank = 150L)
  expect_output(print.marker_lamp(m), "variant at the FIP 3' end")
  expect_output(print.marker_lamp(m), "B3:")
  expect_output(print.marker_lamp(m), "BIP:")
})

test_that("scheme_def_lamp declares the LAMP stage consistently", {
  s <- scheme_def_lamp()
  expect_equal(s$scheme, "lamp")
  expect_equal(s$family, "lamp")
  expect_equal(s$loci, 1L)
  expect_true(all(s$gates %in% names(s$gate_short)))
  expect_true(all(s$gates %in% names(s$gate_limit)))
  expect_true(all(s$relax %in% s$gates))
  expect_true(all(s$gates %in% names(s$gate_floor_label) | !(s$gates %in% s$relax)))
  expect_identical(unname(s$tiers), c("variant discrimination", "inner Tm set", "primer structure", "amplicon length"))
})

test_that("the lamp DAG obeys the contract", expect_dag_invariants("lamp"))
