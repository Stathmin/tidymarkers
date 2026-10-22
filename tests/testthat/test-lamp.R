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

test_that("LAMP cross-dimer dG is evaluated at the LAMP temperature, not 37 C (Z3)", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  m <- design_lamp_marker(f$genome, f$snp, flank = 150L)
  lamp <- scheme_conditions("lamp")
  set <- unlist(m$primers)
  # the reported joint screen is the scheme-temperature evaluation
  expect_equal(m$compatibility$dg, primer_set_compatibility(set, conditions = lamp)$dg)
  # pair for pair, the same value dimer_thermo() returns at 65 C
  dg <- dimer_thermo(set, conditions = lamp)
  expect_equal(
    m$compatibility$dg,
    purrr::map2_dbl(m$compatibility$a, m$compatibility$b, function(x, y) {
      dg$dg[(dg$a == x & dg$b == y) | (dg$a == y & dg$b == x)][[1L]]
    })
  )
  cold <- primer_set_compatibility(set)
  # 65 C is the top of the LAMP window: every duplex is weaker than at 37 C, so
  # the screen is not failing a set on a 37 C-only number
  expect_true(all(m$compatibility$dg > cold$dg))
  # ... and the verdict is a base-pair-run screen: it does not read `conditions`,
  # so it (and so the ladder rung and the returned set) cannot move.  Compare
  # against the rung the set was accepted at: its `slack` widens the run rules
  # only, never the temperature.
  slack <- if ("gate_lamp_cross_dimer" %in% m$relaxed) LAMP_DIMER_SLACK else 0L
  expect_equal(m$compatibility$pass, lamp_pool_compatibility(m$primers, slack, conditions = NULL)$pass)
  # the helper's default is the LAMP scheme, and the ladder's slack only widens
  # the run rules (never the temperature)
  expect_equal(lamp_pool_compatibility(m$primers)$dg, m$compatibility$dg)
  relaxed <- lamp_pool_compatibility(m$primers, slack = 3L)
  expect_equal(relaxed$dg, m$compatibility$dg)
})

test_that("lamp_geometry accepts a well-formed inner-primer pair", {
  # plus-strand template coordinates for the fixture design: F2 ends at the
  # variant (300), B1/B2 sit downstream, and both 3' termini are G/C clamps.
  fip <- list(
    primer = "AGACCGCTAGAGTGCGCCAACGACAATATCCGCGAGTGAGATTC",
    arm5 = c(238L, 257L), arm3 = c(277L, 300L)
  )
  bip <- list(
    primer = "AGGTGCGGACATAGCTGTAGGATGCTCAACACAAATTCCCGCTGC",
    arm5 = c(324L, 346L), arm3 = c(363L, 384L)
  )
  g <- lamp_geometry(fip, bip, max_product = 300L)
  expect_true(g$ok)
  expect_equal(g$f2_b2_distance, 85L)
  expect_equal(g$loop_length, 23L)
  expect_equal(unname(g$loops), c(19L, 23L, 16L))
  expect_true(g$f2_b2_ok)
  expect_true(g$loop_ok)
  expect_true(g$fip_3p_stable)
  expect_true(g$bip_3p_stable)
  expect_true(all(g$checks))
  expect_length(g$failed, 0L)
  expect_length(g$messages, 0L)
})

test_that("lamp_geometry names a failed loop / 3'-end geometry", {
  # B1 starts on the F2 terminus (zero-length loop) and FIP's 3' end folds back
  # on itself (TTTTCAAAA), so both named checks must fail.
  fip <- list(primer = "TTTTCAAAA", arm5 = c(238L, 257L), arm3 = c(277L, 300L))
  bip <- list(primer = "ACGTACGTACGTACGCGG", arm5 = c(301L, 320L), arm3 = c(363L, 384L))
  g <- lamp_geometry(fip, bip, max_product = 300L)
  expect_false(g$ok)
  expect_true(g$f2_b2_ok)
  expect_false(g$loop_ok)
  expect_false(g$fip_3p_stable)
  expect_true(g$bip_3p_stable)
  expect_setequal(g$failed, c("loop_length", "fip_3p_stability"))
  expect_match(paste(g$messages, collapse = " "), "loop")
  expect_match(paste(g$messages, collapse = " "), "FIP 3' end")
})

test_that("lamp_geometry flags an F2-B2 span above the amplicon ceiling", {
  fip <- list(
    primer = "AGACCGCTAGAGTGCGCCAACGACAATATCCGCGAGTGAGATTC",
    arm5 = c(238L, 257L), arm3 = c(277L, 300L)
  )
  bip <- list(
    primer = "AGGTGCGGACATAGCTGTAGGATGCTCAACACAAATTCCCGCTGC",
    arm5 = c(324L, 346L), arm3 = c(363L, 384L)
  )
  g <- lamp_geometry(fip, bip, max_product = 60L)
  expect_false(g$f2_b2_ok)
  expect_true(g$loop_ok)
  expect_match(paste(g$messages, collapse = " "), "F2-B2 distance")
})

test_that("design_lamp_marker surfaces failed inner-primer geometry", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  snp <- f$snp
  # FIP's 3'-terminal base is the ALT by construction; a T ALT is not a G/C
  # clamp, so the set must be returned with a warning and `geometry$ok = FALSE`.
  snp$alt <- "T"
  m <- NULL
  expect_warning(
    m <- design_lamp_marker(f$genome, snp, flank = 150L),
    "geometry check failed"
  )
  expect_s3_class(m, "marker_lamp")
  expect_false(m$geometry$ok)
  expect_true("fip_3p_stability" %in% m$geometry$failed)
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

test_that("the LAMP variant gate re-checks FIP placement (not a decorative TRUE)", {
  fip <- list(arm3 = c(1L, 10L), primer = "AAAAC")
  expect_true(lamp_variant_gate(fip, "C", 10L))
  expect_false(lamp_variant_gate(fip, "A", 10L)) # FIP's 3' base is not ALT
  expect_false(lamp_variant_gate(fip, "C", 9L)) # the arm does not end at the variant
  expect_false(lamp_variant_gate(NULL, "C", 10L))
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

test_that("a LAMP design carries a paired REF set, not only an ALT one (Z9.1)", {
  # A REF-only template and a failed reaction both give no amplification, so an
  # ALT-only set cannot tell them apart.  The pair is the control construct: the
  # same primer positions with FIP's 3' base reading REF.
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  m <- design_lamp_marker(f$genome, f$snp, flank = 150L)
  p <- m$paired
  expect_type(p, "list")
  expect_identical(p$status, "complete")
  # the two FIPs are the same primer but for the single 3' base
  alt <- m$primers$FIP
  ref <- p$ref_fip
  expect_false(is.na(ref))
  expect_equal(nchar(ref), nchar(alt))
  expect_identical(substr(ref, 1L, nchar(ref) - 1L), substr(alt, 1L, nchar(alt) - 1L))
  expect_false(identical(substr(ref, nchar(ref), nchar(ref)), substr(alt, nchar(alt), nchar(alt))))
  # the discriminating base pairs the REF and ALT alleles
  expect_identical(p$discriminating_base, paste0(f$snp$ref, "/", f$snp$alt))
  # a pair is a construct, not a discrimination claim: nothing predicts a gap
  expect_false(any(grepl("gap|delta|ct", names(p), ignore.case = TRUE)))
})

test_that("the pair says why it could not be built instead of failing silently", {
  # an indel has no single 3' base to flip, so the helper reports the reason
  # instead of returning an empty pair
  snp <- list(ref = "AT", alt = "A")
  p <- lamp_paired_set(
    "ACGTACGT", 4L, list(primer = "ACGT"), snp, 20:24, c(55, 65), c(0.4, 0.6), 58
  )
  expect_identical(p$status, "not_applicable")
  expect_match(p$note, "single-base")
  expect_true(is.na(p$ref_fip))
})
