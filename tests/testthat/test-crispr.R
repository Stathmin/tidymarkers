# CRISPR-Dx (SHERLOCK / DETECTR): the guide design and the classed marker.

test_that("design_crispr_marker builds a Cas12a guide with a TTTV PAM and an amplicon", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 7L)
  substr(W, 293L, 296L) <- "TTTA" # 5'-TTTV PAM immediately 5' of the protospacer
  substr(W, 297L, 316L) <- "ACGTGCTAGCTAGCATGCA" # GC-balanced 20 nt spacer
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "T", alt = "A")
  m <- design_crispr_marker(genome, snp, flank = 200L, cas = "cas12a")
  expect_s3_class(m, "marker_crispr")
  expect_equal(m$cas, "cas12a")
  expect_match(m$pam, "^TTT[ACG]$")
  expect_gte(nchar(m$guide), 20L)
  expect_lte(nchar(m$guide), 23L)
  expect_lt(m$seed_offset, 8L) # the variant sits in the PAM-proximal seed
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  expect_gte(m$amplicon$product, 80L)
  expect_lte(m$amplicon$product, 250L)
  expect_output(print(m), "marker_crispr")
})

test_that("design_crispr_marker builds a PAM-free Cas13 guide", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  genome <- f$genome
  snp <- f$snp
  m <- design_crispr_marker(genome, snp, flank = 200L, cas = "cas13")
  expect_s3_class(m, "marker_crispr")
  expect_equal(m$cas, "cas13")
  expect_true(is.na(m$pam)) # Cas13 targets RNA: no PAM
  expect_true(is.na(m$gates[["gate_crispr_pam"]])) # PAM-free: the gate is NA
  expect_false(is.na(m$amplicon$product)) # ... and the marker still designs
  expect_gte(nchar(m$guide), 22L)
  expect_lte(nchar(m$guide), 28L)
  expect_output(print(m), "marker_crispr")
})

test_that("design_crispr_marker informs and returns NULL when no guide passes", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("AT", 300L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "C")
  expect_message(m <- design_crispr_marker(genome, snp, flank = 200L, cas = "cas12a"), "guide")
  expect_null(m)
})

test_that("crispr_guide_candidates honours the strand restriction", {
  W <- bseq(300L, 3L)
  args <- list(seq = W, snp_i = 150L, cas = "cas13", spacer_lengths = 22:28, gc_range = c(0.35, 0.65))
  both <- do.call(crispr_guide_candidates, args)
  plus <- do.call(crispr_guide_candidates, c(args, list(strands = "+")))
  minus <- do.call(crispr_guide_candidates, c(args, list(strands = "-")))
  expect_setequal(both$strand, c("+", "-"))
  expect_true(all(plus$strand == "+"))
  expect_true(all(minus$strand == "-"))
  expect_equal(nrow(plus) + nrow(minus), nrow(both))
})

test_that("scheme_def_crispr declares its prefixed gates, tiers and lookups", {
  def <- scheme_def_crispr()
  expect_equal(def$scheme, "crispr")
  expect_equal(def$family, "crispr")
  expect_true(all(grepl("^gate_crispr_", def$gates)))
  expect_setequal(def$gates, names(def$gate_short))
  expect_setequal(def$gates, names(def$gate_limit))
  expect_true(all(def$relax %in% def$gates))
  expect_setequal(names(def$tiers), c("guide", "seed", "specificity", "amplicon"))
  expect_true(all(names(def$gate_floor_label) %in% def$relax))
  expect_true(all(names(def$gate_hard_label) %in% setdiff(def$gates, def$relax)))
})

test_that("guide specificity is a hard gate", {
  # an off-target protospacer could license the reporter on the wrong template,
  # so this gate is never relaxed (unlike the amplicon / GC / spacer gates).
  expect_true(
    "gate_crispr_specificity" %in% setdiff(scheme_def_crispr()$gates, scheme_def_crispr()$relax)
  )
})

test_that("design_crispr_marker supports only the RPA pre-amplification", {
  # `lamp` was accepted but only reused the flanking-pair designer, which is not
  # a LAMP reaction (no FIP/BIP): the choice is gone until a LAMP engine exists.
  expect_false(any(grepl("lamp", deparse(formals(design_crispr_marker)$preamp))))
  expect_identical(formals(design_crispr_marker)$preamp, "rpa")
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("AT", 300L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "C")
  expect_error(design_crispr_marker(genome, snp, preamp = "lamp"), "rpa")
})

test_that("the crispr DAG obeys the contract", expect_dag_invariants("crispr"))
