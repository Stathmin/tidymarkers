# Indel engine: allele templates, the junction scan, and end-to-end assembly.

test_that("allele_templates keeps the left flank and shifts the right one", {
  W <- Biostrings::DNAString(bseq(120, 7))
  snp_i <- 60L
  ref <- as.character(Biostrings::subseq(W, snp_i, snp_i + 2L)) # e.g. deletion block
  alt <- substr(ref, 1L, 1L)
  at <- allele_templates(W, snp_i, ref, alt)
  expect_equal(length(at$alt), length(W) + at$delta)
  expect_equal(
    as.character(Biostrings::subseq(at$alt, 1L, snp_i - 1L)),
    as.character(Biostrings::subseq(W, 1L, snp_i - 1L))
  )
  # the ALT block now sits where the REF block was
  expect_equal(as.character(Biostrings::subseq(at$alt, snp_i, snp_i)), alt)
})

test_that("the junction scan proposes discriminating primers for an insertion", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(bseq(160, 11))
  snp_i <- 80L
  ref <- as.character(Biostrings::subseq(W, snp_i, snp_i))
  alt <- paste0(ref, "TTG") # 2-nt insertion, VCF-anchored
  d <- design_as_primers_indel(W, snp_i, ref, alt, lengths = 20:22)
  expect_setequal(unique(d$allele), c("ref", "alt"))
  expect_gt(max(d$discrimination), 0)
  expect_true(all(c("primer", "on_eff", "off_eff", "as_start", "as_end") %in% names(d)))
})

test_that("the junction scan proposes discriminating primers for a deletion", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(bseq(160, 13))
  snp_i <- 80L
  ref <- as.character(Biostrings::subseq(W, snp_i, snp_i + 2L))
  alt <- substr(ref, 1L, 1L)
  d <- design_as_primers_indel(W, snp_i, ref, alt, lengths = 20:22)
  expect_gt(max(d$discrimination), 0)
})

test_that("design_as_primers_indel refuses equal-length alleles", {
  withr::local_envvar(LANGUAGE = "en")
  W <- Biostrings::DNAString(bseq(160, 17))
  expect_error(design_as_primers_indel(W, 80L, "A", "G"), "unequal")
})

test_that("assemble_kasp_assay designs an indel end-to-end", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(bseq(400, 4))
  # pick a locus whose neighbourhood has no long homopolymer (a real designer
  # would too - a run-rich junction is legitimately rejected by the checks)
  cand <- seq(160L, 240L, by = 1L)
  ok <- vapply(
    cand, function(p) !has_run(as.character(Biostrings::subseq(W, p - 35L, p + 35L))),
    logical(1)
  )
  pos <- cand[which(ok)[1L]]
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- as.character(Biostrings::subseq(W, pos, pos))
  snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = ref, alt = paste0(ref, "GGG"))
  a <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = NULL, flank = 150L))
  expect_equal(a$snp$ref, ref)
  expect_false(nchar(a$snp$ref) == nchar(a$snp$alt)) # indel => alleles differ in length
  expect_false(a$assay$as_ref[1] == a$assay$as_alt[1]) # the two AS primers differ
  expect_true(all(c("as_start", "as_alt_start") %in% names(a$assay)))
  chk <- marker_design_checks(a)
  expect_true(all(c("criterion", "observed", "limit", "pass") %in% names(chk)))
  expect_true(is.logical(chk$pass))
})

# ---- the indel DAG and its (formerly vacuous) wobble gate -------------------
test_that("the indel DAG obeys the contract and shares the KASP core stages", {
  expect_dag_invariants("indel")
  # KASP InDel discrimination is the SAME 3'-terminal mechanism as a SNP: the AS
  # primers end at the variant site (EasyKASP 2025), so the stage is deliberately
  # the shared KASP core, not a bespoke length-based one
  expect_kasp_core_dag("indel")
})

test_that("an indel is genuinely evaluated by the wobble gate (C7)", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(bseq(400, 4))
  ref <- as.character(Biostrings::subseq(W, 200L, 200L))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  snp <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = ref, alt = paste0(ref, "GGG"))
  a <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = NULL, flank = 150L))
  d <- a$as_candidates
  # the declared soft gate must carry a real flag, not NA (the old defect)
  expect_true(is.logical(d$wobble))
  expect_false(all(is.na(d$wobble)))
  # and it must match the documented rule: a substitution-like (equal-length)
  # junction is a real G:T / A:C wobble question; a length-changing bulge is not
  want <- nchar(d$ref_primer) == nchar(d$alt_primer) &
    wobble_mismatch(plus_3p_base(d$alt_primer, d$strand), plus_3p_base(d$ref_primer, d$strand))
  expect_identical(d$wobble, want)
})

test_that("the junction wobble asks about equal-length primers only", {
  # equal-length: the two 3' bases really do differ by one substitution
  expect_true(wobble_mismatch(plus_3p_base("ACGTACGTG", "+"), plus_3p_base("ACGTACGTA", "+")))
  # a '-' primer carries the complement, so the same pair is read on the plus strand
  expect_true(wobble_mismatch(plus_3p_base("ACGTACGTG", "-"), plus_3p_base("ACGTACGTA", "-")))
  # no wobble for a same-base terminal pair
  expect_false(wobble_mismatch(plus_3p_base("ACGTACGTG", "+"), plus_3p_base("ACGTACGTG", "+")))
})
