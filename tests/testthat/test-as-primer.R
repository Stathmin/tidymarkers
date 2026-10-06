# Invariants for the KASP core: the SNP must land on the 3'-terminal base, and
# the DECIPHER "target" must be the reverse complement of the primer.  A silent
# swap here would invert the discrimination metric, so it is tested explicitly.

W <- Biostrings::DNAString(paste(rep("ACGT", 19L), collapse = "")) # 76 bp
SNP_I <- 30L

seg_with_allele <- function(L, strand, allele) {
  seg <- as.character(as_window(W, SNP_I, L, strand))
  idx <- if (identical(strand, "+")) L else 1L
  substr(seg, idx, idx) <- allele
  Biostrings::DNAString(seg)
}

test_that("make_as_primer puts the requested allele at the 3' end on both strands", {
  L <- 22L
  for (st in c("+", "-")) {
    for (allele in c("A", "C")) {
      p <- make_as_primer(W, SNP_I, L, allele, strand = st)
      # "+": the primer carries the plus-strand base itself;
      # "-": the primer binds the plus strand, so its 3' base is complementary.
      expect_equal(nchar(p), L)
      expect_equal(
        substr(p, L, L),
        if (identical(st, "+")) {
          allele
        } else {
          as.character(Biostrings::complement(Biostrings::DNAString(allele)))
        }
      )
    }
  }
})

test_that("the DECIPHER target is the reverse complement of the primer (no stabiliser)", {
  L <- 21L
  for (st in c("+", "-")) {
    for (allele in c("A", "C", "G")) {
      p <- make_as_primer(W, SNP_I, L, allele, strand = st)
      tgt <- as.character(as_target(seg_with_allele(L, st, allele), st))
      expect_equal(as.character(Biostrings::reverseComplement(Biostrings::DNAString(p))), tgt)
    }
  }
})

test_that("stabilising mismatch replaces the intended base but keeps the 3' allele", {
  L <- 24L
  p0 <- make_as_primer(W, SNP_I, L, "A", strand = "+")
  b0 <- substr(p0, L - 3L + 1L, L - 3L + 1L)
  alt <- setdiff(c("A", "C", "G", "T"), b0)[1L]
  p1 <- make_as_primer(W, SNP_I, L, "A", strand = "+", stab_at = 3L, stab_base = alt)
  expect_equal(nchar(p1), L)
  expect_equal(substr(p1, L - 3L + 1L, L - 3L + 1L), alt)
  expect_equal(substr(p1, L, L), "A")
})

test_that("several ALTs are scored conservatively and carry a degenerate code", {
  skip_if_no_oligoarrayaux()
  wm <- Biostrings::DNAString(bseq(300L, 5L))
  v <- 150L
  ref <- as.character(Biostrings::subseq(wm, v, v))
  alts <- setdiff(c("A", "C", "G", "T"), ref)[1:2]
  multi <- design_as_primers(wm, v,
    ref = ref, alt = alts, lengths = 20:21, stab_ats = NA, strands = "+"
  )
  one <- design_as_primers(wm, v,
    ref = ref, alt = alts[1], lengths = 20:21, stab_ats = NA, strands = "+"
  )
  expect_gt(nrow(multi), 0L)
  # min-on / max-off over the ALTs is never more optimistic than a single ALT
  m <- multi$discrimination[order(multi$as_length)]
  o <- one$discrimination[order(one$as_length)]
  expect_true(all(m <= o + 1e-9))
  # the degenerate consensus carries the IUPAC code at the 3'-terminal block
  code <- iupac_code(alts)
  last <- nchar(multi$alt_degen)
  expect_true(all(substr(multi$alt_degen, last, last) == code))
  # a single ALT is unchanged: no degeneracy
  expect_equal(one$alt_degen, one$alt_primer)
})

test_that("a multi-allelic SNP is designed as one degenerate ALT primer", {
  skip_if_no_oligoarrayaux()
  withr::local_envvar(LANGUAGE = "en")
  wm <- Biostrings::DNAString(bseq(600L, 2L))
  genome <- Biostrings::DNAStringSet(c(Chr1 = wm))
  ref <- as.character(Biostrings::subseq(wm, 200L, 200L))
  others <- setdiff(c("A", "C", "G", "T"), ref)
  alts <- others[1:2]
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 200L, ref = ref, alt = alts[1], alt2 = alts[2]
  )
  a <- suppressWarnings(assemble_kasp_assay(genome, snp, flank = 150L, cores = 1L))
  expect_false(inherits(a, "kasp_doomed"))
  expect_equal(a$site$class[1], "degenerate")
  expect_equal(a$site$code[1], iupac_code(alts))
  # the ordered AS-alt oligo is the degenerate consensus (IUPAC 3' base)...
  last <- nchar(a$assay$as_alt_degen[1])
  expect_equal(substr(a$assay$as_alt_degen[1], last, last), a$site$code[1])
  # ...while the internal (concrete) AS-alt keeps a real base, so the
  # specificity/thermo machinery never sees an IUPAC letter
  last_c <- nchar(a$assay$as_alt[1])
  expect_true(substr(a$assay$as_alt[1], last_c, last_c) %in% alts)
  # export carries the degenerate oligo, not the concrete member
  tbl <- kasp_primer_table(a, list(symbol = "X"))
  expect_equal(tbl$seq_no_tail[tbl$oligo_id == "AS2_HEX"], a$assay$as_alt_degen[1])
  # the ALT primer is recorded as covering every ALT (a first-class criterion)
  chk <- kasp_design_checks(a)
  expect_true(chk$pass[chk$criterion == "ALT primer covers all ALTs"])
})
