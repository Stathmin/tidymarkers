# Package helpers: class predicate, all-pass predicate, env helper.

mk_assay <- function() {
  p <- kasp_primers()
  structure(list(
    snp = tibble::tibble(chrom = "Chr4", pos = 1L, ref = "G", alt = "A"),
    strand = "+",
    assay = tibble::tibble(
      amplicon = 90L, as_ref = p[["AS_ref"]], as_alt = p[["AS_alt"]],
      com = p[["COM"]], com_tm = 60, as_gc = 0.5, com_gc = 0.5,
      n_var = 0L, as_n_3p = 0L, as_n_var = 0L, as_dropout = 0, as_foot_dropout = 0, com_dropout = 0
    ),
    thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), gc = 0.5, hairpin_dg = -1),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3),
      b = rep(names(p), times = 3),
      dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
    ),
    amplicons = tibble::tibble(
      primer = c("AS_ref", "AS_alt"),
      n_products = c(1L, 0L), width = c(90L, NA_integer_)
    ),
    efficiency = tibble::tibble(discrimination = 0.9),
    specificity = tibble::tibble(primer = names(p), seed = "ACGT", hits = 1L)
  ), class = "marker_assay")
}

test_that("is_marker_assay recognises the class", {
  expect_true(is_marker_assay(mk_assay()))
  expect_false(is_marker_assay(list()))
  expect_false(is_marker_assay(NULL))
})

test_that("passes_all reflects the acceptance checks", {
  a <- mk_assay()
  expect_true(passes_all(a))
  a$assay$amplicon <- 500L # break the amplicon criterion
  expect_false(passes_all(a))
})

test_that("with_unafloddat sets and restores UNAFOLDDAT", {
  withr::local_envvar(LANGUAGE = "en")
  d <- withr::local_tempdir(pattern = "dat")
  old <- Sys.getenv("UNAFOLDDAT", unset = NA)
  seen <- with_unafloddat(d, Sys.getenv("UNAFOLDDAT"))
  expect_equal(seen, d)
  expect_equal(Sys.getenv("UNAFOLDDAT", unset = NA), old) # restored afterwards
  expect_error(with_unafloddat(file.path(d, "nope"), 1), "not found")
})

test_that("so_rank matches canonical SO terms, comma lists and VEP impact", {
  # literal expected ranks (not the same dict the function reads)
  expect_equal(so_rank("stop_gained"), 1L)
  expect_equal(so_rank("frameshift_variant"), 1L)
  expect_equal(so_rank("missense_variant"), 3L)
  expect_equal(so_rank("5_prime_UTR_variant"), 6L) # Ensembl spelling
  expect_equal(so_rank("synonymous_variant,stop_gained"), 1L) # most severe wins
  expect_equal(so_rank("weird_term"), SO_RANK_OTHER) # fallback, no error
  # impact ladder wins where present; NA impact falls back to the SO term
  expect_equal(
    so_rank(c("synonymous_variant", "missense_variant"), c("HIGH", NA)),
    c(1L, 3L)
  )
})

test_that("pick_target_snp prefers the locus with a unique 3'-seed", {
  rep30 <- "ACGTTGCAAGTCCTAGGCATCGATTCGGAT"
  pad <- "GCGCGCGCGCGCGCGCGCGC"
  uniq <- withr::with_seed(1, paste(sample(c("A", "C", "G", "T"), 80L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = paste0(pad, rep30, rep30, rep30, uniq)))
  transition <- c(A = "G", G = "A", C = "T", T = "C")
  a_pos <- 20L + 60L # end of the 2nd repeat: both strands' seeds repeat
  b_pos <- 20L + 90L + 40L # inside the unique stretch
  refa <- as.character(Biostrings::subseq(genome[["Chr1"]], a_pos, a_pos))
  refb <- as.character(Biostrings::subseq(genome[["Chr1"]], b_pos, b_pos))
  variants <- tibble::tibble(
    id = c("A", "B"), chrom = "Chr1", pos = c(a_pos, b_pos),
    ref_genome = list(c(refa, unname(transition[[refa]])), c(refb, unname(transition[[refb]]))),
    consequence = "missense_variant"
  )
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = 200L, strand = 1L)

  tgt <- pick_target_snp(variants, genome, gene)
  expect_equal(tgt$id, "B") # A is leftmost but its seed repeats
  expect_equal(tgt$seed_hits, 1L)
})

test_that("pick_target_snp is invariant to the seed-scan backend", {
  rep30 <- "ACGTTGCAAGTCCTAGGCATCGATTCGGAT"
  uniq <- withr::with_seed(1, paste(sample(c("A", "C", "G", "T"), 80L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = paste0(strrep("GC", 10L), rep30, rep30, uniq)))
  transition <- c(A = "G", G = "A", C = "T", T = "C")
  refs <- vapply(c(45L, 85L), function(p) as.character(Biostrings::subseq(genome[["Chr1"]], p, p)), character(1))
  variants <- tibble::tibble(
    id = c("x", "y"), chrom = "Chr1", pos = c(45L, 85L),
    ref_genome = list(
      c(refs[1], unname(transition[[refs[1]]])),
      c(refs[2], unname(transition[[refs[2]]]))
    ),
    consequence = "missense_variant"
  )
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = 200L, strand = 1L)
  a <- pick_target_snp(variants, genome, gene)
  b <- pick_target_snp(variants, genome, gene, specificity_candidates = 0L)
  expect_equal(a$id, b$id) # skip vs scan agree on the winner here
  skip_if_not_installed("BiocParallel")
  d <- pick_target_snp(variants, genome, gene, BPPARAM = BiocParallel::SerialParam())
  expect_equal(d, a)
})

test_that("pick_target_snp prefers the higher VEP impact", {
  seq <- strrep("A", 40L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = seq))
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L)
  variants <- tibble::tibble(
    id = c("lo", "hi"), chrom = "Chr1", pos = c(10L, 20L),
    ref_genome = list(c("A", "G"), c("A", "G")),
    consequence = c("missense_variant", "missense_variant"),
    impact = c("LOW", "HIGH")
  )
  tgt <- pick_target_snp(variants, genome, gene)
  expect_equal(tgt$id, "hi")
})
