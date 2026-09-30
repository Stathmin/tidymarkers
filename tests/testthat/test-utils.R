# Package helpers: class predicate, all-pass predicate, env helper.

mk_assay <- function() {
  p <- c(AS_ref = "TGACTGCAAGTCCGATCAGTA",
         AS_alt = "TGACTGCAAGTCCGATCAGTT",
         COM    = "GCTTAGCCATGGTCAATGCTA")
  structure(list(
    snp = tibble::tibble(chrom = "Chr4", pos = 1L, ref = "G", alt = "A"),
    strand = "+",
    assay = tibble::tibble(amplicon = 90L, as_ref = p[["AS_ref"]], as_alt = p[["AS_alt"]],
                           com = p[["COM"]], com_tm = 60, as_gc = 0.5, com_gc = 0.5),
    thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), gc = 0.5, hairpin_dg = -1),
    complementarity = tibble::tibble(a = rep(names(p), each = 3),
                                     b = rep(names(p), times = 3),
                                     dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)),
    amplicons = tibble::tibble(primer = c("AS_ref", "AS_alt"),
                               n_products = c(1L, 0L), width = c(90L, NA_integer_)),
    efficiency = tibble::tibble(discrimination = 0.9),
    specificity = tibble::tibble(primer = names(p), seed = "ACGT", hits = 1L)
  ), class = "kasp_assay")
}

test_that("is_kasp_assay recognises the class", {
  expect_true(is_kasp_assay(mk_assay()))
  expect_false(is_kasp_assay(list()))
  expect_false(is_kasp_assay(NULL))
})

test_that("passes_all reflects the acceptance checks", {
  a <- mk_assay()
  expect_true(passes_all(a))
  a$assay$amplicon <- 500L          # break the amplicon criterion
  expect_false(passes_all(a))
})

test_that("with_unafloddat sets and restores UNAFOLDDAT", {
  d <- tempfile("dat"); dir.create(d)
  old <- Sys.getenv("UNAFOLDDAT", unset = NA)
  seen <- with_unafloddat(d, Sys.getenv("UNAFOLDDAT"))
  expect_equal(seen, d)
  expect_equal(Sys.getenv("UNAFOLDDAT", unset = NA), old)   # restored afterwards
  expect_error(with_unafloddat(file.path(d, "nope"), 1), "not found")
})

test_that("so_rank matches canonical SO terms, comma lists and VEP impact", {
  # literal expected ranks (not the same dict the function reads)
  expect_equal(so_rank("stop_gained"), 1L)
  expect_equal(so_rank("frameshift_variant"), 1L)
  expect_equal(so_rank("missense_variant"), 3L)
  expect_equal(so_rank("5_prime_UTR_variant"), 6L)          # Ensembl spelling
  expect_equal(so_rank("synonymous_variant,stop_gained"), 1L)  # most severe wins
  expect_equal(so_rank("weird_term"), SO_RANK_OTHER)        # fallback, no error
  # impact ladder wins where present; NA impact falls back to the SO term
  expect_equal(so_rank(c("synonymous_variant", "missense_variant"), c("HIGH", NA)),
               c(1L, 3L))
})

test_that("seed_hits_at counts a repeated seed, and is NA off the start", {
  rep30 <- "ACGTTGCAAGTCCTAGGCATCGATTCGGAT"
  genome <- Biostrings::DNAStringSet(c(Chr1 = paste0("GCGCGCGCGCGCGCGCGCGC", rep30, rep30, rep30)))
  expect_gt(seed_hits_at(genome, "Chr1", 80L, 1L, "+"), 1L)   # last 15 nt of rep30 -> 3 copies
  expect_gt(seed_hits_at(genome, "Chr1", 80L, 1L, "-"), 1L)
  expect_true(is.na(seed_hits_at(genome, "Chr1", 5L, 1L, "+")))  # runs off the start
})

test_that("pick_target_snp prefers the locus with a unique 3'-seed", {
  rep30 <- "ACGTTGCAAGTCCTAGGCATCGATTCGGAT"
  pad   <- "GCGCGCGCGCGCGCGCGCGC"
  set.seed(1)
  uniq  <- paste(sample(c("A", "C", "G", "T"), 80L, TRUE), collapse = "")
  genome <- Biostrings::DNAStringSet(c(Chr1 = paste0(pad, rep30, rep30, rep30, uniq)))
  transition <- c(A = "G", G = "A", C = "T", T = "C")
  a_pos <- 20L + 60L              # end of the 2nd repeat: both strands' seeds repeat
  b_pos <- 20L + 90L + 40L        # inside the unique stretch
  refa <- as.character(Biostrings::subseq(genome[["Chr1"]], a_pos, a_pos))
  refb <- as.character(Biostrings::subseq(genome[["Chr1"]], b_pos, b_pos))
  variants <- tibble::tibble(
    id = c("A", "B"), chrom = "Chr1", pos = c(a_pos, b_pos),
    ref_genome = list(c(refa, unname(transition[[refa]])), c(refb, unname(transition[[refb]]))),
    consequence = "missense_variant")
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = 200L, strand = 1L)

  tgt <- pick_target_snp(variants, genome, gene)
  expect_equal(tgt$id, "B")        # A is leftmost but its seed repeats
  expect_equal(tgt$seed_hits, 1L)
})

test_that("pick_target_snp prefers the higher VEP impact", {
  seq <- strrep("A", 40L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = seq))
  gene <- list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L)
  variants <- tibble::tibble(
    id = c("lo", "hi"), chrom = "Chr1", pos = c(10L, 20L),
    ref_genome = list(c("A", "G"), c("A", "G")),
    consequence = c("missense_variant", "missense_variant"),
    impact = c("LOW", "HIGH"))
  tgt <- pick_target_snp(variants, genome, gene)
  expect_equal(tgt$id, "hi")
})
