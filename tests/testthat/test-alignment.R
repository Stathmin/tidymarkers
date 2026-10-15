# A gene's transcripts as a haplotype set: FASTA -> alignment -> haplotype table.

# The shipped SAT2 seqdump (11 Triticum accessions), aligned once for the file.
sat2_path <- system.file("extdata", "sat2_transcripts.fa", package = "tidykasp")

sat2_hp <- local({
  if (!nzchar(sat2_path)) {
    return(NULL)
  }
  haplotypes_from_alignment(read_transcripts(sat2_path), chrom = "SAT2")
})

test_that("read_transcripts reads one sequence per accession and understands ranges", {
  skip_if(!nzchar(sat2_path), "bundled SAT2 FASTA not found")
  seqs <- read_transcripts(sat2_path)
  expect_s4_class(seqs, "DNAStringSet")
  expect_length(seqs, 11L)
  expect_identical(names(seqs)[1L], "JX896657.1")
  # the coordinate-range records of one accession are joined, not kept apart
  expect_false(any(grepl(":", names(seqs), fixed = TRUE)))
  # keeping them apart gives back every record
  expect_length(read_transcripts(sat2_path, merge_accessions = FALSE), 18L)
})

test_that("haplotypes_from_alignment returns a reference genome and a tidy table", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  expect_length(sat2_hp$reference, 1209L)
  expect_s4_class(sat2_hp$genome, "DNAStringSet")
  expect_identical(names(sat2_hp$genome), "SAT2")
  expect_s4_class(sat2_hp$alignment, "DNAMultipleAlignment")
  expect_equal(nrow(sat2_hp$alignment), 11L)

  h <- sat2_hp$haplotypes
  expect_true(all(c("haplotype", "freq", "chrom", "pos", "ref", "alt") %in% names(h)))
  expect_gt(nrow(h), 50L)
  expect_true(all(h$chrom == "SAT2"))
  expect_true(all(h$pos >= 1L & h$pos <= length(sat2_hp$reference)))
  expect_setequal(unique(h$freq), 1 / 11) # equal frequencies by default
  # both kinds of variant are present (at least one indel)
  expect_true(any(nchar(h$ref) != nchar(h$alt)))
  expect_true(any(nchar(h$ref) == 1L & nchar(h$alt) == 1L))
})

test_that("haplotypes_from_alignment works on a tiny inline example and picks the reference", {
  seqs <- c(
    ref = "ACGTACGTACGTACGTACGT", hap1 = "ACGTACGAACGTACGTACGT",
    hap2 = "ACGTACGTACGTACGTACGA"
  )
  hp <- haplotypes_from_alignment(seqs, chrom = "GENE")
  expect_identical(hp$reference_name, "ref")
  expect_equal(as.character(hp$reference), seqs[["ref"]])
  expect_setequal(unique(hp$haplotypes$haplotype), c("hap1", "hap2"))
  # named reference selection
  hp2 <- haplotypes_from_alignment(seqs, reference = "hap1", chrom = "GENE")
  expect_identical(hp2$reference_name, "hap1")
})

test_that("a haplotype set designs with the haplotype engine", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  skip_if_no_oligoarrayaux()
  h <- sat2_hp$haplotypes
  # a substitution comfortably inside the reference (room for a 150 bp flank)
  v <- h[nchar(h$ref) == 1L & nchar(h$alt) == 1L & h$pos > 200L & h$pos < 900L, , drop = FALSE]
  skip_if(nrow(v) == 0L, "no interior substitution in the SAT2 haplotypes")
  snp <- tibble::tibble(chrom = "SAT2", pos = v$pos[1L], ref = v$ref[1L], alt = v$alt[1L])
  a <- suppressWarnings(assemble_kasp_assay(
    sat2_hp$genome, snp,
    haplotypes = sat2_hp$haplotypes, flank = 150L, cores = 1L
  ))
  # it designs or dooms - never errors - and carries the haplotype verdicts when it designs
  if (inherits(a, "kasp_assay")) {
    expect_true(all(c("haplotypes", "haplotype_hits") %in% names(a)))
    expect_equal(nrow(a$haplotype_hits), 11L)
  } else {
    expect_s3_class(a, "kasp_doomed")
  }
})

test_that("the plan layer carries haplotype context", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  plan <- kasp_plan(sat2_hp$genome) +
    kasp_targets(chrom = "SAT2", pos = 356L, ref = "T", alt = "C") +
    kasp_context(haplotypes = sat2_hp$haplotypes)
  expect_s3_class(plan, "kasp_plan")
  expect_equal(nrow(plan$targets), 1L)
  expect_equal(nrow(plan$haplotypes), nrow(sat2_hp$haplotypes))
})

test_that("variants carry pop_freq, and a MAF table overrides it", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  v <- sat2_hp$variants
  expect_true(all(c("chrom", "pos", "ref", "alt", "pop_freq", "consequence") %in% names(v)))
  expect_true(all(v$pop_freq > 0 & v$pop_freq <= 1))
  # the reference haplotype carries nothing, so equal 1/11 freqs give k/11 weights
  expect_true(any(abs(v$pop_freq - round(v$pop_freq * 11) / 11) < 1e-9))
  # an explicit MAF wins for the matching variant
  target <- v[!duplicated(v[, c("pos", "ref", "alt")]), ][1L, ]
  hp2 <- haplotypes_from_alignment(
    read_transcripts(sat2_path),
    chrom = "SAT2",
    maf = tibble::tibble(pos = target$pos, ref = target$ref, alt = target$alt, maf = 0.02)
  )
  hit <- hp2$variants[hp2$variants$pos == target$pos & hp2$variants$alt == target$alt, ]
  expect_equal(hit$pop_freq, 0.02, tolerance = 1e-9)
})

test_that("haplotype frequencies may be given as a table", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  nms <- names(read_transcripts(sat2_path))
  tab <- tibble::tibble(haplotype = nms, freq = seq_along(nms) / sum(seq_along(nms)))
  hp3 <- haplotypes_from_alignment(read_transcripts(sat2_path), chrom = "SAT2", freq = tab)
  expect_true(all(unique(hp3$haplotypes$freq) %in% tab$freq))
})

test_that("plot_alignment returns a ggplot", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot_alignment(sat2_hp), "ggplot")
  expect_s3_class(plot_alignment(sat2_hp, type = "all"), "ggplot")
})

test_that("plot_alignment draws a feature track", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  skip_if_not_installed("ggplot2")
  feats <- tibble::tibble(
    start = c(1L, 54L), end = c(1209L, 1067L), strand = "+", name = c("SAT2", "CDS")
  )
  expect_s3_class(plot_alignment(sat2_hp, features = feats), "ggplot")
  expect_s3_class(plot_alignment(sat2_hp, type = "all", features = feats), "ggplot")
})

test_that("an AS primer's genomic interval matches the genome (no off-by-one)", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  skip_if_no_oligoarrayaux()
  snp <- tibble::tibble(chrom = "SAT2", pos = 431L, ref = "A", alt = "G")
  a <- suppressWarnings(assemble_kasp_assay(sat2_hp$genome, snp, flank = 150L, cores = 1L))
  skip_if(!inherits(a, "kasp_assay"), "SAT2:431 did not design here")
  aa <- a$assay[1L, ]
  ref <- as.character(sat2_hp$reference)
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
  # the AS primer spells exactly the genome it covers (its strand), and its 3' base
  # sits on the variant - the interval must not be shifted by one
  seg <- substr(ref, aa$as_start, aa$as_end)
  expect_equal(if (identical(aa$strand, "+")) seg else rc(seg), aa$as_ref)
  expect_equal(nchar(aa$as_ref), aa$as_end - aa$as_start + 1L)
  three <- if (identical(aa$strand, "+")) aa$as_end else aa$as_start
  expect_equal(three, snp$pos)
})

test_that("kasp_chunk slices a haplotype set to a local, window-local chunk", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  ch <- kasp_chunk(sat2_hp, from = 281L, to = 581L)
  expect_s3_class(ch, "kasp_chunk")
  expect_equal(ch$reference_length, 301L)
  expect_identical(names(chunk_genome(ch)), "SAT2")
  expect_equal(
    as.character(chunk_reference(ch)),
    as.character(Biostrings::subseq(sat2_hp$genome[[1L]], 281L, 581L))
  )
  v <- chunk_variants(ch)
  expect_true(all(v$pos >= 1L & v$pos <= 301L)) # window-local
  expect_s4_class(chunk_alignment(ch), "DNAMultipleAlignment")
  expect_output(print(ch), "kasp_chunk")
})

test_that("a chunk round-trips between the representations", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  ch <- kasp_chunk(sat2_hp, from = 281L, to = 581L)
  # alignment -> (variants, haplotypes): both derived, and consistent
  h <- chunk_haplotypes(ch)
  expect_gte(nrow(h), 1L)
  expect_equal(nrow(chunk_variants(ch)), nrow(variants_from_haplotypes(h)))
  # haplotypes -> a variants-only chunk: its haplotypes are honestly absent
  ch2 <- kasp_chunk(as.data.frame(chunk_variants(ch)), reference = chunk_reference(ch), chrom = "SAT2")
  expect_null(chunk_haplotypes(ch2))
  expect_equal(nrow(chunk_variants(ch2)), nrow(chunk_variants(ch)))
})

test_that("chunk_weight applies MAF-style weights", {
  skip_if(is.null(sat2_hp), "bundled SAT2 FASTA not found")
  ch <- kasp_chunk(sat2_hp, from = 281L, to = 581L)
  t1 <- chunk_variants(ch)[1L, ]
  ch2 <- chunk_weight(ch, tibble::tibble(pos = t1$pos, ref = t1$ref, alt = t1$alt, maf = 0.02))
  got <- chunk_variants(ch2)
  expect_equal(got$pop_freq[got$pos == t1$pos & got$alt == t1$alt][[1L]], 0.02)
})
