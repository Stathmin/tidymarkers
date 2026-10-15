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
