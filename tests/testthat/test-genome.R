# Indexed genome source: regions on demand, same results as an in-memory set.

make_fa <- function(fa) {
  writeLines(c(
    ">1 desc", "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT",
    ">2 desc", "TTTTGGGGCCCCAAAATTTT"
  ), fa)
  fa
}

test_that("a marker_genome reads the same regions as an in-memory DNAStringSet", {
  skip_if_not_installed("Rsamtools")
  fa <- make_fa(withr::local_tempfile(fileext = ".fa"))
  gs <- marker_genome(fa, chrom_map = ATH_CHROM_MAP)
  expect_s3_class(gs, "marker_genome")
  expect_setequal(genome_chroms(gs), c("Chr1", "Chr2"))
  expect_equal(genome_length(gs, "Chr1"), 40L)
  expect_equal(as.character(genome_seq(gs, "Chr1", 5L, 12L)), "ACGTACGT")

  mem <- read_genome_fasta(fa, chrom_map = ATH_CHROM_MAP, use_index = FALSE)
  expect_equal(as.character(genome_seq(gs, "Chr1")), as.character(mem[["Chr1"]]))
  expect_equal(
    as.character(genome_seq(gs, "Chr2", 3L, 6L)),
    as.character(genome_seq(mem, "Chr2", 3L, 6L))
  )
})

test_that("a whole-sequence read is cached on the object", {
  skip_if_not_installed("Rsamtools")
  fa <- make_fa(withr::local_tempfile(fileext = ".fa"))
  gs <- marker_genome(fa, chrom_map = ATH_CHROM_MAP)
  first <- as.character(genome_seq(gs, "Chr1"))

  # a second whole read must not touch the FASTA
  gs$fa <- "does-not-exist.fa"
  expect_equal(as.character(genome_seq(gs, "Chr1")), first)
  # a partial read still does, so the cache is not a blanket lie
  expect_error(genome_seq(gs, "Chr1", 5L, 12L))

  withr::local_options(tidymarkers.genome_cache = FALSE)
  gs_off <- marker_genome(fa, chrom_map = ATH_CHROM_MAP)
  genome_seq(gs_off, "Chr1")
  gs_off$fa <- "does-not-exist.fa"
  expect_error(genome_seq(gs_off, "Chr1"))
})

test_that("design helpers accept a marker_genome (no whole genome in RAM)", {
  skip_if_not_installed("Rsamtools")
  fa <- make_fa(withr::local_tempfile(fileext = ".fa"))
  gs <- marker_genome(fa, chrom_map = ATH_CHROM_MAP)
  tmpl <- snp_template(gs, "Chr1", 20L, flank = 10L)
  expect_equal(length(tmpl$seq), 21L)
  expect_equal(tmpl$snp_i, 11L)
  hint <- strand_gc_hint(gs, "Chr1", 20L, flank = 5L)
  expect_equal(nrow(hint), 2L)
})

test_that("snp_template refuses a fractional flank before the window is cast", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(100L, 1L)))
  # as.integer(10.9) is 10: the window would be one base narrower
  expect_error(snp_template(genome, "Chr1", 20L, flank = 10.9), "positive whole number of base pairs")
  expect_equal(length(snp_template(genome, "Chr1", 20L, flank = 10L)$seq), 21L)
})

test_that("in_silico_pcr works through an indexed marker_genome", {
  skip_if_not_installed("Rsamtools")
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">chr1", paste0(strrep("A", 30L), "ACGTACGTACGTACGTACGT", strrep("T", 30L))), fa)
  gs <- marker_genome(fa)
  pr <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", gs)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$width, 20L)
})

test_that("a BSgenome-backed genome has no FASTA, so the index builders say so (№8)", {
  # built by hand: the point is the missing FASTA, not the BSgenome package
  bs <- structure(list(path = NULL, bsgenome = "any"), class = "marker_genome")
  for (build in list(marker_blast_index, marker_bowtie_index, marker_fm_index)) {
    expect_error(build(bs), "BSgenome-backed")
    expect_error(build(bs), "marker_genome")
  }
  # the escape hatch names the argument that builder actually takes
  expect_error(marker_blast_index(bs), "blast_db =")
  expect_error(marker_bowtie_index(bs), "bowtie_index =")
  expect_error(marker_fm_index(bs), "fm_index =")

  # a FASTA-backed genome is indexed exactly as before
  skip_if_not_installed("Rsamtools")
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">chr1", strrep("ACGT", 25L)), fa)
  gs <- marker_genome(fa)
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  skip_if_not_installed("rBLAST")
  expect_identical(marker_blast_index(gs), gs$path)
})

test_that("genome_seq refuses a coordinate it would have truncated", {
  skip_if_not_installed("Rsamtools")
  fa <- make_fa(withr::local_tempfile(fileext = ".fa"))
  gs <- marker_genome(fa, chrom_map = ATH_CHROM_MAP)
  mem <- read_genome_fasta(fa, chrom_map = ATH_CHROM_MAP, use_index = FALSE)

  for (g in list(gs, mem)) {
    # 2.9 used to be read as base 2, so the caller got a span nobody asked for
    expect_error(genome_seq(g, "Chr1", 2.9, 6L), "positive whole number")
    expect_error(genome_seq(g, "Chr1", 2L, 6.9), "positive whole number")
    expect_error(genome_seq(g, "Chr1", 0L, 6L), "positive whole number")
    expect_error(genome_seq(g, "Chr1", NA, 6L), "positive whole number")
    expect_error(genome_seq(g, "Chr1", 1L, Inf), "positive whole number")
    expect_error(genome_seq(g, "Chr1", c(1L, 2L), 6L), "positive whole number")
    expect_error(genome_seq(g, "Chr1", 6L, 2L), "must not be past")
    # the same region the caller named still reads
    expect_equal(as.character(genome_seq(g, "Chr1", 2L, 6L)), "CGTAC")
  }
})

test_that("a marker_genome can be backed by a BSgenome (getSeq)", {
  skip_if_not_full()
  skip_if_not_installed_noload("BSgenome.Athaliana.TAIR.TAIR9")
  gs <- marker_genome(bsgenome = "BSgenome.Athaliana.TAIR.TAIR9")
  expect_s3_class(gs, "marker_genome")
  expect_true("Chr4" %in% genome_chroms(gs))
  expect_equal(genome_length(gs, "Chr4"), 18585056L)
  # the FRI stop variant at Chr4:269745 is a G on TAIR9
  expect_equal(as.character(genome_seq(gs, "Chr4", 269745L, 269745L)), "G")
})

# One planted primer and one absent primer, for both indexed backends.  The
# absent one is the regression guard of D7: an unaligned record must not be
# counted, so a backend that reports it as one hit fails here.
make_planted <- function(fa) {
  present <- "ACGTACGTACGTACGTACGT"
  writeLines(c(">chr1", paste0(strrep("A", 300L), present, strrep("T", 300L))), fa)
  list(
    primers = c(present = present, absent = "GGGGGGGGGGGGGGGGGGGG"),
    present = present
  )
}

test_that("specificity_bowtie finds a planted primer and none that is absent", {
  skip_if_not_installed("Rbowtie")
  skip_if_not_installed("Rsamtools")
  fa <- withr::local_tempfile(fileext = ".fa")
  planted <- make_planted(fa)
  gs <- marker_genome(fa)
  idx <- marker_bowtie_index(gs, basename = withr::local_tempfile(pattern = "bt"))
  hits <- specificity_bowtie(planted$primers, idx)
  # bowtie reports both strands, so one genomic site is two alignments
  expect_equal(hits$hits[hits$primer == "present"], 2L)
  expect_gte(hits$hits[hits$primer == "present"], 1L)
  expect_equal(hits$hits[hits$primer == "absent"], 0L)
})

test_that("specificity_fm finds a planted primer and none that is absent", {
  skip_if_not_installed("Rsubread")
  skip_if_not_installed("Rsamtools")
  fa <- withr::local_tempfile(fileext = ".fa")
  planted <- make_planted(fa)
  gs <- marker_genome(fa)
  idx <- marker_fm_index(gs, basename = withr::local_tempfile(pattern = "fm"))
  hits <- specificity_fm(planted$primers, idx)
  # Subread maps the forward strand by default, so the one planted site is one hit
  expect_equal(hits$hits[hits$primer == "present"], 1L)
  expect_gte(hits$hits[hits$primer == "present"], 1L)
  # D7: an absent primer aligns nowhere, so an unaligned record is not a hit
  expect_equal(hits$hits[hits$primer == "absent"], 0L)
})
