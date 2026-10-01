# Indexed genome source: regions on demand, same results as an in-memory set.

make_fa <- function() {
  fa <- tempfile(fileext = ".fa")
  writeLines(c(
    ">1 desc", "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT",
    ">2 desc", "TTTTGGGGCCCCAAAATTTT"
  ), fa)
  fa
}

test_that("a kasp_genome reads the same regions as an in-memory DNAStringSet", {
  skip_if_not_installed("Rsamtools")
  fa <- make_fa()
  gs <- kasp_genome(fa, chrom_map = ATH_CHROM_MAP)
  expect_s3_class(gs, "kasp_genome")
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

test_that("design helpers accept a kasp_genome (no whole genome in RAM)", {
  skip_if_not_installed("Rsamtools")
  fa <- make_fa()
  gs <- kasp_genome(fa, chrom_map = ATH_CHROM_MAP)
  tmpl <- snp_template(gs, "Chr1", 20L, flank = 10L)
  expect_equal(length(tmpl$seq), 21L)
  expect_equal(tmpl$snp_i, 11L)
  hint <- strand_gc_hint(gs, "Chr1", 20L, flank = 5L)
  expect_equal(nrow(hint), 2L)
})

test_that("in_silico_pcr works through an indexed kasp_genome", {
  skip_if_not_installed("Rsamtools")
  fa <- tempfile(fileext = ".fa")
  writeLines(c(">chr1", paste0(strrep("A", 30L), "ACGTACGTACGTACGTACGT", strrep("T", 30L))), fa)
  gs <- kasp_genome(fa)
  pr <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", gs)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$width, 20L)
})

test_that("a kasp_genome can be backed by a BSgenome (getSeq)", {
  skip_if_not_installed("BSgenome.Athaliana.TAIR.TAIR9")
  gs <- kasp_genome(bsgenome = "BSgenome.Athaliana.TAIR.TAIR9")
  expect_s3_class(gs, "kasp_genome")
  expect_true("Chr4" %in% genome_chroms(gs))
  expect_equal(genome_length(gs, "Chr4"), 18585056L)
  # the FRI stop variant at Chr4:269745 is a G on TAIR9
  expect_equal(as.character(genome_seq(gs, "Chr4", 269745L, 269745L)), "G")
})
