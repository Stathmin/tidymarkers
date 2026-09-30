# Input readers: VCF (SNP/MNP/InDel/multiallelic) and arbitrary FASTA.

test_that("read_variants_vcf keeps SNP, InDel and multiallelic records", {
  skip_if_not_installed("vroom")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "4\t100\t.\tA\tG\t.\t.\t.",
    "4\t200\t.\tAT\tA\t.\t.\t.",          # deletion
    "4\t300\t.\tA\tG,T\t.\t.\t.",          # multi-allelic
    "4\t400\t.\tA\t<DEL>\t.\t.\t."         # symbolic -> dropped
  ), vcf)
  v <- read_variants_vcf(vcf)
  expect_setequal(v$pos, c(100L, 200L, 300L, 300L))
  expect_true(any(v$ref == "AT" & v$alt == "A"))
  expect_equal(sum(v$pos == 300L), 2L)
  expect_true(all(v$chrom == "Chr4"))          # names normalised
})

test_that("read_variants_vcf filters by region", {
  skip_if_not_installed("vroom")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t10\t.\tA\tG\t.\t.\t.",
    "1\t50\t.\tA\tC\t.\t.\t.",
    "1\t90\t.\tA\tT\t.\t.\t."
  ), vcf)
  expect_setequal(read_variants_vcf(vcf, chrom = "Chr1", start = 20L, end = 60L)$pos, 50L)
})

test_that("normalise_allele folds IUPAC codes and case to ACGTN", {
  expect_equal(normalise_allele(c("a", "g")), c("A", "G"))
  expect_equal(normalise_allele(c("R", "y", "N")), c("N", "N", "N"))
  expect_equal(normalise_allele("ACRT"), "ACNT")          # ambiguity inside a longer allele
})

test_that("read_variants_vcf normalises IUPAC/lowercase instead of dropping", {
  skip_if_not_installed("vroom")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "4\t500\t.\tA\tR\t.\t.\t.",          # IUPAC ambiguity -> N
    "4\t600\t.\ta\tg\t.\t.\t."           # lowercase -> uppercase
  ), vcf)
  v <- read_variants_vcf(vcf)
  expect_equal(v$alt[v$pos == 500L], "N")
  expect_equal(v$ref[v$pos == 600L], "A")
  expect_equal(v$alt[v$pos == 600L], "G")
})

test_that("the vroom and VariantAnnotation backends agree", {
  skip_if_not_installed("vroom")
  # VariantAnnotation's namespace load emits a benign import-replacement warning
  if (!suppressWarnings(suppressMessages(requireNamespace("VariantAnnotation", quietly = TRUE)))) {
    testthat::skip("VariantAnnotation not installed")
  }
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "4\t100\t.\tA\tG\t.\t.\t.",
    "4\t200\t.\tAT\tA\t.\t.\t.",          # deletion
    "4\t300\t.\tA\tG,T\t.\t.\t.",          # multi-allelic
    "4\t400\t.\tA\t<DEL>\t.\t.\t."         # symbolic -> dropped
  ), vcf)
  a <- read_variants_vcf(vcf, backend = "vroom")
  b <- suppressWarnings(read_variants_vcf(vcf, backend = "VariantAnnotation"))
  key <- function(x) sort(paste(x$chrom, x$pos, x$ref, x$alt))
  expect_equal(key(a), key(b))
  expect_setequal(a$pos, c(100L, 200L, 300L, 300L))
})

test_that("the VariantAnnotation backend filters by region too", {
  skip_if_not_installed("VariantAnnotation")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
               "1\t10\t.\tA\tG\t.\t.\t.",
               "1\t50\t.\tA\tC\t.\t.\t.",
               "1\t90\t.\tA\tT\t.\t.\t."), vcf)
  v <- suppressWarnings(
    read_variants_vcf(vcf, chrom = "Chr1", start = 20L, end = 60L,
                      backend = "VariantAnnotation"))
  expect_setequal(v$pos, 50L)
  expect_true(all(v$chrom == "Chr1"))
})

test_that("read_genome_fasta normalises and filters sequence names", {
  fa <- tempfile(fileext = ".fa")
  writeLines(c(">1 dna:chromosome chromosome:TAIR10:1:1:100:1 REF", "ACGTACGT",
               ">Mt dna:chromosome", "ACGT", ">scaffold_9 x", "TTTT"), fa)
  g <- read_genome_fasta(fa, chrom_map = ATH_CHROM_MAP, use_index = FALSE)
  expect_setequal(names(g), c("Chr1", "ChrM"))     # scaffold dropped
  expect_equal(as.character(g[["Chr1"]]), "ACGTACGT")
  expect_setequal(names(read_genome_fasta(fa, seqids = "1", use_index = FALSE)), "1")
})

test_that("read_genome_fasta reads only requested ids via the FASTA index", {
  skip_if_not_installed("Rsamtools")
  fa <- tempfile(fileext = ".fa")
  writeLines(c(">chr1", "ACGTACGTACGT", ">chr2", "TTTTGGGGCCCC", ">chr3", "GGGGCCCCAAAA"), fa)
  g <- read_genome_fasta(fa, seqids = c("chr2", "chr3"))
  expect_setequal(names(g), c("chr2", "chr3"))     # chr1 not read
  expect_equal(as.character(g[["chr2"]]), "TTTTGGGGCCCC")
  expect_true(file.exists(paste0(fa, ".fai")))       # index created
})

test_that("read_genome_fasta map+index keeps only mapped, requested chroms", {
  skip_if_not_installed("Rsamtools")
  fa <- tempfile(fileext = ".fa")
  writeLines(c(">1 desc", "ACGTACGT", ">2 desc", "TTTTGGGG"), fa)
  g <- read_genome_fasta(fa, chrom_map = ATH_CHROM_MAP, seqids = "Chr2")
  expect_equal(names(g), "Chr2")
  expect_equal(as.character(g[["Chr2"]]), "TTTTGGGG")
})
