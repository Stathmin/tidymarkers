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

test_that("read_genome_fasta normalises and filters sequence names", {
  fa <- tempfile(fileext = ".fa")
  writeLines(c(">1 dna:chromosome chromosome:TAIR10:1:1:100:1 REF", "ACGTACGT",
               ">Mt dna:chromosome", "ACGT", ">scaffold_9 x", "TTTT"), fa)
  g <- read_genome_fasta(fa, chrom_map = ATH_CHROM_MAP)
  expect_setequal(names(g), c("Chr1", "ChrM"))     # scaffold dropped
  expect_equal(as.character(g[["Chr1"]]), "ACGTACGT")
  expect_setequal(names(read_genome_fasta(fa, seqids = "1")), "1")
})
