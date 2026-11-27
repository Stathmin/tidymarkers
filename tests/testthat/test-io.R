# Variant interchange (roadmap Z10.3): variants tibble <-> GRanges / VCF.

.interchange_variants <- function() {
  tibble::tibble(
    chrom = c("Chr4", "Chr4", "Chr4", "Chr1"),
    pos = c(100L, 200L, 300L, 10L),
    ref = c("A", "AT", "A", "G"), # SNP, deletion, multi-allelic, insertion
    alt = c("G", "A", "T", "GAT"),
    depth = c(30L, 12L, NA, 8L)
  )
}

test_that("marker_as_granges round-trips the variants tibble and keeps depth", {
  skip_if_not_installed("GenomicRanges")
  v <- .interchange_variants()
  gr <- marker_as_granges(v)
  expect_s4_class(gr, "GRanges")
  expect_equal(as.character(gr), c("Chr4:100", "Chr4:200-201", "Chr4:300", "Chr1:10"))
  expect_equal(gr$ref, v$ref) # REF stays a string, InDel keeps its anchor base
  expect_equal(gr$alt, v$alt)
  expect_equal(gr$depth, v$depth)
  expect_equal(marker_from_granges(gr), v)
})

test_that("marker_as_granges carries depth from read_variant_depth as an mcol", {
  skip_if_not_installed("GenomicRanges")
  v <- .interchange_variants()[1:2, c("chrom", "pos", "ref", "alt")]
  dep <- tibble::tibble(depth = c(30L, 12L), genotype = c("0/1", "1/1"))
  back <- marker_from_granges(marker_as_granges(v, extra = dep))
  expect_equal(back$depth, c(30L, 12L))
  expect_equal(back$genotype, c("0/1", "1/1"))
  expect_named(back, c("chrom", "pos", "ref", "alt", "depth", "genotype"))
})

test_that("write_variants_vcf round-trips SNP, InDel and multi-allelic through VCF", {
  skip_if_not_installed("vroom")
  v <- .interchange_variants()
  f <- withr::local_tempfile(fileext = ".vcf")
  expect_invisible(write_variants_vcf(v, f))
  # one ALT per record, so compare on the allele key; depth comes from INFO/DP
  got <- merge(
    read_variants_vcf(f),
    read_variant_depth(f)[, c("chrom", "pos", "ref", "alt", "depth")],
    by = c("chrom", "pos", "ref", "alt")
  )
  key <- function(x) paste(x$chrom, x$pos, x$ref, x$alt)
  expect_setequal(key(got), key(v))
  got <- got[match(key(v), key(got)), ]
  expect_equal(got$depth, v$depth) # DP survives; NA stays NA
  # the caller's own temp file, not a file that disappears with a frame
  expect_true(file.exists(f))
  expect_match(readLines(f)[1L], "^##fileformat=VCFv4.2$")
})

# Input readers: VCF (SNP/MNP/InDel/multiallelic) and arbitrary FASTA.

test_that("a reader names a path that is not one existing file", {
  readers <- list(read_variants_vcf, read_variant_genotypes, read_genome_fasta)
  for (reader in readers) {
    # `file_exists(NULL)` is a base R error, and a directory passes it: both
    # must be the reader's own named refusal, not "argument is not a character
    # vector" or an empty tibble with a warning.
    expect_error(reader(NULL), "must be one", class = "rlang_error")
    expect_error(reader(character()), "must be one", class = "rlang_error")
    expect_error(reader(tempdir()), "names a directory", class = "rlang_error")
  }
})

test_that("read_variants_vcf keeps SNP, InDel and multiallelic records", {
  skip_if_not_installed("vroom")
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "4\t100\t.\tA\tG\t.\t.\t.",
    "4\t200\t.\tAT\tA\t.\t.\t.", # deletion
    "4\t300\t.\tA\tG,T\t.\t.\t.", # multi-allelic
    "4\t400\t.\tA\t<DEL>\t.\t.\t." # symbolic alleles are dropped
  ), vcf)
  v <- read_variants_vcf(vcf)
  expect_setequal(v$pos, c(100L, 200L, 300L, 300L))
  expect_true(any(v$ref == "AT" & v$alt == "A"))
  expect_equal(sum(v$pos == 300L), 2L)
  expect_true(all(v$chrom == "Chr4")) # names normalised
})

test_that("read_variants_vcf names the symbolic alleles it drops (1602)", {
  skip_if_not_installed("vroom")
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "4\t100\t.\tA\tG\t.\t.\t.",
    "4\t400\t.\tA\t<DEL>\t.\t.\t.",
    "4\t500\t.\tA\t*\t.\t.\t."
  ), vcf)
  expect_message(v <- read_variants_vcf(vcf), "2 VCF records")
  expect_setequal(v$pos, 100L) # the designable record survives untouched

  clean <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "4\t100\t.\tA\tG\t.\t.\t."
  ), clean)
  expect_message(read_variants_vcf(clean), NA) # quiet when nothing is dropped
})

test_that("read_variants_vcf handles empty, mixed and long-header VCFs (Y5)", {
  skip_if_not_installed("vroom")
  hdr <- c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  # a header with no records returns a typed zero-row frame, not an error
  empty <- withr::local_tempfile(fileext = ".vcf")
  writeLines(hdr, empty)
  v <- read_variants_vcf(empty)
  expect_equal(nrow(v), 0L)
  expect_named(v, c("chrom", "pos", "ref", "alt"))
  # `.` (absent) and `*` (spanning) alleles are dropped, a real SNV is kept
  mixed <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    hdr,
    "1\t250\t.\tA\tG\t.\t.\t.",
    "1\t260\t.\tC\t.\t.\t.\t.",
    "1\t270\t.\tT\t*\t.\t.\t."
  ), mixed)
  expect_setequal(read_variants_vcf(mixed)$pos, 250L)
  # more `##` metadata lines than a fixed 1000-line scan: the header is still found
  long <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste0("##INFO=<ID=T", seq_len(1500), ",Number=1,Type=Integer,Description=\"x\">"),
    hdr[-1L],
    "1\t250\t.\tA\tG\t.\t.\t."
  ), long)
  expect_equal(read_variants_vcf(long)$pos, 250L)
  # a header without ALT is refused with a readable reason, not a vroom internal
  noalt <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tQUAL\tFILTER\tINFO",
    "1\t250\t.\tA\t.\t.\t."
  ), noalt)
  expect_error(read_variants_vcf(noalt), "missing")
  # a file with no header at all errors instead of silently returning nothing
  bare <- withr::local_tempfile(fileext = ".vcf")
  writeLines("1\t250\t.\tA\tG", bare)
  expect_error(read_variants_vcf(bare), "#CHROM")
})

test_that("read_variants_vcf filters by region", {
  skip_if_not_installed("vroom")
  vcf <- withr::local_tempfile(fileext = ".vcf")
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
  expect_equal(normalise_allele("ACRT"), "ACNT") # ambiguity inside a longer allele
})

test_that("read_variants_vcf normalises IUPAC/lowercase instead of dropping", {
  skip_if_not_installed("vroom")
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "4\t500\t.\tA\tR\t.\t.\t.", # IUPAC ambiguity -> N
    "4\t600\t.\ta\tg\t.\t.\t." # lowercase is normalised
  ), vcf)
  v <- read_variants_vcf(vcf)
  expect_equal(v$alt[v$pos == 500L], "N")
  expect_equal(v$ref[v$pos == 600L], "A")
  expect_equal(v$alt[v$pos == 600L], "G")
})

test_that("the vroom and VariantAnnotation backends agree", {
  skip_if_not_full()
  skip_if_not_installed("vroom")
  # VariantAnnotation's namespace load emits a benign import-replacement warning
  if (!suppressWarnings(suppressMessages(requireNamespace("VariantAnnotation", quietly = TRUE)))) {
    testthat::skip("VariantAnnotation not installed")
  }
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "4\t100\t.\tA\tG\t.\t.\t.",
    "4\t200\t.\tAT\tA\t.\t.\t.", # deletion
    "4\t300\t.\tA\tG,T\t.\t.\t.", # multi-allelic
    "4\t400\t.\tA\t<DEL>\t.\t.\t." # symbolic alleles are dropped
  ), vcf)
  a <- read_variants_vcf(vcf, backend = "vroom")
  b <- suppressWarnings(read_variants_vcf(vcf, backend = "VariantAnnotation"))
  key <- function(x) sort(paste(x$chrom, x$pos, x$ref, x$alt))
  expect_equal(key(a), key(b))
  expect_setequal(a$pos, c(100L, 200L, 300L, 300L))
})

test_that("the VariantAnnotation backend filters by region too", {
  skip_if_not_full()
  skip_if_not_installed_noload("VariantAnnotation")
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t10\t.\tA\tG\t.\t.\t.",
    "1\t50\t.\tA\tC\t.\t.\t.",
    "1\t90\t.\tA\tT\t.\t.\t."
  ), vcf)
  v <- suppressWarnings(
    read_variants_vcf(vcf,
      chrom = "Chr1", start = 20L, end = 60L,
      backend = "VariantAnnotation"
    )
  )
  expect_setequal(v$pos, 50L)
  expect_true(all(v$chrom == "Chr1"))
})

test_that("the VariantAnnotation backend streams a region like filtering", {
  skip_if_not_full()
  skip_if_not_installed_noload("VariantAnnotation")
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t10\t.\tA\tG\t.\t.\t.",
    "1\t50\t.\tA\tC\t.\t.\t.",
    "1\t90\t.\tA\tT\t.\t.\t."
  ), vcf)
  reg <- suppressWarnings(read_variants_vcf(vcf,
    chrom = "Chr1", start = 20L, end = 60L,
    backend = "VariantAnnotation"
  ))
  all <- suppressWarnings(read_variants_vcf(vcf, backend = "VariantAnnotation"))
  ref <- dplyr::filter(all, .data$chrom == "Chr1", .data$pos >= 20L, .data$pos <= 60L)
  expect_equal(reg, ref)
})

test_that("read_genome_fasta normalises and filters sequence names", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(
    ">1 dna:chromosome chromosome:TAIR10:1:1:100:1 REF", "ACGTACGT",
    ">Mt dna:chromosome", "ACGT", ">scaffold_9 x", "TTTT"
  ), fa)
  g <- read_genome_fasta(fa, chrom_map = ATH_CHROM_MAP, use_index = FALSE)
  expect_setequal(names(g), c("Chr1", "ChrM")) # scaffold dropped
  expect_equal(as.character(g[["Chr1"]]), "ACGTACGT")
  expect_setequal(names(read_genome_fasta(fa, seqids = "1", use_index = FALSE)), "1")
})

test_that("a FASTA sequence with no name is refused", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">", "ACGTACGT", ">Chr2", "TTTTGGGG"), fa)
  expect_error(
    read_genome_fasta(fa, use_index = FALSE),
    "sequence with no name \\(record 1\\)"
  )
})

test_that("read_genome_fasta reads only requested ids via the FASTA index", {
  skip_if_not_installed("Rsamtools")
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">chr1", "ACGTACGTACGT", ">chr2", "TTTTGGGGCCCC", ">chr3", "GGGGCCCCAAAA"), fa)
  g <- read_genome_fasta(fa, seqids = c("chr2", "chr3"))
  expect_setequal(names(g), c("chr2", "chr3")) # chr1 not read
  expect_equal(as.character(g[["chr2"]]), "TTTTGGGGCCCC")
  expect_true(file.exists(paste0(fa, ".fai"))) # index created
})

test_that("read_genome_fasta map+index keeps only mapped, requested chroms", {
  skip_if_not_installed("Rsamtools")
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">1 desc", "ACGTACGT", ">2 desc", "TTTTGGGG"), fa)
  g <- read_genome_fasta(fa, chrom_map = ATH_CHROM_MAP, seqids = "Chr2")
  expect_equal(names(g), "Chr2")
  expect_equal(as.character(g[["Chr2"]]), "TTTTGGGG")
})

test_that("read_genome_fasta refuses a base outside the DNA alphabet", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">Chr1", "ACGTNXYZ"), fa)
  expect_error(
    read_genome_fasta(fa, use_index = FALSE),
    "outside the DNA/IUPAC alphabet"
  )
  # the indexed path reaches the same guard through scanFa
  expect_error(
    suppressWarnings(read_genome_fasta(fa, seqids = "Chr1")),
    "outside the DNA/IUPAC alphabet"
  )
})

test_that("a FASTA with one name twice is refused on both backends", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">Chr1", "ACGTACGTAC", ">Chr1", "TTTTGGGGCC"), fa)
  for (use_index in c(FALSE, TRUE)) {
    expect_error(
      read_genome_fasta(fa, seqids = "Chr1", use_index = use_index),
      "names sequence Chr1 more than once",
      class = "tidymarkers_duplicate_sequence_error"
    )
  }
  # two source names mapping to one target name is the same ambiguity
  two <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">1", "ACGT", ">2", "ACGT"), two)
  expect_error(
    read_genome_fasta(two, chrom_map = c("1" = "Chr1", "2" = "Chr1")),
    "more than once after mapping",
    class = "tidymarkers_duplicate_sequence_error"
  )
})

test_that("marker_example carries one of each target type", {
  ex <- marker_example()
  expect_setequal(unique(ex$type), c("SNP", "MNP", "insertion", "deletion"))
  expect_true(any(!is.na(ex$alt2))) # the multi-allelic example
  expect_equal(nrow(marker_example("deletion")), 1L)
  # each target names the gene it sits in (the panel fetches its model from this)
  expect_true("gene" %in% names(ex))
  expect_false(any(is.na(ex$gene)))
  # the labelled type matches variant_type() for every non-multi-allelic row
  one <- ex[is.na(ex$alt2), ]
  expect_equal(variant_type(one$ref, one$alt), one$type)
  # the two ALTs of the multi-allelic row are degenerate-coverable
  ma <- ex[!is.na(ex$alt2), ]
  expect_true(marker_site_class(ma$ref, c(ma$alt, ma$alt2))$disc)
})

# The example locus table on disk is what the exported examples are cut from, so
# it is checked as data: the transition column must be the base chemistry, not
# a hand-set flag.
test_that("example_targets: transition is purine-for-purine, and indels are not", {
  ex <- utils::read.csv(
    tm_extdata("example_targets.csv"),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
  expect_equal(nrow(ex), 9L)
  expect_setequal(unique(ex$type), c("SNP", "MNP", "insertion", "deletion"))
  purine <- c("A", "G")
  snp <- ex[ex$type == "SNP", ]
  expect_gt(nrow(snp), 0L)
  implied <- (snp$ref %in% purine) == (snp$alt %in% purine)
  expect_equal(implied, snp$transition)
  # a transition is a single-base concept: an indel or an MNP is never one
  expect_false(any(ex$transition[ex$type != "SNP"] %in% TRUE))
  # exactly one row is multi-allelic, and its second ALT sits on a substitution
  expect_equal(sum(!is.na(ex$alt2)), 1L)
  expect_true(ex$type[!is.na(ex$alt2)] %in% c("SNP", "MNP"))
})

test_that("read_genbank parses features, caches an index and reuses it", {
  gb <- withr::local_tempfile(fileext = ".gbk")
  unlink(paste0(gb, ".gbidx"))
  write_genbank <- function(extra = NULL) {
    writeLines(c(
      "LOCUS       TESTLOCUS         120 bp    DNA     linear   UNK 01-JAN-2026",
      "DEFINITION  Test record.",
      "ACCESSION   TESTACC",
      "VERSION     TESTACC.1",
      "FEATURES             Location/Qualifiers",
      "     source          1..120",
      '                     /organism="Arabidopsis thaliana"',
      '                     /mol_type="genomic DNA"',
      "     gene            complement(10..60)",
      '                     /gene="GENE1"',
      '                     /locus_tag="TAG1"',
      "     CDS             join(10..30,40..60)",
      '                     /gene="GENE1"',
      '                     /locus_tag="TAG1"',
      '                     /product="protein one"',
      "     misc_feature    100..110",
      '                     /note="a feature"',
      extra,
      "ORIGIN",
      "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt",
      "       61 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt",
      "//"
    ), gb)
  }

  write_genbank()
  x <- read_genbank(gb)
  expect_equal(nrow(x), 4L)
  expect_true(attr(x, "file_reads") >= 1L)
  expect_true(all(x$chrom == "TESTLOCUS"))

  cds <- x[x$type == "CDS", ]
  gene <- x[x$type == "gene", ]
  source <- x[x$type == "source", ]
  misc <- x[x$type == "misc_feature", ]
  expect_equal(cds$start, 10L)
  expect_equal(cds$end, 60L)
  expect_equal(cds$strand, "+")
  expect_equal(gene$start, 10L)
  expect_equal(gene$end, 60L)
  expect_equal(gene$strand, "-")
  expect_equal(source$start, 1L)
  expect_equal(source$end, 120L)
  expect_equal(source$strand, "+")
  expect_equal(misc$start, 100L)
  expect_equal(misc$end, 110L)

  expect_equal(cds$gene, "GENE1")
  expect_equal(cds$product, "protein one")
  expect_true(is.na(gene$product))
  expect_equal(gene$locus_tag, "TAG1")

  y <- read_genbank(gb, region = list(chrom = "TESTLOCUS", start = 1, end = 30))
  expect_equal(attr(y, "file_reads"), 0L)
  expect_equal(nrow(y), 3L)
  expect_setequal(y$type, c("source", "gene", "CDS"))

  write_genbank(extra = c("     misc_feature    80..90", '                     /note="added"'))
  z <- read_genbank(gb)
  expect_true(attr(z, "file_reads") >= 1L)
  expect_equal(nrow(z), 5L)
  expect_true(any(z$type == "misc_feature" & z$start == 80L))

  # a same-byte-size edit keeps `size` unchanged, so only a newer mtime can
  # invalidate the cache; bump mtime explicitly to make the test deterministic
  txt <- readLines(gb)
  txt[grepl("DEFINITION", txt)] <- sub("Test record", "Test rec0rd", txt[grepl("DEFINITION", txt)])
  writeLines(txt, gb)
  Sys.setFileTime(gb, Sys.time() + 5)
  w <- read_genbank(gb)
  expect_true(attr(w, "file_reads") >= 1L)
  expect_equal(nrow(w), 5L)
})

test_that("read_genbank keeps a nested complement inside join as mixed strand", {
  gb <- withr::local_tempfile(fileext = ".gbk")
  unlink(paste0(gb, ".gbidx"))
  writeLines(c(
    "LOCUS       TESTLOCUS         60 bp     DNA     linear   UNK 01-JAN-2026",
    "FEATURES             Location/Qualifiers",
    "     source          1..60",
    "     CDS             join(complement(10..30),40..60)",
    "ORIGIN",
    "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt",
    "//"
  ), gb)
  x <- read_genbank(gb)
  cds <- x[x$type == "CDS", ]
  expect_equal(cds$start, 10L)
  expect_equal(cds$end, 60L)
  expect_equal(cds$strand, "*")
})

test_that("read_genbank folds a wrapped qualifier value into one feature", {
  gb <- withr::local_tempfile(fileext = ".gbk")
  unlink(paste0(gb, ".gbidx"))
  writeLines(c(
    "LOCUS       TESTLOCUS         60 bp     DNA     linear   UNK 01-JAN-2026",
    "FEATURES             Location/Qualifiers",
    "     misc_feature    10..20",
    '                     /note="a long note that wraps',
    "                     across two lines",
    "ORIGIN",
    "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt",
    "//"
  ), gb)
  x <- read_genbank(gb)
  expect_equal(nrow(x), 1L)
  expect_equal(x$qualifiers[[1]]$note, "a long note that wraps across two lines")
})

test_that("read_genbank folds a wrapped join() location into full start/end", {
  gb <- withr::local_tempfile(fileext = ".gbk")
  unlink(paste0(gb, ".gbidx"))
  writeLines(c(
    "LOCUS       TESTLOCUS        300 bp    DNA     linear   UNK 01-JAN-2026",
    "FEATURES             Location/Qualifiers",
    "     CDS             join(1..100,",
    "                     150..300)",
    "ORIGIN",
    "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt",
    "//"
  ), gb)
  x <- read_genbank(gb)
  cds <- x[x$type == "CDS", ]
  expect_equal(cds$start, 1L)
  expect_equal(cds$end, 300L)
})

test_that("read_variant_genotypes reads GT/AD into long dosage rows", {
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tL1\tL2\tL3",
    "1\t250\t.\tA\tG\t.\t.\t.\tGT:AD\t0/0:10,0\t0/1:6,4\t1/1:0,9",
    "1\t250\t.\tA\tG,T\t.\t.\t.\tGT:AD\t1/2:5,5\t0/1:5,5\t0/0:9,0",
    "1\t300\t.\tC\tT\t.\t.\t.\tAD\t5,5\t0,10\t.",
    "1\t400\t.\tG\tA\t.\t.\t.\tGT\t./.\t0|1\t1|1"
  ), vcf)

  expect_message(
    g <- read_variant_genotypes(vcf),
    "multi-allelic"
  )
  expect_s3_class(g, "genotypes")
  expect_identical(
    names(g),
    c("line", "chrom", "pos", "ref", "alt", "genotype", "depth")
  )
  expect_equal(nrow(g), 9L)
  expect_identical(unique(g$chrom), "Chr1")

  at <- function(pos, line) g$genotype[g$pos == pos & g$line == line]
  expect_equal(at(250L, "L1")[[1L]], 0)
  expect_equal(at(250L, "L2")[[1L]], 1)
  expect_equal(at(250L, "L3")[[1L]], 2)
  # AD decides the dosage when GT is absent; `.` is no call
  expect_equal(at(300L, "L1"), 1)
  expect_equal(at(300L, "L2"), 2)
  expect_true(is.na(at(300L, "L3")))
  # a phased call uses the same dosage, a missing call stays NA
  expect_equal(at(400L, "L2"), 1)
  expect_equal(at(400L, "L3"), 2)
  expect_true(is.na(at(400L, "L1")))
  expect_equal(g$depth[g$pos == 250L & g$line == "L3"], 9)

  one <- read_variant_genotypes(vcf, samples = "L2")
  expect_equal(nrow(one), 3L)
  expect_identical(unique(one$line), "L2")
  expect_error(read_variant_genotypes(vcf, samples = "nope"), "no sample")
})

test_that("genotypes feed marker_informativeness", {
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tL1\tL2\tL3\tL4",
    "1\t250\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t0/1",
    "1\t300\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t0/0\t1/1"
  ), vcf)
  g <- read_variant_genotypes(vcf)
  targets <- tibble::tibble(
    chrom = "Chr1", pos = c(250L, 300L), ref = c("A", "C"), alt = c("G", "T")
  )
  info <- marker_informativeness(targets, genotypes = g, min_lines = 4L)
  expect_equal(nrow(info), 2L)
  expect_true(all(info$n_lines >= 4L))
})

test_that("as_variant_row gives exactly one typed row", {
  v <- marker_example("SNP")
  row <- as_variant_row(v, chrom = "Chr1")
  expect_s3_class(row, "variant")
  expect_equal(nrow(row), 1L)
  expect_identical(names(row), c("chrom", "pos", "ref", "alt"))

  two <- tibble::tibble(
    id = c("a", "b"), chrom = c("1", "1"), pos = c(10L, 20L),
    ref = c("A", "C"), alt = c("G", "T")
  )
  expect_identical(as_variant_row(two, id = "b")$pos, 20L)
  expect_identical(as_variant_row(two, pos = 10L)$ref, "A")
  expect_error(as_variant_row(two), "exactly one variant")
  expect_error(as_variant_row(two, pos = 99L), "exactly one variant")

  # the selected row is admitted before ath_chrom() can turn "" into "ChrNA"
  bad <- tibble::tibble(chrom = "", pos = 10L, ref = "A", alt = "G")
  expect_error(as_variant_row(bad), "chromosome name")
  bad$chrom <- "Chr1"
  bad$pos <- 0L
  expect_error(as_variant_row(bad), "positive integer `pos`")
  bad$pos <- -1L
  expect_error(as_variant_row(bad), "positive integer `pos`")

  # a symbolic or absent allele is refused, not renamed to a fake N (ADR 0018)
  bad$pos <- 10L
  for (alt in c("<DEL>", "*", ".", NA_character_)) {
    bad$alt <- alt
    expect_error(as_variant_row(bad), "concrete REF and ALT")
  }
})

test_that("a VCF POS outside integer range is refused by every reader", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "Chr1\t0\t.\tA\tC\t.\t.\t.\tGT\t0/1"
  ), path)

  expect_error(read_variants_vcf(path), "positive integer range")

  # the streaming reader once dropped the record and said why; a record one
  # reader refuses must not become an empty table in another (ADR 0016, 0018)
  expect_error(
    read_variant_genotypes(path),
    "record 1 has a POS outside the positive integer range"
  )
  expect_error(read_variant_depth(path), "positive integer range")
})

test_that("a VCF record with no chromosome is refused", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "\t300\t.\tA\tG\t.\tPASS\t."
  ), path)

  # ath_chrom("") would otherwise produce the plausible "ChrNA" contig
  expect_error(read_variants_vcf(path), "no chromosome name")

  # the genotype reader is the second reader of the same file and must not
  # invent a contig the variant reader refuses (ZV15 class, chromosome field)
  gpath <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"
  ), gpath)
  expect_error(read_variant_genotypes(gpath), "record 1 has no chromosome")
})

test_that("a VCF header without ALT names the missing column", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF"
  ), path)

  expect_error(read_variants_vcf(path), "missing column.*ALT")
})

test_that("read_variant_genotypes drops a symbolic or absent ALT like the variant reader", {
  skip_if_not_installed("vroom")
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tL1",
    "1\t100\t.\tA\t<DEL>\t.\t.\t.\tGT:AD\t0/1:5,5",
    "1\t200\t.\tA\t.\t.\t.\t.\tGT:AD\t0/1:5,5",
    "1\t300\t.\tA\tG\t.\t.\t.\tGT:AD\t0/1:5,5"
  ), vcf)

  expect_message(
    g <- read_variant_genotypes(vcf),
    "symbolic or absent ALT"
  )
  # no fake `N` allele survives, and both readers keep the same records
  expect_identical(g$pos, 300L)
  expect_identical(g$alt, "G")
  v <- suppressMessages(read_variants_vcf(vcf))
  expect_identical(v$pos, g$pos)
})

test_that("an absent REF or ALT is dropped by every reader, not a raw failure", {
  skip_if_not_installed("vroom")
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "1\t100\t.\tA\t\t.\t.\t.\tGT:AD\t0/1:5,5",
    "1\t200\t.\t\tG\t.\t.\t.\tGT:AD\t0/1:5,5",
    "1\t300\t.\tA\tG\t.\t.\t.\tGT:AD\t0/1:5,5"
  ), vcf)

  # an empty field used to reach `if (NA > 0L)` in the variant reader and
  # raise a base-R error; the other two readers disagreed (one kept the row)
  expect_message(v <- read_variants_vcf(vcf), "2 VCF records")
  expect_message(g <- read_variant_genotypes(vcf), "2 VCF records")
  d <- suppressMessages(read_variant_depth(vcf))

  expect_identical(v$pos, 300L)
  expect_identical(g$pos, 300L)
  expect_identical(d$pos, 300L)
  expect_identical(d$ref, "A")
  expect_false(anyNA(v$ref) || anyNA(v$alt) || any(!nzchar(v$alt)))
})

test_that("a truncated VCF record is counted as short, not as missing FORMAT", {
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tL1",
    "1\t100\t.\tA\tG\t.\t.\t.\tGT:AD\t0/1:5,5",
    "1\t200\t.\tA\tG\t.\t.\t."
  ), vcf)
  expect_message(g <- read_variant_genotypes(vcf), "fewer columns")
  expect_identical(g$pos, 100L)
})
