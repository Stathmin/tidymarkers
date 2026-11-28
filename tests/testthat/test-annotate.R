# Structured annotation: a TxDb gives the consequence instead of the input string.

# Close the SQLite connections the annotation objects (TxDb/EnsDb) open, so
# RSQLite does not warn ("call dbDisconnect()...") when R shuts down afterwards.
withr::defer(
  {
    try(
      {
        for (dbc in DBI::dbListConnections(RSQLite::SQLite())) DBI::dbDisconnect(dbc)
        suppressWarnings(gc())
      },
      silent = TRUE
    )
  },
  teardown_env()
)

test_that("annotate_variants reads the gene model (locateVariants)", {
  skip_if_not_full()
  skip_if_not_installed_noload("TxDb.Athaliana.BioMart.plantsmart28")
  skip_if_not_installed_noload("VariantAnnotation")
  v <- tibble::tibble(
    id = "fri", chrom = "Chr4", pos = 269745L,
    ref_genome = list(c("G", "A")), consequence = NA_character_
  )
  a <- annotate_variants(v)
  expect_equal(a$gene_id, "AT4G00650") # FRI
  expect_equal(a$consequence, "coding_sequence_variant")
})

test_that("a BSgenome reference gives codon-level consequences", {
  skip_if_not_full()
  skip_if_not_installed_noload("TxDb.Athaliana.BioMart.plantsmart28")
  skip_if_not_installed_noload("BSgenome.Athaliana.TAIR.TAIR9")
  skip_if_not_installed_noload("VariantAnnotation")
  v <- tibble::tibble(
    id = "fri", chrom = "Chr4", pos = 269745L,
    ref_genome = list(c("G", "A")), consequence = NA_character_
  )
  a <- annotate_variants(v, seqsource = marker_bsgenome())
  expect_equal(a$consequence, "stop_gained") # TGG becomes TGA
  expect_equal(a$codon_change, "TGG>TGA")
  # the codon call outranks the region-level term
  expect_equal(a$gene_id, "AT4G00650")
})

test_that("marker_txdb / marker_bsgenome load the default annotation packages", {
  skip_if_not_full()
  skip_if_not_installed_noload("TxDb.Athaliana.BioMart.plantsmart28")
  skip_if_not_installed_noload("BSgenome.Athaliana.TAIR.TAIR9")
  expect_s4_class(marker_txdb(), "TxDb")
  expect_s4_class(marker_bsgenome(), "BSgenome")
})

test_that("an empty or already-annotated table is returned unchanged", {
  # the other tests here guard the same way: VariantAnnotation is a Suggest, and
  # an empty table still goes through the loader before it can be short-circuited
  skip_if_not_installed_noload("VariantAnnotation")
  empty <- tibble::tibble(
    id = character(), chrom = character(), pos = integer(),
    ref_genome = list(), consequence = character()
  )
  expect_equal(nrow(annotate_variants(empty)), 0L)
})

test_that("codon_strings() reads codons from every column shape", {
  # `predictCoding()` hands the codons over as an XStringSet on release and as a
  # list-like column on newer Biostrings; both must become one string per row.
  expect_equal(tidymarkers:::codon_strings(c("TGG", "TGA")), c("TGG", "TGA"))
  set <- Biostrings::DNAStringSet(c("TGG", "TGA"))
  expect_equal(tidymarkers:::codon_strings(set), c("TGG", "TGA"))
  plain <- list(Biostrings::DNAString("TGG"), Biostrings::DNAString("TGA"))
  expect_equal(tidymarkers:::codon_strings(plain), c("TGG", "TGA"))
})

test_that("marker_ensdb builds an EnsDb from a GTF", {
  skip_if_not_full()
  skip_if_not_installed("ensembldb")
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("processx")
  gtf <- withr::local_tempfile(fileext = ".gtf")
  row <- function(feature, start, end, frame, extra) {
    paste("1", "src", feature, start, end, ".", "+", frame, extra, sep = "\t")
  }
  writeLines(c(
    row("gene", 100L, 200L, ".", "gene_id \"G1\"; gene_name \"FRI\"; gene_biotype \"protein_coding\";"),
    # a real Ensembl GTF always carries the transcript feature; without it
    # ensembldb takes its derived `range()` branch and leaves the CDS bounds as
    # logical `NA` columns, which fails its own integer-column check
    row("transcript", 100L, 200L, ".", paste0(
      "gene_id \"G1\"; transcript_id \"T1\"; ",
      "transcript_biotype \"protein_coding\";"
    )),
    row("exon", 100L, 200L, ".", "gene_id \"G1\"; transcript_id \"T1\"; exon_number \"1\";"),
    row("CDS", 110L, 190L, "0", "gene_id \"G1\"; transcript_id \"T1\"; exon_number \"1\";")
  ), gtf)
  # Build in a subprocess: ensDbFromGtf() leaves an RSQLite connection that
  # otherwise leaks a benign "call dbDisconnect()" note at this session's exit.
  cache <- withr::local_tempfile()
  # the code string is parsed by the child R, so a Windows path must not carry
  # single backslashes into it
  gtf_arg <- normalizePath(gtf, winslash = "/", mustWork = FALSE)
  cache_arg <- normalizePath(cache, winslash = "/", mustWork = FALSE)
  code <- sprintf(
    'edb <- tidymarkers::marker_ensdb(gtf = "%s", cache = "%s"); cat("GENES", length(GenomicFeatures::genes(edb)))',
    gtf_arg, cache_arg
  )
  res <- processx::run(file.path(R.home("bin"), "Rscript"), c("-e", code),
    error_on_status = FALSE
  )
  expect_equal(res$status, 0L, info = res$stderr)
  expect_match(res$stdout, "GENES 1")
})

# the cache names are refused before the engine is even looked for, so this
# runs without `ensembldb` installed
test_that("marker_ensdb refuses a link at either cache name", {
  cache <- withr::local_tempdir()
  outside <- fs::path(withr::local_tempdir(), "foreign")
  writeLines("foreign", outside)
  skip_if(
    !file.symlink(outside, fs::path(cache, "Arabidopsis_thaliana.EnsDb.sqlite")),
    "this filesystem cannot make a symlink"
  )
  gtf <- withr::local_tempfile(fileext = ".gtf")
  writeLines("foreign", gtf)
  expect_error(marker_ensdb(gtf = gtf, cache = cache), "which is a link")

  fs::file_delete(fs::path(cache, "Arabidopsis_thaliana.EnsDb.sqlite"))
  skip_if(
    !file.symlink(outside, fs::path(cache, "Arabidopsis_thaliana.TAIR10.63.gtf.gz")),
    "this filesystem cannot make a symlink"
  )
  expect_error(marker_ensdb(cache = cache), "which is a link")
})
