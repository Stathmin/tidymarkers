# Lax KASP: one ALT primer per distinct ALT haplotype, on a shared backbone.
# Exercised on a real local chunk of the SAT2 transcripts (no synthetic fixture).

lax_hp <- local({
  f <- system.file("extdata", "sat2_transcripts.fa", package = "tidykasp")
  if (!nzchar(f)) {
    return(NULL)
  }
  haplotypes_from_alignment(read_transcripts(f), chrom = "SAT2")
})

# a 301 bp chunk around SAT2:431 (window-local position 151)
lax_chunk <- function() {
  if (is.null(lax_hp)) {
    return(NULL)
  }
  kasp_chunk(lax_hp, from = 281L, to = 581L)
}

test_that("design_lax_marker groups matching ALT haplotypes into one primer each", {
  skip_if(is.null(lax_hp), "bundled SAT2 FASTA not found")
  skip_if_no_oligoarrayaux()
  ch <- lax_chunk()
  snp <- tibble::tibble(chrom = "SAT2", pos = 151L, ref = "A", alt = "G")
  x <- suppressWarnings(design_lax_marker(
    chunk_genome(ch), snp, chunk_haplotypes(ch),
    flank = 150L, cores = 1L
  ))
  expect_s3_class(x, "kasp_lax")
  # every ALT haplotype is covered, but haplotypes with the same linked variants
  # collapse to a single ALT primer group
  expect_gte(length(x$alt), 5L)
  n_groups <- length(unique(x$haplotype_groups$group))
  expect_lt(n_groups, length(x$alt))
  expect_equal(nrow(x$haplotype_groups), length(x$alt))
  # the tube: one shared REF primer and one shared COM, plus the ALT primers
  expect_equal(sum(names(x$pool) == "AS_ref"), 1L)
  expect_equal(sum(names(x$pool) == "COM"), 1L)
  expect_gte(sum(!is.na(x$alt_primers)), length(x$alt) - 1L)
  # the *joint* pool is screened, and it is clean here
  expect_true(is.data.frame(x$compatibility))
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(x$compatibility)))
  expect_true(all(x$compatibility$pass))
  expect_equal(x$verdict, "pass")
  expect_output(print(x), "kasp_lax")
})

test_that("design_lax_marker degenerates to one plain ALT primer without haplotypes", {
  skip_if(is.null(lax_hp), "bundled SAT2 FASTA not found")
  skip_if_no_oligoarrayaux()
  ch <- lax_chunk()
  snp <- tibble::tibble(chrom = "SAT2", pos = 151L, ref = "A", alt = "G")
  x <- suppressWarnings(design_lax_marker(chunk_genome(ch), snp, flank = 150L, cores = 1L))
  expect_s3_class(x, "kasp_lax")
  expect_equal(x$alt, "ALT")
  expect_length(x$assays, 1L)
  expect_output(print(x), "plain ALT primer")
})

test_that("the lax scheme is registered and has its own decision stage", {
  expect_true("lax" %in% kasp_schemes()$scheme)
  s <- kasp_decision_schema("lax")
  expect_true("lax" %in% s$nodes$stage)
  expect_true(all(c("as", "pair") %in% s$nodes$stage))
})
