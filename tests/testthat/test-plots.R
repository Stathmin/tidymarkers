# The two-level plot vocabulary: feature level (one oligo) vs alignment level.

test_that("check_level splits criteria into feature and alignment", {
  expect_equal(check_level("AS Tm (C)"), "feature")
  expect_equal(check_level("worst hairpin dG37"), "feature")
  expect_equal(check_level("amplicon length (bp)"), "alignment")
  expect_equal(check_level("in-silico PCR (REF 1 / ALT 0)"), "alignment")
})

test_that("clamp_window follows the primer's 3' end on either strand", {
  expect_equal(clamp_window(10L, 19L, forward = TRUE, n = 5L), c(15L, 19L))
  expect_equal(clamp_window(10L, 19L, forward = FALSE, n = 5L), c(10L, 14L))
})

test_that("check values are monospaced, <= 4 chars, formatted per criterion", {
  expect_equal(format_check_value("AS GC", 0.444), "44%")
  expect_equal(format_check_value("COM GC", 0.571), "57%")
  expect_equal(format_check_value("AS Tm", 63.526), "63.5")
  expect_equal(format_check_value("dimer", -3.26), "-3.3")
  expect_equal(format_check_value("homodimer", -1.63), "-1.6")
  expect_equal(format_check_value("hairpin", -2.294), "-2.3")
  expect_equal(format_check_value("hairpin", 0.544), "0.5")
  expect_equal(format_check_value("COM\u2265AS Tm", 0.409), "0.4")
  expect_equal(format_check_value("COM\u2265AS Tm", -3.47), "-3.5")
  expect_equal(format_check_value("discrim.", 0.932), "93%")
  expect_equal(format_check_value("run", 3), "3")
  expect_equal(format_check_value("amplicon", 67), "67")
  expect_equal(format_check_value("PCR", 1), "1")
  expect_equal(format_check_value("PCR", NA_real_), "-")
  # extreme counts compact to dK / dM
  expect_equal(format_check_value("amplicon", 1400), "1K")
  expect_equal(format_check_value("amplicon", 2.6e6), "3M")

  names <- c("AS GC", "AS Tm", "dimer", "homodimer", "hairpin", "run", "amplicon", "COM\u2265AS Tm", "discrim.", "PCR")
  vals <- c(0.44, 63.5, -3.3, -1.6, -2.29, 3, 67, 0.4, 0.93, 1)
  expect_true(all(nchar(format_check_value(names, vals)) <= 4L))
})

# A tiny real assay on a random genome (skips without OligoArrayAux).
make_plot_assay <- function(seed = 4L, n = 500L, pos = 250L) {
  for (s in seed:(seed + 19L)) {
    withr::with_seed(s, {
      w <- paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
      genome <- Biostrings::DNAStringSet(c(Chr1 = w))
    })
    ref <- as.character(Biostrings::subseq(genome[["Chr1"]], pos, pos))
    snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
    vs <- tibble::tibble(chrom = "Chr1", pos = pos, consequence = "missense_variant", pop_freq = 0.3)
    a <- try(suppressWarnings(assemble_kasp_assay(genome, snp, variants = vs, flank = 150L)), silent = TRUE)
    if (!inherits(a, "try-error")) {
      return(list(assay = a, snp = snp, variants = vs, genome = genome))
    }
  }
  testthat::skip("no random genome assembled into an assay")
}

test_that("the assay plots render and carry the two-level vocabulary", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  m <- make_plot_assay()

  feats <- primer_features(m$assay)
  expect_setequal(feats$role, c("AS-ref", "AS-alt", "COM"))
  expect_true(all(feats$status %in% c("ok", "broken")))
  expect_true(all(feats$clamp %in% c("ok", "broken")))

  # styled letters: clamp bases are bold, motif bases are red, rest neutral
  let <- primer_letters(m$assay)
  expect_setequal(unique(let$face), c("plain", "bold"))
  expect_true("black" %in% let$colour)
  # the 3' GC-clamp is the LAST letters of the string, on either strand
  for (r in unique(let$role)) {
    sub <- let[let$role == r, ]
    expect_true(all(utils::head(sub$face, -5L) == "plain"), info = r)
    expect_true(all(utils::tail(sub$face, 5L) == "bold"), info = r)
  }

  expect_s3_class(plot_checks(m$assay), "ggplot")
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp), "patchwork")
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp, genome = m$genome), "patchwork")
  expect_s3_class(plot_alignment(m$assay), "ggplot")

  # an optional pop_freq weights the footprint load and sizes the markers
  expect_true("w_var" %in% names(m$assay$assay))
})

# the genome-wide ideogram (exact primer hits + in-silico PCR products)
test_that("assay_ideogram marks primer hits and in-silico products", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  m <- make_plot_assay()
  expect_null(assay_ideogram(m$assay, m$snp, NULL))
  expect_s3_class(assay_ideogram(m$assay, m$snp, m$genome), "ggplot")

  hits <- primer_hits(assay_primers(m$assay), m$genome)
  expect_setequal(names(hits), c("primer", "chrom", "pos"))
  expect_true(all(hits$primer %in% c("AS_ref", "AS_alt", "COM")))
  # on the reference the REF/common primers land at the target locus; the ALT
  # primer carries the alternate base, so it has no exact hit there (by design)
  at_target <- hits$chrom == m$snp$chrom & abs(hits$pos - m$snp$pos) < 150L
  expect_true(all(c("AS_ref", "COM") %in% hits$primer[at_target]))

  # the shared pairing helper is the one assay_amplicons() summarises
  pairs <- as_pcr_pairs(m$assay)
  expect_setequal(names(pairs), c("AS_ref", "AS_alt"))
  expect_true(all(vapply(pairs, function(p) all(c("F", "R") %in% names(p)), logical(1))))
})

test_that("the feature row draws annotation, or notes its absence", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  m <- make_plot_assay()
  # no annotation: the row says "no feature annotation", still a patchwork
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp), "patchwork")
  feat <- tibble::tibble(start = m$snp$pos - 40L, end = m$snp$pos + 40L, strand = "+", name = "GENE")
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp, features = feat), "patchwork")
  # a gene record (start/end/strand/symbol) works the same way
  gene <- list(symbol = "FRI", start = m$snp$pos - 30L, end = m$snp$pos + 30L, strand = "-")
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp, features = gene), "patchwork")
})
