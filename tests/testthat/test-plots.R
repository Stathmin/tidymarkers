test_that("clamp_window follows the primer's 3' end on either strand", {
  expect_equal(clamp_window(10L, 19L, forward = TRUE, n = 5L), c(15L, 19L))
  expect_equal(clamp_window(10L, 19L, forward = FALSE, n = 5L), c(10L, 14L))
})

test_that("one letter size rule covers primers and the amplicon they sit on", {
  # a short KASP product and a long Sanger read share the same scale within a
  # plot, but a long product shrinks the whole stack
  expect_equal(letter_size_for_span(c(1L, 100L)), LETTER_SIZE_DEFAULT)
  expect_lt(letter_size_for_span(c(1L, 600L)), letter_size_for_span(c(1L, 100L)))
  expect_equal(letter_size_for_span(c(1L, 600L)), 0.9)
})

test_that("every registered check has a brief tumbler name", {
  crit <- marker_criteria()
  chk <- crit$label[crit$role == "check"]
  short <- check_short(chk)
  expect_false(anyNA(short))
  expect_true(all(nzchar(short)))
  expect_true(all(nchar(short) <= 16L))
  # the label -> short table is the single, backtraceable source
  expect_equal(check_short(names(CHECK_SHORT_MAP)), unname(CHECK_SHORT_MAP))
  expect_setequal(chk, names(CHECK_SHORT_MAP))
})

test_that("map_spec is the only amplicon/primer row constructor", {
  m <- list(
    amplicon = tibble::tibble(start = 10L, end = 100L),
    primers = list(forward = "ACGT", reverse = "TGCA")
  )
  spec <- map_spec(
    m,
    features = map_pair_rows(m, 10L, 100L),
    sequences = list(forward = "ACGT", reverse = "TGCA")
  )
  expect_equal(spec$feature, c("amplicon", "forward", "reverse"))
  expect_equal(spec$seq, c(NA_character_, "ACGT", "TGCA"))
  expect_equal(spec$start, c(10L, 10L, 97L))
})

test_that("only primers inherit the 3' clamp", {
  expect_equal(
    feature_clamp_n(c("primer", "seq", "probe", "guide", "site", "blocker")),
    c(5L, 5L, 0L, 0L, 0L, 0L)
  )
})

test_that("primer nose carries Tm text and a status dot", {
  spec <- tibble::tibble(
    feature = c("forward", "reverse"),
    start = c(10L, 90L), end = c(29L, 109L),
    strand = c("+", "-"), kind = "primer",
    seq = c("ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"),
    y = c(2L, 1L)
  )
  nose <- primer_nose_spec(spec, "caps")
  expect_equal(nrow(nose), 2L)
  expect_true(all(grepl("\u00b0C", nose$label)))
  expect_true(all(c("label_x", "dot_x", "hjust", "status") %in% names(nose)))
})

test_that("a generic map replaces the SNP strip with a haplotype track", {
  x <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "G"),
      amplicon = tibble::tibble(start = 250L, end = 350L),
      primers = list(
        forward = "ACGTACGTACGTACGTACGT", reverse = "TGCATGCATGCATGCATGCA"
      ),
      site = NULL,
      haplotypes = tibble::tibble(
        haplotype = c("H1", "H2", "H1"), freq = c(0.6, 0.4, 0.6),
        chrom = "Chr1", pos = c(300L, 300L, 310L),
        ref = c("A", "A", "C"), alt = c("G", "G", "T")
      )
    ),
    class = "marker_caps"
  )
  expect_s3_class(plot_marker_map(x, "caps"), "patchwork")
})

test_that("gel_spec owns band colour and Poisson width", {
  g <- gel_spec(c("REF", "ALT"), c(100, 120), c(TRUE, FALSE))
  expect_equal(g$col, c(GEL_DIAG, GEL_SHARED))
  expect_true(all(g$width > 0.2 & g$width <= 0.4))
  expect_equal(g$x, c(1, 2))
})

test_that("clamp_n = 0 keeps every base plain", {
  d <- oligo_letters("ACGTACGT", 1L, 8L, TRUE, "probe", clamp_n = 0L)
  expect_true(all(d$face == "plain"))
})

test_that("oligo_letters displays a non-DNA mark without DNA structure lookup", {
  d <- oligo_letters("ACGTXACGT", 1L, 9L, TRUE, "probe")
  expect_equal(nrow(d), 9L)
  expect_equal(d$letter[5L], "X")
})

test_that("oligo_letters accepts reusable span and base marks", {
  d <- oligo_letters("ACGTACGT", 10L, 17L, TRUE, "p",
    spans = tibble::tibble(from = 2L, to = 3L, colour = "#ff0000"),
    bases = tibble::tibble(pos = 2L, face = "bold")
  )
  expect_equal(d$colour[d$x == 11L], "#ff0000")
  expect_equal(d$colour[d$x == 12L], "#ff0000")
  expect_equal(d$face[d$x == 11L], "bold")
  # minus strand reads right-to-left; local sequence indices do not change
  m <- oligo_letters("ACGT", 17L, 20L, FALSE, "p",
    bases = tibble::tibble(pos = 1L, colour = "#00ff00")
  )
  expect_equal(m$x[1L], 20L)
  expect_equal(m$colour[1L], "#00ff00")
})

test_that("spec_row carries the shared status/label highlight contract", {
  r <- spec_row("probe", 10L, 20L, kind = "probe", status = "broken", label = "P")
  expect_equal(r$status, "broken")
  expect_equal(r$label, "P")
  e <- spec_bind(list())
  expect_named(e, c("feature", "start", "end", "strand", "kind", "status", "label"))
})

test_that("check values are monospaced, <= 4 chars, formatted per criterion", {
  expect_equal(format_check_value("AS GC", 0.444), "44%")
  expect_equal(format_check_value("COM GC", 0.571), "57%")
  expect_equal(format_check_value("AS Tm", 63.526), "63.5")
  expect_equal(format_check_value("Tm", 60.24), "60.2")
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

test_that("the tumbler draws NA for neutral rows and codes for categorical ones", {
  # prose in `detail` must not reach the figure: it overruns the narrow column
  long <- "not evaluated: no population variant table was supplied"
  ev <- tibble::tibble(
    level = c("site", "feature"), criterion = c("ALT freq", "amplicon length (bp)"),
    status = c("neutral", "pass"), detail = c(long, NA_character_),
    observed = c(NA_real_, 80), value = c(long, "80")
  )
  p <- tumbler_plot(ev)
  expect_equal(p$data$value, c("NA", "80"))
  expect_equal(ev$detail[[1L]], long)

  # categorical rows: the site class is one word, the Kwok class its first word
  fake <- list(
    population = NULL,
    evaluation = tibble::tibble(
      level = c("feature", "feature", "feature", "site"),
      criterion = c(rep("AS 3' mismatch base pair", 3L), "target site class"),
      observed = NA_real_,
      status = c("pass", "pass", "unknown", "pass"),
      detail = c(
        "severe terminal mismatch (Kwok 1990); discrimination is kinetics, not a binary",
        "moderate terminal mismatch (Kwok 1990); discrimination is kinetics, not a binary",
        "substitution pair not in Kwok 1990's table (no published class)",
        "biallelic"
      )
    )
  )
  pc <- plot_checks(fake)
  expect_equal(pc$data$value, c("sev", "mod", "NA", "bial"))
  expect_false(any(grepl("not a binary", pc$data$value, fixed = TRUE)))
  expect_true(all(nchar(pc$data$value) <= 4L))
})

# A tiny real assay on a random genome (skips without OligoArrayAux).  The
# assembly is deterministic and costs seconds, so build each fixture once and let
# the tests share it.
make_plot_assay <- local({
  cache <- new.env(parent = emptyenv())
  function(seed = 4L, n = 500L, pos = 250L) {
    key <- paste(seed, n, pos, sep = ":")
    if (!is.null(cache[[key]])) {
      return(cache[[key]])
    }
    for (s in seed:(seed + 19L)) {
      withr::with_seed(s, {
        w <- paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
        genome <- Biostrings::DNAStringSet(c(Chr1 = w))
      })
      ref <- as.character(Biostrings::subseq(genome[["Chr1"]], pos, pos))
      snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
      vs <- tibble::tibble(chrom = "Chr1", pos = pos, consequence = "missense_variant", pop_freq = 0.3)
      a <- try(suppressWarnings(assemble_kasp_assay(genome, snp, variants = vs, flank = 150L)), silent = TRUE)
      if (!inherits(a, "try-error") && !inherits(a, "marker_doomed")) {
        cache[[key]] <- list(assay = a, snp = snp, variants = vs, genome = genome)
        return(cache[[key]])
      }
    }
    testthat::skip("no random genome assembled into an assay")
  }
})

test_that("the assay plots render and carry the two-level vocabulary", {
  skip_if_no_oligoarrayaux()
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

  # an optional pop_freq weights the footprint load and sizes the markers
  expect_true("w_var" %in% names(m$assay$assay))
  # the ideogram inputs and the selection outcome are carried on the object, so
  # plotting only reads it (no genome rescan, no re-running the selection)
  expect_true(!is.null(m$assay$ideogram))
  expect_true(all(c("relaxed", "violated") %in% names(m$assay$selection)))
})

# the genome-wide ideogram (exact primer hits + in-silico PCR products)
test_that("assay_ideogram marks primer hits and in-silico products", {
  skip_if_no_oligoarrayaux()
  m <- make_plot_assay()
  # the inputs are computed once and stored on the assay, so plotting needs no genome
  expect_s3_class(assay_ideogram(m$assay, m$snp), "ggplot")
  expect_s3_class(assay_ideogram(m$assay, m$snp, m$genome), "ggplot")
  # without stored data and without a genome there is nothing to draw
  bare <- m$assay
  bare$ideogram <- NULL
  expect_null(assay_ideogram(bare, m$snp, NULL))

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
  skip_if_no_oligoarrayaux()
  m <- make_plot_assay()
  # no annotation: the row says "no feature annotation", still a patchwork
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp), "patchwork")
  feat <- tibble::tibble(start = m$snp$pos - 40L, end = m$snp$pos + 40L, strand = "+", name = "GENE")
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp, features = feat), "patchwork")
  # a gene record (start/end/strand/symbol) works the same way
  gene <- list(symbol = "FRI", start = m$snp$pos - 30L, end = m$snp$pos + 30L, strand = "-")
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp, features = gene), "patchwork")
})

test_that("the variant landscape draws a bare target without an `id` column", {
  skip_if_no_oligoarrayaux()
  m <- make_plot_assay()
  # a bare (chrom, pos, ref, alt) target carries no `id` - the label must come
  # from the variant itself, not a column the caller never supplies
  gene <- list(symbol = "GENE", chrom = "Chr1", start = 1L, end = 500L, strand = "+")
  p <- plot_variant_landscape(gene, m$variants, m$snp, m$assay)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetZoom") # the zoom is facet_zoom(), not a manual inset
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("a doomed target is a classed object and still plots its context", {
  skip_if_no_oligoarrayaux()
  m <- make_plot_assay()
  # force the hard footprint gate to fail for every placement
  variants <- tibble::tibble(
    chrom = "Chr1", pos = c(m$snp$pos, m$snp$pos - 6L, m$snp$pos + 6L),
    ref = c(m$snp$ref, "A", "A"), alt = c(m$snp$alt, "C", "T")
  )
  d <- suppressWarnings(assemble_kasp_assay(
    m$genome, m$snp,
    variants = variants, flank = 150L, max_footprint_variants = 0L
  ))
  expect_s3_class(d, "marker_doomed")
  expect_true(nzchar(d$reason))
  expect_output(print(d), "marker_doomed")
  # plot_primers() dispatches to the context plot (a plain ggplot, not patchwork)
  p <- plot_primers(d, variants, m$snp)
  expect_s3_class(p, "ggplot")
})

test_that("plot_oligo_structure renders the primers' hairpin/SSR structure", {
  skip_if_no_oligoarrayaux()
  m <- make_plot_assay()
  p <- plot_oligo_structure(m$assay)
  expect_s3_class(p, "ggplot")
  # one facet per oligo (its strip labels the hairpin/SSR counts)
  expect_equal(nlevels(p$data$panel), 3L)
  # it must actually build (a 0-column layer would fail only at build time)
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("plot_self_complementarity marks every Watson-Crick pairing", {
  s <- "GAAGGTGACCAAGTTCATGCT"
  p <- plot_self_complementarity(s)
  expect_s3_class(p, "ggplot")
  b <- strsplit(s, "")[[1L]]
  expect_equal(nrow(p$data), sum(outer(b, chartr("ACGT", "TGCA", b), "==")))
})
