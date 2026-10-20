# Every marker object represents itself: `plot(marker)` dispatches on its class.
# The KASP assay is tested on a real design; the other chemistries use a minimal
# but structurally faithful stub, so this file stays cheap.

stub_snp <- function() {
  tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "G")
}

stub_base <- function() {
  list(
    snp = stub_snp(),
    amplicon = tibble::tibble(start = 250L, end = 350L),
    primers = list(
      forward = "ACGTACGTACGTACGTACGT", reverse = "TGCATGCATGCATGCATGCA",
      product = 100L, start = 250L, end = 350L
    ),
    selection = list(relaxed = character(), violated = character())
  )
}

test_that("the KASP assay plots itself", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  a <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp, flank = 150L, cores = 1L))
  expect_s3_class(plot(a), "patchwork")
})

test_that("each non-KASP marker class plots its own panel", {
  snp <- stub_snp()
  base <- stub_base()
  amp <- base$amplicon
  prim <- base$primers

  stub <- function(cls, ...) {
    out <- base
    fields <- list(...)
    for (nm in names(fields)) out[[nm]] <- fields[[nm]]
    structure(out, class = cls)
  }

  assay_stub <- stub("marker_assay",
    assay = tibble::tibble(
      as_ref = "ACGTACGTACGTACGTACGT", as_alt = "ACGTACGTACGTACGTACGA",
      com = "TGCATGCATGCATGCATGCA",
      as_start = 281L, as_end = 300L, as_alt_start = 281L, as_alt_end = 300L,
      com_start = 320L, com_end = 340L, strand = "+"
    ),
    scheme = "snp"
  )

  ssr_obj <- stub("marker_ssr")
  ssr_obj[["repeat"]] <- tibble::tibble(start = 295L, end = 305L, unit = 2L, copies = 5L)

  objects <- list(
    taqman = stub("marker_taqman", probes = tibble::tibble(start = 290L, end = 310L)),
    dpcr = stub("marker_dpcr", probes = tibble::tibble(start = 290L, end = 310L)),
    snapshot = stub("marker_snapshot", extension = list(start = 281L, end = 300L, strand = "+")),
    lamp = stub("marker_lamp",
      primers = list(
        F3 = "ACGTACGTACGTACGTACGTA", B3 = "TGCATGCATGCATGCATGCAT",
        FIP = strrep("ACGT", 10L), BIP = strrep("TGCA", 10L)
      ),
      loop = tibble::tibble(start = 305L, end = 335L)
    ),
    rpa = stub("marker_rpa", probe = NULL),
    crispr = stub("marker_crispr",
      cas = "cas12a", strand = "+", guide = strrep("ACGT", 5L),
      pam = "TTTA", guide_start = 281L, guide_end = 300L, off_target = 0L
    ),
    sanger = stub("marker_sanger",
      seq_primers = list(forward = "AAAACCCCGGGGTTTT", reverse = "TTTTGGGGCCCCAAAA", nested = TRUE)
    ),
    ssr = ssr_obj,
    indel_size = stub("marker_indel_size",
      amplicon = tibble::tibble(start = 250L, end = 350L, ref_size = 100L, alt_size = 110L),
      gap = 10L
    ),
    duplex = stub("marker_duplex", assays = list(assay_stub, assay_stub), readout = "co-dominant")
  )

  for (nm in names(objects)) {
    p <- plot(objects[[nm]])
    expect_true(inherits(p, c("patchwork", "ggplot")), info = nm)
  }
})

test_that("marker_design_checks evaluates a non-KASP marker", {
  x <- structure(
    list(
      snp = stub_snp(),
      amplicon = tibble::tibble(start = 250L, end = 350L),
      primers = list(
        forward = "ACGTACGTACGTACGTACGT", reverse = "TGCATGCATGCATGCATGCA"
      ),
      selection = list(relaxed = character(), violated = character())
    ),
    class = "marker_caps"
  )
  ch <- marker_design_checks(x, scheme = "caps")
  expect_equal(nrow(ch), 13L)
  expect_true(all(c("criterion", "observed", "limit", "pass", "level") %in% names(ch)))
  expect_true("amplicon length (bp)" %in% ch$criterion)
  expect_s3_class(marker_decision_graph(x, scheme = "caps"), "tbl_graph")
})

test_that("Tm guidelines are method-specific, not KASP gates", {
  expect_equal(marker_tm_window("taqman", "probe_REF"), c(65, 72))
  expect_equal(marker_tm_window("taqman", "forward"), c(58, 62))
  expect_equal(marker_tm_window("lamp", "FIP"), c(58, 65))
  expect_null(marker_tm_window("rpa", "forward"))
  expect_null(marker_tm_window("crispr", "guide"))
})

test_that("non-KASP checks report one observed Tm per primer", {
  skip_if_no_oligoarrayaux()
  x <- structure(
    list(
      snp = stub_snp(),
      amplicon = tibble::tibble(start = 250L, end = 350L),
      primers = list(
        forward = "ACGTACGTACGTACGTACGT", reverse = "TGCATGCATGCATGCATGCA"
      )
    ),
    class = "marker_caps"
  )
  checks <- marker_primer_checks(x, "caps")
  tm <- checks[grepl("^Tm ", checks$name), ]
  expect_equal(nrow(tm), 2L)
  expect_setequal(tm$name, c("Tm F", "Tm R"))
  expect_true(all(tm$status %in% c("ok", "broken")))
  expect_true(all(nchar(tm$value) <= 4L))
})

test_that("the check column shows unknown and not_applicable, never ok", {
  x <- structure(
    list(
      snp = stub_snp(),
      amplicon = tibble::tibble(start = 250L, end = 350L),
      primers = list(
        forward = "ACGTACGTACGTACGTACGT", reverse = "TGCATGCATGCATGCATGCA"
      )
    ),
    class = "marker_caps"
  )
  checks <- marker_primer_checks(x, "caps")
  # no genome/ideogram -> the off-target count was not computed: unknown
  off <- checks[checks$name == "off-target", ]
  expect_equal(nrow(off), 1L)
  expect_equal(off$status, "unknown")
  expect_false(any(checks$status == "ok" & checks$name == "off-target"))
  # RPA has no Tm guideline, so its Tm rows are not_applicable, not neutral-ok
  class(x) <- "marker_rpa"
  skip_if_no_oligoarrayaux()
  rpa <- marker_primer_checks(x, "rpa")
  tm <- rpa[grepl("^Tm ", rpa$name), ]
  expect_true(nrow(tm) >= 1L)
  expect_true(all(tm$status == "not_applicable"))
  # the plot adapter maps both to neutral (never green) and still draws
  p <- plot_marker_checks(x, "rpa")
  expect_s3_class(p, "ggplot")
  expect_setequal(p$data$status[p$data$name %in% tm$name], "neutral")
})

test_that("a fragment readout names the repeat formula and the observed sizes", {
  x <- structure(
    list(
      snp = stub_snp(),
      amplicon = tibble::tibble(start = 250L, end = 349L),
      primers = list(forward = "ACGT", reverse = "TGCA", product = 100L)
    ),
    class = "marker_ssr"
  )
  x[["repeat"]] <- tibble::tibble(
    start = 295L, end = 304L, unit = 2L, copies = 5L,
    motif = "ACACACACAC", unit_seq = "AC"
  )
  p <- plot_fragment_readout(x, "ssr")
  expect_s3_class(p, "ggplot")
  expect_true(any(grepl("(AC)5", p$data$label, fixed = TRUE)))
  expect_true(any(grepl("(AC)6", p$data$label, fixed = TRUE)))
  expect_setequal(p$data$size, c(100L, 102L))
})

test_that("a doomed marker still plots its target context", {
  d <- structure(list(
    snp = stub_snp(),
    template = list(snp_i = 151L, seq = Biostrings::DNAString(strrep("ACGT", 100L))),
    reason = "no candidate", stage = "as",
    candidates = tibble::tibble()
  ), class = "marker_doomed")
  expect_true(inherits(plot(d), c("patchwork", "ggplot")))
})
