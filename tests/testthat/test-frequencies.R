# Frequency-aware footprints: a rare variant under a primer is reported
# (the load is counted) but must not doom the assay; a common one costs more,
# and the cost scales with the population frequency.

# A clean synthetic assay plus one extra variant inside the AS footprint.
freq_fixture <- function(freq, kind = c("snp", "indel"), seed = 4L, n = 500L, pos = 250L) {
  kind <- match.arg(kind)
  W <- withr::with_seed(seed, paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- as.character(Biostrings::subseq(genome[["Chr1"]], pos, pos))
  snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
  base <- suppressWarnings(assemble_kasp_assay(genome, snp, flank = 150L, cores = 1L))
  aa <- base$assay[1L, ]
  vpos <- aa$as_start + 10L # comfortably inside the AS footprint (mid-body)
  b <- substring(W, vpos, vpos)
  alt <- if (kind == "snp") setdiff(c("A", "C", "G", "T"), b)[1L] else paste0(b, "TG")
  extra <- tibble::tibble(
    chrom = "Chr1", pos = vpos, ref = b, alt = alt, pop_freq = freq,
    consequence = "missense_variant"
  )
  variants <- dplyr::bind_rows(tibble::tibble(chrom = "Chr1", pos = pos, consequence = "missense_variant"), extra)
  list(
    assay = suppressWarnings(assemble_kasp_assay(genome, snp, variants = variants, flank = 150L, cores = 1L)),
    vpos = vpos
  )
}

test_that("a rare variant under a primer is reported but not fatal", {
  skip_if_no_oligoarrayaux()
  for (kind in c("snp", "indel")) {
    m <- freq_fixture(0.01, kind)
    expect_false(inherits(m$assay, "marker_doomed"))
    a <- m$assay$assay[1L, ]
    expect_equal(a$as_n_var, 1L) # the variant IS seen under the footprint
    expect_lt(a$as_foot_dropout, KASP_LIMITS$footprint_dropout)
    # ... and it is surfaced in the acceptance checks, not hidden
    ch <- marker_design_checks(m$assay)
    expect_gt(ch$observed[ch$criterion == "AS footprint variants (n)"], 0)
  }
})

test_that("the footprint cost scales with the population frequency", {
  skip_if_no_oligoarrayaux()
  rare <- freq_fixture(0.01)
  common <- freq_fixture(1.00)
  expect_false(inherits(rare$assay, "marker_doomed"))
  expect_false(inherits(common$assay, "marker_doomed")) # a single mid-body variant is still tolerable
  expect_lt(rare$assay$assay$as_foot_dropout[1L], common$assay$assay$as_foot_dropout[1L])
})

# A weight is a proportion before it scores anything (ADR 0086).

test_that("a frequency outside [0, 1] is refused, not used as a weight", {
  expect_error(
    tidymarkers:::check_variant_weights(c(0.5, 5)),
    "\\[0, 1\\]"
  )
  expect_error(
    tidymarkers:::check_variant_weights(c(0.5, -0.2)),
    "\\[0, 1\\]"
  )
  expect_error(tidymarkers:::check_variant_weights("5%"), "\\[0, 1\\]")
  expect_equal(tidymarkers:::check_variant_weights(c(0.5, NA)), c(0.5, 1))
  expect_null(tidymarkers:::check_variant_weights(NULL))
})

test_that("an unknown weight counts the variant once, not as an NA load", {
  vpos <- c(10L, 20L, 30L)
  iv <- IRanges::IRanges(5, 40)
  expect_equal(
    tidymarkers:::variant_load_in(iv, vpos, c(0.5, NA, 1)),
    tidymarkers:::variant_load_in(iv, vpos, c(0.5, 1, 1))
  )
})

test_that("an unknown weight still costs the dropout it would with weight 1", {
  seq <- "ACGTACGTACGTACGTTTGGCCAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGG"
  tmpl <- list(seq = Biostrings::DNAString(seq))
  p <- substr(seq, 1L, 20L)
  zone <- function(w) {
    tidymarkers:::zone_dropout(
      p, 1L, 40L, TRUE, 1L, tmpl, 1L, 30L, integer(0), "A", w
    )
  }
  expect_gt(zone(1), 0)
  expect_equal(zone(NA_real_), zone(1))
})

test_that("a percent-scale pop_freq is refused on the design path", {
  skip_if_no_oligoarrayaux()
  expect_error(freq_fixture(5), "\\[0, 1\\]")
})
