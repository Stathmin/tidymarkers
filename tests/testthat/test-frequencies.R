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
