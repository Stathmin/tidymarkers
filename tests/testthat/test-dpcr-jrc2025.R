# dPCR assay anchor: 10.2807/1560-7917.es.2025.30.33.2500183 (PMC12372895).
#
# `inst/extdata/dpcr_tm_jrc2025.csv` carries the printed Table 1 of a real dPCR
# assay set: six oligo Tms (the only dPCR paper found that prints one) and the
# two amplicons they define.  The Tm column is the design tool's output, not a
# UV melt, so nothing here asserts a Tm model; what is asserted is the paper's
# own arithmetic and the design contrast it exposes.

ext <- function(name) {
  utils::read.csv(tm_extdata(name), stringsAsFactors = FALSE)
}

rc <- function(x) {
  b <- c(A = "T", C = "G", G = "C", T = "A")
  paste(rev(b[strsplit(x, "")[[1L]]]), collapse = "")
}

test_that("jrc2025: the shipped table holds 8 rows, six oligos and two amplicons", {
  d <- ext("dpcr_tm_jrc2025.csv")
  expect_equal(nrow(d), 8L)
  expect_setequal(unique(d$assay_id), c("jrc_ha", "jrc_mp"))
  expect_equal(sum(d$primer_or_probe == "amplicon"), 2L)
  expect_setequal(
    d$primer_or_probe,
    c("primer_fwd", "primer_rev", "probe_fam", "probe_vic", "amplicon")
  )
  expect_equal(sum(is.na(d$tm_measured_C)), 2L) # the amplicon rows measure no Tm
})

test_that("jrc2025: every printed Tm is the one the paper prints", {
  d <- ext("dpcr_tm_jrc2025.csv")
  oligo <- d[d$primer_or_probe != "amplicon", ]
  got <- stats::setNames(oligo$tm_measured_C, paste(oligo$assay_id, oligo$primer_or_probe))
  expect_equal(
    got,
    c(
      "jrc_ha primer_fwd" = 59.2, "jrc_ha primer_rev" = 60.2, "jrc_ha probe_fam" = 60.9,
      "jrc_mp primer_fwd" = 61.0, "jrc_mp primer_rev" = 61.9, "jrc_mp probe_vic" = 57.8
    )
  )
})

test_that("jrc2025: the amplicon is the span its two primers define", {
  d <- ext("dpcr_tm_jrc2025.csv")
  for (id in unique(d$assay_id)) {
    a <- d$sequence[d$assay_id == id & d$primer_or_probe == "amplicon"]
    fwd <- d$sequence[d$assay_id == id & d$primer_or_probe == "primer_fwd"]
    rev <- d$sequence[d$assay_id == id & d$primer_or_probe == "primer_rev"]
    probe <- d$sequence[d$assay_id == id & d$primer_or_probe %in% c("probe_fam", "probe_vic")]
    expect_length(a, 1L)
    expect_equal(nchar(a), unique(d$amplicon_size_bp[d$assay_id == id]))
    expect_true(startsWith(a, fwd))
    expect_true(endsWith(a, rc(rev)))
    expect_true(grepl(probe, a, fixed = TRUE))
  }
  expect_equal(unique(d$amplicon_size_bp[d$assay_id == "jrc_ha"]), 139L)
  expect_equal(unique(d$amplicon_size_bp[d$assay_id == "jrc_mp"]), 114L)
})

test_that("jrc2025: neither probe meets the +8-10 C probe rule, and the two land on opposite sides", {
  d <- ext("dpcr_tm_jrc2025.csv")
  probe_gap <- function(id) {
    sub <- d[d$assay_id == id, ]
    probe <- sub$tm_measured_C[grepl("^probe", sub$primer_or_probe)]
    expect_length(probe, 1L)
    probe - max(sub$tm_measured_C[grepl("^primer", sub$primer_or_probe)])
  }
  # measured: JRC-HA's probe is 0.7 C ABOVE its hotter primer, JRC-MP's is
  # 4.1 C BELOW its hotter primer -- neither is anywhere near +8-10 C
  expect_equal(probe_gap("jrc_ha"), 0.7, tolerance = 1e-9)
  expect_equal(probe_gap("jrc_mp"), -4.1, tolerance = 1e-9)
  expect_true(all(c(probe_gap("jrc_ha"), probe_gap("jrc_mp")) < 8))
})
