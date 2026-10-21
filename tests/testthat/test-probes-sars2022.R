# =============================================================================
#  probes/sars2022  --  a double-mismatch primer, measured as a Ct gap
# =============================================================================
#  Hale et al., J Clin Virol 2022, doi:10.1016/j.jcv.2021.105049, PMC8628602.
#  Two reactions whose allele-specific primers differ ONLY at the 3'-terminal
#  base, each also carrying a deliberate C four bases from the 3' end, read out
#  on 42 clinical samples.
#
#  The number Z7 wants is "how much margin does the second mismatch buy".
#  This paper gives one honest answer and one trap:
#    * honest: in the 4 samples with any cross-reaction the gap is 14.3-18.8
#      cycles, and it is 0 or infinite otherwise -- never a middle value;
#    * trap: "typically around 17 cycles" is not a typical sample. 38 of the 42
#      samples have NO Ct in the non-matching reaction, so their gap is
#      unbounded. The tests below keep the two apart.

ext_probe <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name),
    show_col_types = FALSE,
    na = c("", "NA"), progress = FALSE
  )
}

test_that("probes/sars2022: the two forward primers differ only at the 3' base", {
  p <- ext_probe("sars2022_dmas_primers.csv")
  expect_equal(nrow(p), 5L)
  a <- p$sequence[p$oligo == "non_delta_fwd"]
  b <- p$sequence[p$oligo == "delta_fwd"]
  expect_equal(nchar(a), 24L)
  expect_equal(nchar(b), 24L)
  expect_equal(substr(a, 1L, 23L), substr(b, 1L, 23L))
  expect_false(substr(a, 24L, 24L) == substr(b, 24L, 24L))
  expect_equal(substr(a, 24L, 24L), p$primer_3p_base[p$oligo == "non_delta_fwd"])
})

test_that("probes/sars2022: the artificial mismatch is a C at the 4th base", {
  p <- ext_probe("sars2022_dmas_primers.csv")
  as <- p[!is.na(p$artificial_mismatch_base), ]
  expect_equal(nrow(as), 4L)
  expect_true(all(as$artificial_mismatch_base == "C"))
  expect_true(all(as$artificial_mismatch_offset == -3L))
  # offset -3 means the 4th base counting the 3'-terminal base as 1
  pos <- nchar(as$sequence) - 3L
  expect_equal(substr(as$sequence, pos, pos), as$artificial_mismatch_base)
})

test_that("probes/sars2022: a gap exists only in the four cross-reacting samples", {
  d <- ext_probe("sars2022_dmas_ct.csv")
  expect_equal(nrow(d), 42L)
  with_gap <- d[d$nonmatching_reaction_detected, ]
  expect_setequal(with_gap$sample, c("1", "2", "20", "37"))
  expect_equal(nrow(d) - nrow(with_gap), 38L)
  # the gap is exactly Ct(non-matching) - Ct(matching), never a model value
  expect_equal(with_gap$ct_gap,
    as.numeric(with_gap$ct_nonmatching) - with_gap$ct_matching,
    tolerance = 1e-9
  )
  # and no gap is invented for the other samples
  expect_true(all(is.na(d$ct_gap[!d$nonmatching_reaction_detected])))
})

test_that("probes/sars2022: the paper's 14.3-18.8 range is the recomputed one", {
  d <- ext_probe("sars2022_dmas_ct.csv")
  g <- sort(stats::na.omit(d$ct_gap))
  expect_equal(length(g), 4L)
  expect_equal(min(g), 14.3)
  expect_equal(max(g), 18.8)
  # the sample the paper narrates: B = 20.0, A = 37.4
  s1 <- d[d$sample == "1", ]
  expect_equal(as.numeric(s1$ct_delta_reaction), 20.0)
  expect_equal(as.numeric(s1$ct_non_delta_reaction), 37.4)
  expect_equal(s1$ct_gap, 17.4)
})

test_that("probes/sars2022: a detected cross-reaction only happens at high load", {
  d <- ext_probe("sars2022_dmas_ct.csv")
  x <- d[d$nonmatching_reaction_detected, ]
  # the paper's own condition: E-gene Ct <= 22
  e <- as.numeric(x$e_gene_ct)
  expect_true(all(e <= 22))
  # every sample that did cross-react has a matching-reaction Ct below 22 too
  expect_true(all(x$ct_matching <= 22))
  # and the ones that did not cross-react are not simply the low-load ones
  y <- d[!d$nonmatching_reaction_detected & !is.na(as.numeric(d$e_gene_ct)), ]
  expect_gte(sum(as.numeric(y$e_gene_ct) <= 22), 5L)
})

test_that("probes/sars2022: the sequences are still in the fetched paper", {
  f <- test_path("../../.tmp/lit/sars2022/fulltext.xml")
  skip_if_not(file.exists(f), "sars2022 paper not fetched")
  suppressPackageStartupMessages(library(xml2))
  doc <- read_xml(f)
  txt <- paste(xml_text(xml_find_all(doc, "//table-wrap")), collapse = " ")
  p <- ext_probe("sars2022_dmas_primers.csv")
  for (s in p$sequence) {
    expect_true(grepl(s, gsub("[[:space:]]", "", txt), fixed = TRUE), info = s)
  }
  body <- paste(xml_text(xml_find_all(doc, "//body//p")), collapse = " ")
  expect_true(grepl("range 14.3 to 18.8 cycles", body, fixed = TRUE))
})
