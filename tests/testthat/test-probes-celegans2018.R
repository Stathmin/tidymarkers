# =============================================================================
#  probes/celegans2018  --  a -2 mismatch the sequencer can SEE in the product
# =============================================================================
#  microPublication Biology 2018, PMC7255807. The paper has no numeric
#  measurement at all: it reports that the -2 mismatch base introduced into
#  the allele-specific primers of a tetra-primer ARMS assay shows up in the
#  SANGER-sequenced amplicon of samples that do not carry it genomically.
#
#  That is the only published observation we have that a deliberate primer
#  mismatch can end up inside the product, which is the mechanism by which a
#  designed mismatch could erode discrimination as cycles go by. But the paper
#  gives no rate, no peak heights and no cycle series, so the extract is
#  evidence_tier = "weak" and the tests below refuse to let anything treat it
#  as a measured rate.

ext_probe <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA"), progress = FALSE
  )
}

test_that("probes/celegans2018: four primers, and the mismatch is at -2", {
  p <- ext_probe("celegans2018_arms_primers.csv")
  expect_equal(nrow(p), 4L)
  expect_setequal(p$primer, c(
    "flanking_forward", "flanking_reverse",
    "wt_reverse", "n3441_forward"
  ))
  as <- p[!is.na(p$mismatch_added), ]
  expect_equal(nrow(as), 2L)
  expect_setequal(as$mismatch_added, c("GC", "TG"))
  expect_true(all(as$mismatch_offset == -2L))
  # the mismatch dinucleotide is the last two bases before the 3' base, which
  # is what offset -2 MEANS; a hand-edited offset would break this. The pair
  # is written template-first, so the n3441 primer's "TG" is read as "GT" on
  # the primer strand.
  primer_pair <- tolower(c(GC = "gc", TG = "gt"))[as$mismatch_added]
  expect_equal(
    substr(as$sequence, nchar(as$sequence) - 2L, nchar(as$sequence) - 1L),
    unname(primer_pair)
  )
})

test_that("probes/celegans2018: each allele-specific primer carries its allele", {
  p <- ext_probe("celegans2018_arms_primers.csv")
  wt <- p[p$primer == "wt_reverse", ]
  n3 <- p[p$primer == "n3441_forward", ]
  # the WT primer is a REVERSE primer, so its 3' base is the complement of G
  expect_equal(substr(wt$sequence, nchar(wt$sequence), nchar(wt$sequence)), "c")
  expect_equal(substr(n3$sequence, nchar(n3$sequence), nchar(n3$sequence)), "a")
  expect_equal(wt$primer_3p_base, "c")
  expect_equal(n3$primer_3p_base, "a")
})

test_that("probes/celegans2018: nothing here is presented as a measured rate", {
  p <- ext_probe("celegans2018_arms_primers.csv")
  expect_true(all(p$evidence_tier == "weak"))
  # no numeric column other than the design parameters
  expect_false(any(grepl("rate|fraction|pct|percent", names(p))))
  obs <- unique(stats::na.omit(p$observation))
  expect_equal(length(obs), 1L)
  expect_true(grepl("suggest the possibility", obs, fixed = TRUE))
  expect_true(grepl("incorporating the primer sequence", obs, fixed = TRUE))
  # the outcome definition must refuse the cycle-series reading explicitly
  expect_true(all(grepl("NOT supported", p$outcome_definition, fixed = TRUE)))
})

test_that("probes/celegans2018: the design parameters are the paper's", {
  p <- ext_probe("celegans2018_arms_primers.csv")
  expect_true(all(p$cycles == 30L))
  expect_true(all(p$annealing_c == 61L))
  expect_equal(sort(unique(stats::na.omit(p$product_bp))), c(431L, 534L))
  expect_true(all(p$flanking_product_bp == 927L))
  # the -2 site is one of the -1..-3 positions the authors cite
  expect_true(all(p$mismatch_positions_cited[!is.na(p$mismatch_offset)] == "-1 to -3"))
})

test_that("probes/celegans2018: the sequences are still in the fetched paper", {
  f <- test_path("../../.tmp/lit/celegans2018/fulltext.xml")
  skip_if_not(file.exists(f), "celegans2018 paper not fetched")
  txt <- paste(readLines(f, warn = FALSE), collapse = "")
  p <- ext_probe("celegans2018_arms_primers.csv")
  for (s in p$sequence) {
    expect_true(grepl(s, txt, fixed = TRUE), info = s)
  }
  expect_true(grepl("-2 base", txt, fixed = TRUE))
  expect_true(grepl("30 cycles", txt, fixed = TRUE))
})
