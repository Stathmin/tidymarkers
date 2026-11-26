# =============================================================================
#  Y7  --  the Ye 2001 tetra-primer ARMS-PCR set: 12 primers, 3 assays
# =============================================================================
#  `ye2001_arms_primer_set.csv` is the package's only source of an
#  author-printed Tm for a primer whose allele-specific 3' base AND deliberate
#  second mismatch are both known.  It is the worked example behind the rule
#  already carried in `mismatch_discrimination.csv`.  This test is the
#  consumer gate for data-raw/ye2001_arms_primer_set.R.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

ye <- function() ext("ye2001_arms_primer_set.csv")

test_that("Y7: twelve primers, four roles per assay, three assays", {
  d <- ye()
  expect_equal(nrow(d), 12L)
  expect_equal(ncol(d), 10L)
  expect_equal(
    d$polymorphism,
    rep(c("TNF -308G>A", "IL6 -174G>C", "AGTR1 1166A>C"), each = 4L)
  )
  expect_equal(
    d$primer_role,
    rep(
      c("forward_inner", "reverse_inner", "forward_outer", "reverse_outer"),
      3L
    )
  )
})

test_that("Y7: sequences are real oligos, all distinct, inner >= 26 nt", {
  d <- ye()
  expect_true(all(grepl("^[ACGT]{20,40}$", d$sequence)))
  expect_false(anyDuplicated(d$sequence) > 0L)
  inner <- grepl("inner", d$primer_role)
  expect_true(all(nchar(d$sequence[inner]) >= 26L))
})

test_that("Y7: outer primers carry no allele, inner ones do", {
  d <- ye()
  expect_true(all(is.na(d$allele[grepl("outer", d$primer_role)])))
  expect_false(anyNA(d$allele[grepl("inner", d$primer_role)]))
  expect_setequal(d$allele[grepl("inner", d$primer_role)], c("A", "C", "G"))
})

test_that("Y7: per-assay conditions are inherited from the first row", {
  d <- ye()
  first <- d[!duplicated(d$polymorphism), ]
  expect_equal(first$touchdown, c(TRUE, TRUE, FALSE))
  expect_equal(first$annealing_c, c(63, 63, 58))
  expect_equal(first$mgcl2_mm, c(2.5, 2.5, 3.5))
  # every row of an assay repeats its assay's conditions
  expect_equal(length(unique(d[c("polymorphism", "annealing_c")])$annealing_c), 3L)
})

test_that("Y7: the two inner amplicons of an assay are different sizes", {
  d <- ye()
  # this is the whole point of a tetra-primer design: the two allele-specific
  # inner products must be separable on a gel
  for (p in unique(d$polymorphism)) {
    inner_bp <- d$amplicon_bp[d$polymorphism == p][1:2]
    expect_false(anyNA(inner_bp), info = p)
    expect_equal(length(unique(inner_bp)), 2L, info = p)
  }
  expect_equal(
    d$amplicon_bp,
    c(154, 224, 323, NA, 205, 176, 326, NA, 251, 224, 421, NA)
  )
  expect_equal(
    d$product,
    c(
      "A allele", "G allele", "outer pair", NA,
      "G allele", "C allele", "outer pair", NA,
      "A allele", "C allele", "outer pair", NA
    )
  )
})

test_that("Y7: the author Tm column is the paper's six values", {
  d <- ye()
  expect_equal(sort(unique(d$tm_c)), c(62, 64, 65, 67, 68, 72))
})

test_that("Y7: our NN Tm tracks the author's Tm to within 5 C", {
  d <- ye()
  skip_if_not_installed("TmCalculator")
  ours <- mapply(
    function(s, mg) {
      primer_tm(
        s,
        conditions = reaction_conditions(Na = 50, Mg = mg, dNTPs = 200),
        nn_table = "DNA_NN_Breslauer_1986"
      )
    },
    d$sequence,
    d$mgcl2_mm
  )
  delta <- ours - d$tm_c
  # measured 2026-11-26 with DNA_NN_Breslauer_1986 + Schildkraut2010 salt
  # at 50 mM Na / 200 uM dNTP / the assay's own Mg:
  #   delta range -3.28 .. +4.51, mean +0.65, median +0.85, mean |delta| 2.06,
  #   Pearson r = 0.81
  # The paper used Breslauer NN + Rychlik's formula and printed whole degrees,
  # so a ~2 C agreement is the honest expectation; a systematic bias would
  # show up as a large mean delta, which is what this pins.
  expect_true(all(abs(delta) < 5))
  expect_lt(abs(mean(delta)), 1.5)
  expect_lt(mean(abs(delta)), 2.5)
  expect_gt(stats::cor(ours, d$tm_c), 0.7)
})
