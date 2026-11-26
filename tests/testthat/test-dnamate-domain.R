# =============================================================================
#  dnaMATE  --  the domain its Tm model is valid over, and the salt term we share
# =============================================================================
#  `dnamate_operating_domain.csv` is the consumer gate for
#  data-raw/dnamate_operating_domain.R.  dnaMATE publishes no Tm table, so the
#  only checkable claim it supports is the boundary of its own model: 16-30 nt,
#  50-600 mM Na+, 89 % of predictions within 5 C inside that window, the
#  Schildkraut-Lifson salt factor at 16.6 log[Na+], the b strand-count term, and
#  the three nearest-neighbour tables it averages.  The algebraic test is the
#  one that ties the paper to `primer_tm()`: the package's `Schildkraut2010`
#  salt method must move a Tm by exactly 16.6 log10([Na+]2/[Na+]1).

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

dnamate <- function() ext("dnamate_operating_domain.csv")

test_that("dnaMATE: 13 operating-domain rows, unique items, all cite the paper", {
  d <- dnamate()
  expect_equal(nrow(d), 13L)
  expect_equal(ncol(d), 7L)
  expect_false(anyDuplicated(d$item) > 0L)
  expect_true(all(nzchar(d$reported)))
  expect_true(all(d$doi == "10.1093/nar/gki379"))
  expect_true(all(d$pmcid == "PMC1160140"))
  expect_setequal(unique(d$section), c("Abstract", "Methods", "Results"))
})

test_that("dnaMATE: the numeric domain parses to the values the paper prints", {
  d <- dnamate()
  num <- function(item) as.numeric(d$value[d$item == item])
  expect_equal(num("sequence_length_min"), 16)
  expect_equal(num("sequence_length_max"), 30)
  expect_equal(num("monovalent_salt_min"), 50)
  expect_equal(num("monovalent_salt_max"), 600)
  expect_equal(num("accuracy_within_5C"), 89)
  expect_equal(num("server_batch_limit"), 5000)
  expect_equal(num("salt_adjustment_coefficient"), 16.6)
  expect_equal(num("gas_constant"), 1.987)
  expect_equal(num("b_nonself_complementary"), 4)
  expect_equal(num("b_self_complementary_or_excess"), 1)
  expect_equal(num("consensus_zone_error_bound"), 5)
})

test_that("dnaMATE: the three NN tables it averages are all reachable here", {
  d <- dnamate()
  tables <- strsplit(d$value[d$item == "nn_tables"], ";", fixed = TRUE)[[1L]]
  expect_setequal(tables, c("Breslauer", "SantaLucia", "Sugimoto"))
  available <- eval(formals(TmCalculator::tm_nn)$nn_table)
  expect_true(all(c(
    "DNA_NN_Breslauer_1986", "DNA_NN_SantaLucia_2004", "DNA_NN_Sugimoto_1996"
  ) %in% available))
  # the package's own default is one of them
  expect_equal(NN_TABLE, "DNA_NN_SantaLucia_2004")
})

test_that("dnaMATE: the package applies the paper's salt rule, not another one", {
  expect_equal(TM_SALT_METHOD, "Schildkraut2010")
})

test_that("dnaMATE: 50 -> 1000 mM Na+ shifts Tm by exactly 16.6 log10(0.05)", {
  seq <- "GAAGGTGACCAAGTTCATGCT"
  low <- primer_tm(seq, conditions = reaction_conditions(Na = 50, Mg = 0, dNTPs = 0))
  high <- primer_tm(seq, conditions = reaction_conditions(Na = 1000, Mg = 0, dNTPs = 0))
  expect_equal(high - low, 16.6 * log10(1000 / 50), tolerance = 1e-6)
  # the same delta must survive a Mg/dNTP-containing mix: the salt term is
  # taken on the monovalent equivalent, and Mg must not double-count it
  low_mg <- primer_tm(seq, conditions = reaction_conditions(Na = 50, Mg = 2, dNTPs = 200))
  high_mg <- primer_tm(seq, conditions = reaction_conditions(Na = 1000, Mg = 2, dNTPs = 200))
  expect_equal(high_mg - low_mg, 16.6 * log10(1000 / 50), tolerance = 1e-6)
})

test_that("dnaMATE: outside the 16-30 nt window no accuracy claim is carried", {
  d <- dnamate()
  lo <- as.numeric(d$value[d$item == "sequence_length_min"])
  hi <- as.numeric(d$value[d$item == "sequence_length_max"])
  expect_lt(lo, hi)
  # the paper's own wording is boundary-qualified; the extract must keep it
  acc <- d$reported[d$item == "accuracy_within_5C"]
  expect_match(acc, "89%", fixed = TRUE)
  expect_match(acc, "<5", fixed = TRUE)
  # and the consensus-map caveat must be tied to the no-consensus zone, so a
  # caller cannot quote 89 % as if it held everywhere
  expect_match(d$value[d$item == "black_zone"], "no consensus", fixed = TRUE)
})
