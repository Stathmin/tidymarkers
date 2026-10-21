# An independent Tm oracle for `primer_tm()`: what the anchors are, and what
# they are not.
#
# WHAT IS HERE
#   * Same-lineage agreement: SantaLucia (2004, the package default) vs Allawi
#     (1998), through TmCalculator, with an explicit tolerance per model pair.
#     The numbers live in helper-tm-tolerance.R and are keyed on the pair, not
#     a single constant (see ROADMAP ZW3).
#   * Cross-lineage structure, asserted per table: Sugimoto - SantaLucia is
#     positive on the unremarkable 20-mers (a one-sided floor, min 2.2 C), and
#     Breslauer - SantaLucia is asserted only on the GC-rich subset, where the
#     gap actually is (min 9.5 C).
#   * Structural anchors with NO tolerance: Tm rises with GC at fixed length,
#     with length at fixed pattern, and is symmetric under reverse complement.
#   * GC-edge and homopolymer probes (GC < 35 %, GC > 65 %, 5+ base runs) inside
#     the compared set, so the per-pair tolerance also holds on them.
#   * Primer3 `oligotm()` as a THIRD oracle, when it happens to be installed:
#     a three-way matrix is printed and written to dev/logs/ (gitignored), and
#     the three-way disagreement is reported, never averaged or asserted
#     against a constant.
#
# WHAT IS NOT HERE
#   * No external EXPERIMENTAL Tm anchor.  There is none in this repository:
#     von Ahsen 2001 is closed access, dnaMATE is empty, Huguet 2017's figure is
#     raster-only.  These were checked and must not be re-attempted, and nothing
#     here fetches from the network.  Every anchor below compares one
#     thermodynamic model to another; agreement between models is not evidence
#     about a real duplex, and no dPCR Tm window should be calibrated on it.
#   * No absolute Tm values: nothing asserts "this oligo melts at X C", because
#     such a number would be a property of this machine's TmCalculator.
#   * No OligoArrayAux dG cross-check (a different quantity; covered by
#     test-oligo-dg.R).

# The main oligo set: 24 pseudo-random 20-mers.  Deterministic -- a fixed seed
# and no bare `sample()` at the top level -- so a failure names a stable oligo.
tm_oligos <- function(n = 24L, length = 20L) {
  set.seed(20040301L)
  vapply(seq_len(n), function(i) {
    paste(sample(c("A", "C", "G", "T"), length, replace = TRUE), collapse = "")
  }, character(1))
}

# The compared set: the main set plus the edge probes (defined once, in the
# helper).  Named, so a failure points at the probe that broke.
tm_oracle_set <- function() {
  s <- c(tm_oligos(), tm_edge_oligos())
  rlang::set_names(unname(s), c(paste0("oligo", seq_along(tm_oligos())), names(tm_edge_oligos())))
}

# The longest run of one base in an oligo, and its terminal 3' run.
max_run <- function(s) max(rle(strsplit(s, "")[[1L]])$lengths)
three_prime_run <- function(s) {
  r <- rle(strsplit(s, "")[[1L]])
  r$lengths[[length(r$lengths)]]
}

tm_table <- function(seqs, nn_table) {
  rlang::set_names(primer_tm(seqs, nn_table = nn_table), names(seqs))
}

# Is Primer3's oligotm on PATH?  Not a package dependency and not installed on
# this machine, so the third oracle skips rather than failing a CRAN check.
has_oligotm <- function() nzchar(Sys.which("oligotm"))

skip_if_no_oligotm <- function() {
  if (!has_oligotm()) {
    testthat::skip("Primer3 oligotm not available (binary 'oligotm' not on PATH)")
  }
}

# `oligotm -tp 1 -sc 1 -n 0 <seq>`: Primer3's SantaLucia (2004) NN parameters,
# salt-corrected, 0.05 M monovalent -- the same model family `primer_tm()`
# defaults to, but a separate implementation and salt code.
oligotm_tm <- function(seqs) {
  vapply(seqs, function(s) {
    out <- system2("oligotm", c("-tp", "1", "-sc", "1", "-n", "0", s), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0L) {
      rlang::abort(paste0("oligotm failed on ", s, ": ", paste(out, collapse = " ")))
    }
    as.numeric(out[[1L]])
  }, numeric(1))
}

# --- the tolerance table itself ---------------------------------------------

test_that("the tolerance table is well formed and order-insensitive", {
  tol <- tm_tolerances()
  expect_s3_class(tol, "tbl_df")
  expect_named(tol, c("model_a", "model_b", "max_delta_c", "basis"))
  expect_true(all(is.finite(tol$max_delta_c)))
  expect_true(all(tol$max_delta_c > 0))
  # a number without its provenance is not a tolerance anyone can audit
  expect_true(all(nzchar(tol$basis)))
  expect_equal(nrow(dplyr::distinct(tol, .data$model_a, .data$model_b)), nrow(tol))

  fwd <- tm_tolerance(NN_SANTA_LUCIA, NN_ALLAWI)
  rev <- tm_tolerance(NN_ALLAWI, NN_SANTA_LUCIA)
  expect_equal(fwd$max_delta_c, rev$max_delta_c)
  expect_error(tm_tolerance(NN_SANTA_LUCIA, NN_SANTA_LUCIA), "no Tm tolerance")
  expect_error(tm_tolerance(NN_ALLAWI, NN_BRESLAUER), "no Tm tolerance")
})

test_that("primer_tm()'s default table is the one the oracle pins", {
  # the default is the symbol NN_TABLE, not a literal, so compare by name
  expect_identical(deparse(formals(primer_tm)$nn_table), "NN_TABLE")
  expect_identical(tidymarkers:::NN_TABLE, NN_SANTA_LUCIA)
})

test_that("every pinned model is a table primer_tm() accepts", {
  s <- "ACGTACGTACGTACGTACGT"
  for (model in tm_models()) {
    expect_true(is.finite(primer_tm(s, nn_table = model)), info = model)
  }
})

# --- same lineage: an explicit tolerance ------------------------------------

test_that("SantaLucia vs Allawi agree within the pinned tolerance", {
  seqs <- tm_oracle_set()
  sl <- tm_table(seqs, NN_SANTA_LUCIA)
  al <- tm_table(seqs, NN_ALLAWI)
  delta <- abs(al - sl)
  tol <- tm_tolerance(NN_SANTA_LUCIA, NN_ALLAWI)$max_delta_c
  testthat::expect_lte(max(delta), tol)
  expect_named(delta, names(seqs))
})

test_that("the same tolerance holds on the GC-edge and homopolymer probes", {
  seqs <- tm_edge_oligos()
  # the set must actually reach the edges it claims, or the test is vacuous
  expect_lt(min(gc_content(seqs)), 35)
  expect_gt(max(gc_content(seqs)), 65)
  expect_gte(max(vapply(seqs, max_run, integer(1L))), 5L)
  # a run of the same base at the 3' end, which is what "3' edge" means
  expect_gte(three_prime_run(seqs[["three_prime_a"]]), 5L)
  expect_gte(max_run(seqs[["three_prime_g"]]), 5L)

  delta <- abs(tm_table(seqs, NN_ALLAWI) - tm_table(seqs, NN_SANTA_LUCIA))
  tol <- tm_tolerance(NN_SANTA_LUCIA, NN_ALLAWI)$max_delta_c
  testthat::expect_lte(max(delta), tol)
})

# --- cross lineage: structure per table, never a blanket rule ---------------

test_that("Sugimoto sits above SantaLucia on the unremarkable 20-mers", {
  seqs <- tm_oligos()
  sl <- tm_table(seqs, NN_SANTA_LUCIA)
  su <- tm_table(seqs, NN_SUGIMOTO)
  delta <- su - sl
  floor_c <- tm_tolerance(NN_SUGIMOTO, NN_SANTA_LUCIA)$max_delta_c

  expect_true(all(delta > 0))
  expect_gt(mean(delta), floor_c)
  # the same table, a different order: it must rank the oligos alike
  expect_gt(stats::cor(su, sl, method = "spearman"), 0.6)
})

test_that("Breslauer's gap over SantaLucia is asserted only where it exists", {
  seqs <- tm_oracle_set()
  sl <- tm_table(seqs, NN_SANTA_LUCIA)
  delta <- tm_table(seqs, NN_BRESLAUER) - sl
  tol <- tm_tolerance(NN_BRESLAUER, NN_SANTA_LUCIA)$max_delta_c

  gc_rich <- names(seqs)[gc_content(seqs) >= 65]
  expect_gt(length(gc_rich), 0L)
  expect_true(all(delta[gc_rich] > tol), info = paste(gc_rich, collapse = ", "))
  expect_true(all(delta > 0))

  # A/T-rich: the gap is smaller here than on GC-rich probes, but not
  # uniformly -- an A homopolymer still reaches 10.9 C -- so a claim of the
  # GC-rich size would be false, which is exactly why it is not made.
  at_rich <- names(seqs)[gc_content(seqs) < 35]
  expect_gt(length(at_rich), 0L)
  expect_true(all(delta[at_rich] > 0), info = paste(at_rich, collapse = ", "))
  expect_lt(median(delta[at_rich]), tol)
})

# --- structural anchors: no tolerance ---------------------------------------

test_that("Tm rises with GC at a fixed length", {
  # the same A/G pattern with a swapped C/G every other base: composition
  # changes, length and the alternating skeleton do not
  seqs <- c(
    "ATATATATATATATATATAT",
    "GTGTGTGTGTGTGTGTGTGT",
    "GCGCGCGCGCGCGCGCGCGC"
  )
  tm <- primer_tm(seqs)
  expect_equal(nchar(seqs), rep(20L, 3L))
  expect_true(all(diff(tm) > 0), info = paste(round(tm, 2), collapse = ", "))
})

test_that("Tm rises with length at a fixed composition", {
  # lengthening "GGCATGCATGCA..." by one base preserves only length, so a run
  # of the same base is used: 12..20 A's, then 12..20 C's.  Composition is
  # constant within each run and the two runs disagree on GC, which is the
  # point: it is length, not base identity, that has to move the Tm.
  runs <- c(
    lapply(12:20, function(n) strrep("A", n)),
    lapply(12:20, function(n) strrep("C", n))
  )
  tm_a <- primer_tm(unlist(runs[1:9]))
  tm_c <- primer_tm(unlist(runs[10:18]))
  expect_true(all(diff(tm_a) > 0), info = paste(round(tm_a, 2), collapse = ", "))
  expect_true(all(diff(tm_c) > 0), info = paste(round(tm_c, 2), collapse = ", "))
})

test_that("Tm does not depend on the strand", {
  seqs <- tm_oracle_set()
  rc <- vapply(seqs, function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAString(s))), character(1))
  for (model in tm_models()) {
    expect_equal(tm_table(seqs, model), rlang::set_names(primer_tm(rc, nn_table = model), names(seqs)),
      tolerance = 1e-6, info = model
    )
  }
})

# --- third oracle: Primer3 --------------------------------------------------

test_that("Primer3 oligotm agrees with both TmCalculator tables within tolerance", {
  skip_if_no_oligotm()
  seqs <- tm_oracle_set()
  sl <- tm_table(seqs, NN_SANTA_LUCIA)
  al <- tm_table(seqs, NN_ALLAWI)
  p3 <- rlang::set_names(oligotm_tm(seqs), names(seqs))

  tol_sl <- tm_tolerance(NN_SANTA_LUCIA, NN_ALLAWI)$max_delta_c
  testthat::expect_lte(max(abs(p3 - sl)), tol_sl)
  testthat::expect_lte(max(abs(p3 - al)), tol_sl)
})

test_that("the three oracles are reported as a matrix, not averaged", {
  skip_if_no_oligotm()
  seqs <- tm_oracle_set()
  matrix <- tibble::tibble(
    oligo = names(seqs),
    seq = unname(seqs),
    gc = gc_content(seqs),
    santa_lucia = unname(tm_table(seqs, NN_SANTA_LUCIA)),
    allawi = unname(tm_table(seqs, NN_ALLAWI)),
    primer3 = unname(oligotm_tm(seqs))
  )
  matrix$spread <- apply(matrix[c("santa_lucia", "allawi", "primer3")], 1L, max) -
    apply(matrix[c("santa_lucia", "allawi", "primer3")], 1L, min)

  print(matrix, n = nrow(matrix))
  log_dir <- file.path("dev", "logs")
  if (dir.exists(log_dir)) {
    utils::write.csv(matrix, file.path(log_dir, "tm-oracle-three-way.csv"), row.names = FALSE)
  }
  # the report itself is the assertion: every oligo carries three numbers
  expect_false(anyNA(matrix))
  expect_true(all(matrix$spread >= 0))
})
