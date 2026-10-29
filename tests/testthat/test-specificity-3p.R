# ---------------------------------------------------------------------------
# Y3.3: 3'-end-aware off-target grading (offtarget_profile + exact_hit_table
# profile columns).  One primer, three plus-strand sites laid out with 10-nt
# A-flanks and a 10-nt spacer between sites, so the sites start at 11 (exact),
# 52 (3'-terminal T->A) and 93 (inner pos-10 C->G).

p <- "GAAGGTGACCAAGTTCATGCT" # 21 nt, terminal base T
term <- paste0(substr(p, 1L, 20L), "A") # 3'-terminal T->A
inner <- paste0(substr(p, 1L, 9L), "G", substr(p, 11L, 21L)) # pos 10 C->G
g <- Biostrings::DNAStringSet(c(chr1 = paste0(
  strrep("A", 10L), p, strrep("A", 10L),
  strrep("A", 10L), term, strrep("A", 10L),
  strrep("A", 10L), inner, strrep("A", 10L)
)))

profile_cols <- c(
  "n_mismatch", "mismatch_positions", "three_prime_mismatch", "seed_mismatch",
  "terminal_mismatch", "three_prime_terminal", "penalty"
)
hit_cols <- c(
  "primer", "chrom", "start", "end", "width", "strand", profile_cols
)

test_that("exact_hit_table carries the seven profile fields", {
  h <- exact_hit_table(p, g)
  expect_equal(names(h), hit_cols)
  expect_equal(h$n_mismatch, 0L)
  expect_equal(h$penalty, 0)
})

test_that("offtarget_profile grades an exact, a 3'-terminal and an inner miss", {
  hits <- exact_hit_table(c(AS = p), g, max_mismatch = 1L)
  h <- offtarget_profile(hits, c(AS = p), g)
  expect_equal(nrow(h), 3L)
  expect_true(all(h$strand == "+"))
  expect_equal(h$start, c(11L, 52L, 93L))

  exact <- h[h$start == 11L, , drop = FALSE]
  expect_equal(exact$n_mismatch, 0L)
  expect_true(exact$tolerated)

  term_row <- h[h$start == 52L, , drop = FALSE]
  expect_true(term_row$three_prime_terminal)
  expect_equal(term_row$three_prime_mismatch, 1L)
  expect_equal(term_row$mismatch_positions[[1L]], 21L)
  expect_false(term_row$tolerated)

  inner_row <- h[h$start == 93L, , drop = FALSE]
  expect_equal(inner_row$n_mismatch, 1L)
  expect_equal(inner_row$three_prime_mismatch, 0L)
  expect_false(inner_row$three_prime_terminal)
  expect_equal(inner_row$mismatch_positions[[1L]], 10L)
  expect_true(inner_row$tolerated)
})

test_that("a minus-strand 3'-terminal miss is mirrored into the primer's frame", {
  g_minus <- Biostrings::DNAStringSet(c(chr1 = paste0(
    strrep("A", 10L),
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(term))),
    strrep("A", 10L)
  )))
  h <- offtarget_profile(
    exact_hit_table(c(AS = p), g_minus, max_mismatch = 1L), c(AS = p), g_minus
  )
  expect_equal(nrow(h), 1L)
  expect_equal(h$strand, "-")
  expect_true(h$three_prime_terminal)
  expect_equal(h$mismatch_positions[[1L]], 21L)
})

test_that("offtarget_profile keeps the column contract on zero rows", {
  hits <- exact_hit_table(c(AS = p), g, max_mismatch = 1L)
  h <- offtarget_profile(hits, c(AS = p), g)
  h0 <- offtarget_profile(h[0, ], c(AS = p), g)
  expect_equal(nrow(h0), 0L)
  expect_equal(names(h0), names(h))
})

test_that("offtarget_profile rejects an unknown primer and a width mismatch", {
  hits <- exact_hit_table(c(AS = p), g)
  expect_error(
    offtarget_profile(hits, c(OTHER = p), g),
    "must be a name in `primers`"
  )
  bad <- hits
  bad$end <- bad$end + 1L
  expect_error(
    offtarget_profile(bad, c(AS = p), g),
    "gapless equal-width"
  )
})

test_that("an inner miss is tolerated, a 3'-terminal one never is", {
  hits <- exact_hit_table(c(AS = p), g, max_mismatch = 1L)
  for (scheme in c("pcr", "paralog")) {
    h <- offtarget_profile(hits, c(AS = p), g, scheme = scheme)
    expect_true(h$tolerated[h$start == 93L])
    expect_false(h$tolerated[h$start == 52L])
  }
})

test_that("exact_hit_table with no hits keeps the full column contract", {
  none <- exact_hit_table(c(none = "CCCCCC"), g)
  expect_equal(nrow(none), 0L)
  expect_equal(names(none), hit_cols)
})
