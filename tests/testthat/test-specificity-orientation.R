# ---------------------------------------------------------------------------
# Position-aware off-target in the primer's own orientation (ROADMAP Z5.1).
#
# The two sites a verdict must tell apart:
#   (i)  an EXACT off-target inside the amplicon - a second, perfect copy of the
#        whole product, i.e. the most dangerous off-target there is;
#   (ii) a genome-wide 3'-TERMINAL miss - a site whose *last* base differs from
#        the primer's 3' end, in the primer's own 5'->3' orientation, i.e. a
#        site the polymerase cannot extend.
# Both come out of the same genome-wide scan (`in_silico_pcr()` ->
# `primer_hit_table()` -> `mismatch_profile()`), so the distinction under test
# is a property of the scan, not of a separately hand-built product.

or_f <- "ACGTACGTACGTACGTACGT"
or_r <- "TTTTGGGGAAAACCCCGGGG"

or_rc <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
}

# swap one base at `pos` (1-based from the primer's 5' end)
or_mutate <- function(x, pos) {
  bases <- strsplit(x, "", fixed = TRUE)[[1L]]
  bases[pos] <- if (bases[pos] == "A") "C" else "A"
  paste(bases, collapse = "")
}

# A site is flanked by non-homopolymer fillers long enough to keep a
# mismatch-tolerant scan from inventing extra products, and by a spacer
# (`GGT`) shorter than a primer so no product can span two sites:
#   5' <right> <forward> <left> <rc(reverse)> 3'
# `left` and `right` are filler only, so the primer hit is the only evidence.
or_leads <- "GTACCTGATCAGTCAGTCA"
or_fills <- "TGTACGCTGAAGTACGTCA"

or_site <- function(f) {
  paste0(or_leads, f, or_fills, or_rc(or_r))
}

or_copies <- function(n) {
  paste(rep(or_site(or_f), n), collapse = "")
}

# One contig, three sites: an exact copy, then a second exact copy (the
# in-amplicon off-target), then one whose forward primer is terminally
# substituted at its position 20 (its own 3' end).
or_genome <- function() {
  Biostrings::DNAStringSet(c(chr1 = paste0(
    or_site(or_f), or_site(or_f), or_site(or_mutate(or_f, 20L))
  )))
}

# the same with the substitution at position 10 instead (a middle miss)
or_genome_mid <- function() {
  Biostrings::DNAStringSet(c(chr1 = paste0(
    or_site(or_f), or_site(or_f), or_site(or_mutate(or_f, 10L))
  )))
}

test_that("the genome-wide scan grades every near hit, not only degenerates", {
  genome <- or_genome()
  # the 3'-terminal site is one substitution away: a radius scan must find it,
  # and exact mode does not (it now says so, naming the near match)
  exact <- in_silico_pcr(or_f, or_r, genome, max_mismatch = 0L)
  # two exact sites, each hit by the one forward and the one reverse primer, so
  # `in_silico_pcr()` pairs every forward start with every downstream reverse
  # end: a product over one site (59 bp), over two (59 + 19 + 59 = 137 bp) and
  # over all three (215 bp)
  expect_setequal(unique(exact$width), c(59L, 137L, 215L))
  expect_true(all(exact$n_mismatch == 0L))
  expect_warning(
    in_silico_pcr(or_f, or_r, genome, max_mismatch = 0L),
    "near-match"
  )

  radius <- in_silico_pcr(or_f, or_r, genome, max_mismatch = 1L)
  expect_gt(nrow(radius[radius$n_mismatch > 0L, , drop = FALSE]), 0L)
  # every scanned hit carries the position ladder, not a bare count
  expect_true(all(PCR_MISMATCH_COLUMNS %in% names(radius)))
  expect_true(all(radius$penalty >= radius$n_mismatch))
  # the score is `100 - 10 * penalty` (hit_score()): it is *not* floored at 0,
  # so a 12-weight terminal miss legitimately scores below zero
  expect_equal(radius$score, round(100 - 10 * radius$penalty, 4L))
})

test_that("a 3'-terminal miss is scored as terminal in the primer's own frame", {
  genome <- or_genome()
  prods <- in_silico_pcr(or_f, or_r, genome, max_mismatch = 1L)
  term <- prods[prods$three_prime_terminal, , drop = FALSE]
  expect_equal(nrow(term), 1L)
  expect_equal(term$n_mismatch, 1L)
  # position 20 of the forward primer *as written 5'->3'*, not a plus-strand
  # offset: the profile is mirrored into the primer's own orientation
  expect_equal(term$f_mismatch_positions[[1L]], 20L)
  expect_length(term$r_mismatch_positions[[1L]], 0L)
  expect_equal(term$three_prime_mismatch, 1L)
  expect_gt(term$f_penalty, 0)
  expect_equal(term$r_penalty, 0)
  expect_gt(term$seed_mismatch, 0L) # the last-15 window holds it too

  v <- specificity_verdict(prods, "pcr", genome = genome)
  expect_true(
    "terminal" %in% v$products$violations[[which(v$products$three_prime_terminal)]]
  )
  expect_false(any(v$products$tolerated[v$products$three_prime_terminal]))
})

test_that("an exact in-amplicon off-target and a 3'-terminal miss get different verdicts", {
  genome <- or_genome()
  prods <- in_silico_pcr(or_f, or_r, genome, max_mismatch = 1L)

  # the exact duplicate is a tolerated product: more than one perfect product,
  # so the pair is not specific.  (Counted with `scheme`, not from a hand-built
  # row: the filter below is the scan's own output.)
  perfect <- prods[prods$n_mismatch == 0L, , drop = FALSE]
  dup <- specificity_verdict(
    in_silico_pcr(or_f, or_r, genome, scheme = "pcr"), "pcr",
    genome = genome
  )
  expect_gte(dup$n_products, nrow(perfect))
  expect_gte(dup$n_tolerated, 2L)
  expect_identical(dup$verdict, "fail")
  expect_equal(dup$best$penalty, 0)

  # (ii) keep only the first exact product and the 3'-terminal site: one
  # tolerated product is a pass - the opposite call on the same genome.  The
  # two rows are re-graded by `specificity_verdict()`, which recomputes
  # `tolerated` and `violations` from the scan columns.
  first <- perfect[perfect$start == min(perfect$start), , drop = FALSE][1L, ]
  term <- prods[prods$three_prime_terminal, , drop = FALSE]
  mixed <- rbind(
    as.data.frame(first)[names(empty_pcr_products())],
    as.data.frame(term)[names(empty_pcr_products())]
  )
  mixed_v <- specificity_verdict(mixed, "pcr", genome = genome)
  expect_equal(mixed_v$n_tolerated, 1L)
  expect_identical(mixed_v$verdict, "pass")
  expect_true("terminal" %in% mixed_v$products$violations[[2L]])
  expect_length(mixed_v$products$violations[[1L]], 0L)
})

test_that("an exact in-amplicon hit costs less than a genome-wide 3'-end miss", {
  genome <- or_genome()
  prods <- in_silico_pcr(or_f, or_r, genome, max_mismatch = 1L)
  perfect <- prods[prods$n_mismatch == 0L, , drop = FALSE]
  term <- prods[prods$three_prime_terminal, , drop = FALSE]
  pol <- specificity_policy("pcr")
  expect_equal(unique(perfect$penalty), 0)
  expect_equal(unique(perfect$score), 100)
  expect_gt(min(term$penalty), max(perfect$penalty))
  expect_equal(unique(term$penalty), pol$weight_terminal)
  expect_lt(min(term$score), min(perfect$score))

  # same edit distance, different position: the middle miss is cheaper, so the
  # ordering is positional and not a mismatch count in disguise
  mid <- in_silico_pcr(or_f, or_r, or_genome_mid(), max_mismatch = 1L)
  mid_off <- mid[mid$n_mismatch == 1L, , drop = FALSE]
  expect_equal(nrow(mid_off), 1L)
  expect_false(mid_off$three_prime_terminal)
  expect_equal(mid_off$f_mismatch_positions[[1L]], 10L)
  expect_equal(mid_off$n_mismatch, term$n_mismatch)
  expect_lt(mid_off$penalty, term$penalty)
  expect_gt(mid_off$score, term$score)
})
