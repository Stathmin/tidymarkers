# Virtual PCR with Biostrings::matchProbePair (no external tools needed).

test_that("in_silico_pcr finds the expected product and nothing else", {
  genome <- Biostrings::DNAStringSet(
    c(chr1 = paste0(strrep("A", 30L), "ACGTACGTACGTACGTACGT", strrep("T", 30L)))
  )
  pr <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$width, 20L)
})

test_that("in_silico_pcr reports no product when a primer is absent", {
  genome <- Biostrings::DNAStringSet(c(chr1 = strrep("A", 60L)))
  expect_equal(nrow(in_silico_pcr("ACGTACGTACGTACGTACGT", "TTTTGGGGCCCCAAAATTTT", genome)), 0L)
})

# A primer pair is a product on ONE chromosome.  The scan used to pool the hits
# of every chromosome and pair them inside a per-chromosome loop, so a forward
# site on chr1 and a reverse site on chr2 produced a phantom product - once per
# chromosome in the genome, each labelled with the loop's chromosome.
test_that("in_silico_pcr never pairs hits from two different chromosomes", {
  f <- "ACGTACGTACGTACGTACGT"
  r <- "TTTTGGGGCCCCAAAATTTT"
  r_site <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(r)))
  genome <- Biostrings::DNAStringSet(c(
    chr1 = paste0(strrep("A", 30L), f, strrep("T", 40L)),
    chr2 = paste0(strrep("C", 30L), r_site, strrep("G", 20L))
  ))
  # the two sites are on different chromosomes, so there is no amplicon
  expect_equal(nrow(in_silico_pcr(f, r, genome)), 0L)
  # both on one chromosome: exactly one product, on that chromosome
  together <- Biostrings::DNAStringSet(c(
    chr1 = paste0(strrep("A", 30L), f, strrep("T", 40L), r_site, strrep("G", 20L)),
    chr2 = paste0(strrep("C", 30L), r_site, strrep("G", 20L))
  ))
  pr <- in_silico_pcr(f, r, together)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$chrom, "chr1")
  expect_equal(pr$width, 80L)
})

test_that("in_silico_pcr is invariant to the parallel backend", {
  skip_if_not_installed("BiocParallel")
  amp <- "GGGGGGGGGGACGTACGTACGTACGTACGTGGGGGGGGGG"
  genome <- Biostrings::DNAStringSet(c(
    chr1 = paste0(amp, strrep("C", 40L), amp),
    chr2 = amp
  ))
  serial <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
  par <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_equal(par, serial)
  skip_on_os("windows")
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  mc <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome,
    BPPARAM = bp
  )
  expect_equal(mc, serial)
})

test_that("in_silico_pcr supports degenerate primers", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGGAAAATGCAAAACGT"))
  pr <- in_silico_pcr("ACGS", "TGCA", genome)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$start, 1L)
  expect_equal(pr$end, 12L)
  expect_equal(pr$width, 12L)
})

test_that("in_silico_pcr finds a product with one mismatch in a primer", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGGAAAATGCAAAACGT"))
  # ACGT vs ACGG is one mismatch at the 4th base; the reverse primer stays intact
  pr <- in_silico_pcr("ACGT", "TGCA", genome, max_mismatch = 1L)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$start, 1L)
  expect_equal(pr$end, 12L)
  # exact mode still finds nothing; the near-match warning is asserted below
  expect_equal(nrow(suppressWarnings(in_silico_pcr("ACGT", "TGCA", genome))), 0L)
})

test_that("exact mode warns when a one-mismatch product is hidden", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGGAAAATGCAAAACGT"))
  expect_warning(
    in_silico_pcr("ACGT", "TGCA", genome),
    "near-match"
  )
})

test_that("in_silico_pcr finds off-target products too", {
  amp <- "GGGGGGGGGGACGTACGTACGTACGTACGTGGGGGGGGGG"
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(amp, strrep("C", 40L), amp)))
  pr <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
  # two intra-copy products plus the real cross-copy product (forward in copy 1,
  # reverse in copy 2)
  expect_equal(nrow(pr), 3L)
})

test_that("the specificity engine advertises in-silico PCR without changing the matcher", {
  eng <- marker_engines()$specificity
  expect_true(capabilities(eng, "in_silico_pcr"))
  expect_true(capabilities(eng, "iupac"))
  expect_identical(capabilities(eng, "not_a_capability"), "unknown")
  # metadata only: the matcher keeps its signature (no conditions argument)
  expect_false("conditions" %in% names(formals(in_silico_pcr)))
  expect_identical(eng$pcr, in_silico_pcr)
})

# ---------------------------------------------------------------------------
# Position- and strand-aware off-target scoring (ROADMAP Z5).
#
# The amplicon is a short forward primer, a C/G spacer and the reverse primer's
# reverse complement; single bases are swapped in the genome, so every product
# below is a genuine near-match and not an artefact of the search radius.

z5_snp_f <- "ACGTACGTAC"
z5_snp_r <- "TTTTGGGGAAAACCCCGGGG"

z5_rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))

# swap one base of `x` for one it is not (so a mutation always changes a base)
z5_swap <- function(x, pos, to) {
  bases <- strsplit(x, "", fixed = TRUE)[[1L]]
  if (bases[pos] == to) to <- setdiff(c("A", "C", "G", "T"), to)[1L]
  bases[pos] <- to
  paste(bases, collapse = "")
}

# a list-column element must be a bare integer vector (an AsIs one would nest)
z5_positions <- function(col, i) {
  v <- col[[i]]
  if (is.list(v) && length(v) == 1L) v[[1L]] else v
}

z5_amplicon <- function(f = z5_snp_f, r_site = z5_rc(z5_snp_r)) {
  paste0("ACCAACCAAC", f, "CCAAAGGGTCCAAAAGGGTCC", r_site, "CCAACCAACC")
}

z5_snp_genome <- function(...) Biostrings::DNAStringSet(c(chr1 = z5_amplicon(...)))

# the best-scoring, then longest product: the real product, not an end-padded
# artefact of the search radius
z5_product <- function(prods) {
  keep <- prods[prods$n_mismatch == min(prods$n_mismatch), , drop = FALSE]
  keep[keep$width == max(keep$width), , drop = FALSE][1L, , drop = FALSE]
}

test_that("mismatch positions are strand-corrected into each primer's own frame", {
  # one base of the reverse primer's plus-strand site is swapped; the reverse
  # primer is searched as its reverse complement, so site position 1 is the
  # primer's 5' end (position 20) and site position 20 is its 3' end (position 1)
  cases <- tibble::tribble(
    ~site_pos, ~primer_pos, ~terminal,
    20L, 1L, TRUE,
    10L, 11L, FALSE,
    1L, 20L, FALSE
  )
  for (i in seq_len(nrow(cases))) {
    cs <- cases[i, ]
    r_site <- z5_swap(z5_rc(z5_snp_r), cs$site_pos, "T") # T is no site base here
    pr <- z5_product(suppressWarnings(
      in_silico_pcr(
        z5_snp_f, z5_snp_r, z5_snp_genome(r_site = r_site),
        max_mismatch = 1L
      )
    ))
    expect_equal(z5_positions(pr$r_mismatch_positions, 1L), cs$primer_pos)
    expect_length(z5_positions(pr$f_mismatch_positions, 1L), 0L)
    expect_equal(pr$three_prime_terminal, cs$terminal)
  }
  # the forward primer is searched as itself: no mirroring
  for (pos in c(1L, 5L, 10L)) {
    pr <- z5_product(suppressWarnings(
      in_silico_pcr(
        z5_snp_f, z5_snp_r, z5_snp_genome(f = z5_swap(z5_snp_f, pos, "T")),
        max_mismatch = 1L
      )
    ))
    expect_equal(z5_positions(pr$f_mismatch_positions, 1L), pos)
    expect_length(z5_positions(pr$r_mismatch_positions, 1L), 0L)
  }
})

test_that("a 3'-terminal mismatch dominates the score, on either primer", {
  # the forward primer, one mismatch at each end and in the middle
  cases <- tibble::tribble(
    ~pos, ~penalty, ~terminal,
    1L, 1, FALSE, # 5' end: a plain mismatch
    5L, 1, FALSE, # middle
    6L, 8, FALSE, # last 5 nt, still extendable
    8L, 12, TRUE, # first of the last 3 nt
    10L, 12, TRUE # 3' terminal
  )
  for (i in seq_len(nrow(cases))) {
    pr <- z5_product(suppressWarnings(
      in_silico_pcr(
        z5_snp_f, z5_snp_r, z5_snp_genome(f = z5_swap(z5_snp_f, cases$pos[i], "G")),
        max_mismatch = 1L
      )
    ))
    expect_equal(pr$penalty, cases$penalty[i])
    expect_equal(pr$three_prime_terminal, cases$terminal[i])
  }
  # the reverse primer's 3' end (primer position 20) gets the same treatment
  r_site <- z5_swap(z5_rc(z5_snp_r), 20L, "G")
  term <- z5_product(suppressWarnings(
    in_silico_pcr(z5_snp_f, z5_snp_r, z5_snp_genome(r_site = r_site), max_mismatch = 1L)
  ))
  expect_true(term$three_prime_terminal)
  expect_equal(term$r_penalty, 12L * 1)
  # an exact 3' end is never penalised for a mismatch elsewhere: the middle
  # mismatch costs its plain weight, not the terminal one
  mid <- z5_product(suppressWarnings(
    in_silico_pcr(
      z5_snp_f, z5_snp_r, z5_snp_genome(f = z5_swap(z5_snp_f, 5L, "G")),
      max_mismatch = 1L
    )
  ))
  expect_false(mid$three_prime_terminal)
  expect_lt(mid$penalty, 12)
  expect_gt(mid$score, term$score)
})

test_that("a mid-primer SNP 100 bp from the other primer is a reported near-match", {
  # the user-reported case: one SNP in the middle of the forward primer, the
  # reverse primer far away (an ordinary amplicon length)
  fill <- paste0(rep("CCAAAGGGTC", 10L), collapse = "") # 100 bp spacer
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
    "ACCAACCAAC", z5_swap(z5_snp_f, 5L, "G"), fill, z5_rc(z5_snp_r), "CCAACCAACC"
  )))
  # exact mode: no product, and the near-match is not silently hidden
  expect_warning(exact <- in_silico_pcr(z5_snp_f, z5_snp_r, genome), "near-match")
  expect_equal(nrow(exact), 0L)
  # reported: the full product, with its position and its penalty
  pr <- in_silico_pcr(z5_snp_f, z5_snp_r, genome, max_mismatch = 1L)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$width, 10L + 100L + 20L)
  expect_equal(pr$n_mismatch, 1L)
  expect_equal(z5_positions(pr$f_mismatch_positions, 1L), 5L)
  expect_length(z5_positions(pr$r_mismatch_positions, 1L), 0L)
  expect_false(pr$three_prime_terminal)
  expect_equal(pr$f_penalty, 1L * 1) # its own weight, not the terminal one
  expect_equal(pr$penalty, 1L * 1)
  expect_equal(pr$score, 90)
  # ... and classified by the chemistry: a middle mismatch is tolerated by pcr
  kasp <- in_silico_pcr(z5_snp_f, z5_snp_r, genome, scheme = "kasp")
  expect_equal(kasp$scheme, "pcr")
  expect_equal(nrow(kasp), 1L)
  expect_true(kasp$tolerated)
})

test_that("a degenerate match is not counted as a mismatch, with or without a radius", {
  # ACST against ACCT: S covers C, so this is an exact IUPAC match
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
    "ACCAACCAAC", "ACCT", "CCAAAGGGTCCAAAAGGGTCC", z5_rc(z5_snp_r), "CCAACCAACC"
  )))
  pr <- in_silico_pcr("ACST", z5_snp_r, genome)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$n_mismatch, 0L)
  expect_length(z5_positions(pr$mismatch_positions, 1L), 0L)
  expect_equal(pr$penalty, 0)
  expect_equal(pr$score, 100)
  # With a radius only the exact IUPAC scan is used (it says so): the covered
  # code is still a match, and the one position is counted once, never twice
  short <- Biostrings::DNAStringSet(c(chr1 = paste0("G", "ACCT", "G", z5_rc(z5_snp_r))))
  rad <- withCallingHandlers(
    in_silico_pcr("ACST", z5_snp_r, short, max_mismatch = 1L),
    tidymarkers_degenerate_exact = function(w) invokeRestart("muffleWarning")
  )
  expect_equal(nrow(rad), 1L)
  expect_equal(rad$n_mismatch, 0L)
  expect_equal(rad$penalty, 0)
})

test_that("the scheme's own policy weights the penalty", {
  genome <- z5_snp_genome(f = z5_swap(z5_snp_f, 5L, "G"))
  kasp <- z5_product(in_silico_pcr(z5_snp_f, z5_snp_r, genome, scheme = "kasp"))
  rpa <- z5_product(in_silico_pcr(z5_snp_f, z5_snp_r, genome, scheme = "rpa"))
  # both find the same one-mismatch product, and both tolerate it
  expect_equal(kasp$n_mismatch, rpa$n_mismatch)
  expect_true(kasp$tolerated && rpa$tolerated)
  expect_equal(kasp$penalty, 1L * 1)
  # a 10-nt primer's last-5 window starts at position 6: rpa weighs it lower
  expect_equal(specificity_policy("rpa", length = 10L)$weight[6L], 5)
  expect_equal(specificity_policy("pcr", length = 10L)$weight[6L], 8)
  # the 3'-terminal rule is the same for every chemistry: no scheme rescues it
  term <- z5_snp_genome(f = z5_swap(z5_snp_f, 10L, "G"))
  for (scheme in c("kasp", "rpa")) {
    pr <- z5_product(in_silico_pcr(z5_snp_f, z5_snp_r, term, scheme = scheme))
    expect_true(pr$three_prime_terminal)
    expect_false(pr$tolerated)
    expect_equal(pr$penalty, 12L * 1)
  }
})

test_that("many exact hits pair without an outer product", {
  # A low-complexity primer has thousands of exact sites.  Pairing every forward
  # start with every reverse end via `outer()` is a dense 1e11-cell matrix on a
  # real genome (99.6 GB, allocation failure), although only a reverse end at
  # least one primer width downstream of a start can be a product.  So a dense
  # synthetic sequence must still answer, and the answer must obey that rule.
  f <- "ACGTTGCAAGTCCTAGGTAC"
  r <- "GCTAGGCATTGACCATGCTA"
  # in_silico_pcr searches `f` and the reverse complement of `r` on the plus strand
  rev <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(r)))
  spacer <- paste0(strrep("T", 6L), "GGA")
  unit <- paste0(f, spacer, rev)
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(rep(unit, 40L), collapse = "")))
  pr <- in_silico_pcr(f, r, genome, max_mismatch = 0L)
  expect_gt(nrow(pr), 0L)
  expect_equal(nrow(pr), nrow(dplyr::distinct(pr, .data$start, .data$end)))
  # every product: the forward start is an `f` site and the reverse end an
  # `rev` site, at least one primer width apart
  expect_true(all(pr$end - pr$start + 1L >= nchar(f)))
  chr <- as.character(genome[["chr1"]])
  expect_true(all(unname(substr(chr, pr$start, pr$start + nchar(f) - 1L)) == f))
  expect_true(all(unname(substr(chr, pr$end - nchar(r) + 1L, pr$end)) == rev))
})
