# Specificity seed counting: a known seed must be counted on both strands.

test_that("specificity_seed counts seed occurrences on both strands", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  primers <- c(P = "TTTTTTAAACCCGGG") # 3'-terminal 9 nt = AAACCCGGG
  # forward "AAACCCGGG" x2 + reverse-complement "CCCGGGTTT" x2 = 4
  res <- specificity_seed(primers, genome, k = 9L, max_mismatch = 0L)
  expect_equal(res$hits, 4L)
  expect_equal(res$seed, "AAACCCGGG")
})

test_that("specificity_seed returns zero hits for an absent seed", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  primers <- c(P = "TTTTTTACGACGAC")
  expect_equal(specificity_seed(primers, genome, k = 9L, max_mismatch = 0L)$hits, 0L)
})

test_that("the mismatch scan is invariant to the parallel backend", {
  skip_if_not_installed("BiocParallel")
  genome <- Biostrings::DNAStringSet(c(
    chr1 = "AAACCCGGGTTTAAACCCGGGTTT",
    chr2 = "GGGTTTAAACCC"
  ))
  primers <- c(P = "TTTTTTAAACCCGGG", Q = "TTTTTTGGGAAACCC")
  serial <- specificity_seed(primers, genome, k = 9L, max_mismatch = 1L)
  par <- specificity_seed(primers, genome,
    k = 9L, max_mismatch = 1L,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_equal(par$hits, serial$hits)
  skip_on_os("windows")
  skip_if_not_full()
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  mc <- specificity_seed(primers, genome,
    k = 9L, max_mismatch = 1L,
    BPPARAM = bp
  )
  expect_equal(mc$hits, serial$hits)
})

test_that("specificity counts the reference genome, not the variant/allele list", {
  # the seed is looked up in the reference genome only, so a position that merely
  # *varies* across accessions (a row in the variant table) is one reference locus,
  # not a phantom repeat; a genuine duplication in the genome is what adds a hit.
  seed <- "GTCAGTCC"
  primer <- paste0("AAAA", seed)
  one <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTGTCAGTCCACGTACGT"))
  two <- Biostrings::DNAStringSet(c(chr1 = "GTCAGTCCACGTACGTGTCAGTCC"))
  expect_equal(specificity_seed(c(P = primer), one, k = 8L, max_mismatch = 0L)$hits, 1L)
  expect_equal(specificity_seed(c(P = primer), two, k = 8L, max_mismatch = 0L)$hits, 2L)
  # the specificity API takes a genome, never a variant table - so it cannot
  # double-count the alternative alleles of a multi-variant feature
  expect_false("variants" %in% names(formals(specificity_seed)))
})

test_that("the exact PDict path matches a manual both-strand count", {
  genome <- Biostrings::DNAStringSet(c(
    chr1 = "AAACCCGGGTTTAAACCCGGGTTT",
    chr2 = "GGGAAACCC"
  ))
  primers <- c(P = "TTTTTTAAACCCGGG", Q = "TTTTTTGGGAAACCC")
  res <- specificity_seed(primers, genome, k = 9L, max_mismatch = 0L)
  manual <- vapply(unname(substr(primers, 8L, 16L)), function(s) {
    d <- Biostrings::DNAString(s)
    r <- Biostrings::reverseComplement(d)
    sum(Biostrings::vcountPattern(d, genome)) + sum(Biostrings::vcountPattern(r, genome))
  }, numeric(1))
  expect_equal(res$hits, as.integer(manual))
})

test_that("a BiocParallel BPPARAM gives the same mismatch counts as serial", {
  skip_if_not_installed("BiocParallel")
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  primers <- c(P = "TTTTTTAAACCCGGG")
  serial <- specificity_seed(primers, genome, k = 9L, max_mismatch = 1L)
  par <- specificity_seed(primers, genome,
    k = 9L, max_mismatch = 1L,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_equal(par$hits, serial$hits)
})

test_that("the engines registry reports unknown for an undeclared capability", {
  eng <- marker_engines()
  expect_named(eng, c("thermo", "specificity"))
  expect_identical(capabilities(eng, "not_a_capability"), "unknown")
  expect_identical(capabilities(eng$thermo, "not_a_capability"), "unknown")
  expect_identical(capabilities(eng$specificity, "not_a_capability"), "unknown")
  expect_true(capabilities(eng, "seed_counts"))
  # a non-engine object is unknown too, never a silent pass
  expect_identical(capabilities(list(), "seed_counts"), "unknown")
})

# ---------------------------------------------------------------------------
# Mismatch-tolerant, position-aware specificity (ROADMAP Z5).
#
# A synthetic genome, so none of this needs an external aligner:
#   5' AAAAAAAAAA <forward> GGGGGGGGGGGGGGGGGGGG <rc(reverse)> TTTTTTTTTT 3'
# The product is the forward primer + 20 G/C nt + the reverse primer.

z5_f <- "ACGTACGTACGTACGTACGT"
z5_r <- "TTTTGGGGAAAACCCCGGGG"

z5_rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))

z5_genome <- function(f = z5_f, r = z5_r, fill = "GGGGGGGGGGGGGGGGGGGG") {
  Biostrings::DNAStringSet(c(chr1 = paste0("AAAAAAAAAA", f, fill, z5_rc(r), "TTTTTTTTTT")))
}

# the same, but spacers and flanks avoid long single-base runs, which a
# mismatch-tolerant scan otherwise pads the product ends with
z5_genome_hard <- function(f = z5_f) {
  Biostrings::DNAStringSet(c(chr1 = paste0(
    "ACCAACCAAC", f, "CCAAAGGGTCCAAAAGGGTCC", z5_rc(z5_r), "CCAACCAACC"
  )))
}

# a tandem duplication of the whole amplicon, no single-base runs: exactly two
# products from an exact search
a <- paste0("ACCAACCAAC", z5_f, "CCAAAGGGTCCAAAAGGGTCC", z5_rc(z5_r), "CCAACCAACC")
z5_genome_twice <- function() {
  Biostrings::DNAStringSet(c(chr1 = a, chr2 = a))
}

# swap a base at `pos` (1-based from the primer's 5' end) for a different one
z5_mutate <- function(x, pos) {
  bases <- strsplit(x, "", fixed = TRUE)[[1L]]
  bases[pos] <- if (bases[pos] == "A") "C" else "A"
  paste(bases, collapse = "")
}

test_that("a perfect product reports zero mismatches and the legacy columns", {
  prods <- in_silico_pcr(z5_f, z5_r, z5_genome())
  expect_equal(nrow(prods), 1L)
  expect_true(all(c("chrom", "start", "end", "width") %in% names(prods)))
  expect_true(all(PCR_MISMATCH_COLUMNS %in% names(prods)))
  expect_equal(prods$n_mismatch, 0L)
  expect_length(prods$mismatch_positions[[1L]], 0L)
  expect_equal(prods$three_prime_mismatch, 0L)
  expect_equal(prods$seed_mismatch, 0L)
  expect_false(prods$three_prime_terminal)
  expect_equal(prods$penalty, 0)
  expect_equal(prods$score, 100)
  expect_equal(prods$start, 11L)
  expect_equal(prods$width, 20L + 20L + 20L)
})

test_that("a middle primer mismatch with a far second primer is detected and reported", {
  # the blind spot: exact mode finds nothing, and a 3'-terminal seed cannot see a
  # substitution in the middle of the forward primer
  genome <- z5_genome(f = z5_mutate(z5_f, 10L))
  expect_warning(exact <- in_silico_pcr(z5_f, z5_r, genome), "near-match")
  expect_equal(nrow(exact), 0L)
  # one middle mismatch: found with a mismatch-tolerant scan and fully reported
  prods <- in_silico_pcr(z5_f, z5_r, genome, max_mismatch = 1L)
  expect_equal(nrow(prods), 1L)
  expect_equal(prods$n_mismatch, 1L)
  expect_equal(prods$mismatch_positions[[1L]], 10L)
  expect_equal(prods$three_prime_mismatch, 0L)
  expect_equal(prods$seed_mismatch, 0L)
  expect_false(prods$three_prime_terminal)
  # the far second primer still pairs: the product is the full amplicon
  expect_equal(prods$width, 60L)
})

test_that("a 3'-end mismatch is weighted harder than the same middle mismatch", {
  mid <- in_silico_pcr(z5_f, z5_r, z5_genome(f = z5_mutate(z5_f, 10L)), max_mismatch = 1L)
  term <- in_silico_pcr(z5_f, z5_r, z5_genome(f = z5_mutate(z5_f, 20L)), max_mismatch = 1L)
  expect_equal(nrow(mid), 1L)
  expect_equal(nrow(term), 1L)
  expect_equal(mid$mismatch_positions[[1L]], 10L)
  expect_equal(term$mismatch_positions[[1L]], 20L)
  # a terminal substitution counts as three-prime, and as seed only for pcr
  expect_equal(term$three_prime_mismatch, 1L)
  expect_true(term$three_prime_terminal)
  expect_gt(term$penalty, mid$penalty)
  expect_lt(term$score, mid$score)
  # the per-primer shares: only the forward primer changed
  expect_equal(mid$r_penalty, 0)
  expect_equal(term$r_penalty, 0)
  expect_gt(term$f_penalty, mid$f_penalty)
  # a mismatch four nt from the end is in the last-5 window but still extendable
  near <- in_silico_pcr(z5_f, z5_r, z5_genome(f = z5_mutate(z5_f, 16L)), max_mismatch = 1L)
  expect_equal(near$three_prime_mismatch, 1L)
  expect_false(near$three_prime_terminal)
  expect_gt(near$penalty, mid$penalty)
  expect_lt(near$penalty, term$penalty)
})

test_that("an IUPAC degenerate primer is a match, not a mismatch", {
  # ACST: the subject's 3rd base is C, which S covers - a stale exact-hit
  # PRIMER_MISMATCH would call it a mismatch; the trailing T is a real base
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
    "AAAAACGTGGGGGGGGGGGGGGGGGG", z5_rc(z5_r)
  )))
  prods <- in_silico_pcr("ACST", z5_r, genome)
  expect_equal(nrow(prods), 1L)
  expect_equal(prods$n_mismatch, 0L)
  expect_length(prods$mismatch_positions[[1L]], 0L)
  expect_equal(prods$penalty, 0)
  expect_equal(prods$score, 100)
})

test_that("the pcr policy tolerates a middle mismatch but not a 3'-end one", {
  mid <- in_silico_pcr(z5_f, z5_r, z5_genome(f = z5_mutate(z5_f, 10L)), scheme = "kasp")
  expect_true(all(c("scheme", "tolerated") %in% names(mid)))
  expect_equal(mid$scheme, "pcr")
  expect_equal(nrow(mid), 1L)
  expect_true(mid$tolerated)
  term <- in_silico_pcr(z5_f, z5_r, z5_genome(f = z5_mutate(z5_f, 20L)), scheme = "kasp")
  expect_false(any(term$tolerated))
})

test_that("the rpa policy accepts a tolerated mismatch that pcr rejects", {
  # more mismatches than pcr's radius: tolerated by rpa (outside the 3' end),
  # but invisible to the exact pcr scan - so exact-hit-only is wrong for rpa
  genome <- z5_genome(f = z5_mutate(z5_mutate(z5_f, 4L), 10L))
  expect_equal(nrow(in_silico_pcr(z5_f, z5_r, genome, scheme = "pcr")), 0L)
  rpa_products <- in_silico_pcr(z5_f, z5_r, genome, scheme = "rpa")
  expect_gt(nrow(rpa_products), 1L)
  tolerated <- rpa_products[rpa_products$tolerated, , drop = FALSE]
  expect_equal(nrow(tolerated), 1L)
  expect_equal(tolerated$n_mismatch, 2L)
  expect_false(tolerated$three_prime_terminal)
  # the same product, graded: one tolerated product is a pass
  v <- specificity_verdict(rpa_products, "rpa", genome = genome)
  expect_equal(v$verdict, "pass")
  expect_equal(v$n_tolerated, 1L)
  expect_equal(v$n_mismatch, 2L)
})

test_that("crispr grades a near-cognate off-target, not just exact hits", {
  # a second site differing from the guide by a 5' block, so the exact scan cannot
  # see it; the 10-nt spacer left of it and the poly-C filler right of it are
  # themselves no off-target
  cognate <- paste0("CCAAAGGGTCCAAAAGGGTCC", z5_rc(z5_r), "CCAACCAACC")
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
    "ACCAACCAACACCAACCAAC", "AAAAAAAAAAA", "CAAAAAAAGGGGGGGGGGG", cognate,
    "ACCAACCAAC", z5_f, cognate
  )))
  expect_equal(nrow(in_silico_pcr(z5_f, z5_r, genome)), 1L) # one exact hit only
  near <- in_silico_pcr(z5_f, z5_r, genome, scheme = "cas12a")
  # the near-cognate is a full-length product that the exact scan misses
  off <- near[!near$tolerated & near$width > 63L, , drop = FALSE]
  expect_gt(nrow(off), 0L)
  expect_true(all(off$n_mismatch >= 4L))
  v <- specificity_verdict(near, "cas12a", genome = genome, pam = "TTTA")
  # the near-cognate alone does not break the policy: it is not counted as a
  # tolerated product, so the guide verdict still rests on the one exact hit
  expect_equal(v$n_products, nrow(near))
  expect_gt(v$n_products, 1L)
  expect_equal(v$best$width, 61L) # the complete product, not a padded end
})

test_that("a CRISPR PAM change is a hard loss of target, not a tolerated mismatch", {
  genome <- z5_genome_hard()
  prods <- in_silico_pcr(z5_f, z5_r, genome)
  expect_equal(nrow(prods), 1L)
  expect_equal(prods$n_mismatch, 0L)
  # the crispr scan radius also reports products with a padded end; those are
  # artefacts of the radius, so the best tolerated evidence is the tight product
  scanned <- in_silico_pcr(z5_f, z5_r, genome, scheme = "cas12a")
  intact <- specificity_verdict(scanned, "cas12a", genome = genome, pam = "TTTA")
  expect_equal(intact$pam, "match")
  # the complete product is the evidence, not an end-padded scan artefact
  expect_equal(max(intact$products$width[intact$products$tolerated & intact$products$n_mismatch == 0L]), 61L)
  expect_equal(intact$best$width, 61L)
  expect_equal(intact$best$n_mismatch, 0L)
  # a changed PAM: the guide has lost its target, whatever the mismatches say
  # (TTTG is still a valid TTTV PAM; AATA is not)
  damaged <- specificity_verdict(scanned, "cas12a", genome = genome, pam = "AATA")
  expect_equal(damaged$verdict, "fail")
  expect_equal(damaged$pam, "pam_loss")
  expect_equal(damaged$n_tolerated, 0L)
  expect_true(all(vapply(damaged$products$violations, function(v) "pam" %in% v, logical(1))))
  expect_equal(pam_state("TTTA"), "match")
  expect_equal(pam_state("TTTC"), "match") # V = any of A/C/G
  expect_equal(pam_state("TTTG"), "match") # still TTTV
  expect_equal(pam_state("AATA"), "pam_loss")
  expect_equal(pam_state(NULL), "absent")
  # Cas13 has no DNA PAM but its PFS is a weak preference, never "PAM-free"
  expect_equal(pam_state("PFS"), "match")
  expect_equal(specificity_policy("cas13")$pfs, "weak")
})

test_that("a missing genome or engine gives unknown, never a silent pass", {
  no_genome <- specificity_verdict(empty_pcr_products(), "kasp", genome = NULL)
  expect_equal(no_genome$verdict, "unknown")
  expect_true(nzchar(no_genome$reason))
  expect_true(is.na(no_genome$backend))
  # a genome with no in-memory sequences and no BLAST database is not an engine
  no_engine <- specificity_verdict(empty_pcr_products(), "kasp", genome = "not_a_genome")
  expect_equal(no_engine$verdict, "unknown")
  expect_false(no_engine$backend$engine)
  expect_false(no_engine$backend$capabilities[["exact_scan"]])
  # CRISPR genome-wide off-target search is not implemented: the flag says so
  expect_false(specificity_backend(z5_genome(), "crispr")$capabilities[["crispr_offtarget"]])
  expect_true(specificity_backend(z5_genome(), "crispr")$engine)
})

test_that("a verdict fails when the products are not specific", {
  # the same primer pair on two chromosomes: two products, so neither is unique
  twice <- z5_genome_twice()
  prods <- in_silico_pcr(z5_f, z5_r, twice)
  expect_equal(nrow(prods), 2L)
  v <- specificity_verdict(prods, "kasp", genome = twice)
  expect_equal(v$verdict, "fail")
  expect_equal(v$n_tolerated, 2L)
  # no product at all is a failure of the assay, not a pass
  none <- specificity_verdict(empty_pcr_products(), "kasp", genome = z5_genome())
  expect_equal(none$verdict, "fail")
  expect_equal(none$n_products, 0L)
  one <- specificity_verdict(
    in_silico_pcr(z5_f, z5_r, z5_genome_hard()), "kasp",
    genome = z5_genome_hard()
  )
  expect_equal(one$verdict, "pass")
})

test_that("specificity_policy encodes the chemistry tolerances", {
  pcr <- specificity_policy("pcr")
  rpa <- specificity_policy("rpa")
  crispr <- specificity_policy("crispr")
  paralog <- specificity_policy("paralog")
  expect_false(pcr$calibrated)
  expect_equal(pcr$mode, "pcr")
  expect_equal(pcr$max_mismatch, 1L)
  expect_equal(pcr$three_prime_max, 0L)
  expect_false(pcr$requires_pam)
  # rpa is the tolerant chemistry: several mismatches can still amplify
  expect_equal(rpa$mode, "rpa")
  expect_gte(rpa$max_mismatch, 3L)
  expect_gt(rpa$max_mismatch, pcr$max_mismatch)
  expect_gt(rpa$three_prime_max, pcr$three_prime_max)
  # crispr: a PAM-proximal seed (~8-12 nt) and a licensing PAM
  expect_true(crispr$requires_pam)
  expect_equal(crispr$pam, "TTTV")
  expect_gte(crispr$seed_len, 8L)
  expect_lte(crispr$seed_len, 12L)
  expect_equal(crispr$seed_max, 1L)
  expect_equal(paralog$mode, "paralog")
  expect_gt(paralog$max_mismatch, 0L)
  expect_false(paralog$requires_pam)
  # only crispr requires a PAM
  expect_identical(
    vapply(list(pcr, rpa, crispr, paralog), `[[`, logical(1), "requires_pam"),
    c(FALSE, FALSE, TRUE, FALSE)
  )
  # aliases resolve to the same policy as the canonical name
  expect_equal(specificity_policy("KASP")$mode, "pcr")
  expect_equal(specificity_policy("lamp")$mode, "rpa")
  expect_equal(specificity_policy("homeolog")$mode, "paralog")
  expect_equal(specificity_policy("Cas12a")$mode, "crispr")
  expect_error(specificity_policy("not_a_scheme"), "Unknown specificity scheme")
  # the per-position ladder: terminal > last-5 > middle
  w <- specificity_policy("pcr", length = 20L)$weight
  expect_length(w, 20L)
  expect_gt(w[20L], w[16L])
  expect_gt(w[16L], w[10L])
})

test_that("the policy table carries the per-scheme position ladder", {
  fields <- c("weight_middle", "weight_seed", "weight_three_prime", "weight_terminal")
  pcr <- specificity_policy("pcr")
  rpa <- specificity_policy("rpa")
  crispr <- specificity_policy("crispr")
  # pcr is the historical scoring: its weights must not drift
  expect_equal(unname(unlist(pcr[fields])), c(1, 3, 8, 12))
  expect_equal(unname(unlist(specificity_policy("kasp")[fields])), c(1, 3, 8, 12))
  for (pol in list(pcr, rpa, crispr)) {
    expect_true(all(fields %in% names(pol)))
    expect_true(all(unlist(pol[fields]) > 0))
    # a terminal mismatch is never the cheap position, whatever the chemistry
    expect_gt(pol$weight_terminal, pol$weight_three_prime)
    expect_gt(pol$weight_three_prime, pol$weight_middle)
  }
  # RPA tolerates mismatches: the middle and the seed cost less than for PCR,
  # and its seed weighs less than its last-5 window
  expect_lt(rpa$weight_seed, pcr$weight_seed)
  expect_lt(rpa$weight_three_prime, pcr$weight_three_prime)
  expect_lt(rpa$weight_seed, rpa$weight_three_prime)
  # the same terminal rule everywhere (no chemistry extends a mismatched 3' end)
  expect_equal(rpa$weight_terminal, pcr$weight_terminal)
  # crispr weighs its PAM-proximal seed heavier than pcr does
  expect_gt(crispr$weight_seed, pcr$weight_seed)
  # the ladder vector agrees with the fields
  w <- specificity_policy("rpa", length = 20L)$weight
  expect_equal(w[5L], rpa$weight_middle)
  expect_equal(w[15L], rpa$weight_middle)
  expect_equal(w[16L], rpa$weight_three_prime)
  expect_equal(w[20L], rpa$weight_terminal)
})

test_that("a degenerate primer is scored exactly, never rescuing a 3'-terminal mismatch", {
  # S covers C: an exact IUPAC match, and the reverse primer's site runs 3'->5'
  # on the plus strand, so the site's LAST base is the primer's 3'-terminal one
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
    "ACCAACCAAC", "ACCTACGTAC", "CCAAAGGGTCCAAAAGGGTCC", z5_rc(z5_r), "CCAACCAACC"
  )))
  pr <- in_silico_pcr("ACSTACGTAC", z5_r, genome)
  expect_equal(nrow(pr), 1L)
  expect_equal(pr$n_mismatch, 0L)
  expect_equal(pr$penalty, 0)
  # the reverse primer's site, with its 3'-terminal base swapped
  r_site <- z5_rc(z5_r)
  bases <- strsplit(r_site, "", fixed = TRUE)[[1L]]
  bases[length(bases)] <- "G"
  g <- Biostrings::DNAStringSet(c(chr1 = paste0(
    "ACCAACCAAC", "ACCTACGTAC", "CCAAAGGGTCCAAAAGGGTCC", paste(bases, collapse = ""),
    "CCAACCAACC"
  )))
  # the degenerate primer is scanned exactly, so the hit is the exact degenerate
  # site (the forward primer is a match, not a mismatch) and the product is the
  # one carrying the reverse primer's 3'-terminal substitution
  term <- suppressWarnings(in_silico_pcr("ACSTACGTAC", z5_r, g, max_mismatch = 1L))
  term <- term[term$width == 51L, , drop = FALSE]
  expect_equal(nrow(term), 1L)
  expect_length(unlist(term$f_mismatch_positions), 0L)
  expect_true(term$three_prime_terminal)
  expect_true(1L %in% unlist(term$r_mismatch_positions))
  expect_equal(term$r_penalty, 12L * 1)
})
