# =============================================================================
#  marker_specificity_pcr.R  --  in-silico PCR and product pairing
# =============================================================================
#  Split out of marker_specificity.R (T7, cycle2/tidy-20261126).  The virtual
#  PCR half of the specificity engine: hit tables in, paired products out, the
#  chemistry filter applied, plus the pooled (near) pass and the assay-level
#  amplicon helpers.  Function names are unchanged by the split.
# =============================================================================

#' In-silico PCR: every product a primer pair would amplify
#'
#' Uses IUPAC-aware `Biostrings::matchPattern()` to find forward and reverse
#' primer matches, then pairs every forward start with a downstream reverse end
#' into a product.  This supports degenerate primers and `max_mismatch`; a KASP
#' test must amplify exactly one product per allele.
#'
#' Every product also carries its primer mismatches (`n_mismatch` = the worse of
#' the two primers), where they sit, how close they are to the extendable end,
#' and a weighted penalty.  Positions are 1-based from each primer's **own**
#' 5' end: `f_mismatch_positions` for the forward primer, `r_mismatch_positions`
#' for the reverse one (the reverse primer is searched as its reverse complement,
#' so its plus-strand positions are mirrored back - position 1 there is the
#' reverse primer's 3'-terminal base).  `mismatch_positions` is the union of the
#' two and is therefore only meaningful as a set, not by position.
#'
#' A mismatch in the last `THREE_PRIME_TERMINAL` (3) nt of *either* primer sets
#' `three_prime_terminal` and carries the policy's `weight_terminal`, so it always
#' dominates a mismatch elsewhere and no scheme tolerates it.  A degenerate code
#' that still covers the genomic base is a match, counted once, so neither
#' degeneracy nor a scheme radius can double-count it.  `three_prime_mismatch`
#' (last `THREE_PRIME_WINDOW` = 5 nt) and `seed_mismatch` (last `seed_len` nt)
#' are the sums over both primers.
#'
#' The penalty is the sum over the two primers of the scheme's position ladder
#' (`weight_middle`, `weight_seed`, `weight_three_prime`, `weight_terminal` of
#' [specificity_policy()]), so an RPA product with several middle mismatches
#' scores differently from a PCR one.  A degenerate primer can only be scanned
#' exactly (it is reported as exact matches, with a warning); a one-base radius
#' over an ambiguous code is not a well-defined distance.
#'
#' Passing a `scheme` makes the search *policy-aware*: the scan radius is the
#' policy's `max_mismatch` and a product is kept only when the policy tolerates
#' it (`max` total, seed and three-prime caps; a changed PAM is a hard loss for
#' CRISPR).  Without `scheme` the call behaves exactly as before: a search
#' radius of `max_mismatch`, no chemistry filter.
#'
#' @param f_primer forward primer (5'->3'), binds the plus strand.
#' @param r_primer reverse primer (5'->3').
#' @param genome DNAStringSet.
#' @param max_mismatch mismatches tolerated in the primer match (the search
#'   radius; with `scheme`, defaults to the policy's `max_mismatch`).
#' @param BPPARAM optional `BiocParallel` backend; the virtual PCR is split by
#'   genome sequence (one worker per chromosome with `MulticoreParam()`).
#' @param scheme optional specificity scheme (`"pcr"`, `"rpa"`, `"crispr"`,
#'   `"paralog"`); see [specificity_policy()].  When given, the returned table
#'   gains `scheme` and `tolerated` (whether this chemistry would amplify);
#'   `f_penalty`/`r_penalty` in the returned table are the forward and reverse
#'   primer's share of `penalty`.
#' @param pam optional PAM / protospacer-flanking site to check for `"crispr"`
#'   (a `TTTV` string, e.g. `"TTTA"`, or `"PFS"` for a Cas13 guide).  An exact
#'   case-insensitive match is the only tolerated state; a mismatch is a hard
#'   loss of target (never a tolerated primer mismatch).
#' @return tibble `chrom`, `start`, `end`, `width`, `n_mismatch`,
#'   `mismatch_positions`, `f_mismatch_positions`, `r_mismatch_positions`,
#'   `three_prime_mismatch`, `seed_mismatch`, `three_prime_terminal`, `penalty`,
#'   `f_penalty`, `r_penalty`, `score` (one row per product).
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' f <- as.character(Biostrings::subseq(genome[["Chr4"]], 1000, 1019))
#' r <- as.character(Biostrings::reverseComplement(Biostrings::subseq(genome[["Chr4"]], 1180, 1199)))
#' in_silico_pcr(f, r, genome)
#' in_silico_pcr(f, r, genome, scheme = "rpa")
#' @family specificity
#' @export
in_silico_pcr <- function(f_primer, r_primer, genome, max_mismatch = NULL, BPPARAM = NULL,
                          scheme = NULL, pam = NULL) {
  policy <- if (is.null(scheme)) NULL else specificity_policy(scheme)
  radius <- if (is.null(max_mismatch)) {
    if (is.null(policy)) 0L else policy$max_mismatch
  } else {
    check_scan_radius(max_mismatch, "max_mismatch")
  }
  # the reverse primer binds the minus strand, so its reverse complement is the
  # plus-strand pattern to search
  rev <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(r_primer)))
  # the score ladder is the scheme's; without a scheme it is pcr's (the
  # historical weights), so the bare call keeps its scoring
  scoring <- policy$mode %||% "pcr"
  chroms <- genome_chroms(genome)
  # The whole genome goes in once.  `primer_hit_table()` already loops over
  # `genome_chroms()`, so calling it inside the per-chromosome loop below
  # rescanned every chromosome n_chrom times AND paired hits from different
  # chromosomes into phantom products (each scan iteration relabelled the same
  # pool with its own chromosome).  `dev/specificity-oracle.R` enumerates the
  # answer and caught it: one real product came back as twelve.
  fh_all <- primer_hit_table(f_primer, genome, radius, scoring, side = "f")
  rh_all <- primer_hit_table(rev, genome, radius, scoring, side = "r")
  out <- pcr_products_from_hits(fh_all, rh_all, chroms, f_primer, r_primer, BPPARAM)
  out <- finish_pcr_products(out, policy)
  warn_near_products(f_primer, r_primer, genome, radius, out, BPPARAM)
  out
}

#' Pair forward and reverse hit tables into in-silico PCR products
#'
#' The pairing half of [in_silico_pcr()], factored out so a caller that scanned
#' several primers in one pass can pair each pair from its slice of that one
#' table (`pcr_products_from_hits()`).
#' @param fh_all forward-primer hit table (`primer_hit_table(side = "f")` on the
#'   primer itself).
#' @param rh_all reverse-primer hit table, searched side (`side = "r"` on the
#'   reverse complement).
#' @param chroms genome sequence names to scan.
#' @param f_primer,r_primer the primers, for the minimum product width.
#' @param BPPARAM optional `BiocParallel` backend.
#' @return the product tibble without `score`/scheme columns (see
#'   `finish_pcr_products()`).
#' @keywords internal
#' @noRd
pcr_products_from_hits <- function(fh_all, rh_all, chroms, f_primer, r_primer, BPPARAM = NULL) {
  min_width <- max(nchar(f_primer), nchar(r_primer))
  scan_one <- function(nm) {
    worker_biostrings()
    fh <- fh_all[fh_all$chrom == nm, , drop = FALSE]
    rh <- rh_all[rh_all$chrom == nm, , drop = FALSE]
    if (nrow(fh) == 0L || nrow(rh) == 0L) {
      return(NULL)
    }
    # Pair every forward start with a downstream reverse end.  A whole-genome
    # exact scan returns ~1.6e5 hits per primer, so the full `outer()` is
    # 1e11 cells (99.6 GB) and cannot be built; a product also needs
    # `width >= min_width`, i.e. a reverse end at least `min_width - 1` past the
    # forward start.  So each forward start is matched against the slice of
    # `rh` inside `[start + min_width - 1, last_end]` - binary search, never a
    # cross product.  Distinct end positions are enough: two hits with the same
    # end differ only in their penalty/positions, and the row kept below is the
    # cheapest one anyway.
    r_end <- rh$end
    r_ord <- order(r_end)
    r_end <- r_end[r_ord]
    rh <- rh[r_ord, , drop = FALSE]
    n_r <- nrow(rh)
    # `findInterval(x, v)` is the count of `v <= x`, so `n_of` right per start:
    # `[start + min_width - 1, n_r]` is the slice of rows that can close a
    # product, and it is never the whole matrix
    lo <- findInterval(fh$start + min_width - 2L, r_end) + 1L
    n_of <- pmax(0L, n_r - lo + 1L)
    f_i <- rep(seq_len(nrow(fh)), n_of)
    r_i <- sequence(n_of) + rep(lo - 1L, n_of)
    if (!length(f_i)) {
      return(NULL)
    }
    f <- fh[f_i, , drop = FALSE]
    r <- rh[r_i, , drop = FALSE]
    out <- tibble::tibble(
      chrom = nm,
      start = f$start,
      end = r$end,
      width = r$end - f$start + 1L,
      n_mismatch = pmax(f$n_mismatch, r$n_mismatch),
      mismatch_positions = lapply(seq_len(nrow(f)), function(i) {
        sort(unique(c(
          f$mismatch_positions[[i]], r$mismatch_positions[[i]],
          recursive = TRUE
        )))
      }),
      f_mismatch_positions = f$mismatch_positions,
      r_mismatch_positions = r$mismatch_positions,
      three_prime_mismatch = f$three_prime_mismatch + r$three_prime_mismatch,
      seed_mismatch = f$seed_mismatch + r$seed_mismatch,
      three_prime_terminal = f$three_prime_terminal | r$three_prime_terminal,
      f_penalty = f$penalty,
      r_penalty = r$penalty,
      penalty = f$penalty + r$penalty
    )
    # one product per (forward start, reverse end): a primer can hit the same
    # place with several substitutions, and only the best-aligned hit is real.
    # Order first (cheapest penalty wins), then keep the first of each pair - a
    # `!duplicated()` before the sort would keep the wrong row.
    out <- out[out$width >= min_width, , drop = FALSE]
    out <- dplyr::arrange(out, .data$penalty, .data$start, .data$end)
    out <- dplyr::distinct(out, .data$chrom, .data$start, .data$end, .keep_all = TRUE)
    if (!nrow(out)) NULL else out
  }
  res <- if (is.null(BPPARAM)) {
    lapply(chroms, scan_one)
  } else {
    bp_lapply(chroms, scan_one, BPPARAM)
  }
  purrr::list_rbind(res) %||% empty_pcr_products()
}

#' Virtual PCR for many primer pairs with one genome pass per radius
#'
#' The cost of a scan is the number of passes over the genome, not the number of
#' primers in a pass (see `pattern_hits_batch()`), and a design asks for several
#' pairs at the same radius - the two allele-specific pairings of an assay, the
#' `pair_screen` candidates of an assembly.  So every pattern of the pool is
#' searched once (both plus-strand halves: the forward primers as themselves,
#' the reverse primers as their reverse complements), and each pair is then
#' paired from its own slice of the shared table.
#'
#' The result is `in_silico_pcr()`'s, pair for pair and in pool order.  The search
#' is one pattern, and the miss counts come from `pattern_hits_scan()`'s
#' `unambiguous_blocks()` only.  Two things matter here:
#' * a pair is paired only from the hits of its own two primers, so the pairing
#'   is unchanged (a hit of another pair's primer is not a candidate);
#' * `in_silico_pcr_pool()` does not warn by default.  `near = TRUE` takes one
#'   radius-1 pass over every pattern of the pool and warns, which is opt-in
#'   because the pooled near call is expensive: 130.05 s on TAIR10 before the
#'   radius-1 batch stopped building a `PDict` (now `matchPattern`; 193 s is the
#'   same pooled path on a differently loaded machine) and 1.99 s after, against
#'   the 1.14 s warm exact call it repeats (3.63 s as a process's first pool
#'   call). The per-pair radius-1 pass is ~0.61 s, unchanged by the fix, and the
#'   internal assembly screen does not ask for `near`.
#' @param pairs named list of `c(F = , R = )` primer pairs.
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param max_mismatch search radius.
#' @param scheme chemistry policy name, or `NULL` for none.
#' @param near with `max_mismatch = 0`, warn for every pair whose products miss
#'   a <= 1-mismatch product (`warn_near_products()`), from one radius-1 pass
#'   over the whole pool.  `FALSE` is the historical silent pool.
#' @return a named list per pair, each an `in_silico_pcr()` product table.
#' @keywords internal
#' @noRd
in_silico_pcr_pool <- function(pairs, genome, max_mismatch = 0L, scheme = NULL, near = FALSE) {
  if (!length(pairs)) {
    return(list())
  }
  policy <- if (is.null(scheme)) NULL else specificity_policy(scheme)
  scoring <- policy$mode %||% "pcr"
  radius <- check_scan_radius(max_mismatch, "max_mismatch")
  fwd <- vapply(pairs, function(x) x[["F"]], character(1))
  # the reverse primer binds the minus strand, so its reverse complement is the
  # plus-strand pattern to search (as in `in_silico_pcr()`)
  rev <- vapply(pairs, function(x) x[["R"]], character(1)) |>
    lapply(function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))) |>
    unlist(use.names = FALSE)
  hits <- pattern_hits_batch(c(fwd, rev), genome, radius)
  # one hit table per *unique* primer, then the per-pair slice is its own
  tables <- new.env(parent = emptyenv())
  primer_table <- function(pattern, side) {
    key <- paste(side, pattern, sep = "\r")
    if (!exists(key, envir = tables, inherits = FALSE)) {
      tbl <- primer_hit_table_from_hits(pattern, hits[[pattern]], genome, radius, scoring, side = side)
      assign(key, tbl, envir = tables)
    }
    get(key, envir = tables, inherits = FALSE)
  }
  chroms <- genome_chroms(genome)
  out <- stats::setNames(
    lapply(seq_along(pairs), function(i) {
      fh <- primer_table(fwd[[i]], "f")
      rh <- primer_table(rev[[i]], "r")
      out <- pcr_products_from_hits(fh, rh, chroms, fwd[[i]], pairs[[i]][["R"]])
      finish_pcr_products(out, policy)
    }),
    names(pairs)
  )
  if (near && radius == 0L) {
    warn_near_pool(pairs, out, genome, fwd, rev)
  }
  out
}

#' Add the score (and the scheme columns) to a product table
#' @param out product table from `pcr_products_from_hits()`.
#' @param policy the scheme policy, or `NULL`.
#' @keywords internal
#' @noRd
finish_pcr_products <- function(out, policy) {
  if (is.null(out)) {
    out <- empty_pcr_products()
  }
  out$score <- hit_score(out)
  # an explicit ordering: a zero-row table must keep the same column contract
  # (`list_rbind()` of an all-`NULL` list yields a column-less data frame)
  out <- if (nrow(out)) {
    tibble::as_tibble(out)[c(
      "chrom", "start", "end", "width", PCR_MISMATCH_COLUMNS,
      "f_penalty", "r_penalty", "score"
    )]
  } else {
    empty_pcr_products()[
      c("chrom", "start", "end", "width", PCR_MISMATCH_COLUMNS, "f_penalty", "r_penalty")
    ] |>
      tibble::add_column(score = numeric())
  }
  if (!is.null(policy)) {
    out$scheme <- policy$mode
    out$tolerated <- pcr_products_tolerated(out, policy)
  }
  out
}

#' Warn when exact mode hides a one-mismatch product
#'
#' Exact mode must not silently hide a one-mismatch primer product.  Scan one
#' mismatch and warn; callers who want those products in the table pass
#' `max_mismatch = 1` explicitly.  The near scan is unscored (no scheme): every
#' product of it that is not in `out` is a mismatch product.  An IUPAC primer is
#' searched exactly only, so it has no near scan (and must not reach the
#' unambiguous-primer guard of `primer_hit_table()`).
#' @param f_primer,r_primer the pair.
#' @param genome the genome.
#' @param radius the radius `out` was scanned at.
#' @param out the exact-mode product table.
#' @param BPPARAM optional `BiocParallel` backend.
#' @keywords internal
#' @noRd
warn_near_products <- function(f_primer, r_primer, genome, radius, out, BPPARAM = NULL) {
  if (radius != 0L) {
    return(invisible(NULL))
  }
  near <- near_products(f_primer, r_primer, genome, BPPARAM)
  warn_near_extra(near, out)
  invisible(NULL)
}

#' The near pass of one pair (`NULL` when the pair cannot take one)
#'
#' An IUPAC primer is searched exactly only, so it has no near scan (and must not
#' reach the unambiguous-primer guard of `primer_hit_table()`).
#' @param f_primer,r_primer the pair.
#' @param genome the genome.
#' @param BPPARAM optional `BiocParallel` backend.
#' @return the radius-1 product table, or `NULL` for a degenerate pair.
#' @keywords internal
#' @noRd
near_products <- function(f_primer, r_primer, genome, BPPARAM = NULL) {
  if (grepl("[^ACGT]", f_primer) || grepl("[^ACGT]", r_primer)) {
    return(NULL)
  }
  in_silico_pcr(f_primer, r_primer, genome, max_mismatch = 1L, BPPARAM = BPPARAM)
}

#' Warn when a pair's near pass found a product the exact table does not hold
#' @param near the radius-1 product table, or `NULL`.
#' @param out the exact-mode product table.
#' @keywords internal
#' @noRd
warn_near_extra <- function(near, out) {
  if (is.null(near)) {
    return(invisible(NULL))
  }
  extra <- dplyr::anti_join(near, out, by = c("chrom", "start", "end"))
  if (nrow(extra)) {
    msg <- gettext(paste0(
      "%d near-match in-silico PCR product(s) with <= 1 primer mismatch ",
      "were found; exact mode ignores them. Set max_mismatch = 1 to include them."
    ))
    cli::cli_warn(sprintf(msg, nrow(extra)))
  }
  invisible(NULL)
}

#' The near pass of `warn_near_products()` for every pair of a pool, in one pass
#'
#' One radius-1 scan of the whole pool replaces the per-pair
#' `in_silico_pcr()` each `warn_near_products()` would issue, and the warnings
#' come out pair by pair, in pool order.
#' @param pairs the pool's pairs.
#' @param out the pool's exact-mode product tables.
#' @param genome the genome.
#' @param fwd,rev the pool's plus-strand patterns (see `in_silico_pcr_pool()`).
#' @keywords internal
#' @noRd
warn_near_pool <- function(pairs, out, genome, fwd, rev) {
  searchable <- vapply(seq_along(pairs), function(i) {
    !grepl("[^ACGT]", fwd[[i]]) && !grepl("[^ACGT]", pairs[[i]][["R"]])
  }, logical(1), USE.NAMES = FALSE)
  if (!any(searchable)) {
    return(invisible(NULL))
  }
  near <- in_silico_pcr_pool(pairs[searchable], genome, max_mismatch = 1L)
  for (k in seq_along(which(searchable))) {
    warn_near_extra(near[[k]], out[[which(searchable)[[k]]]])
  }
  invisible(NULL)
}

#' Zero-row in-silico PCR table with the full (stable) column contract
#' @keywords internal
#' @noRd
empty_pcr_products <- function() {
  out <- tibble::tibble(
    chrom = character(), start = integer(), end = integer(), width = integer()
  )
  out[["n_mismatch"]] <- integer()
  out[["mismatch_positions"]] <- list(integer())
  out[["f_mismatch_positions"]] <- list(integer())
  out[["r_mismatch_positions"]] <- list(integer())
  out[["three_prime_mismatch"]] <- integer()
  out[["seed_mismatch"]] <- integer()
  out[["three_prime_terminal"]] <- logical()
  out[["penalty"]] <- numeric()
  out[["f_penalty"]] <- numeric()
  out[["r_penalty"]] <- numeric()
  out
}

#' Did a pair of mismatched primer hits still amplify under a policy?
#'
#' The chemistry's caps, applied per product: at most `max_mismatch` in total,
#' at most `seed_max` in either seed, at most `three_prime_max` in the last
#' `THREE_PRIME_WINDOW` nt of either primer.  A 3'-terminal mismatch is a hard
#' failure (the primer cannot be extended) for every mode.  The PAM is *not*
#' graded here: it is a licensing element, not a mismatch, so
#' `specificity_verdict()` applies it once to the whole guide.  A policy that
#' carries a central window (a CRISPR guide length whose middle is intolerant)
#' makes a mismatch inside that window the `"central"` violation.
#' @keywords internal
#' @noRd
pcr_products_tolerated <- function(products, policy) {
  if (!nrow(products)) {
    return(logical())
  }
  ok <- products$n_mismatch <= policy$max_mismatch &
    products$seed_mismatch <= policy$seed_max &
    products$three_prime_mismatch <= policy$three_prime_max &
    !products$three_prime_terminal
  ok <- ok & !central_violated(products, policy)
  ok
}

#' Do a product's mismatches fall in the policy's intolerant central window?
#'
#' Some CRISPR guide lengths tolerate a mismatch only away from the middle of
#' the spacer (`spacer_len_central`).  The window is graded on
#' `f_mismatch_positions` -- the guide's own 5'->3' frame, 1-based from the
#' spacer's 5' end, exactly what `mismatch_profile()` reports.  The untyped
#' `mismatch_positions` column is the union of the forward and reverse primer
#' positions and mirrors one guide mismatch `p` to `21 - p`, so grading it
#' would also flag the guide positions past the middle; it is only a fallback
#' for a table without the per-primer columns.  A policy without a resolved
#' window (any non-CRISPR scheme, an unlisted guide length, or the unresolved
#' policy table) tolerates every position.
#' @keywords internal
#' @noRd
central_violated <- function(products, policy) {
  window <- policy[["spacer_len_central"]]
  has_window <- is.numeric(window) && length(window) == 2L
  guide_pos <- if ("f_mismatch_positions" %in% names(products)) {
    products$f_mismatch_positions
  } else {
    products[["mismatch_positions"]]
  }
  if (!has_window || is.null(guide_pos)) {
    return(rep(FALSE, nrow(products)))
  }
  vapply(
    guide_pos,
    function(pos) any(pos >= window[[1L]] & pos <= window[[2L]]),
    logical(1)
  )
}

#' Forward/reverse primer pairs for the two allele-specific virtual PCRs
#'
#' Uses the template-matched AS primers (no deliberate stabiliser), so the pair
#' answers "is the target site itself amplifiable?": on the plus strand the AS
#' primer is the forward primer, on the minus strand it is the reverse one.  A
#' single source of truth shared by [assay_amplicons()] and the genome ideogram.
#' @keywords internal
#' @noRd
as_pcr_pairs <- function(assay) {
  aa <- assay$assay[1L, ]
  ref_len <- nchar(assay$snp$ref)
  # ASQ puts the variant at the penultimate base, so the shipped AS primers and
  # their anchor carry the assay's offset (KASP's is 0).
  voff <- as.integer(assay$as_variant_offset %||% 0L)
  anchor <- as_anchor(assay$template$snp_i, ref_len, aa$strand, voff)
  seed <- c(
    AS_ref = make_as_primer(assay$template$seq, anchor, aa$as_len, assay$snp$ref,
      strand = aa$strand, variant_offset = voff
    ),
    AS_alt = make_as_primer(assay$template$seq, anchor, aa$as_len, assay$snp$alt,
      strand = aa$strand, variant_offset = voff
    )
  )
  purrr::imap(seed, function(as_p, nm) {
    if (identical(aa$strand, "+")) c(F = as_p, R = aa$com) else c(F = aa$com, R = as_p)
  })
}

#' Products predicted for each allele-specific primer against the common one
#'
#' On the reference genome the REF primer must give exactly one product; the
#' ALT primer (the other base / the indel) should give none.  Uses the
#' template-matched primers (no deliberate stabiliser), like [specificity_seed()].
#' @param assay assembled assay (from [assemble_kasp_assay()]).
#' @param genome DNAStringSet.
#' @param max_mismatch mismatches tolerated in the primer match.
#' @param BPPARAM optional `BiocParallel` backend (see [in_silico_pcr()]).
#' @param products precomputed products from `assay_pcr_products()`; when given,
#'   the genome is not rescanned (the ideogram and the amplicon counts then share
#'   one virtual PCR).
#' @return a tibble `primer`, `n_products`, `width`, one row per AS primer.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci[1, ])
#' assay_amplicons(assay, genome)
#' }
#' @family specificity
#' @export
assay_amplicons <- function(assay, genome, max_mismatch = 0L, BPPARAM = NULL, products = NULL) {
  products <- products %||%
    assay_pcr_products(assay, genome, max_mismatch = max_mismatch, BPPARAM = BPPARAM)
  assay_amplicon_counts(products)
}

#' Every in-silico PCR product of an assay's two allele-specific pairings
#'
#' One virtual-PCR pass per AS primer (`as_pcr_pairs()`); computed once at design
#' time and reused for both the amplicon counts and the genome ideogram, so the
#' genome is scanned once instead of twice.
#' @keywords internal
assay_pcr_products <- function(assay, genome, max_mismatch = 0L, BPPARAM = NULL, near = FALSE) {
  # both allele-specific pairings share the COM primer and differ only in the AS
  # one, so one grouped scan answers both (two genome passes instead of four)
  in_silico_pcr_pool(as_pcr_pairs(assay), genome, max_mismatch = max_mismatch, near = near)
}

#' Turn per-primer product tables into the amplicon summary
#' @keywords internal
#' @noRd
assay_amplicon_counts <- function(products) {
  purrr::imap(products, function(prods, nm) {
    tibble::tibble(
      primer = nm, n_products = nrow(prods),
      width = if (nrow(prods)) max(prods$width) else NA_integer_
    )
  }) |>
    purrr::list_rbind()
}
