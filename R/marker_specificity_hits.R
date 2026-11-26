# =============================================================================
#  marker_specificity_hits.R  --  coordinate hit tables and off-targets
# =============================================================================
#  Split out of marker_specificity.R (T7, cycle2/tidy-20261126).  Turns a
#  pattern into per-chromosome hits (`primer_hit_table()`, the batch half
#  `primer_hit_table_from_hits()`), exposes the coordinate-level table
#  (`exact_hit_table()`) and grades an off-target profile.  Function names are
#  unchanged by the split.
# =============================================================================

#' Per-primer hit table for one primer against every chromosome
#'
#' One `matchPattern()` pass per primer (plain-ACGT patterns are scanned with
#' `fixed = TRUE`, IUPAC ones with `fixed = FALSE`) with
#' the mismatch profile of every hit computed against its genomic footprint.
#' The search is exact-or-radius as asked: a mismatched replacement character in
#' the pattern is the only way a terminal mismatch can be scanned here, because
#' a base-space radius over the pattern cannot see the pattern's own last bases.
#' Mismatch degree therefore comes from `mismatch_profile()` on the footprint,
#' never from the search window.
#' `pattern` is the plus-strand pattern searched; `side` says which primer it
#' is, so the reported positions and the 3'-terminal flag are in the primer's
#' own 5'->3' orientation (the reverse primer is searched as its reverse
#' complement, whose plus-strand positions are mirrored back).
#' @param pattern the plus-strand pattern to search (IUPAC-aware).
#' @param genome DNAStringSet.
#' @param max_mismatch search radius.
#' @param scheme chemistry.
#' @param seed_len seed length override (`NA` = the policy's).
#' @param side `"f"` (forward) or `"r"` (reverse complement searched).
#' @keywords internal
#' @noRd
primer_hit_table <- function(pattern, genome, max_mismatch, scheme, seed_len = NA_integer_,
                             side = c("f", "r")) {
  # An IUPAC primer cannot be scanned with a mismatch radius (the radius is a
  # base-space distance, not a code distance); scoring it exactly is a degraded,
  # not a faithful, answer, so say so rather than only raising the error.
  if (max_mismatch > 0L && grepl("[^ACGT]", pattern)) {
    rlang::warn(gettext(paste(
      "A degenerate primer cannot be scanned with `max_mismatch > 0`;",
      "only its exact (IUPAC-aware) matches are scored."
    )), class = "tidymarkers_degenerate_exact")
  }
  side <- match.arg(side)
  hit_tables_per_chrom(
    pattern, genome, genome_chroms(genome), max_mismatch, scheme, seed_len, side,
    hits = NULL
  )
}

#' Hit table of one primer from hits already found (the unsearched half of
#' `primer_hit_table()`)
#'
#' `primer_hit_table()` is "one search per chromosome, then profile every hit".
#' A caller that searched several primers in one pass (`pattern_hits_batch()`)
#' has the second half to do only: this profiles the hits it was handed.
#' @param pattern the plus-strand pattern searched (IUPAC-aware).
#' @param hits the per-chromosome hits from `pattern_hits_batch()`, named by
#'   chromosome.
#' @param genome `DNAStringSet` / `marker_genome`, read for the hit footprints.
#' @param max_mismatch search radius the hits were found at (for the degenerate
#'   warning, so the batch and the per-primer path complain alike).
#' @param scheme chemistry.
#' @param seed_len seed length override (`NA` = the policy's).
#' @param side `"f"` or `"r"` (see `primer_hit_table()`).
#' @keywords internal
#' @noRd
primer_hit_table_from_hits <- function(pattern, hits, genome, max_mismatch, scheme,
                                       seed_len = NA_integer_, side = c("f", "r")) {
  if (max_mismatch > 0L && grepl("[^ACGT]", pattern)) {
    rlang::warn(gettext(paste(
      "A degenerate primer cannot be scanned with `max_mismatch > 0`;",
      "only its exact (IUPAC-aware) matches are scored."
    )), class = "tidymarkers_degenerate_exact")
  }
  side <- match.arg(side)
  hit_tables_per_chrom(
    pattern, genome, names(hits), max_mismatch, scheme, seed_len, side,
    hits = hits
  )
}

#' Build one primer's hit table, searching only the chromosomes it was not given
#' @param pattern the plus-strand pattern.
#' @param genome `DNAStringSet` / `marker_genome`.
#' @param chroms chromosomes to report, in order.
#' @param max_mismatch search radius.
#' @param scheme chemistry.
#' @param seed_len seed length override.
#' @param side `"f"` or `"r"`.
#' @param hits per-chromosome hits already found, or `NULL` to search them.
#' @keywords internal
#' @noRd
hit_tables_per_chrom <- function(pattern, genome, chroms, max_mismatch, scheme, seed_len, side,
                                 hits = NULL) {
  scan_one <- function(nm) {
    worker_biostrings()
    subject <- genome_seq(genome, nm)
    found <- if (is.null(hits)) {
      exact_pattern_hits(pattern, subject, max_mismatch)
    } else {
      hits[[nm]] %||% IRanges::IRanges()
    }
    if (length(found) == 0L) {
      return(NULL)
    }
    # `exact_pattern_hits()` returns an `XStringViews` on the no-N path and a
    # plain `IRanges` on the N path, so `as.character(hits)` was the DNA in one
    # case and "start-end" in the other; take the bases from the subject.
    footprint <- as.character(Biostrings::extractAt(subject, found))
    prof <- hit_profile(pattern, footprint, scheme, seed_len, side == "r")
    tibble::tibble(
      chrom = nm,
      start = IRanges::start(found),
      end = IRanges::end(found),
      n_mismatch = prof$n_mismatch,
      mismatch_positions = prof$mismatch_positions,
      three_prime_mismatch = prof$three_prime_mismatch,
      seed_mismatch = prof$seed_mismatch,
      terminal_mismatch = prof$terminal_mismatch,
      three_prime_terminal = prof$three_prime_terminal,
      penalty = prof$penalty
    )
  }
  res <- lapply(chroms, scan_one)
  out <- purrr::list_rbind(res)
  if (is.null(out)) {
    out <- tibble::tibble(
      chrom = character(), start = integer(), end = integer(),
      n_mismatch = integer(), mismatch_positions = list(integer()),
      three_prime_mismatch = integer(), seed_mismatch = integer(),
      terminal_mismatch = integer(), three_prime_terminal = logical(), penalty = numeric()
    )
  }
  out
}

#' Exact genome hits with coordinates
#'
#' The coordinate-level counterpart of [specificity_seed()]: instead of a count,
#' it returns where every exact (or radius-bounded) occurrence sits.  Both
#' strands are searched - the oligo itself on the plus strand and its reverse
#' complement on the minus strand - and every hit keeps the primer's own name,
#' so a table of oligos comes back as one tidy set of loci.
#'
#' Patterns are grouped by width and searched with `Biostrings::PDict` /
#' `matchPDict()` (the precomputed-dictionary fast path, one pass per dictionary
#' per chromosome).  A degenerate (IUPAC) oligo, or a pattern shorter than 6 nt
#' with `max_mismatch > 0` (the mismatch tree's minimum width), is searched with
#' `exact_pattern_hits()` instead; an ambiguous genomic base is never a hit.
#'
#' @param sequences named character vector of oligos (a `DNAStringSet` is also
#'   accepted).  Unnamed input is labelled by its own sequence.
#' @param genome a `marker_genome` (streamed one chromosome at a time) or a
#'   `DNAStringSet`.
#' @param max_mismatch search radius (default 0, exact).
#' @param BPPARAM optional `BiocParallel` backend; the scan is split by genome
#'   sequence.
#' @param scheme chemistry for the mismatch profile (see [specificity_policy()]).
#' @param seed_len seed length override (`NA` = the policy's).
#' @return a tibble `primer`, `sequence` (the oligo as supplied, 5'->3'),
#'   `chrom`, `start`, `end`, `width`, `strand`
#'   (`strand` is `"+"` for a plus-strand occurrence and `"-"` for one on the
#'   minus strand), then the seven mismatch-profile fields of
#'   `mismatch_profile()` in the primer's own 5'->3' frame (`n_mismatch`,
#'   `mismatch_positions`, `three_prime_mismatch`, `seed_mismatch`,
#'   `terminal_mismatch`, `three_prime_terminal`, `penalty`), one row per hit.
#'   `start`/`end` are 1-based inclusive genomic coordinates in the genome's own
#'   forward orientation.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' primer <- as.character(Biostrings::subseq(genome[["Chr4"]], 1000, 1019))
#' exact_hit_table(c(AS = primer), genome)
#' @family specificity
#' @export
exact_hit_table <- function(sequences, genome, max_mismatch = 0L, BPPARAM = NULL,
                            scheme = "pcr", seed_len = NA_integer_) {
  if (inherits(sequences, "XStringSet")) {
    sequences <- as.character(sequences)
  }
  if (!is.character(sequences) || !length(sequences)) {
    abort_msg("`sequences` must be a non-empty character vector of oligos.")
  }
  max_mismatch <- as.integer(max_mismatch)
  if (is.na(max_mismatch) || max_mismatch < 0L) {
    abort_msg("`max_mismatch` must be a non-negative integer.")
  }
  nm <- names(sequences)
  if (is.null(nm)) {
    nm <- unname(sequences)
  }
  # one dictionary per width: PDict needs equal-length patterns.  A degenerate
  # oligo cannot enter a PDict, and its mismatch tree needs width >= 6, so those
  # patterns take the per-pattern scan instead.
  groups <- split(seq_along(sequences), nchar(sequences))
  dicts <- lapply(groups, function(idx) {
    strict <- grepl("^[ACGT]+$", sequences[idx])
    dict_ok <- strict & (max_mismatch == 0L | nchar(sequences[idx][[1L]]) >= 6L)
    pd <- NULL
    if (any(dict_ok)) {
      pos <- idx[dict_ok]
      fwd <- sequences[pos]
      rev <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(fwd)))
      pd <- list(
        pdict = Biostrings::PDict(
          Biostrings::DNAStringSet(unname(c(fwd, rev))),
          max.mismatch = max_mismatch
        ),
        idx = c(pos, pos),
        strand = rep(c("+", "-"), each = length(pos)),
        pat = unname(c(fwd, rev))
      )
    }
    list(pdict = pd, scan = idx[!dict_ok])
  })
  hit_rows <- function(primer, sequence, chrom, ranges, strand, pattern, subject, scheme, seed_len) {
    if (!length(ranges)) {
      return(NULL)
    }
    footprints <- as.character(Biostrings::extractAt(subject, ranges))
    prof <- hit_profile(pattern, footprints, scheme, seed_len, identical(strand, "-"))
    tibble::tibble(
      primer = primer,
      sequence = sequence,
      chrom = chrom,
      start = IRanges::start(ranges),
      end = IRanges::end(ranges),
      width = IRanges::width(ranges),
      strand = strand,
      n_mismatch = prof$n_mismatch,
      mismatch_positions = prof$mismatch_positions,
      three_prime_mismatch = prof$three_prime_mismatch,
      seed_mismatch = prof$seed_mismatch,
      terminal_mismatch = prof$terminal_mismatch,
      three_prime_terminal = prof$three_prime_terminal,
      penalty = prof$penalty
    )
  }
  scan_chrom <- function(chrom) {
    worker_biostrings()
    subject <- genome_seq(genome, chrom)
    rows <- list()
    for (g in dicts) {
      if (!is.null(g$pdict)) {
        found <- Biostrings::matchPDict(
          g$pdict$pdict, subject,
          max.mismatch = max_mismatch
        )
        for (k in seq_along(g$pdict$idx)) {
          i <- g$pdict$idx[[k]]
          rows[[length(rows) + 1L]] <- hit_rows(
            nm[i], sequences[[i]], chrom, found[[k]], g$pdict$strand[[k]],
            g$pdict$pat[[k]], subject, scheme, seed_len
          )
        }
      }
      for (i in g$scan) {
        for (strand in c("+", "-")) {
          pattern <- sequences[[i]]
          if (strand == "-") {
            pattern <- as.character(
              Biostrings::reverseComplement(Biostrings::DNAString(pattern))
            )
          }
          rows[[length(rows) + 1L]] <- hit_rows(
            nm[[i]], sequences[[i]], chrom,
            exact_pattern_hits(pattern, subject, max_mismatch), strand,
            pattern, subject, scheme, seed_len
          )
        }
      }
    }
    out <- purrr::list_rbind(rows)
    if (is.null(out) || !nrow(out)) NULL else out
  }
  res <- if (is.null(BPPARAM)) {
    lapply(genome_chroms(genome), scan_chrom)
  } else {
    bp_lapply(genome_chroms(genome), scan_chrom, BPPARAM)
  }
  out <- purrr::list_rbind(res)
  if (is.null(out) || !ncol(out)) {
    out <- tibble::tibble(
      primer = character(), sequence = character(), chrom = character(),
      start = integer(),
      end = integer(), width = integer(), strand = character(),
      n_mismatch = integer(), mismatch_positions = list(integer()),
      three_prime_mismatch = integer(), seed_mismatch = integer(),
      terminal_mismatch = integer(), three_prime_terminal = logical(),
      penalty = numeric()
    )
  }
  out
}

#' Grade an off-target locus by its 3'-end mismatches
#'
#' The 3'-aware post-filter for an aligner's loci.  `exact_hit_table()` and any
#' external aligner report *where* a hit sits but not *where* its mismatches
#' sit; a short-read aligner cannot tell a 3'-terminal mismatch from an inner
#' one.  This fetches each footprint and profiles it with `mismatch_profile()`,
#' so a 3'-terminal or 3'-end mismatch is graded without re-running a search.
#'
#' The `tolerated` verdict is the single-primer form of
#' `pcr_products_tolerated()`, built from the same policy fields
#' (`max_mismatch`, `seed_max`, `three_prime_max`) and the same hard rule that a
#' 3'-terminal mismatch always fails; the two stay numerically identical.  The
#' thresholds are heuristics (`calibrated = FALSE`), not calibrated quantities.
#'
#' @param hits a data frame/tibble with at least `primer`, `chrom`, `start`,
#'   `end`, `strand`; the shape [exact_hit_table()] returns is accepted and its
#'   profile columns are recomputed, not trusted.
#' @param primers named character vector, name -> oligo 5'->3', or `NULL` to
#'   derive it from a `sequence` column on `hits` (as [exact_hit_table()]
#'   returns).  Supply it explicitly for an external aligner's loci, which do
#'   not carry the oligo.  Every `hits$primer` must be a name.
#' @param genome a `DNAStringSet` or a `marker_genome`.
#' @param scheme chemistry (see [specificity_policy()]).
#' @param seed_len seed length override (`NA` = the policy's).
#' @return `hits` with the seven profile columns of `mismatch_profile()` appended
#'   (in the primer's own 5'->3' frame) and a logical `tolerated`.
#' @examples
#' p <- "GAAGGTGACCAAGTTCATGCT"
#' g <- Biostrings::DNAStringSet(c(chr1 = paste0(strrep("A", 10), p, strrep("A", 10))))
#' h <- exact_hit_table(c(AS = p), g)
#' offtarget_profile(h, genome = g)
#' @family specificity
#' @export
offtarget_profile <- function(hits, primers = NULL, genome, scheme = "pcr", seed_len = NA_integer_) {
  need_cols(hits, c("primer", "chrom", "start", "end", "strand"), arg = "hits")
  if (is.null(primers)) {
    if (!"sequence" %in% names(hits)) {
      abort_msg(
        "`hits` carries no `sequence` column, so `primers` must be supplied."
      )
    }
    key <- as.character(hits$primer)
    seq <- as.character(hits$sequence)
    if (anyNA(seq)) {
      abort_msg("`hits$sequence` must not contain NA.")
    }
    pairs <- unique(data.frame(name = key, seq = seq, stringsAsFactors = FALSE))
    if (anyDuplicated(pairs$name)) {
      abort_msg(
        "`hits$sequence` gives one primer name more than one oligo; supply `primers` explicitly."
      )
    }
    primers <- stats::setNames(pairs$seq, pairs$name)
  }
  if (!is.character(primers) || is.null(names(primers)) || any(is.na(names(primers)) | names(primers) == "")) {
    abort_msg("`primers` must be a named character vector of oligos.")
  }
  primer_name <- as.character(hits$primer)
  strand <- as.character(hits$strand)
  bad_strand <- !is.na(strand) & !(strand %in% c("+", "-"))
  if (any(bad_strand)) {
    abort_msg(
      "`hits$strand` must be {.val {c('+', '-')}}; got {.val {unique(strand[bad_strand])}}."
    )
  }
  start <- hits$start
  if (!is.numeric(start) || anyNA(start) || any(start < 1) || any(start != floor(start))) {
    abort_msg("`hits$start` must be a positive integer.")
  }
  missing_primers <- setdiff(unique(primer_name), names(primers))
  if (length(missing_primers)) {
    abort_msg(
      "Every `hits$primer` must be a name in `primers`; missing: {.val {missing_primers}}."
    )
  }
  primer_seq <- unname(primers[primer_name])
  width <- as.integer(hits$end) - as.integer(hits$start) + 1L
  bad_width <- width != nchar(primer_seq)
  if (any(bad_width)) {
    i <- which(bad_width)[[1L]]
    cli::cli_abort(c(
      gettext("`hits` must describe gapless equal-width alignments."),
      x = gettext(paste0(
        "Row ", i, " (primer {.val {primer_name[[i]]}}) has width ", width[[i]],
        ", but primer length is ", nchar(primer_seq[[i]]), "."
      ))
    ))
  }
  profiles <- vector("list", nrow(hits))
  for (i in seq_len(nrow(hits))) {
    chrom <- as.character(hits$chrom[[i]])
    end_fetch <- min(as.integer(hits$end[[i]]), genome_length(genome, chrom))
    window <- genome_seq(genome, chrom, as.integer(hits$start[[i]]), end_fetch)
    pattern <- primer_seq[[i]]
    if (strand[[i]] == "-") {
      pattern <- as.character(
        Biostrings::reverseComplement(Biostrings::DNAString(pattern))
      )
    }
    profiles[[i]] <- hit_profile(
      pattern, as.character(window), scheme, seed_len, strand[[i]] == "-"
    )
  }
  n_mismatch <- vapply(profiles, function(p) p$n_mismatch, integer(1))
  mismatch_positions <- lapply(profiles, function(p) p$mismatch_positions[[1L]])
  three_prime_mismatch <- vapply(profiles, function(p) p$three_prime_mismatch, integer(1))
  seed_mismatch <- vapply(profiles, function(p) p$seed_mismatch, integer(1))
  terminal_mismatch <- vapply(profiles, function(p) p$terminal_mismatch, integer(1))
  three_prime_terminal <- vapply(profiles, function(p) p$three_prime_terminal, logical(1))
  penalty <- vapply(profiles, function(p) p$penalty, numeric(1))
  pol <- specificity_policy(scheme)
  tolerated <- n_mismatch <= pol$max_mismatch &
    seed_mismatch <= pol$seed_max &
    three_prime_mismatch <= pol$three_prime_max &
    !three_prime_terminal
  base <- hits[setdiff(names(hits), c(SINGLE_PRIMER_PROFILE_COLUMNS, "tolerated"))]
  base[["n_mismatch"]] <- n_mismatch
  base[["mismatch_positions"]] <- mismatch_positions
  base[["three_prime_mismatch"]] <- three_prime_mismatch
  base[["seed_mismatch"]] <- seed_mismatch
  base[["terminal_mismatch"]] <- terminal_mismatch
  base[["three_prime_terminal"]] <- three_prime_terminal
  base[["penalty"]] <- penalty
  base[["tolerated"]] <- tolerated
  base
}
