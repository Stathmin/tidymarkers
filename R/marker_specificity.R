# =============================================================================
#  marker_specificity.R  --  genome-wide specificity, in-silico PCR
# =============================================================================
# =============================================================================
#  6. Genome-wide specificity (pure R, plus optional BLAST)
# =============================================================================
#
#  Specificity is mismatch- and position-aware, not "exact hit + a 15-nt seed".
#  Every product carries its primer mismatches, where they sit, and a penalty
#  that weights the 3' end and the seed above the middle.  What is *tolerated*
#  is chemistry-specific, so `specificity_policy(scheme)` returns that
#  chemistry's tolerance object and `specificity_verdict()` grades hits against
#  it.  The numbers are heuristics (calibrated = FALSE), see the policy table.
# =============================================================================

#' Count genome-wide near-matches of each primer's 3'-terminal k-mer
#'
#' For the exact case (`max_mismatch = 0`) counts come from `Biostrings::PDict`
#' (`vcountPDict`) - the precomputed-dictionary fast path, one call for all
#' primers.  With mismatches Biostrings has no counting index (`PDict()` builds an
#' `MTB_PDict`, but `vcountPDict()` does not support it), so the scan is done with
#' `vcountPattern` **parallelised over the genome's sequences** (and, per
#' sequence, once for all seeds), which is how this scales to large genomes.
#' Results are summed over both strands and over all genome sequences.
#' @param primers named character vector of oligos.
#' @param genome DNAStringSet.
#' @param k 3'-seed length.
#' @param max_mismatch mismatches allowed in the seed.
#' @param BPPARAM optional `BiocParallel` backend; the mismatch scan is split by
#'   genome sequence, so `MulticoreParam()` uses one worker per chromosome.
#' @return a tibble `primer`, `seed`, `hits` (both-strand, genome-wide counts).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' specificity_seed(c(AS = "GAAGGTGACCAAGTTCATGCT"), genome)
#' }
#' @family KASP specificity
#' @export
specificity_seed <- function(primers, genome, k = SEED_K, max_mismatch = SEED_MISMATCH,
                             BPPARAM = NULL) {
  if (inherits(genome, "marker_genome")) {
    cli::cli_abort(c(gettext("specificity_seed() needs an in-memory DNAStringSet."),
      i = gettext("For an indexed genome use specificity_of()/specificity_blast().")
    ))
  }
  seeds <- unname(substr(primers, nchar(primers) - k + 1L, nchar(primers)))
  if (max_mismatch == 0L) {
    hits <- exact_seed_counts(seeds, genome)
  } else {
    # vcountPattern takes a single pattern, so loop over seeds *per chromosome*;
    # parallelism is over the genome's sequences (one worker per chromosome)
    fwd <- Biostrings::DNAStringSet(seeds)
    rev <- Biostrings::reverseComplement(fwd)
    count_chr <- function(i) {
      worker_biostrings()
      chr <- genome[[i]]
      cf <- vapply(seq_along(fwd), function(j) {
        Biostrings::countPattern(fwd[[j]], chr, max.mismatch = max_mismatch)
      }, numeric(1))
      cr <- vapply(seq_along(rev), function(j) {
        Biostrings::countPattern(rev[[j]], chr, max.mismatch = max_mismatch)
      }, numeric(1))
      cf + cr
    }
    per_chr <- if (is.null(BPPARAM)) {
      lapply(seq_along(genome), count_chr)
    } else {
      bp_lapply(seq_along(genome), count_chr, BPPARAM)
    }
    hits <- Reduce(`+`, per_chr)
  }
  tibble::tibble(primer = names(primers), seed = seeds, hits = as.integer(hits))
}

# Mismatch position weighting and the PCR product contract ---------------------
#
# The last 5 nt are the 3' window the PCR/KASP literature singles out; the
# position ladder (middle < seed < last 5 < terminal) is a heuristic, not a
# calibrated free-energy model - see mismatch_weights().  A mismatch in the
# last THREE_PRIME_TERMINAL nt is a 3'-terminal mismatch (no extension).
THREE_PRIME_WINDOW <- 5L
THREE_PRIME_TERMINAL <- 3L

# IUPAC primer code -> the bases it covers (Biostrings' DNA alphabet)
IUPAC_BASES <- list(
  A = "A", C = "C", G = "G", T = "T",
  R = c("A", "G"), Y = c("C", "T"), S = c("G", "C"), W = c("A", "T"),
  K = c("G", "T"), M = c("A", "C"), B = c("C", "G", "T"),
  D = c("A", "G", "T"), H = c("A", "C", "T"), V = c("A", "C", "G"),
  N = c("A", "C", "G", "T")
)

#' Table: which IUPAC code covers which subject base
#'
#' The primer is the row, the subject base the column, so
#' `cover[primer_code, subject_base]` is one logical answer per position, with no
#' per-base call.  `N` in the subject is covered by every code (an ambiguous
#' genomic base is not evidence of a mismatch), and `N` as a primer code covers
#' every base.  A subject base outside the alphabet (an `NA` past a sequence end)
#' is never covered.
#' @return logical matrix, both dimnames the IUPAC codes.
#' @keywords internal
#' @noRd
iupac_cover_table <- function() {
  codes <- names(IUPAC_BASES)
  cover <- matrix(FALSE, length(codes), length(codes), dimnames = list(codes, codes))
  for (code in codes) {
    cover[code, IUPAC_BASES[[code]]] <- TRUE
    cover[code, "N"] <- TRUE
  }
  cover
}

# Columns every in-silico PCR product carries, in order, after the legacy four
PCR_MISMATCH_COLUMNS <- c(
  "n_mismatch", "mismatch_positions", "f_mismatch_positions",
  "r_mismatch_positions", "three_prime_mismatch",
  "seed_mismatch", "three_prime_terminal", "penalty"
)

#' Specificity score of an in-silico PCR product (higher is better)
#'
#' `100 - 10 * penalty`: a perfect product scores 100, and the score is directly
#' comparable across products because it is built only from the position-weighted
#' mismatch penalty.  Products carrying `NA` (no profile) score `NA`.
#' @keywords internal
#' @noRd
hit_score <- function(products) {
  if (!nrow(products)) {
    return(numeric())
  }
  round(100 - 10 * products$penalty, 4L)
}

#' Exact both-strand occurrence counts of seeds via a `PDict` dictionary
#' @keywords internal
#' @noRd
exact_seed_counts <- function(seeds, genome) {
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seeds)))
  pats <- c(seeds, rc)
  keep <- !duplicated(pats)
  # vcountPDict returns one row per dictionary element, in the given order
  m <- Biostrings::vcountPDict(Biostrings::PDict(pats[keep]), genome)
  cnt <- rlang::set_names(as.integer(rowSums(m)), pats[keep])
  as.integer(cnt[seeds] + cnt[rc])
}

#' Per-position mismatch weight for a primer hit (heuristic, not calibrated)
#'
#' A mismatch near the primer's 3' terminus destabilises the annealed end far
#' more than one in the middle, because extension starts there; the seed region
#' (the 3'-terminal `seed_len` nt, the mismatch-critical region of an
#' allele-specific primer) sits just behind it.  The ladder is the scheme's own
#' (`weight_seed`, `weight_three_prime`, `weight_terminal` of
#' [specificity_policy()]), so RPA's tolerant filaments are not scored like a
#' proofreading polymerase.  It is a monotone heuristic, not a free-energy
#' calculation.
#' @param n_nt primer length.
#' @param scheme chemistry (see [specificity_policy()]).
#' @param seed_len seed length (default: the policy's).
#' @return numeric vector of length `n_nt`, weight per position (5'->3').
#' @keywords internal
#' @noRd
mismatch_weights <- function(n_nt, scheme = "pcr", seed_len = NULL) {
  pol <- specificity_policy(scheme)
  seed_len <- seed_len %||% pol$seed_len
  n_nt <- as.integer(n_nt)
  if (!n_nt) {
    return(numeric())
  }
  w <- rep(pol$weight_middle, n_nt)
  seed <- seq.int(max(1L, n_nt - seed_len + 1L), n_nt)
  w[seed] <- pol$weight_seed
  tail5 <- seq.int(max(1L, n_nt - pol$three_prime_window + 1L), n_nt)
  w[tail5] <- pol$weight_three_prime
  w[tail_terminal(n_nt)] <- pol$weight_terminal
  w
}

#' Last `THREE_PRIME_TERMINAL` positions of `n_nt` (5'->3', 1-based)
#' @keywords internal
#' @noRd
tail_terminal <- function(n_nt) {
  seq.int(max(1L, n_nt - THREE_PRIME_TERMINAL + 1L), n_nt)
}

#' Mismatch profile of one primer hit against the genomic footprint
#'
# Each mismatching position is counted once, at its heaviest role, so a
# degenerate code plus `max_mismatch` cannot double-count (a seed position is
# not also added as a generic mismatch):
#   terminal (last 3 nt) > last-5 window > seed > middle.
#' Compares the searched (plus-strand) primer pattern with the subject slice(s) it
#' matched.  IUPAC codes in the primer are expanded, so a degenerate position
#' that still covers the subject base is a *match*, not a mismatch; an N in the
#' subject never mismatches (it matches anything).  Mismatch positions are
#' reported 1-based from the 5' end of the primer, i.e. in the primer's own
#' 5'->3' orientation - but `pattern` here is the *plus-strand* pattern, so the
#' reverse primer's positions are mirrored back by [plus_to_primer()] in
#' [primer_hit_table()].
#'
#' Every mismatching position is counted once, at its heaviest role (the masks
#' are nested: terminal <= last-5 window <= seed <= middle), so a degenerate
#' code plus a search radius can never double-count one position.  A terminal
#' mismatch contributes `weight_terminal` each and dominates the sum; the other
#' mismatches contribute their ladder weight divided by one more than the
#' terminal count.
#' @param pattern searched primer as a character string (plus strand).
#' @param subject the genomic slice(s) of the same width (character); one
#'   slice is a scalar and several are a vector, so a whole hit table costs one
#'   call.  A slice shorter than the pattern has `NA` subject bases past its end,
#'   which are never mismatches - exactly what a search near a sequence edge sees.
#' @param scheme chemistry.
#' @param seed_len seed length override (`NA` = the policy's).
#' @return `n_mismatch`, `mismatch_positions` (list, plus-strand),
#'   `three_prime_mismatch`, `seed_mismatch`, `terminal_mismatch`,
#'   `three_prime_terminal`, `penalty`.
#' @keywords internal
#' @noRd
mismatch_profile <- function(pattern, subject, scheme = "pcr", seed_len = NA_integer_) {
  n_nt <- nchar(pattern)
  if (!n_nt) {
    return(list(
      n_mismatch = 0L, mismatch_positions = list(integer()),
      three_prime_mismatch = 0L, seed_mismatch = 0L, terminal_mismatch = 0L,
      three_prime_terminal = FALSE, penalty = 0
    ))
  }
  n_sub <- length(subject)
  # every slice is exactly the pattern's width, so the bases unlist into a vector
  # in hit-then-position order (row = position for the matrix below).  `strsplit`
  # is the whole cost of a genome-wide table (~0.8 s for 183k hits), so it is
  # replaced by fixed-offset slicing whenever every slice has the pattern's width;
  # only a short slice (a search near a sequence end) needs the padded path, and
  # then `strsplit` runs on that one slice.  A padded base is NA, never a mismatch.
  if (all(lengths(subject) == n_nt)) {
    sub <- substring(unlist(subject, use.names = FALSE), rep(seq_len(n_nt), n_sub))
  } else {
    sub <- strsplit(subject, "", fixed = TRUE)
    sub <- unlist(lapply(sub, function(x) {
      n_pad <- max(0L, n_nt - length(x))
      c(x, rep(NA_character_, n_pad))[seq_len(n_nt)]
    }), use.names = FALSE)
  }
  # which pattern code covers which subject base, as a 15 x 15 table built once
  # per call (a set test, not a per-base `%in%`: `%in%` against a multi-base code
  # compares element-wise and would miss N), and an N in the subject always a
  # match.  This is one matrix index for a whole genome-wide hit table: a
  # low-complexity primer returns ~180k hits, and a per-hit call cost seconds.
  cover <- iupac_cover_table()
  pat_code <- match(strsplit(pattern, "", fixed = TRUE)[[1L]], rownames(cover))
  sub_code <- match(sub, rownames(cover))
  sub_code[is.na(sub_code)] <- match("N", rownames(cover))
  covered <- cover[cbind(rep(pat_code, length.out = length(sub_code)), sub_code)]
  # column-major: the bases are in hit-then-position order already
  mism <- matrix(!covered, nrow = n_nt)
  s_len <- if (is.na(seed_len)) specificity_policy(scheme)$seed_len else seed_len
  w <- mismatch_weights(n_nt, scheme, s_len)
  if (n_sub > 1L) {
    return(mismatch_profile_grid(mism, n_sub, w, s_len))
  }
  # an `NA` in `mism` is never a mismatch (an `N` subject base, or a slice past
  # the sequence end).  This R's `which()` treats `NA` as unset and can return
  # `0`, so the test is written so no `NA` reaches it; a single hit is a
  # one-column matrix, so the positions are its TRUE rows
  pos <- which(mism == TRUE)
  n_mm <- length(pos)
  # A mismatch in the last THREE_PRIME_TERMINAL nt stops extension, so its
  # penalty is the scheme's terminal weight times the number of terminal
  # mismatches.  Ordinary mismatches then add their ladder weight *discounted*
  # by the terminal ones: same-role mismatches add up, but no number of
  # elsewhere mismatches can outscore a 3'-terminal one.  A terminal mismatch
  # therefore always dominates, and an exact 3' end is never penalised for a
  # mismatch elsewhere.
  term_w <- max(max(w), 0)
  # positional in the primer's frame here; the masks are logical, so a mismatch
  # is in a region earlier in the ladder unless a later region claims it
  terminal <- pos %in% tail_terminal(n_nt)
  penalty <- sum(terminal) * term_w +
    sum(w[pos[!terminal]]) / (sum(terminal) + 1L)
  list(
    n_mismatch = n_mm,
    mismatch_positions = list(pos),
    three_prime_mismatch = sum(pos %in% seq.int(max(1L, n_nt - THREE_PRIME_WINDOW + 1L), n_nt)),
    seed_mismatch = sum(pos %in% seq.int(max(1L, n_nt - s_len + 1L), n_nt)),
    terminal_mismatch = sum(terminal),
    three_prime_terminal = sum(terminal) > 0L,
    penalty = penalty
  )
}

#' The same profile for many hit slices, in one vectorised pass
#'
#' A genome-wide one-mismatch scan of a low-complexity primer returns tens of
#' thousands of hits, and profiling them one at a time cost ~0.4 s per 200
#' (a projected 400-500 s for a single primer).  The mismatch matrix is n_nt x
#' n_hits, so every per-hit number is a column sum of one masked matrix, and the
#' positions are one `row()`/`col()` split.  Behaviour is identical to the
#' per-hit loop; it is the same masks, not a re-derivation.
#' @param mism logical matrix, `n_nt` rows x one column per hit.
#' @param n_hits number of hit columns.
#' @param w the scheme's per-position weights (5'->3', length `n_nt`).
#' @param s_len seed length.
#' @keywords internal
#' @noRd
mismatch_profile_grid <- function(mism, n_hits, w, s_len) {
  n_nt <- nrow(mism)
  term_w <- max(max(w), 0)
  # `NA` is a covered base (an `N` subject base, or past a sequence end), never a
  # mismatch; this R's `row()`/`col()`/`which()` treat `NA` as unset, so it is
  # dropped once, up front, and the matrix is a plain logical from here on
  mism[is.na(mism)] <- FALSE
  terminal <- row(mism) %in% tail_terminal(n_nt) & mism
  n_term <- colSums(terminal)
  # ordinary mismatches add their ladder weight discounted by the terminal count;
  # the matrix product sums the weight of every non-terminal mismatch per hit
  penalty <- n_term * term_w + (w %*% (mism & !terminal))[1L, ] / (n_term + 1L)
  pos <- row(mism)[mism]
  cols <- col(mism)[mism]
  list(
    n_mismatch = as.integer(colSums(mism)),
    mismatch_positions = split(pos, factor(cols, levels = seq_len(n_hits))),
    three_prime_mismatch = as.integer(colSums(mism[
      seq.int(max(1L, n_nt - THREE_PRIME_WINDOW + 1L), n_nt), ,
      drop = FALSE
    ])),
    seed_mismatch = as.integer(colSums(mism[
      seq.int(max(1L, n_nt - s_len + 1L), n_nt), ,
      drop = FALSE
    ])),
    terminal_mismatch = as.integer(n_term),
    three_prime_terminal = n_term > 0L,
    penalty = as.numeric(penalty)
  )
}

#' Convert a plus-strand mismatch position to the primer's own 5'->3' position
#'
#' The reverse primer is searched as its reverse complement (the plus strand),
#' so a position counted there runs 3'->5' of the primer: 1 on the plus strand
#' is the reverse primer's 3'-terminal base.  Mirroring maps position `p` back
#' to `n_nt - p + 1`, otherwise the report would call the 5' end the 3' end.
#' @param pos integer positions on the searched (plus-strand) pattern.
#' @param n_nt primer length.
#' @param reverse is the searched pattern the reverse complement of the primer?
#' @return integer positions 1-based from the primer's 5' end.
#' @keywords internal
#' @noRd
plus_to_primer <- function(pos, n_nt, reverse) {
  if (reverse) n_nt - pos + 1L else pos
}

#' Does an IUPAC primer code cover a subject base?
#'
#' `N` in the subject matches anything (an ambiguous genomic base is not evidence
#' of a mismatch); a subject base outside `IUPAC_BASES` is always a mismatch.
#' @keywords internal
#' @noRd
iupac_covers <- function(code, base) {
  pat <- IUPAC_BASES[[code]]
  if (is.null(pat) || !base %in% IUPAC_BASES) {
    return(FALSE)
  }
  if (identical(base, "N")) {
    return(TRUE)
  }
  base %in% pat
}

#' Runs of unambiguous bases in a genomic sequence (cached per sequence)
#'
#' A subject `N` must never satisfy a pattern (`Biostrings::matchPattern()` with
#' `fixed = FALSE` treats it as covering every base, so a 27-mer returns ~1.6e5
#' hits on TAIR10 and an in-silico PCR product table of ~1e10 rows).  So the
#' search runs on the unambiguous blocks, and this finds them once per sequence.
#'
#' The ambiguous positions come from the byte codes: `as.integer()` on an XString
#' is one C-level pass and N is 15 in Biostrings' DNA alphabet (A 1, C 2, G 4,
#' T 8).  `letterFrequencyInSlidingView(subject, 1L, "N")` answers the same
#' question but spends ~30 s per genome in S4 dispatch (measured), so it is not
#' used.  The runs are the gaps in `diff()` of the ambiguous positions, which is
#' 0.05 s on 30 Mb.
#'
#' `as.integer()` still allocates a 60 MB integer vector, and every hit table and
#' every `primer_hits()` call on the same genome would pay it again, so the
#' answer is memoised per process.  The key is the sequence's content signature
#' (`digest::digest` on the bytes is not free either, so the key is the string
#' hash of the length and the first/last bases plus the N count) - a *changed*
#' sequence must not reuse a stale block list.  The cache lives in the package
#' namespace and is bounded.
#'
#' @param subject the genomic sequence as an `XString`.
#' @return a list of integer vectors `start`, `end` (1-based, inclusive).
#' @keywords internal
#' @noRd
unambiguous_blocks <- function(subject) {
  n_amb <- Biostrings::letterFrequency(subject, "N")
  if (n_amb == 0L) {
    return(list(start = 1L, end = length(subject)))
  }
  len <- length(subject)
  key <- paste(len, n_amb, as.character(Biostrings::subseq(subject, 1L, 32L)),
    as.character(Biostrings::subseq(subject, max(1L, len - 31L), len)),
    sep = "|"
  )
  cached <- get0(key, envir = .marker_state$blocks, inherits = FALSE)
  if (!is.null(cached)) {
    return(cached)
  }
  ap <- which(as.integer(subject) == 15L)
  gap <- diff(ap) > 1L
  a_start <- ap[c(TRUE, gap)]
  a_end <- ap[c(gap, TRUE)]
  # the clean blocks are the complement of the ambiguous runs
  b_start <- c(1L, a_end + 1L)
  b_end <- c(a_start - 1L, len)
  ok <- b_end >= b_start
  out <- list(start = b_start[ok], end = b_end[ok])
  # a simple bound: the block list of the current genome is what matters
  if (length(.marker_state$blocks) > 16L) {
    rm(list = ls(.marker_state$blocks, all.names = TRUE), envir = .marker_state$blocks)
  }
  assign(key, out, envir = .marker_state$blocks)
  out
}

#' Primer hits that never match an ambiguous genomic base
#'
#' `Biostrings::matchPattern(fixed = FALSE)` treats the pattern as a set of IUPAC
#' codes and the *subject* N as covering every base, so a run of N in the genome
#' returns a hit at every offset in it.  On TAIR10 that is 163,958 spurious hits
#' for a 27-mer, and a 27-mer in a N-run pair with a reverse hit in the same run:
#' the in-silico PCR product table then holds ~1e10 rows and cannot be allocated.
#'
#' An `N` in the subject is not evidence of a binding site, so it must never be a
#' match ([mismatch_profile()] already reports it as a non-mismatch *within* a
#' hit, which is a different question: whether an aligned base counts against the
#' primer).  This searches the unambiguous blocks of the subject separately and
#' shifts each hit back to genomic coordinates, so a hit can never span an
#' ambiguous base.  Blocks shorter than the pattern are skipped, so this is also
#' cheaper than the whole-chromosome scan on a fragmented genome.
#'
#' @param pattern the plus-strand pattern to search (IUPAC-aware).
#' @param subject the genomic sequence as an `XString`.
#' @param max_mismatch search radius.
#' @return an `IRanges` of hit starts/ends in `subject` coordinates.
#' @keywords internal
#' @noRd
exact_pattern_hits <- function(pattern, subject, max_mismatch) {
  width <- nchar(pattern)
  if (Biostrings::letterFrequency(subject, "N") == 0L) {
    return(Biostrings::matchPattern(pattern, subject, max.mismatch = max_mismatch, fixed = FALSE))
  }
  bl <- unambiguous_blocks(subject)
  len <- bl$end - bl$start + 1L
  keep <- len >= width
  if (!any(keep)) {
    return(IRanges::IRanges())
  }
  blocks <- lapply(which(keep), function(i) {
    s <- bl$start[[i]]
    sub <- Biostrings::subseq(subject, s, bl$end[[i]])
    h <- Biostrings::matchPattern(pattern, sub, max.mismatch = max_mismatch, fixed = FALSE)
    if (!length(h)) {
      return(NULL)
    }
    IRanges::IRanges(IRanges::start(h) + s - 1L, width = width)
  })
  do.call(c, c(blocks, list(IRanges::IRanges())))
}

#' Per-primer hit table for one primer against every chromosome
#'
#' One `matchPattern()` pass per primer (`fixed = FALSE`, so IUPAC-aware) with
#' the mismatch profile of every hit computed against its genomic footprint.
#' The search is exact-or-radius as asked: a mismatched replacement character in
#' the pattern is the only way a terminal mismatch can be scanned here, because
#' a base-space radius over the pattern cannot see the pattern's own last bases.
#' Mismatch degree therefore comes from [mismatch_profile()] on the footprint,
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
  orient <- function(prof) {
    prof$mismatch_positions <- lapply(prof$mismatch_positions, function(pos) {
      # an AsIs element keeps the column a vector of integer vectors (a bare list
      # of length 1 would be treated as a column and nested).  A forward hit is
      # already ascending, so only a mirrored one needs sorting
      pos <- plus_to_primer(pos, nchar(pattern), side == "r")
      if (side == "r") pos <- sort(pos)
      structure(as.integer(pos), class = "AsIs")
    })
    prof
  }
  scan_one <- function(nm) {
    worker_biostrings()
    subject <- genome_seq(genome, nm)
    hits <- exact_pattern_hits(pattern, subject, max_mismatch)
    if (length(hits) == 0L) {
      return(NULL)
    }
    footprint <- as.character(hits)
    prof <- orient(mismatch_profile(pattern, footprint, scheme, seed_len))
    tibble::tibble(
      chrom = nm,
      start = IRanges::start(hits),
      end = IRanges::end(hits),
      n_mismatch = prof$n_mismatch,
      mismatch_positions = prof$mismatch_positions,
      three_prime_mismatch = prof$three_prime_mismatch,
      seed_mismatch = prof$seed_mismatch,
      terminal_mismatch = prof$terminal_mismatch,
      three_prime_terminal = prof$three_prime_terminal,
      penalty = prof$penalty
    )
  }
  res <- lapply(genome_chroms(genome), scan_one)
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
#' \dontrun{
#' genome <- ath_genome()
#' in_silico_pcr("GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC", genome)
#' in_silico_pcr("GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC", genome,
#'   scheme = "rpa"
#' )
#' }
#' @family KASP specificity
#' @export
in_silico_pcr <- function(f_primer, r_primer, genome, max_mismatch = NULL, BPPARAM = NULL,
                          scheme = NULL, pam = NULL) {
  policy <- if (is.null(scheme)) NULL else specificity_policy(scheme)
  radius <- if (is.null(max_mismatch)) {
    if (is.null(policy)) 0L else policy$max_mismatch
  } else {
    as.integer(max_mismatch)
  }
  # the reverse primer binds the minus strand, so its reverse complement is the
  # plus-strand pattern to search
  rev <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(r_primer)))
  # the score ladder is the scheme's; without a scheme it is pcr's (the
  # historical weights), so the bare call keeps its scoring
  scoring <- policy$mode %||% "pcr"
  min_width <- max(nchar(f_primer), nchar(r_primer))
  chroms <- genome_chroms(genome)
  # The whole genome goes in once.  `primer_hit_table()` already loops over
  # `genome_chroms()`, so calling it inside the per-chromosome loop below
  # rescanned every chromosome n_chrom times AND paired hits from different
  # chromosomes into phantom products (each scan iteration relabelled the same
  # pool with its own chromosome).  `dev/specificity-oracle.R` enumerates the
  # answer and caught it: one real product came back as twelve.
  fh_all <- primer_hit_table(f_primer, genome, radius, scoring, side = "f")
  rh_all <- primer_hit_table(rev, genome, radius, scoring, side = "r")
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
  out <- purrr::list_rbind(res)
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
  # Exact mode must not silently hide a one-mismatch primer product.  Scan one
  # mismatch and warn; callers who want those products in the table pass
  # max_mismatch = 1 explicitly.  The near scan is unscored (no scheme): every
  # product of it that is not in `out` is a mismatch product.  An IUPAC primer is
  # searched exactly only, so it has no near scan (and must not reach the
  # unambiguous-primer guard of `primer_hit_table()`).
  if (radius == 0L && !grepl("[^ACGT]", f_primer) && !grepl("[^ACGT]", r_primer)) {
    near <- in_silico_pcr(
      f_primer, r_primer, genome,
      max_mismatch = 1L, BPPARAM = BPPARAM
    )
    extra <- dplyr::anti_join(
      near, out,
      by = c("chrom", "start", "end")
    )
    if (nrow(extra)) {
      msg <- gettext(paste0(
        "%d near-match in-silico PCR product(s) with <= 1 primer mismatch ",
        "were found; exact mode ignores them. Set max_mismatch = 1 to include them."
      ))
      cli::cli_warn(sprintf(msg, nrow(extra)))
    }
  }
  out
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
#' `specificity_verdict()` applies it once to the whole guide.
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
  ok
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
  anchor <- as_anchor(assay$template$snp_i, ref_len, aa$strand)
  seed <- c(
    AS_ref = make_as_primer(assay$template$seq, anchor, aa$as_len, assay$snp$ref, strand = aa$strand),
    AS_alt = make_as_primer(assay$template$seq, anchor, aa$as_len, assay$snp$alt, strand = aa$strand)
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
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' assay_amplicons(assay, genome)
#' }
#' @family KASP specificity
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
assay_pcr_products <- function(assay, genome, max_mismatch = 0L, BPPARAM = NULL) {
  purrr::imap(as_pcr_pairs(assay), function(fr, nm) {
    in_silico_pcr(fr[["F"]], fr[["R"]], genome,
      max_mismatch = max_mismatch, BPPARAM = BPPARAM
    )
  })
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

#' BLAST every primer against a local genome database (rBLAST, optional)
#' @param primers named character vector of oligos.
#' @param db_path BLAST database path.
#' @return a tibble of BLAST hits, or `NULL` when the database is missing or the
#'   query fails.
#' @examples
#' \dontrun{
#' blast_primers(c(AS = "GAAGGTGACCAAGTTCATGCT"), "tair10")
#' }
#' @family KASP specificity
#' @export
blast_primers <- function(primers, db_path) {
  if (is.null(db_path) || !fs::file_exists(db_path)) {
    return(NULL)
  }
  rlang::check_installed("rBLAST")
  query <- Biostrings::DNAStringSet(rlang::set_names(primers, names(primers)))
  out <- tryCatch(
    {
      bl <- rBLAST::blast(db_path, type = "blastn")
      stats::predict(bl,
        newdata = query,
        BLAST_args = c("-task", "blastn-short", "-word_size", "7", "-dust", "no")
      )
    },
    error = function(e) {
      cli::cli_warn(gettext("BLAST failed: {conditionMessage(e)}"))
      NULL
    }
  )
  if (is.null(out) || nrow(out) == 0L) {
    return(out)
  }
  tibble::as_tibble(out)
}

#' Build (once) a BLAST nucleotide database from the TAIR10 genome
#' @param genome DNAStringSet.
#' @param path FASTA path to write (and index).
#' @return the `path` to the FASTA that was written and indexed.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' build_blast_db(genome)
#' }
#' @family KASP specificity
#' @export
build_blast_db <- function(genome, path = fs::path(marker_cache_dir(), "tair10.fa")) {
  fs::dir_create(fs::path_dir(path))
  if (!fs::file_exists(path)) Biostrings::writeXStringSet(genome, path)
  rBLAST::makeblastdb(path, dbtype = "nucl", verbose = FALSE)
  path
}

# =============================================================================
#  7. Scheme-specific mismatch tolerance (heuristic, position-aware)
# =============================================================================

# The tolerance table.  One row per chemistry; `specificity_policy()` returns
# the row as an object.  Every number is a heuristic (calibrated = FALSE) taken
# from the thermodynamics the markers rely on, not fitted to data:
#
#  * pcr/KASP: polymerase extension starts at the 3' terminus, so a 3'-terminal
#    mismatch is not a "mismatch" but loss of extension.  A mismatch in the last
#    ~5 nt is destabilising enough to swallow ordinary allele discrimination, so
#    it is counted separately and capped at 0.  One mismatch elsewhere (the
#    deliberate-mismatch and 5'/middle positions) is tolerated.
#    KASP allele-specific PCR relies on the AS primer's 3' base, so an
#    exact-hit-only screen and a 15-nt 3' seed both miss productive off-targets
#    when the mismatch sits in the middle - hence `max_mismatch = 1` here rather
#    than 0.
#  * rpa: recombinase filaments tolerate mismatches and strand exchange is not
#    proofreading, so several mismatches (up to ~3-5) can still amplify - an
#    exact-hit-only in-silico RPA is wrong.  The 3' end still matters (the
#    strand-displacement polymerase needs a paired terminus) but less strictly
#    than for PCR: one mismatch in the last 5 nt, no terminal mismatch.
#  * crispr: a Cas12a guide's PAM-proximal ~8-12 nt seed is the mismatch-critical
#    region (0-1 seed mismatch is tolerated); Cas13's critical region is central,
#    reported in `critical_from`/`critical_to` for clarity.  A changed PAM is a
#    HARD loss of target (a mismatch in the PAM is not a tolerated primer
#    mismatch and there is no partial credit).  Cas13 has no DNA PAM but a weak
#    PFS preference, carried in `pfs = "weak"` - Cas13 is not PAM-free.
#  * paralog/homeolog: no exact-hit assumption at all - a homeolog may differ by
#    several bases yet still be a source of cross-hybridisation - so it is the
#    most permissive DNA mode, with the 3' end and seed still capped.
SPECIFICITY_POLICIES <- list(
  pcr = list(
    mode = "pcr", max_mismatch = 1L, seed_len = 5L, seed_max = 0L,
    three_prime_max = 0L, requires_pam = FALSE, pam = NA_character_, pfs = NA_character_,
    weight_middle = 1, weight_seed = 3, weight_three_prime = 8, weight_terminal = 12
  ),
  rpa = list(
    mode = "rpa", max_mismatch = 4L, seed_len = 5L, seed_max = 2L,
    three_prime_max = 1L, requires_pam = FALSE, pam = NA_character_, pfs = NA_character_,
    weight_middle = 1, weight_seed = 2, weight_three_prime = 5, weight_terminal = 12
  ),
  crispr = list(
    mode = "crispr", max_mismatch = 5L, seed_len = 8L, seed_max = 1L,
    three_prime_max = 2L, requires_pam = TRUE, pam = "TTTV", pfs = "weak",
    critical_from = 1L, critical_to = 8L,
    weight_middle = 1, weight_seed = 4, weight_three_prime = 8, weight_terminal = 12
  ),
  paralog = list(
    mode = "paralog", max_mismatch = 4L, seed_len = 8L, seed_max = 1L,
    three_prime_max = 1L, requires_pam = FALSE, pam = NA_character_, pfs = NA_character_,
    weight_middle = 1, weight_seed = 2, weight_three_prime = 5, weight_terminal = 12
  )
)

# Chemistry aliases -> a policy name (the scheme a marker object would carry)
SPECIFICITY_SCHEME_ALIASES <- c(
  kap = "pcr", kas = "pcr", kasp = "pcr", tas = "pcr", arms = "pcr",
  caps = "pcr", scar = "pcr", ssr = "pcr", dcap = "pcr", dpcr = "pcr",
  as = "pcr", hrm = "pcr", pcrr = "pcr",
  rpa = "rpa", lamp = "rpa", aslamp = "rpa",
  cas9 = "crispr", cas12 = "crispr", cas12a = "crispr", cas13 = "crispr",
  cas13a = "crispr", cas13b = "crispr", sherlock = "crispr", detectr = "crispr",
  paralog = "paralog", homeolog = "paralog", homeologue = "paralog",
  # a homoeolog is the same concept: a duplicated subgenome copy, not a gene
  # family member.  Spellings follow the literature (wheat papers use
  # "homoeolog"); the tables keep the shorter `homeolog` column name.
  homoeolog = "paralog", homoeologue = "paralog",
  homolog = "paralog", homologe = "paralog"
)

#' Mismatch tolerance for a specificity scheme
#'
#' Encodes the *known biology* each chemistry imposes on a target site: how many
#' mismatches can still amplify, where they may sit, and whether a licensing
#' element (a CRISPR PAM) must be intact.  The tolerance numbers are
#' **heuristics** (`calibrated = FALSE`): they follow the qualitative rule the
#' literature states, they are not fitted to an assay dataset.
#'
#' | scheme | max_mismatch | seed_len | seed_max | three_prime_max | requires_pam |
#' |---|---|---|---|---|---|
#' | `pcr` (KASP/ARMS/CAPS/SSR/dPCR) | 1 | 5 | 0 | 0 | no |
#' | `rpa` (RPA/LAMP) | 4 | 5 | 2 | 1 | no |
#' | `crispr` (Cas9/Cas12/Cas13) | 5 | 8 | 1 | 2 | yes (Cas12a `TTTV`) |
#' | `paralog` (homeolog/homolog) | 4 | 8 | 1 | 1 | no |
#'
#' Position ladder (per-position weight of a mismatch, 5'->3' of the primer):
#'
#' | scheme | `weight_middle` | `weight_seed` | `weight_three_prime` | `weight_terminal` |
#' |---|---|---|---|---|
#' | `pcr` | 1 | 3 | 8 | 12 |
#' | `rpa` | 1 | 2 | 5 | 12 |
#' | `crispr` | 1 | 4 | 8 | 12 |
#' | `paralog` | 1 | 2 | 5 | 12 |
#'
#' `weight_seed` applies to the 3'-terminal `seed_len` nt, `weight_three_prime`
#' to the last `three_prime_window` (5) nt and `weight_terminal` to the last 3 nt
#' - the 3'-terminal mismatch that stops extension.  RPA/LAMP tolerate several
#' mismatches, so a middle mismatch costs less there than for a proofreading
#' polymerase; the terminal weight is the same everywhere, because no chemistry
#' extends from a mismatched 3' end.  The pcr row is the historical scoring, so
#' `pcr`/`kasp` callers' scores are unchanged.
#'
#' Biological basis:
#' * **PCR/KASP** - extension begins at the 3' terminus: a 3'-terminal mismatch
#'   is a hard loss (never a tolerated mismatch), 0 mismatches are allowed in the
#'   last ~5 nt, and 0-1 outside them.  An exact-hit-only screen is not
#'   sufficient, because a middle mismatch can still amplify.
#' * **RPA** - the recombinase filament and strand-displacement polymerase are
#'   mismatch-tolerant: several mismatches (up to ~3-5) can still amplify, so
#'   exact-hit-only is wrong.  The 3' end still matters, so at most one mismatch
#'   in the last 5 nt and none terminal.
#' * **CRISPR** - the Cas12a seed is PAM-proximal (~8-12 nt from the PAM): 0-1
#'   seed mismatch, and a changed PAM is a hard loss of target, not a tolerated
#'   mismatch.  Cas13 has no DNA PAM, but its PFS is a weak preference, carried
#'   in `pfs = "weak"` - Cas13 is never PAM-free.
#' * **paralog/homeolog** - mismatch-aware by construction: a homeolog differs
#'   from the assayed copy by several bases and still cross-hybridises, so an
#'   exact-hit assumption understates off-targets.
#'
#' @param scheme a scheme name (`"pcr"`, `"rpa"`, `"crispr"`, `"paralog"`) or a
#'   chemistry alias (`"kasp"`, `"lamp"`, `"cas12a"`, `"homeolog"`, ...).
#' @param length optional target-site length; when given, the returned object
#'   carries a named `weight` vector with the per-position ladder (5'->3').
#' @return a `specificity_policy` list: `mode`, `max_mismatch`, `seed_len`,
#'   `seed_max`, `three_prime_max`, `three_prime_window`, `requires_pam`, `pam`,
#'   `pfs`, `critical_from`/`critical_to` (CRISPR), the position ladder
#'   (`weight_middle`, `weight_seed`, `weight_three_prime`, `weight_terminal`),
#'   `priors` (the tolerance row as a one-row tibble) and `calibrated` (`FALSE`).
#' @examples
#' specificity_policy("kasp")$max_mismatch
#' specificity_policy("rpa")$max_mismatch
#' specificity_policy("cas12a")$requires_pam
#' specificity_policy("homeolog")$mode
#' @family KASP specificity
#' @export
specificity_policy <- function(scheme = "pcr", length = NULL) {
  if (is.null(scheme)) {
    scheme <- "pcr"
  }
  if (!is.character(scheme) || length(scheme) != 1L) {
    cli::cli_abort(gettext("`scheme` must be a single scheme name."))
  }
  key <- gsub("[-_ ]", "", tolower(scheme))
  alias <- match(key, names(SPECIFICITY_SCHEME_ALIASES))
  key <- if (!is.na(alias)) {
    unname(SPECIFICITY_SCHEME_ALIASES[[alias]])
  } else if (key %in% names(SPECIFICITY_POLICIES)) {
    key
  } else {
    cli::cli_abort(c(
      gettext("Unknown specificity scheme {.val {scheme}}."),
      i = gettext(paste(
        "Use one of {.val {names(SPECIFICITY_POLICIES)}} or an alias",
        "({.val {names(SPECIFICITY_SCHEME_ALIASES)}})."
      ))
    ))
  }
  pol <- SPECIFICITY_POLICIES[[key]]
  pol$scheme <- scheme
  pol$three_prime_window <- THREE_PRIME_WINDOW
  pol$calibrated <- FALSE
  pol$priors <- tibble::tibble(
    mode = pol$mode, max_mismatch = pol$max_mismatch, seed_len = pol$seed_len,
    seed_max = pol$seed_max, three_prime_max = pol$three_prime_max,
    requires_pam = pol$requires_pam, calibrated = FALSE
  )
  if (!is.null(length)) {
    pol$weight <- mismatch_weights(as.integer(length), pol$mode, pol$seed_len)
  }
  structure(pol, class = c("specificity_policy", "list"))
}

#' @export
print.specificity_policy <- function(x, ...) {
  cat("<specificity_policy:", x$mode, "(heuristic, calibrated =", x$calibrated, ")>\n")
  cat("  max_mismatch", x$max_mismatch, "| seed", x$seed_max, "of", x$seed_len,
    "| last", x$three_prime_window, "nt", x$three_prime_max, "\n",
    sep = " "
  )
  cat("  weights: middle", x$weight_middle, "| seed", x$weight_seed,
    "| last", x$three_prime_window, "nt", x$weight_three_prime,
    "| terminal", x$weight_terminal, "\n",
    sep = " "
  )
  if (x$requires_pam) cat("  PAM", x$pam, "required (a changed PAM = loss of target)\n")
  if (!is.na(x$pfs)) cat("  PFS preference:", x$pfs, "\n")
  invisible(x)
}

#' Is there a genome/specificity engine able to answer for this genome?
#'
#' The in-memory scan needs a `DNAStringSet` (or an indexed `marker_genome` read
#' through Rsamtools).  An indexed genome without a BLAST database on this
#' machine has **no** genome-wide engine, and that is reported as missing - the
#' caller must not read it as a pass.  A Bowtie/BWA/crisprVerse backend is not
#' implemented here; when one is needed (CRISPR genome-wide), it is reported
#' missing too, with its capability flag.
#' @param genome a `DNAStringSet`, a `marker_genome`, or `NULL`.
#' @param scheme specificity scheme ([specificity_policy()]).
#' @return `TRUE`/`FALSE` for `engine` (a usable in-memory/indexed scan) and
#'   named capability flags (`exact_scan`, `mismatch_scan`, `blast`,
#'   `crispr_offtarget`, `bowtie`, `bwa`, `crisprverse`).
#' @keywords internal
#' @noRd
specificity_backend <- function(genome = NULL, scheme = "pcr") {
  pol <- specificity_policy(scheme)
  in_memory <- inherits(genome, "DNAStringSet")
  blast_ok <- has_blast() && isTRUE(requireNamespace("rBLAST", quietly = TRUE))
  if (inherits(genome, "marker_genome")) {
    in_memory <- !is.null(genome$blast_db) && blast_ok
  }
  caps <- c(
    exact_scan = in_memory,
    mismatch_scan = in_memory,
    blast = blast_ok,
    crispr_offtarget = FALSE,
    bowtie = nzchar(Sys.which("bowtie")),
    bwa = nzchar(Sys.which("bwa")),
    crisprverse = isTRUE(requireNamespace("crisprVerse", quietly = TRUE))
  )
  list(
    engine = in_memory,
    capabilities = caps,
    scheme = pol$mode,
    reason = if (in_memory) {
      NA_character_
    } else if (is.null(genome)) {
      gettext("No genome was supplied.")
    } else {
      gettext(paste(
        "No genome-wide engine for this genome: it is not an in-memory",
        "DNAStringSet and it has no usable BLAST database."
      ))
    }
  )
}

# Six-frame translation need not be rigorous here: a protein PAM is a 3-codon
# string.  (The guide-level PAM check uses `pam_state()` below.)

#' Is a CRISPR PAM/PFS intact?
#'
#' `"match"` when a PAM matching the scheme's pattern (case-insensitive, IUPAC
#' `V`/`D` expanded) was observed; `"pam_loss"` when it does not - a hard loss of
#' target, never a tolerated primer mismatch; `"absent"` when no PAM was
#' supplied.  `"PFS"` (a Cas13 guide) is checked only for its weak documented
#' rule: a `G` immediately 3' of the protospacer is disfavoured, any other base
#' is `"match"`.  That weak PFS is a preference, not a PAM, and is never
#' reported as "PAM-free".
#' @param observed the PAM/PFS sequence read from the target, or `"PFS"`.
#' @param pattern the licensed PAM for the effector (default Cas12a's `TTTV`).
#' @keywords internal
#' @noRd
pam_state <- function(observed, pattern = "TTTV") {
  if (is.null(observed) || is.na(observed) || !nzchar(observed)) {
    return("absent")
  }
  obs <- toupper(observed)
  if (identical(obs, "PFS")) {
    return("match")
  }
  if (nchar(obs) != nchar(pattern)) {
    return("pam_loss")
  }
  ok <- mapply(function(code, base) iupac_covers(code, base),
    strsplit(pattern, "", fixed = TRUE)[[1L]],
    strsplit(obs, "", fixed = TRUE)[[1L]],
    USE.NAMES = FALSE
  )
  if (all(ok)) "match" else "pam_loss"
}

#' Grade in-silico PCR products against a scheme's mismatch policy
#'
#' The structured answer to "is this specific enough?": `pass` when at least one
#' product is tolerated *and* it is the only tolerated product on its chromosome;
#' `fail` when every product breaks the policy (including a changed PAM), or
#' more than one tolerated product exists; `unknown` when there is no genome-wide
#' engine - a missing backend is never a silent pass.  The result carries the
#' products, the tolerated flag, the best product and the policy it was graded
#' against, so the evidence travels with the verdict.
#' @param hits in-silico PCR products ([in_silico_pcr()], as a tibble: one row
#'   per product, `chrom`/`start`/`end`/`width` and the mismatch evidence;
#'   extra columns - e.g. a `guide` id from a multi-guide table - are ignored.
#'   An optional `scheme`/`tolerated` pair is recomputed here, so a table from
#'   `in_silico_pcr(scheme = )` is graded exactly as a bare one.
#' @param scheme specificity scheme ([specificity_policy()]).
#' @param genome optional genome; when supplied, a missing engine turns the
#'   answer into `unknown` (see [specificity_backend()]).
#' @param pam optional PAM/PFS observation for `"crispr"` (see [pam_state()]).
#' @return a `specificity_verdict` list: `verdict` (`"pass"`/`"fail"`/
#'   `"unknown"`), `reason` (an English sentence), `scheme`, `policy`,
#'   `n_products`, `n_tolerated`, `n_mismatch` (best product), `best` (one-row
#'   product tibble), `products` (the graded tibble), `pam` (state) and
#'   `backend` (the capability list, or `NA` when no genome was supplied).
#' @examples
#' genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
#'   strrep("A", 20L), "ACGTACGTACGTACGTACGT", strrep("T", 20L),
#'   "ACGTACGTACGTACGTACGT", strrep("A", 20L)
#' )))
#' hits <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
#' # `unknown` when the genome is missing: a missing engine is never a pass
#' specificity_verdict(hits, "kasp", genome = NULL)$verdict
#' @family KASP specificity
#' @export
specificity_verdict <- function(hits, scheme = "pcr", genome = NULL, pam = NULL) {
  pol <- specificity_policy(scheme)
  backend <- if (is.null(genome)) NA else specificity_backend(genome, scheme)
  products <- if (is.null(hits)) {
    empty_pcr_products()
  } else {
    tibble::as_tibble(hits)[intersect(names(empty_pcr_products()), names(hits))]
  }
  pam_given <- !is.null(pam)
  if (!all(c("n_mismatch", "seed_mismatch", "three_prime_mismatch") %in% names(products))) {
    # a legacy exact table: it carries no evidence, so it cannot be graded
    products$n_mismatch <- 0L
    products$seed_mismatch <- 0L
    products$three_prime_mismatch <- 0L
    products$three_prime_terminal <- FALSE
  }
  if (!"penalty" %in% names(products)) products$penalty <- 0
  # A CRISPR PAM is licensing, not a mismatch: when it was observed, a changed
  # PAM removes the tolerated status of every product, even a perfect one.
  pam_val <- if (pam_given && identical(pol$mode, "crispr")) {
    pam_state(pam, pol$pam %||% "TTTV")
  } else {
    NA_character_
  }
  products$tolerated <- pcr_products_tolerated(products, pol)
  products$violations <- mapply(pcr_violations,
    n_mismatch = products$n_mismatch, seed_mismatch = products$seed_mismatch,
    three_prime_mismatch = products$three_prime_mismatch,
    three_prime_terminal = products$three_prime_terminal,
    MoreArgs = list(policy = pol), USE.NAMES = FALSE
  )
  if (identical(pam_val, "pam_loss")) {
    products$tolerated <- rep(FALSE, nrow(products))
    products$violations <- lapply(products$violations, function(v) union(v, "pam"))
  }
  # The best evidence: a tolerated product first, then the smallest penalty.  At
  # equal penalty the complete (widest) product wins over an end-padded one, so
  # the scan radius cannot manufacture the headline number.
  best <- if (nrow(products)) {
    products[order(
      !products$tolerated, products$penalty,
      -products$width, products$start
    ), ][1L, ]
  } else {
    NULL
  }
  n_tol <- sum(products$tolerated)
  # a backend that is not a list at all (no genome given) carries no engine
  if (!is.list(backend) || !isTRUE(backend$engine)) {
    out <- list(
      verdict = "unknown",
      reason = if (is.null(genome)) {
        gettext("No genome-wide engine was consulted; a missing engine is never a pass.")
      } else if (nzchar(backend$reason)) {
        backend$reason
      } else {
        gettext("No genome-wide engine for this genome; a missing engine is never a pass.")
      },
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = n_tol, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else if (n_tol == 0L) {
    # a CRISPR guide whose PAM changed has lost its target: there is no
    # tolerated product to be specific about, so this is not a specificity pass
    reason <- if (!nrow(products)) {
      gettext("No in-silico PCR product at all: the primers do not amplify.")
    } else if (identical(pam_val, "pam_loss")) {
      gettext("The guide's PAM changed: a hard loss of target, not a tolerated mismatch.")
    } else {
      gettext("Every product breaks the scheme's mismatch policy.")
    }
    out <- list(
      verdict = "fail",
      reason = reason,
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = 0L, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else if (n_tol > 1L) {
    out <- list(
      verdict = "fail",
      reason = gettextf(
        "%d tolerated in-silico PCR products: not specific.", n_tol
      ),
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = n_tol, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else {
    out <- list(
      verdict = "pass",
      reason = gettext("Exactly one tolerated in-silico PCR product."),
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = 1L, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  }
  structure(out, class = c("specificity_verdict", "list"))
}

#' Which caps of a policy a product breaks
#' @keywords internal
#' @noRd
pcr_violations <- function(n_mismatch, seed_mismatch, three_prime_mismatch,
                           three_prime_terminal, policy) {
  v <- character()
  if (n_mismatch > policy$max_mismatch) v <- c(v, "total")
  if (seed_mismatch > policy$seed_max) v <- c(v, "seed")
  if (three_prime_mismatch > policy$three_prime_max) v <- c(v, "three_prime")
  if (isTRUE(three_prime_terminal)) v <- c(v, "terminal")
  v
}

#' @export
print.specificity_verdict <- function(x, ...) {
  cat("<specificity_verdict:", x$verdict, ">\n")
  cat("  scheme:", x$scheme, "| products:", x$n_products,
    "| tolerated:", x$n_tolerated, "\n",
    sep = " "
  )
  cat("  reason:", x$reason, "\n")
  invisible(x)
}
