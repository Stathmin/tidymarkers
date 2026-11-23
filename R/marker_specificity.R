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
#
#  Evidence is graded, never assumed.  `specificity_seed()` is a cheap
#  *prefilter* (3'-seed counts, fixed false-negative rate) for use while a
#  candidate is built; the verdict `specificity_verdict()` gives prefers
#  pair-level in-silico PCR with product-size bounds, then a BLAST database or
#  FM index (the only tiers that can fail a primer genome-wide), and is
#  `"unknown"` when none of them was consulted.  The result names the tier.
# =============================================================================

# The in-memory count is O(genome) per primer, so past this size the "prefilter"
# is the cost.  Basis: 60 seeds over TAIR10 Chr1 (30.4 Mb) measured 0.58 s.
SEED_SCAN_WARN_BP <- 1e9
SEED_SCAN_S_PER_PRIMER_PER_MB <- 3.2e-4

#' Warn once per session when the in-memory scan pays O(genome) per primer
#'
#' The seed tier holds no index: every call re-reads the genome once per primer.
#' That is invisible while the genome fits in a few hundred Mb and dominates the
#' design at wheat scale, so say the price once instead of 15 minutes into a
#' locus.
#' @param bp genome size in bases.
#' @param n_primers number of seeds scanned.
#' @keywords internal
#' @noRd
seed_scan_cost_warn <- function(bp, n_primers) {
  if (!is.finite(bp) || bp < SEED_SCAN_WARN_BP) {
    return(invisible(NULL))
  }
  session_once(paste0("seed-scan-", bp), {
    size <- format(bp, big.mark = ",", scientific = FALSE)
    eta <- format(
      round(bp / 1e6 * n_primers * SEED_SCAN_S_PER_PRIMER_PER_MB),
      big.mark = ",", scientific = FALSE
    )
    cli::cli_warn(c(
      gettextf("This seed scan reads %s bp in memory, O(genome) per primer.", size),
      i = gettextf("~%s s for %d primer(s); only an index is flat in genome size.", eta, n_primers),
      i = gettext("Build marker_bowtie_index(genome) and pass bowtie_index= (~1.34 bytes/bp, ~0.05 s per primer set).")
    ))
  })
  invisible(NULL)
}

#' Count genome-wide near-matches of each primer's 3'-terminal k-mer (prefilter)
#'
#' This is a **prefilter, not a genome-wide guarantee**.  It counts only the
#' 3'-terminal `k`-mer, so it has a fixed false-negative rate: with `k = 15` a
#' single mismatch at offset -5 (i.e. inside the last `k` nt), or a second
#' mismatch in the seed, makes a real off-target site invisible to the count.
#' A zero count therefore means "no site shares this seed", not "no site binds".
#' Genome-wide answers come from pair-level [in_silico_pcr()] with product-size
#' bounds, or from an indexed scan ([specificity_blast()], [specificity_fm()]);
#' [specificity_verdict()] grades against those and reports `"unknown"` when
#' only the seed tier was available.
#'
#' For the exact case (`max_mismatch = 0`) counts come from `Biostrings::PDict`
#' (`vcountPDict`) - the precomputed-dictionary fast path, one call for all
#' primers.  With mismatches Biostrings has no counting index (`PDict()` builds an
#' `MTB_PDict`, but `vcountPDict()` does not support it), so the scan is done with
#' `vcountPattern` **parallelised over the genome's sequences** (and, per
#' sequence, once for all seeds), which is how this scales to large genomes.
#' Results are summed over both strands and over all genome sequences.
#'
#' A `marker_genome` is dispatched to its index instead of being refused: the
#' seed itself is the query, so a whole-primer scan of the seed answers exactly
#' this question.  With a `bowtie_index` the count is a seed within
#' `max_mismatch` substitutions; with only a `blast_db` it is an **exact** seed
#' count (that path keeps full-length exact hits only), i.e. a stricter screen
#' than the in-memory one.
#' @param primers named character vector of oligos.
#' @param genome a `DNAStringSet`, or a `marker_genome` carrying a
#'   `bowtie_index` or a `blast_db`.
#' @param k 3'-seed length.
#' @param max_mismatch mismatches allowed in the seed (in-memory and bowtie
#'   paths; the BLAST path counts exact matches only).
#' @param BPPARAM optional `BiocParallel` backend; the mismatch scan is split by
#'   genome sequence, so `MulticoreParam()` uses one worker per chromosome.
#' @return a tibble `primer`, `seed`, `hits` (both-strand counts of the seed,
#'   over every genome sequence - a prefilter, see Details), `engine` (the
#'   adapter that answered: `"Biostrings"` in memory, else the indexed engine
#'   the call delegates to) and an `engine_stamp` attribute (see
#'   [marker_engine_stamp()]).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' specificity_seed(c(AS = "GAAGGTGACCAAGTTCATGCT"), genome)
#' }
#' @family specificity
#' @export
specificity_seed <- function(primers, genome, k = SEED_K, max_mismatch = SEED_MISMATCH,
                             BPPARAM = NULL) {
  seeds <- unname(substr(primers, nchar(primers) - k + 1L, nchar(primers)))
  if (inherits(genome, "marker_genome")) {
    # the seed is the whole query, so an indexed whole-primer scan answers this
    # question; refusing here silently dropped the AS 3'-seed tie-break
    query <- rlang::set_names(seeds, names(primers))
    if (!is.null(genome$bowtie_index)) {
      return(specificity_bowtie(query, genome$bowtie_index,
        max_mismatches = max_mismatch,
        nthreads = if (is.null(BPPARAM)) 4L else 1L
      ))
    }
    if (!is.null(genome$blast_db)) {
      return(specificity_blast(query, genome$blast_db,
        k = NA_integer_, BPPARAM = BPPARAM
      ))
    }
    cli::cli_abort(c(
      gettext("This indexed genome has neither a bowtie index nor a BLAST database."),
      i = gettext("Build marker_bowtie_index(genome) (bowtie_index=) or marker_blast_index(genome) (blast_db=).")
    ))
  }
  if (inherits(genome, "DNAStringSet")) {
    seed_scan_cost_warn(sum(Biostrings::width(genome)), length(seeds))
  }
  if (max_mismatch == 0L) {
    hits <- exact_seed_counts(seeds, genome)
  } else if (max_mismatch == 1L && all(grepl("^[ACGT]+$", seeds))) {
    # A site within one substitution of a seed is exactly one of the seed's 46
    # Hamming neighbours, so a single exact dictionary pass over the 46 * n
    # neighbours answers every seed.  Measured on TAIR10 Chr1, 60 seeds:
    # 4.40 s for one `countPattern(max.mismatch = 1)` pass per seed vs 0.58 s
    # for the batched exact pass, identical counts.  Degenerate or longer-radius
    # input takes the per-seed loop below.
    hits <- near_seed_counts(seeds, genome)
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
  marker_with_engine(
    tibble::tibble(primer = names(primers), seed = seeds, hits = as.integer(hits)),
    "Biostrings"
  )
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

# Bases of one chromosome scanned per grouped `matchPattern()` call.  A wheat
# chromosome is 700 Mb and a batch holds every primer of a candidate pool, so an
# unbounded pass would copy the whole chromosome once per width group; windowing
# it bounds the per-copy memory at this many bases.  The genome passes become
# `n_chrom * n_widths * ceil(chrom_len / HIT_BATCH_MAX)` rather than one per
# primer.  (Declared here rather than next to the hit cache in `marker_core.R`
# because the grouped search is the only reader.)
HIT_BATCH_MAX <- 1e7L

# The seven profile fields a single-primer hit table carries (a locus, not a
# product): one `mismatch_positions` and a `terminal_mismatch`, unlike the
# per-primer f/r split of PCR_MISMATCH_COLUMNS.
SINGLE_PRIMER_PROFILE_COLUMNS <- c(
  "n_mismatch", "mismatch_positions", "three_prime_mismatch", "seed_mismatch",
  "terminal_mismatch", "three_prime_terminal", "penalty"
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

#' Every string within one substitution of a seed
#'
#' A genomic site matches a seed with at most one mismatch exactly when the site
#' equals one of these strings, so counting the Hamming ball exactly is the same
#' question as counting with `max.mismatch = 1`.
#' @param seed a plain-ACGT string.
#' @return a character vector: the seed itself, then its 3 * nchar substitutions.
#' @keywords internal
#' @noRd
seed_ball <- function(seed) {
  bases <- c("A", "C", "G", "T")
  sp <- strsplit(seed, "")[[1L]]
  out <- seed
  for (i in seq_along(sp)) {
    for (b in bases) {
      if (b != sp[[i]]) {
        v <- sp
        v[[i]] <- b
        out <- c(out, paste(v, collapse = ""))
      }
    }
  }
  out
}

#' Both-strand seed counts allowing one mismatch, by exact dictionary counting
#'
#' `PDict`/`vcountPDict()` cannot count with a mismatch radius, so the radius is
#' moved into the pattern set: the exact counts of a seed's 46 Hamming neighbours
#' sum to its one-mismatch count, and one pass over all neighbours of all seeds
#' replaces one pass per seed.  Same counts, one genome scan.
#' @param seeds plain-ACGT seeds (all the same width).
#' @param genome a `DNAStringSet`.
#' @return integer vector, one count per seed.
#' @keywords internal
#' @noRd
near_seed_counts <- function(seeds, genome) {
  ball <- lapply(seeds, seed_ball)
  n <- lengths(ball)
  pats <- unlist(ball, use.names = FALSE)
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(pats)))
  dict <- unique(c(pats, rc))
  cnt <- as.integer(rowSums(Biostrings::vcountPDict(Biostrings::PDict(dict), genome)))
  names(cnt) <- dict
  as.integer(rowsum(cnt[pats] + cnt[rc], rep(seq_along(seeds), n)))
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
#' reported 1-based from the 5' end of `pattern`, i.e. the plus strand that was
#' searched.  A reverse primer is searched as its reverse complement, so its
#' caller re-orients the pattern and footprint into the primer's own 5'->3'
#' frame before profiling (see `hit_profile()`), rather than mirroring positions
#' afterwards.
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
#' answer is memoised per process - the N count included, because that O(length)
#' pass is the cost and it too is a property of the sequence, not of the pattern.
#' The key is the sequence's content signature (`digest::digest` on the bytes is
#' not free either, so the key is the string of the length plus the first, middle
#' and last 32 bases) - a *changed* sequence must not reuse a stale block list.
#' The cache lives in the package namespace and is bounded.
#'
#' @param subject the genomic sequence as an `XString`.
#' @return a list of integer vectors `start`, `end` (1-based, inclusive).
#' @keywords internal
#' @noRd
unambiguous_blocks <- function(subject) {
  len <- length(subject)
  key <- subject_signature(subject)
  cached <- get0(key, envir = .marker_state$blocks, inherits = FALSE)
  if (!is.null(cached)) {
    return(cached)
  }
  out <- if (Biostrings::letterFrequency(subject, "N") == 0L) {
    list(start = 1L, end = len)
  } else {
    ap <- which(as.integer(subject) == 15L)
    gap <- diff(ap) > 1L
    a_start <- ap[c(TRUE, gap)]
    a_end <- ap[c(gap, TRUE)]
    # the clean blocks are the complement of the ambiguous runs
    b_start <- c(1L, a_end + 1L)
    b_end <- c(a_start - 1L, len)
    ok <- b_end >= b_start
    list(start = b_start[ok], end = b_end[ok])
  }
  # memory bound: keep only the current genome's block list
  if (length(.marker_state$blocks) > 16L) {
    rm(list = ls(.marker_state$blocks, all.names = TRUE), envir = .marker_state$blocks)
  }
  assign(key, out, envir = .marker_state$blocks)
  out
}

#' Content signature of a genomic sequence
#'
#' The key both genome-derived caches use: length plus the first, middle and
#' last 32 bases.  It deliberately does *not* contain the N count - counting N
#' is the O(length) cost (`letterFrequency` on a 30 Mb chromosome), and it used
#' to be paid on every pattern of every primer pair on that chromosome.
#'
#' @param subject the genomic sequence as an `XString`.
#' @return a single string.
#' @keywords internal
#' @noRd
subject_signature <- function(subject) {
  len <- length(subject)
  mid <- max(1L, min(len, len %/% 2L))
  paste0(len, "|", paste(
    as.character(Biostrings::subseq(subject, 1L, min(len, 32L))),
    as.character(Biostrings::subseq(subject, mid, min(len, mid + 31L))),
    as.character(Biostrings::subseq(subject, max(1L, len - 31L), len)),
    sep = "|"
  ))
}

#' Short content-derived name for a genome's cache files
#'
#' Two genomes never share a derived file name by construction: an in-memory
#' `DNAStringSet` is keyed by `subject_signature()` of every sequence, an indexed
#' `marker_genome` by its sequence names and lengths.
#'
#' A `BSgenome` never reaches this function: `genome_fasta_source()` materialises
#' it as a FASTA-backed `marker_genome` first, and that FASTA is keyed by its DNA
#' (`bsgenome_fasta_key()`) rather than by names and lengths, which cannot
#' separate two same-sized genomes.
#'
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @return a 12-character hex string.
#' @keywords internal
#' @noRd
genome_signature <- function(genome) {
  sig <- if (inherits(genome, "marker_genome")) {
    paste0(names(genome$seqnames), ":", genome$lengths, collapse = "|")
  } else {
    paste(vapply(genome, subject_signature, character(1)), collapse = "|")
  }
  substr(rlang::hash(sig), 1L, 12L)
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
#' match (`mismatch_profile()` already reports it as a non-mismatch *within* a
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
  # One design run asks for the same (pattern, chromosome, radius) several
  # times - once per candidate pairing it is paired in, once more by the
  # `max_mismatch = 1` near-match pass of `in_silico_pcr()`, once for the final
  # assay - and a whole-chromosome scan is the expensive part (82 s of a 97 s
  # 6-SNP design).  The answer is a pure function of the three arguments, so it
  # is cached per subject content signature.  A pathological result (a long run
  # of hits, e.g. a primer matching an N-free repeat) is not cached: the entry
  # would cost more than the rescan.  The signature is the one
  # `unambiguous_blocks()` uses (length + first/middle/last 32 bases), so two
  # subjects that differ only in the middle collide - the same assumption that
  # cache already makes, and one genome per session is the intended use.
  key <- paste(max_mismatch, subject_signature(subject), pattern, sep = "\r")
  cached <- get0(key, envir = .marker_state$hits, inherits = FALSE)
  if (!is.null(cached)) {
    return(cached)
  }
  out <- pattern_hits_scan(pattern, subject, max_mismatch)
  if (length(out) <= HIT_CACHE_MAX) {
    keys <- ls(.marker_state$hits, all.names = TRUE)
    # Evict half the cache, not all of it: flushing it whole would drop the
    # answers the current design is still re-asking for.  `ls()` sorts rather
    # than listing insertion order (an environment has none to expose), so which
    # half goes is arbitrary; the point is that the cache survives the eviction.
    if (length(keys) >= HIT_CACHE_SLOTS) {
      rm(list = keys[seq_len(length(keys) %/% 2L)], envir = .marker_state$hits)
    }
    assign(key, out, envir = .marker_state$hits)
  }
  out
}

#' One pattern search over one sequence (the uncached body of `exact_pattern_hits()`)
#' @param pattern the plus-strand pattern to search (IUPAC-aware).
#' @param subject the genomic sequence as an `XString`.
#' @param max_mismatch search radius.
#' @return an `IRanges` of hit starts/ends in `subject` coordinates.
#' @keywords internal
#' @noRd
pattern_hits_scan <- function(pattern, subject, max_mismatch) {
  width <- nchar(pattern)
  if (!grepl("[^ACGT]", pattern)) {
    # A plain-ACGT pattern needs no block split: `fixed = TRUE` matches the
    # subject byte-for-byte, so an `N` (or any other ambiguity code) in the
    # subject can never satisfy a pattern base - the same invariant the block
    # path buys by searching each unambiguous block on its own.  Measured on
    # TAIR10: identical starts, 3.20 s -> 0.27 s per pattern at `max_mismatch =
    # 0`.  With a mismatch radius a hit *can* straddle an ambiguity code (it is
    # then one mismatch), so those hits are dropped explicitly - that is the
    # one thing the block split used to enforce.
    hits <- Biostrings::matchPattern(pattern, subject,
      max.mismatch = max_mismatch, fixed = TRUE
    )
    if (length(hits) && max_mismatch > 0L) {
      bl <- unambiguous_blocks(subject)
      s <- IRanges::start(hits)
      i <- findInterval(s, bl$start)
      inside <- i >= 1L & s + width - 1L <= bl$end[pmax(i, 1L)]
      hits <- hits[inside]
    }
    return(hits)
  }
  bl <- unambiguous_blocks(subject)
  if (length(bl$start) == 1L && bl$start[[1L]] == 1L && bl$end[[1L]] == length(subject)) {
    return(Biostrings::matchPattern(pattern, subject, max.mismatch = max_mismatch, fixed = FALSE))
  }
  len <- bl$end - bl$start + 1L
  keep <- len >= width
  if (!any(keep)) {
    return(IRanges::IRanges())
  }
  starts <- lapply(which(keep), function(i) {
    s <- bl$start[[i]]
    sub <- Biostrings::subseq(subject, s, bl$end[[i]])
    h <- Biostrings::matchPattern(pattern, sub, max.mismatch = max_mismatch, fixed = FALSE)
    if (!length(h)) {
      return(integer())
    }
    IRanges::start(h) + s - 1L
  })
  starts <- unlist(starts, use.names = FALSE)
  # Build one IRanges from the starts.  `do.call(c, ...)` over the per-block
  # ranges used to work only when the *first* block had a hit: a leading `NULL`
  # made `c()` fall back to list semantics, and `start()` of that list returned
  # the index of each block instead of the genomic coordinate.
  if (!length(starts)) {
    return(IRanges::IRanges())
  }
  IRanges::IRanges(starts, width = width)
}

#' Profile one primer's hit footprints, oriented into the primer's own frame
#'
#' `mismatch_profile()` scores the searched (plus-strand) pattern, so for a
#' minus-strand search the pattern is the primer's reverse complement and the
#' primer's 3'-terminal base sits at plus-strand position 1, not at the tail
#' where the 3'-window masks look.  Rather than re-derive the weight ladder, the
#' reverse case is profiled in the primer's own 5'->3' frame by reverse-
#' complementing the pattern and the footprints before delegating, so every
#' field - positions, the 3' windows and the penalty - lands on the primer.
#' @param pattern searched plus-strand pattern (character; for a minus-strand
#'   search this is already the primer's reverse complement).
#' @param footprints character vector of plus-strand genomic slices of exactly
#'   `nchar(pattern)`.
#' @param scheme chemistry.
#' @param seed_len seed length override (`NA` = the policy's).
#' @param reverse is the searched pattern the reverse complement of the primer?
#' @return the seven profile vectors (see `mismatch_profile()`), with positions
#'   in the primer's own 5'->3' frame and `AsIs`-wrapped.
#' @keywords internal
#' @noRd
hit_profile <- function(pattern, footprints, scheme, seed_len = NA_integer_,
                        reverse = FALSE) {
  if (reverse) {
    pattern <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(pattern))
    )
    footprints <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAStringSet(footprints))
    )
  }
  prof <- mismatch_profile(pattern, footprints, scheme, seed_len)
  prof$mismatch_positions <- lapply(prof$mismatch_positions, function(pos) {
    # an AsIs element keeps the column a vector of integer vectors (a bare list
    # of length 1 would be treated as a column and nested)
    structure(as.integer(pos), class = "AsIs")
  })
  prof
}


#' Scan several plus-strand patterns over several chromosomes in one pass each
#'
#' The cost of a specificity scan is the number of passes over the genome, not the
#' number of patterns in a pass: one `matchPattern()` call on a whole 30 Mb
#' chromosome is 0.27 s whether it is handed one 27-mer or forty of them.  A
#' caller that needs the hits of several primers (the two allele-specific
#' pairings of an assay, the `pair_screen` candidates of an assembly) therefore
#' collects every pattern it will ask for, and this answers all of them with one
#' grouped search per chromosome per radius.
#'
#' The result is identical to calling `pattern_hits_scan()` once per pattern:
#' * equal-width patterns are joined into one string and searched in one
#'   `matchPattern()` call, and a hit of the joined string is `nchar(joined)`
#'   wide, so its pattern and its offset are arithmetic, not a second search;
#' * the join is plain-ACGT exactly when every part is, and `fixed` follows that
#'   (the same `fixed = TRUE` fast path `pattern_hits_scan()` takes);
#' * the unambiguous-block split is the same one, with the same per-pattern
#'   width filter, so no hit straddles an ambiguity code.
#'
#' The radius stays base-space.  Moving it into the pattern set (the 46 Hamming
#' neighbours `near_seed_counts()` uses) is 15-19x slower for a primer of this
#' length, and a `PDict` refuses a mixed width - both measured, see ROADMAP.
#'
#' A wheat chromosome is 700 Mb, so a batch is scanned in blocks bounded by
#' `HIT_BATCH_MAX` bases rather than one copy of the whole chromosome per width
#' group.
#' @param patterns plus-strand patterns to search (IUPAC-aware); a repeat is
#'   scanned once.
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param max_mismatch search radius.
#' @return a named list per *unique* pattern, each a list named by chromosome
#'   holding the `IRanges` of hits in that chromosome's coordinates, so a caller
#'   reads `hits[[pattern]][[chrom]]`.
#' @keywords internal
#' @noRd
pattern_hits_batch <- function(patterns, genome, max_mismatch) {
  pats <- unique(patterns)
  out <- stats::setNames(
    lapply(pats, function(p) stats::setNames(vector("list", length(genome_chroms(genome))), genome_chroms(genome))),
    pats
  )
  if (!length(pats)) {
    return(out)
  }
  # one search per (chromosome, width): `matchPattern()` needs equal widths
  widths <- sort(unique(nchar(pats)))
  for (chrom in genome_chroms(genome)) {
    worker_biostrings()
    subject <- genome_seq(genome, chrom)
    bl <- unambiguous_blocks(subject)
    bl_len <- bl$end - bl$start + 1L
    for (w in widths) {
      group <- pats[nchar(pats) == w]
      found <- batch_probe(group, subject, bl, bl_len, w, max_mismatch)
      for (p in group) {
        out[[p]][[chrom]] <- found[[p]]
      }
    }
  }
  out
}

#' Search one width group over one chromosome in one grouped `PDict()` pass
#'
#' The blocks of the subject are walked in windows of `HIT_BATCH_MAX` bases (the
#' window is a whole count of blocks), which is one search per window per group
#' and keeps the memory of a wheat chromosome flat.  An exact plain-ACGT group is
#' searched with one `matchPDict()` over the whole group (which is what makes a
#' batch cost less than a scan per pattern); at a mismatch radius a `PDict()` is
#' no help - its trusted-band tree is rebuilt per block and Biostrings spends
#' most of the search collecting garbage (92.7 % of a 193 s radius-1 pool pass
#' was `gc()` inside `.match.MTB_PDict`) - so every pattern is searched with the
#' `matchPattern()` engine of `pattern_hits_scan()`, whose hit set is identical
#' (measured on TAIR10: same starts, 0.07 s vs 1.2 s per pattern per chromosome).
#' A group holding an IUPAC pattern has no `PDict`, so each pattern is searched
#' on its own with the same block walk.
#' @param group equal-width patterns.
#' @param subject the genomic sequence.
#' @param bl,bl_len the unambiguous blocks (`start`, `end`, length).
#' @param w the pattern width.
#' @param max_mismatch search radius.
#' @return a named list of `IRanges` per pattern, in `subject` coordinates.
#' @keywords internal
#' @noRd
batch_probe <- function(group, subject, bl, bl_len, w, max_mismatch) {
  n <- length(group)
  fixed <- !any(grepl("[^ACGT]", group))
  # `PDict()` only pays on an exact pass; a mismatch radius is searched with
  # `matchPattern()` below (see the function's note).  A pattern shorter than 6
  # nt is searched the same way on either engine.
  pdict <- if (fixed && max_mismatch == 0L) {
    Biostrings::PDict(Biostrings::DNAStringSet(group), max.mismatch = 0L)
  } else {
    NULL
  }
  single_block <- length(bl$start) == 1L && bl$start[[1L]] == 1L &&
    bl$end[[1L]] == length(subject)
  if (!is.null(pdict) && max_mismatch == 0L && single_block) {
    # An exact search of a gapless subject needs no block walk.  The walk exists so
    # a hit cannot straddle an ambiguity code or make the subseq copy of a whole
    # wheat chromosome; with one block holding the whole subject neither applies,
    # and a `PDict` already answers the group, so ONE pass replaces one per
    # `HIT_BATCH_MAX` window (measured on a 30 Mb subject, five primers: 0.24 s
    # for the walk vs 0.03 s for this pass, identical hits).
    found <- lapply(
      Biostrings::matchPDict(pdict, subject, max.mismatch = 0L),
      function(h) if (length(h)) IRanges::IRanges(IRanges::start(h), width = w) else IRanges::IRanges()
    )
    names(found) <- group
    return(found)
  }
  # a block must hold a whole pattern, or no hit can start in it
  keep <- which(bl_len >= w)
  starts <- vector("list", n)
  if (length(keep)) {
    # Every pattern is searched inside ONE unambiguous block.  The blocks are
    # packed into `HIT_BATCH_MAX` windows to bound the memory of a subseq copy, and
    # each block of a window is searched on its own: a search over a whole window
    # would let a hit start in the last `w - 1` bases of a block and run on into
    # the ambiguity code after it, which the per-block search of
    # `pattern_hits_scan()` never reports.  (Measured: 6,046 vs 4,595 hits for an
    # IUPAC pattern on TAIR10 Chr2.)
    grp <- cumsum(c(TRUE, diff(cumsum(bl_len[keep])) > HIT_BATCH_MAX))
    for (win in split(keep, grp[seq_along(keep)])) {
      s <- bl$start[[win[[1L]]]]
      e <- bl$end[[win[[length(win)]]]]
      sub <- if (s == 1L && e == length(subject)) subject else Biostrings::subseq(subject, s, e)
      off <- s - 1L
      for (i in win) {
        bs <- bl$start[[i]]
        be <- bl$end[[i]]
        if (bs == s && be == e) {
          piece <- sub
          rel0 <- 0L
        } else {
          piece <- Biostrings::subseq(sub, bs - off, be - off)
          rel0 <- bs - s
        }
        # one hit set per pattern, in the block's own coordinates
        v <- if (!is.null(pdict)) {
          Biostrings::matchPDict(pdict, piece, max.mismatch = max_mismatch)
        } else {
          lapply(group, function(p) Biostrings::matchPattern(p, piece, max.mismatch = max_mismatch, fixed = fixed))
        }
        for (j in seq_len(n)) {
          if (length(v[[j]])) {
            starts[[j]] <- c(starts[[j]], IRanges::start(v[[j]]) + rel0 + off)
          }
        }
      }
    }
  }
  found <- lapply(seq_len(n), function(i) {
    s <- starts[[i]]
    if (length(s)) IRanges::IRanges(s, width = w) else IRanges::IRanges()
  })
  names(found) <- group
  found
}

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
#' \dontrun{
#' genome <- ath_genome()
#' exact_hit_table(c(AS = "GAAGGTGACCAAGTTCATGCT"), genome)
#' }
#' @family specificity
#' @export
exact_hit_table <- function(sequences, genome, max_mismatch = 0L, BPPARAM = NULL,
                            scheme = "pcr", seed_len = NA_integer_) {
  if (inherits(sequences, "XStringSet")) {
    sequences <- as.character(sequences)
  }
  if (!is.character(sequences) || !length(sequences)) {
    cli::cli_abort(gettext("`sequences` must be a non-empty character vector of oligos."))
  }
  max_mismatch <- as.integer(max_mismatch)
  if (is.na(max_mismatch) || max_mismatch < 0L) {
    cli::cli_abort(gettext("`max_mismatch` must be a non-negative integer."))
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
  required <- c("primer", "chrom", "start", "end", "strand")
  missing_cols <- setdiff(required, names(hits))
  if (length(missing_cols)) {
    cli::cli_abort(gettext("`hits` is missing required columns: {.val {missing_cols}}."))
  }
  if (is.null(primers)) {
    if (!"sequence" %in% names(hits)) {
      cli::cli_abort(gettext(
        "`hits` carries no `sequence` column, so `primers` must be supplied."
      ))
    }
    key <- as.character(hits$primer)
    seq <- as.character(hits$sequence)
    if (anyNA(seq)) {
      cli::cli_abort(gettext("`hits$sequence` must not contain NA."))
    }
    pairs <- unique(data.frame(name = key, seq = seq, stringsAsFactors = FALSE))
    if (anyDuplicated(pairs$name)) {
      cli::cli_abort(gettext(
        "`hits$sequence` gives one primer name more than one oligo; supply `primers` explicitly."
      ))
    }
    primers <- stats::setNames(pairs$seq, pairs$name)
  }
  if (!is.character(primers) || is.null(names(primers)) || any(is.na(names(primers)) | names(primers) == "")) {
    cli::cli_abort(gettext("`primers` must be a named character vector of oligos."))
  }
  primer_name <- as.character(hits$primer)
  strand <- as.character(hits$strand)
  bad_strand <- !is.na(strand) & !(strand %in% c("+", "-"))
  if (any(bad_strand)) {
    cli::cli_abort(gettext(
      "`hits$strand` must be {.val {c('+', '-')}}; got {.val {unique(strand[bad_strand])}}."
    ))
  }
  start <- hits$start
  if (!is.numeric(start) || anyNA(start) || any(start < 1) || any(start != floor(start))) {
    cli::cli_abort(gettext("`hits$start` must be a positive integer."))
  }
  missing_primers <- setdiff(unique(primer_name), names(primers))
  if (length(missing_primers)) {
    cli::cli_abort(gettext(
      "Every `hits$primer` must be a name in `primers`; missing: {.val {missing_primers}}."
    ))
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
#' @family specificity
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
#' table ([pcr_products_batch()]).
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
#'   because that pass used to cost an order of magnitude more than the exact one
#'   (measured 4.6 s exact vs 193 s radius-1 for one pair on a 120 Mb genome,
#'   before the radius-1 batch stopped building a `PDict`; with that fix the near
#'   pass measured 2.0 s against 3.6 s exact, same hits), and the internal
#'   assembly screen does not ask for it.
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
  radius <- as.integer(max_mismatch)
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
#' \dontrun{
#' genome <- ath_genome()
#' snp <- marker_example("SNP")[1, c("chrom", "pos", "ref", "alt")]
#' assay <- assemble_kasp_assay(genome, snp)
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
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' snp <- marker_example("SNP")[1, c("chrom", "pos", "ref", "alt")]
#' assay <- assemble_kasp_assay(genome, snp)
#' assay_pcr_products(assay, genome)
#' }
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

#' BLAST every primer against a local genome database (rBLAST, optional)
#' @param primers named character vector of oligos.
#' @param db_path BLAST database path.
#' @return a tibble of BLAST hits, or `NULL` when the database is missing or the
#'   query fails.
#' @examples
#' \dontrun{
#' blast_primers(c(AS = "GAAGGTGACCAAGTTCATGCT"), "tair10")
#' }
#' @family specificity
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

#' Build (once) a BLAST nucleotide database from a genome
#'
#' Kept as an alias of [marker_blast_index()] so there is one implementation of
#' the build.  The index (`<path>.nhr/.nin/.nsq`) is written next to the FASTA
#' and reused when present, so a second call does not repave a genome-sized
#' file.  This used to write under a fixed `tair10.fa` name for every genome.
#'
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param path FASTA path to write (and index); `NULL` derives a content-keyed
#'   name in the cache.  Ignored for a `marker_genome`, whose own FASTA is
#'   indexed in place.
#' @return the FASTA path the database was built from (the database prefix).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' build_blast_db(genome)
#' }
#' @family specificity
#' @export
build_blast_db <- function(genome, path = NULL) {
  marker_blast_index(genome, path = path)
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
    weight_middle = 1, weight_seed = 4, weight_three_prime = 8, weight_terminal = 12,
    spacer_len_max_mismatch = c(`16` = 0L, `20` = 1L, `23` = 2L, `28` = 0L),
    spacer_len_central = list(`20` = c(7L, 9L))
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
#' ## CRISPR spacer length
#'
#' The `crispr` count above is spacer-independent and stays the fallback: an
#' effector's mismatch budget also depends on the guide's length, and the
#' published series gives an endpoint per length.
#'
#' | spacer_len (nt) | max_mismatch | evidence |
#' |---|---|---|
#' | 16 | 0 | a 16-nt LbuCas13a spacer still required a perfect match (PMC10810210) |
#' | 20 | 1 | a single middle mismatch (positions +7 to +9) costs up to 70 % activity (PMC10810210) |
#' | 23 | 2 | 2 central mismatches still silence the K-RAS G12D transcript 54.8 %; 3 leave 11 % (PMC11654686) |
#' | 28 | 0 | a 28-nt LbuCas13a spacer has no single-mismatch discriminatory capacity (PMC10810210) |
#'
#' The 20-nt row is stricter than "anywhere": its one tolerated mismatch must
#' sit outside spacer positions 7-9 (1-based from the 5' end), the middle window
#' the paper measured.  These numbers stay heuristics (`calibrated = FALSE`):
#' they are published endpoints, not a fit.
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
#' @param spacer_len spacer (guide) length in nt for `"crispr"`; one of the
#'   published series (`16`, `20`, `23`, `28`) sets the mismatch budget and, for
#'   `20`, an intolerant middle window.  An unlisted length keeps the scheme's
#'   `max_mismatch` and reports `spacer_len_known = FALSE`; `NULL` keeps the
#'   length-independent policy.
#' @return a `specificity_policy` list: `mode`, `max_mismatch`, `seed_len`,
#'   `seed_max`, `three_prime_max`, `three_prime_window`, `requires_pam`, `pam`,
#'   `pfs`, `critical_from`/`critical_to` (CRISPR), the position ladder
#'   (`weight_middle`, `weight_seed`, `weight_three_prime`, `weight_terminal`),
#'   `priors` (the tolerance row as a one-row tibble) and `calibrated` (`FALSE`).
#'   A CRISPR policy given a `spacer_len` also carries `spacer_len`,
#'   `spacer_len_known` and `spacer_len_central` (the intolerant middle window,
#'   or `NULL` when that length has none).
#' @examples
#' specificity_policy("kasp")$max_mismatch
#' specificity_policy("rpa")$max_mismatch
#' specificity_policy("cas12a")$requires_pam
#' specificity_policy("homeolog")$mode
#' @family specificity
#' @export
specificity_policy <- function(scheme = "pcr", length = NULL, spacer_len = NULL) {
  if (is.null(scheme)) {
    scheme <- "pcr"
  }
  if (!is.character(scheme) || length(scheme) != 1L) {
    cli::cli_abort(gettext("`scheme` must be a single scheme name."))
  }
  if (!is.null(spacer_len)) {
    n <- if (is.numeric(spacer_len)) spacer_len else NA_real_
    whole <- length(n) == 1L && is.finite(n) && n >= 1 &&
      n <= .Machine$integer.max && n == round(n)
    if (!whole) {
      cli::cli_abort(gettext(
        "`spacer_len` must be a single positive whole number of nucleotides."
      ))
    }
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
  # a published spacer length overrides the scheme budget; an unlisted length
  # keeps it and says so, so a length the evidence does not cover is never
  # guessed at, and `priors` below then carries the overriding number
  if (identical(key, "crispr") && !is.null(spacer_len)) {
    pol$spacer_len <- as.integer(spacer_len)
    name <- as.character(pol$spacer_len)
    if (name %in% names(pol[["spacer_len_max_mismatch"]])) {
      pol$max_mismatch <- pol[["spacer_len_max_mismatch"]][[name]]
      pol$spacer_len_known <- TRUE
    } else {
      pol$spacer_len_known <- FALSE
    }
    pol$spacer_len_central <- if (!is.null(pol[["spacer_len_central"]])) {
      pol[["spacer_len_central"]][[name]]
    } else {
      NULL
    }
  }
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
  if (!is.null(x[["spacer_len"]])) {
    cat("  spacer", x[["spacer_len"]], "nt | max mismatch", x$max_mismatch,
      if (!is.null(x[["spacer_len_central"]])) {
        c(
          "| no mismatch in positions",
          paste0(x[["spacer_len_central"]], collapse = "-")
        )
      },
      "\n",
      sep = " "
    )
  }
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
    crisprverse = namespace_available("crisprVerse")
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
#' product is tolerated *and* it is the only tolerated product;
#' `fail` when every product breaks the policy (including a changed PAM),
#' more than one tolerated product exists, or the backend is primer-level; `unknown`
#' when there is no genome-wide engine - a missing backend is never a silent pass.
#' The result carries the products, the tolerated flag, the best product and the
#' policy it was graded against, so the evidence travels with the verdict, plus
#' `tier`, the evidence level behind the verdict:
#'
#' 1. `"pair"` - pair-level in-silico PCR with product-size bounds, graded here;
#' 2. `"blast"` (`"fm"`) - a primer-level BLAST database / FM-index scan.  The
#'    per-product `n_mismatch` is this pair's primer evidence, so at primer level
#'    it cannot stand in for a genome-wide off-target count: the verdict grades
#'    product formation, and a changed PAM still fails;
#' 3. `"none"` - nothing genome-wide was consulted, so the verdict is `unknown`.
#'
#' [specificity_seed()] never produces a verdict: it is a prefilter, and its
#' fixed false-negative rate is not genome-wide evidence.
#' @param hits in-silico PCR products ([in_silico_pcr()], as a tibble: one row
#'   per product, `chrom`/`start`/`end`/`width` and the mismatch evidence;
#'   extra columns - e.g. a `guide` id from a multi-guide table - are ignored.
#'   An optional `scheme`/`tolerated` pair is recomputed here, so a table from
#'   `in_silico_pcr(scheme = )` is graded exactly as a bare one.
#' @param scheme specificity scheme ([specificity_policy()]).
#' @param genome optional genome; when supplied, a missing engine turns the
#'   answer into `unknown` (see `specificity_backend()`).
#' @param pam optional PAM/PFS observation for `"crispr"` (see `pam_state()`).
#' @param spacer_len optional spacer (guide) length in nt for `"crispr"`; a
#'   length in the published series sets the mismatch budget and, for `20`, the
#'   intolerant middle window (see [specificity_policy()]).
#' @return a `specificity_verdict` list: `verdict` (`"pass"`/`"fail"`/
#'   `"unknown"`), `reason` (an English sentence), `scheme`, `policy` (the policy
#'   graded against, carrying `spacer_len`/`spacer_len_known` for a CRISPR
#'   guide), `n_products`, `n_tolerated`, `n_mismatch` (best product), `best`
#'   (one-row product tibble), `products` (the graded tibble), `pam` (state),
#'   `backend` (the capability list, or `NA` when no genome was supplied) and
#'   `tier` (`"pair"`/`"blast"`/`"fm"`/`"none"`).
#' @examples
#' genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
#'   strrep("A", 20L), "ACGTACGTACGTACGTACGT", strrep("T", 20L),
#'   "ACGTACGTACGTACGTACGT", strrep("A", 20L)
#' )))
#' hits <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
#' # `unknown` when the genome is missing: a missing engine is never a pass
#' specificity_verdict(hits, "kasp", genome = NULL)$verdict
#' specificity_verdict(hits, "kasp", genome = NULL)$tier
#' @family specificity
#' @export
specificity_verdict <- function(hits, scheme = "pcr", genome = NULL, pam = NULL,
                                spacer_len = NULL) {
  pol <- specificity_policy(scheme, spacer_len = spacer_len)
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
  central_bad <- central_violated(products, pol)
  products$tolerated <- pcr_products_tolerated(products, pol) & !central_bad
  products$violations <- mapply(pcr_violations,
    n_mismatch = products$n_mismatch, seed_mismatch = products$seed_mismatch,
    three_prime_mismatch = products$three_prime_mismatch,
    three_prime_terminal = products$three_prime_terminal,
    MoreArgs = list(policy = pol), USE.NAMES = FALSE
  )
  if (any(central_bad)) {
    products$violations[central_bad] <- lapply(
      products$violations[central_bad], function(v) union(v, "central")
    )
  }
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
  # the evidence level: a BLAST/FM index is primer-level (it can fail a primer
  # genome-wide), an in-memory genome is pair-level here, no engine is no tier
  tier <- if (!is.list(backend) || !isTRUE(backend$engine)) {
    "none"
  } else if (identical(backend$engine, "blast")) {
    "blast"
  } else if (identical(backend$engine, "fm")) {
    "fm"
  } else {
    "pair"
  }
  primer_level <- tier %in% c("blast", "fm")
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
      tier = tier,
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
      tier = tier,
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
      tier = tier,
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = n_tol, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else if (primer_level) {
    # the product data is this pair's primer evidence, not an off-target count:
    # promote it to a pass only if a genome-wide primer scan was consulted
    out <- list(
      verdict = "pass",
      reason = gettext(
        "One tolerated product and a genome-wide primer scan: no other binding site."
      ),
      tier = tier,
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = 1L, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else {
    out <- list(
      verdict = "pass",
      reason = gettext(paste(
        "Exactly one tolerated in-silico PCR product; pair-level evidence only,",
        "not a genome-wide primer scan."
      )),
      tier = tier,
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
  cat("  scheme:", x$scheme, "| tier:", x$tier, "| products:", x$n_products,
    "| tolerated:", x$n_tolerated, "\n",
    sep = " "
  )
  cat("  reason:", x$reason, "\n")
  invisible(x)
}
