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
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' primer <- as.character(Biostrings::subseq(genome[["Chr4"]], 1000, 1019))
#' specificity_seed(c(AS = primer), genome)
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
    abort_msg(
      "This indexed genome has neither a bowtie index nor a BLAST database.",
      i = "Build marker_bowtie_index(genome) (bowtie_index=) or marker_blast_index(genome) (blast_db=)."
    )
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
#' no help - its trusted-band tree is rebuilt per block - so every pattern is
#' searched with the `matchPattern()` engine of `pattern_hits_scan()`, whose hit
#' set is identical, which `tests/testthat/test-specificity.R` locks ("the
#' radius-1 batch search keeps the per-pattern hit set on a block-walked
#' subject").
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
  if (!is.null(pdict) && single_block) {
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
