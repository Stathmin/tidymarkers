# =============================================================================
#  dev/specificity-oracle.R  --  the specificity screens against ground truth
# =============================================================================
#  PolyMarker is not a ground truth: it is another implementation with its own
#  heuristics.  On TAIR10 the right answer is computable, so this script
#  enumerates it and scores the package's three specificity paths against the
#  enumeration:
#
#    ORACLE 1  3'-seed screen.  `specificity_seed()` counts, per primer, genome
#              occurrences of the last 15 nt (SEED_K) with <= 1 substitution,
#              on both strands.  The truth is the same definition evaluated
#              exhaustively: the seed itself plus all 45 single-substitution
#              variants, counted exactly on both strands.  Equality is asserted
#              base for base, not within a tolerance.
#
#    ORACLE 2  in-silico PCR, exact mode.  `in_silico_pcr(max_mismatch = 0)`
#              pairs every exact forward hit with every downstream exact
#              reverse-complement hit at least `max(len(F), len(R))` bases
#              away and dedupes on (chrom, start, end).  The truth is that
#              pairing done by enumeration with `vmatchPattern()`.
#
#    ORACLE 3  BLAST+ as a BASELINE, not an oracle.  The same primers go
#              through `blastn -task blastn-short` on an index of the same
#              genome; what matters is recall of the enumerated hits and the
#              wall clock beside the in-memory scan.
#
#  Panel: 8 deterministic 20-mers taken from TAIR10 itself (seeded starts,
#  N-free), so the counts are non-trivial and reproducible; plus 2 designed
#  primer pairs (a real assay's forward and reverse) for ORACLE 2.
#
#  Run:  Rscript dev/specificity-oracle.R
#  Needs: makeblastdb, blastn (BLAST+); writes .tmp/ only (git-ignored).
#  Results: results/specificity-oracle.csv (git-ignored).
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

need <- c("makeblastdb", "blastn")
have <- Sys.which(need)
if (any(!nzchar(have))) stop("missing binary/binaries: ", paste(need[!nzchar(have)], collapse = ", "))

suppressMessages(pkgload::load_all(root, quiet = TRUE))
res_dir <- fs::path(root, "results")
fs::dir_create(res_dir)
tmp <- fs::path(root, ".tmp")
fs::dir_create(tmp)

genome <- ath_genome()
cat(sprintf("genome: %d sequences, %s Mb\n", length(genome),
  formatC(sum(Biostrings::width(genome)) / 1e6, format = "f", digits = 1)))

fmt <- function(x, d = 2L) formatC(x, format = "f", digits = d)
rc <- function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))

# All overlapping occurrences of each pattern on one chromosome.
count_in <- function(chr, pats) {
  sum(vapply(seq_along(pats), function(k) {
    length(Biostrings::matchPattern(pats[[k]], chr))
  }, integer(1)))
}

# --- the panel ---------------------------------------------------------------
set.seed(20261023L)
N_PRIMERS <- 8L
panel <- character(0)
i <- 0L
while (length(panel) < N_PRIMERS && i < 500L) {
  i <- i + 1L
  ch <- names(genome)[(i %% 5L) + 1L]
  st <- 1e5 + sample(seq_len(1e5), 1L)
  s <- as.character(Biostrings::subseq(genome[[ch]], st, st + 19L))
  if (!grepl("[^ACGT]", s) && !s %in% panel) panel <- c(panel, s)
}
panel <- rlang::set_names(panel, paste0("p", seq_along(panel)))

# --- ORACLE 1: the 3'-seed screen, enumerated --------------------------------
single_sub <- function(s) {
  b <- strsplit(s, "")[[1L]]
  out <- character(0)
  for (i in seq_along(b)) {
    for (alt in setdiff(c("A", "C", "G", "T"), b[[i]])) {
      x <- b
      x[[i]] <- alt
      out <- c(out, paste(x, collapse = ""))
    }
  }
  unique(out)
}

seed_truth <- function(seed) {
  pats <- Biostrings::DNAStringSet(unique(c(seed, single_sub(seed))))
  sum(vapply(seq_along(genome), function(i) {
    chr <- genome[[i]]
    count_in(chr, pats) + count_in(chr, Biostrings::reverseComplement(pats))
  }, numeric(1)))
}

t0 <- proc.time()[["elapsed"]]
truth_hits <- vapply(substr(panel, nchar(panel) - SEED_K + 1L, nchar(panel)), seed_truth, numeric(1))
t_truth <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
fast <- specificity_seed(panel, genome)
t_fast <- proc.time()[["elapsed"]] - t0

cat("\n== ORACLE 1: 3'-seed screen (k = 15, <= 1 substitution, both strands) ==\n")
cat(sprintf("primers %d; enumerated %s s, screen %s s\n", length(panel),
  fmt(t_truth, 2L), fmt(t_fast, 2L)))
cmp <- tibble::tibble(
  primer = names(panel), seq = unname(panel), truth = as.numeric(truth_hits),
  screen = as.numeric(fast$hits)
)
cmp$equal <- cmp$truth == cmp$screen
print(as.data.frame(cmp), row.names = FALSE)
cat(sprintf("exact agreement: %d/%d\n", sum(cmp$equal), nrow(cmp)))

# --- ORACLE 2: in-silico PCR, exact mode, enumerated -------------------------
# Two real primer pairs: a 20-mer from the panel against its own reverse
# complement's neighbourhood is not a product, so the pairs are the panel
# primer and a second panel primer (a mismatch-free pair is rare in a 120 Mb
# genome, which is exactly what makes the zero row informative), plus one
# designed pair from Tidymarkers on Chr4.
pcr_truth <- function(f, r) {
  min_width <- max(nchar(f), nchar(r))
  revc <- rc(r)
  total <- 0L
  for (i in seq_along(genome)) {
    chr <- genome[[i]]
    fpos <- Biostrings::start(Biostrings::matchPattern(f, chr))
    rpos <- Biostrings::start(Biostrings::matchPattern(revc, chr))
    if (!length(fpos) || !length(rpos)) next
    rpos <- sort(unique(rpos))
    fpos <- sort(unique(fpos))
    # one product per (forward start, reverse end) with width >= min_width
    total <- total + sum(vapply(fpos, function(s) sum(rpos >= s + min_width - 1L), numeric(1)))
  }
  as.integer(total)
}

pairs <- list(
  
  chr1_span120 = c(
    as.character(Biostrings::subseq(genome[["Chr1"]], 200000L, 200019L)),
    rc(as.character(Biostrings::subseq(genome[["Chr1"]], 200120L, 200139L)))
  ),
  chr3_span500 = c(
    as.character(Biostrings::subseq(genome[["Chr3"]], 1500000L, 1500019L)),
    rc(as.character(Biostrings::subseq(genome[["Chr3"]], 1500500L, 1500519L)))
  ),
  chr5_span60 = c(
    as.character(Biostrings::subseq(genome[["Chr5"]], 3000000L, 3000019L)),
    rc(as.character(Biostrings::subseq(genome[["Chr5"]], 3000060L, 3000079L)))
  )
)

t0 <- proc.time()[["elapsed"]]
pcr <- purrr::imap_dfr(pairs, function(pr, nm) {
  truth <- pcr_truth(pr[[1L]], pr[[2L]])
  t_a <- proc.time()[["elapsed"]]
  got <- suppressWarnings(in_silico_pcr(pr[[1L]], pr[[2L]], genome, max_mismatch = 0L))
  tibble::tibble(
    pair = nm, f = pr[[1L]], r = pr[[2L]],
    truth = truth, found = nrow(got), equal = truth == nrow(got)
  )
})
t_pcr <- proc.time()[["elapsed"]] - t0
cat("\n== ORACLE 2: exact in-silico PCR, enumerated ==\n")
print(as.data.frame(pcr), row.names = FALSE)
cat(sprintf("exact agreement: %d/%d; %s s for the three pairs\n",
  sum(pcr$equal), nrow(pcr), fmt(t_pcr, 2L)))

# --- ORACLE 3: BLAST+ baseline ----------------------------------------------
fa <- fs::path(tmp, "tair10.fa")
Biostrings::writeXStringSet(genome, fa)
db <- fs::path(tmp, "tair10db")
if (any(!fs::file_exists(paste0(db, c(".nhr", ".nin", ".nsq"))))) {
  system2("makeblastdb", c("-in", fa, "-dbtype", "nucl", "-out", db),
    stdout = TRUE, stderr = TRUE)
}
qfa <- fs::path(tmp, "specificity-primers.fa")
Biostrings::writeXStringSet(
  Biostrings::DNAStringSet(rlang::set_names(unname(panel), names(panel))), qfa, width = 60L
)

t0 <- proc.time()[["elapsed"]]
blast <- system2("blastn", c(
  "-task", "blastn-short", "-word_size", "7", "-evalue", "10", "-dust", "no",
  "-db", db, "-query", qfa,
  "-outfmt", shQuote("6 qseqid sseqid sstart send length nident mismatch gaps")
), stdout = TRUE, stderr = FALSE)
t_blast <- proc.time()[["elapsed"]] - t0
b <- utils::read.delim(text = blast, header = FALSE, stringsAsFactors = FALSE,
  col.names = c("q", "sseqid", "sstart", "send", "length", "nident", "mismatch", "gaps"))
# a full-length exact hit is what the enumeration counts
exact <- dplyr::filter(b, .data$length == 20L, .data$nident == 20L, .data$gaps == 0L)
qid <- sub(" .*$", "", exact$q)
blast_counts <- table(factor(qid, levels = names(panel)))
cat("\n== ORACLE 3: BLAST+ baseline (blastn-short, word 7, evalue 10) ==\n")
cat(sprintf("HSPs %d, full-length exact %d, %s s\n", nrow(b), nrow(exact), fmt(t_blast, 2L)))
blast_cmp <- tibble::tibble(
  primer = names(panel),
  truth_full = vapply(unname(panel), function(s) {
    sum(vapply(seq_along(genome), function(i) {
      chr <- genome[[i]]
      count_in(chr, Biostrings::DNAStringSet(s)) +
        count_in(chr, Biostrings::DNAStringSet(rc(s)))
    }, numeric(1)))
  }, numeric(1)),
  blast = as.integer(blast_counts)
)
blast_cmp$recall <- blast_cmp$blast / pmax(blast_cmp$truth_full, 1)
print(as.data.frame(blast_cmp), row.names = FALSE)
cat(sprintf("BLAST recall over the panel: min %s, mean %s; %s s vs %s s for the seed screen\n",
  fmt(min(blast_cmp$recall), 3L), fmt(mean(blast_cmp$recall), 3L),
  fmt(t_blast, 2L), fmt(t_fast, 2L)))

# --- the numbers -------------------------------------------------------------
out <- dplyr::bind_rows(
  tibble::tibble(test = "seed", item = cmp$primer, truth = cmp$truth,
    found = cmp$screen, equal = cmp$equal, seconds = NA_real_),
  tibble::tibble(test = "pcr", item = pcr$pair, truth = pcr$truth,
    found = pcr$found, equal = pcr$equal, seconds = NA_real_),
  tibble::tibble(test = "blast_full_length", item = blast_cmp$primer,
    truth = blast_cmp$truth_full, found = blast_cmp$blast,
    equal = NA, seconds = NA_real_)
)
out$seconds[[1L]] <- t_fast
out$seconds[[which(out$test == "pcr")[1L]]] <- t_pcr
out$seconds[[which(out$test == "blast_full_length")[1L]]] <- t_blast
utils::write.csv(out, file.path(res_dir, "specificity-oracle.csv"), row.names = FALSE)

cat("\n===== SUMMARY =====\n")
print(round(c(
  n_primers = length(panel),
  seed_exact = sum(cmp$equal),
  seed_n = nrow(cmp),
  seed_seconds = t_fast,
  seed_truth_seconds = t_truth,
  pcr_exact = sum(pcr$equal),
  pcr_n = nrow(pcr),
  blast_min_recall = min(blast_cmp$recall),
  blast_seconds = t_blast
), 3L))
cat("wrote results/specificity-oracle.csv\n")
