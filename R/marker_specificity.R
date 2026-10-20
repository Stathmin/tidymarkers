# =============================================================================
#  marker_specificity.R  --  genome-wide specificity, in-silico PCR
# =============================================================================
# =============================================================================
#  6. Genome-wide specificity (pure R, plus optional BLAST)
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

#' In-silico PCR: every product a primer pair would amplify
#'
#' Uses IUPAC-aware `Biostrings::matchPattern()` to find forward and reverse
#' primer matches, then pairs every forward start with a downstream reverse end
#' into a product.  This supports degenerate primers and `max_mismatch`; a KASP
#' test must amplify exactly one product per allele.
#'
#' @param f_primer forward primer (5'->3'), binds the plus strand.
#' @param r_primer reverse primer (5'->3').
#' @param genome DNAStringSet.
#' @param max_mismatch mismatches tolerated in the primer match.
#' @param BPPARAM optional `BiocParallel` backend; the virtual PCR is split by
#'   genome sequence (one worker per chromosome with `MulticoreParam()`).
#' @return tibble `chrom`, `start`, `end`, `width` (one row per product).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' in_silico_pcr("GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC", genome)
#' }
#' @family KASP specificity
#' @export
in_silico_pcr <- function(f_primer, r_primer, genome, max_mismatch = 0L, BPPARAM = NULL) {
  fwd <- Biostrings::DNAString(f_primer)
  # the reverse primer binds the minus strand, so its reverse complement is the
  # plus-strand pattern to search
  rev <- Biostrings::reverseComplement(Biostrings::DNAString(r_primer))
  min_width <- max(nchar(f_primer), nchar(r_primer))
  chroms <- genome_chroms(genome)
  scan_one <- function(nm) {
    worker_biostrings()
    subject <- genome_seq(genome, nm)
    # matchPattern(fixed = FALSE) is IUPAC-aware and supports max.mismatch;
    # matchProbePair() returns individual probe hits, not paired products.
    fh <- Biostrings::matchPattern(fwd, subject, max.mismatch = max_mismatch, fixed = FALSE)
    rh <- Biostrings::matchPattern(rev, subject, max.mismatch = max_mismatch, fixed = FALSE)
    if (length(fh) == 0L || length(rh) == 0L) {
      return(NULL)
    }
    fs <- IRanges::start(fh)
    re <- IRanges::end(rh)
    idx <- which(outer(fs, re, "<="), arr.ind = TRUE)
    if (!nrow(idx)) {
      return(NULL)
    }
    out <- tibble::tibble(
      chrom = nm, start = fs[idx[, 1L]], end = re[idx[, 2L]]
    )
    out$width <- out$end - out$start + 1L
    out <- out[out$width >= min_width, , drop = FALSE]
    if (!nrow(out)) NULL else dplyr::distinct(out)
  }
  res <- if (is.null(BPPARAM)) {
    lapply(chroms, scan_one)
  } else {
    bp_lapply(chroms, scan_one, BPPARAM)
  }
  out <- purrr::list_rbind(res)
  if (!nrow(out)) {
    out <- tibble::tibble(
      chrom = character(), start = integer(), end = integer(), width = integer()
    )
  }
  # Exact mode must not silently hide a one-mismatch primer product.  Scan one
  # mismatch and warn; callers who want those products in the table pass
  # max_mismatch = 1 explicitly.
  if (max_mismatch == 0L) {
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
