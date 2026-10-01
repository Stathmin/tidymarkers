# =============================================================================
#  kasp_scan.R  --  indexed, genome-scale specificity (BLAST index / FM-style)
# =============================================================================
#  An Rsamtools `.fai` index lets you *read regions* fast, but a whole-genome
#  specificity scan still touches every base - so for large genomes the right
#  tool is a search index.  `specificity_seed()` counts with `Biostrings` (exact
#  `PDict` path, or a chromosome-parallel `countPattern` scan) and is perfect up
#  to genome sizes you can hold in RAM; beyond that, build a BLAST database once
#  (`kasp_blast_index()`, NCBI BLAST+ makes its own index - the same idea as an
#  FM-index) and query it with `specificity_blast()`.
#
#  `specificity_of()` picks the right path for the genome object it is given, so
#  the rest of the pipeline stays unaware.
# =============================================================================

#' Is NCBI BLAST+ available (for the indexed genome scan)?
#' @family KASP specificity
#' @export
has_blast <- function() nzchar(Sys.which("blastn")) && nzchar(Sys.which("makeblastdb"))

#' Build (once) a BLAST database for a genome
#'
#' For an in-memory `DNAStringSet` the FASTA is written to `path` first; for a
#' `kasp_genome` its indexed FASTA is used directly.  The returned path is the
#' database to pass to [specificity_blast()]; attach it to a `kasp_genome` with
#' `kasp_genome(..., blast_db = <path>)` to let [specificity_of()] use it.
#' @param genome a `kasp_genome` or a `DNAStringSet`.
#' @param path FASTA path to write a `DNAStringSet` to (ignored for `kasp_genome`).
#' @return the BLAST database path (invisibly used by [specificity_blast()]).
#' @family KASP specificity
#' @export
kasp_blast_index <- function(genome, path = fs::path(fs::path_temp(), "kasp_genome.fa")) {
  rlang::check_installed("rBLAST", reason = gettext("for the indexed genome scan"))
  if (!has_blast()) {
    rlang::abort(c(gettext("NCBI BLAST+ ('blastn'/'makeblastdb') not found on the PATH."),
                   i = gettext("install NCBI BLAST+ (e.g. `apt-get install ncbi-blast+`).")))
  }
  fasta <- if (inherits(genome, "kasp_genome")) genome$path else {
    fs::dir_create(fs::path_dir(path))
    Biostrings::writeXStringSet(genome, path)
    path
  }
  rBLAST::makeblastdb(fasta, dbtype = "nucl", verbose = FALSE)
  fasta
}

#' Genome-wide primer hit counts via a BLAST database (the indexed path)
#'
#' `blastn -task blastn-short` searches both strands, so the per-primer hit count
#' is directly comparable to [specificity_seed()]'s both-strand count.  The
#' thresholds are BLAST's (word size 7, no `dust`), i.e. this is a heuristic
#' index-based screen, not an exact count.
#'
#' @param primers named character vector of oligos.
#' @param db_path BLAST database path (from [kasp_blast_index()]).
#' @param k if not `NA`, count only the 3'-terminal `k` bases of each primer
#'   (matching the seed convention); `NA` uses the whole primer.
#' @param BPPARAM optional `BiocParallel` backend; queries are split across workers.
#' @return tibble `primer`, `seed`, `hits`.
#' @family KASP specificity
#' @export
specificity_blast <- function(primers, db_path, k = NA_integer_, BPPARAM = NULL) {
  rlang::check_installed("rBLAST", reason = gettext("for the indexed genome scan"))
  if (!has_blast()) {
    rlang::abort(c(
      gettext("No NCBI BLAST+ on PATH (has_blast() is FALSE)."),
      i = gettext("Install BLAST+, or use specificity_seed() on an in-memory genome.")))
  }
  seeds <- if (is.na(k)) as.character(primers)
           else unname(substr(primers, nchar(primers) - k + 1L, nchar(primers)))
  query <- Biostrings::DNAStringSet(stats::setNames(seeds, names(primers)))
  # rBLAST::blast() validates the database itself (globs the index files and runs
  # `blastdbcmd -info`), so a missing or corrupt database surfaces there; we only
  # re-raise it with the pipeline's own hint instead of a hand-rolled .nsq probe.
  bl <- rlang::try_fetch(
    rBLAST::blast(db_path, type = "blastn"),
    error = function(e) rlang::abort(
      glue::glue(gettext(
        "Could not open a BLAST database at {db_path} (build it with kasp_blast_index()).")),
      parent = e))
  run <- function(q) {
    # a SOCK worker is a fresh process: load rBLAST so predict.BLAST() dispatches
    requireNamespace("rBLAST", quietly = TRUE)
    stats::predict(bl, newdata = q,
                   BLAST_args = c("-task", "blastn-short", "-word_size", "7", "-dust", "no"))
  }
  res <- if (is.null(BPPARAM)) {
    run(query)
  } else {
    nw <- max(1L, BiocParallel::bpnworkers(BPPARAM))
    idx <- split(seq_along(query), cut(seq_along(query), breaks = nw))
    purrr::list_rbind(BiocParallel::bplapply(idx, function(ii) run(query[ii]), BPPARAM = BPPARAM))
  }
  cnt <- if (is.null(res) || nrow(res) == 0L) integer(0) else table(res$qseqid)
  hits <- as.integer(cnt[names(primers)]); hits[is.na(hits)] <- 0L
  tibble::tibble(primer = names(primers), seed = seeds, hits = hits)
}

#' Specificity by the right path for the genome at hand
#'
#' In-memory `DNAStringSet` -> [specificity_seed()] (exact `PDict` path, or a
#' chromosome-parallel scan).  Indexed `kasp_genome` -> [specificity_blast()]
#' through its BLAST database.  This is what the rest of the pipeline calls.
#'
#' @param primers named character vector of oligos.
#' @param genome a `kasp_genome` (with `blast_db`) or a `DNAStringSet`.
#' @param k 3'-seed length.
#' @param max_mismatch mismatches allowed (in-memory path only).
#' @param BPPARAM optional `BiocParallel` backend.
#' @family KASP specificity
#' @export
specificity_of <- function(primers, genome, k = SEED_K, max_mismatch = SEED_MISMATCH,
                           BPPARAM = NULL) {
  if (!inherits(genome, "kasp_genome")) {
    return(specificity_seed(primers, genome, k = k, max_mismatch = max_mismatch,
                            BPPARAM = BPPARAM))
  }
  db <- genome$blast_db
  if (is.null(db)) {
    rlang::abort(c(
      gettext("This indexed genome has no BLAST database attached."),
      i = gettext("Build one with kasp_blast_index(genome) and pass it as `blast_db=`.")))
  }
  specificity_blast(primers, db, k = k, BPPARAM = BPPARAM)
}
