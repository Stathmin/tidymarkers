# =============================================================================
#  marker_scan.R  --  indexed, genome-scale specificity (BLAST index / FM-style)
# =============================================================================
#  An Rsamtools `.fai` index lets you *read regions* fast, but a whole-genome
#  specificity scan still touches every base - so for large genomes the right
#  tool is a search index.  `specificity_seed()` counts with `Biostrings` (exact
#  `PDict` path, or a chromosome-parallel `countPattern` scan) and is perfect up
#  to genome sizes you can hold in RAM; beyond that, build a BLAST database once
#  (`marker_blast_index()`, NCBI BLAST+ makes its own index - the same idea as an
#  FM-index) and query it with `specificity_blast()`.
#
#  `specificity_of()` picks the right path for the genome object it is given, so
#  the rest of the pipeline stays unaware.
# =============================================================================

#' Is NCBI BLAST+ available (for the indexed genome scan)?
#'
#' Looks for `blastn`/`makeblastdb` in `options(tidymarkers.blast)`, a project
#' `.tools/`, the user cache and finally the `PATH`.
#' @return `TRUE` when both `blastn` and `makeblastdb` are found.
#' @examples
#' has_blast()
#' @family KASP specificity
#' @export
has_blast <- function() {
  !is.null(blast_bin()) || (nzchar(Sys.which("blastn")) && nzchar(Sys.which("makeblastdb")))
}

#' Build (once) a BLAST database for a genome
#'
#' For an in-memory `DNAStringSet` the FASTA is written to `path` first; for a
#' `marker_genome` its indexed FASTA is used directly.  The returned path is the
#' database to pass to [specificity_blast()]; attach it to a `marker_genome` with
#' `marker_genome(..., blast_db = <path>)` to let [specificity_of()] use it.
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param path FASTA path to write a `DNAStringSet` to (ignored for `marker_genome`).
#' @return the BLAST database path (invisibly used by [specificity_blast()]).
#' @examples
#' \dontrun{
#' db <- marker_blast_index(marker_genome("genome.fa"))
#' }
#' @family KASP specificity
#' @export
marker_blast_index <- function(genome, path = fs::path(fs::path_temp(), "marker_genome.fa")) {
  rlang::check_installed("rBLAST", reason = gettext("for the indexed genome scan"))
  if (!has_blast()) {
    cli::cli_abort(c(gettext("NCBI BLAST+ ('blastn'/'makeblastdb') was not found."),
      i = gettext("Install it: apt-get install ncbi-blast+ (Debian/Ubuntu), brew install blast (macOS)."),
      i = gettext("On Windows run the NCBI BLAST+ installer."),
      i = gettext("Or point the package at an existing build with options(tidymarkers.blast = ).")
    ))
  }
  fasta <- if (inherits(genome, "marker_genome")) {
    genome$path
  } else {
    fs::dir_create(fs::path_dir(path))
    Biostrings::writeXStringSet(genome, path)
    path
  }
  with_blast(rBLAST::makeblastdb(fasta, dbtype = "nucl", verbose = FALSE))
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
#' @param db_path BLAST database path (from [marker_blast_index()]).
#' @param k if not `NA`, count only the 3'-terminal `k` bases of each primer
#'   (matching the seed convention); `NA` uses the whole primer.
#' @param BPPARAM optional `BiocParallel` backend; queries are split across workers.
#' @return tibble `primer`, `seed`, `hits`.
#' @examples
#' \dontrun{
#' db <- marker_blast_index(marker_genome("genome.fa"))
#' specificity_blast(c(AS = "ACGTACGTACGTACGT"), db)
#' }
#' @family KASP specificity
#' @export
specificity_blast <- function(primers, db_path, k = NA_integer_, BPPARAM = NULL) {
  rlang::check_installed("rBLAST", reason = gettext("for the indexed genome scan"))
  if (!has_blast()) {
    cli::cli_abort(c(
      gettext("No NCBI BLAST+ found (has_blast() is FALSE)."),
      i = gettext("Install it: apt-get install ncbi-blast+ (Debian/Ubuntu), brew install blast (macOS)."),
      i = gettext("On Windows run the NCBI BLAST+ installer."),
      i = gettext("Or use specificity_seed() on an in-memory genome."),
      i = gettext("Or point the package at an existing build with options(tidymarkers.blast = ).")
    ))
  }
  seeds <- if (is.na(k)) {
    as.character(primers)
  } else {
    unname(substr(primers, nchar(primers) - k + 1L, nchar(primers)))
  }
  query <- Biostrings::DNAStringSet(rlang::set_names(seeds, names(primers)))
  # rBLAST::blast() validates the database itself (globs the index files and runs
  # `blastdbcmd -info`), so a missing or corrupt database surfaces there; we only
  # re-raise it with the pipeline's own hint instead of a hand-rolled .nsq probe.
  bl <- rlang::try_fetch(
    with_blast(rBLAST::blast(db_path, type = "blastn")),
    error = function(e) {
      cli::cli_abort(
        gettext(
          "Could not open a BLAST database at {db_path} (build it with marker_blast_index())."
        ),
        parent = e
      )
    }
  )
  run <- function(q) {
    # a SOCK worker is a fresh process: load rBLAST so predict.BLAST() dispatches
    requireNamespace("rBLAST", quietly = TRUE)
    with_blast(stats::predict(bl,
      newdata = q,
      BLAST_args = c("-task", "blastn-short", "-word_size", "7", "-dust", "no")
    ))
  }
  res <- if (is.null(BPPARAM)) {
    run(query)
  } else {
    # let BiocParallel partition the queries across workers
    purrr::list_rbind(bp_lapply(seq_along(query), function(i) run(query[i]), BPPARAM))
  }
  cnt <- if (is.null(res) || nrow(res) == 0L) integer(0) else table(res$qseqid)
  hits <- as.integer(cnt[names(primers)])
  hits[is.na(hits)] <- 0L
  tibble::tibble(primer = names(primers), seed = seeds, hits = hits)
}

#' Specificity by the right path for the genome at hand
#'
#' In-memory `DNAStringSet` -> [specificity_seed()] (exact `PDict` path, or a
#' chromosome-parallel scan).  Indexed `marker_genome` -> [specificity_blast()]
#' through its BLAST database.  This is what the rest of the pipeline calls.
#'
#' @param primers named character vector of oligos.
#' @param genome a `marker_genome` (with `blast_db`) or a `DNAStringSet`.
#' @param k 3'-seed length.
#' @param max_mismatch mismatches allowed (in-memory path only).
#' @param BPPARAM optional `BiocParallel` backend.
#' @return a tibble with one row per primer (`primer`, `seed`, `hits`); the
#'   in-memory path delegates to [specificity_seed()], the indexed path to
#'   [specificity_blast()].
#' @examples
#' genome <- Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGTACGT"))
#' specificity_of(c(AS = "ACGTACGT"), genome)
#' @family KASP specificity
#' @export
specificity_of <- function(primers, genome, k = SEED_K, max_mismatch = SEED_MISMATCH,
                           BPPARAM = NULL) {
  if (!inherits(genome, "marker_genome")) {
    return(specificity_seed(primers, genome,
      k = k, max_mismatch = max_mismatch,
      BPPARAM = BPPARAM
    ))
  }
  db <- genome$blast_db
  if (is.null(db)) {
    cli::cli_abort(c(
      gettext("This indexed genome has no BLAST database attached."),
      i = gettext("Build one with marker_blast_index(genome) and pass it as `blast_db=`.")
    ))
  }
  specificity_blast(primers, db, k = k, BPPARAM = BPPARAM)
}
