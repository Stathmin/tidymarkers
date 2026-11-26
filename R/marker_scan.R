# =============================================================================
#  marker_scan.R  --  indexed, genome-scale specificity (BLAST index / FM-style)
# =============================================================================
#  An Rsamtools `.fai` index lets you *read regions* fast, but a whole-genome
#  specificity scan still touches every base - so for large genomes the right
#  tool is a search index.  `specificity_seed()` counts with `Biostrings` (exact
#  `PDict` path, or a chromosome-parallel `countPattern` scan) and is perfect up
#  to genome sizes you can hold in RAM; beyond that, build a search index once
#  and query it: bowtie1 (`marker_bowtie_index()` / `specificity_bowtie()`, the
#  default for a large genome) or an NCBI BLAST database (`marker_blast_index()`
#  / `specificity_blast()`, a heuristic).
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
#' @family specificity
#' @export
has_blast <- function() {
  !is.null(blast_bin())
}

#' Build (once) a BLAST database for a genome
#'
#' For an in-memory `DNAStringSet` or a Bioconductor `BSgenome` the FASTA is written
#' to `path` first (a `BSgenome` goes to a name-keyed file in the session's temp
#' directory, rewritten only when that file is gone); for a `marker_genome` its
#' indexed FASTA is used directly.  A `marker_genome()` built *from* a `BSgenome`
#' has no FASTA, so it gets a named error; attach a database it
#' built elsewhere with `marker_genome(..., blast_db = <path>)` instead.  The
#' returned path is the database to pass to [specificity_blast()]; attach it to a
#' `marker_genome` with `marker_genome(..., blast_db = <path>)` to let
#' [specificity_of()] use it.
#' @param genome a `marker_genome`, a `DNAStringSet` or a `BSgenome`.
#' @param path FASTA path to write an in-memory `DNAStringSet` to; the BLAST
#'   index files are created alongside it and reused when they exist, so delete
#'   the `.nsq` to force a rebuild.  `NULL` derives a content-keyed name in the
#'   session's temp directory.  Ignored for a `marker_genome`, whose own FASTA is
#'   indexed in place, and for a `BSgenome`, which is keyed by its name and
#'   lengths.
#' @return the FASTA path the database was built from.  This is the database
#'   *prefix* (NCBI BLAST+ writes `<path>.nhr/.nin/.nsq` next to it), the value
#'   [specificity_blast()] expects; pass it as-is, do not append an extension.
#' @examples
#' \dontrun{
#' db <- marker_blast_index(marker_genome("genome.fa"))
#' }
#' @family specificity
#' @export
marker_blast_index <- function(genome, path = NULL) {
  # the genome check comes first, so without BLAST+ a BSgenome genome still says so
  src <- genome_fasta_source(genome, "marker_blast_index")
  genome <- src$genome
  fasta <- src$fasta
  rlang::check_installed("rBLAST", reason = gettext("for the indexed genome scan"))
  if (!has_blast()) {
    abort_msg("NCBI BLAST+ ('blastn'/'makeblastdb') was not found.",
      i = "Install it: apt-get install ncbi-blast+ (Debian/Ubuntu), brew install blast (macOS).",
      i = "On Windows run the NCBI BLAST+ installer.",
      i = "Or point the package at an existing build with options(tidymarkers.blast = )."
    )
  }
  if (inherits(genome, "marker_genome")) {
    # build once: BLAST+ writes <fasta>.nhr/.nin/.nsq beside the FASTA, and a
    # 14.8 GB genome must not be re-indexed on every call.  Delete the index to
    # force a rebuild.
    if (!file.exists(paste0(fasta, ".nsq"))) {
      with_blast(rBLAST::makeblastdb(fasta, dbtype = "nucl", verbose = FALSE))
    }
  } else {
    # a content-keyed name, so two in-memory genomes never share a database
    if (is.null(path)) {
      path <- fs::path(fs::path_temp(), paste0("genome-", genome_signature(genome), ".fa"))
    }
    fs::dir_create(fs::path_dir(path))
    if (!file.exists(paste0(path, ".nsq"))) {
      Biostrings::writeXStringSet(genome, path)
      with_blast(rBLAST::makeblastdb(path, dbtype = "nucl", verbose = FALSE))
    }
    fasta <- path
  }
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
#' @return a tibble (`primer`, `seed`, `hits`, `engine`) with the `engine` column
#'   naming the adapter that answered (`"BLAST+"`) and an `engine_stamp`
#'   attribute (see [marker_engine_versions()]).
#' @examples
#' \dontrun{
#' db <- marker_blast_index(marker_genome("genome.fa"))
#' specificity_blast(c(AS = "ACGTACGTACGTACGT"), db)
#' }
#' @family specificity
#' @export
specificity_blast <- function(primers, db_path, k = NA_integer_, BPPARAM = NULL) {
  rlang::check_installed("rBLAST", reason = gettext("for the indexed genome scan"))
  if (!has_blast()) {
    abort_msg(
      "No NCBI BLAST+ found (has_blast() is FALSE).",
      i = "Install it: apt-get install ncbi-blast+ (Debian/Ubuntu), brew install blast (macOS).",
      i = "On Windows run the NCBI BLAST+ installer.",
      i = "Or use specificity_seed() on an in-memory genome.",
      i = "Or point the package at an existing build with options(tidymarkers.blast = )."
    )
  }
  # a per-chromosome set is queried shard by shard and summed, like bowtie
  if (length(db_path) > 1L) {
    rows <- lapply(unname(db_path), function(one) {
      specificity_blast(primers, one, k = k, BPPARAM = BPPARAM)
    })
    return(marker_with_engine(
      tibble::tibble(
        primer = names(primers), seed = rows[[1L]]$seed,
        hits = Reduce(`+`, lapply(rows, `[[`, "hits"))
      ),
      "BLAST+"
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
      abort_msg(
        "Could not open a BLAST database at {db_path} (build it with marker_blast_index()).",
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
  # Count the full-length exact matches only.  `predict.BLAST` returns every
  # HSP, so a 20-mer reports one row per 7-mer seed window (106,924 rows for 16
  # primers whose true total is 6); without this filter the count is the seed
  # multiplicity, not the primer's.
  qlen <- stats::setNames(nchar(as.character(seeds)), names(primers))
  res <- if (is.null(res) || nrow(res) == 0L) {
    res
  } else {
    want <- qlen[as.character(res$qseqid)]
    full <- res$mismatch == 0L & res$length == want &
      res$qstart == 1L & res$qend == want
    res[full, , drop = FALSE]
  }
  cnt <- if (is.null(res) || nrow(res) == 0L) integer(0) else table(res$qseqid)
  hits <- as.integer(cnt[names(primers)])
  hits[is.na(hits)] <- 0L
  marker_with_engine(
    tibble::tibble(primer = names(primers), seed = seeds, hits = hits),
    "BLAST+"
  )
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
      warn_msg("BLAST failed: {conditionMessage(e)}")
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

#' Specificity by the right path for the genome at hand
#'
#' In-memory `DNAStringSet` -> [specificity_seed()] (exact `PDict` path, or a
#' chromosome-parallel scan).  A `marker_genome` with a `bowtie_index` ->
#' [specificity_bowtie()] (the indexed path built for large genomes); with only a
#' `blast_db` -> [specificity_blast()] (a heuristic BLAST index).  This is what
#' the rest of the pipeline calls.
#'
#' The indexed paths answer **different questions**, and the answer's meaning
#' follows the path: [specificity_seed()] counts 3'-seed matches, whereas the
#' indexed paths scan the whole primer (bowtie: both strands, exact by default;
#' BLAST: `blastn-short`).  `k` and `max_mismatch` are honoured only on the
#' in-memory path; a non-default value on a bowtie-index genome warns rather than
#' being silently ignored.
#'
#' A genome carrying only an `fm_index` is refused instead of dispatched: on
#' TAIR10 Chr1 the Subread index costs 31.7 bytes per genome base (wheat 16 Gb
#' would be ~490 GB), its query cost does not fall with `nthreads`, and its BAM
#' carries no mismatch tags to verify a count against; the same site measured
#' the bowtie1 index at 1.34 bytes/base and 5000 primers in 0.1 s, at a cost that
#' does not grow with the genome.
#' Build one with [marker_bowtie_index()], or call [specificity_fm()] explicitly
#' if you want the read mapper.
#'
#' @param primers named character vector of oligos.
#' @param genome a `marker_genome` (with `bowtie_index` or `blast_db`), a
#'   [genome_index()] carrying a declared strategy, or a `DNAStringSet`, or a
#'   **named** list of such genomes (a pangenome or a homoeolog set), answered
#'   one genome at a time.
#' @param k 3'-seed length.
#' @param max_mismatch mismatches allowed (in-memory path only).
#' @param BPPARAM optional `BiocParallel` backend.
#' @return a tibble with one row per primer (`primer`, `seed`, `hits`, `engine`)
#'   and a `scan` attribute naming the question the counts answer: `"3'-seed"` on
#'   the in-memory and BLAST paths, `"whole-primer"` on the bowtie path.  The
#'   `engine` column names the adapter that answered (`"Biostrings"`,
#'   `"BLAST+"`, `"Rbowtie"`), and the `engine_stamp` attribute records the
#'   versions (see [marker_engine_versions()]).  Dispatch is in-memory ->
#'   [specificity_seed()], `bowtie_index` -> [specificity_bowtie()], `blast_db`
#'   -> [specificity_blast()].  A [genome_index()] dispatches on its declared
#'   `strategy` instead of a slot, so `attr(x, "scan")` is always declared.
#'
#'   Given a **named** list, one row per primer **per genome**, with a `genome`
#'   column first and the rows in list order then primer order; `scan` is a named
#'   vector with one entry per genome, so a mixed set (some seeds, some whole
#'   primers) stays readable instead of being averaged.  `engine_stamp` is still
#'   present but holds one stamp carried over from the bound parts — they all
#'   share the package-wide engine set — so it is not a per-genome record; the
#'   `engine` column names the adapter that answered each row.
#' @examples
#' genome <- Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGTACGT"))
#' specificity_of(c(AS = "ACGTACGT"), genome)
#' pangenome <- list(
#'   A = Biostrings::DNAStringSet(c(Chr1 = "AAAACCCCGGGGTTTTAAAACCCC")),
#'   B = Biostrings::DNAStringSet(c(Chr1 = "AAAACCCCGGGGAAAACCCCAAA"))
#' )
#' specificity_of(c(AS = "AAAACCCCGGGGTTTT"), pangenome)
#' @family specificity
#' @export
specificity_of <- function(primers, genome, k = SEED_K, max_mismatch = SEED_MISMATCH,
                           BPPARAM = NULL) {
  # a plain list is a genome *set*; `is.object()` keeps the S4 genomes (and a
  # `marker_genome`, itself a classed list) out of this branch
  if (is.list(genome) && !is.object(genome)) {
    if (is.null(names(genome)) || any(!nzchar(names(genome)))) {
      abort_msg(
        "`genome` must be a named list of genomes, one name per element.",
        i = "Name every genome, e.g. list(A = , B = , D = )."
      )
    }
    parts <- vector("list", length(genome))
    scans <- character(length(genome))
    for (i in seq_along(genome)) {
      part <- specificity_of(primers, genome[[i]],
        k = k, max_mismatch = max_mismatch, BPPARAM = BPPARAM
      )
      scans[[i]] <- attr(part, "scan")
      attr(part, "scan") <- NULL
      parts[[i]] <- part
    }
    names(parts) <- names(genome)
    names(scans) <- names(genome)
    out <- dplyr::bind_rows(parts, .id = "genome")
    attr(out, "scan") <- scans
    class(out) <- c("specificity_scan", class(out))
    return(out)
  }
  if (inherits(genome, "genome_index")) {
    return(genome_index_scan(primers, genome,
      k = k, max_mismatch = max_mismatch, BPPARAM = BPPARAM
    ))
  }
  if (!inherits(genome, "marker_genome")) {
    return(mark_scan(
      specificity_seed(primers, genome,
        k = k, max_mismatch = max_mismatch,
        BPPARAM = BPPARAM
      ),
      "3'-seed"
    ))
  }
  if (!is.null(genome$bowtie_index)) {
    seed_default <- identical(as.integer(k), as.integer(SEED_K)) &&
      identical(as.integer(max_mismatch), as.integer(SEED_MISMATCH))
    if (!seed_default) {
      cli::cli_warn(c(
        gettext("The bowtie path scans the whole primer and ignores `k` / `max_mismatch`."),
        i = gettext("Use specificity_seed() on an in-memory genome for the 3'-seed question.")
      ))
    }
    # a BPPARAM means the caller is already parallel across loci; 4 bowtie
    # threads per worker oversubscribed a 6-worker run by 2.75x
    return(mark_scan(
      specificity_bowtie(primers, genome$bowtie_index,
        nthreads = if (is.null(BPPARAM)) 4L else 1L
      ),
      "whole-primer"
    ))
  }
  if (!is.null(genome$fm_index)) {
    abort_msg(
      "This genome has only an FM index; the whole-genome scan no longer dispatches to it.",
      i = "On TAIR10 the Subread index cost 31.7 bytes/bp and 53 s/query; wheat (16 Gb) needs ~490 GB.",
      i = "Build a bowtie index instead: marker_bowtie_index(genome) (bowtie_index=).",
      i = "Or call specificity_fm() explicitly if you want the read mapper."
    )
  }
  db <- genome$blast_db
  if (is.null(db)) {
    abort_msg(
      "This indexed genome has neither a bowtie index nor a BLAST database attached.",
      i = "Build marker_bowtie_index(genome) (bowtie_index=) or marker_blast_index(genome) (blast_db=)."
    )
  }
  mark_scan(specificity_blast(primers, db, k = k, BPPARAM = BPPARAM), "3'-seed")
}

