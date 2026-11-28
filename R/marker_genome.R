# =============================================================================
#  marker_genome.R  --  an indexed genome source (no whole-genome DNA in RAM)
# =============================================================================
#  The design window (a few hundred bp around a variant) is all the *designer*
#  needs; the genome-wide verification steps are the only ones that must see the
#  whole genome, and those have their own indexed path (see marker_scan.R: a
#  bowtie1 index by default, BLAST as the heuristic fallback).
#
#  A `marker_genome` is an Rsamtools-indexed FASTA you read regions from on demand
#  (`genome_seq()`), so `snp_template()`, `strand_gc_hint()` and the reference
#  lookup in `pick_target_snp()` work without loading chromosomes into memory.
#  Sequence names and lengths come from the library index (`seqinfo()`), not from
#  hand-parsing the `.fai`.  Every accessor also accepts a plain in-memory
#  `DNAStringSet`, so the two representations are interchangeable.
# =============================================================================

#' Indexed genome source backed by an Rsamtools FASTA index
#'
#' Wraps an (optionally compressed) FASTA so regions can be read on demand with
#' [Rsamtools::scanFa()], instead of holding whole chromosomes as a `DNAStringSet`
#' in memory.  The accessors [genome_chroms()], [genome_length()] and
#' [genome_seq()] accept both a `marker_genome` and a plain `DNAStringSet`, so the
#' design code does not care which one it is given.
#'
#' A `marker_genome` can also be built from a Bioconductor `BSgenome` object
#' (`bsgenome =`), in which case regions are read with `BSgenome::getSeq()` and no
#' FASTA needs to live on disk.
#'
#' Whole-sequence reads are cached in the object the first time they are asked
#' for: the genome-wide verification steps ask for a whole chromosome once per
#' primer pair, and re-reading it from the FASTA (not the scan) was the cost (40
#' of 58 CPU-seconds in a measured 6-locus design on TAIR10).  The cache costs
#' one chromosome of RAM per chromosome actually scanned; set
#' `options(tidymarkers.genome_cache = FALSE)` before building the genome to
#' turn it off on a memory-tight machine.
#'
#' @param path FASTA path (gzipped is decompressed once, see [fasta_index()]).
#'   Must be `NULL` when `bsgenome` is given: a BSgenome genome has no FASTA on
#'   disk, so [marker_blast_index()], [marker_bowtie_index()] and
#'   [marker_fm_index()] raise an error for it, and the index must be built
#'   elsewhere and attached with `blast_db=`, `bowtie_index=` or `fm_index=`.
#' @param chrom_map optional named character vector `defline -> name`
#'   (e.g. `ATH_CHROM_MAP`); when given, only mapped sequences are kept.
#' @param blast_db optional prebuilt BLAST database path for the indexed
#'   genome-wide scan (see [marker_blast_index()] / [specificity_blast()]).
#' @param fm_index optional prebuilt FM-index basename for the indexed
#'   genome-wide scan (see [marker_fm_index()] / [specificity_fm()]).
#' @param bowtie_index optional prebuilt bowtie1 index basename for the indexed
#'   genome-wide scan (see [marker_bowtie_index()] / [specificity_bowtie()]);
#'   this is the indexed path [specificity_of()] prefers.
#' @param build optional release id (`"A"`, `"B"`, `"D"`, `"TAIR10"`).  One
#'   genome may omit it; a pangenome component must name it, so a scan of that
#'   component carries the build and builds cannot be pooled by accident.
#' @param bsgenome a `BSgenome` object (or a package name) to read from instead
#'   of a FASTA (e.g. `"BSgenome.Athaliana.TAIR.TAIR9"`).  Sequencing works, but
#'   a BSgenome genome has no FASTA on disk, so the index builders
#'   [marker_blast_index()], [marker_bowtie_index()] and [marker_fm_index()]
#'   raise an error for it; pass an index you built elsewhere with
#'   `blast_db=`, `bowtie_index=` or `fm_index=`.
#' @return a `marker_genome` object.
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' marker_genome(fa, chrom_map = c("1" = "Chr1"))
#' @family input
#' @export
marker_genome <- function(path = NULL, chrom_map = NULL, blast_db = NULL,
                          bsgenome = NULL, fm_index = NULL, bowtie_index = NULL,
                          build = NULL) {
  build <- check_genome_build(build)
  if (!is.null(bsgenome)) {
    return(bsgenome_kasp_genome(
      bsgenome, chrom_map, blast_db, fm_index, bowtie_index, build
    ))
  }
  plain <- fasta_index(path)
  fa <- Rsamtools::FaFile(plain)
  si <- GenomeInfoDb::seqinfo(fa) # names + lengths from the .fai
  raw <- as.character(GenomeInfoDb::seqnames(si))
  len <- rlang::set_names(as.integer(GenomeInfoDb::seqlengths(si)), raw)
  norm <- if (!is.null(chrom_map)) unname(chrom_map[raw]) else raw
  keep <- if (!is.null(chrom_map)) !is.na(norm) else rep(TRUE, length(raw))
  seqnames <- rlang::set_names(raw[keep], (if (!is.null(chrom_map)) norm else raw)[keep])
  structure(
    list(
      path = plain, fa = fa, seqnames = seqnames,
      lengths = rlang::set_names(unname(len[raw[keep]]), names(seqnames)),
      blast_db = blast_db, fm_index = fm_index, bowtie_index = bowtie_index,
      build = build,
      cache = if (isTRUE(getOption("tidymarkers.genome_cache", TRUE))) new.env(parent = emptyenv()) else NULL
    ),
    class = "marker_genome"
  )
}

# The declared capability matrix (ADR 0005).  It is data so a caller can ask
# what a genome can answer instead of inferring it from populated index slots.
GENOME_INDEX_CAPABILITIES <- list(
  memory = list(
    seed_scan = TRUE, whole_primer_scan = FALSE,
    mismatch_tolerant = TRUE, max_mismatch = Inf, engine = "Biostrings"
  ),
  bowtie = list(
    seed_scan = TRUE, whole_primer_scan = TRUE,
    mismatch_tolerant = TRUE, max_mismatch = 3, engine = "Rbowtie"
  ),
  blast = list(
    seed_scan = TRUE, whole_primer_scan = TRUE,
    mismatch_tolerant = FALSE, max_mismatch = 0, engine = "BLAST+"
  ),
  fm = list(
    seed_scan = FALSE, whole_primer_scan = TRUE,
    mismatch_tolerant = TRUE, max_mismatch = 2, engine = "Rsubread"
  )
)

#' Attach a declared genome-access strategy to a genome
#'
#' ADR 0005 makes the genome a strategy, not a set of optional index slots.
#' `genome_index()` records which strategy will answer a genome-wide scan and
#' what that strategy can do; [specificity_of()] dispatches on `strategy`
#' instead of guessing from a non-`NULL` slot, and tags the result with
#' `attr(x, "scan")`.  A missing index is built lazily on the first query, in
#' `dir`, so constructing the object never writes a multi-gigabyte index.
#'
#' `manifest` is one row per shard (`shard`, `chrom`, `path`, `bytes`,
#' `sha256`, `tool`, `tool_version`, `build_id`, `genome_build`,
#' `source_sha256`, `created`).  The constructor fills the paths and the tool
#' version but not the hashes, so opening a 16 Gb genome does not read it;
#' Chromosome shards of a bowtie index with `chroms`; above
#' `BOWTIE_MONOLITH_LIMIT` the shards are chosen even when `chroms` is `NULL`.
#' `index_manifest()` refreshes the measured columns from disk.  Above
#' `BOWTIE_MONOLITH_LIMIT` (21.5 GB, 1.34 bytes per base) a bowtie strategy with
#' `chroms = NULL` is sharded per chromosome, and each shard is built from one
#' chromosome at a time (`index_chrom_seq()`), never from the whole FASTA.
#'
#' The capability matrix is:
#'
#' | strategy | `seed_scan` | `whole_primer_scan` | `mismatch_tolerant` | `max_mismatch` | engine |
#' | -------- | ----------- | ------------------- | ------------------- | -------------- | ------ |
#' | memory   | yes         | no                  | yes                 | any            | Biostrings |
#' | bowtie   | yes         | yes                 | yes                 | 3              | Rbowtie |
#' | blast    | yes (exact) | yes (exact)         | no                  | 0              | BLAST+ |
#' | fm       | no          | yes                 | yes                 | 2              | Rsubread |
#'
#' @param genome a `DNAStringSet` for `strategy = "memory"`, or a
#'   `marker_genome`, `DNAStringSet` or `BSgenome` for an index strategy.
#' @param strategy the declared strategy; never inferred at query time.
#' @param chroms optional sequence names for a per-chromosome build (bowtie and
#'   BLAST); `NULL` leaves the whole-genome default, which a bowtie genome above
#'   `BOWTIE_MONOLITH_LIMIT` is sharded into per-chromosome shards anyway.
#' @param dir directory for a lazily built index; defaults to [tempdir()].
#' @param threads threads for an index build.
#' @param build optional release id; when given it is written into every
#'   manifest row (and into the shard paths), and a `genome_index` over a
#'   pangenome component must have one.  A plain list is not an index: build one
#'   index per element.
#' @return a `genome_index` object with `strategy`, `capabilities`, `manifest`,
#'   `build` and the wrapped `genome`.
#' @examples
#' genome <- Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGTACGT"))
#' idx <- genome_index(genome, strategy = "memory")
#' idx$strategy
#' capabilities(idx)
#' @family input
#' @export
genome_index <- function(genome, strategy = c("memory", "bowtie", "blast", "fm"),
                         chroms = NULL, dir = NULL, threads = 1L, build = NULL) {
  threads <- check_worker_count(threads, "threads")
  strategy <- match.arg(strategy)
  build <- check_genome_build(build)
  if (is.list(genome) && !is.object(genome)) {
    abort_msg(
      "A pangenome is one index per build, not one index over a list.",
      i = "Call genome_index() once per element, each with its own `build =`."
    )
  }
  if (!is.null(build)) {
    if (inherits(genome, "marker_genome")) {
      genome$build <- build
    } else if (inherits(genome, "DNAStringSet")) {
      md <- S4Vectors::metadata(genome) %||% list()
      md$genome_build <- build
      S4Vectors::metadata(genome) <- md
    } else {
      abort_msg("`build` needs a {.cls marker_genome} or a {.cls DNAStringSet}.")
    }
  }
  if (identical(strategy, "memory") && !inherits(genome, "DNAStringSet")) {
    abort_msg(
      "The memory strategy needs an in-memory {.cls DNAStringSet}.",
      i = "Use an index strategy for a {.cls marker_genome} or a {.cls BSgenome}."
    )
  }
  # a monolithic wheat index does not fit here, so shard it by chromosome and
  # remember that in the manifest the caller can read
  if (is.null(chroms) && strategy %in% c("bowtie", "blast") && index_should_shard(genome, strategy)) {
    chroms <- names(if (inherits(genome, "DNAStringSet")) genome else genome$seqnames)
  }
  structure(
    list(
      genome = genome, strategy = strategy,
      capabilities = GENOME_INDEX_CAPABILITIES[[strategy]],
      chroms = chroms, dir = dir, threads = threads,
      build = index_genome_build(genome),
      manifest = index_shards(genome, strategy, chroms, dir)
    ),
    class = "genome_index"
  )
}

#' @rdname capabilities
#' @export
capabilities.genome_index <- function(engine, capability = NULL, ...) {
  capability_value(engine, capability)
}

# Resolve (and, for an index strategy, lazily build) the handle a scan queries.
# Shards come from the manifest, so a per-chromosome set is built one
# chromosome at a time and an existing shard is reused.
genome_index_handle <- function(index) {
  index_build(index)
}

# Build a lazy index once, before an outer parallel pool starts.  A worker then
# only opens the shards the manifest names; it never races another worker to
# build the same shard.  Returns the input unchanged for a plain genome.
index_prebuild <- function(genome) {
  if (inherits(genome, "genome_index")) {
    invisible(index_build(genome))
  }
  genome
}

# The dispatcher: choose the declared strategy and announce the question the
# returned counts answer.  A capability the request exceeds is an error, never a
# warning followed by a different question (D2).
genome_index_scan <- function(primers, index, k = SEED_K,
                              max_mismatch = SEED_MISMATCH, BPPARAM = NULL) {
  caps <- index$capabilities
  if (!is.na(max_mismatch) && max_mismatch > caps$max_mismatch) {
    abort_msg(
      "The {index$strategy} strategy supports `max_mismatch <= {caps$max_mismatch}`, not {max_mismatch}.",
      i = "Its declared capability is `mismatch_tolerant = {caps$mismatch_tolerant}`.",
      class = "tidymarkers_capability_error"
    )
  }
  out <- switch(index$strategy,
    memory = mark_scan(
      specificity_seed(primers, index$genome,
        k = k, max_mismatch = max_mismatch, BPPARAM = BPPARAM
      ),
      "3'-seed"
    ),
    bowtie = {
      if (!identical(as.integer(k), as.integer(SEED_K))) {
        abort_msg(
          "The bowtie strategy scans the whole primer: it cannot honour `k = {k}`.",
          i = "Its declared capabilities are `seed_scan = {caps$seed_scan}` and `whole_primer_scan`.",
          i = "The whole primer is scanned, so there is no seed window to shorten.",
          i = "Use strategy = \"memory\" for the 3'-seed question.",
          class = "tidymarkers_capability_error"
        )
      }
      mark_scan(
        specificity_bowtie(primers, genome_index_handle(index),
          max_mismatches = max_mismatch,
          nthreads = if (is.null(BPPARAM)) 4L else 1L
        ),
        "whole-primer"
      )
    },
    blast = mark_scan(
      specificity_blast(primers, genome_index_handle(index),
        k = k, BPPARAM = BPPARAM
      ),
      if (is.na(k)) "whole-primer" else "3'-seed"
    ),
    fm = mark_scan(
      specificity_fm(primers, genome_index_handle(index),
        max_mismatches = max_mismatch, nthreads = 1L
      ),
      "whole-primer"
    )
  )
  # every row names the build it answers for; `NA` is an unnamed single genome
  out$build <- index_result_build(index)
  attr(out, "genome_build") <- index_result_build(index)
  out
}

# Tag a specificity table with the question it answers.  The attribute is the
# wire format (`attr(x, "scan")`); the class is the type, so a caller can
# require it; `scan_question()` is the only sanctioned reader.
mark_scan <- function(x, question) {
  stopifnot(question %in% c("3'-seed", "whole-primer"))
  attr(x, "scan") <- question
  class(x) <- c("specificity_scan", class(x))
  x
}

#' The question a specificity table answers
#'
#' Every genome-wide scan answers one of two questions: `"3'-seed"` (a seed
#' match anywhere in the genome) or `"whole-primer"` (the whole oligo binds
#' somewhere else).  The counts are not comparable across the two, so the
#' question travels with the table as `attr(x, "scan")` and this is the only
#' sanctioned reader.  A named vector is returned for a pangenome scan, one
#' entry per genome.
#'
#' @param x a specificity table (as returned by [specificity_of()]).
#' @return `"3'-seed"` or `"whole-primer"`; a named vector for a pangenome.
#'   An absent or unknown value is an error, never `NULL`.
#' @examples
#' specificity_of(
#'   c(AS = "ACGTACGTACGTACGTACGT"),
#'   Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGTACGTACGTACGT"))
#' ) |>
#'   scan_question()
#' @family specificity
#' @export
scan_question <- function(x) {
  q <- attr(x, "scan")
  if (is.null(q) || !length(q)) {
    abort_msg(
      "This table carries no `scan` question.",
      i = "Only a result of specificity_of() names the question its counts answer.",
      class = "tidymarkers_scan_error"
    )
  }
  bad <- unname(q)[is.na(match(unname(q), c("3'-seed", "whole-primer")))]
  if (length(bad)) {
    abort_msg(
      "Unknown scan question{?s} {.val {bad}}.",
      class = "tidymarkers_scan_error"
    )
  }
  q
}

#' Build a `marker_genome` from a `BSgenome` object
#'
#' `getSeq()` on a `BSgenome` is the Bioconductor-native way to read reference
#' sequence, so a `marker_genome` can be backed by one instead of a FASTA file.
#' Such a genome has no FASTA, so it is built with an attached index
#' (`blast_db=`/`bowtie_index=`/`fm_index=`) rather than by indexing it here.
#' @keywords internal
#' @noRd
bsgenome_kasp_genome <- function(bsgenome, chrom_map, blast_db, fm_index = NULL,
                                 bowtie_index = NULL, build = NULL) {
  load_namespace_quietly("BSgenome", reason = gettext("for a BSgenome genome source"))
  if (is.character(bsgenome)) {
    load_namespace_quietly(bsgenome, reason = gettext("for the BSgenome genome source"))
    bsgenome <- getExportedValue(bsgenome, bsgenome)
  }
  si <- GenomeInfoDb::seqinfo(bsgenome) # names + lengths, no DNA
  raw <- as.character(GenomeInfoDb::seqnames(si))
  len <- rlang::set_names(as.integer(GenomeInfoDb::seqlengths(si)), raw)
  norm <- if (!is.null(chrom_map)) unname(chrom_map[raw]) else raw
  keep <- if (!is.null(chrom_map)) !is.na(norm) else rep(TRUE, length(raw))
  seqnames <- rlang::set_names(raw[keep], (if (!is.null(chrom_map)) norm else raw)[keep])
  structure(
    list(
      path = NULL, fa = NULL, bsgenome = bsgenome, seqnames = seqnames,
      lengths = rlang::set_names(unname(len[raw[keep]]), names(seqnames)),
      blast_db = blast_db, fm_index = fm_index, bowtie_index = bowtie_index,
      build = build
    ),
    class = "marker_genome"
  )
}

# One release id per genome. A pangenome names each component; an unnamed
# build is `NA` and can still be queried alone, never pooled with another.
check_genome_build <- function(build) {
  if (is.null(build)) {
    return(NULL)
  }
  if (!is.character(build) || length(build) != 1L || is.na(build) || !nzchar(build)) {
    abort_msg("`build` must be one short release id, e.g. \"A\" or \"TAIR10\".")
  }
  if (grepl("[/\\\\]", build) || build %in% c(".", "..")) {
    abort_msg("`build` must be one short release id, e.g. \"A\" or \"TAIR10\".")
  }
  build
}

# `genome_signature()` keys a BSgenome by names and lengths only (its class is
# always "BSgenome", so even the package is absent), so two same-sized genomes
# would share one FASTA.  The FASTA *is* whole-genome DNA, so key it by the DNA.
# That costs one `getSeq()` pass (measured 2.3 s cold / 1.5 s warm on TAIR9,
# 120 Mbp, so ~40 s on a 3 Gbp genome, once per index build); a cheaper base probe
# is no cheaper, because BSgenome decompresses a whole per-chromosome rda for any
# subsequence.  One chromosome at a time caps peak memory.
bsgenome_fasta_key <- function(bsgenome) {
  si <- BSgenome::seqinfo(bsgenome)
  parts <- vapply(as.character(GenomeInfoDb::seqnames(si)), function(n) {
    # hash the character, not the DNAString: hashing the latter is unstable
    paste0(n, ":", rlang::hash(as.character(BSgenome::getSeq(bsgenome, n))))
  }, character(1))
  substr(rlang::hash(paste(parts, collapse = "|")), 1L, 12L)
}

# an index build needs a FASTA file; a bare BSgenome is written out once (№8, №10)
genome_fasta_path <- function(genome, call) {
  # `getSeq()` names the DNAStringSet, so the FASTA carries the genome's names
  if (inherits(genome, "BSgenome")) {
    path <- fs::path(fs::path_temp(), paste0("genome-", bsgenome_fasta_key(genome), ".fa"))
    # the name is content-keyed, and a link at it is not a cache hit (ADR 0033)
    if (!cache_present(path)) {
      declared_path(path)
      Biostrings::writeXStringSet(BSgenome::getSeq(genome), path)
    }
    return(path)
  }
  # a BSgenome-backed `marker_genome` keeps no path, so it still says why
  if (inherits(genome, "marker_genome") && is.null(genome$path)) {
    abort_msg(
      "{.fn {call}} needs a DNAStringSet or a FASTA-backed genome, and this genome is BSgenome-backed.",
      i = "Pass the {.cls BSgenome} itself instead of {.code marker_genome(bsgenome = )}.",
      i = "Or attach a prebuilt index: {.code marker_genome(bsgenome = , {index_arg(call)} = <path>)}.",
      i = "{.fn specificity_of} uses that index as-is."
    )
  }
  genome$path
}

# the `marker_genome()` argument that attaches each builder's index
index_arg <- function(call) {
  switch(call,
    marker_blast_index = "blast_db",
    marker_bowtie_index = "bowtie_index",
    marker_fm_index = "fm_index",
    abort_msg("Unknown index call: {call}")
  )
}

# the builders branch on `inherits(genome, "marker_genome")` to use the path
# returned here, so a materialised BSgenome is handed on as one
genome_fasta_source <- function(genome, call) {
  path <- genome_fasta_path(genome, call)
  if (inherits(genome, "BSgenome")) {
    genome <- bsgenome_kasp_genome(genome, NULL, NULL)
    genome$path <- path
  }
  list(genome = genome, fasta = path)
}

#' Print a KASP genome source
#'
#' Compact summary: FASTA file, sequences and whether a BLAST db is attached.
#'
#' @param x a `marker_genome`.
#' @param ... unused.
#' @return `x`, invisibly (called for its side effect of printing a summary).
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' print(marker_genome(fa, chrom_map = c("1" = "Chr1")))
#' @method print marker_genome
#' @family input
#' @export
print.marker_genome <- function(x, ...) {
  source <- if (!is.null(x$bsgenome)) paste0("<BSgenome> ", class(x$bsgenome)[1L]) else fs::path_file(x$path)
  cat("<marker_genome>", source, "\n")
  cat(
    "  sequences:", length(x$seqnames),
    paste0("(", paste(utils::head(names(x$seqnames), 6L), collapse = ", "), ")"), "\n"
  )
  cat("  BLAST db:", if (is.null(x$blast_db)) "none (see marker_blast_index())" else x$blast_db, "\n")
  cat("  bowtie index:", if (is.null(x$bowtie_index)) "none (see marker_bowtie_index())" else x$bowtie_index, "\n")
  cat("  FM index:", if (is.null(x$fm_index)) "none (see marker_fm_index())" else x$fm_index, "\n")
  invisible(x)
}

#' Sequence names of a genome (indexed source or DNAStringSet)
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @return a character vector of sequence names.
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' genome_chroms(marker_genome(fa, chrom_map = c("1" = "Chr1")))
#' @family input
#' @export
genome_chroms <- function(genome) {
  if (inherits(genome, "marker_genome")) names(genome$seqnames) else names(genome)
}

#' Length of one sequence of a genome
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param chrom sequence name.
#' @return the sequence length as a single integer.
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' genome_length(marker_genome(fa, chrom_map = c("1" = "Chr1")), "Chr1")
#' @family input
#' @export
genome_length <- function(genome, chrom) {
  if (!inherits(genome, "marker_genome")) {
    return(length(genome[[chrom]]))
  }
  len <- genome$lengths[[chrom]]
  if (is.null(len)) abort_msg("Unknown sequence: {chrom}")
  len
}

#' Read a region as a DNAString, from an indexed source or a DNAStringSet
#'
#' This is the single accessor the design code uses, so an indexed genome and an
#' in-memory one behave identically.
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param chrom sequence name (as returned by [genome_chroms()]).
#' @param start,end 1-based inclusive coordinates; default to the whole sequence.
#' @return a `DNAString` for the requested interval.
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' genome_seq(marker_genome(fa, chrom_map = c("1" = "Chr1")), "Chr1", 2, 6)
#' @family input
#' @export
genome_seq <- function(genome, chrom, start = 1L, end = genome_length(genome, chrom)) {
  start <- as.integer(start)
  end <- as.integer(end)
  if (!inherits(genome, "marker_genome")) {
    # the whole sequence needs no copy: `subseq()` of a full range memcpys the
    # chromosome, and a genome-scale scan asks for it once per primer pair
    if (start == 1L && end == length(genome[[chrom]])) {
      return(genome[[chrom]])
    }
    return(Biostrings::subseq(genome[[chrom]], start = start, end = end))
  }
  raw <- genome$seqnames[[chrom]]
  if (is.null(raw)) abort_msg("Unknown sequence: {chrom}")
  if (!is.null(genome$bsgenome)) {
    gr <- GenomicRanges::GRanges(raw, IRanges::IRanges(start, end))
    return(BSgenome::getSeq(genome$bsgenome, gr)[[1L]])
  }
  # A genome-wide step asks for the whole chromosome once per primer pair, so
  # the read (not the scan) repeats; keep the first read per sequence.
  whole <- start == 1L && end == as.integer(genome$lengths[[chrom]])
  cache <- genome$cache
  if (whole && !is.null(cache)) {
    cached <- cache[[chrom]]
    if (!is.null(cached)) {
      return(cached)
    }
  }
  gr <- GenomicRanges::GRanges(raw, IRanges::IRanges(start, end))
  seq <- Rsamtools::scanFa(genome$fa, gr)[[1L]]
  if (whole && !is.null(cache)) {
    cache[[chrom]] <- seq
  }
  seq
}
