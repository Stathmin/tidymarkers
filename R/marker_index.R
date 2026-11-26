# =============================================================================
#  marker_index.R  --  the manifest and shard layout behind a genome_index
# =============================================================================
#  ADR 0005 makes the genome a declared strategy.  The strategy announces its
#  capabilities (marker_genome.R) and this file records where its shards live:
#  one row per shard with the tool, version, source and content hash, so a
#  stale or missing shard is detectable without opening a sequence handle.
#
#  A file-backed index is never loaded into RAM: the object stores paths only,
#  and a per-chromosome build reads one chromosome at a time (see
#  `index_chrom_seq()`), writes it as a one-sequence FASTA, indexes it and
#  removes the copy.  A monolithic build is the fallback while the estimate fits
#  `BOWTIE_MONOLITH_LIMIT`; above it a bowtie genome is always sharded.
# =============================================================================

# the ADR 0005 manifest columns, in order
INDEX_MANIFEST_COLUMNS <- c(
  "shard", "chrom", "strategy", "path", "bytes", "sha256", "tool",
  "tool_version", "build_id", "genome_build", "source_sha256", "created"
)

# 1.34 bytes per genome base measured on TAIR10 (see marker_bowtie.R).  Wheat
# (16 Gb) is ~21.5 GB; a host that cannot hold that build gets per-chromosome
# shards by default.
BOWTIE_MONOLITH_LIMIT <- 21.5 * 1024^3

# The tool and its version for a strategy's manifest rows.
index_pkg_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

index_tool <- function(strategy) {
  switch(strategy,
    memory = c("Biostrings", index_pkg_version("Biostrings")),
    bowtie = c("Rbowtie", index_pkg_version("Rbowtie")),
    blast = c("BLAST+", marker_engine_versions()$version[
      match("BLAST+", marker_engine_versions()$engine)
    ]),
    fm = c("Rsubread", index_pkg_version("Rsubread"))
  )
}

# The identity of the *source*: a path for a FASTA genome, a sequence hash for
# an in-memory one.  Cheap on purpose - the full file hash is filled by
# `index_manifest()` / `index_verify()`, not when the genome is opened.
index_source_id <- function(genome) {
  if (inherits(genome, "DNAStringSet")) {
    paste0("DNAStringSet:", hash_object(paste(as.character(genome), collapse = "|")))
  } else if (inherits(genome, "marker_genome")) {
    genome$path %||% "<BSgenome>"
  } else {
    "<BSgenome>"
  }
}

index_source_sha256 <- function(genome) {
  if (inherits(genome, "DNAStringSet")) {
    return(hash_object(paste(as.character(genome), collapse = "|")))
  }
  path <- if (inherits(genome, "marker_genome")) genome$path else NULL
  if (!is.null(path) && fs::file_exists(path)) {
    return(hash_file(path))
  }
  NA_character_
}

index_genome_build <- function(genome) {
  if (inherits(genome, "marker_genome")) genome$build %||% NA_character_ else NA_character_
}

# an attached legacy index slot, if this genome carries one
index_attached <- function(genome, slot) {
  if (inherits(genome, "marker_genome")) genome[[slot]] else NULL
}

# Planned shards: one row per index the strategy would write.  No file is read
# and no handle is opened; `bytes`, `sha256`, `source_sha256` and `created` stay
# NA until `index_manifest()` refreshes them from disk.
index_shards <- function(genome, strategy, chroms, dir) {
  dir <- dir %||% fs::path_temp()
  build <- index_genome_build(genome)
  id <- hash_object(list(index_source_id(genome), strategy, chroms))
  id <- substr(id, 1L, 12L)
  rows <- switch(strategy,
    memory = data.frame(
      shard = "memory", chrom = NA_character_, strategy = "memory",
      path = index_source_id(genome), stringsAsFactors = FALSE
    ),
    bowtie = if (is.null(chroms)) {
      data.frame(
        shard = "bowtie", chrom = NA_character_, strategy = "bowtie",
        path = index_attached(genome, "bowtie_index") %||%
          as.character(fs::path(dir, "genome_index_bowtie")),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        shard = paste0("bowtie_", chroms), chrom = chroms, strategy = "bowtie",
        path = as.character(fs::path(dir, paste0("genome_index_bowtie_", chroms))),
        stringsAsFactors = FALSE
      )
    },
    blast = if (is.null(chroms)) {
      data.frame(
        shard = "blast", chrom = NA_character_, strategy = "blast",
        path = index_attached(genome, "blast_db") %||%
          index_attached(genome, "path") %||%
          as.character(fs::path(dir, "genome_index_blast", "genome.fa")),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        shard = paste0("blast_", chroms), chrom = chroms, strategy = "blast",
        path = as.character(fs::path(dir, paste0("genome_index_blast_", chroms), "genome.fa")),
        stringsAsFactors = FALSE
      )
    },
    # Rsubread has no per-chromosome build, so the FM index is one shard
    fm = data.frame(
      shard = "fm", chrom = NA_character_, strategy = "fm",
      path = as.character(fs::path(dir, "genome_index_fm")),
      stringsAsFactors = FALSE
    )
  )
  tool <- index_tool(strategy)
  tool_name <- tool[[1L]]
  tool_version <- tool[[2L]]
  out <- tibble::tibble(
    shard = rows$shard, chrom = rows$chrom, strategy = rows$strategy,
    path = rows$path, bytes = NA_real_, sha256 = NA_character_,
    tool = tool_name, tool_version = tool_version, build_id = id,
    genome_build = build, source_sha256 = NA_character_,
    created = as.POSIXct(NA)
  )
  out[, INDEX_MANIFEST_COLUMNS]
}

# The files one shard owns.  Returned sorted, so the aggregate hash is stable.
index_shard_files <- function(row) {
  if (identical(row$strategy[[1L]], "memory")) {
    return(character())
  }
  if (identical(row$strategy[[1L]], "bowtie")) {
    return(sort(Sys.glob(paste0(row$path, ".*"))))
  }
  if (identical(row$strategy[[1L]], "fm")) {
    return(sort(Sys.glob(paste0(row$path, ".*"))))
  }
  # BLAST writes `<prefix>.nhr/.nin/.nsq` beside the FASTA we pass
  sort(Sys.glob(paste0(row$path, ".*")))
}

# One shard's row, with bytes/sha/created filled from the files it owns.
index_manifest_row <- function(row, source_sha256) {
  files <- index_shard_files(row)
  row$source_sha256 <- source_sha256
  if (!length(files)) {
    if (!identical(row$strategy[[1L]], "memory")) {
      row$bytes <- NA_real_
      row$sha256 <- NA_character_
      row$created <- as.POSIXct(NA)
    }
    return(row)
  }
  # hash the artefact set, not one file: bowtie shards are six files and BLAST
  # shards three, and a mismatch in any of them is a mismatch in the shard
  parts <- vapply(files, function(f) {
    paste0(fs::path_file(f), ":", hash_file(f))
  }, character(1))
  row$bytes <- sum(as.numeric(fs::file_size(files)))
  row$sha256 <- hash_object(paste(parts, collapse = "|"))
  row$created <- as.POSIXct(max(fs::file_info(files)$modification_time))
  row
}

#' Content manifest of a genome index
#'
#' Fills the ADR 0005 manifest columns from the shards on disk: byte size, a
#' SHA-256 over the artefact set, the source hash and the creation time.  Rows
#' for a shard that has not been built keep `NA` bytes/hash.  Called by tests,
#' [index_verify()] and anyone who wants to record what a run used; opening the
#' index never reads the shards.
#' @param index a `genome_index`.
#' @return a tibble with `INDEX_MANIFEST_COLUMNS`.
#' @keywords internal
#' @noRd
index_manifest <- function(index) {
  source_sha256 <- index_source_sha256(index$genome)
  rows <- lapply(seq_len(nrow(index$manifest)), function(i) {
    index_manifest_row(index$manifest[i, , drop = FALSE], source_sha256)
  })
  out <- dplyr::bind_rows(rows)
  if (identical(index$strategy, "memory")) {
    # the memory "shard" is the sequence itself: no file, but its identity is
    # the content hash, which is what a caller records
    out$sha256 <- index_source_sha256(index$genome)
    out$bytes <- NA_real_
  }
  out[, INDEX_MANIFEST_COLUMNS]
}

# A missing shard names the chromosome; a changed one names the shard.
index_abort_missing <- function(row) {
  if (is.na(row$chrom[[1L]])) {
    abort_msg(
      "The {row$strategy} index shard {row$shard} is missing.",
      i = "Build it again with `genome_index(strategy = \"{row$strategy}\")`.",
      class = "tidymarkers_index_missing_error"
    )
  }
  abort_msg(
    "The {row$strategy} index shard {row$shard} is missing: chromosome {row$chrom}.",
    i = "Build it again with `genome_index(strategy = \"{row$strategy}\")`.",
    class = "tidymarkers_index_missing_error"
  )
}

#' Verify a genome index against its manifest
#'
#' Re-hashes every existing shard and compares it with `x$manifest`; a missing
#' shard stops with `tidymarkers_index_missing_error` naming the chromosome and
#' a changed one with `tidymarkers_index_integrity_error` naming the shard.
#' @param x a `genome_index`.
#' @return the refreshed manifest, invisibly.
#' @keywords internal
#' @noRd
index_verify <- function(x) {
  now <- index_manifest(x)
  recorded <- x$manifest
  for (i in seq_len(nrow(now))) {
    if (is.na(now$path[[i]]) || identical(now$strategy[[i]], "memory")) {
      next
    }
    if (is.na(now$sha256[[i]])) {
      index_abort_missing(now[i, , drop = FALSE])
    }
    if (!is.na(recorded$sha256[[i]]) && !identical(recorded$sha256[[i]], now$sha256[[i]])) {
      abort_msg(
        "The {now$strategy[[i]]} index shard {now$shard[[i]]} changed since the manifest was recorded.",
        i = "Delete it and rebuild, or attach the index the manifest names.",
        class = "tidymarkers_index_integrity_error"
      )
    }
  }
  invisible(now)
}

# One chromosome as a DNAString, without going through the genome cache: a
# shard build must hold one chromosome, not the whole genome, and must not
# populate the query cache the design path relies on.
index_chrom_seq <- function(genome, chrom) {
  if (inherits(genome, "DNAStringSet")) {
    return(genome[[chrom]])
  }
  if (inherits(genome, "BSgenome")) {
    return(BSgenome::getSeq(genome, chrom))
  }
  raw <- genome$seqnames[[chrom]]
  if (is.null(raw)) abort_msg("Unknown sequence: {chrom}")
  len <- as.integer(genome$lengths[[chrom]])
  gr <- GenomicRanges::GRanges(raw, IRanges::IRanges(1L, len))
  if (!is.null(genome$bsgenome)) {
    return(BSgenome::getSeq(genome$bsgenome, gr)[[1L]])
  }
  Rsamtools::scanFa(genome$fa, gr)[[1L]]
}

# Write one shard's sequence, index it, and drop the copy: peak disk is one
# chromosome plus its index, not a sliced whole genome.
index_write_chrom_fasta <- function(seq, chrom, path) {
  gs <- Biostrings::DNAStringSet(seq)
  names(gs) <- chrom
  fs::dir_create(fs::path_dir(path))
  Biostrings::writeXStringSet(gs, path)
  path
}

# Should a bowtie build be sharded although `chroms` was not given?  Only when
# the genome is large enough that one index cannot be built here.
index_should_shard <- function(genome, strategy = "bowtie") {
  if (!identical(strategy, "bowtie")) {
    return(FALSE)
  }
  bases <- if (inherits(genome, "DNAStringSet")) {
    sum(Biostrings::width(genome))
  } else if (inherits(genome, "marker_genome")) {
    sum(as.numeric(genome$lengths))
  } else {
    return(FALSE)
  }
  bases * 1.34 > BOWTIE_MONOLITH_LIMIT
}

# Build any missing shard and return the value the scan query takes: a basename
# for a monolithic bowtie index, a named set for per-chromosome shards, a BLAST
# prefix, or an FM basename.  Existing files are reused, so a second run only
# verifies presence.
index_build <- function(index) {
  rows <- index$manifest
  switch(index$strategy,
    memory = index$genome,
    bowtie = {
      out <- vapply(seq_len(nrow(rows)), function(i) {
        index_build_bowtie_shard(index$genome, rows[i, , drop = FALSE], index$threads)
      }, character(1))
      if (length(out) > 1L) {
        names(out) <- rows$chrom
        structure(out,
          class = c("marker_bowtie_index_set", "character"), chroms = rows$chrom
        )
      } else {
        out[[1L]]
      }
    },
    blast = {
      out <- vapply(seq_len(nrow(rows)), function(i) {
        index_build_blast_shard(index$genome, rows[i, , drop = FALSE])
      }, character(1))
      if (length(out) > 1L) structure(out, class = "marker_blast_index_set") else out[[1L]]
    },
    fm = index_build_fm_shard(index$genome, rows[1L, , drop = FALSE])
  )
}

index_build_bowtie_shard <- function(genome, row, threads) {
  if (fs::file_exists(paste0(row$path, ".1.ebwt"))) {
    return(row$path)
  }
  if (is.na(row$chrom)) {
    return(marker_bowtie_index(genome, basename = row$path, threads = threads))
  }
  # one chromosome at a time: peak memory is that chromosome, and the sliced
  # FASTA lives only for the build
  fa <- paste0(row$path, ".fa")
  on.exit(unlink(fa), add = TRUE)
  index_write_chrom_fasta(index_chrom_seq(genome, row$chrom), row$chrom, fa)
  bowtie_build_fasta(fa, row$path, threads)
}

index_build_blast_shard <- function(genome, row) {
  if (fs::file_exists(paste0(row$path, ".nsq"))) {
    return(row$path)
  }
  if (is.na(row$chrom) && inherits(genome, "marker_genome") && !is.null(genome$path)) {
    # index the genome's own FASTA in place, as marker_blast_index() already does
    return(marker_blast_index(genome, path = genome$path))
  }
  seq <- if (is.na(row$chrom)) {
    if (inherits(genome, "marker_genome") && !is.null(genome$bsgenome)) {
      # a BSgenome is already in RAM; materialise it once to index it
      Biostrings::DNAStringSet(BSgenome::getSeq(genome$bsgenome))
    } else {
      genome
    }
  } else {
    index_chrom_seq(genome, row$chrom)
  }
  fs::dir_create(fs::path_dir(row$path))
  db <- marker_blast_index(Biostrings::DNAStringSet(seq), path = row$path)
  # the database carries the sequence; the FASTA copy is redundant
  unlink(db)
  db
}

index_build_fm_shard <- function(genome, row) {
  if (fs::file_exists(paste0(row$path, ".files"))) {
    return(row$path)
  }
  marker_fm_index(genome, basename = row$path)
}
