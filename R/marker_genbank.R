# =============================================================================
#  marker_genbank.R  --  GenBank flat-file reader with a cached feature index
# =============================================================================
#  Roadmap Y5.1: a region (Y2) can come from GenBank, not only from VCF.  The
#  reader parses LOCUS / FEATURES / ORIGIN / // once and caches the feature rows
#  beside the file as `<path>.gbidx`; later calls only read that index, so a
#  region query never re-reads the whole flat file.  The cache is keyed on the
#  source file's size and mtime, and rebuilt when either differs.
#
#  The location parser handles plain `123..456`, `complement(...)`,
#  `join(...)`, `complement(join(...))`, `<123..>456`, single `123` and
#  `123^124` (a between-position, folded to `start = end = 123`); `start`/`end`
#  are the min/max of the joined parts.
# =============================================================================

#' Read a GenBank flat file through a cached feature index
#'
#' Parses a GenBank flat file (`LOCUS` / `FEATURES` / `ORIGIN` / `//`) into one
#' tidy row per feature.  The parse result is cached beside the file as
#' `<path>.gbidx` (an RDS); a later call only reads that index when the source
#' file's size and mtime are unchanged, so a region query never re-reads the
#' whole flat file.  A gzipped file is rejected, because indexing a compressed
#' file needs seek and this reader does not require it.
#'
#' @param path plain-text GenBank path; a `.gz` file is rejected with a message
#'   to gunzip it first.
#' @param region optional single region to keep: a `GRanges` of length one, or a
#'   data.frame / list with `chrom`, `start`, `end` (1-based, closed, as
#'   everywhere in the package).  Only features overlapping the closed interval
#'   (`start <= region_end` and `end >= region_start`) are returned; an empty
#'   result is a zero-row tibble, not an error.  `NULL` returns every feature.
#' @return a tibble, one row per feature, with columns `chrom` (the record's
#'   `LOCUS` name, falling back to `ACCESSION`), `start` (int, 1-based), `end`
#'   (int), `strand` (`"+"`/`"-"`, `"*"` for an `order()` location, i.e. an
#'   unknown orientation), `type` (the FEATURES key, e.g. `gene`, `CDS`),
#'   `gene`, `locus_tag`, `product` (from the matching qualifiers, `NA` when
#'   absent) and `qualifiers` (a list-column of all `/key=value` pairs, so
#'   nothing is lost).  The tibble carries a `file_reads` integer attribute as a
#'   diagnostic: how many times the GenBank file itself was opened/read by this
#'   call — `0` on a cache hit, `>= 1` when the index was (re)built.  Reads of
#'   the `.gbidx` file do not count.
#' @examples
#' gb <- withr::local_tempfile(fileext = ".gbk")
#' writeLines(c(
#'   "LOCUS       TESTLOCUS         120 bp    DNA     linear   UNK 01-JAN-2026",
#'   "FEATURES             Location/Qualifiers",
#'   "     source          1..120",
#'   '                     /organism="Arabidopsis thaliana"',
#'   "     CDS             join(10..30,40..60)",
#'   '                     /gene="GENE1"',
#'   '                     /product="protein one"',
#'   "ORIGIN",
#'   "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac",
#'   "//"
#' ), gb)
#' read_genbank(gb)
#' @family KASP input
#' @export
read_genbank <- function(path, region = NULL) {
  if (grepl("\\.gz$", path)) {
    cli::cli_abort(gettext("{.file {path}} is gzipped; gunzip it first"))
  }
  if (!fs::file_exists(path)) cli::cli_abort(gettext("GenBank not found: {path}"))

  reg <- if (is.null(region)) NULL else normalise_region(region)

  idx <- genbank_index(path)
  features <- idx$features

  if (!is.null(reg)) {
    features <- dplyr::filter(
      features,
      .data$chrom == reg$chrom,
      .data$start <= reg$end,
      .data$end >= reg$start
    )
  }

  attr(features, "file_reads") <- as.integer(idx$file_reads)
  features
}

#' Load (or build) the `.gbidx` feature index for a GenBank file (internal)
#' @keywords internal
#' @noRd
genbank_index <- function(path) {
  idx_path <- paste0(path, ".gbidx")
  info <- file.info(path)
  size <- as.numeric(info$size)
  mtime <- as.numeric(info$mtime)

  cached <- if (fs::file_exists(idx_path)) {
    tryCatch(readRDS(idx_path), error = function(e) NULL)
  } else {
    NULL
  }

  if (!is.null(cached) && isTRUE(cached$source_size == size) && isTRUE(cached$source_mtime == mtime)) {
    return(list(features = cached$features, file_reads = 0L))
  }

  parsed <- parse_genbank_file(path)
  idx <- list(
    source_size = size,
    source_mtime = mtime,
    records = parsed$records,
    features = parsed$features,
    origin_offset = parsed$origin_offset
  )
  # a read-only directory must not fail the read: fall back to the parsed result
  tryCatch(saveRDS(idx, idx_path), error = function(e) NULL)

  list(features = parsed$features, file_reads = parsed$file_reads)
}

#' Read a file into lines with their byte offsets (internal)
#'
#' Reads the whole file as raw bytes and splits on `\n`, so every line's offset
#' is its exact byte position (GenBank is ASCII, so a character position is a
#' byte position).  The whole-file read is only done on a cache build; a cache
#' hit never touches the source file.
#' @keywords internal
#' @noRd
read_lines_with_offsets <- function(path) {
  raw <- readBin(path, "raw", n = file.info(path)$size)
  txt <- rawToChar(raw)
  if (!nzchar(txt)) {
    return(list(lines = character(0), offsets = integer(0), file_reads = 1L))
  }
  nl <- as.integer(gregexpr("\n", txt, fixed = TRUE)[[1]])
  nl <- nl[nl > 0L]
  line_start <- c(1L, nl + 1L)
  content <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  content <- sub("\r$", "", content)
  # `strsplit` drops a trailing empty element when the file ends with a newline,
  # so `line_start` can hold one phantom offset more than `content`: trim it to
  # the actual lines so the two vectors always align.
  if (length(line_start) > length(content)) {
    line_start <- line_start[seq_along(content)]
  }
  list(lines = content, offsets = as.integer(line_start), file_reads = 1L)
}

#' Parse a GenBank file into records, features and origin offsets (internal)
#' @keywords internal
#' @noRd
parse_genbank_file <- function(path) {
  lines <- read_lines_with_offsets(path)
  if (!length(lines$lines)) {
    return(list(
      records = empty_records(),
      features = empty_features(),
      origin_offset = integer(0),
      file_reads = lines$file_reads
    ))
  }

  term <- trimws(lines$lines) == "//"
  keep <- !term
  rec_id <- cumsum(term)[keep] + 1L
  rlines <- split(lines$lines[keep], rec_id)
  roffsets <- split(lines$offsets[keep], rec_id)
  parsed <- purrr::map2(rlines, roffsets, parse_record)

  if (!length(parsed)) {
    return(list(
      records = empty_records(),
      features = empty_features(),
      origin_offset = integer(0),
      file_reads = lines$file_reads
    ))
  }

  records <- purrr::map_dfr(parsed, function(p) {
    tibble::tibble(chrom = p$chrom, accession = p$accession, length_bp = p$length_bp)
  })
  features <- purrr::map_dfr(parsed, function(p) {
    build_features_tibble(p$features, p$chrom)
  })
  origin_offset <- vapply(parsed, function(p) p$origin_offset, integer(1))

  list(
    records = records,
    features = features,
    origin_offset = origin_offset,
    file_reads = lines$file_reads
  )
}

#' Parse one GenBank record's header and feature block (internal)
#' @keywords internal
#' @noRd
parse_record <- function(rlines, roffsets) {
  chrom <- record_chrom(rlines)
  f_idx <- which(startsWith(trimws(rlines), "FEATURES"))
  o_idx <- which(trimws(rlines) == "ORIGIN")

  from <- if (length(f_idx)) f_idx[1] + 1L else length(rlines) + 1L
  to <- if (length(o_idx) && length(f_idx) && o_idx[1] > f_idx[1]) {
    o_idx[1] - 1L
  } else {
    length(rlines)
  }
  feat_lines <- if (from <= to) rlines[from:to] else character(0)

  origin_offset <- if (length(o_idx) && o_idx[1] < length(roffsets)) {
    roffsets[o_idx[1] + 1L]
  } else {
    NA_integer_
  }

  list(
    chrom = chrom,
    accession = record_accession(rlines),
    length_bp = record_length(rlines),
    origin_offset = origin_offset,
    features = parse_features_block(feat_lines)
  )
}

#' First whitespace-separated token after a fixed tag, else `NA` (internal)
#' @keywords internal
#' @noRd
first_token_of_line <- function(rlines, tag) {
  idx <- which(startsWith(trimws(rlines), tag))
  if (!length(idx)) {
    return(NA_character_)
  }
  m <- stringr::str_match(trimws(rlines[idx[1]]), paste0("^", tag, "\\s+(\\S+)"))
  m[1, 2]
}

#' Record chromosome name: LOCUS name, falling back to ACCESSION (internal)
#' @keywords internal
#' @noRd
record_chrom <- function(rlines) {
  locus <- first_token_of_line(rlines, "LOCUS")
  if (!is.na(locus)) {
    return(locus)
  }
  first_token_of_line(rlines, "ACCESSION")
}

#' Record accession (internal)
#' @keywords internal
#' @noRd
record_accession <- function(rlines) first_token_of_line(rlines, "ACCESSION")

#' Record length in bp from the LOCUS line (internal)
#' @keywords internal
#' @noRd
record_length <- function(rlines) {
  idx <- which(startsWith(trimws(rlines), "LOCUS"))
  if (!length(idx)) {
    return(NA_integer_)
  }
  m <- stringr::str_match(trimws(rlines[idx[1]]), "\\s(\\d+)\\s+bp\\b")
  as.integer(m[1, 2])
}

#' Parse a FEATURES block into a list of `type`/`location`/qualifier records (internal)
#' @keywords internal
#' @noRd
parse_features_block <- function(block) {
  features <- list()
  current <- NULL
  for (ln in block) {
    t <- trimws(ln)
    if (!nzchar(t)) next
    if (startsWith(t, "/")) {
      if (!is.null(current)) {
        q <- parse_qualifier_line(t)
        current$qual_keys <- c(current$qual_keys, q$key)
        current$qual_vals <- c(current$qual_vals, q$value)
      }
    } else if (is_qualifier_continuation(ln)) {
      # a wrapped qualifier value continues on a line whose key column is blank;
      # before the first qualifier exists, a blank key column is a wrapped
      # *location* line instead, so it is appended to the location (not dropped)
      if (is.null(current)) {
        next
      } else if (length(current$qual_vals)) {
        k <- length(current$qual_vals)
        current$qual_vals[k] <- paste(current$qual_vals[k], trimws(ln), sep = " ")
      } else {
        current$location <- paste(current$location, trimws(ln), sep = " ")
      }
    } else {
      # a new feature line flushes the previous feature (with its qualifiers),
      # so the qualifiers are never attached to a stale copy of `current`
      if (!is.null(current)) features <- c(features, list(current))
      m <- stringr::str_match(ln, "^\\s*(\\S+)\\s+(\\S.*)$")
      if (is.na(m[1, 1])) {
        current <- list(
          type = t, location = "",
          qual_keys = character(0), qual_vals = character(0)
        )
      } else {
        current <- list(
          type = m[1, 2], location = trimws(m[1, 3]),
          qual_keys = character(0), qual_vals = character(0)
        )
      }
    }
  }
  if (!is.null(current)) features <- c(features, list(current))
  features
}

#' Is a FEATURES line a continuation of the previous qualifier? (internal)
#' @keywords internal
#' @noRd
is_qualifier_continuation <- function(ln) {
  key_col <- substr(ln, 1L, 21L)
  !nzchar(trimws(key_col))
}

#' Split a `/key=value` qualifier into key and (unquoted) value (internal)
#' @keywords internal
#' @noRd
parse_qualifier_line <- function(ln) {
  s <- sub("^/", "", trimws(ln))
  m <- stringr::str_match(s, "^([^=]*)=(.*)$")
  if (is.na(m[1, 1])) {
    return(list(key = s, value = ""))
  }
  list(key = m[1, 2], value = gsub('^"|"$', "", m[1, 3]))
}

#' Build the feature tibble for one record (internal)
#' @keywords internal
#' @noRd
build_features_tibble <- function(features, chrom) {
  if (!length(features)) {
    return(empty_features())
  }
  rows <- purrr::map(features, function(f) {
    loc <- parse_location(f$location)
    q <- if (length(f$qual_keys)) {
      as.list(f$qual_vals) |> rlang::set_names(f$qual_keys)
    } else {
      list()
    }
    tibble::tibble(
      chrom = chrom,
      start = as.integer(loc$start),
      end = as.integer(loc$end),
      strand = loc$strand,
      type = f$type,
      gene = qual_value(q, "gene"),
      locus_tag = qual_value(q, "locus_tag"),
      product = qual_value(q, "product"),
      qualifiers = list(q)
    )
  })
  dplyr::bind_rows(rows)
}

#' Read a qualifier from a named list, `NA` when absent (internal)
#' @keywords internal
#' @noRd
qual_value <- function(q, key) {
  v <- q[[key]]
  if (is.null(v)) NA_character_ else as.character(v)[1]
}

#' Parse a feature location into `start`/`end`/`strand` (internal)
#' @keywords internal
#' @noRd
parse_location <- function(location) {
  loc <- trimws(location)
  outer_complement <- grepl("^complement\\(", loc)
  if (outer_complement) {
    loc <- sub("^complement\\((.+)\\)$", "\\1", loc)
  }

  is_order <- startsWith(loc, "order(")
  is_join <- startsWith(loc, "join(")
  if (is_order || is_join) {
    loc <- sub("^(join|order)\\((.+)\\)$", "\\2", loc)
  }

  parts <- strsplit(loc, ",", fixed = TRUE)[[1]]
  starts <- integer(0)
  ends <- integer(0)
  part_strands <- character(0)

  for (p in parts) {
    p <- trimws(p)
    part_strand <- "+"
    if (grepl("^complement\\(", p)) {
      part_strand <- "-"
      p <- sub("^complement\\((.+)\\)$", "\\1", p)
    }
    p <- gsub("[<>]", "", trimws(p))
    if (grepl("\\^", p)) {
      left <- suppressWarnings(as.integer(sub("\\^.*$", "", p)))
      if (!is.na(left)) {
        starts <- c(starts, left)
        ends <- c(ends, left)
        part_strands <- c(part_strands, part_strand)
      }
    } else if (grepl("\\.\\.", p)) {
      bounds <- strsplit(p, "..", fixed = TRUE)[[1]]
      s <- suppressWarnings(as.integer(bounds[1]))
      e <- suppressWarnings(as.integer(bounds[2]))
      if (!is.na(s)) {
        starts <- c(starts, s)
        ends <- c(ends, e)
        part_strands <- c(part_strands, part_strand)
      }
    } else {
      v <- suppressWarnings(as.integer(p))
      if (!is.na(v)) {
        starts <- c(starts, v)
        ends <- c(ends, v)
        part_strands <- c(part_strands, part_strand)
      }
    }
  }

  strand <- "+"
  if (is_order) {
    strand <- "*"
  } else if (is_join) {
    if (length(part_strands)) {
      strand <- if (length(unique(part_strands)) > 1L) "*" else unique(part_strands)[1]
    }
    if (outer_complement) {
      strand <- switch(strand,
        "+" = "-",
        "-" = "+",
        "*" = "*"
      )
    }
  } else {
    strand <- if (outer_complement) "-" else "+"
  }

  list(
    start = if (length(starts)) min(starts) else NA_integer_,
    end = if (length(ends)) max(ends) else NA_integer_,
    strand = strand
  )
}

#' Empty records tibble with the right schema (internal)
#' @keywords internal
#' @noRd
empty_records <- function() {
  tibble::tibble(
    chrom = character(0),
    accession = character(0),
    length_bp = integer(0)
  )
}

#' Empty features tibble with the right schema (internal)
#' @keywords internal
#' @noRd
empty_features <- function() {
  tibble::tibble(
    chrom = character(0),
    start = integer(0),
    end = integer(0),
    strand = character(0),
    type = character(0),
    gene = character(0),
    locus_tag = character(0),
    product = character(0),
    qualifiers = list()
  )
}
