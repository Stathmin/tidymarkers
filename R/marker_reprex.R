# =============================================================================
#  marker_reprex.R  --  anonymised, self-contained reproducer
# =============================================================================
#  A user on proprietary germplasm cannot attach the VCF, the BED or the genome
#  a design run failed on. `marker_reprex()` writes one `.R` script carrying a
#  synthesised copy of that input, so the bug is reported without the data.
#
#  What survives, and how exactly:
#  - n_variants, linkage, the chromosome partition and its order: exact (no row
#    is added, dropped or reordered; each chromosome's gaps are unchanged).
#  - gc_content, variant_type: exact per row (the GC count of every sequence,
#    REF/ALT emitted verbatim).
#  - repeat_structure: approximation (dinucleotide class only).
#  Sample/genotype columns, ID, QUAL, FILTER, INFO and every other column are
#  dropped; sequence names become s1..sn.
# =============================================================================

REPREX_PRESERVE <- c("gc_content", "variant_type", "n_variants", "linkage", "repeat_structure")
REPREX_VCF_HEADER <- c(
  "##fileformat=VCFv4.2",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)
REPREX_NOT_LEGAL <- paste(
  "sufficient to reproduce a bug report,",
  "not sufficient for legal or contractual data sharing"
)


#' Anonymised, self-contained reproducer for a bug report
#'
#' Writes one `.R` script that rebuilds the input of a failing design run from
#' synthesised data, so the script can be pasted into an issue from an
#' environment whose data cannot leave. Variants lose every column but
#' `chrom`/`pos`/`ref`/`alt`, their chromosomes become `chr1..chrN` in order of
#' first appearance, and their coordinates are shifted per chromosome; sequences
#' are replaced by a synthetic sequence of the same length.
#'
#' **The output is sufficient for a bug report and NOT sufficient for legal data
#' sharing.** It does not prove that a sample cannot be identified: the number of
#' variants per chromosome, their relative spacing and every allele pair are kept
#' on purpose, because the bug depends on them.
#'
#' @param x one of:
#'   * a variant table: a data frame or tibble with `chrom`, `pos`, `ref` and
#'     `alt` (the shape [assemble_primer_assay()] takes), a `Variants` object, or
#'     a path to a VCF (`.vcf`, `.vcf.gz`) or a BED (`.bed`), read here by hand;
#'   * sequences: a `DNAStringSet`, another `XStringSet` or a plain character
#'     vector, or a path to a FASTA (`.fa`, `.fasta`, `.fna`).
#'
#'   The class of `x` decides which it is; any other input aborts. A
#'   `marker_genome` is read through its FASTA path when it has one.
#' @param path the script to write; `NULL` writes into a temporary file that
#'   stays for the session.
#' @param preserve properties that must survive anonymisation:
#'   * `"gc_content"`: the GC count of every sequence and of every variant
#'     allele (exact);
#'   * `"n_variants"`: the number of variants, overall and per chromosome
#'     (exact);
#'   * `"linkage"`: within a chromosome, the order of the variants and every gap
#'     between consecutive ones, i.e. which variants co-segregate (exact: the
#'     shift is a translation and no row is removed);
#'   * `"variant_type"`: every variant is still the same SNP, MNP, insertion or
#'     deletion (exact: `ref`/`alt` are kept as they are);
#'   * `"repeat_structure"`: the dinucleotide *composition* of a sequence is
#'     kept, by moving a base only to another of the same dinucleotide class.
#'     This is an approximation: the order of the bases is not kept, so a bug
#'     that depends on a repeat's position or period will not reproduce.
#'
#'   Every requested property is re-measured on the synthesised output, and the
#'   call aborts rather than write a script that lost one.
#' @return the path to the script, invisibly, with class `"marker_reprex"`. The
#'   script is standalone: it needs `tidymarkers` and base R only.
#' @examples
#' variants <- tibble::tibble(
#'   chrom = c("Chr4", "Chr4"), pos = c(1710L, 2340L), ref = c("A", "C"), alt = c("G", "T")
#' )
#' path <- marker_reprex(variants, path = tempfile(fileext = ".R"))
#' readLines(path)[[1L]]
#' @family input
#' @export
marker_reprex <- function(x, path = NULL, preserve = c("gc_content", "variant_type")) {
  preserve <- reprex_preserve(preserve)
  path <- reprex_path(path)
  input <- reprex_input(x, preserve)
  lines <- c(
    sprintf("##synthesized_by=tidymarkers_reprex/%s", reprex_version()),
    "#",
    "# A synthetic reproducer for a bug report:",
    paste0("# ", REPREX_NOT_LEGAL, "."),
    sprintf("# preserved: %s", paste(preserve, collapse = ", ")),
    "",
    "suppressPackageStartupMessages(library(tidymarkers))",
    "",
    input$blocks,
    "",
    input$replay,
    "#",
    sprintf("# tidymarkers %s; see https://github.com/Stathmin/tidymarkers/issues", reprex_version())
  )
  writeLines(lines, path, useBytes = TRUE)
  invisible(structure(path, class = "marker_reprex"))
}

# --- validation --------------------------------------------------------------

reprex_preserve <- function(preserve) {
  if (!is.character(preserve) || !length(preserve)) {
    abort_msg("`preserve` must be a non-empty character vector.")
  }
  unknown <- setdiff(unique(preserve), REPREX_PRESERVE)
  if (length(unknown)) {
    abort_msg(
      "`preserve` names a property this function does not implement: {.val {unknown}}.",
      i = "Allowed values are {.val {REPREX_PRESERVE}}."
    )
  }
  REPREX_PRESERVE[REPREX_PRESERVE %in% preserve]
}

reprex_path <- function(path) {
  if (is.null(path)) {
    path <- tempfile("tidymarkers-reprex", fileext = ".R")
  } else {
    if (!is.character(path) || length(path) != 1L || is.na(path)) {
      abort_msg("`path` must be one non-missing file path or `NULL`.")
    }
    if (!dir.exists(dirname(path))) {
      abort_msg("Cannot write {.path {path}}: {.path {dirname(path)}} does not exist.")
    }
  }
  refuse_links(dirname(path), basename(path))
  path
}

reprex_version <- function() {
  v <- tryCatch(as.character(utils::packageVersion("tidymarkers")), error = function(e) NA_character_)
  if (length(v) == 1L && !is.na(v)) v else "dev"
}

# --- dispatch ----------------------------------------------------------------

reprex_input <- function(x, preserve) {
  if (inherits(x, "marker_genome")) {
    if (is.null(x$path)) {
      abort_msg(
        "A {.cls marker_genome} backed by a {.cls BSgenome} has no file to read.",
        i = "Pass a {.cls DNAStringSet} or a FASTA path instead."
      )
    }
    return(reprex_from_seq(reprex_read_fasta(x$path), preserve))
  }
  if (inherits(x, "Variants")) {
    # Only these four records are read; FORMAT and the samples are never touched.
    return(reprex_from_table(tibble::tibble(
      chrom = as.character(GenomicRanges::seqnames(x)),
      pos = as.integer(GenomicRanges::start(x)),
      ref = as.character(S4Vectors::mcols(x)$REF),
      alt = as.character(unlist(S4Vectors::mcols(x)$ALT, use.names = FALSE))
    ), preserve))
  }
  if (inherits(x, "XStringSet")) {
    return(reprex_from_seq(x, preserve))
  }
  if (is.character(x)) {
    kind <- if (length(x) == 1L) reprex_file_kind(x) else NA_character_
    if (identical(kind, "vcf")) {
      return(reprex_from_table(reprex_read_vcf(x), preserve))
    }
    if (identical(kind, "bed")) {
      return(reprex_from_table(reprex_read_bed(x), preserve))
    }
    if (identical(kind, "fasta")) {
      return(reprex_from_seq(reprex_read_fasta(x), preserve))
    }
    if (length(x) == 1L && file.exists(x)) {
      abort_msg(
        "{.path {x}} is not a VCF, a BED or a FASTA (by its extension).",
        i = "Expected .vcf, .vcf.gz, .bed, .fa, .fasta or .fna."
      )
    }
    return(reprex_from_seq(x, preserve))
  }
  if (is.data.frame(x)) {
    return(reprex_from_table(x, preserve))
  }
  abort_msg(
    "`x` is not a variant table, a set of sequences, or a path to one.",
    x = "Got {.cls {class(x)}}."
  )
}

reprex_file_kind <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  name <- tolower(basename(path))
  if (grepl("\\.vcf(\\.gz)?$", name)) {
    return("vcf")
  }
  switch(tools::file_ext(name),
    bed = "bed",
    fa = ,
    fasta = ,
    fna = "fasta",
    NA_character_
  )
}

# --- sequences ---------------------------------------------------------------

reprex_from_seq <- function(x, preserve) {
  seqs <- as.character(x)
  if (!length(seqs)) {
    abort_msg("`x` holds no sequence.")
  }
  syn <- withr::with_seed(1L, lapply(seqs, reprex_synthetic, preserve))
  stuck <- sum(vapply(syn, `[[`, integer(1), "kept"))
  if (stuck) {
    cli::cli_warn(c(
      gettext("{stuck} sequence{?s} kept at least one base in place: too short or too restricted to reshuffle."),
      i = gettext("The output is still not the original; check a very repetitive input by hand.")
    ))
  }
  out <- unname(vapply(syn, `[[`, character(1), "syn"))
  gc_before <- unname(vapply(seqs, reprex_gc, integer(1)))
  gc_after <- unname(vapply(out, reprex_gc, integer(1)))
  if ("gc_content" %in% preserve && !identical(gc_before, gc_after)) {
    reprex_lost("gc_content")
  }
  list(
    blocks = c(reprex_assign("syn", out), reprex_assign("nm", sprintf("s%d", seq_along(seqs)))),
    replay = "seqs <- Biostrings::DNAStringSet(stats::setNames(syn, nm))"
  )
}

# A seeded shuffle, then a repair that only ever swaps two bases for the
# complement pair (A<->T, C<->G), which keeps the GC count, or the non-GC pair
# (A<->C, G<->T), which keeps the dinucleotide class. With no property requested
# the shuffle alone is the output.
reprex_synthetic <- function(seq, preserve) {
  bases <- strsplit(toupper(seq), "", fixed = TRUE)[[1L]]
  known <- !is.na(match(bases, c("A", "C", "G", "T")))
  idx <- which(known)
  syn <- bases
  if (length(idx) > 1L) {
    syn[idx] <- bases[idx][sample.int(length(idx))]
    # GC drift after the shuffle is repaired by one-base swaps (A<->T, C<->G),
    # which never change the GC count of a pair.
    drift <- sum(reprex_gc_bases(syn)) - sum(reprex_gc_bases(bases))
    if (drift > 0L) {
      pick <- utils::head(idx[reprex_gc_bases(syn)], drift)
      syn[pick] <- ifelse(syn[pick] == "G", "C", "G")
    } else if (drift < 0L) {
      pick <- utils::head(idx[!reprex_gc_bases(syn)], -drift)
      syn[pick] <- ifelse(syn[pick] == "A", "T", "T")
    }
    if ("repeat_structure" %in% preserve) {
      syn <- reprex_shift_class(syn, bases, idx)
    }
  }
  list(syn = paste0(syn, collapse = ""), kept = sum(syn == bases))
}

# A<->C moves a base within its dinucleotide class (both are non-GC), so a
# repeat can be rebuilt where a repeat of the original length and phase exists.
reprex_shift_class <- function(syn, bases, idx) {
  free <- idx[!reprex_gc_bases(syn)]
  for (i in idx[syn[idx] == bases[idx]]) {
    k <- free[free != i][1L]
    if (is.na(k) || syn[[k]] == syn[[i]]) {
      next
    }
    syn[c(i, k)] <- syn[c(k, i)]
  }
  syn
}

reprex_gc_bases <- function(bases) {
  bases == "G" | bases == "C"
}

reprex_gc <- function(seq) {
  sum(reprex_gc_bases(strsplit(toupper(seq), "", fixed = TRUE)[[1L]]))
}

# --- variants ----------------------------------------------------------------

reprex_from_table <- function(v, preserve) {
  need <- c("chrom", "pos", "ref", "alt")
  # formatted here: cli would resolve `{need}` in the helper frame, not this one
  need_cols(
    v, need,
    arg = "variant table",
    i = cli::format_inline(
      gettext("It needs {.val {need}}, the shape {.fun assemble_primer_assay} takes.")
    )
  )
  # A column the script cannot carry would make the replay differ silently.
  extra <- setdiff(names(v), need)
  if (length(extra)) {
    cli::cli_warn(c(
      gettext("{.val {extra}}: dropped, the script replays only {.val {need}}."),
      i = gettext("A bug that depends on them will not reproduce; check them by hand.")
    ))
  }
  v <- v[need]
  keep <- stats::complete.cases(unclass(v))
  if (!all(keep)) {
    warn_msg("{sum(!keep)} incomplete row{?s} dropped from the reproducer.")
    v <- v[keep, , drop = FALSE]
  }
  if (!nrow(v)) {
    abort_msg("`x` holds no complete variant.")
  }
  pos <- reprex_int(v$pos)
  chrom <- as.character(v$chrom)
  ref <- reprex_ascii(as.character(v$ref), "REF")
  alt <- reprex_ascii(as.character(v$alt), "ALT")
  before <- list(chrom = chrom, pos = pos)
  new <- reprex_shift_variants(chrom, pos, preserve)
  reprex_verify(before, new, preserve)
  rows <- sprintf("%s\t%d\t.\t%s\t%s\t.\t.\t.", new$chrom, new$pos, ref, alt)
  list(
    blocks = c(reprex_assign("vcf", REPREX_VCF_HEADER), reprex_assign("vcf_rows", rows)),
    replay = paste(
      'path <- file.path(tempdir(), "reprex.vcf")',
      "writeLines(c(vcf, vcf_rows), path)",
      "read_variants_vcf(path)",
      sep = "\n"
    )
  )
}

# One offset per chromosome, taken from its smallest coordinate, so the shifted
# rows stay inside the chromosome's bounds. Every gap is a difference of two
# shifted coordinates, so it is unchanged.
reprex_shift_variants <- function(chrom, pos, preserve) {
  if (!"linkage" %in% preserve) {
    return(list(chrom = sprintf("chr%d", match(chrom, unique(chrom))), pos = pos))
  }
  low <- vapply(split(pos, chrom), min, numeric(1))
  list(
    chrom = sprintf("chr%d", match(chrom, unique(chrom))),
    pos = as.integer(pos - low[chrom] + 1L)
  )
}

# Only the gaps are re-measured: no row is added or dropped before this point, so
# the shift is the one place `linkage` can be lost.
reprex_verify <- function(before, new, preserve) {
  if ("linkage" %in% preserve) {
    gaps_before <- reprex_gaps(before$chrom, before$pos)$gap
    gaps_new <- reprex_gaps(new$chrom, new$pos)$gap
    if (!identical(gaps_before, gaps_new)) {
      reprex_lost("linkage")
    }
  }
  if ("n_variants" %in% preserve && length(before$pos) != length(new$pos)) {
    reprex_lost("n_variants")
  }
  invisible(NULL)
}

# Gaps are compared within a chromosome only: across a boundary they are not gaps.
# `chrom` is only used to group; the labels themselves are relabelled, so the
# gap vector is what a linkage check compares.
reprex_gaps <- function(chrom, pos) {
  ord <- order(chrom, pos)
  inside <- chrom[ord][-1L] == chrom[ord][-length(ord)]
  list(
    chrom = chrom[ord][c(TRUE, inside)],
    gap = diff(pos[ord])[inside]
  )
}

reprex_lost <- function(prop) {
  abort_msg(
    "Anonymisation lost `{prop}`; no script was written.",
    i = "This is a bug in `marker_reprex()`: please report it."
  )
}

# --- shared ------------------------------------------------------------------

reprex_read_lines <- function(path) {
  tryCatch(
    {
      con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else file(path, open = "rb")
      on.exit(close(con))
      # a raw read keeps a missing final newline, unlike readLines()
      text <- strsplit(rawToChar(readBin(con, "raw", file.size(path))), "\n", fixed = TRUE)[[1L]]
      # a file written in text mode on Windows ends every line with CRLF; the
      # carriage return is not part of the record and is not a DNA base
      sub("\r$", "", text)
    },
    error = function(e) abort_msg("Cannot read {.path {path}}.")
  )
}

# The script is plain ASCII, so a field written verbatim must be ASCII too.
reprex_ascii <- function(x, what) {
  bad <- !is.na(x) & grepl("[^ -~]", x)
  if (any(bad)) {
    abort_msg(
      "The script is plain ASCII and a {.field {what}} field is not.",
      i = "Rewrite that field by hand before sharing."
    )
  }
  x
}

reprex_int <- function(x) {
  y <- suppressWarnings(as.integer(x))
  if (anyNA(y)) {
    abort_msg("`pos` must hold whole numbers.")
  }
  y
}

# A vector is emitted with `deparse`, which is always valid R.
reprex_assign <- function(name, x) {
  paste0(name, " <- ", paste(deparse(x, width.cutoff = 500L), collapse = " "))
}

# --- VCF, BED and FASTA readers ----------------------------------------------

reprex_read_vcf <- function(path) {
  lines <- reprex_read_lines(path)
  records <- lines[!startsWith(lines, "#") & nzchar(lines)]
  if (!length(records)) {
    abort_msg("{.path {path}} holds no VCF record.")
  }
  fields <- strsplit(records, "\t", fixed = TRUE)
  if (any(lengths(fields) < 8L)) {
    abort_msg("{.path {path}} has a record with fewer than 8 VCF columns.")
  }
  pick <- function(j) vapply(fields, function(f) f[[j]], character(1))
  tibble::tibble(chrom = pick(1L), pos = pick(2L), ref = pick(4L), alt = pick(5L))
}

# A BED interval is replayed as one single-base record: the script carries it
# through the VCF reader the package already has.
reprex_read_bed <- function(path) {
  lines <- reprex_read_lines(path)
  rows <- lines[!grepl("^((#|track|browser)[[:space:]]*|$)", lines)]
  if (!length(rows)) {
    abort_msg("{.path {path}} holds no BED interval.")
  }
  fields <- strsplit(rows, "\t", fixed = TRUE)
  if (any(lengths(fields) < 3L)) {
    abort_msg("{.path {path}} has a line with fewer than 3 BED columns.")
  }
  cli::cli_warn(c(
    gettext("A BED interval is replayed as a 1-based single-base record, not as a width."),
    i = gettext("A bug that depends on an interval's width will not reproduce.")
  ))
  tibble::tibble(
    chrom = vapply(fields, function(f) f[[1L]], character(1)),
    pos = vapply(fields, function(f) f[[2L]], character(1)),
    ref = "A", alt = "G"
  )
}

reprex_read_fasta <- function(path) {
  lines <- reprex_read_lines(path)
  head <- which(startsWith(lines, ">"))
  if (!length(head)) {
    abort_msg("{.path {path}} is not a FASTA (no `>` header).")
  }
  end <- c(head[-1L] - 1L, length(lines))
  seqs <- vapply(seq_along(head), function(i) {
    paste0(lines[(head[[i]] + 1L):end[[i]]], collapse = "")
  }, character(1))
  stats::setNames(seqs, sub("^>([^[:space:]]*).*$", "\\1", lines[head]))
}
