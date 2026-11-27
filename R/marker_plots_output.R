# =============================================================================
#  marker_plots_output.R  --  result and oligo-order writers
# =============================================================================
#  Split from marker_plots.R (roadmap ZT T7).
# =============================================================================
# =============================================================================
#  Output helpers
# =============================================================================

# A gene symbol or chromosome becomes part of a file name: keep it inside
# `dir` even when the caller passes a path-like label.
safe_name_part <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]", "_", paste(as.character(x), collapse = "_"))
  if (!nzchar(x) || x %in% c(".", "..")) {
    abort_msg("The output name is empty after sanitising; give a real gene symbol.")
  }
  x
}

#' Write the assay results as a few consolidated files
#'
#' Works for any assembled assay, not only KASP: every `marker_assay` carries the
#' same `assay`/`thermo`/`specificity`/`complementarity` fields (see
#' [assemble_primer_assay()]).  `write_kasp_results()` is the deprecated former
#' name of this function.
#'
#' Instead of many tiny CSVs, the outputs are grouped by what they describe and
#' the per-consumer formats are column subsets of the same data:
#'
#' * `marker_design.csv` -- candidate assays (one row per common-primer placement);
#' * `marker_primer_qc.csv` -- one row per primer: Tm/GC/hairpin/3'-stability +
#'   3'-seed specificity + in-silico PCR products;
#' * `marker_dimer_dg.csv` -- duplex dG37 per primer pair (long);
#' * `{base}_tidy.csv` / `_idt.csv` / `_lgc.csv` -- the order table and its
#'   IDT / LGC column subsets (see [export_primers()]).
#'
#' @return character vector of the written paths (invisibly).
#' @param assay assembled assay.
#' @param gene gene record (used for the file base name).
#' @param dir output directory.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci)
#' gene <- list(symbol = "AT4G05320")
#' basename(write_marker_results(assay, gene, dir = tempdir()))
#' }
#' @family output
#' @export
write_marker_results <- function(assay, gene, dir = tempdir()) {
  fs::dir_create(dir)
  base <- safe_name_part(glue::glue("marker_{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}"))

  design <- assay$assay |>
    dplyr::select(
      "rank", "strand", "as_ref", "as_alt", "com", "discrimination",
      "as_gc", "amplicon", "as_len", "as_n_var", "as_n_3p", "n_var"
    )

  qc <- assay$thermo |>
    dplyr::transmute(
      .data$name, .data$primer, .data$len, .data$tm, .data$gc,
      .data$hairpin_dg, .data$terminal_dg
    ) |>
    dplyr::left_join(dplyr::rename(assay$specificity, name = "primer"), by = "name")
  if (!is.null(assay$amplicons)) {
    qc <- dplyr::left_join(qc, dplyr::rename(assay$amplicons, name = "primer"),
      by = "name"
    )
  }

  paths <- c(
    fs::path(dir, "marker_design.csv"),
    fs::path(dir, "marker_primer_qc.csv"),
    fs::path(dir, "marker_dimer_dg.csv")
  )
  readr::write_csv(design, paths[1L])
  readr::write_csv(qc, paths[2L])
  readr::write_csv(assay$complementarity, paths[3L])
  if (!is.null(assay$blast)) {
    p <- fs::path(dir, "marker_blast_hits.csv")
    readr::write_csv(assay$blast, p)
    paths <- c(paths, p)
  }
  export_primers(assay, gene, dir = dir) # {base}_tidy.csv / _idt.csv / _lgc.csv
  invisible(c(
    paths,
    fs::path(dir, glue::glue("{base}_tidy.csv")),
    fs::path(dir, glue::glue("{base}_idt.csv")),
    fs::path(dir, glue::glue("{base}_lgc.csv"))
  ))
}

#' @rdname write_marker_results
#' @param ... arguments passed on to [write_marker_results()].
#' @export
write_kasp_results <- function(...) {
  .Deprecated("write_marker_results")
  write_marker_results(...)
}

#' Tidy table of the assay oligos (with KASP tails)
#' @param assay assembled assay.
#' @param gene gene record.
#' @return a tibble (`assay`, `role`, `oligo_id`, `sequence`, `seq_no_tail`,
#'   `tail`, `tm`, `gc`, `length`) with the three assay oligos and their KASP
#'   tails.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci)
#' gene <- list(symbol = "AT4G05320")
#' marker_primer_table(assay, gene)
#' }
#' @family output
#' @export
marker_primer_table <- function(assay, gene) {
  nm <- glue::glue("{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}_{assay$snp$ref}to{assay$snp$alt}")
  roles <- tibble::tibble(
    name        = c("AS_ref", "AS_alt", "COM"),
    role        = c("AS-ref (FAM)", "AS-alt (HEX/VIC)", "common (reverse)"),
    oligo_id    = c("AS1_FAM", "AS2_HEX", "COM"),
    tail        = c(KASP_TAIL_FAM, KASP_TAIL_HEX, ""),
    # the AS-alt oligo is the physical (possibly degenerate) primer, so it is the
    # degenerate consensus that is ordered/carried here
    seq_no_tail = c(assay$assay$as_ref[1], as_alt_sequence(assay), assay$assay$com[1])
  )

  assay$thermo |>
    dplyr::select("name", "tm", "gc") |>
    dplyr::inner_join(roles, by = "name") |> # join, not positional match()
    dplyr::mutate(
      assay = nm,
      sequence = paste0(.data$tail, .data$seq_no_tail), # with 5' tail
      gc = round(.data$gc, 2),
      length = nchar(.data$sequence)
    ) |>
    dplyr::select(
      "assay", "role", "oligo_id", "sequence", "seq_no_tail",
      "tail", "tm", "gc", "length"
    )
}

#' Write IDT, LGC and tidy order sheets for the assay primers
#'
#' - `*_idt.csv`  : IDT bulk oligo input (Name, Sequence, Scale, Purification)
#' - `*_lgc.csv`  : LGC/KASP assay order (one row per assay, wide)
#' - `*_tidy.csv` : the long tidy table
#' @param assay assembled assay.
#' @param gene gene record.
#' @param dir output directory.
#' @return the tidy primer table (see [marker_primer_table()]), invisibly; called
#'   for its side effect of writing the IDT/LGC/tidy order sheets.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci)
#' gene <- list(symbol = "AT4G05320")
#' head(export_primers(assay, gene, dir = tempdir()))
#' }
#' @family output
#' @export
export_primers <- function(assay, gene, dir = tempdir()) {
  fs::dir_create(dir)
  tbl <- marker_primer_table(assay, gene)
  base <- safe_name_part(glue::glue("marker_{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}"))

  tbl |>
    dplyr::transmute(
      Name = glue::glue("{assay}_{oligo_id}"),
      Sequence = .data$sequence, Scale = "25nm", Purification = "STD"
    ) |>
    readr::write_csv(fs::path(dir, glue::glue("{base}_idt.csv")))

  tibble::tibble(
    Assay      = unique(tbl$assay),
    AS1_FAM    = tbl$sequence[tbl$oligo_id == "AS1_FAM"],
    AS2_HEX    = tbl$sequence[tbl$oligo_id == "AS2_HEX"],
    Common     = tbl$sequence[tbl$oligo_id == "COM"]
  ) |>
    readr::write_csv(fs::path(dir, glue::glue("{base}_lgc.csv")))

  readr::write_csv(tbl, fs::path(dir, glue::glue("{base}_tidy.csv")))
  invisible(tbl)
}
