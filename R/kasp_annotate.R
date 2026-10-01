# =============================================================================
#  kasp_annotate.R  --  structured consequences from a Bioconductor annotation
# =============================================================================
#  `pick_target_snp()` ranks variants by consequence severity, and until now the
#  `consequence` string came from the VCF/Ensembl record and was matched against
#  a hand-kept dictionary (`SO_TERM_RANK`).  When a `TxDb` is available the
#  transcript model itself can produce that string: `VariantAnnotation` locates
#  the variant in the gene model (`locateVariants()`) and, given a reference
#  source, calls the codon consequence (`predictCoding()`).  Nothing here is
#  reimplemented - this file only wires the Bioconductor annotation into the
#  tidy variant table so the rest of the pipeline keeps working unchanged.
#
#  Everything is optional: with no `TxDb` installed `annotate_variants()` simply
#  explains what to install, and the pipeline still runs off the input
#  consequence.
# =============================================================================

#' Default Arabidopsis transcript annotation (`TxDb`)
#'
#' Loads `TxDb.Athaliana.BioMart.plantsmart28` (a `Suggests` dependency) without
#' attaching it to the search path.
#'
#' @param package TxDb package name.
#' @return a `TxDb` object.
#' @examples
#' \dontrun{
#' txdb <- kasp_txdb()
#' }
#' @family KASP annotation
#' @export
kasp_txdb <- function(package = "TxDb.Athaliana.BioMart.plantsmart28") {
  load_namespace_quietly(package, reason = gettext("for structured variant annotation"))
  getExportedValue(package, package)
}

# The Bioconductor stack announces a benign export collision
# ("replacing previous import ... when loading 'SummarizedExperiment'") whenever
# `VariantAnnotation`/a `TxDb` first pulls `SummarizedExperiment` in.  It is not
# actionable by users, so load the namespace once and muffle only that message.
# `rlang::is_installed()` would load the namespace (and thus emit the warning)
# itself, so gate on `system.file()`, which never loads.
#' @keywords internal
#' @noRd
load_namespace_quietly <- function(pkg, reason = NULL) {
  if (!nzchar(system.file(package = pkg))) rlang::check_installed(pkg, reason = reason)
  withCallingHandlers(
    requireNamespace(pkg, quietly = TRUE),
    warning = function(w) {
      if (stringr::str_detect(conditionMessage(w), "replacing previous import")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  invisible(TRUE)
}

#' Default Arabidopsis reference genome as a `BSgenome`
#'
#' Loads `BSgenome.Athaliana.TAIR.TAIR9` (a `Suggests` dependency), the reference
#' source `predictCoding()` needs for codon-level consequences.
#'
#' @param package BSgenome package name.
#' @return a `BSgenome` object.
#' @examples
#' \dontrun{
#' bsgenome <- kasp_bsgenome()
#' }
#' @family KASP annotation
#' @export
kasp_bsgenome <- function(package = "BSgenome.Athaliana.TAIR.TAIR9") {
  load_namespace_quietly(package, reason = gettext("for codon-level annotation"))
  getExportedValue(package, package)
}

#' Build (or load) an `EnsDb` gene model from the Ensembl Plants GTF
#'
#' `ensembldb` is the Bioconductor way to hold a gene model as an `EnsDb`, built
#' from a GTF.  There is **no** prebuilt Arabidopsis `EnsDb` package, so this
#' builds one (cached) from the current Ensembl Plants GTF and gives an
#' up-to-date model from a single source.
#'
#' Note on consequences: `ensembldb` exposes no `locateVariants()`/`predictCoding()`
#' and an `EnsDb` cannot be coerced to a `TxDb`, so the codon-level annotation in
#' [annotate_variants()] keeps using a `TxDb` (see [kasp_txdb()]).  This builder is
#' for gene/transcript/exon queries from the current model; converting it to a
#' `TxDb` via `txdbmaker::makeTxDbFromGRanges()` is possible but not lossless
#' (CDS rows are dropped), which is exactly why the `TxDb` path is the annotation
#' backend.
#'
#' @param gtf GTF path; when `NULL`, the Ensembl Plants GTF is downloaded (once).
#' @param cache cache directory for the GTF and the built `EnsDb`.
#' @param release Ensembl Plants release number.
#' @return an `EnsDb` object.
#' @examples
#' \dontrun{
#' edb <- kasp_ensdb()
#' }
#' @family KASP annotation
#' @export
kasp_ensdb <- function(gtf = NULL, cache = kasp_cache_dir(), release = 63L) {
  rlang::check_installed("ensembldb", reason = gettext("for an EnsDb gene model"))
  fs::dir_create(cache)
  if (is.null(gtf)) {
    gtf <- fs::path(cache, glue::glue("Arabidopsis_thaliana.TAIR10.{release}.gtf.gz"))
    if (!fs::file_exists(gtf)) {
      url <- glue::glue(
        "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-{release}/gtf/",
        "arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.{release}.gtf.gz"
      )
      cli::cli_alert_info(gettext("Downloading the Ensembl Plants GTF into {.path {cache}}."))
      utils::download.file(url, gtf, mode = "wb", quiet = TRUE)
    }
  }
  out <- fs::path(cache, "Arabidopsis_thaliana.EnsDb.sqlite")
  if (!fs::file_exists(out)) {
    cli::cli_alert_info(gettext("Building the EnsDb (one-off) from {.file {fs::path_file(gtf)}}."))
    ensembldb::ensDbFromGtf(
      gtf = as.character(gtf), outfile = as.character(out),
      organism = "Arabidopsis_thaliana",
      genomeVersion = "TAIR10", version = release
    )
  }
  ensembldb::EnsDb(out)
}

# VariantAnnotation spells `predictCoding()`'s CONSEQUENCE in its own small
# vocabulary; map it onto the canonical Sequence-Ontology terms `SO_TERM_RANK`
# already ranks, so both annotation paths land on the same ladder.
VA_TO_SO <- c(
  "nonsense"       = "stop_gained",
  "nonsynonymous"  = "missense_variant",
  "synonymous"     = "synonymous_variant",
  "frameshift"     = "frameshift_variant",
  "stop codon"     = "stop_lost",
  "not translated" = "5_prime_UTR_variant",
  "spliceSite"     = "splice_region_variant"
)

# `locateVariants()` LOCATION -> SO term, used when no reference source is given
# (i.e. codon-level calls are impossible).
LOCATION_TO_SO <- c(
  coding       = "coding_sequence_variant",
  fiveUTR      = "5_prime_UTR_variant",
  threeUTR     = "3_prime_UTR_variant",
  intron       = "intron_variant",
  spliceSite   = "splice_region_variant",
  promoter     = "upstream_gene_variant",
  intergenic   = "intergenic_variant"
)

#' Align genome sequence names to a `TxDb`'s seqlevels
#'
#' A `TxDb` often names its seqlevels `1`..`5`/`Mt`/`Pt` while a genome uses
#' `Chr1`..`Chr5`/`ChrM`/`ChrC`.  This maps the genome names onto the `TxDb`
#' names using the package's `ATH_CHROM_MAP`, falling back to the names unchanged
#' (for genomes that already share the `TxDb` convention).
#' @keywords internal
#' @noRd
txdb_seqname_map <- function(chrom, txdb) {
  inv <- stats::setNames(names(ATH_CHROM_MAP), ATH_CHROM_MAP) # values back to names
  txlev <- GenomeInfoDb::seqlevels(txdb)
  mapped <- dplyr::coalesce(unname(inv[chrom]), chrom)
  ifelse(mapped %in% txlev, mapped, chrom)
}

#' Annotate variants with a `TxDb` transcript model
#'
#' Fills (or replaces) the `consequence` column with the annotation the
#' transcript model itself gives, instead of the string that came in with the
#' record.  `VariantAnnotation::locateVariants()` always runs (region of the gene
#' model); when a `seqsource` reference is supplied,
#' `VariantAnnotation::predictCoding()` also runs and provides the codon-level
#' term (`stop_gained`, `missense_variant`, `synonymous_variant`, ...), which is
#' the most discriminating thing `pick_target_snp()` can rank on.
#'
#' The input and output are the tidy variant table used across the package
#' (`chrom`, `pos`, `ref_genome` list-column of `c(ref, alt)`), so annotation is
#' just another step in the pipe:
#'
#' ```r
#' read_variants_vcf(vcf) |>
#'   annotate_variants(seqsource = kasp_bsgenome()) |>
#'   pick_target_snp(genome, gene)
#' ```
#'
#' @param variants tidy variant tibble (see [read_variants_vcf()]).
#' @param txdb a `TxDb`; defaults to [kasp_txdb()].
#' @param seqsource optional reference (`BSgenome` or a `DNAStringSet` whose
#'   seqlevels match the `TxDb`) enabling codon-level consequences; without it
#'   only the region of the gene model is derived.
#' @param gene optional gene filter (`$chrom`, `$start`, `$end`) to restrict the
#'   annotation to one locus.
#' @return `variants` with `consequence`, plus `gene_id` from the `TxDb` (and
#'   `codon_change` when `seqsource` is supplied).
#' @examples
#' \dontrun{
#' annotate_variants(
#'   read_variants_vcf(vcf),
#'   txdb = kasp_txdb(), seqsource = kasp_bsgenome()
#' )
#' }
#' @family KASP annotation
#' @export
annotate_variants <- function(variants, txdb = NULL, seqsource = NULL, gene = NULL) {
  load_namespace_quietly("VariantAnnotation", reason = "for structured annotation")
  txdb <- txdb %||% kasp_txdb()
  if (nrow(variants) == 0L) {
    return(variants)
  }

  if (!is.null(gene) && !is.null(gene$chrom)) {
    variants <- dplyr::filter(variants, .data$chrom == ath_chrom(gene$chrom))
  }

  refs <- purrr::map_chr(variants$ref_genome, 1L)
  alts <- purrr::map_chr(variants$ref_genome, 2L)
  gr <- GenomicRanges::GRanges(
    txdb_seqname_map(variants$chrom, txdb),
    IRanges::IRanges(variants$pos, width = nchar(refs)),
    ref = refs, alt = alts
  )

  # VariantAnnotation/AnnotationDbi emit benign chatter here (out-of-bound
  # notes from a renamed seqinfo, `select()` progress); quieten it so the
  # returned table is the only output
  loc <- suppressWarnings(suppressMessages(
    VariantAnnotation::locateVariants(gr, txdb, VariantAnnotation::AllVariants())
  ))
  loc <- tibble::as_tibble(as.data.frame(loc)) |>
    dplyr::mutate(
      query = as.integer(.data$QUERYID),
      LOCATION = as.character(.data$LOCATION),
      GENEID = as.character(.data$GENEID)
    ) |>
    # one variant can sit in several transcripts; keep the most "coding" location
    dplyr::mutate(loc_rank = match(.data$LOCATION, names(LOCATION_TO_SO))) |>
    dplyr::arrange(.data$query, .data$loc_rank) |>
    dplyr::summarise(
      so_location = unname(LOCATION_TO_SO[utils::head(.data$LOCATION, 1L)]),
      gene_id = utils::head(.data$GENEID, 1L),
      .by = "query"
    )

  if (is.null(variants$consequence)) variants$consequence <- NA_character_
  out <- dplyr::mutate(variants, query = dplyr::row_number()) |>
    dplyr::left_join(loc, by = "query")

  # codon-level calls, when we have a reference to read codons from; they win
  # over the region-level term because they are the discriminating ones
  if (!is.null(seqsource)) {
    if (inherits(seqsource, "BSgenome")) {
      seqsource <- GenomeInfoDb::renameSeqlevels(seqsource, txdb_seqname_map(GenomeInfoDb::seqlevels(seqsource), txdb))
    }
    pc <- suppressWarnings(suppressMessages(
      VariantAnnotation::predictCoding(gr, txdb, seqsource,
        varAllele = Biostrings::DNAStringSet(alts)
      )
    ))
    coding <- tibble::as_tibble(as.data.frame(pc)) |>
      dplyr::mutate(
        query = as.integer(.data$QUERYID),
        CONSEQUENCE = as.character(.data$CONSEQUENCE)
      ) |>
      dplyr::arrange(.data$query) |>
      dplyr::summarise(
        so_codon = unname(VA_TO_SO[utils::head(.data$CONSEQUENCE, 1L)]),
        codon_change = paste(unique(paste0(.data$REFCODON, ">", .data$VARCODON)), collapse = ","),
        .by = "query"
      )
    out <- dplyr::left_join(out, coding, by = "query")
  } else {
    out$so_codon <- NA_character_
  }

  out |>
    dplyr::mutate(consequence = dplyr::coalesce(.data$so_codon, .data$so_location, .data$consequence)) |>
    dplyr::select(-"query", -"so_location", -"so_codon")
}
