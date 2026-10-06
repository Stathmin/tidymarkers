# =============================================================================
#  kasp_target.R  --  sequence retrieval, target SNP choice
# =============================================================================
# =============================================================================
#  1. Sequence retrieval
# =============================================================================

#' Perform an Ensembl REST request (JSON, retry, throttle, user agent)
#'
#' Ensembl asks for a descriptive `User-Agent` and enforces roughly 15 requests
#' per second, so requests are throttled and retried with backoff.
#' @noRd
ensembl_perform <- function(req) {
  req |>
    httr2::req_headers(`Content-Type` = "application/json") |>
    httr2::req_user_agent("tidykasp") |>
    httr2::req_throttle(rate = 15) |>
    httr2::req_retry(max_tries = 5L, backoff = ~3) |>
    httr2::req_perform()
}

#' Fetch a gene record from Ensembl Plants REST (cached on disk)
#' @param gene_id Ensembl gene id (e.g. `"AT4G00650"`).
#' @param species Ensembl species name.
#' @param cache cache directory.
#' @return a one-row tibble (`gene_id`, `symbol`, `chrom`, `start`, `end`,
#'   `strand`, `description`).
#' @examples
#' \dontrun{
#' ensembl_gene("AT4G00650")
#' }
#' @family KASP input
#' @export
ensembl_gene <- function(gene_id, species = "arabidopsis_thaliana",
                         cache = kasp_cache_dir()) {
  fs::dir_create(cache)
  path <- fs::path(cache, glue::glue("{species}_{gene_id}.json"))
  if (!fs::file_exists(path)) {
    resp <- httr2::request("https://rest.ensembl.org") |>
      httr2::req_url_path_append("lookup", "id", gene_id) |>
      httr2::req_url_query(expand = 1) |>
      ensembl_perform()
    readr::write_lines(httr2::resp_body_string(resp), path)
  }
  j <- jsonlite::fromJSON(readr::read_file(path), simplifyDataFrame = FALSE)
  need <- c("id", "seq_region_name", "start", "end", "strand")
  missing <- setdiff(need, names(j))
  if (length(missing)) {
    rlang::abort(c(
      glue::glue(gettext(
        "ensembl_gene(): record is missing field(s) {paste(missing, collapse = ', ')}."
      )),
      i = gettext("The Ensembl REST schema may have changed.")
    ))
  }
  tibble::tibble(
    gene_id = j$id, symbol = j$display_name %||% j$id,
    chrom = ath_chrom(j$seq_region_name),
    start = as.integer(j$start), end = as.integer(j$end),
    strand = as.integer(j$strand),
    description = j$description %||% NA_character_
  )
}

#' Natural variants overlapping a region, from Ensembl Plants (1001 Genomes)
#'
#' Returns `id`, `chrom`, `pos`, `ref_genome` (the allele vector) and
#' `consequence` -- the Ensembl/Sequence-Ontology term(s).  The `overlap/region`
#' endpoint does **not** provide a VEP `impact` column, so [pick_target_snp()]
#' ranks on the SO terms; a table you build yourself with an `impact` column is
#' also accepted.
#' @param chrom,start,end region to fetch (chromosome is normalised).
#' @param species,cache as in [ensembl_gene()].
#' @return a tibble of variants (`id`, `chrom`, `pos`, `ref_genome`,
#'   `consequence`), one row per record.
#' @examples
#' \dontrun{
#' ensembl_variants("Chr4", 268901, 271503)
#' }
#' @family KASP input
#' @export
ensembl_variants <- function(chrom, start, end,
                             species = "arabidopsis_thaliana",
                             cache = kasp_cache_dir()) {
  fs::dir_create(cache)
  region <- glue::glue("{sub('^Chr', '', chrom)}:{start}-{end}")
  path <- fs::path(cache, glue::glue("{species}_var_{region}.json"))
  if (!fs::file_exists(path)) {
    resp <- httr2::request("https://rest.ensembl.org") |>
      httr2::req_url_path_append("overlap", "region", species, region) |>
      httr2::req_url_query(feature = "variation") |>
      ensembl_perform()
    readr::write_lines(httr2::resp_body_string(resp), path)
  }
  df <- jsonlite::fromJSON(readr::read_file(path)) |> tibble::as_tibble()
  need <- c("id", "seq_region_name", "start", "alleles")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    rlang::abort(c(
      glue::glue(gettext(
        "ensembl_variants(): response is missing field(s) {paste(missing, collapse = ', ')}."
      )),
      i = gettext("The Ensembl REST schema may have changed.")
    ))
  }
  if (!"consequence_type" %in% names(df)) df$consequence_type <- NA_character_
  df |>
    dplyr::transmute(
      id = .data$id, chrom = ath_chrom(.data$seq_region_name),
      pos = as.integer(.data$start), ref_genome = .data$alleles,
      consequence = .data$consequence_type
    )
}

#' Extract a DNA window centred on a SNP as a DNAString
#'
#' The window is copied into a standalone `DNAString`: `subseq()` of a
#' `DNAStringSet` keeps a reference to the genome-wide shared pool, which would
#' otherwise bloat anything we serialise (e.g. the saved assay).
#' @param genome DNAStringSet.
#' @param chrom chromosome name.
#' @param pos variant position.
#' @param flank window half-width (bp).
#' @return a list with `seq` (the window as a standalone `DNAString`) and
#'   `snp_i` (the variant's coordinate within it).
#' @examples
#' g <- Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGTACGTACGTACGT"))
#' snp_template(g, "Chr1", 10L, flank = 5L)
#' @family KASP input
#' @export
snp_template <- function(genome, chrom, pos, flank = 150L) {
  seq <- get_region(genome, chrom, pos - flank, pos + flank)
  seq <- Biostrings::DNAString(as.character(seq))
  list(seq = seq, snp_i = as.integer(flank) + 1L)
}

#' GC fraction of the flanks on each side of the SNP (EasyKASP direction hint)
#'
#' EasyKASP chooses the read direction by flank GC: it prefers the side whose
#' +/-`flank` bases sit closest to 45%.  Returned here as a diagnostic.
#' @param genome DNAStringSet.
#' @param chrom chromosome name.
#' @param pos variant position.
#' @param flank bases to look at on each side.
#' @return a tibble with `side` (`upstream`/`downstream`), `gc` and
#'   `distance_to_45`.
#' @examples
#' g <- Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGTACGTACGTACGT"))
#' strand_gc_hint(g, "Chr1", 15L, flank = 5L)
#' @family KASP design
#' @export
strand_gc_hint <- function(genome, chrom, pos, flank = 25L) {
  up <- as.character(get_region(genome, chrom, max(1L, pos - flank), pos - 1L))
  dn <- as.character(get_region(genome, chrom, pos + 1L, pos + flank))
  tibble::tibble(side = c("upstream", "downstream"), gc = c(gc_frac(up), gc_frac(dn))) |>
    dplyr::mutate(distance_to_45 = abs(.data$gc - 0.45))
}

#' Map Ensembl consequence/impact to a design priority (lower = preferred)
#'
#' If an Ensembl VEP `impact` is available it wins (the standard HIGH/MODERATE/
#' LOW/MODIFIER severity ladder in `IMPACT_RANK`).  Otherwise the consequence
#' string is split on commas and matched against the canonical SO terms in
#' `SO_TERM_RANK` by *full token* (not substring), taking the most severe; terms
#' with no known rank (and `NA`s) fall back to `SO_RANK_OTHER`.
#' @keywords internal
#' @noRd
so_rank <- function(consequence, impact = NULL) {
  term_rank <- purrr::map_int(consequence, function(x) {
    toks <- stringr::str_trim(stringr::str_split(x, ",")[[1L]])
    r <- SO_TERM_RANK[toks]
    r <- r[!is.na(r)]
    if (length(r)) min(unname(r)) else SO_RANK_OTHER
  })
  if (is.null(impact)) {
    return(term_rank)
  }
  impact_rank <- unname(IMPACT_RANK[match(toupper(impact), names(IMPACT_RANK))])
  dplyr::coalesce(impact_rank, term_rank)
}

#' 3'-seed sequence of a would-be AS primer anchored at a variant
#'
#' `NA` if the seed would run off the sequence.
#' @keywords internal
#' @noRd
seed_seq <- function(genome, chrom, pos, ref_len, strand, k = SEED_K) {
  n <- genome_length(genome, chrom)
  anchor <- as_anchor(pos, ref_len, strand)
  lo <- if (identical(strand, "+")) anchor - k + 1L else anchor
  hi <- if (identical(strand, "+")) anchor else anchor + k - 1L
  if (lo < 1L || hi > n) {
    return(NA_character_)
  }
  seg <- get_region(genome, chrom, lo, hi)
  if (identical(strand, "+")) {
    as.character(seg)
  } else {
    as.character(Biostrings::reverseComplement(seg))
  }
}

#' 3'-seed specificity of a would-be AS primer anchored at a variant
#'
#' The seed is the `k` 3'-terminal bases of the AS primer: `[anchor-k+1, anchor]`
#' on the plus strand, or the reverse complement of `[anchor, anchor+k-1]` on the
#' minus strand (see [as_anchor()]).  Counted over the whole genome, both
#' strands, like [specificity_seed()].  `NA` if the seed runs off the sequence.
#' @keywords internal
#' @noRd
seed_hits_at <- function(genome, chrom, pos, ref_len, strand,
                         k = SEED_K, max_mismatch = SEED_MISMATCH, BPPARAM = NULL) {
  s <- seed_seq(genome, chrom, pos, ref_len, strand, k)
  if (is.na(s)) {
    return(NA_integer_)
  }
  specificity_of(c(seed = s), genome,
    k = k, max_mismatch = max_mismatch,
    BPPARAM = BPPARAM
  )$hits
}

#' Choose a designable target variant from a variant table
#'
#' Keeps variants we can actually design on: a biallelic substitution block whose
#' REF and ALT are ACGT-only strings of *equal* length (SNP or MNP).  Indels and
#' multi-allelic records are left for a dedicated engine and are not silently
#' forced into a substitution.
#'
#' @section Selection policy:
#' Among the designable substitutions the best one is picked by, in order:
#' 1. the cleanest AS 3' window (`n_3p`, best of the two orientations);
#' 2. consequence severity -- the Ensembl VEP `impact` column if the table has it,
#'    otherwise the canonical SO terms in `SO_TERM_RANK` (`missense_variant`, ...);
#' 3. a transition (A>G / C>T) over a transversion;
#' 4. 3'-seed specificity (`seed_hits`, best of the two strands);
#' 5. the leftmost position.
#' SNPs with complementary alleles (A/T, C/G) are skipped -- a 3'-terminal base
#' cannot tell them apart.  The genome-wide seed scan runs only on the best
#' `specificity_candidates` loci.
#'
#' @param variants variant table (Ensembl REST shape).
#' @param genome DNAStringSet.
#' @param gene one-row gene record.
#' @param variant_id force a specific variant by id.
#' @param clean_3p number of flanking bases of the AS primer that should be free
#'   of *other* variants.  Because the AS primer may be placed on either strand,
#'   the "clean" count is the best (smallest) of the two orientations: the chosen
#'   variant just has to be clean in *at least one* direction.
#' @param specificity_candidates how many top-ranked loci to check for 3'-seed
#'   specificity.  The seed scan is genome-wide, so only this many of the best
#'   candidates (by the other criteria) are checked; set `0` to skip it.
#' @param BPPARAM optional `BiocParallel` backend for the seed scan (split by
#'   chromosome, see [specificity_seed()]).
#' @return a one-row tibble for the chosen variant (`id`, `chrom`, `pos`,
#'   `consequence`, `ref`, `alt`, `ref_len`, `transition`, `n_3p`, `seed_hits`).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' v <- ensembl_variants("Chr4", 268901, 271503)
#' pick_target_snp(v, genome, ensembl_gene("AT4G00650"))
#' }
#' @family KASP design
#' @export
pick_target_snp <- function(variants, genome, gene, variant_id = NULL, clean_3p = 5L,
                            specificity_candidates = 10L, BPPARAM = NULL) {
  if (!is.null(variant_id)) variants <- dplyr::filter(variants, .data$id == variant_id)

  # REF length comes from the VCF (its first allele); the reference bases are
  # then read from the genome, which is the truth we design against.
  v <- variants |>
    dplyr::mutate(
      ref_len = purrr::map_int(.data$ref_genome, ~ nchar(.x[1L])),
      ref = purrr::pmap_chr(
        list(.data$chrom, .data$pos, .data$ref_len),
        function(ch, p, n) as.character(get_region(genome, ch, p, p + n - 1L))
      ),
      # exactly one non-reference allele of the SAME length, else NA (dropped
      # below): a silent "first of several" choice would pick an arbitrary ALT,
      # and unequal lengths are indels, handled separately.
      alt = purrr::map2_chr(.data$ref_genome, .data$ref, function(alleles, ref) {
        other <- setdiff(alleles, ref)
        if (length(other) == 1L && nchar(other) == nchar(ref)) other else NA_character_
      })
    ) |>
    dplyr::filter(
      !is.na(.data$alt),
      stringr::str_detect(.data$ref, "^[ACGT]+$"),
      stringr::str_detect(.data$alt, "^[ACGT]+$"),
      # a SNP whose alleles are complementary (A/T, C/G) cannot be told
      # apart by a 3'-terminal base (one allele's primer pairs the other
      # allele's template), so KASP discrimination collapses - skip them
      !(.data$ref_len == 1L & paste0(pmin(.data$ref, .data$alt), pmax(.data$ref, .data$alt)) %in% c("AT", "CG"))
    )
  if (nrow(v) == 0L) rlang::abort(gettext("No usable substitution variant (SNP/MNP) left after filtering."))

  # every variant of the gene (incl. indels): used to penalise targets whose AS
  # 3' neighbourhood is itself variable (breaks allele specificity across lines)
  allpos <- variants$pos[variants$chrom == gene$chrom]

  # Ensembl VEP IMPACT, when the table carries it, is the library-defined
  # severity ladder; otherwise so_rank() falls back to the SO terms.  (`v[[..]]`
  # is NULL when the column is absent.)
  impact <- v[["impact"]]
  # optional population frequency: prefer the more informative (higher-MAF) locus
  # as a tie-break; absent, `maf` is all `NA` and changes nothing
  pf <- v[["pop_freq"]]

  ranked <- v |>
    dplyr::mutate(
      transition = .data$ref_len == 1L &
        paste0(.data$ref, .data$alt) %in% c("AG", "GA", "CT", "TC"),
      rank = so_rank(.data$consequence, impact),
      maf = if (is.null(pf)) NA_real_ else pmin(pf, 1 - pf),
      # best-case (either strand) count of other variants in the AS 3' window,
      # vectorised over the whole variant table (IRanges).  The windows sit
      # strictly outside the variant block, so the target cannot count itself.
      n_3p = pmin(
        IRanges::countOverlaps(
          IRanges::IRanges(.data$pos - clean_3p + 1L, .data$pos - 1L),
          IRanges::IRanges(allpos, allpos)
        ),
        IRanges::countOverlaps(
          IRanges::IRanges(.data$pos + .data$ref_len, .data$pos + .data$ref_len + clean_3p - 2L),
          IRanges::IRanges(allpos, allpos)
        )
      )
    ) |>
    dplyr::arrange(.data$n_3p, .data$rank, dplyr::desc(.data$transition), dplyr::desc(.data$maf), .data$pos)

  if (specificity_candidates <= 0L) {
    return(ranked |> dplyr::slice(1L) |> dplyr::select(-"rank", -"ref_genome", -"maf"))
  }

  # among the best candidates, prefer the most specific AS 3'-seed (best of the
  # two strands).  The genome-wide scan is one [specificity_seed()] call over
  # both strands' seeds of the whole shortlist, so it can fan out by chromosome.
  short <- dplyr::slice_head(ranked, n = specificity_candidates) |>
    dplyr::mutate(candidate = dplyr::row_number())
  # one seed per candidate per strand; a stable key carries it through the scan
  seed_tbl <- tidyr::crossing(candidate = short$candidate, strand = c("+", "-")) |>
    dplyr::mutate(
      seed = purrr::pmap_chr(
        list(
          short$chrom[.data$candidate], short$pos[.data$candidate],
          short$ref_len[.data$candidate], .data$strand
        ),
        function(ch, p, rl, st) seed_seq(genome, ch, p, rl, st)
      ),
      key = paste0(.data$candidate, .data$strand)
    ) |>
    dplyr::filter(!is.na(.data$seed))
  if (nrow(seed_tbl) > 0L) {
    seed_tbl <- specificity_of(stats::setNames(seed_tbl$seed, seed_tbl$key), genome,
      BPPARAM = BPPARAM
    ) |>
      dplyr::select(key = "primer", "hits") |>
      dplyr::right_join(seed_tbl, by = "key")
  } else {
    seed_tbl$hits <- integer()
  }
  # seed specificity breaks ties among otherwise comparable loci (it does not
  # outrank the biological/chemical preferences above); best of the two strands
  short <- seed_tbl |>
    dplyr::summarise(seed_hits = min(.data$hits), .by = "candidate") |>
    dplyr::right_join(short, by = "candidate")

  short |>
    dplyr::arrange(
      .data$n_3p, .data$rank, dplyr::desc(.data$transition),
      dplyr::desc(.data$maf), .data$seed_hits, .data$pos
    ) |>
    dplyr::slice(1L) |>
    dplyr::select(-"rank", -"ref_genome", -"maf", -"candidate")
}
