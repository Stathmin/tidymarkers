# =============================================================================
#  marker_target.R  --  sequence retrieval, target SNP choice
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
    httr2::req_user_agent("tidymarkers") |>
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
                         cache = marker_cache_dir()) {
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
    cli::cli_abort(c(
      gettext(
        "ensembl_gene(): record is missing field(s) {paste(missing, collapse = ', ')}."
      ),
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
                             cache = marker_cache_dir()) {
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
    cli::cli_abort(c(
      gettext(
        "ensembl_variants(): response is missing field(s) {paste(missing, collapse = ', ')}."
      ),
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
#' Every substitution is separable by its 3'-terminal base, including the
#' transversions A/T and C/G (their wrong-template mismatch is a strong
#' *same-base* pair).  The weak case is the G:T / A:C *wobble* (transitions), a
#' soft penalty carried by [wobble_mismatch()], not a reason to skip a site.
#' The genome-wide seed scan runs only on the best
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
#' @param scheme optional marker scheme to require (see [marker_schemes()]).  When
#'   given, ranked candidates are designed one by one with that scheme's engine
#'   (`marker_design_at()`) and the first *designable* one is returned - "doomed"
#'   depends on the marker type, so a KASP-doomed locus may still make a CAPS or
#'   ARMS marker.  `NULL` keeps the plain best-by-criteria pick.
#' @param max_tries how many ranked candidates to try when `scheme` is given.
#' @return a one-row tibble for the chosen variant (`id`, `chrom`, `pos`,
#'   `consequence`, `ref`, `alt`, `ref_len`, `transition`, `n_3p`, `seed_hits`);
#'   with `scheme`, the designed marker is attached as the `"marker"` attribute.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' v <- ensembl_variants("Chr4", 268901, 271503)
#' pick_target_snp(v, genome, ensembl_gene("AT4G00650"))
#' # a target that designs specifically as a CAPS marker:
#' pick_target_snp(v, genome, ensembl_gene("AT4G00650"), scheme = "caps")
#' }
#' @family KASP design
#' @export
pick_target_snp <- function(variants, genome, gene, variant_id = NULL, clean_3p = 5L,
                            specificity_candidates = 10L, BPPARAM = NULL,
                            scheme = NULL, max_tries = 10L) {
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
      stringr::str_detect(.data$alt, "^[ACGT]+$")
    )
  if (nrow(v) == 0L) cli::cli_abort(gettext("No usable substitution variant (SNP/MNP) left after filtering."))

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

  # among the best candidates, prefer the most specific AS 3'-seed (best of the
  # two strands).  The genome-wide scan is one [specificity_seed()] call over
  # both strands' seeds of the whole shortlist, so it can fan out by chromosome.
  short <- dplyr::slice_head(ranked, n = max(1L, specificity_candidates))
  if (specificity_candidates > 0L) {
    short <- dplyr::mutate(short, candidate = dplyr::row_number())
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
      seed_tbl <- specificity_of(rlang::set_names(seed_tbl$seed, seed_tbl$key), genome,
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
  } else {
    short$candidate <- seq_len(nrow(short))
    short$seed_hits <- NA_integer_
  }
  short <- dplyr::arrange(
    short,
    .data$n_3p, .data$rank, dplyr::desc(.data$transition),
    dplyr::desc(.data$maf), .data$seed_hits, .data$pos
  )
  best <- short |>
    dplyr::slice(1L) |>
    dplyr::select(-"rank", -"ref_genome", -"maf", -"candidate")
  if (is.null(scheme)) {
    return(best)
  }
  # "doomed" depends on the marker type: walk the ranked candidates until the
  # scheme can actually design one (an assay, a CAPS or an ARMS marker each have
  # their own way of failing)
  cand <- dplyr::slice_head(short, n = max(1L, as.integer(max_tries)))
  for (i in seq_len(nrow(cand))) {
    snp <- tibble::tibble(
      chrom = cand$chrom[i], pos = cand$pos[i], ref = cand$ref[i], alt = cand$alt[i]
    )
    marker <- tryCatch(
      suppressWarnings(marker_design_at(genome, snp, scheme, variants = variants, BPPARAM = BPPARAM)),
      error = function(e) NULL
    )
    if (is_target_designable(marker)) {
      out <- dplyr::select(cand[i, ], -"rank", -"ref_genome", -"maf", -"candidate")
      attr(out, "marker") <- marker
      return(out)
    }
  }
  cli::cli_warn(c(
    gettext("No ranked target is designable for scheme {scheme} ({nrow(cand)} tried)."),
    i = gettext("returning the best-ranked target, whose design may be doomed.")
  ))
  best
}

# =============================================================================
#  MarkerContext seam (ROADMAP ZA2)
# =============================================================================
#  `marker_run_context()` is the single object the design seam reads: the genome
#  the target is anchored to, plus the optional variant table,
#  population/haplotype table, reaction conditions and engines.  The public
#  `marker_context()` name remains the plan-layer grammar constructor in
#  marker_plan.R; this is deliberately a different interface.
#
#  Invariants:
#  * `variants`/`haplotypes` are stored exactly as supplied (no normalisation
#    until `marker_attach_haplotypes()` aligns them to a target);
#  * a missing population stays `NULL` - its frequency is `NA`, never 1.
# =============================================================================

#' A marker run context
#'
#' Bundles what a design run needs besides the target itself.  `population` is
#' the tidy haplotype table (`haplotype`, `freq`, `chrom`, `pos`, `ref`, `alt`)
#' used by [marker_attach_haplotypes()]; `conditions` and `engines` are carried
#' for the reaction/engine seams.
#' @param genome a `marker_genome` or `DNAStringSet`, or `NULL`.
#' @param variants optional variant table.
#' @param population optional population/haplotype table.
#' @param haplotypes alias for a haplotype table when `population` is absent.
#' @param conditions optional [reaction_conditions()].
#' @param engines optional [marker_engines()].
#' @param options free-form run options.
#' @return a classed `marker_run_context` list.
#' @family KASP plan
#' @export
marker_run_context <- function(genome = NULL, variants = NULL, population = NULL,
                               haplotypes = NULL, conditions = NULL, engines = NULL,
                               options = list()) {
  # a haplotype table is the population in the current vocabulary; keep the
  # supplied table unchanged and let marker_attach_haplotypes() validate it.
  if (is.null(haplotypes) && is_haplotype_table(population)) {
    haplotypes <- population
  }
  structure(
    list(
      genome = genome,
      variants = variants,
      population = population,
      haplotypes = haplotypes,
      conditions = conditions,
      engines = engines,
      options = options
    ),
    class = "marker_run_context"
  )
}

is_haplotype_table <- function(x) {
  is.data.frame(x) && all(c("haplotype", "freq") %in% names(x))
}

#' Frequency carried by a marker context (NA when the population is unknown)
#' @keywords internal
#' @noRd
marker_run_context_frequency <- function(context) {
  if (is.null(context) || !inherits(context, "marker_run_context")) {
    return(NA_real_)
  }
  pop <- context$population
  if (is.null(pop)) {
    return(NA_real_)
  }
  if (is.data.frame(pop) && "alt_freq" %in% names(pop) && nrow(pop)) {
    return(as.numeric(pop$alt_freq[[1L]]))
  }
  NA_real_
}

#' Resolve a context's variant/haplotype inputs for one design call
#' @keywords internal
#' @noRd
marker_run_context_args <- function(context, variants = NULL, haplotypes = NULL) {
  if (is.null(context)) {
    return(list(variants = variants, haplotypes = haplotypes))
  }
  if (inherits(context, "marker_layer")) {
    cli::cli_abort(gettext(
      "`context` is a marker_plan() layer; pass marker_run_context(genome, ...) instead."
    ))
  }
  if (!inherits(context, "marker_run_context")) {
    cli::cli_abort(gettext("`context` must be a marker_run_context() object."))
  }
  list(
    variants = variants %||% context$variants,
    haplotypes = haplotypes %||% context$haplotypes
  )
}

#' Design a target with one scheme's engine (for designability checks)
#'
#' "Doomed" is not one thing: a SNP may be undesignable as a KASP but fine as a
#' CAPS or ARMS marker.  This dispatches to the scheme's designer so
#' [pick_target_snp()] can walk candidates until one actually designs.
#' @param context optional [marker_run_context()] supplying `variants`/`haplotypes`;
#'   an explicit argument wins over the context.
#' @keywords internal
#' @noRd
marker_design_at <- function(genome, snp, scheme, variants = NULL, haplotypes = NULL,
                             BPPARAM = NULL, ..., context = NULL) {
  ctx <- marker_run_context_args(context, variants, haplotypes)
  variants <- ctx$variants
  haplotypes <- ctx$haplotypes
  # `primer` is deliberately absent: it is not a genotyping chemistry (no allele
  # is read) but the marker-agnostic "any primer pair" QC subset, so it has no
  # locus-level designer and is documentation-only (C13).  Its criteria remain
  # evaluable - marker_design_checks(assay, scheme = "primer") and
  # marker_decision_graph(scheme = "primer") work on any assembled pair - they just
  # are not a marker the rescue ladder can return.  The one check it still lacks
  # is a pair-level "in-silico PCR: exactly one product" node; registering it
  # belongs in marker_criteria.R/marker_schemes.R (a requested serial slot, not this
  # file).
  out <- switch(scheme,
    snp = ,
    indel = ,
    multi = assemble_kasp_assay(
      genome, snp,
      variants = variants, scheme = scheme,
      haplotypes = haplotypes, BPPARAM = BPPARAM, ...
    ),
    caps = design_caps_marker(genome, snp, haplotypes = haplotypes, ...),
    arms = design_arms_marker(genome, snp, ...),
    hrm = design_hrm_marker(genome, snp, ...),
    {
      # a discovered scheme (R/marker_<scheme>.R) supplies its own designer; a
      # composed technology dispatches to its base designer
      # (extra_scheme_designer).  Pass the shared context only when the designer
      # accepts it, so an existing scheme's designer never sees an extra argument.
      f <- extra_scheme_designer(scheme)
      if (is.null(f)) {
        NULL
      } else {
        fm <- names(formals(f))
        ctx <- list()
        if ("variants" %in% fm) ctx$variants <- variants
        if ("haplotypes" %in% fm) ctx$haplotypes <- haplotypes
        if ("BPPARAM" %in% fm) ctx$BPPARAM <- BPPARAM
        rlang::exec(f, genome, snp, !!!ctx, ...)
      }
    }
  )
  # attach the population layer to every chemistry: haplotype-aware engines use
  # it during design; the rest use it for plotting and the check column
  if (!is.null(out) && !is.null(haplotypes)) {
    out <- marker_attach_haplotypes(out, haplotypes, snp)
  }
  out
}

#' Did a scheme's designer return a usable marker? (see `marker_design_at()`)
#' @keywords internal
#' @noRd
is_target_designable <- function(marker) {
  if (inherits(marker, "marker_doomed")) {
    return(FALSE)
  }
  # marker_assay / marker_caps / marker_arms / marker_hrm and any discovered
  # scheme's marker_<scheme> class are designable
  any(startsWith(class(marker), "marker_"))
}

#' Design a marker with the first chemistry that works
#'
#' A locus KASP dooms is not a dead end: it may still make a CAPS/dCAPS, ARMS or
#' HRM marker (see [marker_schemes()]).  This walks the schemes in order, designing
#' each in turn, and returns the first usable one together with a tidy verdict for
#' every scheme it tried - "which chemistry takes this target", in one call.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`; an `alt2` list-column
#'   is honoured by the KASP scheme).
#' @param strategy a [rescue_strategy()] object; its `tried` rows give the order.
#'   `NULL` (with `schemes = NULL`) uses the default per-locus ladder.
#' @param schemes scheme ids to try, in order (see [marker_schemes()]).
#' @param ... passed to the per-scheme designer (e.g. `flank`).
#' @param context optional [marker_run_context()] supplying `variants`/`haplotypes`
#'   (and, for later seams, `population`, `conditions`, `engines`).  Must be
#'   named; supplied through `...` before this argument it would bind there.
#' @return a list: `snp`, `scheme` (the winner, or `NA`), `marker` (the winning
#'   object, or `NULL`), and `tried` (a tibble `scheme`, `ok`, `marker`).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_marker(genome, marker_example("SNP")[1, ])
#' }
#' @family KASP design
#' @export
design_marker <- function(genome, snp, strategy = NULL, schemes = NULL, ...,
                          context = NULL) {
  # the order to try is the rescue policy (see rescue_strategy()); a plain vector
  # keeps the historical KASP -> CAPS -> ARMS -> HRM ladder
  if (is.null(schemes)) {
    schemes <- if (is.null(strategy)) c("snp", "caps", "arms", "hrm") else strategy$scheme[strategy$tried]
  }
  tried <- purrr::map(schemes, function(sc) {
    m <- tryCatch(
      suppressWarnings(marker_design_at(genome, snp, sc, ..., context = context)),
      error = function(e) NULL
    )
    tibble::tibble(scheme = sc, ok = is_target_designable(m), marker = list(m))
  }) |>
    purrr::list_rbind()
  win <- tried[tried$ok, , drop = FALSE]
  scheme <- if (nrow(win)) win$scheme[1L] else NA_character_
  list(
    snp = snp, scheme = scheme,
    marker = if (is.na(scheme)) NULL else win$marker[[1L]],
    tried = tried
  )
}
