# =============================================================================
#  kasp_assemble.R  --  assemble and score one complete assay
# =============================================================================
# =============================================================================
#  7. Assemble + score one complete assay
# =============================================================================

#' A target no assay clears the gates for (a partial design record)
#'
#' Returned by [assemble_kasp_assay()] instead of aborting when the strict gates
#' leave no AS primer / no pairing.  It carries the target, the template and the
#' reason (and the candidate table when one was built), so a caller can report or
#' plot the doomed target rather than lose it to an error.
#' @keywords internal
#' @noRd
kasp_doomed <- function(reason, stage, snp, template, as_candidates = NULL) {
  structure(
    list(snp = snp, template = template, reason = reason, stage = stage, as_candidates = as_candidates),
    class = "kasp_doomed"
  )
}

#' @keywords internal
#' @noRd
is_kasp_doomed <- function(x) inherits(x, "kasp_doomed")

#' Drop the grader's helper columns from an AS candidate table
#' @keywords internal
#' @noRd
tidy_as_candidates <- function(x) {
  drop <- grep("^\\.neg_|^m_|^gate_|^tier|^as_row$", names(x), value = TRUE)
  x[, setdiff(names(x), drop), drop = FALSE]
}

#' Assemble and score the best KASP assay for a target SNP
#'
#' When `variants` is supplied the selection is *variant-aware*: among primers
#' that pass the physico-chemical filters we minimise the number of other
#' variants falling under the primer footprint (the 3'-terminal 5 nt are handled
#' earlier by [pick_target_snp()]), then maximise discrimination.
#'
#' @section Selection policy:
#' The AS placement is chosen lexicographically by:
#' 1. cleanest 3' window (`n_3p`, best of the two orientations);
#' 2. fewest other variants under the primer footprint (`n_var`);
#' 3. highest discrimination (DECIPHER);
#' 4. not a weak G:T / A:C wobble mismatch;
#' 5. quieter 3'-seed genome-wide (`seed_hits`, best of the two strands) -- a
#'    specificity tie-break computed only for the top `specificity_candidates`
#'    placements, mirroring [pick_target_snp()];
#' 6. then `|Tm - target|`, 3'-stability and length.
#' The common primer is then placed on the winning strand and ranked by
#' amplicon length, `COM Tm >= AS Tm`, dimers/hairpin, footprint load and Tm
#' proximity (see [design_common_primers()]).
#' @section KASP design pipeline:
#' Typical use is `ath_genome()` -> `ensembl_gene()`/`ensembl_variants()` ->
#' `pick_target_snp()` -> `assemble_kasp_assay()` -> `kasp_design_checks()`.
#' See `vignette("kasp-intro")`.
#'
#' @param genome DNAStringSet (from [ath_genome()]).
#' @param snp    one-row tibble (chrom, pos, ref, alt) from [pick_target_snp()].
#' @param variants variant table (for the same chromosome) to count overlaps.
#' @param strands which AS-primer orientations to search (see [design_as_primers()]).
#' @param tm_range accept AS primers whose Tm falls in this window (C).
#' @param clean_3p size of the "clean" AS 3' window used to count other variants.
#' @param max_footprint_variants hard cap on *other* variants under the AS primer
#'   footprint ("N-masking" as a constraint): `0` allows none, the default `Inf`
#'   only penalises them in the ranking.  Use [mask_positions()] to blank such
#'   positions in your own template if you need the literal `N`s.
#' @param min_disc minimum discrimination (on- minus off-target efficiency).
#' @param specificity_candidates how many leading AS placements to screen for
#'   3'-seed genome-wide uniqueness when breaking ties (mirrors
#'   [pick_target_snp()]; `0` disables the screen).
#' @param pair_screen how many leading pairings (by the genome-free ranking) to
#'   genome-screen for off-target exact primer hits, a second in-silico PCR
#'   product and an inverted repeat; the first clean one wins.  These need a
#'   genome/PCR pass, so only the leaders are screened (`0` disables; the ranked
#'   winner is then kept).
#' @param run_blast also build/query a genome BLAST db (needs NCBI BLAST+).
#' @param cores processors for the per-SNP subprocess work (the DECIPHER PCR
#'   efficiencies and the OligoArrayAux dG maps); `1` (the default, from
#'   `getOption("tidykasp.cores")`) keeps CRAN happy, and larger values give a
#'   near-linear speed-up of a single assay.  When `BPPARAM` is set it wins.
#' @param BPPARAM optional `BiocParallel` backend for the genome-wide steps
#'   (seed specificity and in-silico PCR; split by chromosome) and for the
#'   subprocess work.  Use it to avoid nesting: inside a locus-level backend pass
#'   `cores = 1` / a serial `BPPARAM`.
#' @return a `kasp_assay` object: a list holding the `snp`, the template, the
#'   chosen `strand`, ranked `assay`/`as_candidates` tibbles, `efficiency`,
#'   `thermo`, `complementarity`, `specificity`, `amplicons`, `blast` results and
#'   the tailed `tails` primers (plus `haplotypes`/`haplotype_hits` when a
#'   haplotype table was given).  When no candidate clears the gates it returns a
#'   `kasp_doomed` object instead (the target, template, `reason`, `stage` and the
#'   candidate table), so a batch or a plot can report the target rather than
#'   abort.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' v <- ensembl_variants("Chr4", 268901, 271503)
#' snp <- pick_target_snp(v, genome, ensembl_gene("AT4G00650"))
#' assemble_kasp_assay(genome, snp)
#' }
#' @section External tool:
#' The thermodynamics engine (DECIPHER via OligoArrayAux) needs `hybrid-min` on
#' the `PATH`; if the data tables are not found automatically, point them out:
#' `with_unafloddat("/path/to/share/oligoarrayaux", assemble_kasp_assay(...))`.
#' @param flank template half-width (bp).
#' @param linked optional two-column table (`target`, `linked` positions) of
#'   co-segregating substitutions; those riding with `snp$pos` are baked into the
#'   ALT primer so it mismatches REF twice (see [design_as_primers()]).
#' @param haplotypes optional tidy haplotype table
#'   (`haplotype`, `freq`, `chrom`, `pos`, `ref`, `alt`).  Its per-position
#'   frequencies become the `pop_freq` weights and its co-segregating positions
#'   the `linked` set (both derived, see `haplotype_variants()`); the object also
#'   carries `haplotypes` and a per-haplotype `haplotype_hits` table (which AS
#'   primer each is expected to bind).
#' @param min_gc,max_gc accepted AS GC window.
#' @family KASP design
#' @export
assemble_kasp_assay <- function(genome, snp, variants = NULL, flank = 150L,
                                strands = c("+", "-"),
                                min_gc = 0.35, max_gc = 0.70,
                                tm_range = KASP_LIMITS$as_tm, clean_3p = 5L,
                                max_footprint_variants = Inf,
                                min_disc = 0.5, specificity_candidates = 10L,
                                pair_screen = 3L,
                                linked = NULL, haplotypes = NULL, run_blast = FALSE,
                                cores = getOption("tidykasp.cores", 1L), BPPARAM = NULL) {
  if (is.null(BPPARAM)) BPPARAM <- kasp_bpparam(cores)
  # one-line install UX: a missing tool is a hint (once per session), not a silent NA
  if (!has_oligoarrayaux() && !isTRUE(.kasp_state$tool_hint_shown)) {
    .kasp_state$tool_hint_shown <- TRUE
    cli::cli_inform(c(
      gettext("OligoArrayAux was not found: thermodynamic and structural checks will be NA."),
      i = gettext("run tidykasp::kasp_install_tools() once to build it.")
    ))
  }
  # a target site may carry more than one ALT (a list-column `alt`, or an `alt2`
  # column): the *site class* must see them all (below), while the design builds
  # one AS primer - for the first ALT until the degenerate-primer engine lands.
  alts_in <- unique(c(
    as.character(unlist(snp$alt)),
    if ("alt2" %in% names(snp)) as.character(snp$alt2) else character(0)
  ))
  alts_in <- alts_in[!is.na(alts_in) & nzchar(alts_in)]
  # a purely SNP row with >1 ALT is a degenerate site: one ALT primer whose 3'
  # base is the IUPAC consensus covers them all, so every ALT goes to the engine.
  # A mixed or non-SNP ALT keeps the single-allele path (first ALT) for now.
  snp_alts <- unique(alts_in[nchar(alts_in) == 1L & grepl("^[ACGT]$", alts_in)])
  degenerate <- length(snp_alts) > 1L && length(snp_alts) == length(alts_in)
  alt_pass <- if (degenerate) snp_alts else alts_in[1L]
  if (length(alts_in) > 1L) snp$alt <- alts_in[1L]
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  hap <- NULL
  if (!is.null(haplotypes)) {
    hap <- haplotype_variants(haplotypes, snp)
    if (is.null(variants)) variants <- hap$variants
    if (is.null(linked)) linked <- hap$linked
  }
  vpos <- if (is.null(variants)) NULL else variants$pos[variants$chrom == snp$chrom]
  # optional per-variant weights (population allele frequency); absent -> all 1
  vweight <- if (!is.null(variants) && "pop_freq" %in% names(variants)) {
    variants$pop_freq[variants$chrom == snp$chrom]
  } else {
    NULL
  }
  off <- snp$pos - flank # template coordinate -> genomic coordinate
  valt <- NULL
  if (!is.null(variants)) {
    vsub <- variants[variants$chrom == snp$chrom, , drop = FALSE]
    if ("alt" %in% names(vsub)) {
      valt <- as.character(vsub$alt)
    } else if ("ref_genome" %in% names(vsub)) {
      # Ensembl shape: alleles as a vector; the ALT allele is the one that is not
      # the (genome) reference base at that position
      p <- vsub$pos
      bi <- p - off + 1L
      ok <- !is.na(bi) & bi >= 1L & bi <= length(tmpl$seq)
      base <- rep(NA_character_, nrow(vsub))
      base[ok] <- substring(as.character(tmpl$seq), bi[ok], bi[ok])
      valt <- purrr::map2_chr(vsub$ref_genome, base, function(al, b) {
        a <- setdiff(as.character(al), b)
        a <- a[nzchar(a) & grepl("^[ACGTNacgtn]+$", a)]
        if (length(a)) a[[1L]] else NA_character_
      })
    }
  }
  # DECIPHER (via OligoArrayAux) is what turns a footprint variant into an
  # expected dropout; without it we fall back to the plain weighted count
  use_eff <- has_oligoarrayaux()
  ref_len <- nchar(snp$ref)
  alt_len <- nchar(snp$alt)
  is_indel <- ref_len != alt_len
  vblock <- seq.int(snp$pos, snp$pos + max(ref_len, alt_len) - 1L) # the target's own bases

  # one service for both critical zones (AS 3' window, COM footprint): the
  # expected dropout over the design template's variants; NA if DECIPHER chokes
  dropout_at <- function(pr, a, b, forward, cnt, exclude = vblock) {
    zone_dropout(pr, a, b, forward, cnt, tmpl, off, vpos, exclude, valt, vweight)
  }
  safe_dropout <- purrr::possibly(dropout_at, otherwise = NA_real_)
  # pool every candidate zone's DECIPHER pairs into one call; if that batch
  # fails, fall back to the per-candidate path (same NA semantics)
  batch_or_safe <- function(primers, starts, ends, forwards, nvars, exclude) {
    if (!use_eff) {
      return(rep(NA_real_, length(primers)))
    }
    tryCatch(
      zone_dropouts(primers, starts, ends, forwards, nvars, tmpl, off, vpos, exclude, valt, vweight,
        cores = kasp_ncores(BPPARAM)
      ),
      error = function(e) {
        purrr::pmap_dbl(
          list(primers, starts, ends, forwards, nvars),
          function(pr, a, b, fw, cnt) safe_dropout(pr, a, b, fw, cnt, exclude = exclude)
        )
      }
    )
  }

  # co-segregating substitutions to bake into the ALT primer (see design_as_primers):
  # positions in `linked$linked` that ride with `snp$pos`, mapped to template coords
  linked_off <- NULL
  linked_skipped <- 0L
  if (!is.null(linked) && nrow(linked) && all(c("target", "linked") %in% names(linked))) {
    lnk <- linked[linked$target == snp$pos, , drop = FALSE]
    # a grouped (per-haplotype) linked table bakes only the positions shared by
    # every ALT haplotype - a single ALT primer must bind them all; a flat
    # (target, linked) table is taken as-is
    if ("shared" %in% names(lnk)) lnk <- lnk[!is.na(lnk$shared) & lnk$shared, , drop = FALSE]
    lp <- unique(lnk$linked)
    li <- match(lp, vpos)
    ok <- !is.na(li) & !is.na(valt[li])
    if (any(ok)) {
      offset <- lp[ok] - off + 1L
      # the full REF allele (1 base for a substitution, longer for an indel) so
      # the ALT window can be rebuilt length-aware
      ref_allele <- if (!is.null(vsub) && "ref" %in% names(vsub)) {
        as.character(vsub$ref[li[ok]])
      } else {
        substring(as.character(tmpl$seq), offset, offset)
      }
      cand <- tibble::tibble(offset = offset, ref = ref_allele, alt = valt[li[ok]])
      cand <- cand[
        cand$offset >= 1L & (cand$offset + nchar(cand$ref) - 1L) <= length(tmpl$seq), ,
        drop = FALSE
      ]
      if (nrow(cand)) linked_off <- cand
    }
    linked_skipped <- length(lp) - (if (is.null(linked_off)) 0L else nrow(linked_off))
  }

  # co-segregating positions ride with the ALT haplotype: the AS primer is built
  # to match them (REF matches the REF base), so they are not "other variants" in
  # the AS zone - a linked variant in the AS 3' window is not a specificity
  # threat.  Under COM they still hurt (COM is built on the REF template).
  link_pos <- if (is.null(linked_off)) integer(0) else snp$pos - tmpl$snp_i + linked_off$offset
  as_exclude <- c(vblock, link_pos)

  # all primer geometry comes from the shared helpers (plus-strand coordinates);
  # the REF AS primer's 3'-terminal base is the variant anchor for its strand
  anchor_of <- function(strand) as_anchor(snp$pos, ref_len, strand)
  com_iv <- function(strand, L, gap) primer_interval(anchor_of(strand), strand, L, gap = gap, role = "COM")

  # candidate AS primers, unified to one schema for either variant kind:
  #   ref_primer/alt_primer, ref_on/ref_off/alt_on/alt_off, discrimination,
  #   tm, gc, as_length, strand, and genomic footprints as_* and as_alt_*
  if (!is_indel) {
    raw_as <- design_as_primers(
      tmpl$seq, tmpl$snp_i, snp$ref, alt_pass,
      strands = strands, linked = linked_off,
      gc_range = c(min_gc, max_gc), tm_range = tm_range, BPPARAM = BPPARAM
    )
    if (nrow(raw_as) == 0L) {
      return(kasp_doomed(gettext("No AS primer passed the GC/Tm filter."), "as", snp, tmpl))
    }
    # the GC/Tm window is already applied inside design_as_primers(); only the
    # sequence-motif filters are left for here
    as_cand <- raw_as |>
      dplyr::filter(!has_run(.data$ref_primer), !has_repeat(.data$ref_primer)) |>
      dplyr::mutate(
        as_start = .data$tmpl_start + off,
        as_end = .data$tmpl_end + off,
        as_alt_start = .data$as_start,
        as_alt_end = .data$as_end,
        wobble = wobble_mismatch(snp$ref, snp$alt)
      )
    if (nrow(as_cand) == 0L) {
      return(kasp_doomed(gettext("No AS primer passed the homopolymer/SSR filter."), "as", snp, tmpl))
    }
  } else {
    ind <- design_as_primers_indel(
      tmpl$seq, tmpl$snp_i, snp$ref, snp$alt,
      strands = strands, gc_range = c(min_gc, max_gc), tm_range = tm_range, BPPARAM = BPPARAM
    )
    if (nrow(ind) == 0L) {
      return(kasp_doomed(gettext("No indel AS primer passed the GC/Tm filter."), "as", snp, tmpl))
    }
    ind <- dplyr::filter(ind, !has_run(.data$primer), !has_repeat(.data$primer))
    if (nrow(ind) == 0L) {
      return(kasp_doomed(gettext("No indel AS primer passed the homopolymer/SSR filter."), "as", snp, tmpl))
    }
    # alt-template coordinate -> REF/genomic coordinate (block clamps to REF)
    alt_to_ref <- function(x) {
      dplyr::case_when(
        x < tmpl$snp_i ~ x,
        x <= tmpl$snp_i + alt_len - 1L ~ pmin(x, tmpl$snp_i + ref_len - 1L),
        TRUE ~ x - (alt_len - ref_len)
      )
    }
    top <- function(a) {
      ind |>
        dplyr::filter(.data$allele == a) |>
        dplyr::group_by(.data$strand) |>
        dplyr::slice_max(.data$discrimination, n = 1L, with_ties = FALSE) |>
        dplyr::ungroup()
    }
    ref_p <- top("ref")
    alt_p <- top("alt")
    as_cand <- dplyr::inner_join(ref_p, alt_p, by = "strand", suffix = c("_ref", "_alt")) |>
      dplyr::transmute(
        strand = .data$strand,
        as_length = .data$as_length_ref,
        ref_primer = .data$primer_ref,
        alt_primer = .data$primer_alt,
        alt_degen = .data$primer_alt, # length-changing ALT: no IUPAC consensus
        ref_on = .data$on_eff_ref, ref_off = .data$off_eff_ref,
        alt_on = .data$on_eff_alt, alt_off = .data$off_eff_alt,
        discrimination = (.data$discrimination_ref + .data$discrimination_alt) / 2,
        tm = .data$tm_ref, gc = .data$gc_ref,
        dg_3p = terminal_dg(.data$primer_ref),
        wobble = NA, # undefined for length-changing alleles
        as_start = .data$as_start_ref + off,
        as_end = .data$as_end_ref + off,
        as_alt_start = alt_to_ref(.data$as_start_alt) + off,
        as_alt_end = alt_to_ref(.data$as_end_alt) + off
      )
  }

  # clean-3' window and footprint variant load come from the stored interval
  as_cand <- as_cand |>
    dplyr::mutate(
      win3p_start = dplyr::if_else(.data$strand == "+", .data$as_end - clean_3p + 1L, .data$as_start + 1L),
      win3p_end = dplyr::if_else(.data$strand == "+", .data$as_end - 1L, .data$as_start + clean_3p - 1L),
      n_3p = purrr::map2_int(.data$win3p_start, .data$win3p_end, function(a, b) {
        count_variants_in(IRanges::IRanges(a, b), vpos, as_exclude)
      }),
      n_var = purrr::map2_int(.data$as_start, .data$as_end, function(a, b) {
        count_variants_in(IRanges::IRanges(a, b), vpos, as_exclude)
      }),
      w_var = purrr::map2_dbl(.data$as_start, .data$as_end, function(a, b) {
        variant_load_in(IRanges::IRanges(a, b), vpos, vweight, as_exclude)
      }),
      # expected dropout of the AS 3' window (position- and identity-aware)
      as_dropout = batch_or_safe(
        .data$ref_primer, .data$win3p_start, .data$win3p_end,
        .data$strand == "+", .data$n_3p, as_exclude
      ),
      # expected dropout of the whole AS footprint: a rare mid-body variant is
      # cheap (Taq still extends from a matched 3' end), a 3'-terminal one costly
      as_foot_dropout = batch_or_safe(
        .data$ref_primer, .data$as_start, .data$as_end,
        .data$strand == "+", .data$n_var, as_exclude
      )
    )

  # "N-masking" as a hard constraint: optionally drop candidates sitting on
  # other catalogued variants (their footprint would be ambiguous)
  if (is.finite(max_footprint_variants)) {
    as_cand <- dplyr::filter(as_cand, .data$n_var <= max_footprint_variants)
  }
  if (nrow(as_cand) == 0L) {
    return(kasp_doomed(
      glue::glue(gettext(
        "No AS primer survived the footprint cap (max_footprint_variants = {max_footprint_variants})."
      )),
      "as", snp, tmpl, as_cand
    ))
  }

  # prioritise with the shared grader: gate the grave sins (a wobble mismatch, a
  # variant in the 3' window or the footprint), then rank lexicographically on the
  # margins.  The genome-wide 3'-seed screen is only a tie-break, so only the
  # leading placements are screened; an indexed genome (no in-memory DNAStringSet)
  # skips it.
  graded_as <- grade_as_candidates(as_cand, KASP_LIMITS, min_disc)
  graded_as$data$seed_hits <- NA_integer_
  if (specificity_candidates > 0L && !inherits(genome, "kasp_genome")) {
    short <- rank_candidates(graded_as$data, graded_as$tiers, setdiff(graded_as$ties, "seed_hits")) |>
      dplyr::slice_head(n = specificity_candidates)
    seeds <- rlang::set_names(
      substr(short$ref_primer, nchar(short$ref_primer) - SEED_K + 1L, nchar(short$ref_primer)),
      short$as_row
    )
    hits <- specificity_seed(seeds, genome,
      k = SEED_K, max_mismatch = SEED_MISMATCH,
      BPPARAM = BPPARAM
    )$hits
    graded_as$data$seed_hits[short$as_row] <- hits
  }
  pick_as <- select_with_ladder(
    graded_as$data, graded_as$gates, KASP_RELAX_AS,
    graded_as$tiers, graded_as$ties, graded_as$gate_margin, graded_as$gate_floor
  )
  if (is.null(pick_as$best)) {
    return(kasp_doomed(
      gettext("No AS primer survived even the relaxed constraints."), "as",
      snp, tmpl, tidy_as_candidates(rank_candidates(graded_as$data, graded_as$tiers, graded_as$ties))
    ))
  }
  best_as <- pick_as$best
  # the ranked candidate table for inspection, minus the grader's helper columns
  as_cand <- tidy_as_candidates(rank_candidates(graded_as$data, graded_as$tiers, graded_as$ties))

  if (!is.na(best_as$n_3p) && best_as$n_3p > 0L) {
    cli::cli_warn(c(
      gettext(
        "Target SNP {snp$chrom}:{snp$pos} has {best_as$n_3p} other variant(s) in the AS 3'-terminal {clean_3p} nt."
      ),
      i = gettext("Allele specificity on strand {best_as$strand} may not hold across accessions.")
    ))
  }
  if (!is.na(best_as$n_var) && best_as$n_var > 0L) {
    cli::cli_warn(c(
      gettext("AS primer overlaps {best_as$n_var} other variant(s) in its footprint."),
      i = gettext("No candidate with fewer passed the GC/Tm/discrimination filters.")
    ))
  }
  if (best_as$discrimination < min_disc) {
    cli::cli_warn(c(
      gettext(
        "Best AS primer discrimination = {round(best_as$discrimination, 2)} is below min_disc = {min_disc}."
      ),
      i = gettext("Consider widening `strands` or `tm_range`.")
    ))
  }
  # a design can hold *only* because co-segregating variants were excluded: they
  # ride with the ALT haplotype (baked into the ALT primer), so the AS 3' window
  # and footprint look clean without them.  That is valid only where the linked
  # haplotype is present - say so rather than pass in silence.
  linked_dependent <- FALSE
  if (length(link_pos)) {
    in_3p <- link_pos >= best_as$win3p_start & link_pos <= best_as$win3p_end
    in_foot <- link_pos >= best_as$as_start & link_pos <= best_as$as_end
    n_link <- sum(in_3p | in_foot)
    linked_dependent <- n_link > 0L
    if (linked_dependent) {
      cli::cli_warn(c(
        gettext(
          "Design at {snp$chrom}:{snp$pos} relies on linkage ({n_link} linked variant(s) in the AS zone)."
        ),
        i = gettext("The ALT primer is valid only in accessions carrying the linked haplotype.")
      ))
    }
  }

  # COM: design on the SAME strand as the winning AS primer so the amplicon
  # closes.  design_common_primers() works in TEMPLATE coordinates, so its
  # anchor is the template-local variant anchor (not the genomic one).
  common <- design_common_primers(tmpl$seq, as_anchor(tmpl$snp_i, ref_len, best_as$strand),
    strand = best_as$strand, tm_range = tm_range, n = 60L
  ) |>
    dplyr::mutate(
      com_homodimer = oligo_dg_vec(.data$com, .data$com, BPPARAM = BPPARAM),
      com_hairpin = oligo_dg_vec(.data$com, BPPARAM = BPPARAM)
    )

  as_tm_mean <- mean(best_as$tm)

  assays <- common |>
    dplyr::mutate(
      strand = best_as$strand,
      com_start = purrr::pmap_int(
        list(.data$strand, .data$len, .data$gap),
        function(st, L, g) IRanges::start(com_iv(st, L, g))
      ),
      com_end = purrr::pmap_int(
        list(.data$strand, .data$len, .data$gap),
        function(st, L, g) IRanges::end(com_iv(st, L, g))
      ),
      n_var = purrr::map2_int(.data$gap, .data$len, function(g, L) {
        count_variants_in(com_iv(best_as$strand, L, g), vpos, vblock)
      }),
      w_var = purrr::map2_dbl(.data$gap, .data$len, function(g, L) {
        variant_load_in(com_iv(best_as$strand, L, g), vpos, vweight, vblock)
      }),
      # expected dropout of the COM footprint (position- and identity-aware)
      com_dropout = batch_or_safe(
        .data$com, .data$com_start, .data$com_end,
        .data$strand != "+", .data$n_var, vblock
      ),
      as_ref = best_as$ref_primer, as_alt = best_as$alt_primer,
      as_alt_degen = best_as$alt_degen, as_dropout = best_as$as_dropout,
      as_foot_dropout = best_as$as_foot_dropout,
      discrimination = best_as$discrimination, as_gc = best_as$gc,
      as_len = best_as$as_length, as_n_var = best_as$n_var, as_n_3p = best_as$n_3p,
      # worst-case product over both alleles (indels differ by delta)
      amplicon = pmax(
        pmax(best_as$as_end, .data$com_end) - pmin(best_as$as_start, .data$com_start) + 1L,
        pmax(best_as$as_alt_end, .data$com_end) - pmin(best_as$as_alt_start, .data$com_start) + 1L
      ),
      # store the geometry once, for the plots to read instead of recomputing
      as_start = best_as$as_start, as_end = best_as$as_end,
      as_alt_start = best_as$as_alt_start, as_alt_end = best_as$as_alt_end,
      win3p_start = best_as$win3p_start, win3p_end = best_as$win3p_end
    )

  # choose among the pairings: gate the grave sins, then rank lexicographically
  # on tier margins, relaxing gates one at a time only when nothing clears them
  graded <- grade_assay_pairs(assays, KASP_LIMITS, as_tm_mean)
  pick <- select_with_ladder(
    graded$data, graded$gates, KASP_RELAX_PAIR, graded$tiers, graded$ties,
    graded$gate_margin, graded$gate_floor
  )
  if (is.null(pick$best)) {
    return(kasp_doomed(
      gettext("No candidate assay survived even the relaxed constraints."), "pair", snp, tmpl, as_cand
    ))
  }
  best <- pick$best
  # keep the chosen pairing first: the rest of the package (checks, plots) reads
  # assay[1, ] as "the" assay, so the ranking must be reflected in the row order
  ord <- c(pick$best$.row, setdiff(seq_len(nrow(assays)), pick$best$.row))
  assays <- assays[ord, , drop = FALSE]
  assays$rank <- seq_len(nrow(assays))
  warn_relaxed(pick$violated)

  # The ladder ranks on genome-free margins, so a second in-silico PCR product, an
  # off-target exact primer hit or an inverted repeat can still sit on the winner.
  # Those cost a genome/PCR pass, so only the leading `pair_screen` pairings are
  # screened; the first clean one wins (the ranked winner is kept when none is).
  if (!is.null(genome) && pair_screen > 1L && nrow(assays) > 1L) {
    n_screen <- min(as.integer(pair_screen), nrow(assays))
    # the AS primers are shared by every pairing, so their exact hits are scanned
    # once; only the COM primer (which varies) is rescanned per candidate
    as_primers <- c(AS_ref = assays$as_ref[1L], AS_alt = assays$as_alt[1L])
    as_hits <- tryCatch(primer_hits(as_primers, genome, BPPARAM = BPPARAM), error = function(e) NULL)
    as_bad <- !is.null(as_hits) && nrow(as_hits) && max(table(as_hits$primer)) > 1L
    screen_pair <- function(i) {
      row <- assays[i, ]
      if (as_bad || any(has_palindrome(c(as_primers, row$com)))) {
        return(FALSE)
      }
      com_hits <- tryCatch(
        primer_hits(c(COM = row$com), genome, BPPARAM = BPPARAM),
        error = function(e) NULL
      )
      if (!is.null(com_hits) && nrow(com_hits) && max(table(com_hits$primer)) > 1L) {
        return(FALSE)
      }
      prod <- tryCatch(
        assay_pcr_products(list(assay = row, template = tmpl, snp = snp), genome, BPPARAM = BPPARAM),
        error = function(e) NULL
      )
      if (is.null(prod)) {
        return(TRUE) # unknown -> do not block
      }
      cnt <- assay_amplicon_counts(prod)
      ref <- cnt$n_products[cnt$primer == "AS_ref"]
      alt <- cnt$n_products[cnt$primer == "AS_alt"]
      length(ref) == 1L && length(alt) == 1L && ref == 1L && alt == 0L
    }
    clean <- purrr::map_lgl(seq_len(n_screen), screen_pair)
    win <- which(clean)[1L]
    if (!is.na(win) && win > 1L) {
      cli::cli_warn(c(
        gettext(
          "Top-ranked pairing failed genome-level checks; kept the next clean one ({win} of {n_screen} screened)."
        ),
        i = gettext("Ranking is genome-free; off-target hits, in-silico PCR and inverted repeats are screened here.")
      ))
      ord <- c(win, setdiff(seq_len(nrow(assays)), win))
      assays <- assays[ord, , drop = FALSE]
      assays$rank <- seq_len(nrow(assays))
    }
  }
  best <- assays[1L, ]
  primers <- c(AS_ref = best$as_ref, AS_alt = best$as_alt, COM = best$com)

  # specificity uses the template-matched allele for the AS primers (no
  # deliberate mismatch): that answers "is the target site unique?".
  # make_as_primer() works in TEMPLATE coordinates, hence the template anchor.
  tmpl_anchor <- as_anchor(tmpl$snp_i, ref_len, best$strand)
  seed_primers <- c(
    AS_ref = make_as_primer(tmpl$seq, tmpl_anchor, best$as_len, snp$ref, strand = best$strand),
    AS_alt = make_as_primer(tmpl$seq, tmpl_anchor, best$as_len, snp$alt, strand = best$strand),
    COM    = best$com
  )

  # in-silico PCR: the REF primer must give one product, the ALT primer none.
  # The products are computed once and reused for the amplicon counts, the
  # ideogram below and (via `out$ideogram`) the plots, so the genome is scanned
  # a single time per assay instead of three.
  stub <- list(assay = best, template = tmpl, snp = snp)
  pcr_products <- tryCatch(
    assay_pcr_products(stub, genome, BPPARAM = BPPARAM),
    error = function(e) {
      cli::cli_warn(gettext("in-silico PCR failed: {conditionMessage(e)}"))
      NULL
    }
  )
  amplicons <- if (is.null(pcr_products)) NULL else assay_amplicon_counts(pcr_products)

  # whole-primer genome scan (only when requested): an attached FM-index is
  # queried first (index-based, multithreaded), otherwise the BLAST database
  blast <- NULL
  if (isTRUE(run_blast)) {
    blast <- if (inherits(genome, "kasp_genome") && !is.null(genome$fm_index)) {
      tryCatch(specificity_fm(primers, genome$fm_index),
        error = function(e) {
          cli::cli_warn(gettext("FM scan failed: {conditionMessage(e)}"))
          NULL
        }
      )
    } else {
      db <- tryCatch(build_blast_db(genome), error = function(e) {
        cli::cli_warn(gettext("BLAST db unavailable: {conditionMessage(e)}"))
        NULL
      })
      blast_primers(primers, db)
    }
  }

  out <- list(
    snp = snp, template = tmpl, strand = best$strand, assay = assays,
    as_candidates = as_cand,
    efficiency = best_as |> dplyr::select("ref_on", "ref_off", "alt_on", "alt_off", "discrimination"),
    thermo = purrr::imap(primers, ~ dplyr::mutate(primer_thermo(.x, BPPARAM = BPPARAM), name = .y)) |>
      purrr::list_rbind(),
    complementarity = dimer_thermo(primers, BPPARAM = BPPARAM),
    specificity = specificity_of(seed_primers, genome, BPPARAM = BPPARAM),
    amplicons = amplicons,
    selection = list(
      relaxed = pick$relaxed, violated = pick$violated, as = pick_as$violated,
      linked = linked_off, linked_skipped = linked_skipped,
      linked_dependent = linked_dependent
    ),
    blast = blast,
    haplotypes = if (is.null(hap)) NULL else hap$haplotypes,
    haplotype_hits = if (is.null(hap)) {
      NULL
    } else {
      haplotype_hits(haplotypes, snp, best$as_start, best$as_end, hap$alt)
    },
    tails = c(
      AS_ref = paste0(KASP_TAIL_FAM, best$as_ref),
      AS_alt = paste0(KASP_TAIL_HEX, best$as_alt),
      COM = best$com
    )
  )
  # precompute the whole-genome ideogram inputs once, so plotting only reads the
  # object instead of rescanning the genome
  out$ideogram <- tryCatch(assay_ideogram_data(out, genome, BPPARAM = BPPARAM, products = pcr_products),
    error = function(e) {
      cli::cli_warn(gettext("genome ideogram failed: {conditionMessage(e)}"))
      NULL
    }
  )
  # the target site's KASP class (biallelic / degenerate / complementary / indel),
  # from every ALT seen at the position in the variant table
  site_alts <- alts_in
  if (!is.null(variants) && all(c("chrom", "pos") %in% names(variants))) {
    rows <- variants[variants$chrom == snp$chrom & variants$pos == snp$pos, , drop = FALSE]
    if (nrow(rows)) {
      a <- if ("ref_genome" %in% names(rows)) {
        # `ref_genome` is a LIST column (one allele vector per variant): unlist
        # first - `as.character()` on a list would stringify the whole vector
        as.character(unlist(rows$ref_genome, use.names = FALSE))
      } else if ("alt" %in% names(rows)) {
        as.character(rows$alt)
      } else {
        character(0)
      }
      a <- a[nzchar(a) & grepl("^[ACGT]+$", a)]
      if (length(a)) site_alts <- unique(c(site_alts, a))
    }
  }
  out$site <- kasp_site_class(snp$ref, site_alts)
  # how many of the site's ALTs the single ALT primer actually covers: all of them
  # for a degenerate all-SNP site (one IUPAC primer), otherwise just the first
  out$site$alt_covered <- if (degenerate) out$site$n_alt[[1L]] else 1L
  # a multi-allelic site is reported as such.  When it is all-SNP it is designed
  # as a degenerate primer covering every ALT; otherwise the ALT primer still
  # covers the first ALT only - say so, do not pretend.
  if (out$site$n_alt[1L] > 1L) {
    if (degenerate) {
      cli::cli_inform(c(
        gettext(
          "Target {snp$chrom}:{snp$pos} is multi-allelic ({out$site$n_alt[1L]} ALTs, class {out$site$class[1L]})."
        ),
        i = gettext(
          "The ALT primer is degenerate (3' base {out$site$code[1L]}): it covers every ALT but cannot tell them apart."
        )
      ))
    } else {
      cli::cli_warn(c(
        gettext(
          "Target {snp$chrom}:{snp$pos} is multi-allelic ({out$site$n_alt[1L]} ALTs, class {out$site$class[1L]})."
        ),
        i = gettext(
          "The ALT primer is built for {snp$alt} only; the other ALTs are not covered."
        )
      ))
    }
  }
  # evaluate every registered acceptance check once, so plots/reports read the
  # results instead of recomputing them
  out$evaluation <- tryCatch(kasp_design_checks(out), error = function(e) NULL)
  # structural / specificity checks that the selection gates do not cover (they
  # need the genome or a PCR pass, computed only after the pair is chosen):
  # rather than return an ugly assay in silence, name the failures
  if (!is.null(out$evaluation)) {
    hit <- !out$evaluation$pass & out$evaluation$criterion %in% KASP_UNGATED_CHECKS
    bad <- out$evaluation$criterion[hit]
    if (length(bad)) {
      cli::cli_warn(c(
        gettext(
          "The chosen assay fails structural/specificity checks: {paste(bad, collapse = '; ')}."
        ),
        i = gettext("Selection does not gate these yet - inspect assay$evaluation and consider a trade-off.")
      ))
    }
  }
  # the selection policy as a graph, with this assay's runtime status stamped on
  out$decision <- tryCatch(kasp_decision_graph(out),
    error = function(e) NULL
  )
  structure(out, class = "kasp_assay")
}
