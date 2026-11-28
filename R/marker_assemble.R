# =============================================================================
#  marker_assemble.R  --  assemble and score one complete assay
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
marker_doomed <- function(reason, stage, snp, template, as_candidates = NULL) {
  structure(
    list(snp = snp, template = template, reason = reason, stage = stage, as_candidates = as_candidates),
    class = c("marker_doomed", "marker")
  )
}

#' @keywords internal
#' @noRd
is_kasp_doomed <- function(x) inherits(x, "marker_doomed")

#' Drop the grader's helper columns from an AS candidate table
#' @keywords internal
#' @noRd
tidy_as_candidates <- function(x) {
  drop <- grep("^\\.neg_|^m_|^gate_|^tier|^as_row$", names(x), value = TRUE)
  x[, setdiff(names(x), drop), drop = FALSE]
}

# ---- assemble_kasp_assay() stages ------------------------------------------
# The assembly below is one straight pipeline; these helpers are its stages, so
# each is readable (and testable) on its own without changing the behaviour.

#' Resolve the ALT allele(s) of a target row
#'
#' A target may carry one ALT or several (a list-column `alt`, or an `alt2`
#' column).  A purely SNP row with >1 ALT is *degenerate*: one ALT primer whose
#' 3' base is the IUPAC consensus covers them all, so every ALT goes to the
#' engine.  A mixed or non-SNP ALT keeps the single-allele path (first ALT).
#' @keywords internal
#' @noRd
marker_snp_alts <- function(snp) {
  alts_in <- unique(c(
    as.character(unlist(snp$alt)),
    if ("alt2" %in% names(snp)) as.character(snp$alt2) else character(0)
  ))
  alts_in <- alts_in[!is.na(alts_in) & nzchar(alts_in)]
  # An "ALT" that is the REF base is not an ALT.  If it were kept, REF would enter
  # the degenerate 3' code and the "ALT" primer would also extend the reference
  # template - i.e. discriminate nothing - so it is dropped from the consensus
  # set (marker_site_class() excludes it the same way).  alts_in itself is left
  # untouched, so the single-ALT fallback and `is_indel` are unchanged.
  snp_alts <- setdiff(
    unique(alts_in[nchar(alts_in) == 1L & grepl("^[ACGT]$", alts_in)]),
    snp$ref
  )
  degenerate <- length(snp_alts) > 1L && length(snp_alts) == length(alts_in)
  alt_pass <- if (degenerate) snp_alts else alts_in[1L]
  if (length(alts_in) > 1L) snp$alt <- alts_in[1L]
  list(alts_in = alts_in, degenerate = degenerate, alt_pass = alt_pass, snp = snp)
}

#' Per-position ALT alleles (and optional weights) of the variant table
#' @keywords internal
#' @noRd
marker_variant_alleles <- function(variants, snp, seq, off) {
  vpos <- if (is.null(variants)) NULL else variants$pos[variants$chrom == snp$chrom]
  # optional per-variant weights (population allele frequency); absent -> all 1,
  # a missing value -> the same whole-variant weight, never a zero cost (ADR 0086)
  vweight <- if (!is.null(variants) && "pop_freq" %in% names(variants)) {
    check_variant_weights(variants$pop_freq[variants$chrom == snp$chrom])
  } else {
    NULL
  }
  vsub <- NULL
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
      ok <- !is.na(bi) & bi >= 1L & bi <= length(seq)
      base <- rep(NA_character_, nrow(vsub))
      if (any(ok)) base[ok] <- substring(as.character(seq), bi[ok], bi[ok])
      valt <- purrr::map2_chr(vsub$ref_genome, base, function(al, b) {
        a <- setdiff(as.character(al), b)
        a <- a[nzchar(a) & grepl("^[ACGTNacgtn]+$", a)]
        if (length(a)) a[[1L]] else NA_character_
      })
    }
  }
  list(vpos = vpos, vweight = vweight, vsub = vsub, valt = valt)
}

#' Co-segregating substitutions to bake into the ALT primer
#'
#' Maps the linked positions that ride with the target to template coordinates
#' and returns them plus the excluded positions of the AS zone (the target's own
#' bases and the linked ones - a linked variant in the AS zone is not a threat).
#' @keywords internal
#' @noRd
marker_linked_offsets <- function(linked, snp, vpos, valt, vsub, seq, off, snp_i, vblock) {
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
        substring(as.character(seq), offset, offset)
      }
      cand <- tibble::tibble(offset = offset, ref = ref_allele, alt = valt[li[ok]])
      cand <- cand[
        cand$offset >= 1L & (cand$offset + nchar(cand$ref) - 1L) <= length(seq), ,
        drop = FALSE
      ]
      if (nrow(cand)) linked_off <- cand
    }
    linked_skipped <- length(lp) - (if (is.null(linked_off)) 0L else nrow(linked_off))
  }
  link_pos <- if (is.null(linked_off)) integer(0) else snp$pos - snp_i + linked_off$offset
  list(
    linked_off = linked_off, linked_skipped = linked_skipped,
    link_pos = link_pos, as_exclude = c(vblock, link_pos)
  )
}

#' Build the zone-dropout service for one template
#'
#' Returns `batch_or_safe(primers, starts, ends, forwards, nvars, exclude)`: it
#' pools the candidates' DECIPHER pairs into one call and falls back to a
#' per-candidate pass on error (NA semantics preserved).  With no OligoArrayAux
#' every zone is `NA`.
#' @keywords internal
#' @noRd
marker_dropout_batch <- function(tmpl, off, vpos, valt, vweight, use_eff, BPPARAM) {
  safe_dropout <- purrr::possibly(
    function(pr, a, b, forward, cnt, exclude) {
      zone_dropout(pr, a, b, forward, cnt, tmpl, off, vpos, exclude, valt, vweight)
    },
    otherwise = NA_real_
  )
  function(primers, starts, ends, forwards, nvars, exclude) {
    if (!use_eff) {
      return(rep(NA_real_, length(primers)))
    }
    tryCatch(
      zone_dropouts(primers, starts, ends, forwards, nvars, tmpl, off, vpos, exclude, valt, vweight,
        cores = marker_ncores(BPPARAM)
      ),
      error = function(e) {
        purrr::pmap_dbl(
          list(primers, starts, ends, forwards, nvars),
          function(pr, a, b, fw, cnt) safe_dropout(pr, a, b, fw, cnt, exclude = exclude)
        )
      }
    )
  }
}

#' Candidate AS primers, unified to one schema for either variant kind
#'
#' Returns the table (`ref_primer`/`alt_primer`, efficiencies, `discrimination`,
#' `tm`, `gc`, `as_length`, `strand`, genomic footprints `as_*`/`as_alt_*`), or a
#' `marker_doomed` object when no candidate clears the filters.
#' @param strands which AS-primer orientations to search (see [design_as_primers()]).
#' @param conditions reaction conditions passed to [design_as_primers()]
#'   (`NULL` keeps the historical buffer).
#' @keywords internal
#' @noRd
marker_as_candidates <- function(tmpl, snp, alt_pass, is_indel, off, strands, tm_range,
                                 min_gc, max_gc, linked_off, ref_len, alt_len, BPPARAM,
                                 conditions = NULL, variant_offset = 0L,
                                 as_lengths = KASP_AS_LENGTHS) {
  if (!is_indel) {
    # cheap pre-gate: when no length can clear the GC/Tm window the whole grid is
    # built and thrown away, so answer the doomed case before paying for it.  The
    # reason/stage are the ones the empty table produced below, unchanged.
    if (!as_grid_can_pass(
      tmpl$seq, tmpl$snp_i, snp$ref, strands, as_lengths,
      gc_range = c(min_gc, max_gc), tm_range = tm_range, conditions = conditions,
      variant_offset = variant_offset
    )) {
      return(marker_doomed(gettext("No AS primer passed the GC/Tm filter."), "as", snp, tmpl))
    }
    raw_as <- design_as_primers(
      tmpl$seq, tmpl$snp_i, snp$ref, alt_pass,
      strands = strands, linked = linked_off, lengths = as_lengths,
      gc_range = c(min_gc, max_gc), tm_range = tm_range, BPPARAM = BPPARAM,
      conditions = conditions, variant_offset = variant_offset
    )
    if (nrow(raw_as) == 0L) {
      return(marker_doomed(gettext("No AS primer passed the GC/Tm filter."), "as", snp, tmpl))
    }
    # the GC/Tm window is already applied inside design_as_primers(); only the
    # sequence-motif filters are left for here
    as_cand <- raw_as |>
      dplyr::filter(!has_run(.data$ref_primer), !has_repeat(.data$ref_primer)) |>
      dplyr::mutate(
        as_start = .data$tmpl_start + off - 1L,
        as_end = .data$tmpl_end + off - 1L,
        as_alt_start = .data$as_start,
        as_alt_end = .data$as_end,
        # the wobble is a property of the 3'-terminal allele pair, so at a
        # degenerate (multi-allelic) site it is the WORST case over every ALT:
        # if any ALT the primer must bind forms a weak G:T / A:C terminal pair
        # with REF, the ALT primer discriminates poorly there.
        wobble = any(wobble_mismatch(snp$ref, alt_pass))
      )
    if (nrow(as_cand) == 0L) {
      return(marker_doomed(gettext("No AS primer passed the homopolymer/SSR filter."), "as", snp, tmpl))
    }
    as_cand <- marker_degenerate_alt_gate(as_cand, alt_pass, snp, tmpl)
    if (is_kasp_doomed(as_cand)) {
      return(as_cand)
    }
  } else {
    ind <- design_as_primers_indel(
      tmpl$seq, tmpl$snp_i, snp$ref, snp$alt,
      strands = strands, gc_range = c(min_gc, max_gc), tm_range = tm_range, BPPARAM = BPPARAM,
      conditions = conditions
    )
    if (nrow(ind) == 0L) {
      return(marker_doomed(gettext("No indel AS primer passed the GC/Tm filter."), "as", snp, tmpl))
    }
    ind <- dplyr::filter(ind, !has_run(.data$primer), !has_repeat(.data$primer))
    if (nrow(ind) == 0L) {
      return(marker_doomed(gettext("No indel AS primer passed the homopolymer/SSR filter."), "as", snp, tmpl))
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
        # A length change has no substitution wobble, so the substitution flag is
        # only defined when the two AS primers are the SAME length: then the
        # junction differs by a single 3'-terminal base (a substitution-like
        # junction) and the usual G:T / A:C rule applies.  A length-differing pair
        # mismatches by a bulge, not a wobble, so it is not flagged.
        # `wobble_mismatch()` reads plus-strand alleles, so a '-' primer's 3' base
        # is normalised to its plus-strand allele first.
        wobble = dplyr::if_else(
          nchar(.data$primer_ref) == nchar(.data$primer_alt),
          wobble_mismatch(
            plus_3p_base(.data$primer_alt, .data$strand),
            plus_3p_base(.data$primer_ref, .data$strand)
          ),
          FALSE
        ),
        as_start = .data$as_start_ref + off - 1L,
        as_end = .data$as_end_ref + off - 1L,
        as_alt_start = alt_to_ref(.data$as_start_alt) + off - 1L,
        as_alt_end = alt_to_ref(.data$as_end_alt) + off - 1L
      )
  }
  as_cand
}

#' Doom a multi-allelic design whose ALT primer is not the ALT-set consensus
#'
#' At a degenerate (all-SNP multi-allelic) site one physical ALT primer must bind
#' every ALT, so its 3'-terminal base has to be the IUPAC code over the whole ALT
#' set (`iupac_code()`) and must not also match REF (an "ALT" primer that extends
#' the reference template discriminates nothing).  The engine carries that
#' consensus in `alt_degen`; this verifies it instead of trusting it, so a site
#' whose ALT set no single primer can cover is *doomed* rather than returned as a
#' single-ALT primer that claims full coverage.
#' @param as_cand the AS candidate table (see `marker_as_candidates()`).
#' @param alt_pass the ALT allele(s) the primer must bind.
#' @param snp,tmpl the target and template (for the doomed record).
#' @return `as_cand`, or a `marker_doomed` object when no candidate carries the code.
#' @keywords internal
#' @noRd
marker_degenerate_alt_gate <- function(as_cand, alt_pass, snp, tmpl) {
  if (length(alt_pass) < 2L) {
    return(as_cand) # a biallelic site has one ALT: nothing to cover
  }
  code <- iupac_code(alt_pass)
  members <- if (is.na(code)) NA_character_ else Biostrings::IUPAC_CODE_MAP[[code]]
  last <- substr(as_cand$alt_degen, nchar(as_cand$alt_degen), nchar(as_cand$alt_degen))
  ok <- !is.na(members) & last == code & !grepl(snp$ref, members, fixed = TRUE)
  if (!any(ok)) {
    return(marker_doomed(
      gettext("No ALT primer's 3' code is the consensus of every ALT."),
      "as", snp, tmpl, as_cand
    ))
  }
  as_cand[ok, , drop = FALSE]
}
#'
#' Adds `n_3p`/`n_var`/`w_var` and the expected `as_dropout`/`as_foot_dropout`,
#' then applies the optional hard `max_footprint_variants` cap.  Returns the
#' enriched table or a `marker_doomed` object when the cap removes everything.
#' @keywords internal
#' @noRd
marker_as_zone_metrics <- function(as_cand, clean_3p, vpos, vweight, as_exclude,
                                   batch_or_safe, max_footprint_variants, snp, tmpl,
                                   variant_offset = 0L) {
  # clean-3' window and footprint variant load come from the stored interval.
  # The window sits just inside the 3' terminus; with `variant_offset` bases
  # between the variant and that terminus (ASQ), it starts past them.
  as_cand <- as_cand |>
    dplyr::mutate(
      win3p_start = dplyr::if_else(
        .data$strand == "+",
        .data$as_end - clean_3p + 1L - variant_offset,
        .data$as_start + 1L + variant_offset
      ),
      win3p_end = dplyr::if_else(
        .data$strand == "+",
        .data$as_end - 1L - variant_offset,
        .data$as_start + clean_3p - 1L + variant_offset
      ),
      n_3p = purrr::map2_int(.data$win3p_start, .data$win3p_end, function(a, b) {
        # clean_3p = 0 collapses the window; an empty range would abort IRanges
        if (is.na(a) || is.na(b) || a > b) {
          return(0L)
        }
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
    return(marker_doomed(
      glue::glue(gettext(
        "No AS primer survived the footprint cap (max_footprint_variants = {max_footprint_variants})."
      )),
      "as", snp, tmpl, as_cand
    ))
  }
  as_cand
}

#' Grade and pick the AS placement (with the genome 3'-seed tie-break)
#'
#' Gates the grave sins, screens the leading placements for 3'-seed uniqueness,
#' then selects via the relaxation ladder.  Returns `best_as`, the raw `pick` and
#' the tidied candidate table, or a `marker_doomed` object.
#' @keywords internal
#' @noRd
marker_pick_as <- function(as_cand, genome, min_disc, specificity_candidates, BPPARAM, snp, tmpl,
                           limits = KASP_LIMITS, tm_target = AS_TM_TARGET) {
  graded_as <- grade_as_candidates(as_cand, limits, min_disc, tm_target = tm_target)
  graded_as$data$seed_hits <- NA_integer_
  if (specificity_candidates > 0L) {
    # Only the leading placements are seed-screened, so as_candidates carries a
    # PARTIAL column: NA means "not screened", never "not specific".  The ladder
    # sorts NA last, so an unscreened row can never outrank a screened one.
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
    return(marker_doomed(
      gettext("No AS primer survived even the relaxed constraints."), "as",
      snp, tmpl, tidy_as_candidates(rank_candidates(graded_as$data, graded_as$tiers, graded_as$ties))
    ))
  }
  list(
    best_as = pick_as$best, pick_as = pick_as,
    as_cand = tidy_as_candidates(rank_candidates(graded_as$data, graded_as$tiers, graded_as$ties))
  )
}

#' Warn about the chosen AS placement and report linkage dependence
#'
#' Names the soft failures (variants in the 3' window / footprint, discrimination
#' below `min_disc`) and the linkage dependence - a design that holds only because
#' co-segregating variants were excluded.  Returns `linked_dependent`.
#' @keywords internal
#' @noRd
marker_as_warnings <- function(best_as, snp, clean_3p, min_disc, link_pos) {
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
  linked_dependent
}
