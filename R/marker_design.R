# =============================================================================
#  marker_design.R  --  common and allele-specific primer design
# =============================================================================
# =============================================================================
#  3. Common (reverse) primer: a small library-scored candidate scan
# =============================================================================

#' Design candidate common primers
#'
#' Enumerates the partner primer and keeps the ones whose length, GC, Tm and
#' complexity look like a good KASP common primer.  `strand` is the orientation
#' of the *allele-specific* primer: the common primer must sit on the opposite
#' side of the SNP and point back at it, otherwise the amplicon would not close.
#'
#'   * `strand = "+"`: AS extends downstream, so the common primer is a reverse
#'     primer *downstream* of the SNP (sequence = reverse complement).
#'   * `strand = "-"`: AS extends upstream, so the common primer is a forward
#'     primer *upstream* of the SNP (sequence = the plus strand itself).
#'
#' Metrics come from packages (TmCalculator, Biostrings); enumeration/filtering
#' is ours.
#' @param window plus-strand template window (from [snp_template()]).
#' @param snp_i template coordinate of the variant anchor.
#' @param strand assay orientation.
#' @param lengths candidate primer lengths.
#' @param gaps distance from the anchor to the common primer.
#' @param target_tm preferred common-primer Tm.
#' @param min_gc,max_gc GC bounds.
#' @param tm_range accepted Tm window (C).
#' @param n number of candidates to return.
#' @param conditions reaction conditions for [primer_tm()] (`NULL` = the
#'   historical 50 mM Na / 37 C buffer).  Only the buffer fields of `conditions`
#'   enter the nearest-neighbour Tm, so no scheme changes a number today.
#' @return a tibble of candidates (`com`, `len`, `gap`, `start`, `end`, `com_gc`,
#'   `com_tm`, `run`, `ssr`, `rank`).
#' @examples
#' set.seed(1)
#' w <- Biostrings::DNAString(paste(sample(c("A", "C", "G", "T"), 160L, TRUE), collapse = ""))
#' design_common_primers(w, 80L, n = 3L)
#' @family KASP design
#' @export
design_common_primers <- function(window, snp_i, strand = "+", lengths = 18:25,
                                  gaps = 18:80, target_tm = COM_TM_TARGET,
                                  min_gc = 0.40, max_gc = 0.60,
                                  tm_range = c(57, 64), n = 5L, conditions = NULL) {
  wlen <- length(window)
  plus <- identical(strand, "+")
  grid <- tidyr::crossing(L = lengths, gap = gaps)

  # vectorised COM geometry (the `COM` branch of primer_interval()).  The
  # scalar-IRanges + one-row-tibble loop was a top profile entry: hundreds of
  # allocations per assembly, all of it pure bookkeeping.
  if (plus) {
    start <- snp_i + grid$gap
    end <- start + grid$L - 1L
  } else {
    end <- snp_i - grid$gap
    start <- end - grid$L + 1L
  }
  keep <- start >= 1L & end <= wlen
  cand <- tibble::tibble(
    com = NA_character_, len = grid$L[keep], gap = grid$gap[keep],
    start = start[keep], end = end[keep]
  )
  if (nrow(cand)) {
    # substring(), not subseq(): a single DNAString does not vectorise start/end
    seq_plus <- substring(as.character(window), cand$start, cand$end)
    cand$com <- if (plus) {
      as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq_plus)))
    } else {
      seq_plus
    }
  }
  cand <- cand |>
    dplyr::mutate(
      com_gc = gc_frac(.data$com),
      com_tm = primer_tm(.data$com, conditions = conditions),
      run = has_run(.data$com),
      ssr = has_repeat(.data$com)
    ) |>
    dplyr::filter(
      !.data$run, !.data$ssr,
      dplyr::between(.data$com_gc, min_gc, max_gc),
      dplyr::between(.data$com_tm, tm_range[1], tm_range[2])
    ) |>
    dplyr::arrange(abs(.data$com_tm - target_tm))

  if (nrow(cand) == 0L) {
    cli::cli_abort(c(gettext("No common primer passed the filters."),
      i = gettext("Widen `tm_range`/GC or increase `flank` in snp_template().")
    ))
  }
  cand |>
    dplyr::slice_head(n = n) |>
    dplyr::mutate(rank = dplyr::row_number())
}


# =============================================================================
#  4. Allele-specific primers (the KASP-specific part)
# =============================================================================

#' Plus-strand binding window of an AS primer ending at `anchor`
#' @keywords internal
#' @noRd
as_window <- function(window, anchor, length, strand) {
  iv <- primer_interval(anchor, strand, length, role = "AS")
  Biostrings::subseq(window, IRanges::start(iv), IRanges::end(iv))
}

#' DECIPHER target (= strand bound by the primer) for a binding window
#' @keywords internal
#' @noRd
as_target <- function(seg, strand) {
  if (identical(strand, "+")) Biostrings::reverseComplement(seg) else seg
}

#' Plus-strand allele at a primer's 3' terminus
#'
#' [wobble_mismatch()] reads plus-strand alleles, but a primer on the `"-"`
#' strand carries the complement of the template base at its 3' end, so both must
#' be put on one strand before the pair can be asked whether it is a weak
#' G:T / A:C mismatch.  Used for the indel junction's wobble analogue.
#' @param primer oligo(s) 5'->3'.
#' @param strand the primer's strand (`"+"`/`"-"`).
#' @return the 3'-terminal base as a plus-strand allele.
#' @keywords internal
#' @noRd
plus_3p_base <- function(primer, strand) {
  b <- substr(primer, nchar(primer), nchar(primer))
  minus <- !is.na(strand) & strand == "-"
  if (any(minus)) {
    b[minus] <- as.character(Biostrings::complement(Biostrings::DNAStringSet(b[minus])))
  }
  b
}

#' Build one allele-specific primer on a given strand
#'
#' `anchor` is the plus-strand coordinate of the primer's 3'-terminal base (see
#' [as_anchor()]).  The allele's block (1 nt for a SNP, several nt for an MNP) is
#' inserted at the 3' end of the binding window, so the primer ends on the
#' variant.  If `stab_at` is given, the base `stab_at` positions from the 3' end
#' is replaced by `stab_base` (the deliberate stabilising mismatch).
#'
#' `strand = "+"`: primer matches the plus strand, 3' base at the anchor.
#' `strand = "-"`: primer is the reverse complement, 3' base at the anchor.
#' @keywords internal
#' @noRd
make_as_primer <- function(window, anchor, length, allele, strand = "+",
                           stab_at = NULL, stab_base = NULL) {
  seg <- as_window(window, anchor, length, strand) # DNAString, contains the variant
  blk <- block_in_seg(nchar(allele), length, strand) # variant block inside `seg`
  p <- Biostrings::replaceAt(seg, blk, allele)
  if (!identical(strand, "+")) p <- Biostrings::reverseComplement(p)
  if (!is.null(stab_at)) {
    i <- length - stab_at + 1L # 3'-relative -> 5'-index
    p <- Biostrings::replaceAt(p, IRanges::IRanges(i, i), stab_base)
  }
  as.character(p)
}

#' Grid-search allele-specific primers maximising Taq-based discrimination
#'
#' For every strand, length, secondary-mismatch position and mismatch base we
#' build the REF and ALT allele-anchored primers and predict the PCR efficiency of
#' each on both alleles (DECIPHER).  Discrimination = on-target minus off-target
#' efficiency, averaged over the two primers.  `ref`/`alt` are strings, so the
#' same engine handles SNP (`nchar == 1`) and substitution blocks / MNP (longer).
#'
#' When `alt` holds several alleles (a multi-allelic site), one physical ALT
#' primer has to bind them all: the configuration is scored against every ALT and
#' reduced conservatively - on-target efficiency is the *minimum* over the ALTs,
#' off-target the *maximum*, and the REF primer likewise sees its worst ALT
#' template.  `alt_primer` is reported as the worst-ALT member (a concrete
#' sequence, feeding the thermostability/specificity machinery) and `alt_degen`
#' as the IUPAC consensus over every ALT (the physical degenerate oligo, for
#' display and ordering).
#'
#' @param v_start plus-strand start of the variant block (a single base for a SNP).
#' @param ref,alt the alleles as ACGT strings of equal length.  `alt` may hold
#'   more than one allele: a multi-allelic site is then scored against every ALT
#'   and reduced conservatively (see Details).
#' @param strands which strands to consider for the AS primer; `"+"` = forward
#'   primer with the variant as its 3'-terminal block extending downstream, `"-"` =
#'   the mirror-image primer extending upstream.  Both are valid KASP layouts.
#' @examples
#' W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))
#' \dontrun{
#' design_as_primers(W, 40L, ref = "T", alt = "A", lengths = 20:22)
#' }
#' @param window plus-strand template window.
#' @param lengths AS-primer lengths to try.
#' @param stab_ats positions (from the 3' end) for the stabilising mismatch;
#'   `NA` = none.  The decided policy ranks [ARMS_EXTRA_MISMATCH_POSITION] first
#'   and [ARMS_EXTRA_MISMATCH_FALLBACK] second (`stab_policy` records which row
#'   carries which); a position outside `stab_ats` is never invented.
#' @param linked optional two-column table (`offset` 1-based within `window`,
#'   `alt`) of co-segregating substitutions to bake into the ALT haplotype and
#'   ALT primer (makes it mismatch REF twice).
#' @param gc_range,tm_range optional accepted windows; configurations outside
#'   them are dropped **before** the DECIPHER batch (the caller would filter them
#'   afterwards anyway).  `NULL` keeps every configuration.
#' @param temp,P,ions as in [pcr_efficiency()].
#' @param BPPARAM optional `BiocParallel` backend; the batched DECIPHER call and
#'   the 3'-stability map are run across cores (`NULL` = serial).
#' @param conditions reaction conditions for [primer_tm()] (`NULL` = the
#'   historical 50 mM Na / 37 C buffer).  Only the buffer fields of `conditions`
#'   enter the nearest-neighbour Tm.
#' @return a tibble, one row per primer configuration, with `strand`,
#'   `as_length`, `stab_at`, `stab_base`, `ref_primer`, `alt_primer`,
#'   `ref_on`/`ref_off`/`alt_on`/`alt_off`, `discrimination`, `tm`, `gc`,
#'   `dg_3p`, `tmpl_start`, `tmpl_end`.
#' @family KASP design
#' @export
design_as_primers <- function(window, v_start, ref, alt, strands = c("+", "-"),
                              lengths = KASP_AS_LENGTHS, stab_ats = c(NA, 2L, 3L),
                              temp = PCR_EFF$temp, P = PCR_EFF$P, ions = PCR_EFF$ions,
                              linked = NULL, gc_range = NULL, tm_range = NULL, BPPARAM = NULL,
                              conditions = NULL) {
  ref_len <- nchar(ref)
  alts <- unique(as.character(alt))
  alts <- alts[!is.na(alts) & nzchar(alts)]
  if (!length(alts)) {
    return(empty_as_primers())
  }
  if (any(nchar(alts) != ref_len)) {
    cli::cli_abort(gettext("every `alt` allele must be as long as `ref`"))
  }
  # linked variants (substitutions and indels co-segregating with the ALT
  # haplotype) are baked into the ALT template and the ALT primer, so it
  # mismatches REF twice; REF is left untouched.  `ref`/`alt` are the template
  # alleles, so a length change (indel) shifts every coordinate downstream of it
  # - including the variant anchor when the edit is upstream of the block.
  if (is.null(linked) || nrow(linked) == 0L) {
    window_alt <- window
    v_alt <- v_start
  } else {
    # `ref` is optional: a lone `(offset, alt)` row is a 1-base substitution
    ref_w <- if ("ref" %in% names(linked)) nchar(as.character(linked$ref)) else rep(1L, nrow(linked))
    ref_w[is.na(ref_w)] <- 1L
    window_alt <- Biostrings::replaceAt(
      window,
      IRanges::IRanges(as.integer(linked$offset), width = ref_w),
      as.character(linked$alt)
    )
    delta <- nchar(as.character(linked$alt)) - ref_w
    delta[is.na(delta)] <- 0L
    before <- as.integer(linked$offset) < v_start
    v_alt <- v_start + sum(delta[before])
  }
  anchor_of <- function(strand) as_anchor(v_start, ref_len, strand)
  anchor_alt <- function(strand) as_anchor(v_alt, ref_len, strand)
  # candidate mismatch bases at `k` bases from the 3' end: any base other than
  # the one the (perfectly matching) primer already carries there.  The base is
  # READ off the window - the same geometry `make_as_primer()` uses, index
  # `L - k + 1` of the 5'->3' primer - rather than built into a whole primer and
  # then `substr()`ed, which cost one S4 object per (length, position, strand).
  mismatch_bases <- function(L, k, strand) {
    iv <- primer_interval(anchor_of(strand), strand, L, role = "AS")
    i <- L - k + 1L # 3'-relative -> 5' index of the primer
    plus_i <- if (identical(strand, "+")) {
      IRanges::start(iv) + i - 1L
    } else {
      # reverse-complementary primer: its 5' base is the window's last base
      IRanges::end(iv) - i + 1L
    }
    base <- as.character(Biostrings::subseq(window, plus_i, plus_i))
    if (!identical(strand, "+")) base <- as.character(Biostrings::complement(Biostrings::DNAString(base)))
    setdiff(BASES, base)
  }

  # enumerate configurations (strand x length x mismatch position x mismatch base)
  cfgs <- tidyr::crossing(strand = strands, as_length = lengths, stab_at = stab_ats) |>
    dplyr::mutate(stab_base = purrr::pmap(
      list(.data$as_length, .data$stab_at, .data$strand),
      function(L, k, st) if (is.na(k)) NA_character_ else mismatch_bases(L, k, st)
    )) |>
    tidyr::unnest_longer("stab_base")

  # build each configuration's REF primer ONCE, apply the caller's GC/Tm window,
  # and only then build the ALT primer and score the four targets - the caller
  # drops out-of-window rows anyway, so nothing is built or scored twice
  cfgs <- cfgs |>
    dplyr::mutate(
      ref_primer = purrr::pmap_chr(
        list(.data$as_length, .data$strand, .data$stab_at, .data$stab_base),
        function(L, st, k, b) {
          make_as_primer(window, anchor_of(st), L, ref,
            strand = st,
            stab_at = if (is.na(k)) NULL else k, stab_base = if (is.na(b)) NULL else b
          )
        }
      ),
      tm = primer_tm(.data$ref_primer, conditions = conditions),
      gc = gc_frac(.data$ref_primer)
    )
  # the package's soft mismatch index: a G:T / A:C 3' wobble is the class
  # `specificity_policy()` already marks, and the one a deliberate extra mismatch
  # exists to rescue (Ye 2001's weak terminal class).  The pair is the AS primer's
  # own 3' base against the REF allele - `make_as_primer()` can run longer than
  # `as_length`, so the primer's full length is the right index.
  term_base <- substr(cfgs$ref_primer, nchar(cfgs$ref_primer), nchar(cfgs$ref_primer))
  on <- rep(as.character(ref)[1L], nrow(cfgs))
  # the package's soft mismatch index, by the same four pairs `wobble_mismatch()`
  # classifies (C:T / A:G / G:T / A:C) but as a direct element-wise test: the
  # helper sorts a whole vector when given more than one pair
  # a 3' wobble on Ts: a T opposite G or C at the primer's terminal base, which
  # is exactly Kwok's `minimal` class (T:G, T:C, T:T).  Read off the pair rather
  # than through `wobble_mismatch()`, which sorts the WHOLE vector when given
  # more than one pair (checked: it answers FALSE for a batch whose single-pair
  # calls are TRUE).
  wobble <- c("G:T", "C:T", "T:T")
  term_pair <- vapply(seq_along(term_base), function(i) {
    paste(sort(c(term_base[[i]], on[[i]])), collapse = ":")
  }, character(1))
  weak_pairs <- term_pair %in% wobble
  cfgs$weak_terminal <- !is.na(cfgs$stab_at) & weak_pairs
  cfgs$weak_terminal[is.na(cfgs$weak_terminal)] <- FALSE
  # `pmap()` below sees every remaining column, so the policy flag is kept on
  # its own table and joined back by key once the configurations survive
  # kept on its own table and joined back, because purrr::pmap() takes its
  # argument names from the column ORDER of `cfgs`
  policy <- tibble::tibble(
    strand = cfgs[["strand"]], as_length = cfgs[["as_length"]],
    stab_at = cfgs[["stab_at"]], stab_base = cfgs[["stab_base"]],
    weak_terminal = cfgs[["weak_terminal"]]
  )
  cfgs <- dplyr::select(cfgs, -"weak_terminal")
  if (!is.null(gc_range)) {
    cfgs <- dplyr::filter(cfgs, dplyr::between(.data$gc, gc_range[1L], gc_range[2L]))
  }
  if (!is.null(tm_range)) {
    cfgs <- cfgs |>
      dplyr::filter(dplyr::between(.data$tm, tm_range[1L], tm_range[2L]))
  }
  if (nrow(cfgs) == 0L) {
    return(empty_as_primers())
  }
  cfgs$.row_key <- paste(
    cfgs[["strand"]], cfgs[["as_length"]], cfgs[["stab_at"]], cfgs[["stab_base"]],
    sep = "\r"
  )
  policy$.row_key <- paste(
    policy[["strand"]], policy[["as_length"]], policy[["stab_at"]], policy[["stab_base"]],
    sep = "\r"
  )
  # `purrr::pmap()` names its arguments from the column order, so it is given an
  # explicit list; the policy flag joins the scoring by the same key
  cfg_list <- list(
    strand = cfgs[["strand"]], as_length = cfgs[["as_length"]],
    stab_at = cfgs[["stab_at"]], stab_base = cfgs[["stab_base"]],
    ref_primer = cfgs[["ref_primer"]], tm = cfgs[["tm"]], gc = cfgs[["gc"]],
    weak_terminal = policy$weak_terminal[match(cfgs$.row_key, policy$.row_key)]
  )

  # the on/off targets per configuration and per ALT allele, scored in ONE
  # batched DECIPHER call; a multi-allelic site is scored against every ALT and
  # reduced conservatively below
  pairs <- purrr::pmap(cfg_list, function(strand, as_length, stab_at, stab_base, ref_primer, tm, gc, weak_terminal) {
    L <- as_length
    # binding window on the plus strand (REF), and each ALT haplotype carrying
    # the linked substitutions as well
    seg_ref <- as_window(window, anchor_of(strand), L, strand)
    t_ref <- as.character(as_target(seg_ref, strand))
    purrr::map(alts, function(a) {
      alt_primer <- make_as_primer(window_alt, anchor_alt(strand), L, a,
        strand = strand,
        stab_at = if (is.na(stab_at)) NULL else stab_at,
        stab_base = if (is.na(stab_base)) NULL else stab_base
      )
      seg_alt <- Biostrings::replaceAt(
        as_window(window_alt, anchor_alt(strand), L, strand),
        block_in_seg(ref_len, L, strand), a
      )
      t_alt <- as.character(as_target(seg_alt, strand))
      # explicit on/off pairing: the "on" target carries the primer's own allele,
      # the "off" target carries the other one.  Written out (not c()-positional)
      # so a swap cannot silently invert the discrimination metric.
      tibble::tibble(
        which = c("ref_on", "ref_off", "alt_on", "alt_off"),
        primer = c(ref_primer, ref_primer, alt_primer, alt_primer),
        target = c(t_ref, t_alt, t_alt, t_ref),
        strand = strand, as_length = L, stab_at = stab_at, stab_base = stab_base,
        ref_primer = ref_primer, alt_primer = alt_primer, tm = tm, gc = gc,
        weak_terminal = weak_terminal
      )
    }) |>
      purrr::list_rbind()
  }) |>
    purrr::list_rbind()
  pairs$eff <- pcr_efficiency(pairs$primer, pairs$target,
    temp = temp, P = P, ions = ions, cores = marker_ncores(BPPARAM)
  )

  # reduce each configuration over its ALT alleles: one physical primer must bind
  # every ALT, so take the worst on-target (min) and the worst off-target (max);
  # the REF primer likewise sees its worst ALT template.  A single-ALT site is
  # unchanged.  `alt_primer` is the worst-ALT member (a concrete sequence, for
  # the thermostability/specificity machinery); `alt_degen` is the IUPAC
  # consensus over every ALT member - the physical degenerate oligo.
  wide <- pairs |>
    tidyr::pivot_wider(
      id_cols = c(
        "strand", "as_length", "stab_at", "stab_base",
        "ref_primer", "alt_primer", "tm", "gc", "weak_terminal"
      ),
      names_from = "which", values_from = "eff"
    ) |>
    dplyr::group_by(.data$strand, .data$as_length, .data$stab_at, .data$stab_base) |>
    dplyr::summarise(
      ref_primer = dplyr::first(.data$ref_primer),
      # `alt_degen` must merge the ORIGINAL group column, so it is computed before
      # `alt_primer` is (re)assigned below - later expressions see the new column
      alt_degen = merge_iupac(.data$alt_primer),
      alt_primer = .data$alt_primer[which.min(.data$alt_on)],
      tm = dplyr::first(.data$tm),
      gc = dplyr::first(.data$gc),
      weak_terminal = dplyr::first(.data$weak_terminal),
      ref_on = dplyr::first(.data$ref_on),
      ref_off = max(.data$ref_off),
      alt_on = min(.data$alt_on),
      alt_off = max(.data$alt_off),
      .groups = "drop"
    ) |>
    dplyr::mutate(discrimination = ((.data$ref_on - .data$ref_off) + (.data$alt_on - .data$alt_off)) / 2)

  # template-local footprint of the AS primer (genomic coords are added by the
  # caller, which knows the window offset)
  wide <- wide |>
    dplyr::rowwise() |>
    dplyr::mutate(
      tmpl_start = IRanges::start(primer_interval(anchor_of(.data$strand),
        .data$strand, .data$as_length,
        role = "AS"
      )),
      tmpl_end = IRanges::end(primer_interval(anchor_of(.data$strand),
        .data$strand, .data$as_length,
        role = "AS"
      ))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(dg_3p = terminal_dg(.data$ref_primer, BPPARAM = BPPARAM))
  # rank by the DECIDED position policy before discrimination: the measured
  # default (-3) ahead of the rule-of-thumb fallback (-2), then no mismatch, then
  # anything else the caller asked for.  The pool is unchanged - only the order -
  # so no primer appears or disappears by this policy.  A -2 row is the fallback
  # ONLY when the terminal mismatch is weak; otherwise it is just another
  # position the caller asked for.  Positions are compared by ABSOLUTE value:
  # on this R (4.6.1) the scalar literal `3L == -3L` is FALSE while
  # `c(3L) == -3L` is TRUE, so a signed comparison against a constant is not
  # dependable; position runs -1..-4, so abs() is unambiguous.
  rank_default <- !is.na(wide$stab_at) & abs(wide$stab_at) == abs(ARMS_EXTRA_MISMATCH_POSITION)
  rank_fallback <- !is.na(wide$stab_at) &
    abs(wide$stab_at) == abs(ARMS_EXTRA_MISMATCH_FALLBACK) &
    wide$weak_terminal
  wide$.policy_rank <- ifelse(
    rank_default, 1L,
    ifelse(rank_fallback, 2L, ifelse(is.na(wide$stab_at), 3L, 4L))
  )
  wide$stab_policy <- c("default", "fallback", "none", "other")[wide$.policy_rank]
  wide |>
    dplyr::select(dplyr::all_of(AS_PRIMER_COLS), ".policy_rank") |>
    dplyr::arrange(.data$.policy_rank, dplyr::desc(.data$discrimination)) |>
    dplyr::select(-".policy_rank")
}

#' Can any AS primer clear the caller's GC/Tm window at all?
#'
#' A cheap pre-gate for [design_as_primers()].  The grid is `strands x lengths x
#' mismatch positions x mismatch bases`; a configuration whose REF primer falls
#' outside `gc_range`/`tm_range` is dropped *before* the DECIPHER batch, so when
#' EVERY length fails the window the whole grid is built and discarded for
#' nothing.  This builds one ref primer per length - no stabilising mismatch,
#' which can only shift one base - and answers whether any of them can pass.
#'
#' It is a *necessary* condition, not a decision: it says nothing about whether
#' a survivor also clears the homopolymer/SSR or per-site gates.  A `TRUE` means
#' the caller must run the full engine; a `FALSE` means the full engine cannot
#' return a row, which is the case it exists to skip.
#' @param window,v_start,ref,strands,lengths as in [design_as_primers()].
#' @param gc_range,tm_range the accepted windows (`NULL` = no window).
#' @param conditions reaction conditions for [primer_tm()].
#' @return `TRUE` when at least one length could pass, else `FALSE`.
#' @keywords internal
#' @noRd
as_grid_can_pass <- function(window, v_start, ref, strands, lengths,
                             gc_range = NULL, tm_range = NULL, conditions = NULL) {
  if (is.null(gc_range) && is.null(tm_range)) {
    return(TRUE)
  }
  ref_len <- nchar(ref)
  cfg <- tidyr::crossing(strand = strands, as_length = lengths) |>
    dplyr::mutate(
      ref_primer = purrr::pmap_chr(
        list(.data$as_length, .data$strand),
        function(L, st) make_as_primer(window, as_anchor(v_start, ref_len, st), L, ref, strand = st)
      ),
      tm = primer_tm(.data$ref_primer, conditions = conditions),
      gc = gc_frac(.data$ref_primer)
    )
  ok <- rep(TRUE, nrow(cfg))
  if (!is.null(gc_range)) ok <- ok & dplyr::between(cfg$gc, gc_range[1L], gc_range[2L])
  if (!is.null(tm_range)) ok <- ok & dplyr::between(cfg$tm, tm_range[1L], tm_range[2L])
  any(ok, na.rm = TRUE)
}

# The AS-candidate schema, one definition shared by the real table and its
# zero-row twin (so a column added here cannot drift).
AS_PRIMER_COLS <- c(
  "strand", "as_length", "stab_at", "stab_base", "ref_primer", "alt_primer", "alt_degen",
  "ref_on", "ref_off", "alt_on", "alt_off", "discrimination", "tm", "gc",
  "dg_3p", "tmpl_start", "tmpl_end", "weak_terminal", "stab_policy"
)

#' Zero-row `design_as_primers()` result, keeping the schema
#' @keywords internal
#' @noRd
empty_as_primers <- function() {
  tibble::tibble(
    strand = character(), as_length = integer(), stab_at = integer(),
    stab_base = character(), ref_primer = character(), alt_primer = character(),
    alt_degen = character(),
    ref_on = numeric(), ref_off = numeric(), alt_on = numeric(), alt_off = numeric(),
    discrimination = numeric(), tm = numeric(), gc = numeric(), dg_3p = numeric(),
    tmpl_start = integer(), tmpl_end = integer(), weak_terminal = logical(),
    stab_policy = character()
  )[AS_PRIMER_COLS]
}


# =============================================================================
#  4b. Indels: allele templates + junction-scan allele-specific primers
# =============================================================================

#' Plus-strand REF and ALT allele templates for an indel
#'
#' Both templates share the left flank (up to `snp_i - 1`).  The REF template is
#' the genomic window; the ALT template replaces the REF block with the ALT
#' allele.  Coordinates left of the block are identical, and to the right they
#' differ by `delta = nchar(alt) - nchar(ref)`.
#' @param ref_tmpl plus-strand REF template.
#' @param snp_i template coordinate of the first REF base.
#' @param ref,alt alleles as ACGT strings.
#' @return a list with `ref` and `alt` (`DNAString` templates), `delta` (the
#'   length difference) and the `block_ref`/`block_alt` alleles.
#' @examples
#' allele_templates(Biostrings::DNAString("ACGTACGTACGT"), 5L, "G", "TT")
#' @family KASP design
#' @export
allele_templates <- function(ref_tmpl, snp_i, ref, alt) {
  ref_len <- nchar(ref)
  left <- if (snp_i > 1L) Biostrings::subseq(ref_tmpl, 1L, snp_i - 1L) else Biostrings::DNAString("")
  right <- if (snp_i + ref_len <= length(ref_tmpl)) {
    Biostrings::subseq(ref_tmpl, snp_i + ref_len, length(ref_tmpl))
  } else {
    Biostrings::DNAString("")
  }
  list(
    ref = ref_tmpl,
    alt = Biostrings::xscat(left, Biostrings::DNAString(alt), right),
    delta = nchar(alt) - ref_len,
    block_ref = Biostrings::DNAString(ref),
    block_alt = Biostrings::DNAString(alt)
  )
}

#' Allele-specific primer candidates for an insertion/deletion (junction scan)
#'
#' For indels the two alleles differ in *length*, so a fixed-offset substitution
#' no longer applies.  Instead each primer is built from its own allele template
#' and slid along the junction: every candidate is a perfect match to its own
#' allele, and its PCR efficiency against the OTHER allele is scored with
#' DECIPHER (the same Taq-aware model used elsewhere).  The placement that
#' actually discriminates is therefore chosen by the physics, not by a hand rule.
#'
#' @param ref_tmpl plus-strand template (small design window).
#' @param snp_i template coordinate of the first REF base.
#' @param ref,alt alleles as ACGT strings with `nchar(ref) != nchar(alt)`.
#' @param lengths AS-primer lengths to try.
#' @return tidy tibble of candidates: `allele` (`ref`/`alt`), `strand`,
#'   `as_length`, `as_start`/`as_end` (this allele template's coords), `primer`,
#'   `on_eff`, `off_eff`, `discrimination`, `tm`, `gc`.
#' @param strands strands to scan.
#' @param temp,P,ions as in [pcr_efficiency()].
#' @param gc_range,tm_range optional accepted windows; candidates outside them are
#'   dropped **before** the DECIPHER calls.  `NULL` keeps every candidate.
#' @param BPPARAM optional `BiocParallel` backend for the batched DECIPHER calls.
#' @param conditions reaction conditions for [primer_tm()] (`NULL` = the
#'   historical 50 mM Na / 37 C buffer).  Only the buffer fields of `conditions`
#'   enter the nearest-neighbour Tm.
#' @examples
#' \dontrun{
#' tmpl <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))
#' design_as_primers_indel(tmpl, 40L, ref = "G", alt = "GTT")
#' }
#' @family KASP design
#' @export
design_as_primers_indel <- function(ref_tmpl, snp_i, ref, alt,
                                    strands = c("+", "-"), lengths = 18:32,
                                    temp = PCR_EFF$temp, P = PCR_EFF$P, ions = PCR_EFF$ions,
                                    gc_range = NULL, tm_range = NULL, BPPARAM = NULL,
                                    conditions = NULL) {
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)
  if (ref_len == alt_len) {
    cli::cli_abort(gettext(
      "design_as_primers_indel() needs unequal alleles; use design_as_primers() for a substitution."
    ))
  }
  at <- allele_templates(ref_tmpl, snp_i, ref, alt)

  # The OTHER allele's version of a window on THIS template: the window's block
  # overlap is replaced by the other allele, flanks kept.  Windows that miss the
  # block are identical on both alleles (no discrimination).
  substitute_block <- function(tmpl, other_block, block_len, s, e) {
    block_end <- snp_i + block_len - 1L
    if (s > block_end || e < snp_i) {
      return(Biostrings::subseq(tmpl, s, e))
    }
    left <- if (s < snp_i) Biostrings::subseq(tmpl, s, snp_i - 1L) else Biostrings::DNAString("")
    right <- if (e > block_end) Biostrings::subseq(tmpl, block_end + 1L, e) else Biostrings::DNAString("")
    Biostrings::xscat(left, other_block, right)
  }

  build <- function(allele, tmpl, other_block, block_len, strand, L, e) {
    if (identical(strand, "+")) {
      s <- e - L + 1L
      ee <- e
    } else {
      s <- e
      ee <- e + L - 1L
    }
    if (s < 1L || ee > length(tmpl)) {
      return(NULL)
    }
    seg <- Biostrings::subseq(tmpl, s, ee)
    primer <- if (identical(strand, "+")) {
      as.character(seg)
    } else {
      as.character(Biostrings::reverseComplement(seg))
    }
    off_seg <- substitute_block(tmpl, other_block, block_len, s, ee)
    tibble::tibble(
      allele = allele, strand = strand, as_length = L,
      as_start = s, as_end = ee, primer = primer,
      on = as.character(as_target(seg, strand)),
      off = as.character(as_target(off_seg, strand))
    )
  }

  cand <- purrr::pmap(
    tidyr::crossing(allele = c("ref", "alt"), strand = strands, as_length = lengths),
    function(allele, strand, as_length) {
      tmpl <- if (allele == "ref") at$ref else at$alt
      oblk <- if (allele == "ref") at$block_alt else at$block_ref
      blen <- if (allele == "ref") ref_len else alt_len
      # 3'-end placements: from one base before the block to one past its end
      ends <- seq.int(max(1L, snp_i - 1L), min(length(tmpl), snp_i + blen))
      purrr::map(ends, function(e) build(allele, tmpl, oblk, blen, strand, as_length, e)) |>
        purrr::list_rbind()
    }
  ) |>
    purrr::list_rbind()
  if (nrow(cand) == 0L) cli::cli_abort(gettext("No indel AS-primer placement fit inside the design window."))

  # apply the caller's GC/Tm window before the DECIPHER calls (the caller drops
  # those rows afterwards anyway); an empty result is left for the caller
  cand$tm <- primer_tm(cand$primer, conditions = conditions)
  cand$gc <- gc_frac(cand$primer)
  if (!is.null(gc_range)) {
    cand <- dplyr::filter(cand, dplyr::between(.data$gc, gc_range[1L], gc_range[2L]))
  }
  if (!is.null(tm_range)) {
    cand <- dplyr::filter(cand, dplyr::between(.data$tm, tm_range[1L], tm_range[2L]))
  }
  if (nrow(cand) == 0L) {
    return(empty_indel_primers())
  }

  cores <- marker_ncores(BPPARAM)
  cand$on_eff <- pcr_efficiency(cand$primer, cand$on, temp = temp, P = P, ions = ions, cores = cores)
  cand$off_eff <- pcr_efficiency(cand$primer, cand$off, temp = temp, P = P, ions = ions, cores = cores)
  cand |>
    dplyr::mutate(discrimination = .data$on_eff - .data$off_eff) |>
    dplyr::select(
      "allele", "strand", "as_length", "as_start", "as_end", "primer",
      "on_eff", "off_eff", "discrimination", "tm", "gc"
    ) |>
    dplyr::arrange(dplyr::desc(.data$discrimination))
}

#' Zero-row `design_as_primers_indel()` result, keeping the schema
#' @keywords internal
#' @noRd
empty_indel_primers <- function() {
  tibble::tibble(
    allele = character(), strand = character(), as_length = integer(),
    as_start = integer(), as_end = integer(), primer = character(),
    on_eff = numeric(), off_eff = numeric(), discrimination = numeric(),
    tm = numeric(), gc = numeric()
  )
}
