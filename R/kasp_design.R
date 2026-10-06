# =============================================================================
#  kasp_design.R  --  common and allele-specific primer design
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
                                  tm_range = c(57, 64), n = 5L) {
  wlen <- length(window)
  plus <- identical(strand, "+")
  grid <- tidyr::crossing(L = lengths, gap = gaps)

  cand <- purrr::pmap(grid, function(L, gap) {
    iv <- primer_interval(snp_i, strand, L, gap = gap, role = "COM") # template coords
    start <- IRanges::start(iv)
    end <- IRanges::end(iv)
    if (start < 1L || end > wlen) {
      return(NULL)
    }
    seq_plus <- Biostrings::subseq(window, start, end)
    com <- if (plus) {
      as.character(Biostrings::reverseComplement(seq_plus))
    } else {
      as.character(seq_plus)
    }
    tibble::tibble(com = com, len = L, gap = gap, start = start, end = end)
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      com_gc = gc_frac(.data$com),
      com_tm = primer_tm(.data$com),
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
    rlang::abort(c(gettext("No common primer passed the filters."),
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
#' @param v_start plus-strand start of the variant block (a single base for a SNP).
#' @param ref,alt the two alleles as ACGT strings of equal length.
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
#'   `NA` = none.
#' @param linked optional two-column table (`offset` 1-based within `window`,
#'   `alt`) of co-segregating substitutions to bake into the ALT haplotype and
#'   ALT primer (makes it mismatch REF twice).
#' @param gc_range,tm_range optional accepted windows; configurations outside
#'   them are dropped **before** the DECIPHER batch (the caller would filter them
#'   afterwards anyway).  `NULL` keeps every configuration.
#' @param temp,P,ions as in [pcr_efficiency()].
#' @param BPPARAM optional `BiocParallel` backend; the batched DECIPHER call and
#'   the 3'-stability map are run across cores (`NULL` = serial).
#' @return a tibble, one row per primer configuration, with `strand`,
#'   `as_length`, `stab_at`, `stab_base`, `ref_primer`, `alt_primer`,
#'   `ref_on`/`ref_off`/`alt_on`/`alt_off`, `discrimination`, `tm`, `gc`,
#'   `dg_3p`, `tmpl_start`, `tmpl_end`.
#' @family KASP design
#' @export
design_as_primers <- function(window, v_start, ref, alt, strands = c("+", "-"),
                              lengths = 20:27, stab_ats = c(NA, 2L, 3L),
                              temp = PCR_EFF$temp, P = PCR_EFF$P, ions = PCR_EFF$ions,
                              linked = NULL, gc_range = NULL, tm_range = NULL, BPPARAM = NULL) {
  ref_len <- nchar(ref)
  # linked variants (substitutions co-segregating with the ALT haplotype) are
  # baked into the ALT template and the ALT primer, so it mismatches REF twice;
  # REF is left untouched
  window_alt <- if (is.null(linked) || nrow(linked) == 0L) {
    window
  } else {
    Biostrings::replaceAt(
      window,
      IRanges::IRanges(as.integer(linked$offset), as.integer(linked$offset)),
      as.character(linked$alt)
    )
  }
  anchor_of <- function(strand) as_anchor(v_start, ref_len, strand)
  # candidate mismatch bases at `k` bases from the 3' end: any base other than
  # the one the (perfectly matching) primer already carries there.
  mismatch_bases <- function(L, k, strand) {
    base_p <- make_as_primer(window, anchor_of(strand), L, ref, strand = strand)
    cur <- substr(base_p, L - k + 1L, L - k + 1L)
    setdiff(c("A", "C", "G", "T"), cur)
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
      tm = primer_tm(.data$ref_primer),
      gc = gc_frac(.data$ref_primer)
    )
  if (!is.null(gc_range)) {
    cfgs <- dplyr::filter(cfgs, dplyr::between(.data$gc, gc_range[1L], gc_range[2L]))
  }
  if (!is.null(tm_range)) {
    cfgs <- dplyr::filter(cfgs, dplyr::between(.data$tm, tm_range[1L], tm_range[2L]))
  }
  if (nrow(cfgs) == 0L) {
    return(empty_as_primers())
  }
  cfgs$alt_primer <- purrr::pmap_chr(
    list(cfgs$as_length, cfgs$strand, cfgs$stab_at, cfgs$stab_base),
    function(L, st, k, b) {
      make_as_primer(window_alt, anchor_of(st), L, alt,
        strand = st,
        stab_at = if (is.na(k)) NULL else k, stab_base = if (is.na(b)) NULL else b
      )
    }
  )

  # the four on/off targets per configuration, scored in ONE batched DECIPHER call
  pairs <- purrr::pmap(cfgs, function(strand, as_length, stab_at, stab_base, ref_primer, alt_primer, tm, gc) {
    L <- as_length
    # binding window on the plus strand (REF), and the ALT haplotype carrying the
    # linked substitutions as well
    seg_ref <- as_window(window, anchor_of(strand), L, strand)
    seg_alt <- Biostrings::replaceAt(
      as_window(window_alt, anchor_of(strand), L, strand),
      block_in_seg(ref_len, L, strand), alt
    )
    t_ref <- as.character(as_target(seg_ref, strand))
    t_alt <- as.character(as_target(seg_alt, strand))
    # explicit on/off pairing: the "on" target carries the primer's own allele,
    # the "off" target carries the other one.  Written out (not c()-positional)
    # so a swap cannot silently invert the discrimination metric.
    tibble::tibble(
      which = c("ref_on", "ref_off", "alt_on", "alt_off"),
      primer = c(ref_primer, ref_primer, alt_primer, alt_primer),
      target = c(t_ref, t_alt, t_alt, t_ref),
      strand = strand, as_length = L, stab_at = stab_at, stab_base = stab_base,
      ref_primer = ref_primer, alt_primer = alt_primer, tm = tm, gc = gc
    )
  }) |>
    purrr::list_rbind()
  pairs$eff <- pcr_efficiency(pairs$primer, pairs$target,
    temp = temp, P = P, ions = ions, cores = kasp_ncores(BPPARAM)
  )

  wide <- pairs |>
    dplyr::select(
      "strand", "as_length", "stab_at", "stab_base", "ref_primer", "alt_primer",
      "tm", "gc", "which", "eff"
    ) |>
    tidyr::pivot_wider(names_from = "which", values_from = "eff") |>
    dplyr::mutate(discrimination = ((.data$ref_on - .data$ref_off) + (.data$alt_on - .data$alt_off)) / 2)

  # template-local footprint of the AS primer (genomic coords are added by the
  # caller, which knows the window offset)
  wide |>
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
    dplyr::mutate(dg_3p = terminal_dg(.data$ref_primer, BPPARAM = BPPARAM)) |>
    dplyr::select(dplyr::all_of(AS_PRIMER_COLS)) |>
    dplyr::arrange(dplyr::desc(.data$discrimination))
}

# The AS-candidate schema, one definition shared by the real table and its
# zero-row twin (so a column added here cannot drift).
AS_PRIMER_COLS <- c(
  "strand", "as_length", "stab_at", "stab_base", "ref_primer", "alt_primer",
  "ref_on", "ref_off", "alt_on", "alt_off", "discrimination", "tm", "gc",
  "dg_3p", "tmpl_start", "tmpl_end"
)

#' Zero-row `design_as_primers()` result, keeping the schema
#' @keywords internal
#' @noRd
empty_as_primers <- function() {
  tibble::tibble(
    strand = character(), as_length = integer(), stab_at = integer(),
    stab_base = character(), ref_primer = character(), alt_primer = character(),
    ref_on = numeric(), ref_off = numeric(), alt_on = numeric(), alt_off = numeric(),
    discrimination = numeric(), tm = numeric(), gc = numeric(), dg_3p = numeric(),
    tmpl_start = integer(), tmpl_end = integer()
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
#' @examples
#' \dontrun{
#' tmpl <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))
#' design_as_primers_indel(tmpl, 40L, ref = "G", alt = "GTT")
#' }
#' @family KASP design
#' @export
design_as_primers_indel <- function(ref_tmpl, snp_i, ref, alt,
                                    strands = c("+", "-"), lengths = 20:27,
                                    temp = PCR_EFF$temp, P = PCR_EFF$P, ions = PCR_EFF$ions,
                                    gc_range = NULL, tm_range = NULL, BPPARAM = NULL) {
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)
  if (ref_len == alt_len) {
    rlang::abort(gettext(
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

  cand <- purrr::pmap_dfr(
    tidyr::crossing(allele = c("ref", "alt"), strand = strands, as_length = lengths),
    function(allele, strand, as_length) {
      tmpl <- if (allele == "ref") at$ref else at$alt
      oblk <- if (allele == "ref") at$block_alt else at$block_ref
      blen <- if (allele == "ref") ref_len else alt_len
      # 3'-end placements: from one base before the block to one past its end
      ends <- seq.int(max(1L, snp_i - 1L), min(length(tmpl), snp_i + blen))
      purrr::map_dfr(ends, function(e) build(allele, tmpl, oblk, blen, strand, as_length, e))
    }
  )
  if (nrow(cand) == 0L) rlang::abort(gettext("No indel AS-primer placement fit inside the design window."))

  # apply the caller's GC/Tm window before the DECIPHER calls (the caller drops
  # those rows afterwards anyway); an empty result is left for the caller
  cand$tm <- primer_tm(cand$primer)
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

  cores <- kasp_ncores(BPPARAM)
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
