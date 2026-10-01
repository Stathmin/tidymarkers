# =============================================================================
#  kasp_plots.R  --  library-driven visualisations & oligo ordering helpers
# =============================================================================
#  No hand-rolled geometry: every panel is drawn by a dedicated library.
#
#    - ggmsa (Bioconductor)  : DNA multiple-sequence alignments (+ seq logo)
#    - gggenes               : strand-aware arrows (primers on the template,
#                              gene model + variant track)
#    - ggrepel               : non-overlapping labels
#    - OligoArrayAux         : DNA duplex / hairpin free energies (via oligo_dg())
#
#  Hairpins / self-complementarity are shown the DNA way: an alignment of each
#  primer against its own reverse complement (ggmsa). RNA-only tools (R4RNA /
#  ggmsa::geom_helix) are deliberately NOT used - we work with DNA.
#
#  export_primers() writes tidy order sheets for IDT (bulk oligo input) and
#  LGC/KASP (assay order).
# =============================================================================

#' The three assay primers as a named vector
#' @keywords internal
#' @noRd
assay_primers <- function(assay) {
  c(AS_ref = assay$assay$as_ref[1], AS_alt = assay$assay$as_alt[1], COM = assay$assay$com[1])
}

#' Design window (plus-strand genomic coords) as c(start, end)
#' @keywords internal
#' @noRd
design_window <- function(assay, snp) {
  flank <- (length(assay$template$seq) - 1L) / 2L
  c(snp$pos - flank, snp$pos + flank)
}

#' Variants inside a region, as a tidy tibble
#' @keywords internal
#' @noRd
variants_in_window <- function(variants, chrom, win) {
  variants |>
    dplyr::filter(
      .data$chrom == !!chrom,
      dplyr::between(.data$pos, win[1], win[2])
    ) |>
    dplyr::mutate(consequence = forcats::fct_lump_n(.data$consequence, 6L))
}


# =============================================================================
#  One vocabulary: an assay is judged at TWO levels
# =============================================================================
#  * feature level  - one oligo on its own (sequence + structure): Tm, GC, the
#    3' GC-clamp, homopolymer runs, SSRs, hairpins, dimers.  Drawn *on* the
#    primer arrow (a sub-gene span for the clamp, a pin for the mismatch).
#  * alignment level - how the two allele-specific primers and the common primer
#    sit together on the template: amplicon size, COM Tm >= AS Tm, the
#    discrimination the placement buys, the in-silico PCR products.  Drawn on the
#    template and on the REF-vs-ALT ggmsa alignment.
#  Both levels are coloured the same way: ok (green) vs broken (red).

KASP_STATUS <- c(ok = "#1b7837", broken = "#d7301f", neutral = "#4477aa")

#' Feature-level acceptance criteria (everything else is alignment-level)
#' @keywords internal
#' @noRd
check_level <- function(criterion) {
  feature <- c(
    "AS Tm (C)", "COM Tm (C)", "AS GC", "COM GC",
    "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
    "inverted repeat (arms>=4, loop<=6)", "worst hairpin dG37",
    "worst cross-dimer dG37", "worst homodimer dG37",
    "AS 3' window variants (n)", "COM footprint variants (n)"
  )
  ifelse(criterion %in% feature, "feature", "alignment")
}

#' Short labels for the acceptance criteria (for the tumbler column)
#' @keywords internal
#' @noRd
check_short <- function(criterion) {
  dplyr::recode(criterion,
    "amplicon length (bp)" = "amplicon",
    "AS Tm (C)" = "AS Tm",
    "COM Tm (C)" = "COM Tm",
    "COM Tm >= AS Tm" = "COM\u2265AS Tm",
    "longest homopolymer (bp)" = "run",
    "SSR (2-4 nt, >=3 copies)" = "SSR",
    "inverted repeat (arms>=4, loop<=6)" = "inv. repeat",
    "worst hairpin dG37" = "hairpin",
    "worst cross-dimer dG37" = "dimer",
    "worst homodimer dG37" = "homodimer",
    "AS 3' window variants (n)" = "AS 3' var",
    "COM footprint variants (n)" = "COM var",
    "discrimination" = "discrim.",
    "in-silico PCR (REF 1 / ALT 0)" = "PCR"
  )
}

#' The 3' terminal window of a primer span (plus-strand coordinates)
#' @keywords internal
#' @noRd
clamp_window <- function(start, end, forward, n = 5L) {
  if (forward) c(max(start, end - n + 1L), end) else c(start, min(end, start + n - 1L))
}

#' Feature-level track of the assay primers (arrows + 3' clamp + 3' pin)
#'
#' Pure geometry/sequence work, so the plot stays a plot.  `status` is whether the
#' oligo on its own passes its Tm/GC/run/SSR/palindrome tests; `clamp_ok` is
#' whether the 3'-terminal base is a G/C (the "GC clamp" a polymerase needs).
#' @keywords internal
#' @noRd
primer_features <- function(assay, limits = KASP_LIMITS) {
  aa <- assay$assay[1L, ]
  th <- assay$thermo
  plus <- identical(aa$strand %||% "+", "+")
  p <- assay_primers(assay)
  # AS Tm lives in the per-primer thermo table (the assay row only carries COM's)
  as_tm <- mean(th$tm[th$name %in% c("AS_ref", "AS_alt")])
  d <- tibble::tibble(
    role = factor(c("AS-ref", "AS-alt", "COM"), levels = c("COM", "AS-alt", "AS-ref")),
    start = c(aa$as_start, aa$as_alt_start %||% aa$as_start, aa$com_start),
    end = c(aa$as_end, aa$as_alt_end %||% aa$as_end, aa$com_end),
    forward = c(plus, plus, !plus),
    seq = c(p[["AS_ref"]], p[["AS_alt"]], p[["COM"]]),
    tm = c(as_tm, as_tm, aa$com_tm),
    gc = c(aa$as_gc, aa$as_gc, aa$com_gc)
  )
  tm_ok <- c(
    dplyr::between(as_tm, limits$as_tm[1], limits$as_tm[2]),
    dplyr::between(as_tm, limits$as_tm[1], limits$as_tm[2]),
    dplyr::between(aa$com_tm, limits$com_tm[1], limits$com_tm[2])
  )
  clean <- !vapply(d$seq, has_run, logical(1)) &
    !vapply(d$seq, has_repeat, logical(1)) &
    !vapply(d$seq, has_palindrome, logical(1))
  d$status <- ifelse(tm_ok & dplyr::between(d$gc, limits$gc[1], limits$gc[2]) & clean, "ok", "broken")
  win <- Map(clamp_window, d$start, d$end, d$forward)
  d$cstart <- vapply(win, `[`, numeric(1), 1L)
  d$cend <- vapply(win, `[`, numeric(1), 2L)
  # the primer string is always 5'->3'; the 3' base is its last letter
  last <- nchar(d$seq)
  d$three_x <- ifelse(d$forward, d$end, d$start)
  d$clamp <- ifelse(substr(d$seq, last, last) %in% c("G", "C"), "ok", "broken")
  d$label <- glue::glue("{d$role} \u00b7 {round(d$tm)}C {round(100 * d$gc)}%")
  # the injected destabilising mismatch is not in the assay row; recover it from
  # the AS candidates (the chosen primer's own row)
  d$stab <- NA_integer_
  ac <- assay$as_candidates
  if (!is.null(ac) && all(c("stab_at", "ref_primer") %in% names(ac))) {
    hit <- ac$stab_at[!is.na(ac$stab_at) & ac$ref_primer == p[["AS_ref"]]]
    if (length(hit)) d$stab[1:2] <- hit[[1L]]
  }
  d
}

#' SSR spans (start, end) within one sequence
#' @keywords internal
#' @noRd
repeat_spans <- function(s, units = 2:4, min_copies = 3L) {
  spans <- matrix(integer(0), ncol = 2L)
  for (u in units) {
    pat <- glue::glue("([ACGT]{{{u}}})\\1{{{min_copies - 1L},}}")
    m <- stringr::str_locate_all(s, pat)[[1L]]
    if (nrow(m)) spans <- rbind(spans, m)
  }
  spans
}

#' Per-base letters of the assay primers, styled for KASP
#'
#' The primer sequence as one row of letters at genomic coordinates: neutral
#' black by default, **bold** over the 3' GC-clamp window, and bright red on the
#' bases of an SSR or a hairpin stem.  Coordinates follow the strand, so a
#' minus-strand primer reads right-to-left.
#' @keywords internal
#' @noRd
primer_letters <- function(assay, clamp_n = 5L) {
  d <- primer_features(assay)
  rows <- lapply(seq_len(nrow(d)), function(i) {
    s <- d$seq[i]
    L <- nchar(s)
    xs <- if (d$forward[i]) seq.int(d$start[i], length.out = L) else seq.int(d$end[i], by = -1L, length.out = L)
    face <- rep("plain", L)
    # sequence is grey by default; only the KASP-relevant bases stand out (GC-clamp
    # in bold black, SSR/hairpin bases in red)
    colour <- rep("grey55", L)
    # `s` is always 5'->3', so the 3' GC-clamp is the last `clamp_n` letters,
    # whatever the genomic strand
    clamp <- (L - clamp_n + 1L):L
    clamp <- clamp[clamp >= 1L & clamp <= L]
    face[clamp] <- "bold"
    colour[clamp] <- "black"
    sp <- repeat_spans(s)
    pal <- find_palindromes(s)[[1L]]
    if (nrow(pal)) sp <- rbind(sp, as.matrix(pal[, c("start", "end")]))
    if (nrow(sp)) {
      idx <- unique(unlist(lapply(seq_len(nrow(sp)), function(j) sp[j, 1L]:sp[j, 2L])))
      colour[idx[idx >= 1L & idx <= L]] <- "#d7191c"
    }
    tibble::tibble(
      role = as.character(d$role[i]), x = xs, letter = strsplit(s, "")[[1L]],
      face = face, colour = colour
    )
  })
  dplyr::bind_rows(rows)
}

#' Where every primer matches the genome (positions), for the ideogram
#' @keywords internal
#' @noRd
primer_hits <- function(primers, genome) {
  purrr::list_rbind(lapply(genome_chroms(genome), function(ch) {
    subject <- Biostrings::DNAStringSet(genome_seq(genome, ch))
    purrr::imap_dfr(primers, function(p, nm) {
      # a primer binds either strand: search the plus strand for the primer and
      # for its reverse complement, and report the plus-strand match starts
      pats <- unique(c(p, as.character(Biostrings::reverseComplement(Biostrings::DNAString(p)))))
      pos <- unlist(lapply(pats, function(pt) {
        IRanges::start(Biostrings::vmatchPattern(Biostrings::DNAString(pt), subject, max.mismatch = 0L)[[1L]])
      }))
      if (!length(pos)) {
        return(NULL)
      }
      tibble::tibble(primer = nm, chrom = ch, pos = sort(unique(pos)))
    })
  }))
}

#' Genome ideogram: exact primer hits + in-silico PCR products
#'
#' A thin, schematic view of the whole genome: every chromosome is drawn as a
#' grey bar on its own row (scaled to fill the width, karyotype-style) and two
#' signals are overlaid - where each assay primer matches exactly
#' (`primer_hits()`) and where the designed allele-specific pairs actually
#' amplify ([in_silico_pcr()], the same pairings as [assay_amplicons()]).  The
#' target locus is marked.  Returns `NULL` when no genome is supplied, so
#' [plot_primers()] simply drops the panel.
#' @keywords internal
#' @noRd
assay_ideogram <- function(assay, snp, genome) {
  if (is.null(genome)) {
    return(NULL)
  }
  chroms <- genome_chroms(genome)
  bars <- tibble::tibble(
    chrom = chroms,
    length = vapply(chroms, function(ch) genome_length(genome, ch), numeric(1))
  )
  roles <- c(AS_ref = "AS-ref", AS_alt = "AS-alt", COM = "COM")
  sig_cols <- c("AS-ref" = "#4477aa", "AS-alt" = "#66c2a5", "COM" = "#8da0cb")
  # adaptive axis unit: Mb for real genomes, kb/bp for small ones
  span <- max(bars$length)
  div <- if (span >= 1e6) 1e6 else if (span >= 1e3) 1e3 else 1
  suf <- if (span >= 1e6) " Mb" else if (span >= 1e3) " kb" else " bp"

  hits <- primer_hits(assay_primers(assay), genome)
  if (nrow(hits)) hits$signal <- roles[hits$primer]
  prods <- purrr::imap_dfr(as_pcr_pairs(assay), function(fr, nm) {
    in_silico_pcr(fr[["F"]], fr[["R"]], genome) |> dplyr::mutate(signal = roles[[nm]])
  })

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = bars,
      ggplot2::aes(x = 1, xend = .data$length, y = 1, yend = 1),
      linewidth = 2.8, colour = "grey72", lineend = "round"
    ) +
    # the target locus, as a thin cursor on its chromosome
    ggplot2::geom_segment(
      data = snp,
      ggplot2::aes(x = .data$pos, xend = .data$pos, y = 0.74, yend = 1.26),
      colour = "firebrick", linewidth = 0.4, alpha = 0.8
    )
  # in-silico PCR products (the amplicon), at their midpoint
  if (nrow(prods)) {
    prods$mid <- (prods$start + prods$end) / 2
    p <- p + ggplot2::geom_point(
      data = prods,
      ggplot2::aes(x = .data$mid, y = 1, colour = .data$signal),
      shape = 18, size = 2.1
    )
  }
  # exact primer hits, as small upright ticks
  if (nrow(hits)) {
    p <- p + ggplot2::geom_point(
      data = hits,
      ggplot2::aes(x = .data$pos, y = 1.12, colour = .data$signal),
      shape = 124, size = 2.6
    )
  }
  p +
    ggplot2::scale_colour_manual(values = sig_cols, name = NULL, drop = TRUE) +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(x / div, 1), suf)) +
    ggplot2::scale_y_continuous(limits = c(0.6, 1.4)) +
    # chromosome names on the left, sized for accession-style labels
    ggplot2::facet_grid(chrom ~ ., switch = "y") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = "black"),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 6, colour = "black"),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(6, "pt"),
      legend.text = ggplot2::element_text(size = 6)
    )
}


# =============================================================================
#  Alignments (ggmsa)
# =============================================================================

#' AS primers aligned against the template; SNP + mismatches visible
#'
#' A DNA MSA rendered by ggmsa.  The two template rows are the primer-binding
#' sequences carrying the REF/ALT allele (without the deliberate stabilising
#' mismatch); the two AS rows are the real primers, so a column that differs
#' between a template row and its AS row is the injected secondary mismatch.
#' The AS 3'-terminal column - exactly where REF and ALT disagree - is
#' highlighted (`position_highlight`).  Works for either assay orientation.
#' (ggmsa's add-on geoms are deliberately not used: `geom_seqlogo()` only works
#' when ggmsa is attached, and `geom_GC()` fails to build.)
#' @param assay assembled assay (from [assemble_kasp_assay()]).
#' @return a `ggmsa` alignment plot (a `ggplot` object).
#' @examples
#' \dontrun{
#' plot_alignment(assay)
#' }
#' @family KASP plots
#' @export
plot_alignment <- function(assay) {
  tmpl <- assay$template$seq
  snp_i <- assay$template$snp_i
  p <- assay_primers(assay)
  L <- nchar(p[["AS_ref"]])
  strand <- assay$strand %||% "+"
  anchor <- as_anchor(snp_i, nchar(assay$snp$ref), strand)
  # template rows in the SAME 5'->3' orientation as the primers, one per allele
  t_ref <- make_as_primer(tmpl, anchor, L, assay$snp$ref, strand = strand)
  t_alt <- make_as_primer(tmpl, anchor, L, assay$snp$alt, strand = strand)
  msa <- Biostrings::DNAStringSet(c(
    template_ref = t_ref, template_alt = t_alt,
    AS_ref = p[["AS_ref"]], AS_alt = p[["AS_alt"]]
  ))
  # The AS primer is anchored on the variant, so its 3'-terminal column (column
  # L in this orientation) is exactly where template_ref and template_alt differ:
  # highlight that whole column - it is the base KASP discrimination rests on.
  suppressWarnings(
    ggmsa::ggmsa(msa, 1, L,
      color = "Chemistry_NT", seq_name = TRUE, char_width = 0.6,
      position_highlight = L
    )
  )
}

#' Self-complementarity of each primer (DNA): primer vs its reverse complement
#'
#' Columns that match are self-complementary (hairpin-prone stems).
#' @param assay assembled assay.
#' @param title panel title (or `NULL`).
#' @return a `patchwork` object stacking one primer-vs-reverse-complement
#'   alignment per assay primer.
#' @examples
#' \dontrun{
#' plot_self_complementarity(assay)
#' }
#' @family KASP plots
#' @export
plot_self_complementarity <- function(
  assay,
  title = "Self-complementarity (DNA): matching columns = hairpin-prone stems"
) {
  plots <- purrr::imap(assay_primers(assay), function(s, nm) {
    rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
    suppressWarnings(
      ggmsa::ggmsa(Biostrings::DNAStringSet(c(primer = s, revcomp = rc)),
        1, nchar(s),
        color = "Chemistry_NT", seq_name = TRUE, char_width = 0.5
      )
    ) +
      ggplot2::theme(plot.margin = ggplot2::margin(1, 4, 1, 4))
  })
  p <- patchwork::wrap_plots(plots, ncol = 1)
  if (!is.null(title)) p <- p + patchwork::plot_annotation(title = title)
  p
}


# =============================================================================
#  Gene / primer maps (gggenes)
# =============================================================================

#' Local assay map: styled primers, Tm, and the acceptance checks
#'
#' The local view of the assay on the design window (genomic coordinates).  Each
#' primer is drawn as its sequence of letters: neutral black by default, **bold**
#' over the 3' GC-clamp window, and **bright red** on the bases of an SSR or a
#' hairpin stem; a ring marks an injected mismatch / the 3' discriminating base.
#' The Tm sits just past the 3' end, with its deviation from the target as a
#' green/red tumbler.  Below: the amplicon's own sequence, the local variant
#' track and a feature-annotation row (same x axis), then - when `genome` is
#' given - a thin genome-wide ideogram (`assay_ideogram()`) of exact primer hits
#' and in-silico PCR products.  Right: a thin column of per-criterion tumblers.
#' @param assay assembled assay.
#' @param variants variant table.
#' @param snp target variant (one row).
#' @param genome optional `kasp_genome` or `DNAStringSet`; when given, the
#'   genome-wide ideogram panel is added.
#' @param features optional annotation to draw on the feature row: a data frame
#'   with `start`/`end` (and optional `strand`, `name`), or a gene record from
#'   [ensembl_gene()].  A KASP marker should sit on an annotated feature, so the
#'   row says `"no feature annotation"` when nothing (or no overlap) is given.
#' @return a `patchwork` object: the local map (primers over the variant track,
#'   sharing an x axis) beside the tumbler list from [plot_checks()].
#' @examples
#' \dontrun{
#' plot_primers(assay, variants, snp)
#' }
#' @family KASP plots
#' @export
plot_primers <- function(assay, variants, snp, genome = NULL, features = NULL) {
  win <- design_window(assay, snp)
  aa <- assay$assay[1, ]
  strand <- aa$strand %||% "+"
  plus <- identical(strand, "+")
  d <- primer_features(assay)
  letters <- primer_letters(assay)

  # alignment-level arrow: the amplicon the placement buys, on its own row
  amp <- tibble::tibble(
    role = factor("amplicon", levels = c(levels(d$role), "amplicon")),
    start = min(d$start), end = max(d$end), forward = plus,
    status = ifelse(dplyr::between(aa$amplicon, KASP_LIMITS$amplicon[1], KASP_LIMITS$amplicon[2]), "ok", "broken")
  )
  # the amplicon's own sequence (plus strand) as letters on its row
  ti <- function(p) p - snp$pos + assay$template$snp_i
  amp_seq <- as.character(Biostrings::subseq(assay$template$seq, ti(amp$start), ti(amp$end)))
  letters <- dplyr::bind_rows(letters, tibble::tibble(
    role = "amplicon", x = seq.int(amp$start, amp$end),
    letter = strsplit(amp_seq, "")[[1L]], face = "plain", colour = "grey55"
  ))

  # the Tm as a small number just past the 3' end, and its deviation from the
  # target as a green/red tumbler next to it
  lab <- d
  lab$target <- ifelse(lab$role == "COM", COM_TM_TARGET, AS_TM_TARGET)
  lab$dtm <- lab$tm - lab$target
  lab$meas_x <- ifelse(lab$forward, lab$end + 1, lab$start - 1)
  lab$hj <- ifelse(lab$forward, 0, 1)
  lab$dtm_x <- ifelse(lab$forward, lab$meas_x + 4, lab$meas_x - 4)
  lab$dtm_status <- ifelse(abs(lab$dtm) <= 3, "ok", "broken")

  # pins: the 3' discriminating base of every oligo, and the injected
  # destabilising mismatch of the AS primers (if the placement carries one)
  pins <- tibble::tibble(role = as.character(d$role), x = d$three_x)
  asx <- d[!is.na(d$stab), ]
  if (nrow(asx)) {
    pins <- dplyr::bind_rows(pins, tibble::tibble(
      role = as.character(asx$role),
      x = ifelse(asx$forward, asx$end - (asx$stab - 1L), asx$start + (asx$stab - 1L))
    ))
  }

  span <- range(c(d$start, d$end, lab$meas_x, lab$dtm_x, amp$start, amp$end))

  # the annotated features on the same axis: a KASP marker should sit on a gene,
  # not in intergenic space, so if no (overlapping) annotation is given we say so
  feat <- NULL
  if (!is.null(features)) {
    features <- as.data.frame(features)
    if (nrow(features) && all(c("start", "end") %in% names(features))) {
      nm <- if ("name" %in% names(features)) {
        features$name
      } else if ("symbol" %in% names(features)) {
        features$symbol
      } else {
        NA_character_
      }
      feat <- tibble::tibble(
        start = as.integer(features$start),
        end = as.integer(features$end),
        forward = if ("strand" %in% names(features)) features$strand != "-" else TRUE,
        name = nm
      )
      feat <- feat[feat$end >= span[1] & feat$start <= span[2], ]
      if (nrow(feat) == 0L) feat <- NULL
    }
  }

  p_map <- ggplot2::ggplot() +
    gggenes::geom_gene_arrow(
      data = amp,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$role, forward = .data$forward),
      fill = "grey94", colour = "grey55"
    ) +
    gggenes::geom_gene_arrow(
      data = d,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$role, forward = .data$forward),
      fill = "grey98", colour = "grey75"
    ) +
    ggplot2::geom_text(
      data = letters,
      ggplot2::aes(x = .data$x, y = .data$role, label = .data$letter, fontface = .data$face, colour = .data$colour),
      family = "mono", size = 2.6
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_point(
      data = pins, ggplot2::aes(x = .data$x, y = .data$role),
      shape = 1, size = 3.4, colour = "grey35"
    ) +
    ggplot2::geom_text(
      data = lab,
      ggplot2::aes(x = .data$meas_x, y = .data$role, label = glue::glue("{round(.data$tm)}\u00b0C"), hjust = .data$hj),
      size = 2.6, colour = "black"
    ) +
    ggplot2::geom_point(
      data = lab, ggplot2::aes(x = .data$dtm_x, y = .data$role, fill = .data$dtm_status),
      shape = 21, size = 2.8, colour = "grey40"
    ) +
    ggplot2::geom_vline(xintercept = snp$pos, linetype = "dashed", colour = "grey45") +
    # feature-annotation row (baseline + boxes, or a note when there is none)
    ggplot2::geom_segment(
      data = tibble::tibble(x = span[1], xend = span[2]),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = "feature", yend = "feature"),
      colour = "grey72", linewidth = 0.3
    )
  if (is.null(feat)) {
    p_map <- p_map + ggplot2::annotate("text",
      x = mean(span), y = "feature", label = "no feature annotation",
      size = 2.1, colour = "grey55", fontface = "italic"
    )
  } else {
    p_map <- p_map +
      gggenes::geom_gene_arrow(
        data = feat,
        ggplot2::aes(xmin = .data$start, xmax = .data$end, y = "feature", forward = .data$forward),
        fill = "grey72", colour = "grey40"
      )
    if (any(!is.na(feat$name))) {
      lab_feat <- feat[!is.na(feat$name), ]
      # keep the label inside the visible window even for a feature that spans it
      lab_feat$x <- (pmax(lab_feat$start, span[1]) + pmin(lab_feat$end, span[2])) / 2
      p_map <- p_map + ggplot2::geom_text(
        data = lab_feat,
        ggplot2::aes(x = .data$x, y = "feature", label = .data$name),
        size = 2.1, colour = "black", fontface = "bold", vjust = 0.5
      )
    }
  }
  p_map <- p_map +
    ggplot2::scale_fill_manual(values = KASP_STATUS, guide = "none") +
    ggplot2::scale_y_discrete(limits = c("feature", "amplicon", "COM", "AS-alt", "AS-ref")) +
    ggplot2::coord_cartesian(xlim = span) +
    ggplot2::labs(x = NULL, y = NULL) +
    gggenes::theme_genes() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black")
    )

  vw <- variants_in_window(variants, snp$chrom, win)
  reg3p <- range(c(aa$win3p_start, aa$win3p_end))
  as_span <- c(aa$as_start, aa$as_end)
  com_span <- c(aa$com_start, aa$com_end)
  in_foot <- dplyr::between(vw$pos, as_span[1], as_span[2]) |
    dplyr::between(vw$pos, com_span[1], com_span[2])
  foot <- vw[in_foot, ]
  # the AS 3'-terminal window (the ref/alt side): green when no other known
  # variant sits there - a clean 3' end is what buys allele specificity - red
  # when one does
  win_status <- if (isTRUE(aa$as_n_3p == 0L)) "ok" else "broken"

  # the local variants: a thin strip on the same x axis
  p_var <- ggplot2::ggplot(vw, ggplot2::aes(.data$pos, 1)) +
    ggplot2::annotate("rect",
      xmin = reg3p[1] - 0.5, xmax = reg3p[2] + 0.5,
      ymin = 0.6, ymax = 1.4, fill = KASP_STATUS[[win_status]], alpha = 0.22
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$consequence), size = 1.3) +
    ggplot2::geom_point(data = foot, ggplot2::aes(.data$pos, 1), shape = 1, colour = "black", size = 2.1) +
    ggplot2::geom_point(
      data = snp, ggplot2::aes(.data$pos, 1),
      shape = 21, fill = "firebrick", colour = "white", size = 2.4
    ) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL) +
    ggplot2::scale_y_continuous(breaks = 1, labels = "SNP") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(xlim = span) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black")
    )

  # genome-wide ideogram: where the primers hit and where in-silico PCR amplifies
  p_ideo <- assay_ideogram(assay, snp, genome)

  # primer rows are packed tight; the feature row sits at the bottom
  panels <- if (is.null(p_ideo)) {
    patchwork::wrap_plots(p_map, p_var, ncol = 1, heights = c(2.3, 0.4))
  } else {
    patchwork::wrap_plots(p_map, p_var, p_ideo, ncol = 1, heights = c(2.3, 0.4, 3))
  }

  patchwork::wrap_plots(panels, plot_checks(assay), widths = c(1, 0.12)) +
    patchwork::plot_annotation(
      title = "KASP assay: local map + acceptance checks",
      subtitle = glue::glue(
        "target {snp$chrom}:{snp$pos} ({snp$ref}>{snp$alt}) on strand {strand}"
      )
    )
}

#' Format one acceptance criterion's observed value for the tumbler column
#'
#' Every value is monospaced and at most 4 characters wide, so the column lines
#' up: GC and discrimination as whole percents (`dd%`), Tm and the four dG terms
#' as `dd.d` / `-d.d`, and integer counts plain (compacted to `dK` / `dM` for
#' extreme values).
#' @keywords internal
#' @noRd
format_check_value <- function(name, value) {
  out <- dplyr::case_when(
    name %in% c("AS GC", "COM GC", "discrim.") ~ sprintf("%.0f%%", 100 * value),
    name %in% c("AS Tm", "COM Tm", "dimer", "homodimer", "hairpin", "COM\u2265AS Tm") ~ sprintf("%.1f", value),
    abs(value) >= 1e6 ~ paste0(round(value / 1e6), "M"),
    abs(value) >= 1e3 ~ paste0(round(value / 1e3), "K"),
    TRUE ~ as.character(round(value))
  )
  out[is.na(value)] <- "-"
  out
}

#' Acceptance checks as a thin list of tumblers (feature vs alignment)
#'
#' Runs [kasp_design_checks()] and shows one named tumbler per criterion - a
#' small green (ok) / red (broken) marker with its observed value and limit - in
#' a single narrow column, faceted by the level judged: the *feature* level is
#' one oligo on its own, the *alignment* level is how the primers sit together on
#' the template.  Sits as a thin strip beside [plot_primers()].
#' @param assay assembled assay.
#' @param limits thresholds (see `KASP_LIMITS`).
#' @return a `ggplot` tumbler list.
#' @examples
#' \dontrun{
#' plot_checks(assay)
#' }
#' @family KASP plots
#' @export
plot_checks <- function(assay, limits = KASP_LIMITS) {
  chk <- kasp_design_checks(assay, limits) |>
    dplyr::mutate(
      level = factor(check_level(.data$criterion), levels = c("feature", "alignment")),
      name = check_short(.data$criterion),
      status = ifelse(.data$pass, "ok", "broken"),
      value = format_check_value(.data$name, .data$observed)
    )
  ggplot2::ggplot(chk, ggplot2::aes(x = 0.5, y = forcats::fct_rev(.data$name))) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$status), width = 0.3, height = 0.28) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$value),
      x = 0.92, hjust = 1, size = 2.3, colour = "black", family = "mono"
    ) +
    ggplot2::facet_grid(level ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_manual(values = KASP_STATUS, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0.32, 0.97), breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6.5, colour = "black"),
      strip.text.y = ggplot2::element_blank()
    )
}

#' Gene model with the 1001-Genomes variant track + inset of the target region
#' @param gene gene record.
#' @param variants variant table.
#' @param snp target variant (one row).
#' @param assay assembled assay.
#' @return a `ggplot` object: the gene model and variant track with an inset of
#'   the target region.
#' @examples
#' \dontrun{
#' plot_variant_landscape(gene, variants, snp, assay)
#' }
#' @family KASP plots
#' @export
plot_variant_landscape <- function(gene, variants, snp, assay) {
  gene_df <- tibble::tibble(
    molecule = gene$symbol, start = gene$start, end = gene$end,
    strand = gene$strand
  )
  v <- variants |>
    dplyr::mutate(consequence = forcats::fct_lump_n(.data$consequence, 6L))

  main <- ggplot2::ggplot() +
    gggenes::geom_gene_arrow(
      data = gene_df,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = 1, forward = .data$strand > 0),
      fill = "grey80"
    ) +
    ggplot2::geom_point(
      data = v, ggplot2::aes(.data$pos, 2, colour = .data$consequence),
      size = 1.3, shape = 108
    ) +
    ggplot2::geom_point(
      data = snp, ggplot2::aes(.data$pos, 2),
      shape = 21, fill = "firebrick", colour = "white", size = 3
    ) +
    ggrepel::geom_text_repel(
      data = snp, ggplot2::aes(.data$pos, 2, label = .data$id),
      size = 2.6, nudge_y = 0.6, min.segment.length = 0
    ) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "consequence") +
    ggplot2::scale_y_continuous(
      limits = c(0.5, 3), breaks = c(1, 2),
      labels = c(gene$symbol, "variants")
    ) +
    ggplot2::labs(
      x = glue::glue("{gene$chrom} position (bp)"), y = NULL,
      title = glue::glue("Variation across {gene$symbol} (1001 Genomes)")
    ) +
    ggplot2::theme_minimal()

  win <- design_window(assay, snp)
  vw <- variants_in_window(variants, snp$chrom, win)
  inset <- ggplot2::ggplot(vw, ggplot2::aes(.data$pos, 1)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$consequence), size = 1.1) +
    ggplot2::geom_point(
      data = snp, ggplot2::aes(.data$pos, 1),
      shape = 21, fill = "firebrick", colour = "white", size = 2
    ) +
    ggplot2::scale_colour_brewer(palette = "Dark2", guide = "none") +
    ggplot2::scale_y_continuous(breaks = 1, labels = "variants") +
    ggplot2::labs(
      x = glue::glue("design window ({win[1]}-{win[2]} bp)"), y = NULL,
      title = "target region"
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = "grey70"))

  main + patchwork::inset_element(inset, left = 0.53, bottom = 0.42, right = 1, top = 1)
}


# =============================================================================
#  Cross-complementarity network (ggraph / tidygraph)
# =============================================================================

#' Dimer QC: pairwise duplex dG heatmap + hairpin dG (OligoArrayAux, DNA)
#'
#' Readable replacement for a 3-node network: a symmetric matrix of duplex
#' free energies (kcal/mol, 37 C; more negative = stronger dimer) plus a bar
#' panel of intramolecular folding (hairpin) free energies.
#' @param assay assembled assay.
#' @return a `patchwork` object: the duplex dG heatmap and hairpin dG bars over
#'   the self-complementarity panels.
#' @examples
#' \dontrun{
#' plot_dimer_qc(assay)
#' }
#' @family KASP plots
#' @export
plot_dimer_qc <- function(assay) {
  d <- assay$complementarity # already a full symmetric matrix

  heat <- ggplot2::ggplot(d, ggplot2::aes(.data$a, .data$b, fill = .data$dg)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = round(.data$dg, 1)), size = 4) +
    ggplot2::scale_fill_gradient(
      low = "firebrick", high = "white",
      limits = c(min(d$dg, na.rm = TRUE), 0),
      name = "dG37\n(kcal/mol)"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, title = "Duplex dG37 (kcal/mol)") +
    ggplot2::theme_minimal()

  hair <- ggplot2::ggplot(
    assay$thermo,
    ggplot2::aes(.data$hairpin_dg, stats::reorder(.data$name, .data$hairpin_dg))
  ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_vline(xintercept = HAIRPIN_DG_WARN, linetype = "dashed", colour = "firebrick") +
    ggplot2::labs(
      x = "hairpin dG37 (kcal/mol)", y = NULL,
      title = "Hairpin dG37 (kcal/mol)"
    ) +
    ggplot2::theme_minimal()

  sc <- plot_self_complementarity(assay, title = NULL)
  (heat | hair) / sc +
    patchwork::plot_layout(heights = c(1, 1.1)) +
    patchwork::plot_annotation(title = "Primer QC: duplex/hairpin dG37 + self-complementarity")
}


# =============================================================================
#  Oligo ordering (IDT / LGC)
# =============================================================================

#' Write the assay results as a few consolidated files
#'
#' Instead of many tiny CSVs, the outputs are grouped by what they describe and
#' the per-consumer formats are column subsets of the same data:
#'
#' * `kasp_design.csv` — candidate assays (one row per common-primer placement);
#' * `kasp_primer_qc.csv` — one row per primer: Tm/GC/hairpin/3'-stability +
#'   3'-seed specificity + in-silico PCR products;
#' * `kasp_dimer_dg.csv` — duplex dG37 per primer pair (long);
#' * `{base}_tidy.csv` / `_idt.csv` / `_lgc.csv` — the order table and its
#'   IDT / LGC column subsets (see [export_primers()]).
#'
#' @return character vector of the written paths (invisibly).
#' @param assay assembled assay.
#' @param gene gene record (used for the file base name).
#' @param dir output directory.
#' @examples
#' \dontrun{
#' write_kasp_results(assay, gene, dir = tempdir())
#' }
#' @family KASP output
#' @export
write_kasp_results <- function(assay, gene, dir = tempdir()) {
  fs::dir_create(dir)
  base <- glue::glue("kasp_{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}")

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
    fs::path(dir, "kasp_design.csv"),
    fs::path(dir, "kasp_primer_qc.csv"),
    fs::path(dir, "kasp_dimer_dg.csv")
  )
  readr::write_csv(design, paths[1L])
  readr::write_csv(qc, paths[2L])
  readr::write_csv(assay$complementarity, paths[3L])
  if (!is.null(assay$blast)) {
    p <- fs::path(dir, "kasp_blast_hits.csv")
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

#' Tidy table of the assay oligos (with KASP tails)
#' @param assay assembled assay.
#' @param gene gene record.
#' @return a tibble (`assay`, `role`, `oligo_id`, `sequence`, `seq_no_tail`,
#'   `tail`, `tm`, `gc`, `length`) with the three assay oligos and their KASP
#'   tails.
#' @examples
#' \dontrun{
#' kasp_primer_table(assay, gene)
#' }
#' @family KASP output
#' @export
kasp_primer_table <- function(assay, gene) {
  nm <- glue::glue("{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}_{assay$snp$ref}to{assay$snp$alt}")
  roles <- tibble::tibble(
    name        = c("AS_ref", "AS_alt", "COM"),
    role        = c("AS-ref (FAM)", "AS-alt (HEX/VIC)", "common (reverse)"),
    oligo_id    = c("AS1_FAM", "AS2_HEX", "COM"),
    tail        = c(KASP_TAIL_FAM, KASP_TAIL_HEX, ""),
    seq_no_tail = c(assay$assay$as_ref[1], assay$assay$as_alt[1], assay$assay$com[1])
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
#' @return the tidy primer table (see [kasp_primer_table()]), invisibly; called
#'   for its side effect of writing the IDT/LGC/tidy order sheets.
#' @examples
#' \dontrun{
#' export_primers(assay, gene, dir = tempdir())
#' }
#' @family KASP output
#' @export
export_primers <- function(assay, gene, dir = tempdir()) {
  fs::dir_create(dir)
  tbl <- kasp_primer_table(assay, gene)
  base <- glue::glue("kasp_{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}")

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
