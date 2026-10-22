# =============================================================================
#  marker_plots.R  --  library-driven visualisations & oligo ordering helpers
# =============================================================================
#  No hand-rolled geometry: every panel is drawn by a dedicated library.
#
#    - gggenes               : strand-aware arrows (primers on the template,
#                              gene model + variant track)
#    - ggrepel               : non-overlapping labels
#    - OligoArrayAux         : DNA duplex / hairpin free energies (via oligo_dg())
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

#' A target's allele set as a label (`REF>ALT` or `REF>ALT1/ALT2`)
#'
#' A multi-allelic target carries its extra ALTs in an `alt2` column (or a
#' list-column `alt`), so the label - and the plot - show every allele, not just
#' the one the ALT primer was built for.
#' @keywords internal
#' @noRd
snp_alleles <- function(snp) {
  extra <- if ("alt2" %in% names(snp)) as.character(snp$alt2) else character(0)
  a <- unique(c(as.character(unlist(snp$alt)), extra))
  a[!is.na(a) & nzchar(a)]
}

#' @keywords internal
#' @noRd
snp_label <- function(snp) {
  a <- snp_alleles(snp)
  glue::glue("{snp$ref[1L]}>{paste(a, collapse = '/')}")
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
#    discrimination the placement buys, the in-silico PCR products.
#  Both levels are coloured the same way: ok (green) vs broken (red).

MARKER_STATUS <- c(ok = "#1b7837", broken = "#d7301f", neutral = "#4477aa")
# historical internal name; all schemes now share the same palette
KASP_STATUS <- MARKER_STATUS

#' A two-lane schematic gel (shared by the CAPS / ARMS readouts)
#'
#' Bands are horizontal marks at their fragment size; the caller supplies the
#' per-band colour (diagnostic vs shared) so the two marker families render the
#' same way.  Internal - the exported accessors are `plot_caps_gel()` and
#' `plot_arms_gel()`.
#' @param bands tibble with `x` (numeric lane), `size` (bp) and `col`.
#' @param x_labels lane labels.
#' @param y_label y-axis label.
#' @param title,subtitle plot text.
#' @keywords internal
#' @noRd
# One gel vocabulary for every restriction/ARMS family.  Allele identity lives
# in the lane label; colour is reserved for status, never for REF/ALT.
GEL_DIAG <- "#44403C"
GEL_SHARED <- "#D6D3D1"

# Poisson counting noise (CV = 1/sqrt(lambda)) is an illustrative band-width
# cue, not a gel simulator: loading, staining, diffusion and camera response are
# not modelled.
poisson_band_halfwidth <- function(size) {
  lambda <- pmax(4, size / 25)
  cv <- 1 / sqrt(lambda)
  pmin(0.40, 0.22 + 0.30 * cv)
}

# The only constructor for gel bands.  Callers supply biological facts
# (allele, size, diagnostic); the renderer supplies lane, palette and width.
gel_spec <- function(allele, size, diagnostic) {
  n <- length(size)
  if (length(allele) == 1L && n > 1L) allele <- rep(allele, n)
  if (length(diagnostic) == 1L && n > 1L) diagnostic <- rep(diagnostic, n)
  tibble::tibble(
    allele = as.character(allele), size = as.numeric(size),
    diagnostic = as.logical(diagnostic)
  ) |>
    dplyr::mutate(
      x = as.numeric(factor(.data$allele, levels = unique(.data$allele))),
      width = poisson_band_halfwidth(.data$size),
      col = dplyr::if_else(.data$diagnostic, GEL_DIAG, GEL_SHARED)
    )
}

gel_plot <- function(bands, x_labels, y_label, title, subtitle) {
  # `width` is optional: callers that model band counting set it; the plain
  # two-lane gel keeps the shared default footprint.
  if (!"width" %in% names(bands)) {
    bands$width <- 0.32
  }
  ggplot2::ggplot(bands) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$x - .data$width, xend = .data$x + .data$width,
        y = .data$size, yend = .data$size
      ),
      linewidth = 2.4, lineend = "round", colour = bands$col
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(x_labels), labels = x_labels,
      limits = c(0.4, length(x_labels) + 0.6)
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = NULL, y = y_label, title = title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 10)
}

#' Feature-level acceptance criteria (everything else is alignment-level)
#' @keywords internal
#' @noRd
check_level <- function(criterion) {
  site <- c(
    "target site class", "ALT alleles (n)", "linked variants skipped (n)",
    "ALT primer covers all ALTs", "design relies on linkage"
  )
  specificity <- c("off-target exact primer hits (n)")
  feature <- c(
    "AS Tm (C)", "COM Tm (C)", "AS GC", "COM GC",
    "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
    "inverted repeat (arms>=4, loop<=6)", "worst hairpin dG37",
    "worst cross-dimer dG37", "worst homodimer dG37",
    "AS 3' window variants (n)", "AS footprint variants (n)", "COM footprint variants (n)",
    "AS 3' window dropout", "AS footprint dropout", "COM footprint dropout"
  )
  out <- rep("alignment", length(criterion))
  out[criterion %in% feature] <- "feature"
  out[criterion %in% site] <- "site"
  out[criterion %in% specificity] <- "specificity"
  out
}

# The one label -> short table. Keys are the registry labels from
# marker_criteria(); values are the tumbler text. Tests guard that every
# registered check is present, so a short name is never invented at render time.
CHECK_SHORT_MAP <- c(
  "target site class" = "site class",
  "ALT alleles (n)" = "ALT n",
  "linked variants skipped (n)" = "linked skip",
  "ALT primer covers all ALTs" = "ALT cover",
  "design relies on linkage" = "linkage",
  "off-target exact primer hits (n)" = "off-target",
  "amplicon length (bp)" = "amplicon",
  "AS Tm (C)" = "AS Tm",
  "COM Tm (C)" = "COM Tm",
  "AS GC" = "AS GC",
  "COM GC" = "COM GC",
  "COM Tm >= AS Tm" = "COM\u2265AS Tm",
  "longest homopolymer (bp)" = "run",
  "SSR (2-4 nt, >=3 copies)" = "SSR",
  "inverted repeat (arms>=4, loop<=6)" = "inv. repeat",
  "worst hairpin dG37" = "hairpin",
  "worst cross-dimer dG37" = "dimer",
  "worst homodimer dG37" = "homodimer",
  "AS 3' window variants (n)" = "AS 3' var",
  "AS footprint variants (n)" = "AS foot var",
  "COM footprint variants (n)" = "COM var",
  "AS 3' window dropout" = "AS 3' drop",
  "AS footprint dropout" = "AS foot drop",
  "COM footprint dropout" = "COM drop",
  "discrimination" = "discrim.",
  "in-silico PCR (REF 1 / ALT 0)" = "PCR",
  "restriction site defined" = "site present",
  "restriction site separates alleles" = "site discrim.",
  "band gap (bp)" = "band gap",
  "digest sums to amplicon" = "digest sum",
  "restriction site conserved (population)" = "site conserv.",
  "FIP 3' ALT" = "FIP 3' ALT",
  "inner arm Tm balance (C)" = "arm Tm bal.",
  "inner arms inside amplicon" = "arms inside",
  "variant at FIP 3' end" = "FIP variant",
  "primer length (nt)" = "primer len",
  "THF off variant" = "THF off var",
  "RPA probe spans variant" = "RPA probe",
  "variant under primer/probe" = "variant site",
  "PAM valid" = "PAM",
  "variant in seed" = "seed",
  "guide spans variant" = "guide span",
  "guide/PAM conserved (population)" = "guide cons.",
  "off-target protospacers (n)" = "off-target",
  "probe spans variant" = "probe span",
  "probe Tm - primer Tm (C)" = "probe dTm",
  "ASP 3' discrimination" = "ASP discr.",
  "extension 3' at variant" = "ext at var",
  "variant-read distance (bp)" = "read dist",
  "read primers inside amplicon" = "read inside",
  "read primer pair defined" = "read pair",
  "repeat unit (nt)" = "repeat unit",
  "repeat copies" = "copies",
  "REF/ALT size difference (bp)" = "size diff",
  "fragment size (bp)" = "fragment"
)

#' Short labels for the acceptance criteria (for the tumbler column)
#' @keywords internal
#' @noRd
check_short <- function(criterion) {
  out <- unname(CHECK_SHORT_MAP[criterion])
  miss <- is.na(out) | !nzchar(out)
  out[miss] <- criterion[miss]
  out
}

#' Variant points on a track, sized by population frequency when present
#'
#' Weights are optional; when the variant table carries `pop_freq` the markers
#' grow with it (compressed into a modest range so they never shout), otherwise
#' every variant is a small neutral dot.
#' @keywords internal
#' @noRd
variant_points <- function(variants, y) {
  if (!"consequence" %in% names(variants)) {
    variants$consequence <- "variant"
  }
  variants$consequence <- as.character(variants$consequence)
  variants$consequence[is.na(variants$consequence)] <- "variant"
  if ("pop_freq" %in% names(variants)) {
    list(
      ggplot2::geom_point(
        data = variants,
        ggplot2::aes(x = .data$pos, y = y, colour = .data$consequence, size = .data$pop_freq)
      ),
      ggplot2::scale_size_continuous(range = c(0.9, 2.6), guide = "none")
    )
  } else {
    ggplot2::geom_point(
      data = variants,
      ggplot2::aes(x = .data$pos, y = y, colour = .data$consequence), size = 1.3
    )
  }
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
    # the drawn letters: the AS-alt oligo is the physical (possibly degenerate)
    # primer; `seq` stays concrete so the SSR/palindrome QC below sees real bases
    seq_disp = c(p[["AS_ref"]], as_alt_sequence(assay), p[["COM"]]),
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
  d$status <- dplyr::if_else(tm_ok & dplyr::between(d$gc, limits$gc[1], limits$gc[2]) & clean, "ok", "broken")
  win <- purrr::pmap(list(d$start, d$end, d$forward), clamp_window)
  d$cstart <- vapply(win, `[`, numeric(1), 1L)
  d$cend <- vapply(win, `[`, numeric(1), 2L)
  # the primer string is always 5'->3'; the 3' base is its last letter
  last <- nchar(d$seq)
  d$three_x <- dplyr::if_else(d$forward, d$end, d$start)
  d$clamp <- dplyr::if_else(substr(d$seq, last, last) %in% c("G", "C"), "ok", "broken")
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
  spans <- purrr::map(units, function(u) stringr::str_locate_all(s, repeat_regex(u, min_copies))[[1L]])
  purrr::reduce(spans, rbind, .init = matrix(integer(0), ncol = 2L))
}

#' Per-base letters of the assay primers, styled for KASP
#'
#' The primer sequence as one row of letters at genomic coordinates: neutral
#' black by default, **bold** over the 3' GC-clamp window, and bright red on the
#' bases of an SSR or a hairpin stem.  Coordinates follow the strand, so a
#' minus-strand primer reads right-to-left.
#' @keywords internal
#' @noRd
# The letters layer, reusable by any oligo with a plus-strand interval.  The
# sequence is 5'->3'; a minus-strand oligo is laid right-to-left.
oligo_letters <- function(seq, start, end, forward, role, seq_disp = seq, clamp_n = 5L,
                          spans = NULL, bases = NULL) {
  s <- as.character(seq)[1L]
  s_disp <- as.character(seq_disp)[1L]
  L <- nchar(s)
  xs <- if (isTRUE(forward)) seq.int(start, length.out = L) else seq.int(end, by = -1L, length.out = L)
  face <- rep("plain", L)
  colour <- rep("#292524", L)
  clamp <- if (clamp_n > 0L) (L - clamp_n + 1L):L else integer(0)
  clamp <- clamp[clamp >= 1L & clamp <= L]
  face[clamp] <- "bold"
  colour[clamp] <- "black"
  # A displayed oligo may carry a non-DNA chemistry (e.g. the RPA THF `X`).
  # Show the character, but do not ask the DNA secondary-structure code about it.
  if (grepl("^[ACGT]+$", s)) {
    sp <- repeat_spans(s)
    pal <- find_palindromes(s)[[1L]]
    if (nrow(pal)) sp <- rbind(sp, as.matrix(pal[, c("start", "end")]))
    if (nrow(sp)) {
      idx <- unique(purrr::list_c(purrr::map(seq_len(nrow(sp)), function(j) sp[j, 1L]:sp[j, 2L])))
      colour[idx[idx >= 1L & idx <= L]] <- "#d7191c"
    }
  }
  # explicit highlights win over the defaults: any scheme can mark a PAM, a
  # probe base, a variant or a guide without growing this function
  if (!is.null(spans)) {
    for (i in seq_len(nrow(spans))) {
      idx <- seq.int(spans$from[i], spans$to[i])
      idx <- idx[idx >= 1L & idx <= L]
      if (!length(idx)) next
      if ("colour" %in% names(spans) && !is.na(spans$colour[i])) colour[idx] <- spans$colour[i]
      if ("face" %in% names(spans) && !is.na(spans$face[i])) face[idx] <- spans$face[i]
    }
  }
  if (!is.null(bases)) {
    idx <- as.integer(bases$pos)
    idx <- idx[idx >= 1L & idx <= L]
    if (length(idx)) {
      if ("colour" %in% names(bases)) colour[idx] <- bases$colour[match(idx, as.integer(bases$pos))]
      if ("face" %in% names(bases)) face[idx] <- bases$face[match(idx, as.integer(bases$pos))]
    }
  }
  tibble::tibble(
    role = rep(as.character(role), L), x = xs, letter = strsplit(s_disp, "")[[1L]],
    face = face, colour = colour
  )
}

letter_size_for_feature <- function(width, max_size = LETTER_SIZE_DEFAULT) {
  # a short block (PAM, restriction site) cannot use the amplicon-scale font:
  # shrink the glyphs so the letters stay inside the block
  min(max_size, 0.45 + 0.11 * width)
}

letters_layer <- function(d, size = LETTER_SIZE_DEFAULT) {
  if (!"y" %in% names(d)) d$y <- d$role
  if ("size" %in% names(d)) {
    return(list(
      ggplot2::geom_text(
        data = d,
        ggplot2::aes(
          x = .data$x, y = .data$y, label = .data$letter,
          fontface = .data$face, colour = .data$colour, size = .data$size
        ),
        family = "mono", show.legend = FALSE
      ),
      ggplot2::scale_size_identity()
    ))
  }
  ggplot2::geom_text(
    data = d,
    ggplot2::aes(
      x = .data$x, y = .data$y, label = .data$letter,
      fontface = .data$face, colour = .data$colour
    ),
    family = "mono", size = size, show.legend = FALSE
  )
}

# One sizing rule for every sequence row.  A primer is not drawn at a fixed
# size: it shares the amplicon/reference scale, so a 600 bp Sanger product and
# its 20 bp primers cannot end up in two different typographic worlds.
LETTER_SIZE_DEFAULT <- 2.6
letter_size_for_span <- function(span) {
  n <- diff(range(span)) + 1L
  if (n > 500L) {
    0.9
  } else if (n > 250L) {
    1.1
  } else if (n > 150L) {
    1.5
  } else {
    LETTER_SIZE_DEFAULT
  }
}

primer_letters <- function(assay, clamp_n = 5L) {
  d <- primer_features(assay)
  rows <- lapply(seq_len(nrow(d)), function(i) {
    oligo_letters(
      d$seq[i], d$start[i], d$end[i], d$forward[i], d$role[i],
      seq_disp = d$seq_disp[i], clamp_n = clamp_n
    )
  })
  dplyr::bind_rows(rows)
}

#' Ideogram inputs: chromosome bars + exact primer hits + in-silico PCR products
#'
#' Computed once (at design time) and carried on the assay (`assay$ideogram`), so
#' plotting reads the object instead of rescanning the genome.  When `hits` or
#' `products` are supplied (the design-time virtual PCR, see
#' [assay_pcr_products()]) they are reused instead of rescanned.
#' @keywords internal
#' @noRd
assay_ideogram_data <- function(assay, genome, BPPARAM = NULL, hits = NULL, products = NULL) {
  if (is.null(genome)) {
    return(NULL)
  }
  chroms <- genome_chroms(genome)
  bars <- tibble::tibble(
    chrom = chroms,
    length = vapply(chroms, function(ch) genome_length(genome, ch), numeric(1))
  )
  roles <- c(AS_ref = "AS-ref", AS_alt = "AS-alt", COM = "COM")
  hits <- hits %||% primer_hits(assay_primers(assay), genome, BPPARAM = BPPARAM)
  if (nrow(hits)) hits$signal <- roles[hits$primer]
  products <- products %||%
    assay_pcr_products(assay, genome, BPPARAM = BPPARAM)
  prods <- purrr::imap(products, function(p, nm) {
    if (nrow(p)) dplyr::mutate(p, signal = roles[[nm]]) else NULL
  }) |>
    purrr::list_rbind()
  list(bars = bars, hits = hits, products = prods)
}

#' Where every primer matches the genome (positions), for the ideogram
#'
#' The scan is split by chromosome when a `BPPARAM` is given.
#' @keywords internal
#' @noRd
primer_hits <- function(primers, genome, BPPARAM = NULL) {
  scan_one <- function(ch) {
    worker_biostrings()
    subject <- Biostrings::DNAString(genome_seq(genome, ch))
    purrr::imap(primers, function(p, nm) {
      # a primer binds either strand: search the plus strand for the primer and
      # for its reverse complement, and report the plus-strand match starts.
      # matchPattern(fixed = FALSE) keeps IUPAC degenerate primers usable.
      pats <- unique(c(p, as.character(Biostrings::reverseComplement(Biostrings::DNAString(p)))))
      pos <- purrr::list_c(purrr::map(pats, function(pt) {
        # a genomic N must never satisfy a primer: `fixed = FALSE` lets a subject
        # N match anything, which turns every N run into thousands of hits
        IRanges::start(exact_pattern_hits(pt, subject, 0L))
      }))
      if (!length(pos)) {
        return(NULL)
      }
      tibble::tibble(primer = nm, chrom = ch, pos = sort(unique(pos)))
    }) |>
      purrr::list_rbind()
  }
  chroms <- genome_chroms(genome)
  if (is.null(BPPARAM)) {
    purrr::list_rbind(lapply(chroms, scan_one))
  } else {
    purrr::list_rbind(bp_lapply(chroms, scan_one, BPPARAM))
  }
}

#' Genome ideogram: exact primer hits + in-silico PCR products
#'
#' A thin, schematic view of the whole genome: every chromosome is drawn as a
#' grey bar on its own row (scaled to fill the width, karyotype-style) and two
#' signals are overlaid - where each assay primer matches exactly
#' (`primer_hits()`) and where the designed allele-specific pairs actually
#' amplify ([in_silico_pcr()], the same pairings as [assay_amplicons()]).  The
#' target locus is marked.
#'
#' The underlying data are computed once at design time and carried on the assay
#' (`assay$ideogram`), so plotting only reads the object; when it is absent (an
#' older assay) and a `genome` is supplied they are recomputed.  Returns `NULL`
#' when neither is available, so [plot_primers()] simply drops the panel.
#' @keywords internal
#' @noRd
assay_ideogram <- function(assay, snp, genome = NULL) {
  dat <- assay$ideogram %||% assay_ideogram_data(assay, genome)
  if (is.null(dat)) {
    return(NULL)
  }
  bars <- dat$bars
  hits <- dat$hits
  prods <- dat$products
  sig_cols <- c("AS-ref" = "#4477aa", "AS-alt" = "#66c2a5", "COM" = "#8da0cb")
  extra <- setdiff(unique(hits$signal), names(sig_cols))
  if (length(extra)) {
    more <- c("#e78ac3", "#a6d854", "#ffd92f", "#fc8d62", "#66c2a5")
    sig_cols <- c(sig_cols, stats::setNames(rep(more, length.out = length(extra)), extra))
  }
  # adaptive axis unit: Mb for real genomes, kb/bp for small ones
  span <- max(bars$length)
  div <- if (span >= 1e6) 1e6 else if (span >= 1e3) 1e3 else 1
  suf <- if (span >= 1e6) " Mb" else if (span >= 1e3) " kb" else " bp"

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
#  Gene / primer maps (gggenes)
# =============================================================================

#' Normalise an annotation table to the feature row's shape
#'
#' Accepts a data frame with `start`/`end` (and optional `strand`, `name`/`symbol`)
#' or a gene record, keeps only features overlapping `span`, and returns
#' `start`/`end`/`forward`/`name` (or `NULL`).  One place so the doomed and the
#' designed views cannot drift.
#' @keywords internal
#' @noRd
normalise_features <- function(features, span) {
  if (is.null(features)) {
    return(NULL)
  }
  features <- tibble::as_tibble(features)
  if (!nrow(features) || !all(c("start", "end") %in% names(features))) {
    return(NULL)
  }
  nm <- if ("name" %in% names(features)) {
    features$name
  } else if ("symbol" %in% names(features)) {
    features$symbol
  } else {
    NA_character_
  }
  feat <- tibble::tibble(
    start = as.integer(features$start), end = as.integer(features$end),
    forward = if ("strand" %in% names(features)) features$strand != "-" else TRUE,
    name = nm
  )
  feat <- feat[feat$end >= span[1L] & feat$start <= span[2L], , drop = FALSE]
  if (nrow(feat) == 0L) NULL else feat
}

#' Context plot for a doomed target
#'
#' A target no assay clears the gates for still deserves a picture: the local
#' window, its variants and the feature annotation, on the same x axis as
#' [plot_primers()], with the reason it was doomed.  Drawn by [plot_primers()]
#' when it is handed a `marker_doomed` object.
#' @keywords internal
#' @noRd
plot_doomed_target <- function(assay, variants, snp, features = NULL) {
  win <- design_window(assay, snp)
  vw <- variants[
    variants$chrom == snp$chrom & variants$pos >= win[1L] & variants$pos <= win[2L], ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = tibble::tibble(x = win[1L], xend = win[2L]),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = 1, yend = 1),
      colour = "grey72", linewidth = 0.3
    )
  feat <- normalise_features(features, win)
  if (is.null(feat)) {
    p <- p + ggplot2::annotate("text",
      x = mean(win), y = 1, label = "no feature annotation",
      size = 2.1, colour = "grey55", fontface = "italic"
    )
  } else {
    p <- p + ggplot2::geom_rect(
      data = feat,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = 0.9, ymax = 1.1),
      fill = "grey72", colour = "grey40"
    )
    lab <- feat[!is.na(feat$name), , drop = FALSE]
    if (nrow(lab)) {
      lab$x <- (pmax(lab$start, win[1L]) + pmin(lab$end, win[2L])) / 2
      p <- p + ggplot2::geom_text(
        data = lab, ggplot2::aes(x = .data$x, y = 1, label = .data$name),
        size = 2.1, colour = "black", fontface = "bold"
      )
    }
  }

  p +
    ggplot2::geom_point(data = vw, ggplot2::aes(.data$pos, 0.45), colour = "grey65", size = 0.8) +
    ggplot2::geom_point(
      data = snp, ggplot2::aes(.data$pos, 0.45),
      shape = 21, fill = "firebrick", colour = "white", size = 2.6
    ) +
    ggplot2::geom_vline(xintercept = snp$pos, linetype = "dashed", colour = "grey45") +
    ggplot2::coord_cartesian(xlim = win, ylim = c(0.2, 1.25)) +
    ggplot2::labs(
      x = glue::glue("{snp$chrom} position (bp)"), y = NULL,
      title = glue::glue("Doomed target {snp$chrom}:{snp$pos} {snp_label(snp)}"),
      subtitle = glue::glue("no assay passed the gates - {assay$reason}")
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(colour = "#d7301f", face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
}

#' Layered haplotype track: one row per haplotype, its variants and frequency
#'
#' Drawn only when the assay carries a haplotype table.  Each row is a haplotype,
#' labelled with its population frequency and the AS primer it is expected to
#' bind; its own variants are marked (the target in red), so it is visible which
#' haplotype a primer catches and which variant rides with it.
#' @keywords internal
#' @noRd
haplotype_layers <- function(assay, snp, span) {
  h <- assay$haplotypes
  if (is.null(h) || !nrow(h)) {
    return(NULL)
  }
  haps <- unique(h$haplotype)
  h$y <- match(h$haplotype, haps)
  vari <- h[h$pos != snp$pos, , drop = FALSE]
  tgt <- tibble::tibble(x = snp$pos, y = seq_along(haps))
  hits <- assay$haplotype_hits
  lab <- if (is.null(hits)) {
    haps
  } else {
    hits$y <- match(hits$haplotype, haps)
    glue::glue("{hits$haplotype} {round(hits$freq, 2)} -> {hits$primer}")
  }
  base <- tibble::tibble(x = span[1L], xend = span[2L], y = seq_along(haps))
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = base,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      colour = "grey92"
    ) +
    ggplot2::geom_point(
      data = vari, ggplot2::aes(.data$pos, .data$y),
      shape = 21, fill = "#e6ab02", colour = "white", size = 1.8
    ) +
    ggplot2::geom_text(
      data = vari, ggplot2::aes(.data$pos, .data$y, label = .data$alt),
      vjust = -0.9, size = 1.9, family = "mono", colour = "#7f6000"
    ) +
    ggplot2::geom_point(
      data = tgt, ggplot2::aes(.data$x, .data$y),
      shape = 21, fill = "firebrick", colour = "white", size = 2.0
    ) +
    ggplot2::geom_vline(xintercept = snp$pos, linetype = "dashed", colour = "grey65") +
    ggplot2::scale_y_reverse(breaks = seq_along(haps), labels = lab) +
    ggplot2::coord_cartesian(xlim = span) +
    ggplot2::labs(x = NULL, y = NULL, title = "haplotypes (freq -> AS primer)") +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black")
    )
}

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
#' @param assay assembled assay, or a `marker_doomed` object (then only the target
#'   context and the reason are drawn).
#' @param variants variant table.
#' @param snp target variant (one row).
#' @param genome optional `marker_genome` or `DNAStringSet`; when given, the
#'   genome-wide ideogram panel is added.
#' @param features optional annotation to draw on the feature row: a data frame
#'   with `start`/`end` (and optional `strand`, `name`), or a gene record from
#'   [ensembl_gene()].  A KASP marker should sit on an annotated feature, so the
#'   row says `"no feature annotation"` when nothing (or no overlap) is given.
#' @return a `patchwork` object: the local map (primers over the variant track,
#'   sharing an x axis) beside the tumbler list from [plot_checks()].  For a
#'   `marker_doomed` target, a single context plot with the reason instead.
#' @examples
#' \dontrun{
#' plot_primers(assay, variants, snp)
#' }
#' @family KASP plots
#' @export
plot_primers <- function(assay, variants, snp, genome = NULL, features = NULL) {
  if (is_kasp_doomed(assay)) {
    return(plot_doomed_target(assay, variants, snp, features))
  }
  win <- design_window(assay, snp)
  aa <- assay$assay[1, ]
  strand <- aa$strand %||% "+"
  plus <- identical(strand, "+")
  d <- primer_features(assay)
  letters <- primer_letters(assay)

  # alignment-level arrow: the amplicon the placement buys, on its own row
  amp_ok <- dplyr::between(aa$amplicon, KASP_LIMITS$amplicon[1], KASP_LIMITS$amplicon[2])
  amp <- tibble::tibble(
    role = factor("amplicon", levels = c(levels(d$role), "amplicon")),
    start = min(d$start), end = max(d$end), forward = plus,
    status = dplyr::if_else(amp_ok, "ok", "broken")
  )
  # one size for the whole map: primers and amplicon share the reference scale
  letter_size <- letter_size_for_span(c(amp$start, amp$end))
  # the amplicon's own sequence (plus strand) as letters on its row
  ti <- function(p) p - snp$pos + assay$template$snp_i
  amp_seq <- as.character(Biostrings::subseq(assay$template$seq, ti(amp$start), ti(amp$end)))
  letters <- dplyr::bind_rows(letters, tibble::tibble(
    role = "amplicon", x = seq.int(amp$start, amp$end),
    letter = strsplit(amp_seq, "")[[1L]], face = "plain", colour = "#292524"
  ))

  # the Tm as a small number just past the 3' end, and its deviation from the
  # target as a green/red tumbler next to it
  lab <- d
  lab$target <- dplyr::if_else(lab$role == "COM", COM_TM_TARGET, AS_TM_TARGET)
  lab$dtm <- lab$tm - lab$target
  lab$meas_x <- dplyr::if_else(lab$forward, lab$end + 1, lab$start - 1)
  lab$hj <- dplyr::if_else(lab$forward, 0, 1)
  lab$dtm_x <- dplyr::if_else(lab$forward, lab$meas_x + 4, lab$meas_x - 4)
  lab$dtm_status <- dplyr::if_else(abs(lab$dtm) <= 3, "ok", "broken")

  # pins: the 3' discriminating base of every oligo, and the injected
  # destabilising mismatch of the AS primers (if the placement carries one)
  pins <- tibble::tibble(role = as.character(d$role), x = d$three_x)
  asx <- d[!is.na(d$stab), ]
  if (nrow(asx)) {
    pins <- dplyr::bind_rows(pins, tibble::tibble(
      role = as.character(asx$role),
      x = dplyr::if_else(asx$forward, asx$end - (asx$stab - 1L), asx$start + (asx$stab - 1L))
    ))
  }

  span <- range(c(d$start, d$end, lab$meas_x, lab$dtm_x, amp$start, amp$end))

  # the annotated features on the same axis: a KASP marker should sit on a gene,
  # not in intergenic space, so if no (overlapping) annotation is given we say so
  feat <- normalise_features(features, span)

  p_map <- ggplot2::ggplot() +
    gggenes::geom_gene_arrow(
      data = amp,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$role, forward = .data$forward),
      fill = "grey94", colour = "grey55", linewidth = 0.25
    ) +
    gggenes::geom_gene_arrow(
      data = d,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$role, forward = .data$forward),
      fill = "grey98", colour = "grey75", linewidth = 0.25
    ) +
    letters_layer(letters, size = letter_size) +
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
  # the AS 3'-terminal window and the COM footprint: the two zones a variant
  # must not sit in.  A zone is "clean" when the expected dropout stays within
  # tolerance (DECIPHER), falling back to "no known variant" when the estimate is
  # unavailable.  Same green/red logic for both.
  tol <- KASP_LIMITS$footprint_dropout
  zone_status <- function(drop, n_var) {
    if (length(drop) && !is.na(drop)) {
      if (drop <= tol) "ok" else "broken"
    } else if (isTRUE(n_var == 0L)) {
      "ok"
    } else {
      "broken"
    }
  }
  win_status <- zone_status(aa[["as_dropout"]], aa$as_n_3p)
  com_status <- zone_status(aa[["com_dropout"]], aa$n_var)

  # the local variants: a thin strip on the same x axis.  When the assay carries
  # a linked haplotype, the target + its co-segregating bases are drawn as two
  # parallel tracks (REF vs ALT haplotype), so the linkage is visible.
  lk <- assay$selection$linked
  # when a full haplotype table is present the dedicated layered track carries
  # the linkage, so the two-row REF/ALT view is not drawn twice
  has_link <- !is.null(lk) && nrow(lk) > 0 && is.null(assay$haplotypes)
  if (has_link) {
    gp <- snp$pos - assay$template$snp_i + lk$offset
    hap <- dplyr::bind_rows(
      tibble::tibble(x = snp$pos, y = 1, base = snp$ref, kind = "ref_target"),
      tibble::tibble(x = snp$pos, y = 2, base = snp$alt, kind = "alt_target"),
      tibble::tibble(x = gp, y = 1, base = lk$ref, kind = "ref_linked"),
      tibble::tibble(x = gp, y = 2, base = lk$alt, kind = "alt_linked")
    )
    hap_col <- c(
      ref_target = KASP_STATUS[["neutral"]], alt_target = "#d7301f",
      ref_linked = "grey62", alt_linked = "#e6ab02"
    )
    hap_face <- c(ref_target = "bold", alt_target = "bold", ref_linked = "plain", alt_linked = "bold")
    ybreaks <- c(1, 2)
    ylabs <- c("REF hap", "ALT hap")
    zone <- c(0.5, 2.5)
    tall <- c(0.28, 2.72)
  } else {
    ybreaks <- 1
    ylabs <- "SNP"
    zone <- c(0.6, 1.4)
    tall <- c(0.22, 1.78)
  }
  ys <- ybreaks
  shared <- tibble::tibble(
    pos = rep(vw$pos, times = length(ys)),
    y = rep(ys, each = nrow(vw))
  )
  foot_pts <- tibble::tibble(
    pos = rep(foot$pos, times = length(ys)),
    y = rep(ys, each = nrow(foot))
  )
  p_var <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
      xmin = com_span[1] - 0.5, xmax = com_span[2] + 0.5,
      ymin = zone[1], ymax = zone[2], fill = KASP_STATUS[[com_status]], alpha = 0.22
    ) +
    # the AS 3' window is the discriminating side: drawn as a taller square so it
    # reads as the important zone, coloured by its own dropout
    ggplot2::annotate("rect",
      xmin = reg3p[1] - 0.5, xmax = reg3p[2] + 0.5,
      ymin = tall[1], ymax = tall[2], fill = KASP_STATUS[[win_status]], alpha = 0.30
    ) +
    ggplot2::geom_point(
      data = shared, ggplot2::aes(.data$pos, .data$y),
      colour = "grey70", size = 0.8
    ) +
    ggplot2::geom_point(data = foot_pts, ggplot2::aes(.data$pos, .data$y), shape = 1, colour = "black", size = 2.1)
  if (has_link) {
    p_var <- p_var +
      ggplot2::geom_point(
        data = hap[hap$kind %in% c("ref_target", "alt_target"), ],
        ggplot2::aes(.data$x, .data$y, fill = .data$kind), shape = 21, colour = "white", size = 2.6
      ) +
      ggplot2::geom_text(
        data = hap, ggplot2::aes(.data$x, .data$y, label = .data$base, colour = .data$kind, fontface = .data$kind),
        vjust = -0.9, size = 2.3, family = "mono", show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = hap_col, guide = "none") +
      ggplot2::scale_colour_manual(values = hap_col, guide = "none") +
      ggplot2::scale_discrete_manual("fontface", values = hap_face, guide = "none")
  } else {
    p_var <- p_var +
      variant_points(vw, 1) +
      ggplot2::geom_point(
        data = snp, ggplot2::aes(.data$pos, 1),
        shape = 21, fill = "firebrick", colour = "white", size = 2.4
      ) +
      ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL)
    # a multi-allelic target: name every ALT at the site, so the track reads as
    # one *site* with several alleles, not one allele
    aa <- snp_alleles(snp)
    if (length(aa) > 1L) {
      p_var <- p_var + ggplot2::annotate("text",
        x = snp$pos, y = 1, label = paste(snp$ref[1L], paste(aa, collapse = "/"), sep = ">"),
        vjust = 2.1, size = 2.2, family = "mono", colour = "black"
      )
    }
  }
  p_var <- p_var +
    ggplot2::scale_y_continuous(breaks = ybreaks, labels = ylabs, limits = c(zone[1] - 0.15, zone[2] + 0.15)) +
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
  p_hap <- haplotype_layers(assay, snp, span)

  # primer rows are packed tight; the feature row sits at the bottom
  parts <- list(p_map, p_var)
  heights <- c(2.3, 0.4)
  if (!is.null(p_hap)) {
    parts <- c(parts, list(p_hap))
    heights <- c(heights, 0.22 * length(unique(assay$haplotypes$haplotype)))
  }
  if (!is.null(p_ideo)) {
    parts <- c(parts, list(p_ideo))
    heights <- c(heights, 3)
  }
  panels <- patchwork::wrap_plots(parts, ncol = 1, heights = heights)

  patchwork::wrap_plots(panels, plot_checks(assay), widths = c(1, 0.12)) +
    patchwork::plot_annotation(
      title = "KASP assay: local map + acceptance checks",
      subtitle = glue::glue(
        "target {snp$chrom}:{snp$pos} ({snp_label(snp)}) on strand {strand}"
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
  out <- if (!is.numeric(value)) {
    as.character(value)
  } else {
    dplyr::case_when(
      name %in% c("AS GC", "COM GC", "discrim.", "AS drop", "COM drop") ~ sprintf("%.0f%%", 100 * value),
      name %in% c("AS Tm", "COM Tm", "Tm", "dimer", "homodimer", "hairpin", "COM\u2265AS Tm") ~ sprintf("%.1f", value),
      abs(value) >= 1e6 ~ paste0(round(value / 1e6), "M"),
      abs(value) >= 1e3 ~ paste0(round(value / 1e3), "K"),
      TRUE ~ as.character(round(value))
    )
  }
  out[is.na(value)] <- "-"
  out
}

# One tumbler for KASP acceptance checks and for the per-class gate tables.
# `chk` needs `level`, `name`, `status` (ok/broken/neutral) and `value`.
tumbler_plot <- function(chk) {
  chk$level <- droplevels(factor(chk$level, levels = unique(chk$level)))
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

#' Acceptance checks as a thin list of tumblers (feature vs alignment)
#'
#' Runs [marker_design_checks()] and shows one named tumbler per criterion - a
#' small green (ok) / red (broken) marker with its observed value and limit - in
#' a single narrow column, faceted by the level judged: the *site* level is the
#' target itself (KASP class, ALT count, linked variants), the *feature* level is
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
  chk <- assay$evaluation %||% marker_design_checks(assay, limits)
  if (!"status" %in% names(chk)) {
    chk$status <- dplyr::if_else(chk$pass, "pass", "fail")
  }
  chk <- chk |>
    dplyr::mutate(
      level = factor(.data$level, levels = c("site", "feature", "alignment", "specificity")),
      name = check_short(.data$criterion),
      status = dplyr::case_when(
        .data$status %in% c("unknown", "not_applicable") ~ "neutral",
        .data$status == "pass" ~ "ok",
        TRUE ~ "broken"
      ),
      value = dplyr::coalesce(.data$detail, format_check_value(.data$name, .data$observed))
    )
  # the site class is shown as a 4-char code (monospace column, max 4 wide)
  cls <- c(biallelic = "bial", degenerate = "degn", complementary = "comp", indel = "indl", none = "--")
  is_site <- chk$criterion == "target site class"
  chk$value[is_site] <- unname(cls[chk$detail[is_site]])
  # population facts are chemistry-independent: append ALT frequency and, when
  # available, haplotype coverage to the same check column
  pop <- marker_population_checks(assay)
  if (nrow(pop)) chk <- dplyr::bind_rows(chk, pop)
  tumbler_plot(chk)
}

#' Points of a constant-height arc (bump) spanning `[a, b]` above `y0`
#' @keywords internal
#' @noRd
bump_points <- function(a, b, y0 = 0.7, h = 1, n = 40L) {
  t <- seq(0, pi, length.out = n)
  tibble::tibble(
    x = (a + b) / 2 + (b - a) / 2 * cos(t),
    y = y0 + h * sin(t)
  )
}

#' Hairpins, inverted repeats and SSRs within the assay oligos
#'
#' The sequence-level structure of the three primers side by side: each oligo's
#' bases, its SSR spans underlined in red, and its inverted-repeat / hairpin stems
#' drawn as arcs above the sequence.  The visual twin of the thermodynamic
#' hairpin/SSR screen, so a primer that folds is seen, not merely scored.  Arcs
#' are constant-height bumps (`geom_path()` over computed points) rather than
#' `ggforce::geom_arc()` semicircles: a semicircle's height grows with the stem,
#' which would dominate a compact strip, while the horizontal span already shows
#' the stem length.
#' @param assay assembled assay.
#' @param min_armlength,max_looplength passed to [find_palindromes()].
#' @return a `ggplot` faceted by oligo.
#' @examples
#' \dontrun{
#' plot_oligo_structure(assay)
#' }
#' @family KASP plots
#' @export
plot_oligo_structure <- function(assay, min_armlength = 4L, max_looplength = 6L) {
  p <- assay_primers(assay)
  oligos <- tibble::tibble(oligo = names(p), seq = unname(p), L = nchar(unname(p)))
  pals <- purrr::map(oligos$seq, function(s) {
    find_palindromes(s, min_armlength, max_looplength)[[1L]]
  })
  ssrs <- purrr::map(oligos$seq, repeat_spans)
  lab <- glue::glue(
    "{oligos$oligo}  {purrr::map_int(pals, nrow)} hairpin / {purrr::map_int(ssrs, nrow)} SSR"
  )
  panel_of <- function(oligo) lab[match(oligo, oligos$oligo)]

  # one row of letters per oligo, red where an SSR span sits
  letters <- purrr::pmap(oligos, function(oligo, seq, L) {
    colour <- rep("grey45", L)
    sp <- repeat_spans(seq)
    if (nrow(sp)) {
      idx <- unique(purrr::list_c(purrr::map(seq_len(nrow(sp)), function(j) sp[j, 1L]:sp[j, 2L])))
      colour[idx[idx >= 1L & idx <= L]] <- "#d7191c"
    }
    tibble::tibble(
      panel = panel_of(oligo), oligo = oligo, x = seq_len(L),
      letter = strsplit(seq, "")[[1L]], colour = colour
    )
  }) |>
    purrr::list_rbind()
  letters$panel <- factor(letters$panel, levels = lab)

  ssr <- purrr::imap(ssrs, function(sp, i) {
    if (!nrow(sp)) {
      return(NULL)
    }
    tibble::tibble(panel = lab[i], x = sp[, 1L], xend = sp[, 2L])
  }) |>
    purrr::list_rbind()

  arcs <- purrr::imap(pals, function(h, i) {
    if (!nrow(h)) {
      return(NULL)
    }
    purrr::pmap(h, function(start, end, width) {
      dplyr::mutate(bump_points(start, end), panel = lab[i], arc = paste(lab[i], start, end))
    }) |>
      purrr::list_rbind()
  }) |>
    purrr::list_rbind()

  gg <- ggplot2::ggplot(letters, ggplot2::aes(
    x = .data$x, y = 0, label = .data$letter, colour = .data$colour
  )) +
    ggplot2::geom_text(family = "mono", size = 2.0) +
    ggplot2::scale_colour_identity()
  if (nrow(ssr)) {
    gg <- gg + ggplot2::geom_segment(
      data = ssr,
      ggplot2::aes(x = .data$x - 0.45, xend = .data$xend + 0.45, y = -0.45, yend = -0.45),
      colour = "#d7191c", linewidth = 1.5
    )
  }
  if (nrow(arcs)) {
    gg <- gg + ggplot2::geom_path(
      data = arcs,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$arc),
      colour = "#4477aa", alpha = 0.75, linewidth = 0.4
    )
  }
  gg +
    ggplot2::facet_grid(panel ~ ., scales = "free_x") +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(-0.9, 2.2), breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0, size = 6.5)
    )
}

#' Self-complementarity dot-plot of one oligo
#'
#' Plots every position pair that can Watson-Crick pair inside the same oligo
#' (`i` pairs `j` when base `i` is the complement of base `j`).  A hairpin shows
#' as an off-diagonal block; the closer the block sits to the main diagonal, the
#' tighter its stem/loop.  Sequence-level, so it complements the thermodynamic
#' hairpin free energy with the raw pairing pattern.
#' @param seq one oligo, 5'->3'.
#' @return a `ggplot` dot-plot.
#' @examples
#' plot_self_complementarity("GAAGGTGACCAAGTTCATGCT")
#' @family KASP plots
#' @export
plot_self_complementarity <- function(seq) {
  seq <- as.character(seq)[1L]
  b <- strsplit(seq, "")[[1L]]
  # i pairs j when base i is the complement of base j
  good <- outer(b, chartr("ACGT", "TGCA", b), "==")
  hit <- which(good, arr.ind = TRUE)
  n <- length(b)
  ggplot2::ggplot(tibble::tibble(x = hit[, 1L], y = hit[, 2L]), ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(width = 0.9, height = 0.9, fill = "#4477aa") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.2) +
    ggplot2::coord_equal(xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5)) +
    ggplot2::labs(x = "position", y = "position") +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Gene model with the variant track and a zoom on the design window
#'
#' The context of the target: the gene model and the 1001-Genomes variant track,
#' with `ggforce::facet_zoom()` adding a zoomed panel on the design window and
#' boxing the zoomed area on the full view - the "where is this SNP" picture that
#' sits above the local [plot_primers()] map.
#' @param gene gene record.
#' @param variants variant table.
#' @param snp target variant (one row).
#' @param assay assembled assay.
#' @return a `ggplot` object: the gene model and variant track, faceted by
#'   `ggforce::facet_zoom()` into the full view and a zoom on the design window.
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
  # label the target from the variant itself - a bare (chrom, pos, ref, alt)
  # record carries no `id`, so the plot must not depend on one; a multi-allelic
  # target shows every ALT
  extra_alt <- if ("alt2" %in% names(snp) && !is.na(snp$alt2[1L])) paste0("/", snp$alt2) else ""
  snp <- snp |>
    dplyr::mutate(label = glue::glue("{.data$chrom}:{.data$pos} {.data$ref}>{.data$alt}{extra_alt}"))
  win <- design_window(assay, snp)

  ggplot2::ggplot() +
    gggenes::geom_gene_arrow(
      data = gene_df,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = 1, forward = .data$strand > 0),
      fill = "grey80"
    ) +
    variant_points(v, 2) +
    ggplot2::geom_point(
      data = snp, ggplot2::aes(.data$pos, 2),
      shape = 21, fill = "firebrick", colour = "white", size = 3
    ) +
    ggrepel::geom_text_repel(
      data = snp, ggplot2::aes(.data$pos, 2, label = .data$label),
      size = 2.6, nudge_y = 0.6, min.segment.length = 0
    ) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "consequence") +
    ggplot2::scale_y_continuous(
      limits = c(0.5, 3), breaks = c(1, 2),
      labels = c(gene$symbol, "variants")
    ) +
    ggplot2::labs(
      x = glue::glue("{gene$chrom} position (bp)"), y = NULL,
      title = glue::glue("Variation across {gene$symbol} (1001 Genomes), zoom on the design window")
    ) +
    ggforce::facet_zoom(xlim = win, zoom.size = 1) +
    ggplot2::theme_minimal()
}


# =============================================================================
#  Output helpers
# =============================================================================

#' Write the assay results as a few consolidated files
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
#' \dontrun{
#' write_kasp_results(assay, gene, dir = tempdir())
#' }
#' @family KASP output
#' @export
write_kasp_results <- function(assay, gene, dir = tempdir()) {
  fs::dir_create(dir)
  base <- glue::glue("marker_{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}")

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

#' Tidy table of the assay oligos (with KASP tails)
#' @param assay assembled assay.
#' @param gene gene record.
#' @return a tibble (`assay`, `role`, `oligo_id`, `sequence`, `seq_no_tail`,
#'   `tail`, `tm`, `gc`, `length`) with the three assay oligos and their KASP
#'   tails.
#' @examples
#' \dontrun{
#' marker_primer_table(assay, gene)
#' }
#' @family KASP output
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
#' \dontrun{
#' export_primers(assay, gene, dir = tempdir())
#' }
#' @family KASP output
#' @export
export_primers <- function(assay, gene, dir = tempdir()) {
  fs::dir_create(dir)
  tbl <- marker_primer_table(assay, gene)
  base <- glue::glue("marker_{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}")

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
