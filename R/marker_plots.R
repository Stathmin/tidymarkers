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

#' Keep the `n` most frequent levels, lump the rest into `"Other"`
#'
#' Local stand-in for `forcats::fct_lump_n()`: same level order (kept levels in
#' their original order, `"Other"` last) and the same tie rule (every level
#' tied with the n-th count is kept).  `NA` stays `NA`.
#' @keywords internal
#' @noRd
lump_top_n <- function(x, n, other_level = "Other") {
  f <- factor(x)
  if (nlevels(f) <= n) {
    return(f)
  }
  counts <- as.integer(table(f))
  keep <- if (n <= 0L) {
    character(0)
  } else {
    rank_order <- order(-counts, seq_along(counts))
    levels(f)[counts >= counts[rank_order[[n]]]]
  }
  lumped <- !levels(f) %in% keep
  out <- as.character(f)
  out[!is.na(out) & !out %in% keep] <- other_level
  factor(out, levels = c(levels(f)[!lumped], if (any(lumped)) other_level))
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
    dplyr::mutate(consequence = lump_top_n(.data$consequence, 6L))
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
  "AS 3' mismatch base pair" = "AS 3' pair",
  "AS footprint variants (n)" = "AS foot var",
  "COM footprint variants (n)" = "COM var",
  "AS 3' window dropout" = "AS 3' drop",
  "AS footprint dropout" = "AS foot drop",
  "COM footprint dropout" = "COM drop",
  "discrimination" = "discrim.",
  "in-silico PCR (REF 1 / ALT 0)" = "PCR",
  "in-silico PCR products (n)" = "PCR products",
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
  "fragment size (bp)" = "fragment",
  # §W family check: the scheme applies only to a single-base substitution
  "target is a single-base substitution" = "single base"
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

#' Where every primer matches the genome, from one grouped scan
#'
#' `assay_ideogram_data()` asks `primer_hits()` for the exact hits of every
#' primer.  Searching them one at a time is one genome pass per (primer,
#' strand) - six for a three-primer assay, and the two allele-specific primers
#' are still two patterns, so the scan is one pass per distinct pattern.
#' This searches them in one go (`pattern_hits_batch()`), then shapes the hits
#' as `primer_hits()` returns them: one row per (primer, chromosome) with the
#' sorted, de-duplicated hit starts of that primer on either strand.
#' @param primers named character vector of oligos.
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @return a tibble `primer`, `chrom`, `pos`, or `NULL` when no primer hits.
#' @keywords internal
#' @noRd
batch_primer_hits <- function(primers, genome) {
  if (!length(primers)) {
    return(NULL)
  }
  chroms <- genome_chroms(genome)
  # a primer binds either strand: search every primer and its reverse complement
  # (both strands of both), as `primer_hits()` does
  fwd <- c(primers, stats::setNames(
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(unname(primers)))),
    paste0(names(primers), "_rc")
  ))
  hits <- pattern_hits_batch(unname(fwd), genome, 0L)
  res <- purrr::imap(primers, function(p, nm) {
    back <- fwd[[paste0(nm, "_rc")]]
    purrr::map(chroms, function(ch) {
      pos <- c(
        IRanges::start(hits[[p]][[ch]] %||% IRanges::IRanges()),
        IRanges::start(hits[[back]][[ch]] %||% IRanges::IRanges())
      )
      if (!length(pos)) {
        return(NULL)
      }
      tibble::tibble(primer = nm, chrom = ch, pos = sort(unique(pos)))
    }) |> purrr::list_rbind()
  })
  purrr::list_rbind(res)
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
#' when neither is available, so [plot_primers()] drops the panel.
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
