# =============================================================================
#  marker_plot_methods.R  --  every marker class can draw itself
# =============================================================================
#  `plot(marker)` dispatches on the marker class.  Classes with a dedicated
#  drawing (KASP map, CAPS/ARMS gel, HRM melt) reuse it; the rest share one
#  schematic map engine (`plot_marker_map()`).  Every non-KASP panel carries the
#  class's own gate table on the right, drawn with the same tumbler as
#  `plot_checks()`, so a report never shows another chemistry's checks.
# =============================================================================

#' Plot a marker object
#'
#' Every marker class knows how to draw itself: the KASP assay keeps its local
#' map + acceptance-check panel, the gel chemistries keep their gel, and the
#' remaining classes draw a schematic of their own architecture with the gate
#' table valid for that class on the right.
#'
#' @param x a marker object.
#' @param variants optional variant table (KASP context).
#' @param genome optional genome, used for the KASP ideogram.
#' @param features optional feature table (gene model).
#' @param ... reserved for class-specific options.
#' @return a `ggplot` or `patchwork` object.
#' @name plot-marker
NULL

empty_variants <- function() {
  tibble::tibble(
    chrom = character(), pos = integer(), ref = character(),
    alt = character(), consequence = character()
  )
}

empty_plot <- function() {
  ggplot2::ggplot() +
    ggplot2::theme_void()
}

first_num <- function(...) {
  x <- suppressWarnings(as.numeric(unlist(list(...), use.names = FALSE)))
  x <- x[is.finite(x)]
  if (length(x)) x[1L] else NA_real_
}

# The SSR motif is data, not decoration: read the stored motif when the designer
# recorded it, otherwise recover it from the genome the plot already has.
ssr_motif <- function(rp, genome = NULL, chrom = NULL) {
  motif <- if ("motif" %in% names(rp) && !is.na(rp$motif[1])) rp$motif[1] else NULL
  if (is.null(motif) && !is.null(genome) && !is.null(chrom)) {
    motif <- tryCatch(
      as.character(genome_seq(genome, chrom, rp$start[1L], rp$end[1L])),
      error = function(e) NULL
    )
  }
  motif
}

ssr_unit_seq <- function(rp, motif, unit) {
  if ("unit_seq" %in% names(rp)) {
    return(rp$unit_seq[1L])
  }
  if (!is.null(motif)) substr(motif, 1L, unit) else NULL
}

ssr_formula <- function(unit_seq, copies) {
  if (is.null(unit_seq) || is.na(unit_seq)) {
    return(sprintf("%d copies", copies))
  }
  sprintf("(%s)%d", unit_seq, copies)
}

# One thermal guideline table for the non-KASP families.  A probe has its own
# window (usually 5-10 C above the primers); isothermal methods have no Tm gate,
# so their observed value is shown as neutral.
marker_tm_window <- function(scheme, primer = "") {
  family <- tryCatch(marker_scheme(scheme)$family, error = function(e) NULL)
  probe <- grepl("^probe", primer)
  if (isTRUE(family == "rpa")) {
    return(NULL)
  }
  if (isTRUE(family == "crispr") && identical(primer, "guide")) {
    return(NULL)
  }
  if (probe && family %in% c("taqman", "dpcr", "castpcr")) {
    return(c(65, 72))
  }
  switch(family,
    lamp = c(58, 65),
    taqman = c(58, 62),
    dpcr = c(58, 62),
    castpcr = c(55, 65),
    sanger = c(50, 65),
    snapshot = c(50, 65),
    crispr = c(55, 65),
    c(55, 65)
  )
}

primer_label_short <- function(nm) {
  nm <- sub("^outer_forward$", "oF", nm)
  nm <- sub("^outer_reverse$", "oR", nm)
  nm <- sub("^AS_ref$", "ASr", nm)
  nm <- sub("^AS_alt$", "ASa", nm)
  nm <- sub("^PCR_F$", "PCR F", nm)
  nm <- sub("^PCR_R$", "PCR R", nm)
  nm <- sub("^read_F$", "read F", nm)
  nm <- sub("^read_R$", "read R", nm)
  nm <- sub("^forward$", "F", nm)
  nm <- sub("^reverse$", "R", nm)
  nm <- sub("^probe_", "P", nm)
  nm <- sub("^probe$", "P", nm)
  nm <- sub("^extension$", "ext", nm)
  nm
}

# Full criterion table for a non-KASP marker.  This is the evaluator behind
# `marker_design_checks()` and the decision graph for every non-KASP family.
marker_generic_checks <- function(marker, scheme, genome = NULL, limits = KASP_LIMITS) {
  primers <- tryCatch(marker_primers(marker), error = function(e) character(0))
  primers <- primers[!is.na(primers) & nzchar(primers)]
  dna <- primers[grepl("^[ACGT]+$", primers)]
  crit <- marker_criteria(scheme = scheme, limits = limits)
  checks <- crit[crit$role == "check", , drop = FALSE]
  if (!nrow(checks)) {
    return(tibble::tibble(
      stage = character(), criterion = character(), observed = numeric(),
      limit = character(), pass = logical(), status = marker_check_status(),
      detail = character(), level = character()
    ))
  }
  field <- function(x, nm) if (!is.null(x) && nm %in% names(x)) x[[nm]] else NULL
  amp_span <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    s <- field(x, "start")
    e <- field(x, "end")
    if (is.null(s) || is.null(e)) NULL else e - s + 1L
  }
  a1 <- if (length(marker$assays)) marker$assays[[1L]] else NULL
  if (is.null(a1) && !is.null(marker$base)) a1 <- marker$base
  amp <- first_num(
    field(marker$amplicon, "product"), field(marker$primers, "product"),
    field(marker$amplicon, "ref_size"), field(marker$amplicon, "alt_size"),
    amp_span(marker$amplicon),
    field(a1$amplicon, "product"), field(a1$primers, "product"),
    field(a1$assay, "amplicon"), amp_span(a1$amplicon)
  )
  win <- marker_amplicon_window(marker_scheme(scheme)$family, limits)
  thermo <- if (length(dna) && has_oligoarrayaux()) {
    tryCatch(primer_thermo(dna), error = function(e) NULL)
  } else {
    NULL
  }
  dimers <- if (length(dna) > 1L && has_oligoarrayaux()) {
    tryCatch(dimer_thermo(dna), error = function(e) NULL)
  } else {
    NULL
  }
  hits <- if (!is.null(genome)) {
    tryCatch(primer_hits(dna, genome), error = function(e) NULL)
  } else if (!is.null(marker$ideogram$hits)) {
    marker$ideogram$hits
  } else {
    NULL
  }
  observed <- rep(NA_real_, nrow(checks))
  pass <- rep(NA, nrow(checks))
  family <- marker_scheme(scheme)$family
  native <- identical(class(marker)[1L], paste0("marker_", scheme))
  for (i in seq_len(nrow(checks))) {
    lab <- checks$label[i]
    # family checks only make sense on the matching marker class; a KASP mock
    # graded as e.g. `caps` must get NA, not a fabricated pass/fail.
    if (lab %in% (MARKER_FAMILY_CHECKS[[family]] %||% character(0)) && !native) {
      next
    }
    if (identical(lab, "amplicon length (bp)")) {
      observed[i] <- amp
      pass[i] <- is.finite(amp) && dplyr::between(amp, win[1L], win[2L])
    } else if (identical(lab, "longest homopolymer (bp)") && length(dna)) {
      observed[i] <- max(vapply(dna, longest_run, numeric(1)))
      pass[i] <- observed[i] < limits$max_run
    } else if (identical(lab, "SSR (2-4 nt, >=3 copies)") && length(dna)) {
      observed[i] <- sum(vapply(dna, has_repeat, logical(1)))
      pass[i] <- observed[i] == 0L
    } else if (identical(lab, "inverted repeat (arms>=4, loop<=6)") && length(dna)) {
      observed[i] <- sum(vapply(dna, has_palindrome, logical(1)))
      pass[i] <- observed[i] == 0L
    } else if (!is.null(thermo) && identical(lab, "worst hairpin dG37")) {
      v <- thermo$hairpin_dg[is.finite(thermo$hairpin_dg)]
      if (length(v)) {
        observed[i] <- min(v)
        pass[i] <- observed[i] >= limits$hairpin_dg
      }
    } else if (!is.null(dimers) && identical(lab, "worst cross-dimer dG37")) {
      dg <- dimers$dg[dimers$a != dimers$b]
      dg <- dg[is.finite(dg)]
      if (length(dg)) {
        observed[i] <- min(dg)
        pass[i] <- observed[i] >= limits$duplex_dg
      }
    } else if (!is.null(dimers) && identical(lab, "worst homodimer dG37")) {
      dg <- dimers$dg[dimers$a == dimers$b]
      dg <- dg[is.finite(dg)]
      if (length(dg)) {
        observed[i] <- min(dg)
        pass[i] <- observed[i] >= limits$homodimer_dg
      }
    } else if (!is.null(hits) && identical(lab, "off-target exact primer hits (n)")) {
      observed[i] <- if (nrow(hits)) max(table(hits$primer)) else 0
      pass[i] <- observed[i] <= 1
    } else if (identical(lab, "restriction site defined")) {
      observed[i] <- as.numeric(!is.null(marker$site))
      pass[i] <- !is.null(marker$site)
    } else if (identical(lab, "restriction site separates alleles")) {
      if (!is.null(marker$fragments_ref) && !is.null(marker$fragments_alt)) {
        observed[i] <- as.numeric(!identical(marker$fragments_ref, marker$fragments_alt))
        pass[i] <- !identical(marker$fragments_ref, marker$fragments_alt)
      }
    } else if (identical(lab, "band gap (bp)")) {
      observed[i] <- first_num(marker$gap)
      pass[i] <- is.finite(observed[i]) && observed[i] >= KASP_CAPS_GAP_MIN
    } else if (identical(lab, "digest sums to amplicon")) {
      if (length(marker$fragments_ref) && length(marker$fragments_alt) && is.finite(amp)) {
        observed[i] <- max(
          abs(sum(marker$fragments_ref) - amp),
          abs(sum(marker$fragments_alt) - amp)
        )
        pass[i] <- observed[i] == 0
      }
    } else if (identical(lab, "restriction site conserved (population)")) {
      h <- marker$haplotypes
      site <- marker$site
      if (!is.null(h) && nrow(h) && !is.null(site) && !anyNA(c(site$start[1L], site$end[1L]))) {
        extra <- h[
          h$pos >= site$start[1L] & h$pos <= site$end[1L] &
            h$pos != marker$snp$pos[1L], ,
          drop = FALSE
        ]
        observed[i] <- nrow(extra)
        pass[i] <- observed[i] == 0
      }
    } else if (identical(lab, "FIP 3' ALT")) {
      fip <- marker$primers$FIP
      if (!is.null(fip) && nzchar(fip)) {
        observed[i] <- as.numeric(substr(fip, nchar(fip), nchar(fip)) == marker$snp$alt[1L])
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "inner arm Tm balance (C)")) {
      tm <- marker$tm
      if (!is.null(tm) && all(c("FIP", "BIP") %in% names(tm))) {
        observed[i] <- abs(tm[["FIP"]] - tm[["BIP"]])
        pass[i] <- is.finite(observed[i]) && observed[i] <= 5
      }
    } else if (identical(lab, "inner arms inside amplicon")) {
      arms <- marker$lamp_arms
      lo <- field(marker$amplicon, "start")
      hi <- field(marker$amplicon, "end")
      if (!is.null(arms) && !is.null(lo) && !is.null(hi)) {
        coords <- unlist(lapply(arms, function(a) c(a$arm5, a$arm3)), use.names = FALSE)
        observed[i] <- as.numeric(all(coords >= lo & coords <= hi))
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "variant at FIP 3' end")) {
      fip <- marker$primers$FIP
      if (!is.null(fip) && nzchar(fip)) {
        observed[i] <- as.numeric(substr(fip, nchar(fip), nchar(fip)) == marker$snp$alt[1L])
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "primer length (nt)")) {
      ln <- nchar(primers)
      if (length(ln)) {
        observed[i] <- max(ln)
        pass[i] <- all(ln >= min(RPA_PRIMER_LEN) & ln <= max(RPA_PRIMER_LEN))
      }
    } else if (identical(lab, "THF off variant")) {
      pr <- marker$probe
      if (!is.null(pr) && !is.null(pr$start) && !is.null(pr$thf_index)) {
        voff <- marker$snp$pos[1L] - pr$start[1L] + 1L
        observed[i] <- as.numeric(!identical(as.integer(pr$thf_index), as.integer(voff)))
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "RPA probe spans variant")) {
      pr <- marker$probe
      if (!is.null(pr) && !is.null(pr$start) && !is.null(pr$end)) {
        observed[i] <- as.numeric(pr$start[1L] <= marker$snp$pos[1L] && pr$end[1L] >= marker$snp$pos[1L])
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "variant under primer/probe")) {
      under <- field(marker$primers, "variant_under")
      if (!is.null(under) && length(under) && nzchar(under[1L])) {
        observed[i] <- 1
        pass[i] <- TRUE
      }
    } else if (identical(lab, "PAM valid")) {
      g <- field(marker$gates, "gate_crispr_pam")
      if (!is.null(g) && length(g)) {
        observed[i] <- as.numeric(isTRUE(g[[1L]]))
        pass[i] <- isTRUE(g[[1L]])
      }
    } else if (identical(lab, "variant in seed")) {
      g <- field(marker$gates, "gate_crispr_seed")
      if (!is.null(g) && length(g)) {
        observed[i] <- as.numeric(isTRUE(g[[1L]]))
        pass[i] <- isTRUE(g[[1L]])
      }
    } else if (identical(lab, "guide spans variant")) {
      if (!is.null(marker$guide_start) && !is.null(marker$guide_end)) {
        observed[i] <- as.numeric(
          marker$guide_start[1L] <= marker$snp$pos[1L] &&
            marker$guide_end[1L] >= marker$snp$pos[1L]
        )
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "guide/PAM conserved (population)")) {
      h <- marker$haplotypes
      if (!is.null(h) && nrow(h) && !is.null(marker$guide_start)) {
        lo <- marker$guide_start[1L]
        hi <- marker$guide_end[1L]
        if (!is.null(marker$pam) && !is.na(marker$pam[1L])) {
          lo <- min(lo, marker$guide_start[1L] - 4L)
          hi <- max(hi, marker$guide_start[1L] - 1L)
        }
        extra <- h[
          h$pos >= lo & h$pos <= hi & h$pos != marker$snp$pos[1L], ,
          drop = FALSE
        ]
        observed[i] <- nrow(extra)
        pass[i] <- observed[i] == 0
      }
    } else if (identical(lab, "off-target protospacers (n)")) {
      observed[i] <- first_num(marker$off_target)
      pass[i] <- is.finite(observed[i]) && observed[i] <= 0
    } else if (identical(lab, "probe spans variant")) {
      pr <- marker$probes
      if (!is.null(pr) && all(c("start", "end") %in% names(pr))) {
        observed[i] <- as.numeric(any(
          pr$start <= marker$snp$pos[1L] & pr$end >= marker$snp$pos[1L]
        ))
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "probe Tm - primer Tm (C)")) {
      pr <- marker$probes
      if (!is.null(pr) && "tm" %in% names(pr) && length(dna)) {
        ptm <- primer_tm(dna)
        observed[i] <- min(pr$tm, na.rm = TRUE) - max(ptm, na.rm = TRUE)
        pass[i] <- is.finite(observed[i]) && observed[i] >= 5
      }
    } else if (identical(lab, "ASP 3' discrimination")) {
      observed[i] <- first_num(marker$primers$as_discrimination)
      pass[i] <- is.finite(observed[i]) && observed[i] > 0
    } else if (identical(lab, "extension 3' at variant")) {
      ext <- marker$extension
      if (!is.null(ext) && !is.null(ext$start) && !is.null(ext$end)) {
        d <- min(abs(c(ext$start, ext$end) - marker$snp$pos[1L]))
        observed[i] <- d
        pass[i] <- d <= 1
      }
    } else if (identical(lab, "variant-read distance (bp)")) {
      lo <- field(marker$amplicon, "start")
      hi <- field(marker$amplicon, "end")
      if (!is.null(lo) && !is.null(hi)) {
        observed[i] <- min(marker$snp$pos[1L] - lo, hi - marker$snp$pos[1L])
        pass[i] <- observed[i] >= KASP_SANGER_LIMITS$min_dist
      }
    } else if (identical(lab, "read primers inside amplicon")) {
      sp <- marker$seq_primers
      lo <- field(marker$amplicon, "start")
      hi <- field(marker$amplicon, "end")
      sp_ok <- !is.null(sp) && !is.null(sp$forward_start) && !is.null(sp$reverse_end)
      if (!is.null(lo) && !is.null(hi) && sp_ok) {
        observed[i] <- as.numeric(
          sp$forward_start >= lo && sp$reverse_end <= hi
        )
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "read primer pair defined")) {
      sp <- marker$seq_primers
      if (!is.null(sp)) {
        observed[i] <- as.numeric(!is.null(sp$forward) && !is.null(sp$reverse))
        pass[i] <- observed[i] == 1
      }
    } else if (identical(lab, "repeat unit (nt)")) {
      rp <- marker[["repeat"]]
      if (!is.null(rp) && "unit" %in% names(rp)) {
        observed[i] <- rp$unit[1L]
        pass[i] <- dplyr::between(observed[i], 2, 4)
      }
    } else if (identical(lab, "repeat copies")) {
      rp <- marker[["repeat"]]
      if (!is.null(rp) && "copies" %in% names(rp)) {
        observed[i] <- rp$copies[1L]
        pass[i] <- observed[i] >= 5
      }
    } else if (identical(lab, "REF/ALT size difference (bp)")) {
      rs <- field(marker$amplicon, "ref_size")
      as <- field(marker$amplicon, "alt_size")
      rp <- marker[["repeat"]]
      if (!is.null(rs) && !is.null(as)) {
        observed[i] <- abs(rs[1L] - as[1L])
      } else if (!is.null(rp) && "unit" %in% names(rp)) {
        observed[i] <- rp$unit[1L]
      }
      if (is.finite(observed[i])) pass[i] <- observed[i] >= 3
    } else if (identical(lab, "fragment size (bp)")) {
      observed[i] <- amp
      pass[i] <- is.finite(amp) && dplyr::between(amp, win[1L], win[2L])
    }
  }
  # Family checks belong only to their own marker class; elsewhere they are not
  # applicable, so they neither pass nor doom.  A check with no observed value
  # (missing data or an engine that was not run) is unknown, never a pass.
  status <- rep("unknown", nrow(checks))
  fam_labs <- MARKER_FAMILY_CHECKS[[family]] %||% character(0)
  status[checks$label %in% fam_labs & !native] <- "not_applicable"
  evaluated <- status != "not_applicable" & is.finite(observed) & !is.na(pass)
  status[evaluated] <- ifelse(pass[evaluated], "pass", "fail")
  status <- marker_check_status(status)
  out <- tibble::tibble(
    stage = checks$stage, criterion = checks$label, observed = observed,
    limit = checks$limit, pass = marker_check_pass(status), status = status,
    detail = NA_character_, level = checks$level
  )
  # internal hand-off for the display adapter; not part of the public table
  attr(out, "thermo") <- thermo
  out
}

# The marker-agnostic acceptance checks (the KASP_CHECK_PRIMER class), computed
# from the marker's own primers/amplicon/genome, so the non-KASP tumbler shows
# real observed values + ok/broken colours exactly like plot_checks().
marker_primer_checks <- function(marker, scheme, genome = NULL, limits = KASP_LIMITS) {
  primers <- tryCatch(marker_primers(marker), error = function(e) character(0))
  primers <- primers[!is.na(primers) & nzchar(primers)]
  primers <- primers[grepl("^[ACGT]+$", primers)]
  if (!length(primers)) {
    return(tibble::tibble())
  }
  generic <- marker_generic_checks(marker, scheme, genome, limits)
  thermo <- attr(generic, "thermo")
  base <- generic[generic$criterion %in% KASP_CHECK_PRIMER, , drop = FALSE]
  keep <- !is.na(base$status)
  out <- if (any(keep)) {
    name <- check_short(base$criterion[keep])
    tibble::tibble(
      level = base$level[keep],
      name = name,
      status = dplyr::case_when(
        as.character(base$status[keep]) == "pass" ~ "ok",
        as.character(base$status[keep]) == "fail" ~ "broken",
        as.character(base$status[keep]) == "unknown" ~ "unknown",
        TRUE ~ "not_applicable"
      ),
      value = format_check_value(name, base$observed[keep])
    )
  } else {
    tibble::tibble(level = character(), name = character(), status = character(), value = character())
  }
  if (!is.null(thermo) && nrow(thermo)) {
    nm <- names(primers)
    if (is.null(nm)) nm <- rep("primer", length(primers))
    tm_rows <- purrr::map2_dfr(
      nm, thermo$tm,
      function(primer, tm) {
        win <- marker_tm_window(scheme, primer)
        status <- if (is.null(win)) {
          "not_applicable"
        } else if (!is.finite(tm)) {
          "unknown"
        } else if (dplyr::between(tm, win[1L], win[2L])) {
          "ok"
        } else {
          "broken"
        }
        tibble::tibble(
          level = "feature", name = paste0("Tm ", primer_label_short(primer)),
          status = status, value = format_check_value("Tm", tm)
        )
      }
    )
    out <- dplyr::bind_rows(out, tm_rows)
  }
  dplyr::bind_rows(out, marker_population_checks(marker))
}

plot_marker_checks <- function(marker, scheme, genome = NULL) {
  if (inherits(marker, "marker_assay")) {
    return(plot_checks(marker))
  }
  # match the KASP column: acceptance checks only, each with an observed value
  checks <- marker_primer_checks(marker, scheme, genome)
  if (!nrow(checks)) {
    return(NULL)
  }
  # The tumbler palette has ok / broken / neutral only: an unevaluated or
  # inapplicable check is drawn neutral (never green).
  display <- c(ok = "ok", broken = "broken", unknown = "neutral", not_applicable = "neutral")
  checks$status <- unname(display[checks$status])
  checks$status[is.na(checks$status)] <- "neutral"
  checks$level <- factor(checks$level, levels = c("site", "feature", "alignment", "specificity"))
  tumbler_plot(checks)
}

# Genome-wide off-site panel for any marker class, reusing the KASP ideogram
# brick (bars + exact primer hits + target cursor).
marker_ideogram <- function(marker, genome) {
  if (is.null(genome)) {
    return(NULL)
  }
  primers <- tryCatch(marker_primers(marker), error = function(e) character(0))
  primers <- primers[!is.na(primers) & nzchar(primers)]
  primers <- primers[grepl("^[ACGT]+$", primers)]
  if (!length(primers)) {
    return(NULL)
  }
  chroms <- genome_chroms(genome)
  bars <- tibble::tibble(
    chrom = chroms,
    length = vapply(chroms, function(ch) genome_length(genome, ch), numeric(1))
  )
  hits <- tryCatch(primer_hits(primers, genome), error = function(e) NULL)
  if (is.null(hits)) {
    return(NULL)
  }
  if (!nrow(hits)) {
    hits <- tibble::tibble(primer = character(), chrom = character(), pos = integer())
  }
  roles <- c(AS_ref = "AS-ref", AS_alt = "AS-alt", COM = "COM")
  signal <- unname(roles[hits$primer])
  signal[is.na(signal)] <- hits$primer[is.na(signal)]
  hits$signal <- signal
  # paired markers also get the in-silico PCR products, as in the KASP ideogram
  fwd <- if ("forward" %in% names(primers)) primers[["forward"]] else NULL
  rev <- if ("reverse" %in% names(primers)) primers[["reverse"]] else NULL
  if (is.null(fwd) && "PCR_F" %in% names(primers)) fwd <- primers[["PCR_F"]]
  if (is.null(rev) && "PCR_R" %in% names(primers)) rev <- primers[["PCR_R"]]
  if (is.null(fwd) && "outer F" %in% names(primers)) fwd <- primers[["outer F"]]
  if (is.null(rev) && "outer R" %in% names(primers)) rev <- primers[["outer R"]]
  prods <- tibble::tibble()
  if (!is.null(fwd) && !is.null(rev)) {
    pr <- tryCatch(in_silico_pcr(fwd, rev, genome), error = function(e) NULL)
    if (!is.null(pr) && nrow(pr)) {
      prods <- tibble::tibble(start = pr$start, end = pr$end, signal = "PCR")
    }
  }
  # composite markers (duplex) may not carry a top-level snp; use the first locus
  snp <- marker$snp
  if (is.null(snp) && length(marker$assays)) {
    # duplex: every locus is a target on the ideogram, not just the first
    snp <- dplyr::bind_rows(lapply(marker$assays, function(a) a$snp))
  }
  fake <- list(
    ideogram = list(bars = bars, hits = hits, products = prods),
    snp = snp
  )
  assay_ideogram(fake, snp)
}

# Fragment analysis panel: what the capillary/gel actually reads.
plot_fragment_readout <- function(marker, scheme, genome = NULL) {
  d <- NULL
  if (identical(scheme, "ssr")) {
    amp <- marker$amplicon
    rp <- marker[["repeat"]]
    if (!is.null(amp) && !is.null(rp)) {
      ref <- amp$end[1L] - amp$start[1L] + 1L
      unit <- rp$unit[1L]
      motif <- ssr_motif(rp, genome, marker$snp$chrom[1L])
      unit_seq <- ssr_unit_seq(rp, motif, unit)
      d <- tibble::tibble(
        allele = c("REF", "ALT (+1 unit)"),
        copies = c(rp$copies[1L], rp$copies[1L] + 1L),
        size = c(ref, ref + unit),
        delta = c(0L, unit)
      )
      d$formula <- if (is.null(unit_seq)) d$allele else ssr_formula(unit_seq, d$copies)
    }
  } else if (identical(scheme, "indel_size")) {
    amp <- marker$amplicon
    if (!is.null(amp) && all(c("ref_size", "alt_size") %in% names(amp))) {
      d <- tibble::tibble(
        allele = c("REF", "ALT"),
        size = c(amp$ref_size[1L], amp$alt_size[1L]),
        delta = c(0L, amp$alt_size[1L] - amp$ref_size[1L])
      )
      d$formula <- d$allele
    }
  }
  if (is.null(d) || !nrow(d) || anyNA(d$size)) {
    return(NULL)
  }
  d$label <- paste0(d$formula, " · ", d$size, " bp")
  d$label[d$delta > 0L] <- paste0(d$label[d$delta > 0L], " (+", d$delta[d$delta > 0L], ")")
  d$allele <- factor(d$allele, levels = d$allele)
  fills <- c(REF = "#E7E5E4", "ALT (+1 unit)" = "#FECACA", ALT = "#FECACA")
  lo <- min(d$size) - max(1, diff(range(d$size)) * 0.08)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$size, y = .data$allele)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = lo, xend = .data$size, y = .data$allele, yend = .data$allele),
      colour = "grey70", linewidth = 0.45
    ) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$allele), shape = 21, size = 2.9, colour = "grey35") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), hjust = -0.18, size = 2.2) +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.42))) +
    ggplot2::labs(x = "fragment size (bp)", y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6.5)
    )
}

marker_composite <- function(primary, marker, scheme, title = NULL, genome = NULL,
                             offsite = TRUE) {
  checks <- plot_marker_checks(marker, scheme, genome)
  off <- if (offsite) marker_ideogram(marker, genome) else NULL
  # the off-site ideogram goes BELOW the map, not beside it: the map keeps its
  # full width for sequence letters
  if (!is.null(off)) {
    primary <- patchwork::wrap_plots(primary, off, ncol = 1, heights = c(1, 0.22))
  }
  p <- if (is.null(checks)) primary else patchwork::wrap_plots(primary, checks, widths = c(1, 0.2))
  if (!is.null(title)) {
    p <- p + patchwork::plot_annotation(title = title)
  }
  p
}

# Every `plot.marker_*` method should be this call.  `above`/`below` are the
# only scheme-specific panels (a gel, a fragment readout); the map, ideogram and
# check column are always the shared bricks.
marker_figure <- function(marker, scheme, variants = NULL, genome = NULL,
                          features = NULL, above = NULL, below = NULL,
                          heights = NULL, title = NULL) {
  map <- plot_marker_map(marker, scheme, variants, genome, features)
  panels <- c(
    if (!is.null(above)) list(above),
    list(map),
    if (!is.null(below)) list(below)
  )
  primary <- if (length(panels) == 1L) {
    panels[[1L]]
  } else {
    patchwork::wrap_plots(panels, ncol = 1, heights = heights)
  }
  if (is.null(title)) {
    lab <- tryCatch(marker_scheme(scheme)$label, error = function(e) NULL)
    title <- lab %||% scheme
  }
  marker_composite(primary, marker, scheme, title, genome = genome)
}

spec_row <- function(feature, start, end, strand = "+", kind = "primer",
                     status = "ok", label = NA_character_) {
  bad <- is.null(start) || is.null(end) || length(start) != 1L || length(end) != 1L
  if (bad || is.na(start) || is.na(end)) {
    return(NULL)
  }
  tibble::tibble(
    feature = feature, start = as.integer(start), end = as.integer(end),
    strand = strand, kind = kind, status = status, label = label
  )
}

spec_bind <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(tibble::tibble(
      feature = character(), start = integer(), end = integer(),
      strand = character(), kind = character(), status = character(),
      label = character()
    ))
  }
  purrr::list_rbind(rows)
}

# The map contract.  Coordinates are plus-strand/genomic; primer lengths come
# from the sequences, so the arrows match the real oligos.  Every scheme owns
# its method, in its own file; the engine only knows this generic.
marker_map_spec <- function(marker, scheme = NULL, ...) {
  UseMethod("marker_map_spec")
}

#' @method marker_map_spec default
#' @export
#' @noRd
marker_map_spec.default <- function(marker, scheme = NULL, ...) {
  spec_bind(list())
}

# Common amplicon/paired-primer geometry, shared by the scheme methods.
map_amplicon <- function(marker) {
  amp <- marker$amplicon
  if (!inherits(amp, "data.frame")) {
    amp <- tibble::tibble(start = amp$start, end = amp$end)
  }
  c(start = as.integer(amp$start[1L]), end = as.integer(amp$end[1L]))
}

map_pair_spec <- function(forward, reverse, start, end) {
  fwd <- forward %||% ""
  rev <- reverse %||% ""
  if (!nchar(fwd) || !nchar(rev)) {
    return(list())
  }
  list(
    spec_row("forward", start, start + nchar(fwd) - 1L, "+", "primer"),
    spec_row("reverse", end - nchar(rev) + 1L, end, "-", "primer")
  )
}

map_pair_rows <- function(marker, start, end) {
  map_pair_spec(marker$primers$forward, marker$primers$reverse, start, end)
}

# One constructor for every scheme map.  A scheme supplies feature rows and the
# sequence of each row; the renderer supplies geometry, colour, letter size and
# labels.  No scheme method should rebuild the amplicon row itself.
map_spec_rows <- function(amplicon, features = list(), sequences = list()) {
  rows <- c(
    list(spec_row("amplicon", amplicon[1L], amplicon[2L], "+", "amplicon")),
    features
  )
  spec <- spec_bind(rows)
  spec$seq <- NA_character_
  for (nm in names(sequences)) {
    spec <- map_set_seq(spec, nm, sequences[[nm]])
  }
  spec
}

map_spec <- function(marker, features = list(), sequences = list()) {
  map_spec_rows(map_amplicon(marker), features, sequences)
}

feature_clamp_n <- function(kind) {
  # only an actual primer has a 3' GC clamp; probes, guides, sites and blockers
  # are features, not directional oligos
  ifelse(kind %in% c("primer", "seq"), 5L, 0L)
}

# The KASP "nose" readout for every primer row: Tm text just past the 3' end
# and a status dot one step further out.
primer_nose_spec <- function(spec, scheme, span = range(c(spec$start, spec$end))) {
  primers <- spec[
    spec$kind %in% c("primer", "seq") & !is.na(spec$seq), ,
    drop = FALSE
  ]
  if (!nrow(primers)) {
    return(NULL)
  }
  tm <- rep(NA_real_, nrow(primers))
  dna <- grepl("^[ACGT]+$", primers$seq)
  if (any(dna)) {
    tm[dna] <- suppressWarnings(primer_tm(primers$seq[dna]))
  }
  primers <- primers[is.finite(tm), , drop = FALSE]
  tm <- tm[is.finite(tm)]
  if (!nrow(primers)) {
    return(NULL)
  }
  forward <- primers$strand == "+"
  tip <- ifelse(forward, primers$end, primers$start)
  off <- max(1L, round(diff(span) * 0.008))
  status <- vapply(seq_len(nrow(primers)), function(i) {
    win <- marker_tm_window(scheme, as.character(primers$feature[i]))
    if (is.null(win)) "neutral" else if (dplyr::between(tm[i], win[1L], win[2L])) "ok" else "broken"
  }, character(1))
  tibble::tibble(
    x = tip,
    y = primers$y,
    tm = tm,
    label = paste0(round(tm), "\u00b0C"),
    label_x = ifelse(forward, tip + off, tip - off),
    dot_x = ifelse(forward, tip + 3L * off, tip - 3L * off),
    hjust = ifelse(forward, 0, 1),
    status = status,
    off = off
  )
}

map_probe_rows <- function(probes, prefix = "probe") {
  if (is.null(probes) || !nrow(probes) || !all(c("start", "end") %in% names(probes))) {
    return(list())
  }
  if ("allele" %in% names(probes)) {
    return(lapply(seq_len(nrow(probes)), function(i) {
      spec_row(
        paste(prefix, probes$allele[i]), probes$start[i], probes$end[i], "+", "probe",
        label = paste(prefix, probes$allele[i])
      )
    }))
  }
  list(spec_row(prefix, probes$start[1L], probes$end[1L], "+", "probe", label = prefix))
}

map_set_seq <- function(spec, feature, value) {
  i <- match(feature, spec$feature)
  if (!is.na(i) && length(value) && !is.na(value[1L])) {
    spec$seq[i] <- as.character(value[1L])
  }
  spec
}

#' @method marker_map_spec marker_assay
#' @export
#' @noRd
marker_map_spec.marker_assay <- function(marker, scheme = NULL, ...) {
  aa <- marker$assay[1L, ]
  map_spec_rows(
    amplicon = c(
      min(aa$as_start, aa$as_alt_start, aa$com_start),
      max(aa$as_end, aa$as_alt_end, aa$com_end)
    ),
    features = list(
      spec_row("AS-ref", aa$as_start[1L], aa$as_end[1L], aa$strand[1L], "primer"),
      spec_row("AS-alt", aa$as_alt_start[1L], aa$as_alt_end[1L], aa$strand[1L], "primer"),
      spec_row(
        "COM", aa$com_start[1L], aa$com_end[1L],
        if (aa$strand[1L] == "+") "-" else "+", "primer"
      )
    ),
    sequences = list(
      "AS-ref" = aa$as_ref, "AS-alt" = aa$as_alt, COM = aa$com
    )
  )
}
plot_marker_map <- function(marker, scheme, variants = NULL, genome = NULL, features = NULL) {
  spec <- marker_map_spec(marker, scheme)
  if (!nrow(spec)) {
    return(empty_plot())
  }
  snp_pos <- marker$snp$pos[1L]
  span <- range(c(spec$start, spec$end))
  letter_size <- letter_size_for_span(span)
  pad <- max(1L, round(diff(span) * 0.03))
  has_hap <- !is.null(marker$haplotypes) && nrow(marker$haplotypes) > 0L
  # A haplotype table is richer than a SNP strip: when it exists, replace the
  # single-position track with the haplotype rows.  SSR uses its fragment track.
  if (!identical(scheme, "ssr") && !has_hap) {
    spec <- dplyr::bind_rows(spec, spec_row("variants", snp_pos, snp_pos, "+", "variant"))
  }
  if (!is.null(features)) {
    spec <- dplyr::bind_rows(spec, spec_row("feature", span[1L], span[2L], "+", "feature"))
  }
  spec$feature <- factor(spec$feature, levels = rev(unique(spec$feature)))
  # one numeric y for every layer: mixing numeric rects with discrete text rows
  # was what put letters/labels outside their boxes
  spec$y <- as.integer(spec$feature)
  spec$forward <- spec$strand == "+"
  nose <- primer_nose_spec(spec, scheme, span)
  if (!is.null(nose)) {
    pad <- max(pad, 4L * nose$off[1L])
  }
  # Primers (and sequencing primers) are the only arrowheads.  Every other
  # marker element is the same neutral block: same fill, same outline, same
  # aligned letters; its type is carried only by the superscript tag.
  arrows <- spec[spec$kind %in% c("primer", "seq"), ]
  blocks <- spec[spec$kind %in% c("amplicon", "probe", "site", "tract", "guide", "blocker"), ]
  outline <- function(x) {
    ifelse(x$status == "broken", MARKER_STATUS[["broken"]], "grey45")
  }
  kinds <- unique(spec$kind)
  # a neutral block must still be visible on white, otherwise the letters look
  # unboxed; outline + a slightly darker fill carry the box
  fills <- stats::setNames(rep("#E7E5E4", length(kinds)), kinds)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = tibble::tibble(x = span[1L], xend = span[2L]),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = 1, yend = 1),
      colour = "grey75", linewidth = 0.2
    )
  if (nrow(arrows)) {
    p <- p + gggenes::geom_gene_arrow(
      data = arrows,
      ggplot2::aes(
        xmin = .data$start, xmax = .data$end, y = .data$y,
        forward = .data$forward, fill = .data$kind
      ),
      colour = outline(arrows), linewidth = 0.35
    )
  }
  if (nrow(blocks)) {
    p <- p + ggplot2::geom_rect(
      data = blocks,
      ggplot2::aes(
        xmin = .data$start, xmax = .data$end,
        ymin = .data$y - 0.25,
        ymax = .data$y + 0.25,
        fill = .data$kind
      ),
      colour = outline(blocks), linewidth = 0.4
    )
  }
  # the same letters layer as plot_primers(): grey sequence, bold GC clamp,
  # red SSR / inverted-repeat bases
  lrows <- spec[!is.na(spec$seq), ]
  if (nrow(lrows)) {
    letters <- purrr::list_rbind(lapply(seq_len(nrow(lrows)), function(i) {
      d <- oligo_letters(
        lrows$seq[i], lrows$start[i], lrows$end[i], lrows$strand[i] == "+",
        as.character(lrows$feature[i]),
        clamp_n = feature_clamp_n(lrows$kind[i])
      )
      d$size <- letter_size_for_feature(
        lrows$end[i] - lrows$start[i] + 1L, letter_size
      )
      d$y <- lrows$y[i]
      d
    }))
    p <- p + letters_layer(letters, size = letter_size) +
      ggplot2::scale_colour_identity()
  }
  # the amplicon sequence: the same letters layer and the same size as the
  # primers, so a long product shrinks the whole stack, not just itself
  if (!is.null(genome) && "amplicon" %in% spec$feature) {
    span_all <- range(c(spec$start, spec$end))
    amp_seq <- tryCatch(
      as.character(genome_seq(genome, marker$snp$chrom[1L], span_all[1L], span_all[2L])),
      error = function(e) NULL
    )
    if (!is.null(amp_seq) && length(amp_seq) == 1L && nchar(amp_seq) > 0L) {
      amp_letters <- oligo_letters(
        amp_seq, span_all[1L], span_all[2L], TRUE, "amplicon",
        clamp_n = 0L
      )
      amp_letters$y <- spec$y[match("amplicon", as.character(spec$feature))]
      p <- p + letters_layer(amp_letters, size = letter_size)
    }
  }
  # KASP-style nose readout on every primer/sequencing-primer arrow
  if (!is.null(nose) && nrow(nose)) {
    p <- p + ggplot2::geom_text(
      data = nose,
      ggplot2::aes(x = .data$label_x, y = .data$y, label = .data$label, hjust = .data$hjust),
      inherit.aes = FALSE, size = 2.1, colour = "black"
    )
    for (st in unique(nose$status)) {
      fill <- switch(st,
        ok = MARKER_STATUS[["ok"]],
        broken = MARKER_STATUS[["broken"]],
        "#A8A29E"
      )
      p <- p + ggplot2::geom_point(
        data = nose[nose$status == st, , drop = FALSE],
        ggplot2::aes(x = .data$dot_x, y = .data$y),
        inherit.aes = FALSE, shape = 21, size = 2.2,
        fill = fill, colour = "grey40"
      )
    }
  }
  # the target label: an SSR locus is a repeat, not a REF>ALT SNP
  snp <- marker$snp
  if (!identical(scheme, "ssr") && !has_hap && !is.null(snp) && nrow(snp)) {
    var_y <- spec$y[match("variants", as.character(spec$feature))]
    p <- p + ggplot2::geom_text(
      data = tibble::tibble(
        x = snp$pos[1L], y = var_y,
        label = paste0(snp$ref[1L], ">", snp$alt[1L])
      ),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE, vjust = -0.55, size = 2.2, colour = "firebrick"
    )
  }
  # One superscript hook for every non-primer element: `PAM`, `EcoRI TTS^AA`,
  # `(AC)5`, `probe REF`, `blocker`, ...  The label never sits on the letters.
  labs <- spec[!is.na(spec$label) & !spec$kind %in% c("primer", "seq"), , drop = FALSE]
  if (nrow(labs)) {
    # fixed offset: just above the block, still inside this row
    labs$y <- labs$y + 0.34
    p <- p + ggplot2::geom_text(
      data = labs,
      ggplot2::aes(x = (.data$start + .data$end) / 2, y = .data$y, label = .data$label),
      size = 2.0, colour = "#292524", vjust = 0, fontface = "bold"
    )
  }
  # feature track, same gggenes arrow as plot_primers()
  feat <- if (!is.null(features)) normalise_features(features, span) else NULL
  if (!is.null(feat)) {
    feat_y <- spec$y[match("feature", as.character(spec$feature))]
    p <- p + gggenes::geom_gene_arrow(
      data = feat,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = feat_y, forward = .data$forward),
      fill = "grey72", colour = "grey40"
    )
    lab_feat <- feat[!is.na(feat$name), ]
    if (nrow(lab_feat)) {
      p <- p + ggplot2::geom_text(
        data = lab_feat,
        ggplot2::aes(x = (.data$start + .data$end) / 2, y = feat_y, label = .data$name),
        size = 2.1, fontface = "bold"
      )
    }
  }
  p <- p +
    ggplot2::geom_vline(xintercept = snp_pos, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = spec$y[!duplicated(spec$y)],
      labels = as.character(spec$feature[!duplicated(spec$y)]),
      expand = ggplot2::expansion(add = 0.45)
    ) +
    ggplot2::coord_cartesian(xlim = c(span[1L] - pad, span[2L] + pad)) +
    ggplot2::labs(x = glue::glue("{marker$snp$chrom[1L]} position (bp)"), y = NULL) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(colour = "black", size = 7)
    )
  parts <- list(p)
  heights <- 2.5
  if (has_hap) {
    p_hap <- haplotype_layers(marker, snp, span)
    if (!is.null(p_hap)) {
      parts <- c(parts, list(p_hap))
      heights <- c(heights, 0.22 * length(unique(marker$haplotypes$haplotype)))
    }
  }
  if (!has_hap && !is.null(variants) && nrow(variants)) {
    vw <- variants_in_window(variants, marker$snp$chrom[1L], span)
    if (nrow(vw)) {
      if (!"consequence" %in% names(vw)) vw$consequence <- "variant"
      pv <- ggplot2::ggplot() +
        ggplot2::annotate("rect",
          xmin = min(vw$pos) - 0.5, xmax = max(vw$pos) + 0.5,
          ymin = 0.5, ymax = 1.5, fill = "grey95"
        ) +
        variant_points(vw, 1) +
        ggplot2::geom_point(
          data = marker$snp, ggplot2::aes(.data$pos, 1),
          shape = 21, fill = "firebrick", colour = "white", size = 2.4
        ) +
        ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL) +
        ggplot2::scale_y_continuous(limits = c(0.6, 1.4), breaks = NULL) +
        ggplot2::coord_cartesian(xlim = c(span[1L] - pad, span[2L] + pad)) +
        ggplot2::labs(x = NULL, y = "variants") +
        ggplot2::theme_minimal(base_size = 7) +
        ggplot2::theme(
          legend.position = "none", panel.grid = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()
        )
      parts <- c(parts, list(pv))
      heights <- c(heights, 0.25)
    }
  }
  if (length(parts) > 1L) {
    patchwork::wrap_plots(parts, ncol = 1, heights = heights)
  } else {
    p
  }
}

#' @rdname plot-marker
#' @export
plot.marker_assay <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  # plot_primers() already carries the ideogram AND the check column
  plot_primers(x, variants %||% empty_variants(), x$snp, genome = genome, features = features)
}

#' @rdname plot-marker
#' @export
plot.marker_doomed <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  plot_doomed_target(x, variants %||% empty_variants(), x$snp, features = features)
}

#' @rdname plot-marker
#' @export
plot.marker_caps <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(
    x, "caps", variants, genome, features,
    above = plot_caps_gel(x), heights = c(0.8, 1),
    title = "CAPS / dCAPS: gel, amplicon map, class gates"
  )
}

#' @rdname plot-marker
#' @export
plot.marker_arms <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(
    x, "arms", variants, genome, features,
    above = plot_arms_gel(x), heights = c(0.8, 1),
    title = "ARMS: gel, primer map, class gates"
  )
}

#' @rdname plot-marker
#' @export
plot.marker_hrm <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(x, "hrm", variants, genome, features, title = "HRM: amplicon map, class gates")
}

#' @rdname plot-marker
#' @export
plot.marker_taqman <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(x, "taqman", variants, genome, features, title = "TaqMan: probes on the amplicon")
}

#' @rdname plot-marker
#' @export
plot.marker_castpcr <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(
    x, "castpcr", variants, genome, features,
    title = "castPCR: mutant ASP, WT blocker, locus probe"
  )
}

#' @rdname plot-marker
#' @export
plot.marker_dpcr <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(x, "dpcr", variants, genome, features, title = "dPCR: hydrolysis probes on the amplicon")
}

#' @rdname plot-marker
#' @export
plot.marker_snapshot <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(
    x, "snapshot", variants, genome, features,
    title = "SNaPshot: extension primer at the variant"
  )
}

#' @rdname plot-marker
#' @export
plot.marker_lamp <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(
    x, "lamp", variants, genome, features,
    title = "LAMP: F3/B3 outer, FIP/BIP inner arms"
  )
}

#' @rdname plot-marker
#' @export
plot.marker_rpa <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(x, "rpa", variants, genome, features, title = "RPA: primers and exo probe")
}

#' @rdname plot-marker
#' @export
plot.marker_crispr <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(x, "crispr", variants, genome, features, title = "CRISPR-Dx: guide + PAM")
}

#' @rdname plot-marker
#' @export
plot.marker_sanger <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(x, "sanger", variants, genome, features, title = "Sanger: amplicon and read primers")
}

#' @rdname plot-marker
#' @export
plot.marker_ssr <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(
    x, "ssr", variants, genome, features,
    below = plot_fragment_readout(x, "ssr", genome), heights = c(1, 0.28),
    title = "SSR: repeat tract and fragment sizes"
  )
}

#' @rdname plot-marker
#' @export
plot.marker_indel_size <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  marker_figure(
    x, "indel_size", variants, genome, features,
    below = plot_fragment_readout(x, "indel_size", genome), heights = c(1, 0.28),
    title = "InDel size: REF/ALT products"
  )
}

#' @rdname plot-marker
#' @export
plot.marker_duplex <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  maps <- lapply(x$assays, function(a) {
    if (inherits(a, "marker_assay")) {
      plot_marker_map(a, "snp", variants, genome, features)
    } else {
      plot_doomed_target(a, variants %||% empty_variants(), a$snp)
    }
  })
  marker_composite(
    patchwork::wrap_plots(maps, ncol = 1), x, "duplex", "Duplex: two loci in one tube",
    genome = genome
  )
}

#' @rdname plot-marker
#' @export
plot.marker_lax <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  if (is.null(x$base)) {
    return(empty_plot())
  }
  if (inherits(x$base, "marker_assay")) {
    return(plot_primers(x$base, variants %||% empty_variants(), x$base$snp,
      genome = genome, features = features
    ))
  }
  plot_doomed_target(x$base, variants %||% empty_variants(), x$base$snp)
}
