# =============================================================================
#  marker_checks_generic.R  --  non-KASP criterion table for one marker
# =============================================================================

#' Full criterion table for a non-KASP marker.  This is the evaluator behind
#' `marker_design_checks()` and the decision graph for every non-KASP family.
#' @param engines optional `marker_engines` adapter list (see `engine_has()`);
#'   `NULL` reads the session default; it gates the generic thermo screen.
#' @keywords internal
#' @noRd
marker_generic_checks <- function(marker, scheme, genome = NULL, limits = KASP_LIMITS,
                                  engines = NULL) {
  primers <- tryCatch(marker_primers(marker), error = function(e) character(0))
  primers <- primers[!is.na(primers) & nzchar(primers)]
  dna <- primers[grepl("^[ACGT]+$", primers)]
  crit <- marker_criteria(scheme = scheme, limits = limits)
  checks <- crit[crit$role == "check", , drop = FALSE]
  if (!nrow(checks)) {
    return(marker_evaluation_table(
      stage = character(), criterion = character(), observed = numeric(),
      limit = character(), pass = logical(), level = character()
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
  prod <- first_num(
    field(marker$amplicons, "n_products"), field(a1$amplicons, "n_products")
  )
  win <- marker_amplicon_window(marker_scheme(scheme)$family, limits)
  thermo <- if (length(dna) && engine_has("hairpin_dg", engines)) {
    tryCatch(primer_thermo(dna), error = function(e) NULL)
  } else {
    NULL
  }
  dimers <- if (length(dna) > 1L && engine_has("duplex_dg", engines)) {
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
  # asked once, not once per check row: the labels are a property of the family
  fam_labs <- marker_family_check_labels(family)
  for (i in seq_len(nrow(checks))) {
    lab <- checks$label[i]
    # family checks only make sense on the matching marker class; a KASP mock
    # graded as e.g. `caps` must get NA, not a fabricated pass/fail.
    if (lab %in% fam_labs && !native) {
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
    } else if (identical(lab, "target is a single-base substitution")) {
      # §W: the design function refuses a multi-base target and returns NULL, so
      # such a target was never evaluated against the scheme's bar.  Report it
      # `unknown` (not a pass, and not a `fail` the design never reached).
      ref <- field(marker$snp, "ref")
      alt <- field(marker$snp, "alt")
      if (length(ref) && length(alt) && nzchar(ref[1L]) && nzchar(alt[1L])) {
        single <- nchar(ref[1L]) == 1L && nchar(alt[1L]) == 1L
        observed[i] <- as.numeric(single)
        pass[i] <- if (single) TRUE else NA
      }
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
        pos <- marker$snp$pos[1L]
        fwd_end <- marker$primers$start[1L] + nchar(marker$primers$forward[1L]) - 1L
        rev_start <- marker$primers$end[1L] - nchar(marker$primers$reverse[1L]) + 1L
        observed[i] <- as.numeric(rpa_variant_gate(
          under[1L], pos, fwd_end, rev_start, marker$probe
        ))
        pass[i] <- observed[i] == 1
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
    } else if (identical(lab, "in-silico PCR products (n)")) {
      if (is.finite(prod)) {
        observed[i] <- prod
        pass[i] <- prod <= 1
      }
    }
  }
  # Family checks belong only to their own marker class; elsewhere they are not
  # applicable, so they neither pass nor doom.  A check with no observed value
  # (missing data or an engine that was not run) is unknown, never a pass.
  # A class whose design never paired its primers onto a genome has no product
  # count to grade, so that row is not_applicable - the same reading as the
  # duplex pairing in `R/marker_duplex.R` - not an unknown that poisons the
  # verdict.
  fam_labs <- marker_family_check_labels(family)
  out <- marker_evaluation_table(
    stage = checks$stage, criterion = checks$label, observed = observed,
    limit = checks$limit, pass = pass, level = checks$level,
    not_applicable = (checks$label %in% fam_labs & !native) |
      (checks$label == "in-silico PCR products (n)" & !is.finite(prod))
  )
  # internal hand-off for the display adapter; not part of the public table
  attr(out, "thermo") <- thermo
  # ZW2.2: name the engines that produced these numbers, so a drifted model is
  # visible next to the verdict it moved
  marker_with_stamp(out)
}

#' Checks for every designed marker in a solved object
#' @keywords internal
#' @noRd
marker_solved_checks <- function(solved, limits = KASP_LIMITS, scheme = NULL, engines = NULL) {
  r <- solved$report
  idx <- which(r$ok & !vapply(r$marker, is.null, logical(1)))
  if (!length(idx)) {
    return(tibble::tibble(
      label = character(), stage = character(), criterion = character(),
      observed = numeric(), limit = numeric(), pass = logical(),
      status = factor(), detail = character(), level = character()
    ))
  }
  parts <- lapply(idx, function(i) {
    chk <- marker_design_checks(
      r$marker[[i]],
      limits = limits,
      scheme = scheme %||% r$scheme[[i]], engines = engines
    )
    chk$label <- r$label[[i]]
    chk[, c("label", setdiff(names(chk), "label"))]
  })
  out <- dplyr::bind_rows(parts)
  # bind_rows() drops the engine stamp; all parts share one, so re-stamp
  marker_with_stamp(out, from = parts[[1L]])
}
