# =============================================================================
#  marker_crate.R  --  GWAS loci + a locus dossier -> one marker crate
# =============================================================================
#  The panel orchestrator starts at a target table.  A programme starts one
#  step earlier: a lead variant per GWAS locus, the dossier that supported the
#  locus, and a population reference that says which linked variant is worth
#  assaying.  This file selects one candidate per locus and hands the survivors
#  to marker_design_panel(); it does not design an assay, group a pot or write
#  the design artefact itself.  The crate around that artefact records the
#  GWAS run, the dossier, the rejected candidates and the specificity scan.
# =============================================================================

#' Build one reproducible marker crate from GWAS loci and a locus dossier
#'
#' `marker_locus_crate()` joins a GWAS locus table to its locus dossier, selects
#' one candidate variant per locus from LD, haplotype and MAF context, and hands
#' the selection to [marker_design_panel()].  The panel designs the assays,
#' groups the multiplex and writes its artefact; the crate adds the GWAS run,
#' dossier evidence, rejected candidates, specificity scan and an RO-Crate
#' descriptor around it.
#'
#' @param gwas a path or data frame with one lead row per locus: `locus`,
#'   `chrom`, `pos`, `ref`, `alt`, `population`, `genome_build`, `gwas_run`;
#'   optional `p`, `beta`, `se`, `maf`, `lead`.
#' @param dossier a path or data frame with one or more rows per GWAS locus:
#'   `locus`, `gene`, `source`, `source_version`, `retrieved_at`, `confidence`;
#'   optional `pathway`, `citation`, `doi`.
#' @param genome the reference sequence: a `Biostrings::DNAStringSet` or a
#'   `marker_genome` (see [marker_genome()]).
#' @param variants optional candidate variants keyed by `locus` (`locus`,
#'   `chrom`, `pos`, `ref`, `alt`; optional `maf`, `ld_r2`, `population`).  The
#'   GWAS lead is always a candidate even when the table omits it.
#' @param haplotypes optional phased haplotypes (see [marker_design_panel()]);
#'   they are the dosage source for LD and the population model for MAF.
#' @param chemistries a `KASP_TECH_PRESETS` name or a scheme id vector; the
#'   default `"kasp"` expands to the preset's schemes.
#' @param index a `genome_index` or a named list of indexes (a pangenome);
#'   `NULL` builds `genome_index(genome, strategy = "memory")` for a
#'   fixture-sized `DNAStringSet`.
#' @param maf_min,ld_min optional inclusive lower bounds applied to measured
#'   MAF and LD r2.  A bound on an unmeasured value aborts rather than passing.
#' @param dir the crate directory to write.
#' @param created optional crate timestamp; see [marker_artefact()].
#' @param ... passed to [marker_design_panel()] and then to [marker_solve()].
#' @return the crate directory, a character vector of length one.  The crate is
#'   an RO-Crate: `marker-crate.json` is the crate record and `artefact.json`
#'   is the panel's design record.
#' @examples
#' \dontrun{
#' crate <- marker_locus_crate(
#'   "gwas-loci.tsv", "locus-dossier.tsv",
#'   genome = marker_genome("soybean.fa", build = "Wm82.a4.v1"),
#'   variants = "candidate-variants.tsv",
#'   haplotypes = "population-haplotypes.csv",
#'   chemistries = "kasp"
#' )
#' }
#' @family plan
#' @export
marker_locus_crate <- function(gwas, dossier, genome, variants = NULL,
                               haplotypes = NULL, chemistries = "kasp",
                               index = NULL, maf_min = NULL, ld_min = NULL,
                               dir = tempfile("marker-crate"), created = NULL,
                               ...) {
  if (!inherits(genome, "DNAStringSet") && !inherits(genome, "marker_genome")) {
    abort_msg(
      "`genome` must be a DNAStringSet or a marker_genome (see marker_genome())."
    )
  }
  crate_check_bound(maf_min, "maf_min")
  crate_check_bound(ld_min, "ld_min")
  created <- artefact_created(created)

  gwas_in <- crate_input(gwas, "gwas")
  dossier_in <- crate_input(dossier, "dossier")
  variants_in <- if (is.null(variants)) NULL else crate_input(variants, "variants")
  gwas_tbl <- crate_gwas(gwas_in$data)
  dossier_tbl <- crate_dossier(dossier_in$data, gwas_tbl)
  build <- unique(gwas_tbl$genome_build)
  haplotypes <- panel_haplotypes(haplotypes)
  candidates <- crate_candidates(variants_in, gwas_tbl, haplotypes, maf_min, ld_min)

  selected <- candidates$selected
  selection <- candidates$table

  artefact_dir <- fs::path(dir, "artefact")
  fs::dir_create(artefact_dir)
  readr::write_csv(gwas_tbl, fs::path(artefact_dir, "gwas-loci.csv"))
  readr::write_csv(dossier_tbl, fs::path(artefact_dir, "dossier.csv"))
  readr::write_csv(selection, fs::path(artefact_dir, "selection.csv"))

  schemes <- crate_chemistries(chemistries)
  panel <- marker_design_panel(
    genome,
    variants = candidates$variants,
    haplotypes = haplotypes,
    chemistries = schemes,
    targets = selected,
    created = created,
    dir = artefact_dir,
    ...
  )

  provenance <- crate_provenance(
    gwas_in, gwas_tbl, dossier_in, dossier_tbl,
    candidates$input
  )
  panel$solved$plan$provenance <- provenance
  artefact <- marker_artefact(
    panel$solved,
    dir = artefact_dir,
    created = created,
    scenario = "J3"
  )

  index <- index %||% genome_index(genome, strategy = "memory", build = build)
  crate_check_index(index, build)
  specificity <- crate_specificity(panel, index)
  readr::write_csv(specificity, fs::path(artefact_dir, "specificity.csv"))

  record <- crate_record(
    artefact_dir, artefact, gwas_in, gwas_tbl, dossier_in, dossier_tbl,
    selection, selected, panel, specificity, build, index, maf_min, ld_min,
    created
  )
  writeLines(
    jsonlite::toJSON(
      record,
      auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA
    ),
    fs::path(artefact_dir, "marker-crate.json")
  )

  marker_export_rocrate(
    artefact,
    dir = dir,
    extra = c(
      "gwas-loci.csv", "dossier.csv", "selection.csv", "specificity.csv",
      "marker-crate.json"
    )
  )
  as.character(dir)
}

# ---- input admission --------------------------------------------------------

#' Read a crate input: a delimited path or a data frame
#' @keywords internal
#' @noRd
crate_input <- function(x, arg) {
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    if (!fs::file_exists(x)) {
      abort_msg("No {.arg {arg}} file at {.path {x}}.")
    }
    delim <- if (grepl("\\.csv$", x, ignore.case = TRUE)) "," else "\t"
    data <- readr::read_delim(
      x,
      delim = delim,
      show_col_types = FALSE,
      col_types = readr::cols(.default = readr::col_character())
    )
    return(list(
      data = data, path = as.character(fs::path_abs(x)),
      sha256 = hash_file(x), bytes = as.integer(file.size(x)),
      name = fs::path_file(x)
    ))
  }
  if (!is.data.frame(x)) {
    abort_msg("`{arg}` must be a file path or a data frame.")
  }
  data <- tibble::as_tibble(x)
  list(
    data = data, path = NA_character_, sha256 = hash_object(data),
    bytes = NA_integer_, name = arg
  )
}

#' A value that must be one finite number in [0, 1]
#' @keywords internal
#' @noRd
crate_check_bound <- function(x, arg) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || x > 1) {
    abort_msg("`{arg}` must be NULL or one finite number in [0, 1].")
  }
  invisible(NULL)
}

#' The GWAS locus table, admitted and typed
#' @keywords internal
#' @noRd
crate_gwas <- function(gwas) {
  need_cols(
    gwas,
    c("locus", "chrom", "pos", "ref", "alt", "population", "genome_build", "gwas_run"),
    arg = "gwas"
  )
  if (!nrow(gwas)) {
    abort_msg("`gwas` has no rows; a crate needs at least one locus.")
  }
  out <- tibble::tibble(
    locus = as.character(gwas$locus),
    chrom = as.character(gwas$chrom),
    pos = suppressWarnings(as.integer(gwas$pos)),
    ref = toupper(as.character(gwas$ref)),
    alt = toupper(as.character(gwas$alt)),
    population = as.character(gwas$population),
    genome_build = as.character(gwas$genome_build),
    gwas_run = as.character(gwas$gwas_run)
  )
  for (nm in c("p", "beta", "se", "maf")) {
    if (nm %in% names(gwas)) {
      out[[nm]] <- suppressWarnings(as.numeric(gwas[[nm]]))
    }
  }
  if ("lead" %in% names(gwas)) {
    lead <- toupper(as.character(gwas$lead)) %in% c("TRUE", "T", "1")
    if (anyNA(lead) || !all(lead)) {
      abort_msg(
        "When `gwas$lead` is present every row must be the lead of its locus (`TRUE`)."
      )
    }
    out$lead <- TRUE
  }
  crate_check_locus_rows(out)
  dup <- duplicated(out$locus)
  if (any(dup)) {
    abort_msg(gettextf(
      "The GWAS table repeats locus %s; one lead row per locus is the contract.",
      out$locus[[which(dup)[[1L]]]]
    ))
  }
  builds <- unique(out$genome_build)
  if (length(builds) != 1L) {
    abort_msg(gettextf(
      "One crate means one genome build; the GWAS table names %s.",
      paste(builds, collapse = ", ")
    ))
  }
  if ("maf" %in% names(out)) {
    bad <- !is.na(out$maf) & (!is.finite(out$maf) | out$maf < 0 | out$maf > 1)
    if (any(bad)) {
      abort_msg(gettextf("GWAS row %d has `maf` outside [0, 1].", which(bad)[[1L]]))
    }
  }
  out
}

#' Refuse a lead row that cannot name a variant
#' @keywords internal
#' @noRd
crate_check_locus_rows <- function(x) {
  bad <- is.na(x$locus) | !nzchar(x$locus)
  if (any(bad)) {
    abort_msg(gettextf("GWAS row %d has no locus id.", which(bad)[[1L]]))
  }
  bad <- is.na(x$chrom) | !nzchar(x$chrom) | is.na(x$pos) | x$pos < 1L
  if (any(bad)) {
    abort_msg(gettextf(
      "GWAS row %d needs a chromosome and a positive integer `pos`.", which(bad)[[1L]]
    ))
  }
  dna <- function(v) !is.na(v) & nzchar(v) & grepl("^[ACGTN]+$", v)
  bad <- !dna(x$ref) | !dna(x$alt)
  if (any(bad)) {
    abort_msg(gettextf(
      "GWAS row %d needs DNA alleles in `ref` and `alt` (ACGTN).", which(bad)[[1L]]
    ))
  }
  same <- x$ref == x$alt
  if (any(same)) {
    abort_msg(gettextf(
      "GWAS row %d has `ref` equal to `alt` (%s); that is not a variant.",
      which(same)[[1L]], x$ref[[which(same)[[1L]]]]
    ))
  }
  bad <- is.na(x$population) | !nzchar(x$population) |
    is.na(x$genome_build) | !nzchar(x$genome_build) |
    is.na(x$gwas_run) | !nzchar(x$gwas_run)
  if (any(bad)) {
    abort_msg(gettextf(
      "GWAS row %d needs `population`, `genome_build` and `gwas_run`.", which(bad)[[1L]]
    ))
  }
  invisible(x)
}

#' The locus dossier, admitted and joined to the GWAS loci
#' @keywords internal
#' @noRd
crate_dossier <- function(dossier, gwas) {
  need_cols(
    dossier,
    c("locus", "gene", "source", "source_version", "retrieved_at", "confidence"),
    arg = "dossier"
  )
  out <- tibble::tibble(
    locus = as.character(dossier$locus),
    gene = as.character(dossier$gene),
    source = as.character(dossier$source),
    source_version = as.character(dossier$source_version),
    retrieved_at = as.character(dossier$retrieved_at),
    confidence = suppressWarnings(as.numeric(dossier$confidence))
  )
  for (nm in c("pathway", "citation", "doi")) {
    if (nm %in% names(dossier)) {
      out[[nm]] <- as.character(dossier[[nm]])
    }
  }
  bad <- is.na(out$locus) | !nzchar(out$locus) | is.na(out$gene) | !nzchar(out$gene) |
    is.na(out$source) | !nzchar(out$source) |
    is.na(out$source_version) | !nzchar(out$source_version)
  if (any(bad)) {
    abort_msg(gettextf(
      "Dossier row %d needs a `locus`, `gene`, `source` and `source_version`.",
      which(bad)[[1L]]
    ))
  }
  date_ok <- !is.na(out$retrieved_at) &
    grepl("^\\d{4}-\\d{2}-\\d{2}$", out$retrieved_at) &
    !is.na(suppressWarnings(as.Date(out$retrieved_at, format = "%Y-%m-%d")))
  if (any(!date_ok)) {
    abort_msg(gettextf(
      "Dossier row %d needs an ISO `retrieved_at` date (YYYY-MM-DD).",
      which(!date_ok)[[1L]]
    ))
  }
  bad <- is.na(out$confidence) | !is.finite(out$confidence) |
    out$confidence < 0 | out$confidence > 1
  if (any(bad)) {
    abort_msg(gettextf(
      "Dossier row %d needs `confidence` in [0, 1].", which(bad)[[1L]]
    ))
  }
  missing <- setdiff(gwas$locus, unique(out$locus))
  if (length(missing)) {
    abort_msg(gettextf(
      "The dossier has no row for GWAS locus/loci: %s.",
      paste(missing, collapse = ", ")
    ))
  }
  unknown <- setdiff(unique(out$locus), gwas$locus)
  if (length(unknown)) {
    abort_msg(gettextf(
      "The dossier names locus/loci absent from the GWAS table: %s.",
      paste(unknown, collapse = ", ")
    ))
  }
  out
}

#' Candidate variants keyed by locus, admitted and checked
#' @keywords internal
#' @noRd
crate_candidates_input <- function(variants, gwas) {
  if (is.null(variants)) {
    return(NULL)
  }
  v <- variants$data
  need_cols(v, c("locus", "chrom", "pos", "ref", "alt"), arg = "variants")
  out <- tibble::tibble(
    locus = as.character(v$locus),
    chrom = as.character(v$chrom),
    pos = suppressWarnings(as.integer(v$pos)),
    ref = toupper(as.character(v$ref)),
    alt = toupper(as.character(v$alt)),
    population = if ("population" %in% names(v)) as.character(v$population) else NA_character_,
    maf = if ("maf" %in% names(v)) suppressWarnings(as.numeric(v$maf)) else NA_real_,
    ld_r2 = if ("ld_r2" %in% names(v)) suppressWarnings(as.numeric(v$ld_r2)) else NA_real_
  )
  unknown <- setdiff(unique(out$locus), gwas$locus)
  if (length(unknown)) {
    abort_msg(gettextf(
      "The candidate table names locus/loci absent from the GWAS table: %s.",
      paste(unknown, collapse = ", ")
    ))
  }
  crate_check_variant_rows(out)
  for (nm in c("maf", "ld_r2")) {
    bad <- !is.na(out[[nm]]) & (!is.finite(out[[nm]]) | out[[nm]] < 0 | out[[nm]] > 1)
    if (any(bad)) {
      abort_msg(gettextf("Candidate row %d has `%s` outside [0, 1].", which(bad)[[1L]], nm))
    }
  }
  out
}

#' Refuse a candidate row that cannot name a variant
#' @keywords internal
#' @noRd
crate_check_variant_rows <- function(x) {
  bad <- is.na(x$locus) | !nzchar(x$locus)
  if (any(bad)) {
    abort_msg(gettextf("Candidate row %d has no locus id.", which(bad)[[1L]]))
  }
  bad <- is.na(x$chrom) | !nzchar(x$chrom) | is.na(x$pos) | x$pos < 1L
  if (any(bad)) {
    abort_msg(gettextf(
      "Candidate row %d needs a chromosome and a positive integer `pos`.", which(bad)[[1L]]
    ))
  }
  dna <- function(v) !is.na(v) & nzchar(v) & grepl("^[ACGTN]+$", v)
  bad <- !dna(x$ref) | !dna(x$alt)
  if (any(bad)) {
    abort_msg(gettextf(
      "Candidate row %d needs DNA alleles in `ref` and `alt` (ACGTN).", which(bad)[[1L]]
    ))
  }
  same <- x$ref == x$alt
  if (any(same)) {
    abort_msg(gettextf(
      "Candidate row %d has `ref` equal to `alt` (%s); that is not a variant.",
      which(same)[[1L]], x$ref[[which(same)[[1L]]]]
    ))
  }
  invisible(x)
}

#' Select one candidate per locus and record every candidate's verdict
#' @keywords internal
#' @noRd
crate_candidates <- function(variants, gwas, haplotypes, maf_min, ld_min) {
  cand <- crate_candidates_input(variants, gwas)
  lead <- tibble::tibble(
    locus = gwas$locus, chrom = gwas$chrom, pos = gwas$pos,
    ref = gwas$ref, alt = gwas$alt, population = gwas$population,
    maf = NA_real_, ld_r2 = NA_real_, is_lead = TRUE,
    gwas_maf = if ("maf" %in% names(gwas)) gwas$maf else NA_real_,
    duplicate = FALSE
  )
  if (is.null(cand)) {
    rows <- lead
  } else {
    cand$is_lead <- FALSE
    cand$gwas_maf <- NA_real_
    cand$duplicate <- FALSE
    rows <- dplyr::bind_rows(lead, cand)
  }
  rows <- rows[order(rows$locus, -as.integer(rows$is_lead)), , drop = FALSE]
  key <- paste(rows$locus, rows$chrom, rows$pos, rows$ref, rows$alt, sep = "\r")
  rows$duplicate <- duplicated(key)

  pop <- gwas$population[match(rows$locus, gwas$locus)]
  rows$population <- ifelse(is.na(rows$population), pop, rows$population)
  rows$status <- ifelse(
    rows$duplicate, "duplicate",
    ifelse(rows$population != pop, "population", "")
  )

  measured <- crate_measure(rows, haplotypes)
  rows$ld_measured <- measured$ld_r2
  rows$maf_measured <- measured$maf

  ld_table <- rows$ld_r2
  rows$ld_r2 <- ifelse(
    rows$is_lead, 1,
    ifelse(!is.na(ld_table), ld_table, rows$ld_measured)
  )
  rows$ld_source <- ifelse(
    rows$is_lead, "lead",
    ifelse(!is.na(ld_table), "table",
      ifelse(!is.na(rows$ld_measured), "haplotypes", "missing")
    )
  )
  table_maf <- ifelse(!is.na(rows$maf), rows$maf, NA_real_)
  rows$maf <- ifelse(
    !is.na(table_maf), table_maf,
    ifelse(!is.na(rows$maf_measured), rows$maf_measured,
      ifelse(rows$is_lead, rows$gwas_maf, NA_real_)
    )
  )
  rows$maf_source <- ifelse(
    !is.na(table_maf), "table",
    ifelse(!is.na(rows$maf_measured), "haplotypes",
      ifelse(rows$is_lead & !is.na(rows$gwas_maf), "gwas", "missing")
    )
  )

  rows <- crate_apply_bounds(rows, maf_min, ld_min)
  rows <- crate_rank(rows, maf_min, ld_min)
  selection <- rows[, c(
    "locus", "chrom", "pos", "ref", "alt", "is_lead", "population",
    "ld_r2", "ld_source", "maf", "maf_source", "rank", "status"
  )]
  chosen <- rows[rows$status == "selected", , drop = FALSE]
  if (!nrow(chosen)) {
    abort_msg(
      "No candidate passed the crate selection; every locus would abort. ",
      i = "See selection.csv or relax `maf_min` / `ld_min`."
    )
  }
  targets <- tibble::tibble(
    label = chosen$locus, chrom = chosen$chrom, pos = chosen$pos,
    ref = chosen$ref, alt = chosen$alt
  )
  for (nm in c("p", "beta", "se")) {
    if (nm %in% names(gwas)) {
      targets[[nm]] <- gwas[[nm]][match(chosen$locus, gwas$locus)]
    }
  }
  targets$population <- chosen$population
  targets$genome_build <- gwas$genome_build[match(chosen$locus, gwas$locus)]
  targets$gwas_run <- gwas$gwas_run[match(chosen$locus, gwas$locus)]
  targets$maf <- chosen$maf
  targets$ld_r2 <- chosen$ld_r2
  list(
    table = selection, selected = targets,
    variants = cand, input = variants
  )
}

#' LD r2 against the lead and haplotype MAF for every candidate
#' @keywords internal
#' @noRd
crate_measure <- function(rows, haplotypes) {
  out <- tibble::tibble(ld_r2 = rep(NA_real_, nrow(rows)), maf = rep(NA_real_, nrow(rows)))
  if (is.null(haplotypes) || !nrow(haplotypes)) {
    return(out)
  }
  targets <- tibble::tibble(
    label = paste0("candidate-", seq_len(nrow(rows))),
    chrom = rows$chrom, pos = rows$pos, ref = rows$ref, alt = rows$alt,
    alt2 = NA_character_
  )
  dosage <- panel_dosage(targets, haplotypes)
  r2 <- suppressWarnings(stats::cor(dosage)^2)
  lead_at <- vapply(
    rows$locus,
    function(l) which(rows$is_lead & rows$locus == l)[[1L]],
    integer(1)
  )
  out$ld_r2 <- r2[cbind(seq_len(nrow(rows)), lead_at)]
  info <- marker_informativeness(targets, haplotypes = haplotypes)
  p <- info$alt_freq
  out$maf <- pmin(p, 1 - p)
  out
}

#' Apply the optional MAF and LD lower bounds
#' @keywords internal
#' @noRd
crate_apply_bounds <- function(rows, maf_min, ld_min) {
  if (!is.null(ld_min)) {
    bad <- rows$status == "" & is.na(rows$ld_r2)
    if (any(bad)) {
      i <- which(bad)[[1L]]
      abort_msg(gettextf(
        "Cannot apply `ld_min` to locus %s candidate %s:%d %s>%s: LD is unmeasured.",
        rows$locus[[i]], rows$chrom[[i]], rows$pos[[i]], rows$ref[[i]], rows$alt[[i]]
      ))
    }
    rows$status[rows$status == "" & rows$ld_r2 < ld_min] <- "ld"
  }
  if (!is.null(maf_min)) {
    bad <- rows$status == "" & is.na(rows$maf)
    if (any(bad)) {
      i <- which(bad)[[1L]]
      abort_msg(gettextf(
        "Cannot apply `maf_min` to locus %s candidate %s:%d %s>%s: MAF is unmeasured.",
        rows$locus[[i]], rows$chrom[[i]], rows$pos[[i]], rows$ref[[i]], rows$alt[[i]]
      ))
    }
    rows$status[rows$status == "" & rows$maf < maf_min] <- "maf"
  }
  rows
}

#' Rank the surviving candidates per locus and mark one selected
#' @keywords internal
#' @noRd
crate_rank <- function(rows, maf_min, ld_min) {
  rows$rank <- NA_integer_
  rows$status[rows$status == ""] <- "not_selected"
  for (locus in unique(rows$locus)) {
    at <- which(rows$locus == locus)
    ok <- at[rows$status[at] == "not_selected"]
    if (!length(ok)) {
      filters <- c(
        if (!is.null(ld_min)) "ld_min",
        if (!is.null(maf_min)) "maf_min"
      )
      if (!length(filters)) {
        filters <- "population/duplicate exclusions"
      }
      abort_msg(gettextf(
        "No candidate for locus %s passed the selection; the filters that emptied it were %s.",
        locus, paste(filters, collapse = " and ")
      ))
    }
    ord <- ok[order(
      -rows$ld_r2[ok], -rows$maf[ok], -as.integer(rows$is_lead[ok]),
      rows$pos[ok], rows$alt[ok],
      na.last = TRUE
    )]
    rows$rank[ord] <- seq_along(ord)
    rows$status[ord[[1L]]] <- "selected"
  }
  rows
}

# ---- design, specificity and the crate record -------------------------------

#' Expand a chemistry preset name or pass scheme ids through
#' @keywords internal
#' @noRd
crate_chemistries <- function(chemistries) {
  is_preset <- is.character(chemistries) &&
    length(chemistries) == 1L &&
    chemistries %in% names(KASP_TECH_PRESETS)
  if (is_preset) {
    return(KASP_TECH_PRESETS[[chemistries]])
  }
  chemistries
}

#' The panel's designed primers, scanned against the crate index
#' @keywords internal
#' @noRd
crate_specificity <- function(panel, index) {
  rows <- list()
  for (i in which(panel$targets$ok)) {
    locus <- panel$targets$target[[i]]
    primers <- marker_primers(panel$targets$assay[[i]])
    if (!length(primers)) {
      next
    }
    names(primers) <- paste0(locus, "|", names(primers))
    scan <- specificity_of(primers, index)
    if (!"genome" %in% names(scan)) {
      scan$genome <- NA_character_
    }
    scan$locus <- locus
    scan$target <- locus
    rows[[length(rows) + 1L]] <- scan[, c(
      "locus", "target", "primer", "seed", "hits", "engine", "build", "genome"
    )]
  }
  if (!length(rows)) {
    return(tibble::tibble(
      locus = character(), target = character(), primer = character(),
      seed = character(), hits = integer(), engine = character(),
      build = character(), genome = character()
    ))
  }
  dplyr::bind_rows(rows)
}

#' Refuse an index that cannot answer for the crate's build
#' @keywords internal
#' @noRd
crate_check_index <- function(index, build) {
  builds <- if (is.list(index) && !is.object(index)) {
    vapply(index, index_result_build, character(1))
  } else {
    index_result_build(index)
  }
  known <- builds[!is.na(builds)]
  if (length(known) && !build %in% known) {
    abort_msg(gettextf(
      "The specificity index answers build(s) %s, not the GWAS build %s.",
      paste(known, collapse = ", "), build
    ))
  }
  invisible(index)
}

#' Hashed provenance records for the plan and the crate
#' @keywords internal
#' @noRd
crate_provenance <- function(gwas_in, gwas, dossier_in, dossier, candidates) {
  runs <- unique(gwas$gwas_run)
  out <- list(
    list(
      role = "gwas", source = paste0("gwas:", paste(runs, collapse = ",")),
      name = paste0(nrow(gwas), " locus/loci"),
      path = gwas_in$path, sha256 = gwas_in$sha256, bytes = gwas_in$bytes
    ),
    list(
      role = "dossier",
      source = paste0("dossier:", paste(unique(dossier$source), collapse = ",")),
      name = paste0(nrow(dossier), " dossier row(s)"),
      path = dossier_in$path, sha256 = dossier_in$sha256, bytes = dossier_in$bytes
    )
  )
  if (!is.null(candidates)) {
    out[[length(out) + 1L]] <- list(
      role = "candidates", source = "variants:table",
      name = paste0(nrow(candidates$data), " candidate(s)"),
      path = candidates$path, sha256 = candidates$sha256, bytes = candidates$bytes
    )
  }
  out
}

#' The crate-level record written to marker-crate.json
#' @keywords internal
#' @noRd
crate_record <- function(artefact_dir, artefact, gwas_in, gwas, dossier_in,
                         dossier, selection, selected, panel, specificity,
                         build, index, maf_min, ld_min, created) {
  sidecars <- c(
    "gwas-loci.csv", "dossier.csv", "selection.csv", "specificity.csv", "artefact.json"
  )
  files <- lapply(sidecars, function(file) {
    path <- fs::path(artefact_dir, file)
    list(
      name = file, sha256 = hash_file(path),
      bytes = as.integer(file.size(path))
    )
  })
  selected_rows <- lapply(seq_len(nrow(selection)), function(i) {
    lapply(selection, function(col) col[[i]])
  })
  record <- list(
    schema = "tidymarkers/crate@1",
    contract = list(
      version = "tidymarkers/crate-contract@1",
      panel = "tidymarkers/contract@1"
    ),
    created = created,
    gwas = list(
      runs = I(unique(gwas$gwas_run)), population = I(unique(gwas$population)),
      genome_build = build, rows = nrow(gwas), path = gwas_in$path,
      sha256 = gwas_in$sha256
    ),
    dossier = list(
      rows = nrow(dossier), sources = I(unique(dossier$source)),
      versions = I(unique(dossier$source_version)),
      retrieved_at = I(unique(dossier$retrieved_at)),
      confidence = list(min = min(dossier$confidence), max = max(dossier$confidence)),
      path = dossier_in$path, sha256 = dossier_in$sha256
    ),
    selection = list(
      maf_min = maf_min, ld_min = ld_min,
      candidates = selected_rows,
      selected = I(as.character(selected$label))
    ),
    panel = list(
      artefact_id = artefact$id,
      targets = nrow(panel$targets),
      designed = sum(panel$targets$ok),
      one_pot = list(
        groups = panel$one_pot$n_groups,
        members = lapply(seq_len(nrow(panel$one_pot$groups)), function(i) {
          lapply(panel$one_pot$groups, function(col) col[[i]])
        }),
        unassigned = I(as.character(panel$one_pot$unassigned$label))
      )
    ),
    specificity = list(
      rows = nrow(specificity),
      hits = if (nrow(specificity)) sum(specificity$hits) else 0L,
      max_hits = if (nrow(specificity)) max(specificity$hits) else 0L,
      builds = I(unique(specificity$build)),
      engine = I(unique(specificity$engine)),
      index = if (is.list(index) && !is.object(index)) {
        I(unname(vapply(index, function(x) x$strategy %||% NA_character_, character(1))))
      } else {
        index$strategy %||% NA_character_
      }
    ),
    files = files,
    limitations = I(c(
      "one selected marker per locus",
      "specificity hit counts are advisory; no threshold is applied",
      "a bound on an unmeasured MAF or LD aborts rather than passing",
      "the one-pot grouping is not a validated wet-lab mix"
    ))
  )
  record$id <- hash_object(jsonlite::toJSON(
    record[c("schema", "contract", "created", "gwas", "dossier", "selection", "panel", "specificity")],
    auto_unbox = TRUE, null = "null", digits = NA
  ))
  record
}
