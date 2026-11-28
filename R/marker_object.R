# =============================================================================
#  marker_object.R  --  assay predicates, print method, batch design
# =============================================================================
#' Is `x` a KASP assay object?
#'
#' @param x object to test.
#' @return a single logical.
#' @examples
#' is_marker_assay(list())
#' @family utilities
#' @export
is_marker_assay <- function(x) inherits(x, "marker_assay")

#' Do all acceptance checks pass?
#'
#' @param assay assembled assay.
#' @param limits thresholds (see [marker_design_checks()]).
#' @return a single logical, `TRUE` when every check passes.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' assay <- design_hrm_marker(genome, loci)
#' # FALSE here: the off-target checks are `unknown` without a BLAST index
#' passes_all(assay)
#' @family checks
#' @export
passes_all <- function(assay, limits = KASP_LIMITS) {
  all(marker_design_checks(assay, limits)$pass)
}

#' AS-alt oligo sequence as ordered (the degenerate consensus of a multi-allelic
#' site), falling back to the concrete member for assays built before
#' `as_alt_degen` existed.
#' @param assay assembled assay.
#' @return a single sequence string.
#' @keywords internal
#' @noRd
as_alt_sequence <- function(assay) {
  aa <- assay$assay
  if ("as_alt_degen" %in% names(aa)) aa$as_alt_degen[1] else aa$as_alt[1]
}

#' Print a KASP assay
#'
#' Compact one-screen summary: target, the three oligos, amplicon size,
#' discrimination and the verdict (see `evaluation_verdict()`), so the first thing
#' a user reads already answers "did this design work?".
#'
#' @param x a [assemble_kasp_assay()] result.
#' @param ... unused.
#' @return `x`, invisibly.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci)
#' print(assay)
#' }
#' @family utilities
#' @noRd
marker_print_assay <- function(x, ...) {
  a <- x$assay[1L, ]
  cat(
    "<marker_assay>",
    glue::glue("{x$snp$chrom}:{x$snp$pos} {x$snp$ref}>{x$snp$alt} (strand {x$strand})"), "\n"
  )
  cat("  AS-ref ", a$as_ref, "\n")
  cat("  AS-alt ", as_alt_sequence(x), "\n")
  cat("  COM    ", a$com, "\n")
  cat(" ", glue::glue("amplicon {a$amplicon} bp | discrimination {round(x$efficiency$discrimination[1L], 3)}"), "\n")
  cat(" ", glue::glue("verdict {evaluation_verdict(x)}"), "\n")
  invisible(x)
}

#' Print a doomed target
#' @param x a `marker_doomed` object (see [assemble_kasp_assay()]).
#' @param ... unused.
#' @return `x`, invisibly.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' # no candidate KASP pair survives this AT-rich window, so the target comes
#' # back doomed rather than as an error: the design is returned, not thrown away
#' \dontrun{
#' print(assemble_kasp_assay(genome, loci))
#' }
#' @family utilities
#' @noRd
marker_print_doomed <- function(x, ...) {
  cat(
    "<marker_doomed>",
    glue::glue("{x$snp$chrom}:{x$snp$pos} {x$snp$ref}>{x$snp$alt} (stage {x$stage})"), "\n"
  )
  cat(" ", x$reason, "\n")
  invisible(x)
}


# =============================================================================
#  7b. Batch design over many loci
# =============================================================================

#' Classify a variant by its REF/ALT strings
#' @param ref,alt allele strings (recycled).
#' @return a character vector: `"SNP"`, `"MNP"`, `"insertion"` or `"deletion"`.
#' @examples
#' variant_type("A", "G") # SNP
#' variant_type("ACGT", "TCGA") # MNP
#' variant_type("A", "ACGT") # insertion
#' @family design
#' @export
variant_type <- function(ref, alt) {
  dplyr::if_else(
    nchar(ref) == nchar(alt),
    dplyr::if_else(nchar(ref) == 1L, "SNP", "MNP"),
    dplyr::if_else(nchar(alt) > nchar(ref), "insertion", "deletion")
  )
}

#' Design KASP markers for a table of loci, with a per-locus verdict
#'
#' Runs [assemble_kasp_assay()] independently over every row of `loci`, so one
#' hard locus cannot abort the batch: an error is captured as
#' `verdict = "error"` (a design that ran and failed is `"fail"`), warnings or
#' failing acceptance checks give `verdict = "warn"`, and a fully passing assay
#' gives `verdict = "pass"`.  `loci` needs columns
#' `chrom`, `pos`, `ref`, `alt` (as in [read_variants_vcf()]).
#'
#' This is the KASP batch: it runs the full KASP assembly and acceptance checks,
#' so it returns the per-locus verdict columns (`type`, `strand`, `amplicon`,
#' `discrimination`, `n_fail`, `verdict`, `note`).  For "which chemistry fits this
#' locus", use [design_auto_marker()]; for a thin one-row-per-locus summary across
#' any chemistry, use [design_markers()].
#'
#' @param genome DNAStringSet.
#' @param loci tibble of loci (`chrom`, `pos`, `ref`, `alt`).
#' @param variants variant table for footprint accounting; defaults to `loci`.
#' @param BPPARAM optional `BiocParallel` backend; when given, loci are processed
#'   in parallel (a forked backend avoids copying the genome).  The per-locus
#'   subprocess work is then forced serial, so cores are not oversubscribed; with
#'   no backend, `cores` is forwarded to [assemble_kasp_assay()] instead.
#' @param ... passed on to [assemble_kasp_assay()].
#' @return tidy tibble, one row per locus.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' loci <- tibble::tibble(chrom = "Chr4", pos = 269090L, ref = "G", alt = "A")
#' design_kasp_markers(genome, loci)
#' }
#' @family design
#' @export
design_kasp_markers <- function(genome, loci, variants = loci, BPPARAM = NULL, ...) {
  need_cols(loci, c("chrom", "pos", "ref", "alt"), arg = "loci")
  one <- function(ch, pos, ref, alt, cores) {
    # `locus_row()` judges the position, so a bad row becomes an error verdict
    # below instead of aborting the whole batch (ZV105).
    snp <- NULL
    # one schema for every verdict, so the doomed/ok/error rows cannot drift
    row <- function(strand = NA_character_, amplicon = NA_integer_, discrimination = NA_real_,
                    n_fail = NA_integer_, verdict = "fail", note = NA_character_) {
      chrom_v <- if (is.null(snp)) as.character(ch) else snp$chrom
      pos_v <- if (is.null(snp)) NA_integer_ else snp$pos
      ref_v <- if (is.null(snp)) as.character(ref) else snp$ref
      alt_v <- if (is.null(snp)) as.character(alt) else snp$alt
      tibble::tibble(
        chrom = chrom_v, pos = pos_v, ref = ref_v, alt = alt_v,
        type = variant_type(ref_v, alt_v),
        strand = strand, amplicon = amplicon, discrimination = discrimination,
        n_fail = n_fail, verdict = verdict, note = note
      )
    }
    cap <- new.env(parent = emptyenv())
    cap$warns <- character(0)
    dots <- rlang::list2(...)
    dots$cores <- cores
    res <- withCallingHandlers(
      tryCatch(
        {
          snp <- locus_row(tibble::tibble(chrom = ch, pos = pos, ref = ref, alt = alt))
          a <- rlang::exec(assemble_kasp_assay, genome, snp, variants = variants, !!!dots)
          if (is_kasp_doomed(a)) {
            list(row = row(note = a$reason), marker = NULL)
          } else {
            chk <- marker_design_checks(a)
            list(
              row = row(
                strand = a$strand, amplicon = a$assay$amplicon[1],
                discrimination = a$efficiency$discrimination[1],
                n_fail = sum(!chk$pass),
                verdict = if (all(chk$pass)) "pass" else "warn",
                note = paste(chk$criterion[!chk$pass], collapse = "; ")
              ),
              marker = a
            )
          }
        },
        error = function(e) {
          list(
            row = row(verdict = "error", note = conditionMessage(e)),
            marker = NULL
          )
        }
      ),
      warning = function(w) {
        cap$warns <- c(cap$warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    # surface any assembly warnings on an otherwise built assay
    if (res$row$verdict %in% c("pass", "warn") && length(cap$warns)) {
      res$row$verdict <- if (res$row$verdict == "pass") "warn" else res$row$verdict
      res$row$note <- paste(c(res$row$note, cap$warns), collapse = "; ")
    }
    res
  }

  if (is.null(BPPARAM)) {
    n <- nrow(loci)
    id <- cli::cli_progress_bar(gettext("Designing KASP markers"), total = n)
    parts <- purrr::imap(seq_len(n), function(i, ...) {
      cli::cli_progress_update(id = id)
      one(loci$chrom[i], loci$pos[i], loci$ref[i], loci$alt[i],
        cores = getOption("tidymarkers.cores", 1L)
      )
    })
    cli::cli_progress_done(id = id)
  } else {
    # one parent build, then every worker only opens an existing shard
    index_prebuild(genome)
    parts <- bp_lapply(seq_len(nrow(loci)), function(i) {
      worker_biostrings()
      one(loci$chrom[i], loci$pos[i], loci$ref[i], loci$alt[i], cores = 1L)
    }, BPPARAM)
  }
  # the assays travel as an attribute: the verdict schema stays exactly the one
  # callers already parse, while `as_marker_solved()` can still reach the designs
  out <- purrr::list_rbind(lapply(parts, `[[`, "row"))
  attr(out, "markers") <- lapply(parts, `[[`, "marker")
  class(out) <- c("marker_batch", class(out))
  out
}
