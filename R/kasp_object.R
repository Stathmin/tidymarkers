# =============================================================================
#  kasp_object.R  --  assay predicates, print method, batch design
# =============================================================================
#' Is `x` a KASP assay object?
#'
#' @param x object to test.
#' @return a single logical.
#' @examples
#' is_kasp_assay(list())
#' @family KASP utilities
#' @export
is_kasp_assay <- function(x) inherits(x, "kasp_assay")

#' Do all acceptance checks pass?
#'
#' @param assay assembled assay.
#' @param limits thresholds; see [kasp_design_checks()].
#' @return a single logical, `TRUE` when every check passes.
#' @examples
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' passes_all(assay)
#' }
#' @family KASP checks
#' @export
passes_all <- function(assay, limits = KASP_LIMITS) {
  all(kasp_design_checks(assay, limits)$pass)
}

#' Print a KASP assay
#'
#' Compact one-screen summary: target, the three oligos, amplicon size and
#' discrimination.
#'
#' @param x a [assemble_kasp_assay()] result.
#' @param ... unused.
#' @return `x`, invisibly.
#' @examples
#' \dontrun{
#' print(assay)
#' }
#' @method print kasp_assay
#' @family KASP utilities
#' @export
print.kasp_assay <- function(x, ...) {
  a <- x$assay[1L, ]
  cat(
    "<kasp_assay>",
    glue::glue("{x$snp$chrom}:{x$snp$pos} {x$snp$ref}>{x$snp$alt} (strand {x$strand})"), "\n"
  )
  cat("  AS-ref ", a$as_ref, "\n")
  cat("  AS-alt ", a$as_alt, "\n")
  cat("  COM    ", a$com, "\n")
  cat(" ", glue::glue("amplicon {a$amplicon} bp | discrimination {round(x$efficiency$discrimination[1L], 3)}"), "\n")
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
#' @family KASP design
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
#' hard locus cannot abort the batch: errors are captured as `verdict = "fail"`,
#' warnings or failing acceptance checks give `verdict = "warn"`, and a fully
#' passing assay gives `verdict = "pass"`.  `loci` needs columns
#' `chrom`, `pos`, `ref`, `alt` (as in [read_variants_vcf()]).
#'
#' @param genome DNAStringSet.
#' @param loci tibble of loci (`chrom`, `pos`, `ref`, `alt`).
#' @param variants variant table for footprint accounting; defaults to `loci`.
#' @param BPPARAM optional `BiocParallel` backend; when given, loci are processed
#'   in parallel (a forked backend avoids copying the genome).
#' @param ... passed on to [assemble_kasp_assay()].
#' @return tidy tibble, one row per locus.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' loci <- tibble::tibble(chrom = "Chr4", pos = 269090L, ref = "G", alt = "A")
#' design_kasp_markers(genome, loci)
#' }
#' @family KASP design
#' @export
design_kasp_markers <- function(genome, loci, variants = loci, BPPARAM = NULL, ...) {
  need <- c("chrom", "pos", "ref", "alt")
  missing <- setdiff(need, names(loci))
  if (length(missing)) {
    rlang::abort(glue::glue(gettext("`loci` is missing column(s): {paste(missing, collapse = ', ')}")))
  }
  one <- function(ch, pos, ref, alt) {
    snp <- tibble::tibble(
      chrom = as.character(ch), pos = as.integer(pos),
      ref = as.character(ref), alt = as.character(alt)
    )
    type <- variant_type(snp$ref, snp$alt)
    warns <- character(0)
    res <- withCallingHandlers(
      tryCatch(
        {
          a <- assemble_kasp_assay(genome, snp, variants = variants, ...)
          chk <- kasp_design_checks(a)
          tibble::tibble(
            chrom = snp$chrom, pos = snp$pos, ref = snp$ref, alt = snp$alt, type = type,
            strand = a$strand, amplicon = a$assay$amplicon[1],
            discrimination = a$efficiency$discrimination[1],
            n_fail = sum(!chk$pass),
            verdict = if (all(chk$pass)) "pass" else "warn",
            note = paste(chk$criterion[!chk$pass], collapse = "; ")
          )
        },
        error = function(e) {
          tibble::tibble(
            chrom = snp$chrom, pos = snp$pos, ref = snp$ref, alt = snp$alt,
            type = type, strand = NA_character_, amplicon = NA_integer_,
            discrimination = NA_real_, n_fail = NA_integer_,
            verdict = "fail", note = conditionMessage(e)
          )
        }
      ),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    # surface any assembly warnings on an otherwise built assay
    if (res$verdict != "fail" && length(warns)) {
      res$verdict <- if (res$verdict == "pass") "warn" else res$verdict
      res$note <- paste(c(res$note, warns), collapse = "; ")
    }
    res
  }

  if (is.null(BPPARAM)) {
    n <- nrow(loci)
    id <- cli::cli_progress_bar(gettext("Designing KASP markers"), total = n)
    res <- purrr::imap_dfr(seq_len(n), function(i, ...) {
      cli::cli_progress_update(id = id)
      one(loci$chrom[i], loci$pos[i], loci$ref[i], loci$alt[i])
    })
    cli::cli_progress_done(id = id)
    res
  } else {
    res <- bp_lapply(seq_len(nrow(loci)), function(i) {
      worker_biostrings()
      one(loci$chrom[i], loci$pos[i], loci$ref[i], loci$alt[i])
    }, BPPARAM)
    purrr::list_rbind(res)
  }
}
