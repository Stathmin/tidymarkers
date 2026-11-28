# =============================================================================
#  marker_patent.R  --  WIPO ST.26 sequence listings and draft applications
# =============================================================================
#  Roadmap ZY8: turn the markers a design already produced into the documents a
#  patent filing or a defensive publication needs -- a WIPO ST.26 XML sequence
#  listing and a skeleton application.  Every claim is read from the assay
#  object (`$assay`, `$snp`, `$template`); nothing is re-derived and nothing is
#  invented: applicant, inventor and priority stay out of the draft, and the
#  prose sections are bracketed instructions for a human to fill in.
#
#  The ST.26 XSD does not ship with the package (external data never ships), so
#  `marker_patent_validate()` reads it from the per-user cache directory and
#  aborts with the expected path when it is absent.  The export itself still
#  writes the XML and warns once, so a missing schema costs validation, not the
#  document.
# =============================================================================

# Disclaimer carried by every generated patent document: the draft documents
# this package writes are technical drafts, not legal advice.  Internal on
# purpose -- dev/graphify-health.R only registers roxygen blocks for function
# definitions, so a documented constant reads as an undocumented export.
PATENT_DISCLAIMER <- "draft documents; patent filing requires qualified patent attorney review"

#' Export markers as a WIPO ST.26 sequence listing or a draft application
#'
#' `format = "st26"` writes a WIPO ST.26 XML sequence listing: one `SequenceData`
#' per claimed sequence (allele-specific primers for the reference and the
#' alternative allele, the common primer, then the amplicon, each 5' to 3' and
#' numbered contiguously from 1), inside an `INSDSeq` record with
#' `INSDSeq_moltype` `DNA`, `INSDSeq_division` `PAT` and a `FeatureTable`.
#' `format = "docx"` renders the skeleton application from
#' `inst/templates/patent/application.md`: the sections *Background*, *Summary of
#' invention*, *Detailed description*, *Figures reference*, *Sequence listing
#' reference* and *Claims skeleton*, plus the disclaimer.  Both formats are
#' drafts: the prose sections are bracketed instructions, and no applicant,
#' inventor, priority date or sequence identifier is invented.
#'
#' @section Sequence numbering:
#' Sequences are read from the assay object: `$assay$as_ref`, `$assay$as_alt`,
#' `$assay$com` and the amplicon.  The amplicon is taken from
#' `$assay$amplicon_seq` when present; otherwise it is sliced from
#' `$template$seq` between the outer primer coordinates and, for an assay on the
#' minus strand, reverse-complemented so every sequence is listed 5' to 3' in
#' the orientation of its assay.  Degenerate (IUPAC) primers cannot be listed
#' unambiguously and abort with the offending role.
#'
#' @section Sequence listing schema:
#' WIPO publishes the ST.26 schema; it does not ship with this package (external
#' data never ships).  Put `ST26SequenceListing_V1_3.xsd` in the per-user cache
#' directory, `tools::R_user_dir("tidymarkers", "cache")` (see
#' [marker_cache_dir()]), then validation runs automatically on export and on
#' demand with [marker_patent_validate()].  The standards page is
#' <https://www.wipo.int/standards/en/>; the exact archive URL is supplied by the
#' operator, and the one-line install is
#' `curl -L -o "<cache>/ST26SequenceListing_V1_3.xsd" "<WIPO ST.26 XSD URL>"`.
#' Without the XSD the XML is still written, with a warning that it was not
#' validated.
#'
#' @param assay_or_solved a [marker_assay][assemble_kasp_assay] object, or a
#'   `marker_solved` from [marker_solve()] (only its designed markers are
#'   exported; unsolved targets are skipped).
#' @param format output format, matched with [match.arg()]: `"st26"` for the XML
#'   sequence listing, `"docx"` for the application skeleton.
#' @param output output path.  Defaults to `"sequence_listing.xml"` for ST.26 and
#'   `"application_draft.docx"` for the docx; the parent directory is created.
#' @param ... unused; reserved for future formats.
#' @return the output path, invisibly.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci[1, ])
#' marker_patent_export(assay, format = "st26", output = tempfile(fileext = ".xml"))
#' marker_patent_export(assay, format = "docx", output = tempfile(fileext = ".docx"))
#' }
#' @family patent
#' @export
marker_patent_export <- function(assay_or_solved, format = c("st26", "docx"), output = NULL, ...) {
  formats <- c("st26", "docx")
  format <- tryCatch(
    match.arg(format, formats),
    error = function(e) {
      abort_msg(
        "Unsupported `format` {.val {format}}.",
        i = "Supported formats: {.val {formats}}."
      )
    }
  )
  if (is.null(output)) {
    output <- switch(format,
      st26 = "sequence_listing.xml",
      docx = "application_draft.docx"
    )
  }
  output <- fs::path(output)
  refuse_links(fs::path_dir(output), fs::path_file(output))
  fs::dir_create(fs::path_dir(output))
  markers <- patent_markers(assay_or_solved)

  if (identical(format, "st26")) {
    patent_require("xml2", gettext("writing a WIPO ST.26 sequence listing"))
    patent_st26_doc(patent_claims(markers), output)
    schema <- patent_schema_path()
    if (fs::file_exists(schema)) {
      marker_patent_validate(output, schema)
    } else {
      cli::cli_warn(c(
        gettext("WIPO ST.26 XML written without schema validation."),
        i = gettext(
          gettext(
            "The XSD was not found at {.file {schema}}; ?marker_patent_export gives the filename and curl command."
          )
        )
      ))
    }
    return(invisible(output))
  }

  patent_require("officer", gettext("building the application draft"))
  patent_docx_write(patent_application_text(markers), output)
  invisible(output)
}

#' Validate a WIPO ST.26 sequence listing against the XSD
#'
#' Reads the schema from the per-user cache directory by default (the XSD does
#' not ship with the package) and validates the listing with [xml2::xml_validate()].
#' A missing listing or a missing schema aborts with the path that was looked
#' for, so the operator knows which half is absent.
#'
#' @param path a sequence listing written by [marker_patent_export()].
#' @param schema path to `ST26SequenceListing_V1_3.xsd`; defaults to
#'   `tools::R_user_dir("tidymarkers", "cache")/ST26SequenceListing_V1_3.xsd`.
#' @return `TRUE`, invisibly; a listing that does not validate aborts with the
#'   schema's error messages.
#' @examples
#' \dontrun{
#' path <- marker_patent_export(assay, "st26", tempfile(fileext = ".xml"))
#' marker_patent_validate(path)
#' }
#' @family patent
#' @export
marker_patent_validate <- function(path, schema = patent_schema_path()) {
  if (!fs::file_exists(path)) {
    abort_msg("No sequence listing at {.file {path}}.")
  }
  if (!fs::file_exists(schema)) {
    abort_msg(
      "The WIPO ST.26 XSD was not found at {.file {schema}}.",
      i = "Download ST26SequenceListing_V1_3.xsd into the cache directory {.file {marker_cache_dir()}}.",
      i = "See ?marker_patent_export for the expected filename and the one-line curl command."
    )
  }
  patent_require("xml2", gettext("validating a WIPO ST.26 sequence listing"))
  ok <- xml2::xml_validate(xml2::read_xml(path), xml2::read_xml(schema))
  if (!isTRUE(ok)) {
    cli::cli_abort(c(
      gettext("The sequence listing does not validate against {.file {schema}}."),
      attr(ok, "errors")
    ))
  }
  invisible(TRUE)
}

#' Draft a defensive publication from the designed markers
#'
#' Renders the defensive publication template with the markers' claims: a
#' prior-art disclosure rather than an application, dated by a
#' placeholder, carrying the same `PATENT_DISCLAIMER` and naming the usual
#' destinations (Zenodo, IP.com, Research Disclosure) as text instructions.  The
#' publication identifier is a placeholder -- this function performs no network
#' call and minted nothing.
#'
#' @inheritParams marker_patent_export
#' @param assay_or_solved an assay from [assemble_kasp_assay()], a solved
#'   design, or a list of either.
#' @return a length-one character vector of markdown.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' cat(marker_defensive_publication(design_asq_marker(genome, loci[1, ])))
#' }
#' @family patent
#' @export
marker_defensive_publication <- function(assay_or_solved, ...) {
  claims <- patent_claims(patent_markers(assay_or_solved))
  patent_render_template(
    patent_read_template("defensive_publication.md"),
    patent_render_values(claims)
  )
}

# =============================================================================
#  Internal helpers
# =============================================================================

# The XSD does not ship with the package; this is where validation looks.
patent_schema_path <- function() {
  fs::path(marker_cache_dir(), "ST26SequenceListing_V1_3.xsd")
}

# Suggests-only packages stay behind a guard so R CMD check with
# _R_CHECK_DEPENDS_ONLY_ never reaches an unguarded `::`.
patent_require <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    abort_msg(
      "The {.pkg {pkg}} package is required for {why}.",
      i = "Install it with {.run install.packages(\"{pkg}\")}."
    )
  }
  invisible(TRUE)
}

# One marker from either an assay or a solved plan; unsolved rows are skipped.
patent_markers <- function(assay_or_solved) {
  if (inherits(assay_or_solved, "marker_assay")) {
    return(list(assay_or_solved))
  }
  if (inherits(assay_or_solved, "marker_solved")) {
    markers <- Filter(Negate(is.null), assay_or_solved$report$marker)
    if (!length(markers)) {
      abort_msg("No designed marker to export: every target in this marker_solved is unsolved.")
    }
    return(markers)
  }
  abort_msg(
    "`assay_or_solved` must be a marker_assay or a marker_solved.",
    i = "Got {.cls {class(assay_or_solved)}}."
  )
}

# The one-row assay table, checked for the columns every claim needs.
patent_assay_row <- function(marker, i) {
  row <- marker$assay
  needed <- c("as_ref", "as_alt", "com")
  missing <- setdiff(needed, names(row))
  if (!is.data.frame(row) || nrow(row) < 1L || length(missing)) {
    abort_msg(
      "Marker {i} carries no usable assay table.",
      i = "A marker_assay needs a one-row `$assay` with {.field {needed}}; missing: {.field {missing}}."
    )
  }
  row[1L, , drop = FALSE]
}

patent_locus <- function(marker) {
  snp <- marker$snp
  chrom <- snp$chrom[1L]
  pos <- snp$pos[1L]
  if (is.null(chrom) || is.null(pos) || is.na(chrom) || is.na(pos)) {
    return(NA_character_)
  }
  paste0(chrom, ":", as.integer(pos))
}

# Genomic coordinates to template-window coordinates (the assay stores the
# former, the template the latter).
patent_template_offset <- function(marker, tmpl) {
  snp <- marker$snp
  if (is.null(snp) || is.null(snp$pos) || is.null(tmpl$snp_i)) {
    return(0L)
  }
  as.integer(snp$pos[1L]) - as.integer(tmpl$snp_i)
}

patent_amplicon <- function(marker, assay_row, i) {
  if ("amplicon_seq" %in% names(assay_row)) {
    stored <- assay_row$amplicon_seq[[1L]]
    if (!is.null(stored) && !is.na(stored) && nzchar(stored)) {
      return(as.character(stored))
    }
  }
  tmpl <- marker$template
  if (is.null(tmpl) || is.null(tmpl$seq)) {
    abort_msg(
      "Marker {i} has no template sequence, so the amplicon cannot be reconstructed.",
      i = "Keep `$template$seq` on the marker, or store the amplicon in `$assay$amplicon_seq`."
    )
  }
  cols <- c("as_start", "as_end", "com_start", "com_end")
  missing <- setdiff(cols, names(assay_row))
  if (length(missing)) {
    abort_msg(
      "Marker {i} has no primer coordinates, so the amplicon cannot be reconstructed.",
      i = "Missing `$assay` column(s): {.field {missing}}."
    )
  }
  off <- patent_template_offset(marker, tmpl)
  lo <- min(assay_row$as_start, assay_row$com_start) - off
  hi <- max(assay_row$as_end, assay_row$com_end) - off
  n <- length(tmpl$seq)
  if (lo < 1L || hi > n || hi < lo) {
    abort_msg(
      "Marker {i}: the amplicon falls outside the template.",
      i = "Template positions {lo}..{hi} were asked of a {n}-bp template."
    )
  }
  plus <- substr(as.character(tmpl$seq), lo, hi)
  strand <- if ("strand" %in% names(assay_row)) as.character(assay_row$strand[[1L]]) else "+"
  if (identical(strand, "-")) {
    return(as.character(Biostrings::reverseComplement(Biostrings::DNAString(plus))))
  }
  plus
}

# ST.26 needs an unambiguous DNA sequence; the role names the offender.
patent_clean_sequence <- function(x, marker, role) {
  seq <- toupper(as.character(x[[1L]]))
  chars <- strsplit(seq, "", fixed = TRUE)[[1L]]
  bad <- is.na(match(chars, c("A", "C", "G", "T")))
  if (!nzchar(seq) || any(bad)) {
    locus <- patent_locus(marker)
    at <- if (is.na(locus)) "" else paste0(" at ", locus)
    cli::cli_abort(c(
      gettextf("The %s sequence%s is not plain DNA (non-ACGT).", role, at),
      i = gettext("WIPO ST.26 lists unambiguous DNA; degenerate IUPAC codes cannot be exported.")
    ))
  }
  seq
}

# One row per claimed sequence: AS_ref, AS_alt, COM, amplicon, numbered from 1.
patent_claims <- function(markers) {
  rows <- lapply(seq_along(markers), function(i) {
    marker <- markers[[i]]
    row <- patent_assay_row(marker, i)
    locus <- patent_locus(marker)
    roles <- c("AS_ref", "AS_alt", "COM")
    oligos <- list(row$as_ref[[1L]], row$as_alt[[1L]], row$com[[1L]])
    seqs <- vapply(seq_along(roles), function(k) {
      patent_clean_sequence(oligos[[k]], marker, roles[[k]])
    }, character(1))
    amplicon <- patent_clean_sequence(patent_amplicon(marker, row, i), marker, "amplicon")
    role <- c(roles, "amplicon")
    note <- patent_claim_note(i, role, locus)
    tibble::tibble(
      marker = i,
      role = role,
      locus = locus,
      note = note,
      sequence = c(seqs, amplicon)
    )
  })
  claims <- purrr::list_rbind(rows)
  claims$seq_no <- seq_len(nrow(claims))
  claims[, c("seq_no", "marker", "role", "locus", "note", "sequence")]
}

patent_claim_note <- function(i, role, locus) {
  where <- if (is.na(locus)) "" else paste0(" at ", locus)
  paste0(role, " of marker ", i, where, ", designed by tidymarkers")
}

patent_feature_key <- function(role) {
  if (identical(role, "amplicon")) "misc_feature" else "primer_bind"
}

patent_software_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("tidymarkers")),
    error = function(e) "unknown"
  )
}

# Build the ST.26 tree with xml2 (never paste0) and write it.
patent_st26_doc <- function(claims, path) {
  doc <- xml2::xml_new_root(
    "ST26SequenceListing",
    dtdVersion = "V1_3",
    fileName = basename(path),
    softwareName = "tidymarkers",
    softwareVersion = patent_software_version(),
    productionDate = format(Sys.Date(), "%Y-%m-%d")
  )
  xml2::xml_add_child(doc, "SequenceTotalQuantity", as.character(nrow(claims)))
  for (i in seq_len(nrow(claims))) {
    patent_st26_sequence(doc, claims[i, , drop = FALSE])
  }
  xml2::write_xml(doc, path)
  invisible(path)
}

patent_st26_sequence <- function(doc, claim) {
  seq <- claim$sequence[[1L]]
  len <- nchar(seq)
  data <- xml2::xml_add_child(doc, "SequenceData")
  xml2::xml_set_attrs(data, c(sequenceIDNumber = as.character(claim$seq_no[[1L]])))
  insd <- xml2::xml_add_child(data, "INSDSeq")
  xml2::xml_add_child(insd, "INSDSeq_length", as.character(len))
  xml2::xml_add_child(insd, "INSDSeq_moltype", "DNA")
  xml2::xml_add_child(insd, "INSDSeq_division", "PAT")
  table <- xml2::xml_add_child(xml2::xml_add_child(insd, "INSDSeq_feature-table"), "FeatureTable")
  feature <- xml2::xml_add_child(table, "Feature")
  xml2::xml_add_child(feature, "Feature_key", patent_feature_key(claim$role[[1L]]))
  interval <- xml2::xml_add_child(xml2::xml_add_child(feature, "Location"), "Interval")
  xml2::xml_add_child(interval, "Interval_start", "1")
  xml2::xml_add_child(interval, "Interval_end", as.character(len))
  quals <- xml2::xml_add_child(feature, "Feature_quals")
  qualifier <- xml2::xml_add_child(quals, "Feature_qualifier")
  xml2::xml_add_child(qualifier, "Feature_qualifier_name", "note")
  xml2::xml_add_child(qualifier, "Feature_qualifier_value", claim$note[[1L]])
  xml2::xml_add_child(insd, "INSDSeq_sequence", seq)
  invisible(data)
}

patent_template_dir <- function() {
  fs::path(system.file("templates", "patent", package = "tidymarkers"))
}

patent_read_template <- function(name) {
  path <- fs::path(patent_template_dir(), name)
  if (!fs::file_exists(path)) {
    abort_msg("Template {.file {path}} is missing from the installed package.")
  }
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

# Keep a backslash in a claim literal when it lands in a gsub() replacement.
patent_escape_replacement <- function(v) {
  gsub("\\", "\\\\", as.character(v), fixed = TRUE)
}

# Base-R templating: `{{name}}` substitution into a single string.
patent_render_template <- function(lines, values) {
  text <- paste(lines, collapse = "\n")
  for (name in names(values)) {
    text <- gsub(
      paste0("{{", name, "}}"),
      patent_escape_replacement(values[[name]]),
      text,
      fixed = TRUE
    )
  }
  if (grepl("{{", text, fixed = TRUE)) {
    left <- unique(regmatches(text, gregexpr("\\{\\{[^}]*\\}\\}", text))[[1L]])
    cli::cli_abort(c(
      gettextf("Unrendered template placeholder(s): %s.", paste(left, collapse = ", ")),
      i = gettext("Every `{{name}}` in a template needs a value.")
    ))
  }
  text
}

patent_render_values <- function(claims) {
  loci <- unique(claims$locus[!is.na(claims$locus)])
  list(
    disclaimer = PATENT_DISCLAIMER,
    date = "YYYY-MM-DD",
    publication_id = "PUBLICATION-ID",
    n_markers = as.character(length(unique(claims$marker))),
    n_sequences = as.character(nrow(claims)),
    marker_summary = patent_marker_summary(claims),
    marker_loci = if (length(loci)) paste(loci, collapse = ", ") else "not recorded",
    sequence_table = patent_sequence_table(claims)
  )
}

patent_marker_summary <- function(claims) {
  ids <- unique(claims$marker)
  parts <- vapply(ids, function(id) {
    d <- claims[claims$marker == id, , drop = FALSE]
    ids_seq <- range(d$seq_no)
    roles <- paste(d$role, collapse = ", ")
    locus <- unique(d$locus[!is.na(d$locus)])
    locus <- if (length(locus)) locus[[1L]] else "locus not recorded"
    paste0(
      "- Marker ", id, ": ", locus, "; SEQ ID NO ", ids_seq[[1L]], "-", ids_seq[[2L]],
      " (", roles, ")."
    )
  }, character(1))
  paste(parts, collapse = "\n")
}

patent_sequence_table <- function(claims) {
  head <- c(
    "| SEQ ID NO | Marker | Role | Locus | Length | Sequence |",
    "| --- | --- | --- | --- | --- | --- |"
  )
  rows <- vapply(seq_len(nrow(claims)), function(i) {
    locus <- if (is.na(claims$locus[[i]])) "not recorded" else claims$locus[[i]]
    paste0(
      "| ", claims$seq_no[[i]], " | ", claims$marker[[i]], " | ", claims$role[[i]], " | ",
      locus, " | ", nchar(claims$sequence[[i]]), " | ", claims$sequence[[i]], " |"
    )
  }, character(1))
  paste(c(head, rows), collapse = "\n")
}

patent_application_text <- function(markers) {
  claims <- patent_claims(markers)
  patent_render_template(
    patent_read_template("application.md"),
    patent_render_values(claims)
  )
}

# `# ` and `## ` template lines become Word headings; the rest stay paragraphs.
patent_docx_write <- function(text, path) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  doc <- officer::read_docx()
  for (line in lines) {
    if (startsWith(line, "## ")) {
      doc <- officer::body_add_par(doc, substring(line, 4L), style = "heading 2")
    } else if (startsWith(line, "# ")) {
      doc <- officer::body_add_par(doc, substring(line, 3L), style = "heading 1")
    } else if (nzchar(trimws(line))) {
      doc <- officer::body_add_par(doc, line)
    }
  }
  print(doc, target = path)
  invisible(path)
}
