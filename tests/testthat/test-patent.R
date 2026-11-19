# WIPO ST.26 export, draft application and defensive publication (roadmap ZY8).

patent_rc <- function(x) {
  bases <- rev(strsplit(x, "", fixed = TRUE)[[1L]])
  chartr("ACGT", "TGCA", paste(bases, collapse = ""))
}

patent_toggle <- function(x) {
  last <- substr(x, nchar(x), nchar(x))
  replacement <- if (identical(last, "A")) "C" else "A"
  paste0(substr(x, 1L, nchar(x) - 1L), replacement)
}

# 300 bp window; the SNP template index is 150, the genomic position 1000.
patent_window <- strrep("ACGT", 75L)

patent_row <- function(as_ref, as_alt, com, as_start, as_end, com_start, com_end, strand) {
  tibble::tibble(
    as_ref = as_ref,
    as_alt = as_alt,
    com = com,
    as_start = as_start,
    as_end = as_end,
    com_start = com_start,
    com_end = com_end,
    strand = strand,
    amplicon = 151L
  )
}

patent_marker <- function(row) {
  structure(
    list(
      scheme = "kasp",
      snp = list(chrom = "Chr1", pos = 1000L),
      assay = row,
      template = list(seq = Biostrings::DNAString(patent_window), snp_i = 150L)
    ),
    class = "marker_assay"
  )
}

patent_plus <- patent_marker(patent_row(
  as_ref = substr(patent_window, 41L, 60L),
  as_alt = patent_toggle(substr(patent_window, 41L, 60L)),
  com = patent_rc(substr(patent_window, 171L, 191L)),
  as_start = 891L,
  as_end = 910L,
  com_start = 1021L,
  com_end = 1041L,
  strand = "+"
))

patent_minus <- patent_marker(patent_row(
  as_ref = patent_rc(substr(patent_window, 81L, 100L)),
  as_alt = patent_toggle(patent_rc(substr(patent_window, 81L, 100L))),
  com = substr(patent_window, 211L, 231L),
  as_start = 931L,
  as_end = 950L,
  com_start = 1061L,
  com_end = 1081L,
  strand = "-"
))

patent_solved <- structure(
  list(report = tibble::tibble(marker = list(patent_plus, patent_minus))),
  class = "marker_solved"
)

test_that("marker_patent_export() writes a conforming ST.26 draft", {
  out <- tempfile(fileext = ".xml")
  suppressWarnings(marker_patent_export(patent_plus, format = "st26", output = out))
  expect_true(fs::file_exists(out))

  doc <- xml2::read_xml(out)
  expect_identical(xml2::xml_name(doc), "ST26SequenceListing")
  attrs <- xml2::xml_attrs(doc)
  expect_true(all(
    c("dtdVersion", "fileName", "softwareName", "softwareVersion", "productionDate") %in% names(attrs)
  ))
  expect_identical(attrs[["dtdVersion"]], "V1_3")
  expect_identical(attrs[["softwareName"]], "tidymarkers")
  expect_match(attrs[["productionDate"]], "^\\d{4}-\\d{2}-\\d{2}$")

  expect_identical(xml2::xml_text(xml2::xml_find_first(doc, "//SequenceTotalQuantity")), "4")
  expect_length(xml2::xml_find_all(doc, "//SequenceData"), 4L)
  expect_length(xml2::xml_find_all(doc, "//INSDSeq"), 4L)
  expect_length(xml2::xml_find_all(doc, "//Feature"), 4L)
  expect_length(xml2::xml_find_all(doc, "//INSDSeq_sequence"), 4L)

  expect_identical(unique(xml2::xml_text(xml2::xml_find_all(doc, "//INSDSeq_moltype"))), "DNA")
  expect_identical(unique(xml2::xml_text(xml2::xml_find_all(doc, "//INSDSeq_division"))), "PAT")
  expect_identical(
    xml2::xml_text(xml2::xml_find_all(doc, "//Feature_key")),
    c("primer_bind", "primer_bind", "primer_bind", "misc_feature")
  )

  seqs <- xml2::xml_text(xml2::xml_find_all(doc, "//INSDSeq_sequence"))
  lens <- as.integer(xml2::xml_text(xml2::xml_find_all(doc, "//INSDSeq_length")))
  expect_identical(lens, nchar(seqs))
  expect_identical(
    xml2::xml_attr(xml2::xml_find_all(doc, "//SequenceData"), "sequenceIDNumber"),
    c("1", "2", "3", "4")
  )
})

test_that("the ST.26 sequences are the assay's primers and amplicon", {
  out_plus <- tempfile(fileext = ".xml")
  suppressWarnings(marker_patent_export(patent_plus, "st26", out_plus))
  seqs_plus <- xml2::xml_text(xml2::xml_find_all(xml2::read_xml(out_plus), "//INSDSeq_sequence"))
  expect_identical(
    seqs_plus,
    c(
      patent_plus$assay$as_ref,
      patent_plus$assay$as_alt,
      patent_plus$assay$com,
      substr(patent_window, 41L, 191L)
    )
  )

  out_minus <- tempfile(fileext = ".xml")
  suppressWarnings(marker_patent_export(patent_minus, "st26", out_minus))
  seqs_minus <- xml2::xml_text(xml2::xml_find_all(xml2::read_xml(out_minus), "//INSDSeq_sequence"))
  expect_identical(
    seqs_minus,
    c(
      patent_minus$assay$as_ref,
      patent_minus$assay$as_alt,
      patent_minus$assay$com,
      patent_rc(substr(patent_window, 81L, 231L))
    )
  )
  expect_identical(seqs_minus[[3L]], substr(patent_window, 211L, 231L))
})

test_that("a solved plan numbers every designed marker contiguously", {
  out <- tempfile(fileext = ".xml")
  suppressWarnings(marker_patent_export(patent_solved, "st26", out))
  doc <- xml2::read_xml(out)
  expect_identical(xml2::xml_text(xml2::xml_find_first(doc, "//SequenceTotalQuantity")), "8")
  expect_identical(
    xml2::xml_attr(xml2::xml_find_all(doc, "//SequenceData"), "sequenceIDNumber"),
    as.character(1:8)
  )
  expect_length(xml2::xml_find_all(doc, "//INSDSeq_sequence"), 8L)
})

test_that("marker_patent_export() refuses input it cannot list", {
  no_template <- patent_plus
  no_template$template <- NULL
  expect_error(
    marker_patent_export(no_template, "st26", tempfile(fileext = ".xml")),
    "amplicon"
  )

  degenerate <- patent_plus
  degenerate$assay$as_ref <- paste0(substr(degenerate$assay$as_ref, 1L, 19L), "N")
  expect_error(
    marker_patent_export(degenerate, "st26", tempfile(fileext = ".xml")),
    "AS_ref.*non-ACGT"
  )

  expect_error(
    marker_patent_export(list(), "st26", tempfile(fileext = ".xml")),
    "marker_assay"
  )

  unsolved <- structure(
    list(report = tibble::tibble(marker = list(NULL))),
    class = "marker_solved"
  )
  expect_error(
    marker_patent_export(unsolved, "st26", tempfile(fileext = ".xml")),
    "No designed marker"
  )

  expect_error(marker_patent_export(patent_plus, format = "pdf"), "st26")
})

test_that("export warns once when the ST.26 XSD is not cached", {
  skip_if(fs::file_exists(patent_schema_path()), "the ST.26 XSD is cached")
  warnings <- testthat::capture_warnings(
    marker_patent_export(patent_plus, "st26", tempfile(fileext = ".xml"))
  )
  expect_length(warnings, 1L)
  expect_match(warnings[[1L]], "ST26SequenceListing", fixed = TRUE)
})

test_that("marker_patent_validate() names the schema it looked for", {
  listing <- tempfile(fileext = ".xml")
  writeLines("<ST26SequenceListing/>", listing)
  expect_error(
    marker_patent_validate(listing, schema = tempfile(fileext = ".xsd")),
    "ST26SequenceListing_V1_3.xsd"
  )
  expect_error(marker_patent_validate(tempfile(fileext = ".xml")), "No sequence listing")

  expect_identical(
    as.character(patent_schema_path()),
    as.character(fs::path(marker_cache_dir(), "ST26SequenceListing_V1_3.xsd"))
  )

  skip_if(fs::file_exists(patent_schema_path()), "the ST.26 XSD is cached")
  err <- tryCatch(marker_patent_validate(listing), error = function(e) e)
  msg <- cli::ansi_strip(conditionMessage(err))
  expect_match(msg, "ST26SequenceListing_V1_3.xsd", fixed = TRUE)
  expect_match(msg, as.character(marker_cache_dir()), fixed = TRUE)
})

test_that("a listing validates when the ST.26 XSD is cached", {
  skip_if_not_installed("xml2")
  skip_if(
    !fs::file_exists(patent_schema_path()),
    paste0("ST.26 XSD not cached at ", patent_schema_path())
  )
  out <- tempfile(fileext = ".xml")
  expect_silent(marker_patent_export(patent_plus, "st26", out))
  expect_true(marker_patent_validate(out))
})

test_that("the application skeleton carries every required section", {
  text <- patent_application_text(patent_markers(patent_plus))
  headings <- c(
    "Background", "Summary of invention", "Detailed description",
    "Figures reference", "Sequence listing reference", "Claims skeleton"
  )
  for (heading in headings) {
    expect_match(text, paste0("# ", heading), fixed = TRUE)
  }
  expect_match(text, tidymarkers:::PATENT_DISCLAIMER, fixed = TRUE)
  expect_false(grepl("{{", text, fixed = TRUE))
  expect_match(text, "SEQ ID NO 1-4", fixed = TRUE)
  expect_match(text, patent_plus$assay$as_ref, fixed = TRUE)
  expect_match(text, "Chr1:1000", fixed = TRUE)
})

test_that("the docx export writes the required sections", {
  skip_if_not_installed("officer")
  skip_if_not_installed("xml2")
  out <- tempfile(fileext = ".docx")
  suppressWarnings(marker_patent_export(patent_plus, "docx", out))
  dir <- tempfile("docx")
  dir.create(dir)
  unzip(out, "word/document.xml", exdir = dir)
  text <- xml2::xml_text(xml2::read_xml(file.path(dir, "word/document.xml")))
  expect_match(text, "Background", fixed = TRUE)
  expect_match(text, "Summary of invention", fixed = TRUE)
  expect_match(text, "Detailed description", fixed = TRUE)
  expect_match(text, "Figures reference", fixed = TRUE)
  expect_match(text, "Sequence listing reference", fixed = TRUE)
  expect_match(text, "Claims skeleton", fixed = TRUE)
  expect_match(text, tidymarkers:::PATENT_DISCLAIMER, fixed = TRUE)
})

test_that("the docx export aborts when officer is absent", {
  skip_if(requireNamespace("officer", quietly = TRUE), "officer is installed")
  expect_error(
    marker_patent_export(patent_plus, "docx", tempfile(fileext = ".docx")),
    "officer"
  )
})

test_that("marker_defensive_publication() frames the markers as prior art", {
  for (obj in list(patent_plus, patent_solved)) {
    pub <- marker_defensive_publication(obj)
    expect_type(pub, "character")
    expect_length(pub, 1L)
    expect_match(pub, tidymarkers:::PATENT_DISCLAIMER, fixed = TRUE)
    expect_match(pub, "prior art", ignore.case = TRUE)
    expect_match(pub, "Zenodo", fixed = TRUE)
    expect_match(pub, "IP.com", fixed = TRUE)
    expect_match(pub, "Research Disclosure", fixed = TRUE)
    expect_match(pub, "YYYY-MM-DD", fixed = TRUE)
    expect_match(pub, "PUBLICATION-ID", fixed = TRUE)
    expect_false(grepl("{{", pub, fixed = TRUE))
    expect_false(grepl("http", pub, ignore.case = TRUE))
  }
})
