# =============================================================================
#  data-raw/lib/xlsx-read.R  --  the one minimal xlsx reader
# =============================================================================
#  Extracted from makhoul-supplement.R so every supplement parser shares it.
#  No spreadsheet package is installed in this toolchain, so we read the sheet
#  XML directly. Supports shared strings (t="s"), inline strings
#  (t="inlineStr") and numbers; returns list(sheet_name = list(row_cells)).
#
#  Legacy BIFF .xls is NOT supported -- it is not XML. A .xls path is named
#  explicitly by any script that needs it (e.g. liu2012-supplement.R).
# =============================================================================

xlsx_sheets <- function(path) {
  NS <- c(
    d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
  )
  exdir <- tempfile("xlsx")
  utils::unzip(path, exdir = exdir)
  on.exit(unlink(exdir, recursive = TRUE), add = TRUE)

  sst <- character()
  sst_path <- file.path(exdir, "xl", "sharedStrings.xml")
  if (file.exists(sst_path)) {
    doc <- xml2::read_xml(sst_path)
    sst <- vapply(xml2::xml_find_all(doc, ".//d1:si", NS), function(si) {
      paste(xml2::xml_text(xml2::xml_find_all(si, ".//d1:t", NS)), collapse = "")
    }, character(1))
  }

  wb <- xml2::xml_ns_strip(xml2::read_xml(file.path(exdir, "xl", "workbook.xml")))
  rels <- xml2::xml_ns_strip(
    xml2::read_xml(file.path(exdir, "xl", "_rels", "workbook.xml.rels")))
  rel_nodes <- xml2::xml_find_all(rels, ".//Relationship")
  rel <- stats::setNames(xml2::xml_attr(rel_nodes, "Target"),
                         xml2::xml_attr(rel_nodes, "Id"))

  out <- list()
  for (sh in xml2::xml_find_all(wb, ".//sheet")) {
    nm <- xml2::xml_attr(sh, "name")
    path_sheet <- file.path(exdir, "xl",
                            sub("^/?", "", rel[[xml2::xml_attr(sh, "id")]]))
    doc <- xml2::read_xml(path_sheet)
    out[[nm]] <- lapply(xml2::xml_find_all(doc, ".//d1:row", NS), function(row) {
      cells <- xml2::xml_find_all(row, "./d1:c", NS)
      col <- sub("[0-9]+$", "", xml2::xml_attr(cells, "r"))
      type <- xml2::xml_attr(cells, "t")
      v <- xml2::xml_text(xml2::xml_child(cells, "d1:v", NS))
      is_shared <- !is.na(type) & type == "s"
      v[is_shared] <- sst[as.integer(v[is_shared]) + 1L]
      is_inline <- !is.na(type) & type == "inlineStr"
      if (any(is_inline)) {
        v[is_inline] <- vapply(cells[is_inline], function(cc) {
          paste(xml2::xml_text(xml2::xml_find_all(cc, ".//d1:t", NS)), collapse = "")
        }, character(1))
      }
      stats::setNames(as.list(v), col)
    })
  }
  out
}

# NA-aware cell getter: a missing cell and an empty cell are both NA.
xcell <- function(cells, k) {
  if (is.null(cells[[k]])) return(NA_character_)
  v <- cells[[k]]
  if (!length(v) || is.na(v) || !nzchar(v)) NA_character_ else v
}
