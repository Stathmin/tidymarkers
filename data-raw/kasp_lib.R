# =============================================================================
#  data-raw/kasp_lib.R  --  one xlsx reader + one write-check, for the KASP trio
# =============================================================================
#  `makhoul-supplement.R` already carries a minimal xml2 reader for the sheet
#  XML; there is no readxl/openxlsx in this toolchain.  The three KASP sources
#  need the same reader, so it is factored here once rather than copied three
#  times.  Sourcing this file has no side effects beyond loading xml2.
#
#  `kasp_sheets()` returns a named list of sheets; each sheet is a list of rows;
#  each row is a named character vector keyed by column letter ("A", "B", ...),
#  with "" for an absent cell.  An absent cell inside a merged range is
#  indistinguishable here -- callers must forward-fill explicitly and then
#  assert the totals against the paper, exactly as the Makhoul script does.
# =============================================================================

KASP_NS <- c(
  d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
  r  = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
)

kasp_sheets <- function(path) {
  for (pkg in c("xml2", "readr", "tibble", "fs")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("the KASP generators need ", pkg)
  }
  if (!file.exists(path)) stop("no such workbook: ", path)
  exdir <- tempfile("xlsx")
  on.exit(unlink(exdir, recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = exdir)

  sst <- character()
  sst_path <- file.path(exdir, "xl", "sharedStrings.xml")
  if (file.exists(sst_path)) {
    doc <- xml2::read_xml(sst_path)
    sst <- vapply(xml2::xml_find_all(doc, ".//d1:si", KASP_NS), function(si) {
      paste(xml2::xml_text(xml2::xml_find_all(si, ".//d1:t", KASP_NS)), collapse = "")
    }, character(1))
  }

  wb <- xml2::xml_ns_strip(xml2::read_xml(file.path(exdir, "xl", "workbook.xml")))
  rels <- xml2::xml_ns_strip(
    xml2::read_xml(file.path(exdir, "xl", "_rels", "workbook.xml.rels"))
  )
  rn <- xml2::xml_find_all(rels, ".//Relationship")
  rel <- stats::setNames(xml2::xml_attr(rn, "Target"), xml2::xml_attr(rn, "Id"))

  out <- list()
  for (sh in xml2::xml_find_all(wb, ".//sheet")) {
    nm <- xml2::xml_attr(sh, "name")
    ps <- file.path(exdir, "xl", sub("^/?", "", rel[[xml2::xml_attr(sh, "id")]]))
    doc <- xml2::read_xml(ps)
    out[[nm]] <- lapply(xml2::xml_find_all(doc, ".//d1:row", KASP_NS), function(row) {
      cells <- xml2::xml_find_all(row, "./d1:c", KASP_NS)
      col <- sub("[0-9]+$", "", xml2::xml_attr(cells, "r"))
      type <- xml2::xml_attr(cells, "t")
      v <- xml2::xml_text(xml2::xml_child(cells, "d1:v", KASP_NS))
      is_str <- !is.na(type) & type == "s"
      v[is_str] <- sst[as.integer(v[is_str]) + 1L]
      v[is.na(v)] <- ""
      stats::setNames(as.list(v), col)
    })
  }
  out
}

# one workbook row -> named character vector over A..Z, "" where absent
kasp_row <- function(row, cols = LETTERS[1:26]) {
  vapply(cols, function(k) if (is.null(row[[k]])) "" else row[[k]], character(1))
}

# the fetched copy if present, else the vendored copy (the CSVs must remain
# re-derivable from a clean checkout, where .tmp/lit does not exist)
kasp_input <- function(root, key, rel, vendored, what) {
  live <- fs::path(root, ".tmp", "lit", key, rel)
  if (file.exists(live)) return(live)
  vp <- fs::path(root, "inst", "extdata", vendored)
  if (file.exists(vp)) return(vp)
  stop(
    what, ": neither the fetched copy (", live, ") nor the vendored copy (",
    vp, ") is present. Re-fetch with the recipe in .tmp/sdd/report-kasp.md."
  )
}

# ---------------------------------------------------------------------------
#  a small text anchor for reading numbers out of a cross-validated full text
# ---------------------------------------------------------------------------
#  The regex engine in this toolchain has two failure modes that look like bad
#  data: (a) a space in a MULTI-WORD literal can fail to match ("Thirty out of"
#  does not match text that reads "Thirty-nine out of"), and (b) regmatches()
#  on a pattern that occurs more than once can return the LAST match.  Neither
#  is worth fighting: a bibliographic anchor is located by ONE word (or a
#  pattern with no spaces and no hyphenated suffix), and every number is cut
#  out of the text by position.
kasp_norm <- function(x) gsub("[[:space:]]+", " ", x)

kasp_find <- function(flat, anchor) {
  m <- regexpr(anchor, flat, perl = TRUE)
  if (m[1L] <= 0L) stop("anchor not found in the full text: \"", anchor, "\"")
  list(start = m[1L], anchor = anchor)
}

# the text at the anchor, for exactly `len` characters
kasp_site <- function(loc, flat, len = 60L) {
  substr(flat, loc$start, loc$start + len - 1L)
}

# the first digit run at or after the anchor (searched within `window` chars)
kasp_num_after <- function(loc, flat, window = 120L, what = loc$anchor) {
  s <- kasp_site(loc, flat, window)
  v <- strsplit(s, "")[[1L]]
  d <- grepl("[0-9]", v)
  if (!any(d)) stop("no digits within ", window, " chars of the anchor for ", what,
                    "; the text is: ", s)
  st <- which(d)[1L]
  en <- st
  while (en < length(v) &&
         (grepl("[0-9]", v[en + 1L]) || (v[en + 1L] == "," && en + 2L <= length(v) &&
                                         grepl("[0-9]", v[en + 2L])))) {
    en <- en + 1L
  }
  as.numeric(gsub(",", "", paste(v[st:en], collapse = "")))
}

# the numbers in a sentence fragment, in order (a comma joins two runs only
# when a digit follows it)
kasp_numbers <- function(s) {
  v <- strsplit(s, "")[[1L]]
  out <- numeric()
  i <- 1L
  while (i <= length(v)) {
    if (grepl("[0-9]", v[i])) {
      j <- i
      while (j < length(v) &&
             (grepl("[0-9]", v[j + 1L]) ||
              (v[j + 1L] == "," && j + 2L <= length(v) && grepl("[0-9]", v[j + 2L])))) {
        j <- j + 1L
      }
      out <- c(out, as.numeric(gsub(",", "", paste(v[i:j], collapse = ""))))
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  out
}

# the capture groups of a pattern, cut by the match attributes (never by
# regmatches, which can return a different match than the one that was located)
kasp_groups <- function(flat, pattern, what) {
  m <- regexpr(pattern, flat, perl = TRUE)
  if (m[1L] <= 0L) stop("pattern not found for ", what, ": ", pattern)
  cs <- attr(m, "capture.start")
  cl <- attr(m, "capture.length")
  if (is.null(cs) || !ncol(cs)) stop("the pattern for ", what, " has no capture group")
  vapply(seq_len(ncol(cs)), function(j) {
    substr(flat, cs[1L, j], cs[1L, j] + cl[1L, j] - 1L)
  }, character(1))
}

# every script ends here: the CSV on disk must equal the recomputed table.
# The comparison is on the WRITTEN text, not on a read-back: readr can guess a
# column type differently on the second pass and would report a phantom edit.
kasp_write_checked <- function(tbl, path) {
  path <- as.character(path)
  written <- utils::capture.output(
    readr::write_csv(tbl, stdout(), na = ""), type = "output"
  )
  if (file.exists(path)) {
    on_disk <- readLines(path, warn = FALSE)
    if (!identical(on_disk, written)) {
      nl <- min(length(on_disk), length(written))
      first_bad <- which(on_disk[seq_len(nl)] != written[seq_len(nl)])[1L]
      if (is.na(first_bad)) first_bad <- nl + 1L
      # Report the byte that differs, not just the line: two lines that look
      # identical in a 160-char window differ somewhere past it (an xlsx cell's
      # trailing space, or a hyphen the workbook stores as U+2212).
      a <- on_disk[first_bad]; b <- written[first_bad]
      j <- which(strsplit(a, "", fixed = TRUE)[[1]] !=
                 strsplit(b, "", fixed = TRUE)[[1]])[1]
      if (is.na(j)) j <- min(nchar(a), nchar(b)) + 1L
      stop("the CSV on disk disagrees with the recomputed table; first differing ",
           "line is ", first_bad, " at character ", j, ":\n  on disk: ",
           substr(a, max(1L, j - 40L), j + 60L), "\n  recomputed: ",
           substr(b, max(1L, j - 40L), j + 60L))
    }
  }
  readr::write_csv(tbl, path, na = "")
  invisible(path)
}
