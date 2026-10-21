# =============================================================================
#  data-raw/fix-typing-citations.R
# =============================================================================
#  Three typing sources were merged with the placeholder citation
#  "Anonymised for proposal: <title>" and an EMPTY doi/pmid/year. A registry row
#  whose DOI is blank is a promise, not a provenance record, and the placeholder
#  is exactly the kind of text that survives review because it looks deliberate.
#
#  This script rebuilds those three rows from the fetched JATS XML:
#  first author + year becomes the citation, the DOI is read from
#  <article-id pub-id-type="doi">, the PMID from <article-id pub-id-type="pmid">.
#
#  It stops when the fetched file is absent: run it in a checkout with .tmp/lit.
#  Idempotent: a second run finds nothing to fix.
# =============================================================================

suppressPackageStartupMessages(library(xml2))
root <- normalizePath(".")
lit <- file.path(root, ".tmp", "lit")
path <- file.path(root, "inst/extdata", "scientific_sources.csv")
s <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

# registry key -> fetched directory name (the two disagree because keys were
# renamed after the fetch, and the directory names are the PMCIDs they came from)
map <- c(snapshot_str2026 = "PMC13020337",
         csnp14p2026      = "PMC13644189",
         woerner_stutter2017 = "PMC5704242")

txt <- function(node, xp) {
  v <- xml_text(xml_find_first(node, xp))
  if (length(v) && !is.na(v)) trimws(v) else NA_character_
}
authors <- function(x) {
  n <- xml_find_all(x, "//contrib-group//contrib[not(@contrib-type='editor')]//name")
  if (!length(n)) n <- xml_find_all(x, "//contrib//name")
  sur <- xml_text(xml_find_first(n, ".//surname"))
  data.frame(surname = sur, given = xml_text(xml_find_first(n, ".//given-names")),
             stringsAsFactors = FALSE)
}

fixed <- character(0)
for (key in names(map)) {
  i <- which(s$key == key)
  if (!length(i)) next
  needs <- is.na(s$doi[i]) || !nzchar(s$doi[i]) ||
    grepl("^Anonymised for proposal", s$citation[i])
  if (!needs) next

  f <- file.path(lit, map[[key]], "fulltext.xml")
  if (!file.exists(f)) {
    stop("fetched XML missing for ", key, " (", f, "); cannot rebuild its citation")
  }
  x <- read_xml(f)
  # JATS: <front><article-meta><article-id pub-id-type="doi">..</article-id>
  doi <- txt(x, "//article-id[@pub-id-type='doi']")
  pmid <- txt(x, "//article-id[@pub-id-type='pmid']")
  title <- txt(x, "//article-title")
  journal <- txt(x, "//journal-title")
  year <- txt(x, "//pub-date/year")
  # first author. JATS spells the author block at least three ways seen in these
  # three files: <contrib-group content-type="author"><contrib><name>,
  # <contrib-group><contrib contrib-type="author"><name>, and variants with
  # per-contrib content-type. The XPath selection missed two of three and
  # silently produced NA, so select by STRUCTURE and stop if empty.
  au <- authors(x)
  sur <- if (nrow(au)) au$surname[[1]] else NA_character_
  giv <- if (nrow(au)) au$given[[1]] else NA_character_
  initials <- if (!is.na(giv)) {
    paste0(substr(strsplit(giv, " ")[[1]], 1, 1), collapse = "")
  } else ""
  if (anyNA(c(doi, title, year, sur)) || !nzchar(doi)) {
    stop(key, ": fetched XML lacks a field needed to rebuild the citation ",
         "(doi=", doi, ", title=", is.na(title), ", year=", year, ", surname=", sur, ")")
  }
  s$citation[i] <- paste0(sur, " ", initials, ", et al. ", title, ". ", journal,
                          " ", year, ".")
  s$doi[i] <- doi
  if (!is.na(pmid)) s$pmid[i] <- pmid
  if (!is.na(journal)) s$venue[i] <- journal
  s$year[i] <- as.integer(year)
  s$access_detail[i] <- paste(na.omit(c(s$access_detail[i],
                                        "citation rebuilt from fetched JATS")),
                              collapse = "; ")
  fixed <- c(fixed, sprintf("%s <- %s (%s)", key, doi, journal))
}

if (length(fixed)) {
  utils::write.csv(s, path, row.names = FALSE, na = "NA")
  cat(paste(fixed, collapse = "\n"), "\n")
} else {
  cat("no typing citations needed rebuilding.\n")
}
