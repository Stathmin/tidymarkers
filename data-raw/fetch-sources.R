# Resolve every registry row (scientific_sources.csv + scientific_backlog.csv)
# against Crossref and Europe PMC, then:
#   1. write inst/extdata/tidymarkers.bib  (real metadata, not prose)
#   2. download OA full texts + supplementary files into .tmp/lit/<key>/
#   3. write inst/extdata/source_resolution.csv recording what happened
#
# Nothing here is scientific evidence. Resolution says a paper EXISTS; only a
# data-raw/<paper>.R that asserts a number against the fetched file may put
# that number into inst/extdata. Every parse failure is an ERROR, never a
# quietly-NA field: silent NAs are how a fetch pass "succeeds" with zero bytes.
#
# Usage: Rscript data-raw/fetch-sources.R [--no-download] [--only key1,key2]
#
# PARALLEL WORKERS: each lives in its own git worktree, where .tmp/ is absent
# (it is untracked). Run this script there with --only <its keys>. Resolution
# of a registry row is idempotent, but the CSV/BIB are shared outputs, so only
# the main science worktree writes them.
#
suppressPackageStartupMessages({
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || !length(a) ||
                             (length(a) == 1 && is.na(a))) b else a

options(timeout = 90)
args <- commandArgs(trailingOnly = TRUE)
DO_DOWNLOAD <- !("--no-download" %in% args)
ONLY <- character(0)
if (any(grepl("^--only", args))) {
  v <- sub("^--only=?", "", args[grepl("^--only", args)][1])
  if (!nzchar(v)) {
    i <- which(grepl("^--only", args))[1]
    v <- if (length(args) > i) args[i + 1] else ""
  }
  ONLY <- trimws(strsplit(v, ",")[[1]])
  ONLY <- ONLY[nzchar(ONLY)]
}

root <- normalizePath(".")
lit_dir <- file.path(root, ".tmp", "lit")
dir.create(lit_dir, showWarnings = FALSE, recursive = TRUE)

ua <- "tidymarkers-science-fetch/1.0 (mailto:librarian@tidymarkers.invalid)"
JSON_ACCEPT <- c("-H", shQuote("Accept: application/json"))

get_raw <- function(url, headers = character(0)) {
  tf <- tempfile()
  out <- suppressWarnings(system2(
    "curl", c("-sS", "-L", "-m", "90", "-A", shQuote(ua), headers,
              "-o", tf, "-w", "%{http_code}", url),
    stdout = TRUE, stderr = TRUE
  ))
  # `curl -w` prints the code last; stderr messages are interleaved above.
  code <- sub("^.*\n", "", paste(out, collapse = "\n"))
  ok <- identical(trimws(code), "200") && file.exists(tf) && file.size(tf) > 0
  list(ok = ok, code = code, path = tf,
       size = if (file.exists(tf)) file.size(tf) else 0)
}

get_json <- function(url, headers = JSON_ACCEPT) {
  r <- get_raw(url, headers)
  if (!r$ok) return(NULL)
  out <- try(fromJSON(r$path, simplifyVector = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

jstr <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  # jsonlite returns atoms for scalars (hitCount is 1L, pmcid is "PMC..").
  # `[[1]]` keeps the type, so an integer hitCount stayed integer and every
  # downstream `!is.na(h$pmcid) && nzchar(...)` test was against the wrong type.
  as.character(x[[1]])
}

# --- registry -------------------------------------------------------------

read_registry <- function() {
  src <- read.csv(file.path(root, "inst/extdata/scientific_sources.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  bl <- read.csv(file.path(root, "inst/extdata/scientific_backlog.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
  src <- src[, c("key", "citation", "venue", "year", "doi", "pmid", "pmcid")]
  src$origin <- "sources"
  bl2 <- data.frame(
    key = bl$key, citation = bl$identifier, venue = bl$family,
    year = NA_integer_, doi = NA_character_, pmid = NA_character_,
    pmcid = NA_character_, origin = "backlog", stringsAsFactors = FALSE
  )
  # A backlog identifier that IS a DOI, or carries one in parentheses.
  # NB: registrant "49.6" is real (von Ahsen), so allow single-segment suffixes.
  doi_in <- function(x) {
    m <- regmatches(x, regexpr("10\\.[0-9]{4}/[^ )\"']+", x))
    if (length(m)) m else NA_character_
  }
  bl2$doi <- vapply(bl2$citation, doi_in, character(1))
  # Drop the duplicate ye2001 row (it also lives in scientific_sources.csv).
  bl2 <- bl2[!(bl2$key %in% src$key), , drop = FALSE]
  rbind(src, bl2)
}

# --- resolve --------------------------------------------------------------
# Three routes, tried in order, with the route recorded per row:
#   crossref : DOI -> canonical metadata (title/journal/year/authors)
#   epmc_doi : DOI -> Europe PMC (pmid/pmcid/OA); DOI must be UNQUOTED, the
#              quoted form '10.1093/clinchem/..' returns 0 hits
#   epmc_search : free-text only when no DOI exists. Hits are NOT trusted
#              unless a year filter is present; an unconstrained search
#              returns 25 arbitrary papers.
resolve <- function(row) {
  # NB: EPMC nests hits as resultList$result. Reading `result` gives a list
  # whose sole element is itself `result`, so every field lookup is silently
  # NULL -- that produced a pass that "resolved" 10 papers to zero bytes.
  rec <- list(key = row$key, title = NA_character_, journal = NA_character_,
              year = NA_character_, authors = NA_character_,
              doi = if (is.na(row$doi)) NA_character_ else row$doi,
              pmid = NA_character_, pmcid = NA_character_, oa = FALSE,
              fulltext = NA_character_, source = NA_character_, note = NA_character_)

  if (!is.na(row$doi) && grepl("^https?://", row$doi)) {
    # Software record: there is no Crossref/EPMC entry to resolve against, and
    # reusing the companion paper's DOI made this row silently resolve to that
    # paper. Record the URL verbatim and stop.
    rec$doi <- row$doi
    rec$source <- "url"
    rec$note <- "software repository; provenance only, not evidence"
    return(rec)
  }

  if (!is.na(row$doi) && nzchar(row$doi)) {
    cr <- get_json(paste0("https://api.crossref.org/works/",
                          utils::URLencode(row$doi, reserved = TRUE)))
    if (!is.null(cr$message)) {
      m <- cr$message
      rec$title <- jstr(m$title)
      rec$journal <- jstr(m$`container-title`)
      yr <- jstr(m$`published-print`$`date-parts`[[1]][[1]])
      if (is.na(yr)) yr <- jstr(m$`published-online`$`date-parts`[[1]][[1]])
      rec$year <- as.character(yr)
      fam <- vapply(m$author %||% list(), function(a) jstr(a$family), character(1))
      rec$authors <- paste(fam[!is.na(fam)], collapse = " and ")
      rec$source <- "crossref"

      ep <- get_json(paste0(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=",
        utils::URLencode(paste0("DOI:", row$doi), reserved = TRUE),
        "&format=json&resultType=core&pageSize=1"))
      if (!is.null(ep) && (ep$hitCount %||% 0) > 0) {
        h <- ep$resultList$result[[1]]
        if (is.null(h) || !length(h)) h <- list()
        rec$pmid <- jstr(h$pmid)
        rec$pmcid <- jstr(h$pmcid)
        rec$oa <- identical(jstr(h$isOpenAccess), "Y")
      }
    } else {
      rec$note <- "crossref:no_record"
    }
  }

  if (is.na(rec$source) && !is.na(row$pmid) && nzchar(row$pmid)) {
    ep <- get_json(paste0(
      "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=",
      utils::URLencode(paste0("EXT_ID:", row$pmid), reserved = TRUE),
      "&format=json&resultType=core&pageSize=1"))
    if (!is.null(ep) && (ep$hitCount %||% 0) > 0) {
      h <- ep$resultList$result[[1]]
      if (is.null(h) || !length(h)) h <- list()
      rec$pmid <- jstr(h$pmid)
      rec$pmcid <- jstr(h$pmcid)
      rec$oa <- identical(jstr(h$isOpenAccess), "Y")
      rec$title <- jstr(h$title)
      rec$journal <- jstr(h$journalTitle)
      rec$year <- jstr(h$pubYear)
      rec$source <- "epmc"
    }
  }

  rec
}

# --- bibtex ---------------------------------------------------------------

bib_escape <- function(x) {
  if (is.na(x)) return("")
  x <- gsub("([&%$#_])", "\\\\\\1", x)
  gsub("[{}]", "", x)
}

to_bib <- function(recs) {
  chunks <- character(0)
  for (r in recs) {
    type <- if (identical(r$source, "url")) {
      "misc"
    } else if (!is.na(r$pmid) && nzchar(r$pmid)) "article" else "misc"
    fields <- c(
      title = r$title, author = r$authors, journal = r$journal,
      year = r$year, doi = r$doi, pmid = r$pmid, pmcid = r$pmcid
    )
    if (identical(r$source, "url")) {
      fields <- c(howpublished = r$doi, fields[names(fields) != "doi"])
    }
    fields <- fields[!is.na(fields) & nzchar(fields)]
    body <- paste0(
      "  ", names(fields), " = {", vapply(fields, bib_escape, character(1)), "},",
      collapse = "\n"
    )
    chunks <- c(chunks, sprintf("@%s{%s,\n%s\n}", type, r$key, sub(",$", "", body)))
  }
  chunks
}

# --- download -------------------------------------------------------------

fetch_assets <- function(rec) {
  kdir <- file.path(lit_dir, rec$key)
  dir.create(kdir, showWarnings = FALSE, recursive = TRUE)
  got <- character(0)
  if (!DO_DOWNLOAD || is.na(rec$pmcid) || !nzchar(rec$pmcid)) return(got)
  base <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/", rec$pmcid)

  r <- get_raw(paste0(base, "/fullTextXML"))
  if (r$ok && r$size > 2000) {
    file.copy(r$path, file.path(kdir, "fulltext.xml"), overwrite = TRUE)
    got <- c(got, "fulltext.xml")
  }

  # Endpoint is called supplementaryFiles but returns a ZIP.
  r2 <- get_raw(paste0(base, "/supplementaryFiles"))
  if (r2$ok && r2$size > 1000) {
    sp <- file.path(kdir, "supplements.zip")
    file.copy(r2$path, sp, overwrite = TRUE)
    unz <- try(utils::unzip(sp, exdir = file.path(kdir, "supplements")),
               silent = TRUE)
    if (!inherits(unz, "try-error")) got <- c(got, basename(unz))
  }
  got
}

# --- run ------------------------------------------------------------------

reg <- read_registry()
if (length(ONLY)) {
  missing <- setdiff(ONLY, reg$key)
  if (length(missing)) {
    stop("unknown key(s): ", paste(missing, collapse = ", "))
  }
  reg <- reg[match(ONLY, reg$key), , drop = FALSE]
}
message("registry rows: ", nrow(reg))

recs <- vector("list", nrow(reg))
assets <- character(nrow(reg))
for (i in seq_len(nrow(reg))) {
  rec <- resolve(reg[i, ])
  recs[[i]] <- rec
  assets[i] <- paste(fetch_assets(rec), collapse = ";")
  message(sprintf("  %-24s %-9s oa=%-5s %s", rec$key,
                  if (is.na(rec$source)) "UNRESOLVED" else rec$source,
                  isTRUE(rec$oa), assets[i]))
}

# A subset run must not clobber the shared registry outputs written by the
# full run; it only populates .tmp/lit/ for its own keys.
if (length(ONLY)) {
  message("\n--only run: assets downloaded ", sum(nzchar(assets)), "/", nrow(reg), "; manifest NOT rewritten")
  quit(save = "no")
}

res <- data.frame(
  key = vapply(recs, `[[`, character(1), "key"),
  origin = reg$origin,
  resolved = !vapply(recs, function(r) is.na(r$source), logical(1)),
  source = vapply(recs, function(r) as.character(r$source %||% NA), character(1)),
  title = vapply(recs, function(r) as.character(r$title %||% NA), character(1)),
  journal = vapply(recs, function(r) as.character(r$journal %||% NA), character(1)),
  year = vapply(recs, function(r) as.character(r$year %||% NA), character(1)),
  doi = vapply(recs, function(r) as.character(r$doi %||% NA), character(1)),
  pmid = vapply(recs, function(r) as.character(r$pmid %||% NA), character(1)),
  pmcid = vapply(recs, function(r) as.character(r$pmcid %||% NA), character(1)),
  oa = vapply(recs, function(r) isTRUE(r$oa), logical(1)),
  assets = assets,
  note = vapply(recs, function(r) as.character(r$note %||% NA), character(1)),
  stringsAsFactors = FALSE
)

utils::write.csv(res, file.path(root, "inst/extdata/source_resolution.csv"),
                 row.names = FALSE)

writeLines(c(
  "% Generated by data-raw/fetch-sources.R - do not edit by hand.",
  "% Metadata only: resolution is provenance, not evidence.",
  to_bib(recs)
), file.path(root, "inst/extdata/tidymarkers.bib"))

message("\nresolved: ", sum(res$resolved), "/", nrow(res),
        "  assets downloaded: ", sum(nzchar(res$assets)), "/", nrow(res))
if (any(!res$resolved)) {
  message("unresolved: ", paste(res$key[!res$resolved], collapse = ", "))
}
