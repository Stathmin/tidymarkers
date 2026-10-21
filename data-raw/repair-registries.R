# =============================================================================
#  data-raw/repair-registries.R
# =============================================================================
#  A re-merge after a rename re-imported the OLD PMCID-shaped keys, and the
#  normaliser's `write.csv(na = "NA")` turned every real NA into the two-letter
#  STRING "NA" -- which is not NA to any reader, so the DOI guard saw 58
#  non-DOI values. Both are mechanical damage from running the pipeline twice,
#  not new data. This script repairs them and then re-asserts the invariants the
#  tests care about, so a silent half-repair cannot pass.
#
#  Order matters: repair keys, then NA-strings, then re-run the normalisers.
#  Idempotent, and it stops if a repair would destroy a non-recoverable value.
# =============================================================================

root <- normalizePath(".")
src <- file.path(root, "inst", "extdata", "scientific_sources.csv")
bl  <- file.path(root, "inst", "extdata", "scientific_backlog.csv")
led <- file.path(root, "inst", "extdata", "science_ledger.csv")
journal <- root

read_csv <- function(p) utils::read.csv(p, stringsAsFactors = FALSE,
                                        check.names = FALSE, na.strings = "NA")

## 1. "NA" literal -> real NA, in every character column ----------------------
na_strings <- function(d, label) {
  n <- 0L
  for (j in names(d)) {
    if (is.character(d[[j]])) {
      hit <- !is.na(d[[j]]) & d[[j]] == "NA"
      n <- n + sum(hit)
      d[[j]][hit] <- NA_character_
    }
  }
  if (n) cat(sprintf("  %-24s %d literal \"NA\" cells -> NA\n", label, n))
  d
}

## 2. PMCID keys that came back -> the readable key that replaced them -------
key_map <- c(PMC13020337 = "snapshot_str2026",
             PMC13644189 = "csnp14p2026",
             PMC5704242  = "woerner_stutter2017")

retire <- function(d, label, keycol = "key") {
  k <- d[[keycol]]
  for (old in names(key_map)) {
    new <- key_map[[old]]
    hit <- !is.na(k) & k == old
    if (!any(hit)) next
    if (new %in% d[[keycol]][!hit]) {
      # the renamed row already exists: the PMCID row is the duplicate, drop it
      d <- d[!hit, , drop = FALSE]
      cat(sprintf("  %-24s dropped duplicate %s (already present as %s)\n",
                  label, old, new))
    } else {
      d[[keycol]][hit] <- new
      cat(sprintf("  %-24s renamed %s -> %s\n", label, old, new))
    }
  }
  d
}

s <- read_csv(src)
before_keys <- s$key
# an EMPTY string is not NA to any reader downstream: normalise it first, or the
# DOI-shaped guard reports `belv2024 => ` with nothing after the arrow and the
# message cannot be acted on
for (cn in c("doi", "pmid", "pmcid", "year")) {
  if (cn %in% names(s)) {
    s[[cn]][!is.na(s[[cn]]) & as.character(s[[cn]]) == ""] <- NA
  }
}
# order matters: decide which ROWS exist BEFORE rewriting cells
s <- retire(s, "scientific_sources")
# A `_note` row is a comment about another row, not a source, and is identified
# by having NOTHING to anchor. Test that, never a single empty cell: the
# PolyMarker row keeps its substance in `what_it_anchors` with an empty `claim`,
# and a claim-emptiness test deleted a real source definition.
comment_row <- grepl("_note$", s$key) &
  (is.na(s$what_it_anchors) | !nzchar(s$what_it_anchors))
if (any(comment_row)) {
  cat(sprintf("  %-24s dropping %d comment row(s): %s\n", "scientific_sources",
              sum(comment_row), paste(s$key[comment_row], collapse = ", ")))
  s <- s[!comment_row, , drop = FALSE]
}
s <- na_strings(s, "scientific_sources")
# a row whose every cell is NA is a blank line wrangled by write.csv, not a source
all_na <- apply(is.na(s), 1, all)
if (any(all_na)) {
  cat(sprintf("  %-24s dropping %d all-NA row(s)\n", "scientific_sources", sum(all_na)))
  s <- s[!all_na, , drop = FALSE]
}
# A row must never Vanish between reading and writing. It happened: a source
# landed on one branch while `retire()` derived its drop set from a hard-coded
# key map, and the file came back one row shorter with no message.
no_lost_rows <- function(before, after, label) {
  if (!identical(sort(before), sort(after))) {
    stop(label, ": row keys changed on write; lost: ",
         paste(setdiff(before, after), collapse = ", "), " gained: ",
         paste(setdiff(after, before), collapse = ", "))
  }
  invisible(NULL)
}
no_lost_rows(before_keys, s$key, "scientific_sources")
utils::write.csv(s, src, row.names = FALSE, na = "NA")

b <- read_csv(bl)
b <- retire(b, "scientific_backlog")
b <- na_strings(b, "scientific_backlog")
b <- b[!apply(is.na(b), 1, all), , drop = FALSE]
utils::write.csv(b, bl, row.names = FALSE, na = "NA")

l <- read_csv(led)
l <- retire(l, "science_ledger", keycol = "source_key")
l <- na_strings(l, "science_ledger")
l <- l[!apply(is.na(l), 1, all), , drop = FALSE]
# `blocked` implies a reason and nothing else carries one
l$reason[l$use != "blocked"] <- NA_character_
# A row must never Vanish between reading and writing. It happened: a source
# landed on one branch while `retire()` derived its drop set from a hard-coded
# key map, and the file came back one row shorter with no message. Duplicate
# (source, extract) pairs are dropped for the same reason: the earlier row wins,
# because the later one was seeded by a classifier that did not read the file.
l$reason <- ifelse(is.na(l$reason) | !nzchar(l$reason), NA_character_, l$reason)
l <- l[!duplicated(paste(l$source_key, l$extract)), , drop = FALSE]
utils::write.csv(l, led, row.names = FALSE, na = "NA")

## 3. Re-assert the invariants the tests own ---------------------------------
s <- read_csv(src)
b <- read_csv(bl)
l <- read_csv(led)

# `write.csv(na="NA")` writes the two characters NA; read_csv with
# na.strings="NA" turns them back. But a cell that was genuinely EMPTY comes
# back as "" and no guard can say what it meant. Fill the documented sentinel.
for (nm in c("s", "b")) {
  d <- get(nm)
  if ("doi" %in% names(d)) {
    blank <- is.na(d$doi) | !nzchar(d$doi)
    if (any(blank)) {
      d$doi[blank] <- "DOI-UNRESOLVED"
      assign(nm, d)
      cat(sprintf("  %-24s marked %d blank doi -> DOI-UNRESOLVED\n", nm, sum(blank)))
    }
  }
}
utils::write.csv(s, src, row.names = FALSE, na = "NA")
utils::write.csv(b, bl, row.names = FALSE, na = "NA")

stopifnot(!anyDuplicated(s$key), !anyDuplicated(b$key),
          !anyDuplicated(paste(l$source_key, l$extract)))
stopifnot(all(!is.na(s$key) & nzchar(s$key)))

# a DOI is DOI-shaped, or the explicit sentinel; a sentinel row needs a PMCID
unres <- s$doi == "DOI-UNRESOLVED"
bad_doi <- !unres & !grepl("^10\\.[0-9]{4,9}/", s$doi)
if (any(bad_doi)) {
  stop("DOI not DOI-shaped and not the sentinel:\n",
       paste("  ", s$key[bad_doi], "=>", s$doi[bad_doi], collapse = "\n"))
}
# A row with no DOI must still be reachable: by PMCID, or -- for a documented
# negative search -- by the query log that proves the search happened.
neg_search <- !is.na(s$access) & s$access == "not found"
reachable <- !is.na(s$pmcid) & nzchar(s$pmcid) |
  (neg_search & !is.na(s$access_detail) & nzchar(s$access_detail))
unreach <- unres & !reachable
if (any(unreach)) {
  stop("DOI-UNRESOLVED row(s) with no PMCID and no search log: ",
       paste(s$key[unreach], collapse = ", "))
}

# required-text columns are non-empty, except a row that is a software record
# or that says in its claim why it has no venue
req <- c("key", "citation", "venue", "year", "access", "what_it_anchors", "claim")
for (cn in req) {
  if (!cn %in% names(s)) stop("column missing from scientific_sources: ", cn)
  empty <- is.na(s[[cn]]) | !nzchar(as.character(s[[cn]]))
  if (any(empty)) {
    stop("empty required column '", cn, "' for: ",
         paste(s$key[empty], collapse = ", "))
  }
}
stopifnot(is.logical(s$verified) || all(s$verified %in% c(TRUE, FALSE)))

cat(sprintf("repaired. sources=%d backlog=%d ledger=%d\n",
            nrow(s), nrow(b), nrow(l)))
cat("access:", paste(sort(unique(s$access)), collapse = ", "), "\n")
cat("status:", paste(sort(unique(b$status)), collapse = ", "), "\n")
cat("use   :", paste(sort(unique(l$use)), collapse = ", "), "\n")
