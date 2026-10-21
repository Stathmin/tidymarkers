# =============================================================================
#  data-raw/normalise-sources.R
# =============================================================================
#  Bring the widened scientific_sources.csv back to its documented schema.
#  Workers added free-text labels and PMCID-style keys; both are provenance
#  noise, and one of them was a duplicate article wearing a different key.
#
#  This script is deliberately dumb and loud: every transformation it makes is
#  printed, and anything it cannot decide is a stop(). Run it after merging any
#  batch of proposed rows. Idempotent: a second run prints nothing and writes
#  nothing.
# =============================================================================

root <- normalizePath(".")
path <- file.path(root, "inst/extdata", "scientific_sources.csv")
s <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
before <- s

changed <- character(0)

## 1. Drop rows that duplicate another article under a PMCID-shaped key ------
# `PMC8621402_note` exists only to say "ijms2022 is the same article". That is a
# comment, not a source: its claim is empty and its what_it_anchors says there is
# nothing to extract. Keep the article row, delete the note.
dup_notes <- s$key[grepl("_note$", s$key) & !nzchar(s$claim)]
s <- s[!(s$key %in% dup_notes), , drop = FALSE]
if (length(dup_notes)) {
  changed <- c(changed, paste("dropped duplicate note row(s):", paste(dup_notes, collapse = ", ")))
}

## 2. Rename PMCID-shaped keys to the title's own identity --------------------
# A key is our handle, and `PMC13020337` tells a reader nothing. Derive a
# readable key from the first author + year, falling back to the PMCID when
# that is genuinely all we have. Never silently rename: print the map.
pmc_keys <- grep("^PMC[0-9]+$", s$key, value = TRUE)
if (length(pmc_keys)) {
  rename <- c(
    PMC13020337 = "snapshot_str2026",   # STR/cSNP SNaPshot sensitivity study
    PMC13644189 = "csnp14p2026",        # 14-plex cSNP SNaPshot panel
    PMC5704242  = "woerner_stutter2017" # STR stutter in MPS data
  )
  unmapped <- setdiff(pmc_keys, names(rename))
  if (length(unmapped)) {
    stop("PMCID-shaped key(s) with no mapping: ", paste(unmapped, collapse = ", "),
         "\n  Add them to `rename` after checking which real article they are.")
  }
  for (old in names(rename)) {
    new <- rename[[old]]
    if (new %in% s$key) stop("rename target already present: ", new)
    s$key[s$key == old] <- new
    changed <- c(changed, sprintf("renamed %s -> %s", old, new))
  }
}

## 3. `access` is a closed vocabulary; detail moves to `access_detail` --------
# The test that owns this vocabulary is tests/testthat/test-scientific-data.R.
# Levels mean only ONE thing: how far into the document we actually read.
access_map <- c(
  "fulltext" = "fulltext", "full text" = "fulltext",
  "open" = "fulltext", "full text (Tables 1-2)" = "fulltext",
  "full text (Tables 4, 6, 7)" = "fulltext",
  "full text + supplement" = "supplement",
  "full text + supplements (mmc1-mmc4)" = "supplement",
  "supplement" = "supplement",
  "abstract" = "abstract",
  "stub" = "stub", "blocked" = "stub",
  "transcribed" = "transcribed",
  "source" = "source",
  "closed" = "closed access", "closed access" = "closed access",
  "subscription" = "closed access"
)
bad <- setdiff(unique(s$access), names(access_map))
if (length(bad)) {
  stop("unmapped access label(s): ", paste(bad, collapse = ", "),
       "\n  Map each to one of: ", paste(unique(access_map), collapse = ", "))
}
s$access_detail <- ifelse(
  s$access %in% names(access_map) & access_map[s$access] != s$access,
  s$access, NA_character_)
s$access <- unname(access_map[s$access])
if (any(access_map[unique(before$access)] != unique(before$access))) {
  changed <- c(changed, "normalised `access`; original wording kept in `access_detail`")
}

## 4. A missing DOI must say WHY, not be blank -------------------------------
# The test requires every row to carry a DOI-shaped string. A row with no DOI
# has not been excluded from that rule by accident: it is either still being
# resolved or deliberately transcribed. Make it explicit.
nodoi <- is.na(s$doi) | !nzchar(s$doi)
if (any(nodoi)) {
  s$doi[nodoi] <- "DOI-UNRESOLVED"
  changed <- c(changed, paste("marked no-DOI rows:", paste(s$key[nodoi], collapse = ", ")))
}

## 5. `verified` is logical; prose moves to `access_detail` -------------------
# Lands here as TRUE/FALSE or as prose ('fulltext fetched', 'NOT fetched'). The
# test reads it as logical, so a prose value is a type error, not a label.
if (!is.logical(s$verified)) {
  raw <- as.character(s$verified)
  s$verified <- ifelse(toupper(raw) == "TRUE", TRUE,
                ifelse(toupper(raw) == "FALSE" | grepl("^NOT fetched", raw), FALSE,
                       grepl("fetched", raw)))
  changed <- c(changed, "normalised `verified` to logical; original prose kept in `access_detail`")
  s$access_detail <- ifelse(!is.na(raw) & !toupper(raw) %in% c("TRUE", "FALSE"),
                            paste(na.omit(c(s$access_detail, "verified:", raw)), collapse = " "),
                            s$access_detail)
}

## 6. A closed-access row must SAY it was not read -------------------------
# The test greps for "NOT READ". That is not decoration: a closed row with a
# neutral claim reads as if the paper had been consulted.
cl <- s$access == "closed access" & !grepl("NOT READ", s$claim)
if (any(cl)) {
  s$claim[cl] <- paste0("NOT READ: ", s$claim[cl])
  changed <- c(changed, paste("prefixed NOT READ on closed rows:",
                              paste(s$key[cl], collapse = ", ")))
}

## 7. Reorder columns to the documented order, with access_detail last -------
want <- c("key", "citation", "venue", "year", "doi", "pmid", "pmcid", "access",
          "verified", "verified_on", "what_it_anchors", "claim", "claim_kind",
          "access_detail")
s <- s[, want, drop = FALSE]

if (length(changed)) {
  utils::write.csv(s, path, row.names = FALSE, na = "NA")
  cat(paste(changed, collapse = "\n"), "\n")
  cat(sprintf("scientific_sources.csv: %d -> %d rows\n", nrow(before), nrow(s)))
} else {
  cat("scientific_sources.csv already normalised; nothing written.\n")
}

cat("access levels in use:", paste(sort(unique(s$access)), collapse = ", "), "\n")
