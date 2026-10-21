# =============================================================================
#  data-raw/fix-backlog-doi.R
# =============================================================================
#  The backlog's `doi` column came out of a proposal merge SHIFTED: artemis2024
#  held Ayyadevara's DOI, ayyadevara2000 held Badolo's, and so on. A shifted
#  join is worse than a blank -- every row looks populated and every row is
#  wrong -- so this script does not "repair" the column, it REBUILDS it from
#  the one field a worker writes by hand: `identifier`.
#
#  Order of authority for each key:
#    1. the DOI that occurs verbatim in that row's own `identifier` prose,
#    2. the registry's DOI for the same key, when the prose names none,
#    3. the explicit sentinel DOI-UNRESOLVED.
#  Then it CHECKS the result: a DOI that appears in no row's prose and in no
#  registry row is reported, because that is how a shift looks.
#
#  NB: no `ifelse()` anywhere. In this R build,
#  `ifelse(lengths(m) == 1L, unlist(m), NA_character_)` returns a vector of the
#  UNLISTED length, not the condition's length -- two different list lengths
#  silently recycle. Build vectors element-wise instead.
# =============================================================================

root <- normalizePath(".")
bp <- file.path(root, "inst/extdata/scientific_backlog.csv")
sp <- file.path(root, "inst/extdata/scientific_sources.csv")
read_chr <- function(p) {
  utils::read.csv(p, stringsAsFactors = FALSE, na.strings = c("NA", ""),
                  colClasses = "character")
}
b <- read_chr(bp)
s <- read_chr(sp)

doi_pat <- "10\\.[0-9]{4,9}/[^ )\";,]+"

first_doi <- function(txt) {
  if (is.na(txt) || !nzchar(txt)) return(NA_character_)
  m <- regmatches(txt, regexpr(doi_pat, txt, perl = TRUE))
  if (length(m) == 1L && nzchar(m)) m else NA_character_
}

prose <- vapply(b$identifier, first_doi, character(1), USE.NAMES = FALSE)
reg <- s$doi[match(b$key, s$key)]
bad_reg <- is.na(reg) | !nzchar(reg) | !grepl(paste0("^", doi_pat), reg)
reg[bad_reg] <- NA_character_

new <- prose
need <- is.na(new)
new[need] <- reg[need]
new[is.na(new)] <- "DOI-UNRESOLVED"
# A DOI someone wrote by hand outranks prose: the first version preferred prose
# outright and re-NULLed five DOIs that had been looked up on purpose, which is
# the same class of silent loss this script exists to repair.
keep <- !is.na(b$doi) & nzchar(b$doi) & b$doi != "DOI-UNRESOLVED" &
  grepl(paste0("^", doi_pat), b$doi)
new[keep] <- b$doi[keep]
stopifnot(length(new) == nrow(b))

changed <- which(is.na(b$doi) | b$doi != new)
if (length(changed)) {
  for (i in changed) cat(sprintf("  %-22s %-40s -> %s\n", b$key[i],
                                 if (is.na(b$doi[i])) "(blank)" else b$doi[i], new[i]))
}
b$doi <- new

# the shift test: a DOI that belongs to NO row here and to NO registry row
known <- unique(c(prose, s$doi))
orphan <- b$doi != "DOI-UNRESOLVED" & !(b$doi %in% known)
if (any(orphan)) {
  cat("WARNING: DOI named nowhere in this file or the registry:\n")
  for (i in which(orphan)) cat(sprintf("  %-22s %s\n", b$key[i], b$doi[i]))
}

utils::write.csv(b, bp, row.names = FALSE, na = "NA")
cat(sprintf("backlog: %d rows; %d resolved from prose/registry, %d sentinel, %d row(s) changed\n",
            nrow(b), sum(b$doi != "DOI-UNRESOLVED"), sum(b$doi == "DOI-UNRESOLVED"),
            length(changed)))
