# =============================================================================
#  data-raw/sync-backlog-status.R
# =============================================================================
#  The backlog's `status` column drifted from reality: eight keys still said
#  "unfetched" after their extracts had landed in inst/extdata and their rows
#  had appeared in science_ledger.csv. A backlog that under-reports what exists
#  is as misleading as one that over-reports, and it is the file a reader uses
#  to decide what to work on next.
#
#  The LEDGER is the source of truth for "did something land": a row with a
#  real extract and use in {test, rule, config, report, artefact} means the
#  source was read. This script raises `status` to match, never lowers it:
#  a level is a record of how far the source was read, and only a human who
#  has re-read the paper may withdraw a claim.
# =============================================================================

root <- normalizePath(".")
bp <- file.path(root, "inst/extdata/scientific_backlog.csv")
lp <- file.path(root, "inst/extdata/science_ledger.csv")
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE,
                                  na.strings = c("NA", ""), colClasses = "character")
b <- rd(bp)
l <- rd(lp)

landed <- unique(l$source_key[l$use %in% c("test", "rule", "config", "report",
                                           "artefact")])
# a pure `blocked` ledger row is a read-but-nothing-landed source
blocked <- setdiff(unique(l$source_key[l$use == "blocked"]), landed)

# status rank: how far into the source we got
rank <- c(unfetched = 0L, "fulltext-read-no-assertion" = 1L,
          not_extractable = 2L, extracted = 3L)
old <- b$status
b$status[b$key %in% landed] <- "extracted"
b$status[b$key %in% blocked & rank[b$status] < 2L] <- "not_extractable"

up <- which(b$status != old)
if (length(up)) {
  for (i in up) cat(sprintf("  %-22s %-30s -> %s\n", b$key[i], old[i], b$status[i]))
} else {
  cat("backlog status already agrees with the ledger.\n")
}

# `verified` records that the source artefact was READ, for any status past
# "unfetched" -- the column means what was done, not what was produced.
b$verified <- b$status != "unfetched"
b$verified_on[b$verified & (is.na(b$verified_on) | !nzchar(b$verified_on))] <- "2026-10-21"

utils::write.csv(b, bp, row.names = FALSE, na = "NA")
cat(sprintf("backlog: %d rows; %s\n", nrow(b),
            paste(sprintf("%s=%d", names(table(b$status)), table(b$status)),
                  collapse = " ")))
still <- b$key[b$status == "unfetched"]
if (length(still)) {
  cat("still unfetched (", length(still), "): ", paste(still, collapse = ", "), "\n", sep = "")
}
