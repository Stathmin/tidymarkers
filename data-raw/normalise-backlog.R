# =============================================================================
#  data-raw/normalise-backlog.R
# =============================================================================
#  The backlog's `status` is a CLOSED vocabulary (owned by
#  tests/testthat/test-scientific-data.R). Workers wrote prose into it --
#  "fetched, not extractable in this environment", "fetched + extracted" -- which
#  made `verified` and `status` disagree and would have made any query of the
#  backlog return "no rows".
#
#  Mapping is explicit and printed; an unmapped label is a stop(), never a guess.
#  Idempotent.
# =============================================================================

root <- normalizePath(".")
path <- file.path(root, "inst", "extdata", "scientific_backlog.csv")
b <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
n0 <- nrow(b)

## status is the SHAPE of the outcome; the reason moves to `note` --------------
norm <- c(
  "unfetched" = "unfetched",
  "fetched" = "fetched",
  "extracted" = "extracted",
  "fulltext-read-no-assertion" = "fulltext-read-no-assertion",
  "not_extractable" = "not_extractable"
)
# prose -> canonical, with an optional note that must be preserved
prose <- c(
  "fetched + extracted" = "extracted",
  "fetched; landed as evidence_tier = weak only" = "extracted",
  "fetched, not extractable in this environment" = "not_extractable",
  "fetched, no usable table" = "not_extractable",
  "fetched; not extractable (no validation-rate table)" = "not_extractable",
  "fetched; not extractable (no validation-rate number)" = "not_extractable",
  "blocked: no OA full text" = "unfetched",
  "fetched-not-extractable: the row describes a different paper" =
    "not_extractable",
  "fetched; blocked on licence" = "fulltext-read-no-assertion"
)

changed <- character(0)
for (i in seq_len(nrow(b))) {
  st <- b$status[[i]]
  if (st %in% norm) next
  keep_why <- grepl("^fetched-not-extractable: ", st)
  if (keep_why) {
    b$note[[i]] <- paste(na.omit(c(b$note[[i]], sub("^fetched-not-extractable: ", "", st))),
                         collapse = " | ")
    b$status[[i]] <- "not_extractable"
    changed <- c(changed, sprintf("%s: prose -> not_extractable (reason kept in note)", b$key[[i]]))
    next
  }
  hit <- prose[st]
  if (is.na(hit)) {
    stop("unmapped backlog status for '", b$key[[i]], "': ", st,
         "\n  Add it to `prose` with the canonical level it really means.")
  }
  if (st != hit) {
    b$note[[i]] <- paste(na.omit(c(b$note[[i]], paste0("original status: ", st))),
                         collapse = " | ")
  }
  b$status[[i]] <- unname(hit)
  changed <- c(changed, sprintf("%s: %s -> %s", b$key[[i]], st, hit))
}

## verified must agree with status: it records that the artefact was READ ------
done <- b$status %in% c("extracted", "fulltext-read-no-assertion", "not_extractable")
fix_v <- as.logical(b$verified) != done
b$verified <- done
if (any(fix_v)) {
  changed <- c(changed, paste("corrected `verified` on:", paste(b$key[fix_v], collapse = ", ")))
}

## a row read but not turned into a CSV must say why somewhere ----------------
needs_why <- b$status %in% c("not_extractable", "fulltext-read-no-assertion") &
  (is.na(b$note) | !nzchar(b$note))
if (any(needs_why)) {
  changed <- c(changed, paste("WARNING: read-but-not-extracted rows with no note:",
                              paste(b$key[needs_why], collapse = ", ")))
}

if (length(changed)) {
  utils::write.csv(b, path, row.names = FALSE, na = "NA")
  cat(paste(changed, collapse = "\n"), "\n")
  cat(sprintf("scientific_backlog.csv: %d rows (unchanged count)\n", n0))
} else {
  cat("scientific_backlog.csv already normalised; nothing written.\n")
}
cat("status levels in use:", paste(sort(unique(b$status)), collapse = ", "), "\n")
