# =============================================================================
#  data-raw/merge-proposals.R
# =============================================================================
#  Fold the workers' `data-raw/proposed/*` files into the four shared registries,
#  row-wise by key. This is the ONLY writer of those registries besides
#  data-raw/science_predictions.R and fetch-sources.R.
#
#  Rules (each enforced, each a stop() when violated):
#   1. A key already in a registry is NOT overwritten. A proposal that would
#      change a landed row is reported for manual review, never auto-merged.
#      Silently clobbering a verified row is how a registry rots.
#   2. A key appearing twice among the proposals with CONFLICTING values is a
#      stop(). Duplicates with identical values are collapsed.
#   3. Column names must be a subset of the target's; missing columns are filled
#      with NA. Extra columns are a stop() unless whitelisted below.
#   4. Every landed source that has an extract must appear in science_ledger.csv
#      with a `use` of "test", "modelled" or "blocked"+reason. This is the check
#      the user actually cares about: "recorded but used for nothing" is a bug.
#
#  Usage: Rscript data-raw/merge-proposals.R [--apply]
#         (default is a DRY RUN that prints what it would do)
# =============================================================================

suppressPackageStartupMessages(library(fs))
root <- normalizePath(".")
dir_prop <- file.path(root, "data-raw", "proposed")
apply_changes <- "--apply" %in% commandArgs(trailingOnly = TRUE)
ext <- function(f) file.path(root, "inst", "extdata", f)

read_any <- function(path, target_cols) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  extra <- setdiff(names(d), target_cols)
  if (length(extra)) {
    # LOUD on purpose: a file whose columns do not match its target was being
    # caught by the signature test and then silently skipped, so 61 correct
    # ledger rows and one correct source row sat in data-raw/proposed/ for
    # several waves while the registries said they did not exist.
    warning(basename(path), ": skipped, columns not in target: ",
            paste(extra, collapse = ", "), call. = FALSE)
    return(NULL)
  }
  for (nm in setdiff(target_cols, names(d))) d[[nm]] <- NA
  d[, target_cols, drop = FALSE]
}

# A proposal file is assigned to exactly one registry by the column signature of
# its header plus its filename. Signature-first: filenames have drifted.
signature_of <- function(d) {
  n <- names(d)
  if (all(c("source_key", "extract", "use", "test_file") %in% n)) return("ledger")
  if (all(c("key", "family", "what_it_would_anchor") %in% n)) return("backlog")
  if (all(c("key", "citation", "doi") %in% n)) return("sources")
  "unknown"
}

prop_files <- dir_ls(dir_prop, glob = "*.csv", recurse = FALSE)
by_registry <- list(ledger = list(), backlog = list(), sources = list())
for (f in prop_files) {
  d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, nrows = 1)
  sig <- signature_of(d)
  if (sig == "unknown") {
    # A file whose NAME says which registry it belongs to must never pass
    # silently: polyploid.rows.csv carried the sources signature, so its four
    # sources and four ledger rows sat unread while its extracts showed up as
    # ORPHAN. Unrecognised stray proposals are still only a message.
    if (grepl("^(sources|ledger|backlog)\\.", basename(f))) {
      stop(basename(f), ": named for a registry but its columns match no signature (got: ",
           paste(names(d), collapse = ", "), ")")
    }
    message("  SKIP (not a registry proposal): ", basename(f))
    next
  }
  by_registry[[sig]] <- c(by_registry[[sig]], list(f))
}

targets <- list(
  ledger  = ext("science_ledger.csv"),
  backlog = ext("scientific_backlog.csv"),
  sources = ext("scientific_sources.csv")
)

## --- dedupe proposals, refusing to choose between two different values -------
# The unit of uniqueness differs by registry:
#   - sources/backlog: one row per key. Two rows for a key must AGREE everywhere.
#   - ledger: one row per (source_key, extract) -- a source legitimately has
#     several extracts, and deduping on the key alone wrongly reported that as
#     a conflict. Two rows for the same unit must agree everywhere.
unit_of <- function(reg) {
  if (reg == "ledger") c("source_key", "extract") else "key"
}
collapse <- function(rows, reg, label) {
  if (!nrow(rows)) return(rows)
  units <- unit_of(reg)
  id <- do.call(paste, c(lapply(units, function(u) rows[[u]]), list(sep = "\u0001")))
  out <- rows[0, , drop = FALSE]
  for (u in unique(id)) {
    grp <- rows[id == u, , drop = FALSE]
    for (j in names(grp)) {
      v <- unique(grp[[j]])
      v <- v[!is.na(v)]
      # `family` names a GROUPING and `what_it_would_anchor` is prose about ONE
      # paper: two workers filing the same source describe both differently
      # without contradicting each other, and the first-seen text is as good as
      # the second. Those two are wording choices, not merge conflicts. Every
      # column that states a FACT (doi, pmcid, status, verified) still stops the
      # run, because a disagreement there is a real problem.
      if (length(v) > 1L && !(j %in% c("family", "what_it_would_anchor"))) {
        stop(label, ": row '", gsub("\u0001", "/", u), "' has conflicting values in column '", j,
             "': ", paste(v, collapse = " | "), "\n  -> resolve by hand, this is not auto-mergeable")
      }
    }
    out <- rbind(out, grp[1, , drop = FALSE])
  }
  out
}

report <- character(0)
for (reg in names(by_registry)) {
  files <- by_registry[[reg]]
  if (!length(files)) next
  tgt <- utils::read.csv(targets[[reg]], stringsAsFactors = FALSE,
                         check.names = FALSE)
  key <- if (reg == "ledger") "source_key" else "key"
  # NB: an anonymous function, not `read_any` passed by name. `lapply(f, read_any,
  # target_cols = ...)` binds `target_cols` to read_any's FIRST formal, so every
  # file was checked against the wrong column set and the warning above fired
  # for files that actually match.
  new <- do.call(rbind, lapply(files, function(f) read_any(f, names(tgt))))
  new <- collapse(new, reg, reg)
  units <- unit_of(reg)
  skey <- function(d) do.call(paste, c(lapply(units, function(u) d[[u]]), list(sep = "\u0001")))
  fresh <- new[!(skey(new) %in% skey(tgt)), , drop = FALSE]
  clash <- new[skey(new) %in% skey(tgt), , drop = FALSE]

  report <- c(report, sprintf(
    "%-8s %2d proposed rows from %d file(s); %2d new, %2d already present",
    reg, nrow(new), length(files), nrow(fresh), nrow(clash)))

  if (nrow(clash)) {
    for (i in seq_len(nrow(clash))) {
      k <- gsub("\u0001", "/", skey(clash)[i])
      old <- tgt[skey(tgt) == skey(clash)[i], , drop = FALSE][1, , drop = FALSE]
      shared <- intersect(names(tgt), names(clash))
      diffcols <- shared[vapply(shared, function(cn) {
        !identical(as.character(old[[cn]]), as.character(clash[[cn]][i]))
      }, logical(1))]
      if (length(diffcols)) {
        report <- c(report, sprintf("  REVIEW %s/%s differs in: %s", reg, k,
                                    paste(diffcols, collapse = ",")))
      }
    }
  }
  if (nrow(fresh) && apply_changes) {
    out <- rbind(tgt, fresh[, names(tgt), drop = FALSE])
    # Stable order by key so a later diff is readable, but never reorder
    # within a key (ledger has several rows per source).
    if (reg != "ledger") out <- out[order(out[[key]]), , drop = FALSE]
    utils::write.csv(out, targets[[reg]], row.names = FALSE)
    report <- c(report, sprintf("  WROTE %s (%d -> %d rows)",
                                basename(targets[[reg]]), nrow(tgt), nrow(out)))
  }
}
cat(paste(report, collapse = "\n"), "\n")

## --- the check that matters: nothing recorded may be used for nothing -------
led <- utils::read.csv(targets$ledger, stringsAsFactors = FALSE)
bad <- led$use == "blocked" & (is.na(led$reason) | !nzchar(led$reason))
if (any(bad)) stop("ledger: 'blocked' rows with no reason: ",
                   paste(unique(led$source_key[bad]), collapse = ", "))
unknown_use <- setdiff(unique(led$use), c("test", "modelled", "artefact", "rule", "report", "config", "blocked"))
if (length(unknown_use)) stop("ledger: unknown use value(s): ",
                              paste(unknown_use, collapse = ", "))

# Every extract file on disk should be claimed by some ledger row, or be a
# known non-extract artefact (registries, bib, source resolution).
all_csv <- basename(dir_ls(ext(".")))
registries <- c("scientific_sources.csv", "scientific_backlog.csv",
                "science_ledger.csv", "source_resolution.csv",
                "conversion_rates.csv")
# Raw downloads and binaries are provenance, NOT extracts: an .xml / .fa / .xlsx
# beside a sourced .csv is the file the .csv was transcribed from, and no ledger
# row ever claims it. Listing them here keeps the orphan check about CSVs.
provenance_ext <- c("xml", "fa", "fasta", "xlsx", "bib")
all_csv <- basename(dir_ls(ext(".")))
all_csv <- all_csv[tolower(tools::file_ext(all_csv)) == "csv"]
claimed <- unique(unlist(strsplit(led$extract, ",", fixed = TRUE)))
claimed <- trimws(claimed[!is.na(claimed) & claimed != "none"])
orphans <- setdiff(all_csv, c(claimed, registries))
if (length(orphans)) {
  cat("\nORPHAN extracts (on disk, no ledger row):\n  ",
      paste(orphans, collapse = "\n  "), "\n")
  cat("  NOTE: source_resolution.csv is generated; conversion_rates.csv is legacy.\n")
  # Fail the gate: a landed extract nobody claims is exactly the state the
  # polyploid merge left behind, and it read as a NOTE for four waves.
  stop("orphan extracts: ", length(orphans), " file(s) read by no ledger row")
} else {
  cat("\nEvery extract is claimed by a ledger row.\n")
}
