# =============================================================================
#  data-raw/wave5-ledger-classify.R
# =============================================================================
#  Seeds a ledger row for EVERY registered source that has none, by asking what
#  the corpus can actually say about it. The class is derived from the source's
#  own registry row (access level, claim_kind, whether what_it_anchors says
#  NOT READ / BLOCKED) and from whether any extract on disk cites the source's
#  own DOI. Each row states what the source can and cannot support.
#
#    test    an extract on disk cites this source's DOI and a test reads it
#    config  the extract exists as design INPUT (primers, assay tables)
#    report  an extract exists but carries only what the paper ASSERTS
#    blocked the corpus holds nothing from this source, and the row says why
#
#  Level semantics: tests/testthat/test-science-ledger.R. A `blocked` row must
#  carry a reason and no test file; no other level may carry a reason.
#
#  Additive and idempotent: an existing row is never rewritten. Re-seeding a
#  source's row set happens only when the source has NO row at all.
#
#  Run:  Rscript data-raw/wave5-ledger-classify.R
# =============================================================================

`%or%` <- function(x, y) if (length(x) != 1L || is.na(x)) y else x

a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(a)) dirname(dirname(normalizePath(sub("^--file=", "", a[[1L]])))) else getwd()
setwd(root)
ext_dir <- file.path(root, "inst", "extdata")
led_path <- file.path(ext_dir, "science_ledger.csv")
test_dir <- file.path(root, "tests", "testthat")

s <- utils::read.csv(file.path(ext_dir, "scientific_sources.csv"),
                     stringsAsFactors = FALSE, na.strings = c("NA", ""))
l <- utils::read.csv(led_path, stringsAsFactors = FALSE, na.strings = c("NA", ""),
                     colClasses = "character")
l$reason <- ifelse(is.na(l$reason) | !nzchar(l$reason), NA_character_, l$reason)
files <- list.files(ext_dir, pattern = "\\.csv$")
tests <- list.files(test_dir, pattern = "^test-.*\\.R$")

# every extract that cites a DOI, and every test that names that extract
ext_doi <- lapply(files, function(f) {
  d <- utils::read.csv(file.path(ext_dir, f), stringsAsFactors = FALSE,
                       check.names = FALSE, colClasses = "character")
  if (!"doi" %in% names(d)) return(character(0))
  unique(stats::na.omit(d$doi))
})
names(ext_doi) <- files
test_src <- vapply(tests, function(t) {
  paste(readLines(file.path(test_dir, t), warn = FALSE), collapse = "\n")
}, character(1))
reader <- function(f) {
  hit <- names(test_src)[grepl(f, test_src, fixed = TRUE)]
  hit %or% NA_character_
}
# design-input files never carry a verdict: primers, assay tables, designs
is_config <- function(f) grepl("_primers?\\.csv$|_assays?\\.csv$|_design", f)

mk <- function(k, f, use, test, verdict, reason = NA_character_) {
  data.frame(source_key = k, extract = f, use = use, test_file = test,
             verdict_moved = verdict, reason = reason, stringsAsFactors = FALSE)
}

seed <- list()
for (i in seq_len(nrow(s))) {
  k <- s$key[[i]]
  if (k %in% l$source_key) next
  doi <- s$doi[[i]]
  mine <- if (!is.na(doi) && grepl("^10\\.", doi)) {
    names(ext_doi)[vapply(ext_doi, function(d) doi %in% d, logical(1))]
  } else character(0)
  mine <- setdiff(mine, l$extract)
  anchor <- s$what_it_anchors[[i]] %or% ""
  claim_kind <- s$claim_kind[[i]] %or% ""

  has_test <- vapply(mine, function(f) !is.na(reader(f)), logical(1))
  used <- mine[has_test]
  cfg <- mine[!has_test & vapply(mine, is_config, logical(1))]
  rpt <- setdiff(mine, c(used, cfg))

  for (f in used) {
    seed[[length(seed) + 1L]] <- mk(k, f, "test", reader(f),
      paste0("what ", f, " measures, asserted by ", reader(f)))
  }
  for (f in cfg) {
    seed[[length(seed) + 1L]] <- mk(k, f, "config", reader(f),
      paste0("design INPUT (", f, "): no verdict is asserted from it alone"))
  }
  for (f in rpt) {
    seed[[length(seed) + 1L]] <- mk(k, f, "report", reader(f),
      paste0("REPORTED content of ", f, ": the paper's own claim, kept apart from any measurement"),
      paste0("a reported-claim table read by no test, so ", f,
             " can support a citation but not a verdict"))
  }
  if (!length(mine)) {
    why <- if (grepl("^NOT READ|^BLOCKED", anchor) || claim_kind == "blocked") {
      paste0("nothing was landed from this source: ", sub("^NOT READ: *", "", anchor))
    } else if (s$access[[i]] %in% c("closed access", "stub", "abstract", "not found")) {
      paste0("no full text (", s$access[[i]], "): ", s$access_detail[[i]] %or% anchor,
             ", so no number of this paper could be checked and none was landed")
    } else {
      paste0("the full text was read, but no outcome in it can be asserted by an extract, so nothing was landed: ",
             substr(anchor, 1L, 200L))
    }
    seed[[length(seed) + 1L]] <- mk(k, "none", "blocked", NA_character_, "nothing moved", why)
  }
}

## --- files no ledger row names ---------------------------------------------
# Some belong to no single paper (shared priors, engine probes); most are a
# second table of a source already in the ledger. Either way a row names the
# test that reads the file, so the orphan check in merge-proposals.R stays
# empty for a real reason rather than being silenced.
claimed <- unique(unlist(strsplit(l$extract, ",", fixed = TRUE)))
claimed <- trimws(claimed[!is.na(claimed)])
left <- setdiff(setdiff(files, claimed), c("scientific_sources.csv",
                                           "scientific_backlog.csv",
                                           "science_ledger.csv",
                                           "source_resolution.csv",
                                           "conversion_rates.csv"))
for (f in left) {
  stem <- sub("\\.csv$", "", f)
  cands <- unique(c(s$key[s$key == stem],
                    s$key[s$key == sub("_.*$", "", f)],
                    s$key[s$key == sub("_", "", f, fixed = TRUE)]))
  cands <- cands[nzchar(cands)]
  r <- reader(f)
  use <- if (is.na(r)) "blocked" else "test"
  seed[[length(seed) + 1L]] <- mk(
    if (length(cands)) cands[[1L]] else paste0("shared_", stem), f, use, r,
    if (is.na(r)) "nothing moved" else paste0("what ", f, " measures, asserted by ", r),
    if (is.na(r)) paste0("no test reads ", f, ", and no single paper owns it: it is on disk but moves no verdict") else NA_character_)
}

new <- do.call(rbind, seed)
key <- function(d) paste(d$source_key, d$extract)
fresh <- new[!key(new) %in% key(l), ]
l <- rbind(l, fresh[, names(l), drop = FALSE])
utils::write.csv(l, led_path, row.names = FALSE, na = "")
message("wave5 ledger classify: ", nrow(fresh), " row(s) added; ", nrow(l), " total")
