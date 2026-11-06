#!/usr/bin/env Rscript
# =============================================================================
#  dev/relabel-bonsai.R  --  rename misleading graphify communities, strictly
# =============================================================================
#  `graphify update .` labels each community with its hub node's name.  That is
#  fine for a single-file community and misleading for a multi-file one, and its
#  built-in prompt echoes member symbols and document headings.  This script asks
#  a local chat model for a *purpose* name with a strict prompt, and records the
#  pass in `graphify-out/.graphify_labels.origin.json` so the work is traceable.
#
#  Targets by status (same classifier as dev/graphify-health.R):
#    default      missing + duplicate + echo-stale
#    --all-echo   also echo-honest (a single-file file-echo is honest, not wrong)
#    --dry-run    list the targets and exit (the self-check)
#    --limit=N    cap the batch (smoke-test the model before a long run)
#    --cids=1,2,3 repair specific communities the status classifier cannot see
#
#  The model is OpenAI-compatible and local:
#    OPENAI_BASE_URL (default http://127.0.0.1:8000/v1)
#    OPENAI_MODEL    (default Ternary-Bonsai-2-27B-tight)
#    OPENAI_API_KEY  (default "local")
#
#  Writes `.graphify_labels.json` in place at the end only, and never touches
#  `.graphify_labels.json.sig` (membership-based; editing it renames hubs).
#
#  Run:  Rscript dev/relabel-bonsai.R [--all-echo] [--dry-run] [--limit=N]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

GRAPH <- "graphify-out/graph.json"
LABELS <- "graphify-out/.graphify_labels.json"
ORIGIN <- "graphify-out/.graphify_labels.origin.json"

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
all_echo <- "--all-echo" %in% args
dry_run <- "--dry-run" %in% args
limit_arg <- grep("^--limit=", args, value = TRUE)
limit <- if (length(limit_arg)) as.integer(sub("^--limit=", "", limit_arg[[1L]])) else Inf
cids_arg <- grep("^--cids=", args, value = TRUE)

base <- Sys.getenv("OPENAI_BASE_URL", "http://127.0.0.1:8000/v1")
model <- Sys.getenv("OPENAI_MODEL", "Ternary-Bonsai-2-27B-tight")
key <- Sys.getenv("OPENAI_API_KEY", "local")

if (!file.exists(GRAPH) || !file.exists(LABELS)) {
  message("skip (no ", GRAPH, " or ", LABELS, " - run `graphify update .` first)")
  quit(status = 0L)
}
for (pkg in c("jsonlite", "curl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("skip (", pkg, " not installed)")
    quit(status = 0L)
  }
}

g <- jsonlite::fromJSON(GRAPH, simplifyVector = FALSE)
labels <- jsonlite::fromJSON(LABELS, simplifyVector = TRUE)

nodes <- g$nodes
cid_of <- vapply(nodes, function(n) as.character(n$community %||% NA), character(1))
lab_of <- vapply(nodes, function(n) n$label %||% n$id %||% "", character(1))
file_of <- vapply(nodes, function(n) n$source_file %||% "", character(1))
keep <- !is.na(cid_of)
members <- split(lab_of[keep], cid_of[keep])
files <- lapply(split(file_of[keep], cid_of[keep]), function(x) sort(unique(x[nzchar(x)])))

# `graphify update .` renumbers communities and drops the label entry for each
# one it does not re-attach, so the file can cover fewer communities than the
# graph has.  A missing entry *is* the `missing` status (dev/graphify-health.R
# reads it that way); restore it as "" so this pass can name it.
comm_ids <- sort(unique(cid_of[keep]))
labels <- stats::setNames(
  vapply(comm_ids, function(c) labels[[c]] %||% "", character(1)),
  comm_ids
)

# same echo test as dev/graphify-health.R -- vectorised with `|`, not `||`
is_file_echo <- function(x) {
  grepl(" + ", x, fixed = TRUE) |
    grepl("(^|\\.)(R|r|md|py|yml|yaml|js|json|csv|tsv)$", x)
}

name_v <- unname(labels)
dup_labels <- names(which(table(labels[nzchar(labels)]) > 1L))
nfiles <- vapply(names(labels), function(c) length(files[[c]] %||% character()), integer(1))
echo <- is_file_echo(name_v)
status <- ifelse(
  !nzchar(name_v), "missing",
  ifelse(
    name_v %in% dup_labels, "duplicate",
    ifelse(echo & nfiles >= 2L, "echo-stale",
      ifelse(echo, "echo-honest", "named")
    )
  )
)
names(status) <- names(labels)

wanted <- if (all_echo) status != "named" else status %in% c("missing", "duplicate", "echo-stale")
targets <- sort(as.integer(names(status)[wanted]))
if (length(cids_arg)) {
  # repair lever for a name the classifier now calls `named` but that is wrong
  targets <- sort(intersect(
    as.integer(strsplit(sub("^--cids=", "", cids_arg[[1L]]), ",")[[1L]]),
    as.integer(names(status))
  ))
}
if (is.finite(limit)) targets <- head(targets, limit)

if (!length(targets)) {
  message(
    "no communities need a label (",
    paste(names(table(status)), as.integer(table(status)), sep = "=", collapse = ", "), ")"
  )
  quit(status = 0L)
}

ctx <- function(cid) {
  syms <- members[[as.character(cid)]] %||% character()
  syms <- unique(syms[!is_file_echo(syms) & !grepl("\\(\\)$", syms)])
  list(
    files = head(files[[as.character(cid)]] %||% character(), 14),
    symbols = head(syms, 24)
  )
}

if (dry_run) {
  cat("would relabel ", length(targets), " communit", if (length(targets) == 1L) "y" else "ies", ":\n", sep = "")
  for (cid in targets) {
    cat(sprintf(
      "  %3d  %-10s  %-34s  files=%d\n",
      cid, status[[as.character(cid)]],
      paste0("`", labels[[as.character(cid)]] %||% "", "`"), nfiles[[as.character(cid)]]
    ))
  }
  quit(status = 0L)
}

SYSTEM <- paste0(
  "You name code communities in a plant-marker-design R package. Each community ",
  "is a set of files and symbols that belong together. Give each one a short ",
  "purpose name.\n\n",
  "Rules for every name:\n",
  "- Title Case, 2 to 5 words.\n",
  "- Name the shared job of the cluster, not one member and not a file.\n",
  "- Never include a file name, a symbol name, a document heading, or the id.\n",
  "- Never use '::', '\u00a7', '/', or a file extension.\n",
  "- Return ONLY a JSON object mapping each id (string) to its name.\n\n",
  "Example input:\n",
  '{"31": {"files": ["read_panel.R", "write_panel.R"], ',
  '"symbols": ["read_panel", "write_panel", "panel_cache"]}}\n',
  "Example output:\n",
  '{"31": "Panel Storage Cache"}'
)
DISTINCT <- paste(
  "\nNames must be distinct from each other; use the symbol list to tell them",
  "apart, even when communities share a file."
)

salvage <- function(text) {
  text <- gsub("```(?:json)?", "", text)
  text <- gsub('^["[:space:]`]+|["[:space:]`]+$', "", text)
  obj <- tryCatch(jsonlite::fromJSON(text, simplifyVector = FALSE), error = function(e) NULL)
  if (is.list(obj) && length(obj)) {
    return(stats::setNames(
      vapply(obj, function(x) as.character(x)[[1L]], character(1)),
      names(obj)
    ))
  }
  hits <- regmatches(text, gregexpr('"[0-9]+"[[:space:]]*:[[:space:]]*"[^"]+"', text, perl = TRUE))[[1L]]
  if (!length(hits)) {
    return(stats::setNames(character(), character()))
  }
  caps <- regmatches(hits, regexec('"([0-9]+)"[[:space:]]*:[[:space:]]*"([^"]+)"', hits))
  stats::setNames(
    vapply(caps, `[[`, character(1), 3L),
    vapply(caps, `[[`, character(1), 2L)
  )
}

ask <- function(cids, distinct = FALSE, avoid = character()) {
  payload <- stats::setNames(lapply(cids, ctx), as.character(cids))
  sys <- if (distinct) paste0(SYSTEM, DISTINCT) else SYSTEM
  avoid <- unique(avoid[nzchar(avoid)])
  if (length(avoid)) {
    sys <- paste0(
      sys, "\nDo not reuse these names, they are already taken: ",
      paste(avoid, collapse = "; "), "."
    )
  }
  body <- list(
    model = model,
    messages = list(
      list(role = "system", content = sys),
      list(role = "user", content = jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"))
    ),
    temperature = 0,
    max_tokens = 64L + 24L * length(cids)
  )
  h <- curl::new_handle()
  curl::handle_setopt(
    h,
    post = TRUE,
    postfields = jsonlite::toJSON(body, auto_unbox = TRUE, null = "null"),
    timeout = 300L,
    httpheader = c(
      "Content-Type: application/json",
      paste0("Authorization: Bearer ", key)
    )
  )
  res <- curl::curl_fetch_memory(paste0(sub("/+$", "", base), "/chat/completions"), handle = h)
  if (res$status_code != 200L) {
    stop("HTTP ", res$status_code, ": ", substr(rawToChar(res$content), 1L, 200L))
  }
  parsed <- jsonlite::fromJSON(rawToChar(res$content), simplifyVector = FALSE)
  salvage(parsed$choices[[1L]]$message$content %||% "")
}

clean <- function(x) trimws(gsub('^["[:space:]]+|["[:space:]]+$', "", x))

bad <- function(cid, name, c) {
  if (is.null(name) || is.na(name) || !is.character(name)) {
    return("not a string")
  }
  s <- clean(name)
  if (!nzchar(s)) {
    return("empty")
  }
  if (grepl("::", s, fixed = TRUE) || grepl("\u00a7", s, fixed = TRUE) || grepl("/", s, fixed = TRUE)) {
    return("punctuation")
  }
  if (isTRUE(is_file_echo(s))) {
    return("file-like")
  }
  n_words <- length(strsplit(s, "[[:space:]]+")[[1L]])
  if (n_words < 2L || n_words > 5L) {
    return("word-count")
  }
  low <- tolower(s)
  banned <- c(tolower(c$symbols), tolower(basename(c$files)))
  if (low %in% banned) {
    return("echoes member")
  }
  if (low == as.character(cid) || startsWith(low, "community ") ||
    grepl("^empty ", low)) {
    return("placeholder")
  }
  # a no-op rename is not a rename: it would be counted as relabelled but
  # change nothing (cid 210/217 first `--cids` pass)
  if (low == tolower(labels[[as.character(cid)]] %||% "")) {
    return("unchanged")
  }
  NULL
}

accepted <- list()
rejected <- list()
BATCH <- 12L
for (i in seq(1L, length(targets), by = BATCH)) {
  batch <- targets[i:min(i + BATCH - 1L, length(targets))]
  out <- tryCatch(
    ask(batch),
    error = function(e) {
      message("batch ", (i - 1L) %/% BATCH + 1L, " failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(out)) next
  for (cid in batch) {
    name <- out[[as.character(cid)]]
    why <- bad(cid, name, ctx(cid))
    if (is.null(why)) {
      accepted[[as.character(cid)]] <- clean(name)
    } else {
      rejected[[as.character(cid)]] <- list(name = name, why = why)
    }
  }
  message(sprintf(
    "batch %d/%d: %d accepted, %d rejected",
    (i - 1L) %/% BATCH + 1L, ceiling(length(targets) / BATCH),
    length(accepted), length(rejected)
  ))
}

# second chance: names must be distinct within the smaller batch
if (length(rejected)) {
  rcids <- sort(as.integer(names(rejected)))
  # `unchanged` is rejected only because it matches the current name, so the
  # re-ask has to be told that name is taken or temperature 0 just repeats it
  proposed <- unique(unlist(lapply(rejected, function(r) as.character(r$name %||% ""))))
  for (i in seq(1L, length(rcids), by = 6L)) {
    batch <- rcids[i:min(i + 5L, length(rcids))]
    out <- tryCatch(
      ask(batch, distinct = TRUE, avoid = proposed),
      error = function(e) NULL
    )
    if (is.null(out)) next
    for (cid in batch) {
      name <- out[[as.character(cid)]]
      if (is.null(bad(cid, name, ctx(cid)))) {
        accepted[[as.character(cid)]] <- clean(name)
        rejected[[as.character(cid)]] <- NULL
      }
    }
  }
}

# collision with a label of a community we are not renaming
used <- split(as.integer(names(labels)), unname(labels))
for (cid in names(accepted)) {
  nm <- accepted[[cid]]
  others <- setdiff(used[[nm]] %||% integer(), as.integer(names(accepted)))
  if (length(others)) {
    rejected[[cid]] <- list(name = nm, why = paste("collides with cid", paste(others, collapse = ",")))
  }
}

# duplicate proposals within this pass
final <- accepted[!names(accepted) %in% names(rejected)]
tab <- table(unlist(final, use.names = FALSE))
for (nm in names(tab)[tab > 1L]) {
  for (cid in names(final)[unlist(final) == nm]) {
    rejected[[cid]] <- list(name = nm, why = "duplicate within pass")
  }
}
final <- final[!names(final) %in% names(rejected)]

# names rejected only because two communities proposed the same string: one
# distinct re-ask, since the collision is about the pair, not either name
dup_rej <- names(rejected)[vapply(
  rejected, function(r) identical(r$why, "duplicate within pass"), logical(1)
)]
if (length(dup_rej)) {
  out <- tryCatch(
    ask(
      sort(as.integer(dup_rej)),
      distinct = TRUE,
      avoid = c(unlist(final, use.names = FALSE), vapply(
        rejected[dup_rej], function(r) as.character(r$name %||% ""), character(1)
      ))
    ),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    used_names <- unlist(final, use.names = FALSE)
    for (cid in dup_rej) {
      nm <- out[[cid]]
      if (is.null(bad(cid, nm, ctx(cid))) && !(clean(nm) %in% used_names)) {
        final[[cid]] <- clean(nm)
        used_names <- c(used_names, clean(nm))
      }
    }
  }
}

kept <- setdiff(as.character(targets), names(final))
if (!length(final)) {
  message("no label accepted; nothing written")
  for (cid in names(rejected)) {
    cat(sprintf("  %3s  rejected %s: %s\n", cid, rejected[[cid]]$why, rejected[[cid]]$name))
  }
  quit(status = 0L)
}

old <- lapply(names(final), function(cid) labels[[cid]])
names(old) <- names(final)
for (cid in names(final)) {
  labels[[cid]] <- final[[cid]]
}

jsonlite::write_json(as.list(labels), LABELS, auto_unbox = TRUE, pretty = FALSE)

passes <- if (file.exists(ORIGIN)) {
  prev <- tryCatch(
    jsonlite::fromJSON(ORIGIN, simplifyVector = FALSE),
    error = function(e) list()
  )
  prev$passes %||% list()
} else {
  list()
}
passes[[length(passes) + 1L]] <- list(
  date = format(Sys.Date()),
  model = model,
  prompt = "strict-purpose: Title Case, 2-5 words, no file/symbol/heading echoes",
  requested = length(targets),
  relabelled = length(final),
  kept_previous = length(kept),
  rejected = lapply(rejected[names(rejected) %in% kept], function(r) list(proposal = r$name, reason = r$why)),
  changes = lapply(names(final), function(cid) {
    list(from = old[[cid]], to = final[[cid]])
  }) |>
    stats::setNames(names(final))
)
jsonlite::write_json(list(passes = passes), ORIGIN, auto_unbox = TRUE, pretty = TRUE, null = "null")

cat(sprintf(
  "requested %d, relabelled %d, kept %d, rejected %d\n",
  length(targets), length(final), length(kept), length(rejected)
))
for (cid in names(final)) {
  cat(sprintf("  %3s  %s -> %s\n", cid, old[[cid]] %||% "?", final[[cid]]))
}
for (cid in names(rejected)[names(rejected) %in% kept]) {
  cat(sprintf("  %3s  kept (rejected %s: %s)\n", cid, rejected[[cid]]$why, rejected[[cid]]$name))
}
