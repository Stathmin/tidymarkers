#!/usr/bin/env Rscript
# =============================================================================
#  dev/graphify-faithfulness.R  --  does the roxygen prose describe the body?
# =============================================================================
#  `dev/graphify-health.R` measures doc<->code *distance* deterministically
#  (`@param` names against `formals()`).  It cannot tell whether the prose is
#  true of the body; that needs a model.  The whole-document semantic pass is
#  unusable here (one 18.6 KB doc costs 692 s at this server's 31 tok/s and
#  never finishes), but one export at a time fits: a few hundred prompt tokens
#  and ~40 output tokens per verdict.
#
#  One export per call, `temperature = 0`, thinking off.  Verdicts:
#    match        every claim the prose makes is in the code
#    partial      true but silent about something a caller must know
#    contradicts  the prose states something the code does not do
#
#  Results are appended to graphify-out/.graphify_faithfulness.jsonl, so a
#  re-run resumes where it stopped; the summary goes to
#  graphify-out/faithfulness.md and to stdout.
#
#  The model is OpenAI-compatible and local (same knobs as relabel-bonsai.R):
#    OPENAI_BASE_URL (default http://127.0.0.1:8000/v1)
#    OPENAI_MODEL    (default Ternary-Bonsai-2-27B-tight)
#    OPENAI_API_KEY  (default "local")
#
#  Run:  Rscript dev/graphify-faithfulness.R [--dry-run] [--limit=N]
#                                            [--only=fn] [--fresh]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

OUT_DIR <- "graphify-out"
RESULT <- file.path(OUT_DIR, ".graphify_faithfulness.jsonl")
REPORT <- file.path(OUT_DIR, "faithfulness.md")

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
fresh <- "--fresh" %in% args
limit_arg <- grep("^--limit=", args, value = TRUE)
limit <- if (length(limit_arg)) as.integer(sub("^--limit=", "", limit_arg[[1L]])) else Inf
only_arg <- grep("^--only=", args, value = TRUE)
only <- if (length(only_arg)) sub("^--only=", "", only_arg[[1L]])

base <- Sys.getenv("OPENAI_BASE_URL", "http://127.0.0.1:8000/v1")
model <- Sys.getenv("OPENAI_MODEL", "Ternary-Bonsai-2-27B-tight")
key <- Sys.getenv("OPENAI_API_KEY", "local")

for (pkg in c("jsonlite", "curl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " is required")
  }
}

# --- what a caller can call ---------------------------------------------------
ns <- readLines("NAMESPACE", warn = FALSE)
exports <- unlist(lapply(grep("^export\\(", ns, value = TRUE), function(s) {
  trimws(strsplit(sub("^[^(]+\\((.*)\\)$", "\\1", s), ",", fixed = TRUE)[[1L]])
}))

# --- prose and body, per export -----------------------------------------------
# Same srcref walk as dev/graphify-health.R: a roxygen block is the run of `#'`
# lines directly above the binding.  Prose stops at `@examples`; the code tags
# below it are examples, not claims about the body.
DROP_TAG <- "^@(examples|seealso|family|inherit|inheritParams|inheritSection|aliases|name|keywords|noRd|rdname|export|rawNamespace|importFrom|useDynLib|method|md)\\b"

read_export <- function(file) {
  dlines <- readLines(file, warn = FALSE)
  ex <- parse(file, keep.source = TRUE)
  srs <- attr(ex, "srcref")
  out <- list()
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    if (!is.call(e) || !as.character(e[[1L]]) %in% c("<-", "=")) next
    rhs <- e[[3L]]
    if (!is.call(rhs) || !identical(rhs[[1L]], as.name("function"))) next
    nm <- as.character(e[[2L]])
    sr <- srs[[i]]
    if (is.null(sr)) next
    j <- as.integer(sr[[1L]]) - 1L
    blk <- character()
    while (j >= 1L && grepl("^#'", dlines[[j]])) {
      blk <- c(dlines[[j]], blk)
      j <- j - 1L
    }
    txt <- sub("^#+'\\s?", "", blk)
    rd <- sub("^@rdname\\s+(\\S+).*$", "\\1", grep("^@rdname\\s+", txt, value = TRUE))
    cut <- match(TRUE, grepl(DROP_TAG, txt))
    if (!is.na(cut)) txt <- txt[seq_len(cut - 1L)]
    txt <- trimws(txt[nzchar(trimws(txt))])
    body <- dlines[seq.int(as.integer(sr[[1L]]), as.integer(sr[[3L]]))]
    # a 300-line body costs prompt tokens and buys no verdict a caller needs;
    # the head and the return statement carry the contract
    if (length(body) > 70L) {
      body <- c(body[seq_len(45L)], "  # ... body elided ...", utils::tail(body, 20L))
    }
    out[[nm]] <- list(
      file = file, prose = txt, body = body,
      topic = if (length(rd)) rd[[1L]] else nm
    )
  }
  out
}

pool <- list()
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  pool <- c(pool, read_export(f))
}
targets <- pool[intersect(exports, names(pool))]
# an `@rdname` alias carries only tags; its prose belongs to the topic owner and is
# audited there. Auditing the alias too double-counts one doc block and reports the
# deliberate deprecation forwarder as a contradiction, so aliases are dropped.
aliases <- names(targets)[!vapply(targets, function(x) length(x$prose) > 0L, logical(1))]
if (length(aliases)) {
  cat(sprintf("%d @rdname alias(es) audited via their topic: %s\n",
    length(aliases), paste(aliases, collapse = ", ")))
}
targets <- targets[setdiff(names(targets), aliases)]
if (!is.null(only)) {
  targets <- targets[intersect(only, names(targets))]
}
if (!length(targets)) {
  stop("no target: none of the requested exports has a roxygen block in R/")
}

SYSTEM <- paste0(
  "You audit whether an R function's documentation describes what the function ",
  "body actually does. Read the doc prose and the code, then give one verdict:\n",
  "- \"match\": every claim the prose makes about behaviour, defaults or the ",
  "return value is true of the code.\n",
  "- \"partial\": the prose is true but silent about something the code does ",
  "that a caller must know.\n",
  "- \"contradicts\": the prose states something the code does not do, or the ",
  "code does the opposite.\n",
  "Only say \"contradicts\" when the code shown directly refutes the prose. When the ",
  "body delegates to another function, or carries a `# ... body elided ...` line, the ",
  "evidence may sit in code you were not given: answer \"partial\" instead.\n",
  "Reply with ONE JSON object and nothing else:\n",
  '{"verdict": "match", "evidence": "the prose claim and the code fact"}\n',
  "Keep evidence under 30 words. Judge only the text you are given; do not ",
  "comment on style, tests, performance or missing features."
)

prompt_of <- function(nm) {
  x <- targets[[nm]]
  paste0(
    "Export: ", nm, " (", x$file, ")\n\nDOC PROSE:\n",
    paste(x$prose, collapse = "\n"),
    "\n\nFUNCTION BODY:\n", paste(x$body, collapse = "\n")
  )
}

if (dry_run) {
  cat("targets: ", length(targets), " export(s)\n", sep = "")
  empty <- names(targets)[!vapply(targets, function(x) length(x$prose) > 0L, logical(1))]
  if (length(empty)) {
    stop("no prose for: ", paste(empty, collapse = ", "))
  }
  nm <- names(targets)[[1L]]
  cat("\n--- system ---\n", SYSTEM, "\n--- user (", nm, ") ---\n", prompt_of(nm), "\n", sep = "")
  quit(status = 0L)
}

# --- ask ----------------------------------------------------------------------
verdicts <- c("match", "partial", "contradicts")

salvage <- function(text) {
  text <- gsub("```(?:json)?", "", text)
  v <- regmatches(
    text,
    regexpr("(match|partial|contradicts)", text, ignore.case = TRUE)
  )
  ev <- regmatches(
    text,
    regexpr('"evidence"[[:space:]]*:[[:space:]]*"([^"]*)"', text, perl = TRUE)
  )
  if (length(ev)) {
    ev <- sub('^.*?"([^"]*)"$', "\\1", ev, perl = TRUE)
  } else {
    ev <- ""
  }
  list(
    verdict = if (length(v)) tolower(v) %||% NA_character_ else NA_character_,
    evidence = trimws(ev)
  )
}

ask <- function(nm) {
  body <- list(
    model = model,
    messages = list(
      list(role = "system", content = SYSTEM),
      list(role = "user", content = prompt_of(nm))
    ),
    temperature = 0,
    max_tokens = 120L,
    # measured on this llama.cpp build: the {"thinking":...} shape is ignored,
    # these two together keep output to one JSON object
    chat_template_kwargs = list(enable_thinking = FALSE),
    reasoning_effort = "low"
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

# --- run ----------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE)
if (fresh && file.exists(RESULT)) {
  # `--fresh --only=fn` must not wipe the other exports' records
  prev <- tryCatch(
    jsonlite::stream_in(file(RESULT), verbose = FALSE),
    error = function(e) NULL
  )
  keep <- if (is.data.frame(prev) && !is.null(only)) {
    prev[!prev$export %in% only, , drop = FALSE]
  } else {
    NULL
  }
  if (is.null(keep) || !nrow(keep)) {
    unlink(RESULT)
  } else {
    # write-then-rename, so a crash mid-rewrite cannot lose the kept records
    tmp <- paste0(RESULT, ".tmp")
    con <- file(tmp, "w")
    jsonlite::stream_out(keep, con, verbose = FALSE)
    close(con)
    file.rename(tmp, RESULT)
  }
}
done <- character()
if (file.exists(RESULT)) {
  prev <- tryCatch(
    jsonlite::stream_in(file(RESULT), verbose = FALSE),
    error = function(e) NULL
  )
  if (is.data.frame(prev) && nrow(prev)) done <- prev$export
}
todo <- setdiff(names(targets), done)
todo <- utils::head(todo, limit)

cat(sprintf(
  "auditing %d export(s) with %s (%d already recorded)\n",
  length(todo), model, length(done)
))
failed <- character()
for (k in seq_along(todo)) {
  nm <- todo[[k]]
  t0 <- Sys.time()
  r <- tryCatch(ask(nm), error = function(e) {
    failed <<- c(failed, nm)
    message(sprintf("  %-28s error: %s", nm, conditionMessage(e)))
    NULL
  })
  if (is.null(r)) next
  if (is.na(r$verdict)) {
    failed <- c(failed, nm)
    message(sprintf("  %-28s unparsable answer", nm))
    next
  }
  cat(
    jsonlite::toJSON(
      list(
        export = nm,
        file = targets[[nm]]$file,
        verdict = r$verdict,
        evidence = r$evidence,
        seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1L)
      ),
      auto_unbox = TRUE
    ),
    "\n",
    sep = "",
    file = RESULT,
    append = TRUE
  )
  cat(sprintf("  %-28s %-11s %.1fs\n", nm, r$verdict, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# --- summarise ----------------------------------------------------------------
res <- jsonlite::stream_in(file(RESULT), verbose = FALSE)
bad <- res[res$verdict %in% c("partial", "contradicts"), , drop = FALSE]
bad <- bad[order(match(bad$verdict, c("contradicts", "partial")), bad$export), , drop = FALSE]
tbl <- table(factor(res$verdict, levels = verdicts))
md <- c(
  "# Doc<->code faithfulness (generated)",
  "",
  "Per export: roxygen prose against the function body, one local-model call each.",
  sprintf("Model: `%s`. Export audited: %d. Median %.1f s per call.", model, nrow(res), stats::median(res$seconds)),
  "",
  paste0("- ", names(tbl), ": ", as.integer(tbl)),
  "",
  if (nrow(bad)) c("## Findings", "", paste0("- **", bad$export, "** (", bad$file, ") `", bad$verdict, "` - ", bad$evidence)) else "No finding.",
  if (length(failed)) c("", "## Failed", "", paste0("- ", failed))
)
writeLines(md, REPORT)
cat(sprintf(
  "done: %d match, %d partial, %d contradicts of %d (%s)\n",
  tbl[["match"]], tbl[["partial"]], tbl[["contradicts"]], nrow(res), REPORT
))
if (length(failed)) cat("failed: ", paste(failed, collapse = ", "), "\n", sep = "")
