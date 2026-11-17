#!/usr/bin/env Rscript
# =============================================================================
#  dev/graphify-claims.R  --  does a doc sentence describe what the function does?
# =============================================================================
#  `dev/graphify-faithfulness.R` audits a function's own roxygen prose against
#  its body.  This script audits the *other* direction: a sentence in README,
#  VALIDATION, `docs/*.md` or a vignette that names `fn()`, checked against the
#  prose and body of `fn` itself.  The deterministic half of this question
#  (does the name exist at all?) is already answered by loading the package;
#  the half that needs a model is whether the sentence's claim is true.
#
#  One (sentence, function) pair per call, `temperature = 0`, thinking off.
#  Verdicts:
#    supported     the claim about the named function is true of the code
#    partial       true, but the sentence implies more than the code does
#    contradicts   the code does what the sentence says it does not, or the
#                  sentence states a behaviour/default/return the code lacks
#    unverifiable  the sentence is not about this function's own behaviour
#                  (pipeline prose, a file, a workflow), or the body was elided
#
#  Results append to graphify-out/.graphify_claims.jsonl so a re-run resumes;
#  the report goes to graphify-out/claims.md.
#
#  Same model knobs as relabel-bonsai.R / graphify-faithfulness.R:
#    OPENAI_BASE_URL (default http://127.0.0.1:8000/v1)
#    OPENAI_MODEL    (default Ternary-Bonsai-2-27B-tight)
#    OPENAI_API_KEY  (default "local")
#
#  Run:  Rscript dev/graphify-claims.R [--dry-run] [--limit=N] [--only=fn]
#                                      [--path=docs/x.md] [--fresh]
#
#  The model-call block is a deliberate copy of graphify-faithfulness.R: that
#  file was just measured and works, and a shared helper would change it.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

OUT_DIR <- "graphify-out"
RESULT <- file.path(OUT_DIR, ".graphify_claims.jsonl")
REPORT <- file.path(OUT_DIR, "claims.md")

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
fresh <- "--fresh" %in% args
limit_arg <- grep("^--limit=", args, value = TRUE)
limit <- if (length(limit_arg)) as.integer(sub("^--limit=", "", limit_arg[[1L]])) else Inf
only_arg <- grep("^--only=", args, value = TRUE)
only <- if (length(only_arg)) sub("^--only=", "", only_arg[[1L]])
path_arg <- grep("^--path=", args, value = TRUE)
only_path <- if (length(path_arg)) sub("^--path=", "", path_arg[[1L]])

base <- Sys.getenv("OPENAI_BASE_URL", "http://127.0.0.1:8000/v1")
model <- Sys.getenv("OPENAI_MODEL", "Ternary-Bonsai-2-27B-tight")
key <- Sys.getenv("OPENAI_API_KEY", "local")

for (pkg in c("jsonlite", "curl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " is required")
  }
}

# --- every top-level function in R/, prose and body ---------------------------
# Wider than the export list on purpose: process docs name internal helpers too.
DROP_TAG <- "^@(examples|seealso|family|inherit|inheritParams|inheritSection|aliases|name|keywords|noRd|rdname|export|rawNamespace|importFrom|useDynLib|method|md)\\b"

read_funs <- function(file) {
  dlines <- readLines(file, warn = FALSE)
  ex <- parse(file, keep.source = TRUE)
  srs <- attr(ex, "srcref")
  out <- list()
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    if (!is.call(e) || !is.symbol(e[[1L]]) || !as.character(e[[1L]]) %in% c("<-", "=")) next
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
    cut <- match(TRUE, grepl(DROP_TAG, txt))
    if (!is.na(cut)) txt <- txt[seq_len(cut - 1L)]
    txt <- trimws(txt[nzchar(trimws(txt))])
    body <- dlines[seq.int(as.integer(sr[[1L]]), as.integer(sr[[3L]]))]
    # same threshold as the faithfulness pass: the head and the return statement
    # carry the contract, and a body over 200 lines is not worth prompt tokens
    elided <- length(body) > 200L
    if (elided) {
      body <- c(body[seq_len(45L)], "  # ... body elided ...", utils::tail(body, 40L))
    }
    out[[nm]] <- list(file = file, prose = txt, body = body, elided = elided)
  }
  out
}

pool <- list()
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  pool <- c(pool, read_funs(f))
}

# --- doc units: a paragraph, its sentences, each sentence's line --------------
# Fenced code and tables are dropped: a fenced block is an example, not a claim,
# and a table cell has no verb to audit.  Both are skipped, not mis-verdicts.
sentences <- function(path) {
  lines <- readLines(path, warn = FALSE)
  in_fence <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("^\\s*```", lines[[i]])) {
      in_fence <- !in_fence
      lines[[i]] <- ""
    } else if (in_fence || grepl("^\\s*\\|", lines[[i]])) {
      lines[[i]] <- ""
    }
  }
  # front matter of an .Rmd is metadata, not prose
  if (length(lines) > 2L && grepl("^---\\s*$", lines[[1L]])) {
    end <- match(TRUE, grepl("^---\\s*$", lines[-1L]))
    if (!is.na(end)) lines[seq_len(end + 1L)] <- ""
  }
  blank <- !nzchar(trimws(lines))
  out <- list()
  i <- 1L
  while (i <= length(lines)) {
    if (blank[[i]]) {
      i <- i + 1L
      next
    }
    j <- i
    while (j < length(lines) && !blank[[j + 1L]]) j <- j + 1L
    text <- paste(lines[seq.int(i, j)], collapse = "\n")
    m <- gregexpr("(?<=[.!?])\\s+", text, perl = TRUE)[[1L]]
    starts <- c(1L, m + attr(m, "match.length"))
    ends <- c(if (identical(m[[1L]], -1L)) length(text) else m - 1L, nchar(text))
    for (k in seq_along(starts)) {
      s <- substring(text, starts[[k]], ends[[k]])
      out[[length(out) + 1L]] <- list(
        line = i + sum(strsplit(substr(text, 1L, starts[[k]] - 1L), "")[[1L]] == "\n"),
        text = trimws(s)
      )
    }
    i <- j + 1L
  }
  out
}

REF <- "`([A-Za-z_.][A-Za-z0-9_.]*(?:::{0,2}[A-Za-z0-9_.]+)*)\\(\\)`"

claims <- list()
external <- character()
for (f in files <- c(
  "README.md", "VALIDATION.md",
  setdiff(Sys.glob("docs/*.md"), Sys.glob("docs/generated/*")),
  Sys.glob("vignettes/*.Rmd")
)) {
  if (!is.null(only_path) && !identical(f, only_path)) next
  if (!file.exists(f)) next
  for (s in sentences(f)) {
    m <- regmatches(s$text, gregexpr(REF, s$text, perl = TRUE))[[1L]]
    if (!length(m)) next
    nms <- unique(sub("^.*:{2,3}", "", sub(REF, "\\1", m, perl = TRUE)))
    for (nm in nms) {
      if (!nm %in% names(pool)) {
        if (!nm %in% external) external <- c(external, nm)
        next
      }
      claims[[length(claims) + 1L]] <- list(
        id = paste(f, s$line, nm, sep = "|"),
        path = f, line = s$line, fn = nm, text = s$text
      )
    }
  }
}

if (!length(claims)) {
  stop("no claim: no doc sentence names a function defined in R/")
}
if (!is.null(only)) {
  claims <- claims[vapply(claims, function(x) x$fn == only, logical(1))]
}
if (!length(claims)) {
  stop("no claim naming ", only)
}

SYSTEM <- paste0(
  "You audit one sentence of an R package's documentation against the code of ",
  "the function that sentence names. Read the sentence and the function (its ",
  "doc prose, then its body), then give one verdict:\n",
  "- \"supported\": the sentence's claim about the named function is true of ",
  "the code.\n",
  "- \"partial\": the claim is true, but the sentence says more than the code ",
  "does or omits something a reader must know.\n",
  "- \"contradicts\": the sentence states a behaviour, default, return value or ",
  "guarantee that the code does not have, or the code does the opposite.\n",
  "- \"unverifiable\": the sentence is not about this function's own behaviour ",
  "(it describes a package, a pipeline, a file or a workflow and only names the ",
  "function in passing), or the body you were given is visibly incomplete.\n",
  "Only say \"contradicts\" when the code shown directly refutes the sentence. ",
  "When the body delegates to another function, the evidence may sit in code you ",
  "were not given: answer \"partial\".\n",
  "Reply with ONE JSON object and nothing else:\n",
  '{"verdict": "supported", "evidence": "the sentence claim and the code fact"}\n',
  "Keep evidence under 30 words. Judge only the text you are given; do not ",
  "comment on style, tests, performance or missing features."
)

prompt_of <- function(cl) {
  x <- pool[[cl$fn]]
  paste0(
    "Document: ", cl$path, ":", cl$line, "\n",
    "Named function: ", cl$fn, " (", x$file, ")\n\nSENTENCE:\n", cl$text,
    "\n\nDOC PROSE OF ", cl$fn, ":\n",
    if (length(x$prose)) paste(x$prose, collapse = "\n") else "(none)",
    "\n\nFUNCTION BODY:\n", paste(x$body, collapse = "\n")
  )
}

if (dry_run) {
  cat(sprintf(
    "%d claim(s) over %d function(s); %d doc-named function(s) not defined in R/ (external)\n",
    length(claims), length(unique(vapply(claims, function(x) x$fn, character(1)))),
    length(external)
  ))
  if (length(external)) cat("external: ", paste(external, collapse = ", "), "\n", sep = "")
  cat("\n--- system ---\n", SYSTEM, "\n--- user ---\n", prompt_of(claims[[1L]]), "\n", sep = "")
  quit(status = 0L)
}

# --- ask ----------------------------------------------------------------------
verdicts <- c("supported", "partial", "contradicts", "unverifiable")

salvage <- function(text) {
  text <- gsub("```(?:json)?", "", text)
  v <- regmatches(
    text,
    regexpr("(supported|partial|contradicts|unverifiable)", text, ignore.case = TRUE)
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
    verdict = if (length(v)) tolower(v) else NA_character_,
    evidence = trimws(ev)
  )
}

ask <- function(cl) {
  body <- list(
    model = model,
    messages = list(
      list(role = "system", content = SYSTEM),
      list(role = "user", content = prompt_of(cl))
    ),
    temperature = 0,
    max_tokens = 120L,
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
  prev <- tryCatch(jsonlite::stream_in(file(RESULT), verbose = FALSE), error = function(e) NULL)
  keep <- if (is.data.frame(prev) && !is.null(only)) prev[prev$fn != only, , drop = FALSE] else NULL
  if (is.null(keep) || !nrow(keep)) {
    unlink(RESULT)
  } else {
    tmp <- paste0(RESULT, ".tmp")
    con <- file(tmp, "w")
    jsonlite::stream_out(keep, con, verbose = FALSE)
    close(con)
    file.rename(tmp, RESULT)
  }
}
done <- character()
if (file.exists(RESULT)) {
  prev <- tryCatch(jsonlite::stream_in(file(RESULT), verbose = FALSE), error = function(e) NULL)
  if (is.data.frame(prev) && nrow(prev)) done <- prev$id
}
todo <- claims[!vapply(claims, function(x) x$id %in% done, logical(1))]
todo <- utils::head(todo, limit)

cat(sprintf(
  "auditing %d claim(s) with %s (%d already recorded)\n",
  length(todo), model, length(done)
))
failed <- character()
t_start <- Sys.time()
for (k in seq_along(todo)) {
  cl <- todo[[k]]
  t0 <- Sys.time()
  r <- if (isTRUE(pool[[cl$fn]]$elided)) {
    list(verdict = "unverifiable", evidence = "body elided; not sent to the model")
  } else {
    tryCatch(ask(cl), error = function(e) {
      failed <<- c(failed, cl$id)
      message(sprintf("  %-40s error: %s", cl$id, conditionMessage(e)))
      NULL
    })
  }
  if (is.null(r)) next
  if (is.na(r$verdict)) {
    failed <- c(failed, cl$id)
    message(sprintf("  %-40s unparsable answer", cl$id))
    next
  }
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1L)
  cat(
    jsonlite::toJSON(
      list(
        id = cl$id, path = cl$path, line = cl$line, fn = cl$fn,
        verdict = r$verdict, evidence = r$evidence, seconds = secs
      ),
      auto_unbox = TRUE
    ),
    "\n",
    sep = "",
    file = RESULT,
    append = TRUE
  )
  cat(sprintf("  %-46s %-12s %.1fs\n", cl$id, r$verdict, secs))
}

# --- summarise ----------------------------------------------------------------
res <- jsonlite::stream_in(file(RESULT), verbose = FALSE)
res <- res[order(match(res$verdict, verdicts), res$path, res$line), , drop = FALSE]
tbl <- table(factor(res$verdict, levels = verdicts))
find <- res[res$verdict %in% c("contradicts", "partial"), , drop = FALSE]
md <- c(
  "# Doc claims against the code (generated)",
  "",
  "Per (sentence, function) pair: a doc sentence that names `fn()`, checked against",
  "the prose and body of `fn`. One local-model call each; fenced code and tables",
  "are skipped. Deterministic name resolution is not repeated here.",
  sprintf(
    "Model: `%s`. Claims audited: %d over %d function(s). Median %.1f s per call.",
    model, nrow(res), length(unique(res$fn)), stats::median(res$seconds)
  ),
  "",
  paste0("- ", names(tbl), ": ", as.integer(tbl)),
  "",
  if (nrow(find)) {
    c("## Findings", "", sprintf(
      "- **%s:%d** `%s` - `%s`: %s", find$path, find$line, find$fn, find$verdict, find$evidence
    ))
  } else {
    "No finding."
  },
  if (length(external)) {
    c("", "## Named in docs, not defined in `R/`", "",
      paste0("- ", external), "",
      "Base, suggested packages and user-facing verbs; not auditable here.")
  } else NULL,
  if (length(failed)) c("", "## Failed", "", paste0("- ", failed)) else NULL
)
writeLines(md, REPORT)
cat(sprintf(
  "done: %d supported, %d partial, %d contradicts, %d unverifiable of %d (%s)\n",
  tbl[["supported"]], tbl[["partial"]], tbl[["contradicts"]], tbl[["unverifiable"]],
  nrow(res), REPORT
))
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
if (length(failed)) cat("failed: ", paste(failed, collapse = ", "), "\n", sep = "")
