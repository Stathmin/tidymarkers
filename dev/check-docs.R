# =============================================================================
#  dev/check-docs.R  --  the docs graph gate
# =============================================================================
#  Checks the docs layer built around docs/INDEX.md and docs/artifacts.yml:
#
#    1. no orphan doc -- every document in scope is linked from docs/INDEX.md
#    2. no broken relative link -- every relative link under docs/ resolves
#    3. every include_graphics() in vignettes/*.Rmd has a manifest entry
#    4. backlinks agree -- an artifact's `docs` are exactly the docs that
#       actually reference it, and vice versa
#
#  Scope: README/AGENTS/CONTRIBUTING/ROADMAP/NEWS/INVENTORY, vignettes/*.Rmd,
#  dev/*.md and docs/**/*.md.  report/ and lecture/ are out of scope (their
#  includes are generated at render time and are not part of the product docs).
#
#  Run:  Rscript dev/check-docs.R          (exit 1 on failure)
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)
source(file.path(root, "dev", "artifacts.R"))

failures <- character()
note <- function(...) failures <<- c(failures, paste0(...))
ok <- function(...) message("  ok   ", ...)

doc_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

resolve_rel <- function(doc, target) {
  t <- trimws(target)
  t <- sub("^<", "", t)
  t <- sub(">$", "", t)
  t <- sub("[#?].*$", "", t)
  if (!nzchar(t)) {
    return("")
  }
  # schemes (https:, mailto:, ...) are not relative links
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*:", t)) {
    return("")
  }
  t <- if (grepl("^/", t)) sub("^/+", "", t) else file.path(dirname(doc), t)
  parts <- strsplit(normalize_slashes(t), "/", fixed = TRUE)[[1L]]
  out <- character()
  for (p in parts) {
    if (!nzchar(p) || p == ".") next
    if (p == "..") {
      if (length(out)) out <- out[-length(out)]
    } else {
      out <- c(out, p)
    }
  }
  paste(out, collapse = "/")
}

extract_links <- function(text) {
  md <- regmatches(text, gregexpr("\\]\\(([^)]+)\\)", text, perl = TRUE))[[1L]]
  md <- sub("^\\]\\((.*)\\)$", "\\1", md)
  html <- regmatches(
    text,
    gregexpr("(?:src|href)\\s*=\\s*[\"']([^\"']+)[\"']", text, perl = TRUE)
  )[[1L]]
  html <- sub("^.*=[\"']([^\"']+)[\"']$", "\\1", html)
  unique(c(md, html))
}

extract_include_graphics <- function(text) {
  hits <- regmatches(
    text,
    gregexpr("include_graphics\\(\\s*[\"']([^\"']+)[\"']", text, perl = TRUE)
  )[[1L]]
  sub("^include_graphics\\(\\s*[\"'](.*)[\"']$", "\\1", hits)
}

doc_targets <- function(doc) {
  text <- doc_text(doc)
  targets <- c(extract_links(text), extract_include_graphics(text))
  targets <- vapply(targets, function(t) resolve_rel(doc, t), character(1))
  unique(targets[nzchar(targets)])
}

docs_scope <- function() {
  top <- c(
    "README.md", "AGENTS.md", "CONTRIBUTING.md", "ROADMAP.md", "NEWS.md",
    "INVENTORY.md"
  )
  top <- top[file.exists(top)]
  vign <- list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)
  dev <- list.files("dev", pattern = "[.]md$", full.names = TRUE)
  docs <- list.files("docs", pattern = "[.]md$", recursive = TRUE, full.names = TRUE)
  unique(normalize_slashes(c(top, vign, dev, docs)))
}

# 1. no orphan doc ----------------------------------------------------------
index <- "docs/INDEX.md"
if (!file.exists(index)) {
  note(index, ": missing (the canonical index)")
} else {
  indexed <- doc_targets(index)
  scope <- setdiff(docs_scope(), index)
  orphans <- scope[!scope %in% indexed]
  if (length(orphans)) {
    for (o in orphans) note(o, ": orphan doc (not linked from docs/INDEX.md)")
  } else {
    ok(paste0(length(scope), " documents are all linked from docs/INDEX.md"))
  }
}

# 2. no broken relative link (under docs/) ----------------------------------
docs_md <- list.files("docs", pattern = "[.]md$", recursive = TRUE, full.names = TRUE)
broken <- character()
for (doc in docs_md) {
  for (target in doc_targets(doc)) {
    if (!file.exists(target)) broken <- c(broken, paste0(doc, " -> ", target))
  }
}
if (length(broken)) {
  for (b in broken) note("broken relative link: ", b)
} else {
  ok("no broken relative links in ", length(docs_md), " docs/*.md files")
}

# 3 + 4. include_graphics and backlinks -------------------------------------
artifacts <- read_manifest()
paths <- vapply(artifacts, function(a) a$path, character(1))
by_path <- stats::setNames(artifacts, paths)

includes <- list()
for (v in list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)) {
  v <- normalize_slashes(v)
  for (target in extract_include_graphics(doc_text(v))) {
    resolved <- resolve_rel(v, target)
    includes[[resolved]] <- unique(c(includes[[resolved]], v))
  }
}

for (target in names(includes)) {
  docs <- includes[[target]]
  if (!target %in% paths) {
    note(
      "include_graphics(", target, ") has no docs/artifacts.yml entry (",
      paste(docs, collapse = ", "), ")"
    )
  } else {
    missing <- setdiff(docs, by_path[[target]]$docs)
    if (length(missing)) {
      note(
        target, ": docs not declared in the manifest: ",
        paste(missing, collapse = ", ")
      )
    }
  }
}

for (a in artifacts) {
  for (doc in a$docs) {
    if (!file.exists(doc)) {
      note(a$path, ": declared doc does not exist: ", doc)
    } else if (!a$path %in% doc_targets(doc)) {
      note(a$path, ": declared doc does not reference the artifact: ", doc)
    }
  }
}

if (!length(failures)) {
  ok("include_graphics and manifest backlinks agree")
}

if (length(failures)) {
  message("docs check FAILED:")
  for (f in failures) message("  - ", f)
  quit(status = 1L)
}
message("ok: docs graph is consistent")
