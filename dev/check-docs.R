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
#    5. an Rd topic that is not exported carries no \examples{} -- R CMD check
#       attaches only the exports, so such an example dies at check time
#
#  Scope: README/AGENTS/CONTRIBUTING/ROADMAP/NEWS/VALIDATION, vignettes/*.Rmd,
#  dev/*.md and docs/**/*.md.  dev/lecture/ is out of scope (its includes are
#  generated at render time and are not part of the product docs).
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
source(file.path(root, "dev", "check-reporter.R"))

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
    "VALIDATION.md"
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

# 5. the graph's limit is named (ZA8.2) ---------------------------------------
#  "the graph is empty" and "the code is clean" look the same, so the generated
#  graph report must say so; dev/graphify-health.R --write emits the sentence.
health <- "docs/generated/graphify-health.md"
limit <- paste(
  "The code graph is a **limitation of orientation, not proof of correctness**:",
  "an empty graph and clean code look the same."
)
if (!file.exists(health)) {
  note(health, ": missing (the graph limit must be documented)")
} else if (!any(grepl(limit, readLines(health, warn = FALSE), fixed = TRUE))) {
  note(health, ": does not name the limitation (rerun dev/graphify-health.R --write)")
} else {
  ok("the code graph's limitation is named")
}

# 6. the three-level articulation is stated (ZY16) ---------------------------
#  The README must say what the package is beyond "a primer designer"; the
#  sentence is checked literally so a rewrite cannot drop the framing silently.
levels <- paste(
  "tidymarkers is three things at once: a marker design tool for a single",
  "researcher (Level 1), a unified grammar for marker design across chemistries",
  "(Level 2), and infrastructure for community-wide method refinement through",
  "aggregated design traces (Level 3). Each level stands on its own; together",
  "they define a new way to work."
)
readme <- gsub("[[:space:]]+", " ", doc_text("README.md"))
if (!grepl(levels, readme, fixed = TRUE)) {
  note("README.md: does not state the three-level articulation (ZY16)")
} else {
  ok("README.md states the three-level articulation")
}

# 7. the README quickstart is executed on every build (ZY15) ------------------
#  A snippet that only exists can rot; the same body is an eval-enabled chunk in
#  a vignette, so the build runs it and this check keeps the two copies equal.
fence_body <- function(lines, open) {
  if (open >= length(lines)) {
    return(NULL)
  }
  rest <- lines[seq.int(open + 1L, length(lines))]
  close <- which(grepl("^```\\s*$", rest))
  if (!length(close)) {
    return(NULL)
  }
  body <- rest[seq_len(close[[1L]] - 1L)]
  body <- sub("[[:space:]]+$", "", body)
  while (length(body) && !nzchar(body[[1L]])) body <- body[-1L]
  while (length(body) && !nzchar(body[[length(body)]])) body <- body[-length(body)]
  body
}

vig_lines <- readLines("vignettes/kasp-intro.Rmd", warn = FALSE, encoding = "UTF-8")
quick <- grep("^```\\{r +quickstart[ ,}]", vig_lines)
quick_body <- if (length(quick) == 1L) fence_body(vig_lines, quick) else NULL

readme_lines <- readLines("README.md", warn = FALSE, encoding = "UTF-8")
snippets <- lapply(grep("^```r\\s*$", readme_lines), function(i) {
  fence_body(readme_lines, i)
})
matched <- length(quick_body) > 0L &&
  any(vapply(snippets, function(b) identical(b, quick_body), logical(1)))

if (length(quick) != 1L) {
  note("vignettes/kasp-intro.Rmd: no single chunk labelled 'quickstart' (ZY15)")
} else if (grepl("eval\\s*=\\s*FALSE", vig_lines[[quick]])) {
  note("vignettes/kasp-intro.Rmd: the quickstart chunk is eval = FALSE (ZY15)")
} else if (!matched) {
  note("README.md: no ```r snippet is identical to the vignette quickstart (ZY15)")
} else {
  ok("README.md quickstart is the executed vignette chunk")
}

# 8. the counts quoted in prose are re-derived from the data (ZY17) ---------
#  Each number below is computed from inst/extdata/ and then searched for in the
#  document that quotes it, so a data change that the prose does not follow
#  fails here. It is R and not `grep` because three ledger cells name several
#  tables at once. Rung counts (primary_literature/vendor_guideline/
#  rule_of_thumb) come from criterion_spec() and are gated by
#  dev/check-criteria.R, which loads the package.
read_ext <- function(name) {
  utils::read.csv(file.path("inst/extdata", name), stringsAsFactors = FALSE)
}
flat <- function(path) {
  text <- gsub(intToUtf8(0x2212), "-", doc_text(path), fixed = TRUE)
  gsub("[[:space:]]+", " ", text)
}
sources <- read_ext("scientific_sources.csv")
ledger <- read_ext("science_ledger.csv")
prov <- read_ext("criterion_provenance.csv")
extracts <- list.files("inst/extdata", pattern = "[.]csv$")
extract_names <- function(rows) {
  raw <- rows$extract[!is.na(rows$extract)]
  names <- basename(trimws(unlist(strsplit(raw, ",", fixed = TRUE))))
  unique(names[grepl("[.]csv$", names)])
}
named <- extract_names(ledger)
tested <- extract_names(ledger[ledger$use == "test", , drop = FALSE])
blocked_extract <- trimws(ifelse(is.na(ledger$extract), "", ledger$extract))
blocked_extract <- blocked_extract[ledger$use == "blocked"]
blocked_none <- sum(!nzchar(blocked_extract) | tolower(blocked_extract) == "none")
counts <- c(
  sources = nrow(sources),
  extracts = length(extracts),
  named = length(named),
  tested = length(tested),
  ledger = nrow(ledger),
  test = sum(ledger$use == "test"),
  report = sum(ledger$use == "report"),
  blocked = sum(ledger$use == "blocked"),
  rule = sum(ledger$use == "rule"),
  artefact = sum(ledger$use == "artefact"),
  config = sum(ledger$use == "config"),
  blocked_none = blocked_none,
  blocked_real = sum(ledger$use == "blocked") - blocked_none,
  criteria = nrow(prov),
  named_source = sum(prov$reason_class == "named_source"),
  assumption = sum(prov$source == "package assumption"),
  verified = sum(sources$verified, na.rm = TRUE),
  fulltext = sum(sources$access == "fulltext"),
  supplement = sum(sources$access == "supplement"),
  abstract = sum(sources$access == "abstract"),
  closed = sum(sources$access == "closed access"),
  code = sum(sources$access == "source"),
  stub = sum(sources$access == "stub"),
  transcribed = sum(sources$access == "transcribed"),
  notfound = sum(sources$access == "not found")
)

must_quote <- function(path, snippets) {
  text <- flat(path)
  miss <- snippets[!vapply(snippets, grepl, logical(1), x = text, fixed = TRUE)]
  if (length(miss)) {
    for (m in miss) note(path, ": does not quote the derived count: '", m, "'")
  } else {
    ok(path, ": quotes ", length(snippets), " counts re-derived from inst/extdata/")
  }
}

must_quote("VALIDATION.md", c(
  paste0(counts[["sources"]], " registered sources"),
  paste0(counts[["extracts"]], " extracted tables"),
  paste0(counts[["named"]], " extract tables are named"),
  paste0("the ", counts[["extracts"]], "th CSV"),
  paste0(counts[["tested"]], " of them"),
  paste0(
    "Ledger rows by use: ", counts[["test"]], " `test`, ",
    counts[["report"]], " `report`"
  ),
  paste0(
    "Of the ", counts[["blocked"]], " `blocked` rows, ",
    counts[["blocked_none"]], " name no extract"
  ),
  paste0(counts[["blocked_real"]], " name a real extract"),
  paste0(
    counts[["fulltext"]], " sources were read in full text, ",
    counts[["supplement"]], " through a supplement, ",
    counts[["abstract"]], " as an abstract only, ",
    counts[["closed"]], " are closed access, ",
    counts[["code"]], " are source code, ",
    if (counts[["stub"]] == 1L) "1 is a stub, " else paste0(counts[["stub"]], " are stubs, "),
    counts[["transcribed"]], " was transcribed, ",
    counts[["notfound"]], " were not found"
  ),
  paste0(counts[["verified"]], " of the ", counts[["sources"]], " rows are marked `verified`"),
  paste0(counts[["closed"]], " sources are closed access")
))

must_quote("ROADMAP.md", c(
  paste0(counts[["sources"]], " registered sources"),
  paste0(counts[["extracts"]], " extract tables"),
  paste0(counts[["test"]], " `test` rows"),
  paste0(counts[["ledger"]], " \u0441\u0442\u0440\u043e\u043a ledger"),
  paste0(counts[["criteria"]], " criterion rows"),
  paste0(counts[["named_source"]], " \u0441 \u043d\u0430\u0437\u0432\u0430\u043d\u043d\u044b\u043c")
))

must_quote("docs/decision-audit.md", c(
  paste0("rows of ", counts[["criteria"]]),
  paste0("source registry (", counts[["sources"]], " rows)"),
  paste0("registry (", counts[["sources"]], " rows) | ", counts[["named_source"]]),
  paste0(counts[["assumption"]], " of ", counts[["criteria"]], " rows")
))

must_quote("vignettes/reading-results.Rmd", c(
  paste0(counts[["criteria"]], " criterion rows")
))

# An Rd topic that is not exported must not carry \examples{} -----------------
#  R CMD check runs the examples with only the exports attached, so an example
#  on an internal topic dies there with `could not find function "..."` -- 12
#  minutes into CI, not here.  S3 methods and data topics are reachable
#  without an export of their own, so they are exempt.
ns <- readLines("NAMESPACE", warn = FALSE)
exported <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
s3_parts <- regmatches(ns, regexec("^S3method\\(\\s*([^,]+),\\s*([^)]+)\\)", ns))
s3_names <- vapply(
  s3_parts[vapply(s3_parts, length, integer(1)) == 3L],
  function(m) {
    paste0(gsub('"', "", trimws(m[[2L]])), ".", gsub('"', "", trimws(m[[3L]])))
  },
  character(1)
)
reachable <- c(exported, s3_names)
for (rd in list.files("man", pattern = "[.]Rd$", full.names = TRUE)) {
  text <- doc_text(rd)
  if (!grepl("\\\\examples\\{", text)) next
  if (grepl("\\\\docType\\{data\\}", text)) next
  aliases <- regmatches(text, gregexpr("\\\\alias\\{([^}]*)\\}", text))[[1L]]
  if (length(aliases)) aliases <- gsub("^\\\\alias\\{|\\}$", "", aliases)
  if (length(intersect(aliases, reachable))) next
  note(rd, ": \\examples{} on a topic that is not exported (R CMD check cannot run it)")
}

# The three-step Tm decomposition is one number set quoted in five places -----
#  The steps are a chain (-16.9 -> -9.6 -> ~-5.4 C) and a reader who sees only
#  one of them in one file cannot tell whether the others were updated.
for (path in c(
  "README.md", "VALIDATION.md", "R/marker_dpcr.R",
  "man/design_dpcr_marker.Rd", "data-raw/mm2020.R"
)) {
  body <- sub("^\\s*#'?\\s?", "", readLines(path, warn = FALSE, encoding = "UTF-8"))
  body <- gsub(intToUtf8(0x2212), "-", paste(body, collapse = " "), fixed = TRUE)
  body <- gsub("[[:space:]]+", " ", body)
  chain <- grepl("-16[.]9.*-9[.]6.*-5[.]4", body)
  label <- grepl("three (recorded )?steps", body, ignore.case = TRUE)
  if (chain && label) {
    ok(path, ": carries the three-step Tm decomposition")
  } else {
    note(path, ": does not carry the three-step Tm decomposition (-16.9/-9.6/-5.4)")
  }
}

if (length(failures)) {
  message("docs check FAILED:")
  for (f in failures) message("  - ", f)
  quit(status = 1L)
}
message("ok: docs graph is consistent")
