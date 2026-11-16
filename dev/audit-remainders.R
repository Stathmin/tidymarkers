#!/usr/bin/env Rscript
# Audit ROADMAP.md for closed items that name a deferred part without pointing
# at where that remainer is tracked.
#
# The rule (ROADMAP.md, "## Остатки в закрытых пунктах"): a closed item (`[x]`)
# whose text defers something ("кроме", "пока", "только для", ...) must point at
# either an entry in that section (`№N`) or the parked list.
#
#   Rscript dev/audit-remainders.R           # report only
#   Rscript dev/audit-remainders.R --check   # exit 1 when an item is unpointed
#
# A closed item also defers when its body still carries an open (`[ ]`) or
# partial (`[~]`) sub-item: that is the multiline form of "done, except ...",
# and no wording check on the parent line can see it.
#
# The audit is deliberately word-based: it cannot tell an honest scope note from
# a hidden remainer. Both must be worded so they do not read as a deferral, or
# pointed at. Every hit names what tripped it.

args <- commandArgs(trailingOnly = TRUE)
check_mode <- "--check" %in% args

roadmap <- "ROADMAP.md"
lines <- readLines(roadmap, encoding = "UTF-8", warn = FALSE)

# Only top-level items start a new body; anything indented belongs to the item
# above it, so a nested box cannot hide its parent.
item_re <- "^([ \t]*)(?:-|[0-9]+\\.)[ \t]+\\[(.)\\]"
defer_re <- paste0(
  "кроме|остаё|Остаё|остат|Остат|пока|только для|не покрыт|",
  "отложен|не сделан|TODO|FIXME"
)
pointer_re <- paste0(
  "Остатки в закрытых пунктах|Запаркованные остатки|№[ \t]*[0-9]+|",
  "§[A-Z][A-Z0-9.]*|подпункт[ \t]*:"
)

items <- list()
current <- NULL
for (i in seq_along(lines)) {
  line <- lines[[i]]
  if (grepl("^#{2,3} ", line)) {
    if (!is.null(current)) items[[length(items) + 1L]] <- current
    current <- NULL
    next
  }
  matched <- regmatches(line, regexec(item_re, line))[[1L]]
  if (length(matched) == 3L && !nzchar(matched[[2L]])) {
    if (!is.null(current)) items[[length(items) + 1L]] <- current
    current <- list(line = i, status = matched[[3L]], body = line)
    next
  }
  if (!is.null(current)) current$body <- paste(current$body, line, sep = "\n")
}
if (!is.null(current)) items[[length(items) + 1L]] <- current

closed <- Filter(function(it) identical(it$status, "x"), items)
open_sub_re <- "(?m)^[ \t]+(?:-|[0-9]+\\.)[ \t]+\\[[ ~]\\]"
hit_word <- function(body) {
  if (grepl(open_sub_re, body, perl = TRUE)) return("open sub-item")
  words <- regmatches(body, gregexpr(defer_re, body, perl = TRUE))[[1L]]
  if (identical(words, character(0))) NA_character_ else words[[1L]]
}
closed <- Filter(function(it) !is.na(hit_word(it$body)), closed)

pointed <- Filter(function(it) grepl(pointer_re, it$body), closed)
unpointed <- Filter(function(it) !grepl(pointer_re, it$body), closed)

title_of <- function(it) {
  head <- strsplit(it$body, "\n", fixed = TRUE)[[1L]][[1L]]
  trimws(sub(".*\\[[^]]*\\][ \t]*", "", head))
}

cat("roadmap items:", length(items), "\n")
cat("closed items naming a deferral:", length(closed), "\n")
cat("  pointed at a tracked remainer:", length(pointed), "\n")
cat("  unpointed:", length(unpointed), "\n")

if (length(unpointed)) {
  cat("\nClosed items that read as deferrals without a pointer:\n")
  for (it in unpointed) {
    cat("  l.", it$line, " [", hit_word(it$body), "] ", title_of(it), "\n", sep = "")
  }
}

if (check_mode && length(unpointed)) {
  cat("\nFix: reword the note so it is not a deferral, or add",
      "(`## Остатки в закрытых пунктах`, №N) / (`Запаркованные остатки`, `ID`).\n")
  quit(status = 1L)
}
invisible(NULL)
