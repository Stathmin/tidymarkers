# =============================================================================
#  dev/audit-deps.R  --  the R analogue of deptry/vulture
# =============================================================================
#  Static audit of the dependency and function surface: declared-but-unused
#  `Imports`, unused `importFrom` symbols, duplicate top-level definitions,
#  internal functions nothing calls, and byte-identical bodies. It is a report,
#  not part of `dev/check.R` (CI has no reason to fail on a stale `Suggests`
#  entry), but the first three classes are failures here because they are
#  always mistakes.
#
#  Limits: this reads text, not the loaded namespace. `getExportedValue()`,
#  `requireNamespace()` and NSE (a column name that matches a package) are
#  invisible; S3 methods are skipped because dispatch calls them by class, not
#  by name. A hit therefore means "look here", never "delete this".
#
#  Run:  Rscript dev/audit-deps.R [--strict]
#        exit 1 on unused `Imports`, a duplicate name or a dead internal;
#        `--strict` also exits 1 on identical bodies.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

strict <- "--strict" %in% commandArgs(trailingOnly = TRUE)
failed <- character()
fail <- function(...) failed <<- c(failed, paste0(...))

read_lines <- function(path) {
  if (!file.exists(path)) {
    return(character())
  }
  readLines(path, warn = FALSE)
}

# --- declared dependencies -------------------------------------------------

description <- read.dcf("DESCRIPTION")
field <- function(name) {
  if (!name %in% colnames(description)) {
    return(character())
  }
  parts <- trimws(strsplit(gsub("\n", " ", description[1, name]), ",")[[1]])
  parts <- sub("\\s*\\(.*\\)$", "", parts)
  parts[nzchar(parts) & parts != "R"]
}
imports <- field("Imports")
suggests <- field("Suggests")

r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
r_code <- unlist(lapply(r_files, read_lines))
r_code <- r_code[!grepl("^\\s*#", r_code)]

namespace <- read_lines("NAMESPACE")
imported_symbols <- unique(gsub('"', "", sub(
  "^importFrom\\([^,]+, *\"?([^\"|)]+)\"?.*$", "\\1",
  grep("^importFrom", namespace, value = TRUE)
)))
imported_packages <- unique(gsub('"', "", sub(
  "^import\\(([^,)]+).*$", "\\1", grep("^import\\(", namespace, value = TRUE)
)))

uses_package <- function(pkg) {
  pkg %in% imported_packages ||
    any(grepl(paste0("\\b", gsub("[.]", "\\\\.", pkg), "::"), r_code))
}
unused_imports <- imports[!vapply(imports, uses_package, logical(1))]
if (length(unused_imports)) {
  fail("unused Imports: ", paste(unused_imports, collapse = ", "))
}

tokens_of <- function(lines) {
  strsplit(paste(lines, collapse = "\n"), "[^A-Za-z0-9._]+")[[1]]
}
r_mentions <- table(tokens_of(r_code))
dead_symbols <- imported_symbols[!imported_symbols %in% names(r_mentions)]
if (length(dead_symbols)) {
  fail("importFrom symbols never named in R/: ",
       paste(dead_symbols, collapse = ", "))
}

# Every file that could call a function, including docs and CI scripts.
called_text <- unlist(lapply(c(
  r_files,
  list.files("tests", pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE),
  list.files("data-raw", pattern = "[.]R$", full.names = TRUE),
  list.files("dev", pattern = "[.]R$", full.names = TRUE)
), read_lines))

# --- definitions -----------------------------------------------------------

definitions <- lapply(r_files, function(path) {
  rows <- lapply(parse(path, keep.source = FALSE), function(e) {
    if (!is.call(e) || !identical(e[[1]], as.name("<-")) ||
          !is.call(e[[3]]) || !identical(e[[3]][[1]], as.name("function"))) {
      return(NULL)
    }
    data.frame(
      name = as.character(e[[2]]),
      file = path,
      body = paste(deparse(e[[3]]), collapse = " ")
    )
  })
  do.call(rbind, rows)
})
definitions <- do.call(rbind, definitions)

duplicate_names <- unique(definitions$name[duplicated(definitions$name)])
if (length(duplicate_names)) {
  fail("functions defined more than once: ",
       paste(duplicate_names, collapse = ", "))
}

exports <- gsub('"', "", sub(
  "^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE)
))
methods <- gsub('"', "", sub(
  "^S3method\\(([^,]+), *([^)]+)\\)$", "\\1.\\2",
  grep("^S3method", namespace, value = TRUE)
))
methods <- sub(".*::", "", methods)

# One token count for the whole tree answers "is this name ever mentioned?";
# the definition line itself is the one occurrence every function has.
mentions <- table(tokens_of(called_text))

dead <- character()
for (i in seq_len(nrow(definitions))) {
  name <- definitions$name[i]
  if (name %in% c(exports, methods)) {
    next
  }
  total <- if (name %in% names(mentions)) as.integer(mentions[[name]]) else 0L
  if (total <= sum(definitions$name == name)) {
    dead <- c(dead, paste0(name, " (", definitions$file[i], ")"))
  }
}
if (length(dead)) {
  fail("internal functions never referenced: ", paste(dead, collapse = ", "))
}

twins <- definitions$body[
  duplicated(definitions$body) | duplicated(definitions$body, fromLast = TRUE)
]
if (length(twins)) {
  groups <- vapply(unique(twins), function(body) {
    paste(definitions$name[definitions$body == body], collapse = " == ")
  }, character(1))
  message("identical bodies (delegation is legitimate; check each):")
  message(paste0("  ", groups, collapse = "\n"))
  if (strict) {
    fail("identical function bodies: ", length(unique(twins)), " group(s)")
  }
}

# Never-named `Suggests` is a warning: the package may be needed only by an
# external workflow this scan does not read.
all_text <- unlist(lapply(
  c(
    list.files(".", pattern = "[.](R|Rmd|md|qmd|yml|yaml|Rproj)$", recursive = TRUE),
    list.files(".github", pattern = "[.](yml|yaml|py)$", recursive = TRUE)
  ),
  read_lines
))
unused_suggests <- suggests[!vapply(suggests, function(pkg) {
  any(grepl(paste0("\\b", gsub("[.]", "\\\\.", pkg), "\\b"), all_text))
}, logical(1))]
if (length(unused_suggests)) {
  message("Suggests never named anywhere (check by hand): ",
          paste(unused_suggests, collapse = ", "))
}

# --- summary ---------------------------------------------------------------

message(
  "\n", length(imports), " Imports, ", length(suggests), " Suggests, ",
  nrow(definitions), " top-level functions, ",
  sum(definitions$name %in% exports), " exported"
)
if (length(failed)) {
  message("FAIL")
  message(paste0("  - ", failed, collapse = "\n"))
  quit(status = 1L)
}
message("ok: dependency and function surface is consistent")
