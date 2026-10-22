# =============================================================================
#  dev/artifacts.R  --  maintain the committed-artifact manifest + projection
# =============================================================================
#  docs/artifacts.yml is the source of truth for every committed artifact and
#  docs/generated/artifacts.md is its Markdown projection (the docs/artifacts
#  graph that graphify indexes).  This script owns the hashes and the
#  projection; the metadata (kind, generator, inputs, docs, code) is curated by
#  hand in the manifest.
#
#  Each entry:
#    path         repo-relative artifact
#    kind         figure | decision | table | benchmark
#    generator    the script that writes it
#    inputs       repo-relative source/data paths the generator reads
#    input_hash   md5 over the generator + inputs (content, not mtime)
#    output_hash  md5 of the artifact itself
#    docs         the docs that embed the artifact (the backlinks)
#    code         the package source the artifact is derived from
#    status       verified | unverified  (unverified: hashes are reported, not
#                 enforced -- use it when a generator cannot be justified)
#
#  The manifest is YAML.  `yaml` is a dev-only tool; if it is absent the scripts
#  fall back to `jsonlite`, which reads the JSON subset (`--write` then emits
#  JSON, which is valid YAML).  Neither is added to DESCRIPTION: `dev/` is
#  outside the built package (`.Rbuildignore`).
#
#  Only repository inputs are hashed.  OligoArrayAux (for the gallery) is
#  external and versioned by DESCRIPTION/SystemRequirements; the layout engine
#  of the decision trees is vendored under .tools/elk, but its own sources are
#  not hashed, so a layout-only change does not mark the trees stale.
#
#  Run:
#    Rscript dev/artifacts.R --write   # refresh hashes + docs/generated/artifacts.md
#    Rscript dev/artifacts.R --check   # fail (exit 1) on stale hashes/projection
# =============================================================================

ARTIFACT_MANIFEST <- "docs/artifacts.yml"
ARTIFACT_PROJECTION <- "docs/generated/artifacts.md"
ARTIFACT_HASH_SALT <- "tidymarkers-artifacts-v1"
ARTIFACT_FIELDS <- c(
  "path", "kind", "generator", "inputs", "input_hash", "output_hash",
  "docs", "code", "status"
)

normalize_slashes <- function(x) gsub("\\\\", "/", x)

md5_file <- function(path) {
  if (!file.exists(path)) stop("missing file: ", path, call. = FALSE)
  unname(tools::md5sum(path))
}

# Content hash over files and directories.  Directories are expanded
# recursively and the relative file list is part of the hash, so deleting or
# renaming a source file also invalidates the artifact.
hash_paths <- function(paths) {
  parts <- character()
  for (p in unique(paths)) {
    if (dir.exists(p)) {
      rel <- list.files(p, recursive = TRUE, all.files = FALSE, no.. = TRUE)
      rel <- sort(rel[!dir.exists(file.path(p, rel))])
      if (!length(rel)) parts <- c(parts, paste0(p, "/:<empty>"))
      for (r in rel) {
        f <- file.path(p, r)
        parts <- c(parts, paste0(normalize_slashes(f), ":", md5_file(f)))
      }
    } else {
      parts <- c(parts, paste0(normalize_slashes(p), ":", md5_file(p)))
    }
  }
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(sort(c(ARTIFACT_HASH_SALT, parts)), tmp)
  unname(tools::md5sum(tmp))
}

artifact_inputs <- function(a) as.character(unlist(a$inputs, use.names = FALSE))

artifact_input_hash <- function(a) {
  hash_paths(c(a$generator, artifact_inputs(a)))
}

artifact_output_hash <- function(a) md5_file(a$path)

# YAML when available, JSON (a YAML subset) via jsonlite otherwise.
read_manifest <- function(path = ARTIFACT_MANIFEST) {
  if (!file.exists(path)) stop("no manifest at ", path, call. = FALSE)
  doc <- if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::read_yaml(path, eval.expr = FALSE)
  } else if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  } else {
    stop("need the dev-only `yaml` (or `jsonlite`) package to read ", path, call. = FALSE)
  }
  artifacts <- doc$artifacts
  if (!is.list(artifacts) || !length(artifacts)) {
    stop("manifest has no `artifacts` list: ", path, call. = FALSE)
  }
  lapply(artifacts, function(a) {
    a$inputs <- as.character(unlist(a$inputs, use.names = FALSE))
    a$docs <- as.character(unlist(a$docs, use.names = FALSE))
    a$code <- as.character(unlist(a$code, use.names = FALSE))
    a
  })
}

manifest_header <- c(
  "# Artifact manifest -- committed artifacts only (the source of truth).",
  "#",
  "# Fields: path, kind, generator, inputs, input_hash, output_hash, docs, code,",
  "# status.",
  "#   input_hash  md5 over the generator and its inputs (content, not mtime)",
  "#   output_hash md5 of the artifact itself",
  "#   docs        the docs that embed the artifact (the backlinks)",
  "#   code        the package source the artifact is derived from",
  "#   status      verified | unverified",
  "#",
  "# Hashes cover repository inputs only. OligoArrayAux (for the gallery) is",
  "# external and versioned by DESCRIPTION/SystemRequirements; the decision trees'",
  "# layout engine is vendored under .tools/elk and its sources are not hashed.",
  "#",
  "# vignettes/figures/*.png are generated locally and git-ignored, so they are",
  "# deliberately absent from this manifest. The .txt twins of the decision trees",
  "# are committed next to the SVGs for grepping, but they are not images and so",
  "# are not artifacts here.",
  "#",
  "# Refresh: Rscript dev/artifacts.R --write",
  "# Check:   Rscript dev/check-artifacts.R"
)

write_manifest <- function(artifacts, path = ARTIFACT_MANIFEST) {
  payload <- list(artifacts = artifacts)
  body <- if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::as.yaml(payload)
  } else if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE)
  } else {
    stop("need the dev-only `yaml` (or `jsonlite`) package to write ", path, call. = FALSE)
  }
  writeLines(c(manifest_header, body), path, useBytes = TRUE)
  invisible(path)
}

rel_from_generated <- function(paths) paste0("../../", normalize_slashes(paths))

link_list <- function(paths) {
  if (!length(paths)) {
    return("(none)")
  }
  paste0("[", paths, "](", rel_from_generated(paths), ")", collapse = ", ")
}

projection_lines <- function(artifacts) {
  lines <- c(
    "# Artifact graph (generated)",
    "",
    "> Generated by `Rscript dev/artifacts.R --write` from",
    "> [`docs/artifacts.yml`](../artifacts.yml). Do not edit by hand.",
    "",
    "Every committed artifact is listed once. The heading is the artifact path,",
    "so a plain-text search for the file name lands on its generator, code and",
    "consumers.",
    ""
  )
  for (a in artifacts) {
    lines <- c(
      lines,
      paste0("## ", normalize_slashes(a$path)),
      "",
      paste0("- kind: `", a$kind, "`"),
      paste0("- status: `", a$status, "`"),
      paste0("- artifact: [", a$path, "](", rel_from_generated(a$path), ")"),
      paste0("- generator: [", a$generator, "](", rel_from_generated(a$generator), ")"),
      paste0("- code: ", link_list(a$code)),
      paste0("- inputs: ", link_list(artifact_inputs(a))),
      paste0("- docs: ", link_list(a$docs)),
      paste0("- output_hash: `", a$output_hash, "`"),
      ""
    )
  }
  lines
}

write_projection <- function(artifacts, path = ARTIFACT_PROJECTION) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(projection_lines(artifacts), path, useBytes = TRUE)
  invisible(path)
}

validate_structure <- function(a) {
  missing <- ARTIFACT_FIELDS[!vapply(ARTIFACT_FIELDS, function(f) {
    !is.null(a[[f]])
  }, logical(1))]
  if (length(missing)) {
    stop("manifest entry ", a$path %||% "<no path>", " is missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

refresh_hashes <- function(artifacts) {
  lapply(artifacts, function(a) {
    validate_structure(a)
    a$input_hash <- artifact_input_hash(a)
    a$output_hash <- artifact_output_hash(a)
    a
  })
}

cmd_write <- function() {
  artifacts <- read_manifest()
  artifacts <- refresh_hashes(artifacts)
  write_manifest(artifacts)
  write_projection(artifacts)
  message(
    "wrote ", ARTIFACT_MANIFEST, " and ", ARTIFACT_PROJECTION,
    " (", length(artifacts), " artifacts)"
  )
  invisible(TRUE)
}

cmd_check <- function() {
  artifacts <- read_manifest()
  problems <- character()
  note <- function(...) problems <<- c(problems, paste0(...))

  expected <- projection_lines(artifacts)
  actual <- if (file.exists(ARTIFACT_PROJECTION)) {
    readLines(ARTIFACT_PROJECTION, warn = FALSE)
  } else {
    character()
  }
  if (!identical(expected, actual)) {
    note(ARTIFACT_PROJECTION, " is stale - run Rscript dev/artifacts.R --write")
  }

  for (a in artifacts) {
    validate_structure(a)
    if (!file.exists(a$path)) {
      note(a$path, ": missing artifact")
      next
    }
    if (!file.exists(a$generator)) note(a$path, ": missing generator ", a$generator)
    for (p in c(artifact_inputs(a), a$code, a$docs)) {
      if (!file.exists(p)) note(a$path, ": missing source/doc ", p)
    }
    if (!identical(a$status, "unverified")) {
      if (!identical(artifact_output_hash(a), a$output_hash)) {
        note(a$path, ": output hash changed (artifact edited or lost)")
      }
      if (file.exists(a$generator) && !identical(artifact_input_hash(a), a$input_hash)) {
        note(a$path, ": input hash changed (generator/inputs edited) - stale")
      }
    }
  }

  if (length(problems)) {
    message("artifacts check FAILED:")
    for (p in problems) message("  - ", p)
    quit(status = 1L)
  }
  message("ok: ", length(artifacts), " artifacts fresh, projection current")
  invisible(TRUE)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else {
    getwd()
  }
  setwd(root)
  if ("--write" %in% args) {
    cmd_write()
  } else if ("--check" %in% args) {
    cmd_check()
  } else {
    message("usage: Rscript dev/artifacts.R --write | --check")
    quit(status = 2L)
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

if (sys.nframe() == 0L) main()
