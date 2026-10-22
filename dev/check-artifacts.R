# =============================================================================
#  dev/check-artifacts.R  --  the artifacts freshness gate
# =============================================================================
#  Validates docs/artifacts.yml against the tree:
#
#    1. every manifest path, generator, input, code and doc exists
#    2. every committed artifact is in the manifest (no orphan)
#    3. the input hash is unchanged (generator + inputs)
#    4. the output hash matches the artifact on disk
#    5. docs/generated/artifacts.md is current (the projection is not stale)
#
#  Only tracked files under man/figures and vignettes/figures count as
#  artifacts; vignettes/figures/*.png and the other locally rendered rasters are
#  git-ignored, so they are not orphans.  A stale hash or an orphan exits
#  non-zero.  Entries with `status: unverified` are reported but their hashes
#  are not enforced (no generator could be justified).
#
#  Run:  Rscript dev/check-artifacts.R    (exit 1 on stale/orphan)
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

artifacts <- read_manifest()
paths <- vapply(artifacts, function(a) a$path, character(1))

# 1. structure and paths -----------------------------------------------------
for (a in artifacts) {
  tryCatch(validate_structure(a), error = function(e) note(conditionMessage(e)))
  if (!file.exists(a$path)) note(a$path, ": missing artifact")
  if (!file.exists(a$generator)) note(a$path, ": missing generator ", a$generator)
  for (p in unique(c(artifact_inputs(a), a$code, a$docs))) {
    if (!file.exists(p)) note(a$path, ": missing source/doc ", p)
  }
}

# 2. orphan artifacts (tracked files in the artifact directories) ------------
artifact_ext <- "[.](png|svg|jpe?g|webp|pdf|gif)$"
tracked <- tryCatch(
  suppressWarnings(
    system2(
      "git",
      c("ls-files", "--", "man/figures", "vignettes/figures"),
      stdout = TRUE, stderr = FALSE
    )
  ),
  error = function(e) NULL
)
if (!is.null(tracked) &&
  (is.null(attr(tracked, "status")) || identical(attr(tracked, "status"), 0L))) {
  candidates <- tracked[grepl(artifact_ext, tracked, ignore.case = TRUE)]
} else {
  # no git: fall back to the directory scan, skipping the git-ignored rasters
  candidates <- c(
    list.files("man/figures", pattern = artifact_ext, full.names = TRUE),
    list.files("vignettes/figures", pattern = artifact_ext, full.names = TRUE)
  )
  candidates <- candidates[!grepl("^vignettes/figures/.*[.]png$", candidates)]
}
orphans <- setdiff(normalize_slashes(candidates), paths)
if (length(orphans)) {
  for (o in orphans) note(o, ": orphan artifact (not in docs/artifacts.yml)")
} else {
  ok(paste0(length(candidates), " committed artifacts are all in the manifest"))
}
if (dir.exists("vignettes/figures")) {
  local_only <- list.files("vignettes/figures", pattern = "[.]png$", full.names = TRUE)
  if (length(local_only)) {
    message("  note ", length(local_only), " vignettes/figures/*.png are local-only (git-ignored); not checked")
  }
}

# 3 + 4. hashes --------------------------------------------------------------
for (a in artifacts) {
  if (!file.exists(a$path) || !file.exists(a$generator)) next
  out_now <- artifact_output_hash(a)
  in_now <- artifact_input_hash(a)
  if (identical(a$status, "unverified")) {
    message("  warn ", a$path, ": status unverified - hashes reported, not enforced")
    next
  }
  if (!identical(out_now, a$output_hash)) {
    note(
      a$path, ": output hash mismatch (artifact edited; expected ",
      a$output_hash, ", got ", out_now, ")"
    )
  }
  if (!identical(in_now, a$input_hash)) {
    note(a$path, ": input hash mismatch (generator/inputs changed; stale)")
  }
}
if (!length(failures)) ok(paste0(length(artifacts), " artifacts have matching input and output hashes"))

# 5. projection freshness ----------------------------------------------------
expected <- projection_lines(artifacts)
actual <- if (file.exists(ARTIFACT_PROJECTION)) {
  readLines(ARTIFACT_PROJECTION, warn = FALSE)
} else {
  character()
}
if (!identical(expected, actual)) {
  note(ARTIFACT_PROJECTION, " is stale - run Rscript dev/artifacts.R --write")
} else {
  ok(ARTIFACT_PROJECTION, " is current")
}

if (length(failures)) {
  message("artifacts check FAILED:")
  for (f in failures) message("  - ", f)
  quit(status = 1L)
}
message("ok: artifacts graph is fresh")
