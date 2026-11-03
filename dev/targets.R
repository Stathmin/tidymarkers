# =============================================================================
#  dev/targets.R  --  run the freshness pipeline (ZA8.1)
# =============================================================================
#  A thin front end over dev/targets/_targets.R, so the store path lives in one
#  place.  It is an alternative to `dev/render-all.R`, not a replacement: both
#  call the same renderer, and whichever runs, the figure hashes are the same.
#
#  Modes:
#    --make      render only what changed (tar_make)
#    --outdated  list what a fresh clone would consider stale
#    --check     exit 1 if anything is outdated (the gate's use)
#
#  `targets` is a Suggests: without it this prints a message and exits 0, so
#  the gate falls back to dev/render-all.R on a bare machine.
#
#  Run:  Rscript dev/targets.R --make | --outdated | --check
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

args <- commandArgs(trailingOnly = TRUE)
action <- if ("--make" %in% args) {
  "make"
} else if ("--outdated" %in% args) {
  "outdated"
} else if ("--check" %in% args) {
  "check"
} else {
  NA
}
if (is.na(action)) {
  message("usage: Rscript dev/targets.R --make | --outdated | --check")
  quit(status = 2L)
}

if (!requireNamespace("targets", quietly = TRUE)) {
  message("skip (targets not installed; dev/render-all.R is the fallback)")
  quit(status = 0L)
}

TARGET_SCRIPT <- "dev/targets/_targets.R"
STORE <- "dev/targets/_targets"

if (identical(action, "make")) {
  targets::tar_make(script = TARGET_SCRIPT, store = STORE, callr_function = NULL)
  quit(status = 0L)
}

outdated <- targets::tar_outdated(script = TARGET_SCRIPT, store = STORE, callr_function = NULL)
if (!length(outdated)) {
  message("ok: no target is outdated")
  quit(status = 0L)
}
message("outdated: ", paste(outdated, collapse = ", "))
if (identical(action, "check")) {
  message("FAIL run `Rscript dev/targets.R --make` and commit the result")
  quit(status = 1L)
}
quit(status = 0L)
