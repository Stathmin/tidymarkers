# =============================================================================
#  dev/step-cost.R  --  ZX5.1: what each artifact regeneration step actually costs
# =============================================================================
#  ADR 0002 widened `targets` from "wait for a second consumer" to "run one
#  measured pass over the artifact chain and then decide".  This script is that
#  pass.  It does not change the pipeline; it prices every step the gate and a
#  human already run, so the next widening is argued from numbers.
#
#  Scope is the **artifact chain** only.  The `data-raw` half of ZX5.1 is not
#  measured here on purpose: 56 of the 97 `data-raw/*.R` scripts read an
#  operator-held source outside the tree (`.tmp/...`, an Excel supplement, a
#  fetched genome) and `stop()` with a recipe when it is absent, so a timing run
#  would measure the presence of that source, not the pipeline.  See ROADMAP
#  ZX5.1 for the count and the consequence.
#
#  One child Rscript per step so each gets its own VmHWM (`proc-mem.R`).  The
#  child is a wrapper that shims `commandArgs()` so the step sees exactly what
#  `Rscript <step> <args>` gives it -- `--file=` root discovery included -- and
#  prints its peak through `.Last` as well, because these scripts may `quit()`.
#  `elapsed_s` is the whole child: R startup, `pkgload::load_all()` and the work.
#  That is the number a human pays, not the renderer alone.
#
#  `targets --make` is measured twice: once from a deleted store (a cold clone)
#  and once immediately after (a warm store).  The gap is the value of the
#  content-hash cache.  The store is `dev/targets/_targets/`, git-ignored and
#  fully regenerable; this script deletes it and rebuilds it.
#
#  Writes results/step-cost.csv (results/ is git-ignored; the numbers that must
#  survive in the tree go into ROADMAP.md ZX5.1 and ADR 0002).
#
#  Run:  Rscript dev/step-cost.R [--check]
#    --check  exit non-zero if the CSV is missing, short of a step, or malformed
# =============================================================================

# ---- plumbing ----------------------------------------------------------------

step_root <- function() {
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))))
  }
  getwd()
}

root <- step_root()
setwd(root)
source(file.path(root, "dev", "proc-mem.R"))

argv <- commandArgs(trailingOnly = TRUE)
out_csv <- file.path(root, "results", "step-cost.csv")
store <- file.path(root, "dev", "targets", "_targets")
columns <- c("step", "status", "elapsed_s", "peak_rss_mb", "rerun", "cause")

# `rerun` is how often a human or the gate starts the step, which is what makes
# a step's cost matter: "gate" = every `Rscript dev/check.R`, "change" = after
# touching the engine or a data file, "rare" = release or a deliberate redraw.
steps <- list(
  list(
    name = "targets-make-cold", script = "dev/targets.R", args = "--make",
    rerun = "rare", prepare = "clear_store"
  ),
  list(name = "targets-make-warm", script = "dev/targets.R", args = "--make", rerun = "gate"),
  list(name = "targets-check", script = "dev/targets.R", args = "--check", rerun = "change"),
  list(name = "render-decision", script = "dev/render-decision.R", args = character(),
    rerun = "change"),
  list(name = "render-gallery", script = "dev/render-gallery.R", args = character(),
    rerun = "change"),
  list(name = "render-plots", script = "dev/render-plots.R", args = character(),
    rerun = "change"),
  list(name = "artifacts-check", script = "dev/artifacts.R", args = "--check", rerun = "gate"),
  list(name = "artifacts-write", script = "dev/artifacts.R", args = "--write", rerun = "change")
)

write_csv <- function(tab, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(tab, path, row.names = FALSE, na = "")
}

# ---- check mode --------------------------------------------------------------

if ("--check" %in% argv) {
  if (!file.exists(out_csv)) {
    stop("missing ", out_csv, "; run without --check first", call. = FALSE)
  }
  tab <- utils::read.csv(out_csv, stringsAsFactors = FALSE)
  missing_cols <- setdiff(columns, names(tab))
  if (length(missing_cols)) {
    stop("results CSV is missing columns: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  absent <- setdiff(vapply(steps, `[[`, character(1), "name"), tab$step)
  if (length(absent)) {
    stop("results CSV has no row for: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  measured <- tab$status == "measured"
  if (any(!is.finite(tab$elapsed_s[measured])) || any(!is.finite(tab$peak_rss_mb[measured]))) {
    stop("a measured row has a non-finite elapsed_s or peak_rss_mb", call. = FALSE)
  }
  bad <- tab[!measured & is.na(tab$cause), ]
  if (nrow(bad)) stop(nrow(bad), " non-measured row(s) without a cause", call. = FALSE)
  cat("results/step-cost.csv: ", nrow(tab), " steps, ", sum(measured), " measured\n",
    sep = ""
  )
  quit(save = "no", status = 0L)
}

# ---- child wrapper -----------------------------------------------------------

# The step scripts discover the package root from `--file=` and branch on
# `commandArgs(trailingOnly = TRUE)`.  A wrapper that `sys.source()`s the script
# would give it the wrapper's own argv, so it gets a shim instead: the shim
# answers both calls exactly as `Rscript <step> <args>` would.
wrapper <- file.path(tempdir(), "step-cost-wrapper.R")
writeLines(c(
  "argv <- commandArgs(trailingOnly = TRUE)",
  "proc_mem <- argv[[1L]]",
  "step_script <- argv[[2L]]",
  "step_args <- if (length(argv) > 2L) argv[-c(1L, 2L)] else character()",
  "shim <- function(trailingOnly = FALSE) {",
  "  if (isTRUE(trailingOnly)) return(step_args)",
  "  c('Rscript', paste0('--file=', step_script),",
  "    if (length(step_args)) c('--args', step_args))",
  "}",
  "assign('commandArgs', shim, envir = globalenv())",
  "sys.source(proc_mem, envir = globalenv())",
  "peak <- function() {",
  "  cat('\\nSTEP_PEAK_KB=', format(proc_peak_rss_kb(), scientific = FALSE), '\\n', sep = '')",
  "  flush(stdout())",
  "}",
  "assign('.Last', peak, envir = globalenv())",
  "env <- new.env(parent = globalenv())",
  "sys.source(step_script, envir = env, keep.source = FALSE)",
  "peak()"
), wrapper)

rscript <- file.path(R.home("bin"), "Rscript")
work <- file.path(tempdir(), "step-cost")
dir.create(work, showWarnings = FALSE, recursive = TRUE)

run_step <- function(step) {
  if (identical(step$prepare, "clear_store")) {
    unlink(store, recursive = TRUE)
  }
  log <- file.path(work, paste0(step$name, ".log"))
  err <- file.path(work, paste0(step$name, ".err"))
  t0 <- proc.time()[["elapsed"]]
  status <- system2(
    rscript,
    c(
      "--vanilla", shQuote(wrapper), shQuote(file.path(root, "dev", "proc-mem.R")),
      shQuote(file.path(root, step$script)), step$args
    ),
    stdout = log, stderr = err
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  peaks <- suppressWarnings(as.numeric(sub(
    "^STEP_PEAK_KB=\\s*(-?[0-9]+)$", "\\1",
    grep("^STEP_PEAK_KB=", readLines(log, warn = FALSE), value = TRUE)
  )))
  peaks <- peaks[is.finite(peaks)]
  peak_mb <- if (length(peaks)) max(peaks) / 1024 else NA_real_

  ok <- identical(as.integer(status), 0L) && length(peaks) > 0L
  cause <- NA_character_
  if (!ok) {
    lines <- readLines(err, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    hit <- grep("^Error", lines, value = TRUE)
    cause <- trimws(if (length(hit)) hit[[length(hit)]] else tail(lines, 1L))
    if (!length(cause) || !nzchar(cause)) {
      cause <- paste0("child status ", status, ", no error text")
    }
  }
  data.frame(
    step = step$name,
    status = if (ok) "measured" else "failed",
    elapsed_s = elapsed,
    peak_rss_mb = if (ok) peak_mb else NA_real_,
    rerun = step$rerun,
    cause = cause,
    stringsAsFactors = FALSE
  )
}

# ---- parent: run every step --------------------------------------------------

rows <- lapply(steps, function(step) {
  message("== ", step$name, " ==")
  row <- run_step(step)
  message(
    "  ", row$status, " in ", sprintf("%.1f s", row$elapsed_s),
    if (is.na(row$peak_rss_mb)) "" else paste0(", peak ", sprintf("%.1f MB", row$peak_rss_mb)),
    if (!is.na(row$cause)) paste0(" (", row$cause, ")") else ""
  )
  row
})
tab <- do.call(rbind, rows)
write_csv(tab, out_csv)

cat("\nresults/step-cost.csv\n")
print(tab[, c("step", "status", "elapsed_s", "peak_rss_mb", "rerun")], row.names = FALSE)

measured <- tab[tab$status == "measured", ]
if (nrow(measured)) {
  cat("\nslowest step: ", measured$step[[which.max(measured$elapsed_s)]],
    " at ", sprintf("%.1f s", max(measured$elapsed_s)), "\n",
    sep = ""
  )
  gate <- measured$elapsed_s[measured$rerun == "gate"]
  if (length(gate)) {
    cat("gate-time total (rerun = gate): ", sprintf("%.1f s", sum(gate)), "\n", sep = "")
  }
}
cold <- tab$elapsed_s[tab$step == "targets-make-cold"]
warm <- tab$elapsed_s[tab$step == "targets-make-warm"]
if (length(cold) && length(warm) && all(is.finite(c(cold, warm)))) {
  cat("cold vs warm store: ", sprintf("%.1f s -> %.1f s", cold, warm), "\n", sep = "")
}
