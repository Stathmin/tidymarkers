# =============================================================================
#  dev/check-examples-noengine.R  --  run every Rd example without the engine
# =============================================================================
#  The macOS CI leg installs BLAST+ only, so OligoArrayAux is absent there and
#  `pcr_efficiency()` aborts.  R CMD check on macOS therefore stops on examples,
#  one Rd file at a time; a host with the engine (this one) never sees it, and
#  `R CMD check` cannot emulate the leg.  This script does: it hides
#  `oligoarrayaux_paths()` from the namespace and sources every `man/*.Rd`
#  example exactly as `tools::Rd2ex()` writes them.  An example that reaches the
#  engine is a failure -- wrap that tail in `\dontrun{}` (`\donttest{}` does not
#  help: R CMD check runs it).
#
#  Run:  Rscript dev/check-examples-noengine.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

suppressMessages(pkgload::load_all(".", helpers = FALSE))
ns <- asNamespace("tidymarkers")
orig <- get("oligoarrayaux_paths", envir = ns)
assignInNamespace("oligoarrayaux_paths", function() NULL, ns = "tidymarkers")
on.exit(assignInNamespace("oligoarrayaux_paths", orig, ns = "tidymarkers"), add = TRUE)

out <- file.path(tempdir(), "rd-examples")
dir.create(out, showWarnings = FALSE)
rds <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
fails <- character()
for (rd in rds) {
  ex <- file.path(out, sub("\\.Rd$", ".R", basename(rd)))
  ok <- tryCatch(
    {
      tools::Rd2ex(rd, out = ex)
      TRUE
    },
    error = function(e) FALSE
  )
  if (!ok || !file.exists(ex)) next
  err <- tryCatch(
    {
      suppressMessages(sys.source(ex, envir = new.env(parent = globalenv())))
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (!is.null(err)) {
    fails <- c(fails, sprintf("%s: %s", tools::file_path_sans_ext(basename(rd)), gsub("\n", " ", err)))
  }
}
if (length(fails)) {
  stop(
    length(fails), " of ", length(rds), " examples need the engine and are not in \\dontrun{}:\n  ",
    paste(fails, collapse = "\n  "),
    call. = FALSE
  )
}
message("  ok   all ", length(rds), " examples run without OligoArrayAux")
