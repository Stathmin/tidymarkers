# Activate a resolved OligoArrayAux for the thermodynamics tests, so they run
# instead of skipping.  Same resolution order as the package (see R/kasp_tools.R):
# a project-scoped `.tools/` (dev checkout) or the user cache built by
# `kasp_install_tools()`.  No-op when neither exists - a plain CRAN check, or a
# machine without the tool - and those tests then skip.

local({
  has_exe <- function(bin) {
    file.exists(file.path(bin, "hybrid-min")) || file.exists(file.path(bin, "hybrid-min.exe"))
  }
  prefixes <- function() {
    out <- character()
    d <- getwd()
    repeat {
      out <- c(out, file.path(d, ".tools"))
      up <- dirname(d)
      if (identical(up, d)) break
      d <- up
    }
    cache <- tryCatch(tools::R_user_dir("tidykasp", "cache"), error = function(e) "")
    c(out, if (nzchar(cache)) file.path(cache, "tools"))
  }
  for (p in prefixes()) {
    bin <- file.path(p, "bin")
    if (has_exe(bin)) {
      Sys.setenv(
        PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep),
        UNAFOLDDAT = file.path(p, "share", "oligoarrayaux")
      )
      break
    }
  }
})

# Skip a thermo/QC test when OligoArrayAux cannot be found (options/.tools/cache/PATH).
skip_if_no_oligoarrayaux <- function() {
  if (!has_oligoarrayaux()) testthat::skip("OligoArrayAux not available")
}
