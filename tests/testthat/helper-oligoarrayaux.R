# Activate a locally built OligoArrayAux (see data-raw/install-oligoarrayaux.R)
# so the thermodynamics tests run instead of skipping.  No-op when the tool is
# absent - a plain CRAN check, or a checkout that never built it - and those
# tests then skip exactly as before.
local({
  find_tools <- function(dir) {
    repeat {
      bin <- file.path(dir, ".tools", "bin")
      if (file.exists(file.path(bin, "hybrid-min")) || file.exists(file.path(bin, "hybrid-min.exe"))) {
        return(file.path(dir, ".tools"))
      }
      parent <- dirname(dir)
      if (identical(parent, dir)) {
        return(NULL)
      }
      dir <- parent
    }
  }
  tools <- find_tools(getwd())
  if (is.null(tools)) {
    return(invisible())
  }
  Sys.setenv(
    PATH = paste(file.path(tools, "bin"), Sys.getenv("PATH"), sep = .Platform$path.sep),
    UNAFOLDDAT = file.path(tools, "share", "oligoarrayaux")
  )
})
