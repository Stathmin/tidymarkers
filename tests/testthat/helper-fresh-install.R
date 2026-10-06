# Parallel tests must run the CURRENT source, never a stale build.
#
# Under `test_local()`/`devtools::test()` the package is `pkgload::load_all()`ed,
# so the master process runs the source.  A SnowParam/SOCK worker, however, is a
# fresh R process: it resolves the package by name and loads the *installed*
# `tidykasp`.  If that install is older than the working tree, the SOCK tests
# silently exercise old code and can pass while the source is broken.
#
# So, in a dev session, refresh the install once when the source is newer.  A
# marker file records the source mtime that was last installed (R CMD INSTALL
# does not bump DESCRIPTION's mtime, so we cannot probe the install directly).
# The refresh is serialised with a lock because testthat runs files in parallel.
# If it is not possible (not a dev session, or the library is not writable),
# SOCK tests call `skip_if_stale_install()` and skip loudly rather than test a
# stale build.

kasp_dev_session <- function() {
  requireNamespace("pkgload", quietly = TRUE) &&
    !is.null(tryCatch(pkgload::pkg_path(), error = function(e) NULL))
}

# mtime of the files a build would pick up
kasp_source_mtime <- function(root) {
  files <- c(
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
    file.path(root, "DESCRIPTION"),
    file.path(root, "NAMESPACE")
  )
  files <- files[file.exists(files)]
  if (!length(files)) NA_real_ else as.numeric(max(file.info(files)$mtime, na.rm = TRUE))
}

# the installed copy: `find.package()` returns the *source* under load_all, so
# scan the libraries for a tidykasp that is not the working tree
kasp_installed_path <- function(root) {
  for (lib in .libPaths()) {
    p <- file.path(lib, "tidykasp")
    if (!file.exists(file.path(p, "DESCRIPTION"))) next
    if (identical(normalizePath(p, mustWork = FALSE), normalizePath(root, mustWork = FALSE))) next
    return(p)
  }
  NULL
}

kasp_refresh_install <- function() {
  if (!kasp_dev_session()) {
    # R CMD check installed the package from this very source tree
    options(tidykasp.tests.fresh = TRUE)
    return(invisible(TRUE))
  }
  root <- tryCatch(pkgload::pkg_path(), error = function(e) NULL)
  installed <- if (is.null(root)) NULL else kasp_installed_path(root)
  if (is.null(root) || is.null(installed)) {
    options(tidykasp.tests.fresh = FALSE)
    return(invisible(FALSE))
  }

  cache <- tools::R_user_dir("tidykasp", "cache")
  fs::dir_create(cache)
  marker <- file.path(cache, ".tests-install-mtime")
  fresh <- function() {
    if (!file.exists(marker)) {
      return(FALSE)
    }
    # +1 s absorbs the as.character() rounding of the fractional POSIX mtime
    built <- suppressWarnings(as.numeric(readLines(marker, warn = FALSE)[1L]))
    !is.na(built) && !is.na(kasp_source_mtime(root)) && kasp_source_mtime(root) <= built + 1
  }
  if (fresh()) {
    options(tidykasp.tests.fresh = TRUE)
    return(invisible(TRUE))
  }

  # serialise concurrent test processes (testthat parallelises per file)
  lock <- file.path(cache, ".tests-install.lock")
  for (i in seq_len(120L)) {
    if (dir.create(lock, showWarnings = FALSE)) break
    Sys.sleep(0.5)
    if (fresh()) { # another process refreshed it
      options(tidykasp.tests.fresh = TRUE)
      return(invisible(TRUE))
    }
  }
  if (!dir.exists(lock)) {
    options(tidykasp.tests.fresh = FALSE)
    return(invisible(FALSE))
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)

  if (!fresh()) {
    message("tidykasp: refreshing the installed build so parallel tests use the current source ...")
    # our own lock already serialises installs, so any R lock dir is stale
    unlink(file.path(.libPaths()[1L], "00LOCK-tidykasp"), recursive = TRUE)
    status <- system2(
      file.path(R.home("bin"), "R"),
      c(
        "CMD", "INSTALL", "--no-multiarch", "--no-byte-compile",
        paste0("--library=", .libPaths()[1L]), root
      ),
      stdout = FALSE, stderr = FALSE
    )
    if (identical(as.integer(status), 0L)) {
      writeLines(as.character(kasp_source_mtime(root)), marker)
    }
  }
  options(tidykasp.tests.fresh = fresh())
  invisible(fresh())
}

kasp_refresh_install()

# For tests that spin up a real second process (SnowParam/SOCK).
skip_if_stale_install <- function() {
  if (isTRUE(getOption("tidykasp.tests.fresh"))) {
    return(invisible(TRUE))
  }
  testthat::skip("the installed tidykasp could not be refreshed to the current source")
}
