# =============================================================================
#  kasp_tools.R  --  locating and installing the external command-line tools
# =============================================================================
#  Two optional, non-R dependencies back parts of the pipeline:
#
#    * OligoArrayAux (`hybrid-min`/`hybrid-ss-min`) - the dG engine DECIPHER
#      drives; needs its data tables in `UNAFOLDDAT`.
#    * NCBI BLAST+ (`blastn`/`makeblastdb`) - the indexed genome scan (rBLAST).
#
#  They are found the same way everywhere: an explicit option/path first, then a
#  project-scoped `.tools/` (dev checkouts), then the user cache
#  (`kasp_cache_dir()/tools`, where `kasp_install_tools()` builds them), then the
#  system `PATH`.  Callers that hand the tool name to another library (DECIPHER,
#  rBLAST) wrap the call in `with_oligoarrayaux()`/`with_blast()` so the resolved
#  `PATH`/`UNAFOLDDAT` apply only for the duration of that call.
# =============================================================================

#' Find an executable in a tool directory (as `<dir>/bin/<exe>` or `<dir>/<exe>`)
#' @keywords internal
#' @noRd
find_exe <- function(dir, names) {
  for (sub in c("bin", "")) {
    for (nm in names) {
      p <- fs::path(dir, sub, nm)
      if (fs::file_exists(p)) {
        return(fs::path(p))
      }
    }
  }
  NULL
}

#' Project-scoped `.tools/` directories, from the working directory upward
#' @keywords internal
#' @noRd
project_tools_dirs <- function() {
  out <- character()
  d <- getwd()
  repeat {
    out <- c(out, fs::path(d, ".tools"))
    up <- fs::path_dir(d)
    if (identical(up, d)) break
    d <- up
  }
  out
}

#' User cache `tools/` directory (created by `kasp_install_tools()`)
#' @keywords internal
#' @noRd
cache_tools_dir <- function() fs::path(tools::R_user_dir("tidykasp", "cache"), "tools")

#' Candidate tool directories, in resolution order (option first, `PATH` last)
#' @keywords internal
#' @noRd
tool_dirs <- function(opt) {
  o <- getOption(opt)
  c(if (!is.null(o)) fs::path(o), project_tools_dirs(), cache_tools_dir())
}

#' Resolve OligoArrayAux (executable + data tables)
#'
#' @return `NULL` when not found, else `list(exe, bin, dat)` (`dat` is `NA` when
#'   the `UNAFOLDDAT` tables were not found next to the executable).
#' @keywords internal
#' @noRd
oligoarrayaux_paths <- function() {
  exe <- NULL
  for (d in tool_dirs("tidykasp.oligoarrayaux")) {
    exe <- find_exe(d, c("hybrid-min", "hybrid-min.exe"))
    if (!is.null(exe)) break
  }
  if (is.null(exe)) {
    p <- Sys.which("hybrid-min")
    if (!nzchar(p)) {
      return(NULL)
    }
    exe <- fs::path(p)
  }
  prefix <- fs::path_dir(exe)
  if (identical(fs::path_file(prefix), "bin")) prefix <- fs::path_dir(prefix)
  dat <- fs::path(prefix, "share", "oligoarrayaux")
  list(exe = exe, bin = fs::path_dir(exe), dat = if (fs::dir_exists(dat)) dat else NA_character_)
}

#' Run `code` with OligoArrayAux on `PATH` and `UNAFOLDDAT` set (when found)
#'
#' DECIPHER calls `hybrid-min` by name, so a resolved build must be exposed both
#' ways - only for the duration of the wrapped call.
#' @keywords internal
#' @noRd
with_oligoarrayaux <- function(code) {
  p <- oligoarrayaux_paths()
  if (is.null(p)) {
    return(force(code))
  }
  env <- c(PATH = paste(p$bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
  if (!is.na(p$dat)) env <- c(env, UNAFOLDDAT = p$dat)
  withr::with_envvar(env, code)
}

#' Resolve the NCBI BLAST+ `bin` directory
#' @return `NULL` when both `blastn` and `makeblastdb` are not found together.
#' @keywords internal
#' @noRd
blast_bin <- function() {
  for (d in tool_dirs("tidykasp.blast")) {
    b <- find_exe(d, c("blastn", "blastn.exe"))
    m <- find_exe(d, c("makeblastdb", "makeblastdb.exe"))
    if (!is.null(b) && !is.null(m)) {
      return(fs::path_dir(b))
    }
  }
  NULL
}

#' Run `code` with a resolved BLAST+ `bin` directory on `PATH` (when found)
#' @keywords internal
#' @noRd
with_blast <- function(code) {
  d <- blast_bin()
  if (is.null(d)) {
    return(force(code))
  }
  withr::with_path(d, code)
}

#' sha256 of a file via the first available system tool (`NA` when none)
#' @keywords internal
#' @noRd
sha256_of <- function(path) {
  probes <- list(
    c("sha256sum", path),
    c("shasum", "-a", "256", path),
    c("openssl", "dgst", "-sha256", path)
  )
  for (probe in probes) {
    if (!nzchar(Sys.which(probe[[1L]]))) next
    out <- suppressWarnings(tryCatch(
      system2(probe[[1L]], probe[-1L], stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    ))
    tok <- purrr::list_c(stringr::str_extract_all(out, "[0-9a-f]{64}"))
    if (length(tok)) {
      return(tok[[1L]])
    }
  }
  NA_character_
}

#' Smoke-test an OligoArrayAux `bin`/`dat` pair (one duplex, one number)
#' @keywords internal
#' @noRd
oligoarrayaux_works <- function(bin, dat) {
  exe <- fs::path(bin, if (.Platform$OS.type == "windows") "hybrid-min.exe" else "hybrid-min")
  if (!fs::file_exists(exe)) {
    return(FALSE)
  }
  out <- suppressWarnings(system2(
    exe,
    c(
      "-n", "DNA", "-t", "37", "-T", "37", "-N", "0.05", "-E", "-q",
      "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"
    ),
    env = sprintf("UNAFOLDDAT=%s", shQuote(dat)), stdout = TRUE, stderr = FALSE
  ))
  length(out) == 1L && is.finite(suppressWarnings(as.numeric(out)))
}

#' @keywords internal
#' @noRd
run_step <- function(command, args = character(), dir = NULL) {
  status <- if (is.null(dir)) system2(command, args) else withr::with_dir(dir, system2(command, args))
  if (!identical(as.integer(status), 0L)) {
    rlang::abort(glue::glue(gettext("'{command}' failed (status {status}).")))
  }
  invisible()
}

#' Build OligoArrayAux into the user cache
#'
#' Downloads the OligoArrayAux source release, verifies its checksum, builds it
#' and installs it under `dest` (`kasp_cache_dir()/tools` by default) - the
#' directory [oligoarrayaux_paths()] then searches, so the thermodynamics and QC
#' functions start working without any `PATH`/`UNAFOLDDAT` fiddling.
#'
#' Nothing is downloaded or written unless you call this function; a working
#' build is detected and left alone.  On Windows the build runs through `bash`
#' (Rtools/MSYS2), which must be on the `PATH`.
#'
#' @param dest directory to install into (default: the user cache).
#' @param version OligoArrayAux release to build.
#' @param ask ask for consent before downloading when `interactive()`.
#' @param force rebuild even when a working build already exists.
#' @return `dest`, invisibly.
#' @examples
#' \dontrun{
#' kasp_install_tools()
#' has_oligoarrayaux()
#' }
#' @family KASP utilities
#' @export
kasp_install_tools <- function(dest = fs::path(kasp_cache_dir(), "tools"),
                               version = "3.8.1", ask = interactive(), force = FALSE) {
  dest <- fs::path(dest)
  bin <- fs::path(dest, "bin")
  dat <- fs::path(dest, "share", "oligoarrayaux")
  if (!force && oligoarrayaux_works(bin, dat)) {
    cli::cli_alert_success(gettext("OligoArrayAux is already installed in {.path {dest}}."))
    return(invisible(dest))
  }
  if (ask && !isTRUE(utils::askYesNo(gettext("Download and build OligoArrayAux into the user cache?")))) {
    cli::cli_alert_info(gettext("Cancelled; nothing was downloaded."))
    return(invisible(dest))
  }
  make <- if (.Platform$OS.type == "windows") Sys.which("make") else Sys.which("make")
  if (!nzchar(make)) {
    rlang::abort(c(
      gettext("'make' was not found on the PATH."),
      i = gettext("On Windows install Rtools; on Unix install your build tools.")
    ))
  }
  if (!any(nzchar(Sys.which(c("cc", "gcc", "clang"))))) {
    rlang::abort(c(
      gettext("No C compiler (cc/gcc/clang) was found on the PATH."),
      i = gettext("On Windows install Rtools; on Unix install a compiler.")
    ))
  }

  url <- glue::glue("https://www.unafold.org/download/oligoarrayaux-{version}.tar.gz")
  sha256 <- "b635bcf8c3e8bc154f3b2db8c0147e0ef22bb17bcc4ed8c2c23cc3c6d2d1a5f2"
  fs::dir_create(dest)
  tar <- fs::path(dest, glue::glue("oligoarrayaux-{version}.tar.gz"))
  if (!fs::file_exists(tar)) {
    cli::cli_alert_info(gettext("Downloading OligoArrayAux {version} ..."))
    utils::download.file(url, tar, mode = "wb", quiet = TRUE)
  }
  got <- sha256_of(tar)
  if (is.na(got)) {
    cli::cli_alert_warning(gettext("No sha256 tool found; skipping the checksum."))
  } else if (!identical(tolower(got), sha256)) {
    rlang::abort(glue::glue(gettext("Checksum mismatch for {tar} (expected {sha256}, got {got}).")))
  }

  src <- fs::path(dest, glue::glue("oligoarrayaux-{version}"))
  if (!fs::dir_exists(src)) utils::untar(tar, exdir = dest)

  cli::cli_alert_info(gettext("Configuring and building OligoArrayAux (this takes a minute) ..."))
  if (.Platform$OS.type == "windows") {
    bash <- Sys.which("bash")
    if (!nzchar(bash)) {
      rlang::abort(c(
        gettext("Building on Windows needs 'bash' (Rtools/MSYS2) on the PATH."),
        i = gettext("Install Rtools, or build OligoArrayAux by hand.")
      ))
    }
    cmd <- sprintf(
      "cd %s && ./configure --prefix=%s && make && make install",
      shQuote(src), shQuote(dest)
    )
    run_step(bash, c("-lc", shQuote(cmd)))
  } else {
    run_step("./configure", sprintf("--prefix=%s", shQuote(dest)), dir = src)
    run_step("make", dir = src)
    run_step("make", "install", dir = src)
  }

  if (!oligoarrayaux_works(bin, dat)) {
    rlang::abort(gettext("The build finished but hybrid-min did not return a dG."))
  }
  cli::cli_alert_success(gettext("OligoArrayAux installed in {.path {dest}}."))
  invisible(dest)
}
