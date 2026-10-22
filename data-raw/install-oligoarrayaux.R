# =============================================================================
#  data-raw/install-oligoarrayaux.R  --  build OligoArrayAux from source
# =============================================================================
#  DECIPHER's thermodynamics shells out to OligoArrayAux (`hybrid-min` for
#  duplexes, `hybrid-ss-min` for hairpins).  It is a small C program, so we
#  build it from source into the git-ignored `.tools/`, laid out exactly as
#  `dev/lecture/render.sh` expects:
#
#      .tools/bin/hybrid-min            <- put on PATH
#      .tools/share/oligoarrayaux/      <- the UNAFOLDDAT data tables
#
#  No conda, no distro package.  Deliberately base-R only (no `here`/`fs`), so
#  it can run before `renv::restore()` has produced any library.
#
#  Run:  Rscript data-raw/install-oligoarrayaux.R
#  Re-running is cheap: a working build is smoke-tested and left alone.
# =============================================================================

version <- "3.8.1"
url <- sprintf("https://www.unafold.org/download/oligoarrayaux-%s.tar.gz", version)
sha256 <- "b635bcf8c3e8bc154f3b2db8c0147e0ef22bb17bcc4ed8c2c23cc3c6d2d1a5f2"

# The package root: this file lives in <root>/data-raw/, so go two levels up.
# (Falls back to the working directory when sourced, not run, by Rscript.)
script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()

tools <- file.path(root, ".tools")
src <- file.path(tools, sprintf("oligoarrayaux-%s", version))
tar <- file.path(tools, sprintf("oligoarrayaux-%s.tar.gz", version))
bin <- file.path(tools, "bin")
dat <- file.path(tools, "share", "oligoarrayaux")

say <- function(...) message("oligoarrayaux: ", ...)

# A working install runs hybrid-min on a known duplex and prints one number.
smoke <- function() {
  exe <- file.path(bin, "hybrid-min")
  if (!file.exists(exe)) {
    return(FALSE)
  }
  out <- suppressWarnings(system2(
    exe, c(
      "-n", "DNA", "-t", "37", "-T", "37", "-N", "0.05", "-E", "-q",
      "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"
    ),
    env = sprintf("UNAFOLDDAT=%s", shQuote(dat)), stdout = TRUE, stderr = FALSE
  ))
  length(out) == 1L && is.finite(suppressWarnings(as.numeric(out)))
}

run <- function(command, args = character()) {
  status <- system2(command, args)
  if (!identical(as.integer(status), 0L)) {
    stop(sprintf("'%s' failed (status %s)", command, status), call. = FALSE)
  }
}

with_dir <- function(dir, code) {
  old <- setwd(dir)
  on.exit(setwd(old))
  force(code)
}

# Try the usual sha256 tools: sha256sum (coreutils), shasum (macOS),
# openssl dgst -sha256 (fallback).  Returns NA when none is available.
sha256_of <- function(path) {
  probes <- list(
    c("sha256sum", shQuote(path)),
    c("shasum", "-a", "256", shQuote(path)),
    c("openssl", "dgst", "-sha256", shQuote(path))
  )
  for (probe in probes) {
    if (!nzchar(Sys.which(probe[[1L]]))) next
    out <- suppressWarnings(tryCatch(
      system2(probe[[1L]], probe[-1L], stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    ))
    # sha256sum/shasum print "<hex>  path"; openssl prints "SHA2-256(path)= <hex>"
    tok <- unlist(regmatches(out, gregexpr("[0-9a-f]{64}", out)))
    if (length(tok)) return(tok[[1L]])
  }
  NA_character_
}

main <- function() {
  if (smoke()) {
    say("already built and working")
    return(invisible())
  }

  if (!nzchar(Sys.which("make"))) {
    stop("'make' not found on PATH (Windows: install Rtools)", call. = FALSE)
  }
  if (!any(nzchar(Sys.which(c("cc", "gcc", "clang"))))) {
    stop("no C compiler (cc/gcc/clang) found on PATH (Windows: install Rtools)", call. = FALSE)
  }

  dir.create(tools, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(tar)) {
    say("downloading ", url)
    utils::download.file(url, tar, mode = "wb", quiet = FALSE)
  }
  got <- sha256_of(tar)
  if (is.na(got)) {
    say("no sha256sum available - skipping checksum")
  } else if (!identical(tolower(got), sha256)) {
    stop(sprintf(
      "checksum mismatch for %s\n  expected %s\n  got      %s",
      tar, sha256, got
    ), call. = FALSE)
  }

  if (!dir.exists(src)) {
    say("unpacking ", basename(tar))
    utils::untar(tar, exdir = tools)
  }

  say("configuring + building (this takes a minute)")
  with_dir(src, {
    run("./configure", sprintf("--prefix=%s", shQuote(tools)))
    run("make")
    run("make", "install")
  })

  if (!smoke()) stop("build finished but hybrid-min did not return a dG", call. = FALSE)
  say("done")
}

main()

# Point the session (or the shell) at the result:
#   export PATH="$PWD/.tools/bin:$PATH"
#   export UNAFOLDDAT="$PWD/.tools/share/oligoarrayaux"
say(sprintf("PATH       += %s", bin))
say(sprintf("UNAFOLDDAT  = %s", dat))
