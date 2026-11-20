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
  exe <- file.path(bin, exe_name("hybrid-min"))
  if (!file.exists(exe)) {
    return(FALSE)
  }
  # system2(env=) is not usable: on Windows it appends "NAME=value" to the
  # command line instead of setting the variable (R's system2.Rd), and the
  # program then reads the data path as a sequence.  Set it for the child.
  old <- Sys.getenv("UNAFOLDDAT", unset = NA_character_)
  Sys.setenv(UNAFOLDDAT = dat)
  on.exit(
    if (is.na(old)) Sys.unsetenv("UNAFOLDDAT") else Sys.setenv(UNAFOLDDAT = old),
    add = TRUE
  )
  out <- suppressWarnings(system2(
    exe, c(
      "-n", "DNA", "-t", "37", "-T", "37", "-N", "0.05", "-E", "-q",
      "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"
    ),
    stdout = TRUE, stderr = FALSE
  ))
  dg <- suppressWarnings(as.numeric(out))
  # a perfect 20-mer homoduplex is stabilising; a positive number means the
  # tables never reached the program (this is what hid the Windows bug)
  length(out) == 1L && is.finite(dg) && dg < 0
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

exe_name <- function(program) {
  if (.Platform$OS.type == "windows") paste0(program, ".exe") else program
}

compiler <- function() {
  cc <- Sys.which(c("cc", "gcc", "clang"))
  cc <- unname(cc[nzchar(cc)])
  if (!length(cc)) stop("no C compiler (cc/gcc/clang) found on PATH (Windows: install Rtools)", call. = FALSE)
  cc[[1L]]
}

# A string-valued -D: the quotes are part of the macro, and system2() runs the
# arguments through a shell, so the value is wrapped in a quoted group whose
# inner quotes are escaped for that shell (and for the Windows CRT).
define <- function(name, value) {
  sprintf("-D%s=\"\\\"%s\\\"\"", name, value)
}

# The units and flags are what autotools compiled for each program, plus the
# macros configure wrote into config.h (the tarball ships only config.in):
# without <stdlib.h>/<string.h>, gcc >= 14 rejects the implicit declarations.
build <- function(cc, program) {
  datadir <- normalizePath(dat, winslash = "/", mustWork = FALSE)
  # The tarball is K&R C: `extern int getopt ();` means "any arguments" there.
  # A compiler defaulting to C23 reads `()` as `(void)`, colliding with mingw's
  # real getopt/getenv prototypes.
  run(cc, c(
    "-std=gnu17", "-g", "-O2", "-DENERGY=int", "-DINTEGER", "-DPRECISION=1000",
    "-DHAVE_STDLIB_H=1", "-DHAVE_STRING_H=1", "-DHAVE_FCNTL_H=1",
    define("PACKAGE_NAME", "OligoArrayAux"),
    define("PACKAGE_VERSION", version),
    define("PACKAGE_BUGREPORT", "markham@alum.rpi.edu"),
    # energy.c's openFile() concatenates PKGDATADIR with the file name.
    define("PKGDATADIR", paste0(datadir, "/")),
    "-o", file.path(bin, exe_name(program)),
    "energy.c", paste0(program, ".c"), "getopt.c", "getopt1.c", "xmalloc.c", "-lm"
  ))
}

# The tarball ships no tables: autotools generates the files in data/ from
# data/rules with machine.pl, and data/Makefile.am records the flags for each
# target ($(INFINITY) where a missing pair is forbidden, 0 for enthalpies).
# 999999 is what configure's infinity probe fell back to, and what every
# published OligoArrayAux install carries, so the output does not depend on
# which compiler ran the probe.
tables <- function(perl) {
  makefile <- readLines(file.path(src, "data", "Makefile.am"))
  lists <- grep("^[a-z]+ = ", makefile, value = TRUE)
  targets <- unique(unlist(strsplit(sub("^[a-z]+ = ", "", lists), " +")))
  for (target in targets) {
    at <- which(startsWith(makefile, paste0(target, " : ")))
    if (!length(at)) next
    rule <- makefile[[at[[1L]]]]
    tokens <- strsplit(trimws(makefile[[at[[1L]] + 1L]]), " +")[[1L]]
    args <- tokens[-c(1L, 2L, length(tokens) - 1L, length(tokens))]
    args[args == "$(INFINITY)"] <- "999999"
    input <- sub("^\\$\\(srcdir\\)/", "", sub("^[^ ]+ : ", "", rule))
    run(perl, c(file.path(src, "machine.pl"), args, file.path(src, "data", input), file.path(dat, target)))
  }
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
    if (length(tok)) {
      return(tok[[1L]])
    }
  }
  NA_character_
}

main <- function() {
  if (smoke()) {
    say("already built and working")
    return(invisible())
  }

  cc <- compiler()
  perl <- unname(Sys.which("perl"))
  if (!nzchar(perl)) {
    stop("no perl found on PATH (needed to generate the OligoArrayAux data tables)", call. = FALSE)
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

  dir.create(bin, recursive = TRUE, showWarnings = FALSE)
  dir.create(dat, recursive = TRUE, showWarnings = FALSE)

  say("generating the data tables with machine.pl")
  tables(perl)

  say("compiling with ", cc, " (this takes a minute)")
  with_dir(file.path(src, "src"), {
    build(cc, "hybrid-min")
    build(cc, "hybrid-ss-min")
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
