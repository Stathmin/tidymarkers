#!/usr/bin/env Rscript
# =============================================================================
#  dev/check-licenses.R  --  the Commons licensing gate (ROADMAP ZY18)
# =============================================================================
#  Reads only what DESCRIPTION declares (Imports, Suggests, LinkingTo, Depends),
#  never the installed library: a package we do not declare is not our
#  dependency.  Each installed dependency's `License` field is then decided:
#
#    deny first   AGPL, NonCommercial, NC, Commons Clause, SSPL,
#                 Elastic License, BUSL, UNKNOWN, an empty field      -> fail
#    then allow   MIT, BSD_2_clause, BSD_3_clause, Apache License, GPL, GPL-2,
#                 GPL-3, GPL (>= x), LGPL, Artistic, CC0, CC BY,
#                 file LICENSE                                         -> pass
#    neither      -> fail, the field printed verbatim
#
#  Deny-before-allow is the point, not an implementation detail: `AGPL-3 | MIT`
#  is a failure.  The dual-licence escape hatch is exactly what a
#  Commons-licensing invariant must not accept.  An unrecognised licence is a
#  refusal, not a pass.
#
#  A declared package that is not installed is reported, not failed: Suggests
#  are optional on any given host, and R CMD check's dependency step is what
#  catches a genuinely missing Imports.  The counts in the summary keep a
#  silently shrinking scan visible.
#
#  `R` and the base packages are part of R, not external dependencies, so they
#  are not scanned.
#
#  Run:  Rscript dev/check-licenses.R             (exit 1 on a forbidden licence)
#        Rscript dev/check-licenses.R --self-test  (classifier over literals)
# =============================================================================

root <- local({
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
})
setwd(root)

# The two lists are regexes matched against the `License` field, so a permissive
# combination (`GPL-2 | GPL-3`, `MIT + file LICENSE`, `LGPL (>= 2.1)`) passes on
# its first allowlisted alternative.  `NC` is bounded to a whole token so that it
# catches `CC BY-NC-SA` without matching inside another word.
allowed_licences <- c(
  "MIT",
  "BSD_2_clause",
  "BSD_3_clause",
  "Apache License",
  # `GPL` bare is GPLv2-or-later with the version left unstated (knitr declares
  # it); it is not a weaker grant than `GPL-2`, so it passes.  `AGPL` cannot
  # reach here -- the deny list is checked first.
  "^GPL( |$|\\()",
  "GPL-2",
  "GPL-3",
  "GPL \\(>=\\s*[0-9]",
  "LGPL",
  "Artistic",
  "CC0",
  "CC BY",
  "file LICENSE"
)

forbidden_licences <- c(
  "AGPL",
  "NonCommercial",
  "(^|[^A-Za-z])NC([^A-Za-z]|$)",
  "Commons Clause",
  "SSPL",
  "Elastic License",
  "BUSL",
  "UNKNOWN"
)

# The whole rule, in one place, so the scan and --self-test cannot drift.
# Returns list(ok = <logical>, reason = <character>).
licence_verdict <- function(x) {
  if (length(x) != 1L || is.na(x)) {
    x <- ""
  }
  x <- trimws(x)
  if (!nzchar(x)) {
    return(list(ok = FALSE, reason = "empty license field"))
  }
  if (any(vapply(forbidden_licences, grepl, logical(1), x = x))) {
    return(list(ok = FALSE, reason = "forbidden licence"))
  }
  if (any(vapply(allowed_licences, grepl, logical(1), x = x))) {
    return(list(ok = TRUE, reason = ""))
  }
  list(ok = FALSE, reason = "unrecognised licence")
}

# Named character vector: names are packages, values are the declaring field.
declared_dependencies <- function(file = "DESCRIPTION") {
  wanted <- c("Depends", "Imports", "LinkingTo", "Suggests")
  dcf <- read.dcf(file)
  deps <- character()
  for (field in intersect(wanted, colnames(dcf))) {
    items <- trimws(strsplit(dcf[1L, field], ",", fixed = TRUE)[[1L]])
    items <- trimws(sub("\\s*\\(.*\\)\\s*$", "", items))
    items <- items[nzchar(items)]
    deps <- c(deps, stats::setNames(rep(field, length(items)), items))
  }
  deps <- deps[!duplicated(names(deps))]
  # `R` and the base packages ship with R and are not external dependencies.
  deps[!(names(deps) %in% c("R", rownames(installed.packages(priority = "base"))))]
}

is_installed <- function(pkg) length(find.package(pkg, quiet = TRUE)) > 0L

licence_field <- function(pkg) {
  value <- suppressWarnings(packageDescription(pkg, fields = "License"))
  if (length(value) != 1L) NA_character_ else value
}

scan_licenses <- function(deps) {
  pkgs <- names(deps)
  installed <- vapply(pkgs, is_installed, logical(1))
  problems <- character()
  for (pkg in pkgs[installed]) {
    field <- licence_field(pkg)
    verdict <- licence_verdict(field)
    if (!verdict$ok) {
      shown <- if (is.na(field)) "<missing>" else field
      problems <- c(problems, paste0(shown, "  ", pkg, " (", deps[[pkg]], ")  -- ", verdict$reason))
    }
  }
  list(
    declared = length(pkgs),
    installed = sum(installed),
    not_installed = pkgs[!installed],
    problems = problems
  )
}

self_test <- function() {
  inputs <- c(
    "MIT", "GPL-3", "AGPL-3", "AGPL-3 | MIT", "MIT + file LICENSE",
    "NonCommercial", "Commons Clause", "SSPL", "", "UNKNOWN",
    "LGPL (>= 2.1)", "Apache License (>= 2)", "GPL",
    "GPL-2 | file LICENSE", "CC BY-NC-SA 4.0"
  )
  expected <- c(
    TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
    TRUE, FALSE
  )
  got <- vapply(inputs, function(x) licence_verdict(x)$ok, logical(1))
  for (i in seq_along(inputs)) {
    cat("  ", if (got[[i]]) "pass" else "fail", " ", sQuote(inputs[[i]]), "\n", sep = "")
  }
  if (any(got != expected)) {
    cat("self-test FAILED: ", paste0(sQuote(inputs[got != expected]), collapse = ", "), "\n", sep = "")
    quit(status = 1L)
  }
  cat("self-test ok\n")
  quit(status = 0L)
}

if ("--self-test" %in% commandArgs(trailingOnly = TRUE)) {
  self_test()
}

scan <- scan_licenses(declared_dependencies())
not_installed <- if (length(scan$not_installed)) {
  paste0(" (", length(scan$not_installed), " not installed: ", paste(scan$not_installed, collapse = ", "), ")")
} else {
  ""
}
counts <- paste0(
  scan$declared, " declared, ", scan$installed, " installed", not_installed,
  ", ", scan$installed, " checked"
)

if (length(scan$problems)) {
  message("check-licenses: ", counts)
  for (problem in scan$problems) {
    message(problem)
  }
  message("fail: ", length(scan$problems), " forbidden licence(s)")
  quit(status = 1L)
}
message("check-licenses: ", scan$declared, " declared, ", scan$installed, " installed", not_installed, ", 0 forbidden")
