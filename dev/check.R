# =============================================================================
#  dev/check.R  --  the internal gate: check a change before it lands
# =============================================================================
#  The local gate: CI runs `R CMD check --as-cran` too, and this script adds the
#  checks CI does not render (style, figures, artifact hashes).  Runs the checks
#  a "best package" runs, in order:
#
#    1. NAMESPACE / man/ are in sync with the roxygen comments
#    2. R/ and tests/ are styler-clean and lintr-clean
#    3. every criterion is falsifiable; docs/INDEX.md and the artifact manifest
#       are consistent and fresh
#    4. spelling -- only with `--spelling` (the prose is English now, so spelling
#             is meaningful but needs a WORDLIST before it can gate)
#       dictionary flags every word; the R/roxygen code is checked by lintr)
#    5. figures -- regenerate the figures from the engine and fail if a committed
#       one (man/figures, vignettes/figures/decision-*.svg) changed, so the
#       diagrams never drift from the pipeline (see dev/render-all.R).  By
#       default only the fast, offline decision trees are refreshed; `--render`
#       also redraws the heavy panels (target panel - needs genome
#       data/OligoArrayAux).  Skip the step with `--no-render`.
#    6. install smoke -- R CMD INSTALL must succeed
#    7. R CMD check -- only with `--full` (slow; builds and checks the tarball)
#
#  Modes:
#    --quick  routine loop: steps 1-3 only (roxygen, styler, lintr, criterion
#             fixtures, docs/artifacts); --spelling still adds step 4
#    --no-artifacts  skip step 3 (docs/artifacts).  Use in an isolated worktree:
#             the manifest hashes all of R/, so any edit makes 20 figures stale;
#             the coordinator regenerates it once, after the merge
#    default  the gate before a change lands: steps 1-3 + figures + install smoke
#             (spelling stays opt-in with --spelling)
#    --full   release gate: default + R CMD check --as-cran, which runs the
#             complete test suite (biology/external cross-checks included)
#
#  A missing tool is a message, not a failure, so it runs on a bare machine.
#  Run:  Rscript dev/check.R [--quick] [--full] [--spelling] [--render] [--no-render]
#                            [--no-artifacts]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

full <- "--full" %in% commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
spelling <- "--spelling" %in% commandArgs(trailingOnly = TRUE)
render <- "--render" %in% commandArgs(trailingOnly = TRUE)
no_render <- "--no-render" %in% commandArgs(trailingOnly = TRUE)
no_artifacts <- "--no-artifacts" %in% commandArgs(trailingOnly = TRUE)
if (quick && (full || render)) {
  stop("--quick is the routine loop; use --full/--render for the release gate", call. = FALSE)
}
failed <- character()

message("mode: ", if (quick) {
  "quick (routine loop)"
} else if (full) {
  "full (release gate)"
} else {
  "gate (install smoke)"
})

ok <- function(...) message("  ok   ", ...)

step <- function(name, code) {
  message("\n== ", name, " ==")
  if (!is.null(step_result <- tryCatch({ code(); NULL }, error = function(e) e))) {
    message("  FAIL ", conditionMessage(step_result))
    failed <<- c(failed, name)
  }
}

# 1. roxygen ----------------------------------------------------------------
step("documentation in sync", function() {
  if (!requireNamespace("roxygen2", quietly = TRUE)) return(message("  skip (roxygen2 not installed)"))
  # compare before/after roxygenise (not against git HEAD, which a dirty tree
  # would flag); a new or changed man/ file means the sources drifted
  snapshot <- function() {
    files <- c("NAMESPACE", list.files("man", pattern = "\\.Rd$", full.names = TRUE))
    vapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), character(1))
  }
  before <- snapshot()
  suppressMessages(roxygen2::roxygenise())
  if (!identical(before, snapshot())) stop("stale NAMESPACE/man - run roxygen2::roxygenise() and commit")
  # a quoted export means a non-syntactic name got in - almost always a roxygen
  # tag whose value swallowed prose (the bug that once broke R CMD INSTALL)
  if (any(grepl('export("', readLines("NAMESPACE", warn = FALSE), fixed = TRUE))) {
    stop("NAMESPACE has quoted exports - a roxygen tag likely swallowed text")
  }
  ok("NAMESPACE and man/ match the roxygen comments")
})

# 2. style + lint -----------------------------------------------------------
step("style (styler)", function() {
  if (!requireNamespace("styler", quietly = TRUE)) return(message("  skip (styler not installed)"))
  styler::style_dir("R", dry = "fail")
  styler::style_dir("tests", dry = "fail")
  ok("R/ and tests/ are styler-clean")
})

step("lintr", function() {
  if (!requireNamespace("lintr", quietly = TRUE)) return(message("  skip (lintr not installed)"))
  # object_usage_linter resolves package-internal symbols from a *loaded* build,
  # and this step runs before the install smoke; without loading the source a
  # stale/absent install floods every internal call with "no visible binding"
  if (requireNamespace("pkgload", quietly = TRUE)) {
    suppressMessages(pkgload::load_all(quiet = TRUE))
  }
  l <- lintr::lint_package()
  if (length(l)) {
    print(l)
    stop(length(l), " lint(s)")
  }
  ok("no lints")
})

# 3. criterion fixtures + docs/artifacts graph ------------------------------
step("criterion fixtures", function() {
  # ZW1: a criterion without a fail fixture (or an exemption with a reason) must
  # block the pull request, not just fail a test later.  The ledger lives in
  # tests/testthat/helper-criterion-fixtures.R, which the script sources.
  status <- system2(file.path(R.home("bin"), "Rscript"), "dev/check-criteria.R")
  if (!identical(as.integer(status), 0L)) {
    stop("dev/check-criteria.R failed (status ", status, ")")
  }
  ok("every registered criterion has a fail fixture or a reasoned exemption")
})

step("docs/artifacts", function() {
  if (no_artifacts) {
    return(message("  skip (--no-artifacts; regenerate after the merge)"))
  }
  for (script in c("dev/check-docs.R", "dev/check-artifacts.R")) {
    status <- system2(file.path(R.home("bin"), "Rscript"), script)
    if (!identical(as.integer(status), 0L)) {
      stop(script, " failed (status ", status, ")")
    }
  }
  ok("docs index and artifact manifest are consistent and fresh")
})

# 4. spelling ---------------------------------------------------------------
step("spelling", function() {
  if (quick && !spelling) return(message("  skip (--quick)"))
  if (!spelling) return(message("  skip (pass --spelling)"))
  if (!requireNamespace("spelling", quietly = TRUE)) return(message("  skip (spelling not installed)"))
  words <- spelling::spell_check_package()
  if (nrow(words)) {
    print(words)
    stop(nrow(words), " unknown word(s)")
  }
  ok("no unknown words")
})

# 5. figures ----------------------------------------------------------------
step("figures (auto-render)", function() {
  if (quick) {
    return(message("  skip (--quick)"))
  }
  if (no_render) {
    return(message("  skip (--no-render)"))
  }
  # the layout engine is vendored under .tools/elk (ELK + libavoid), so graphviz
  # `dot` is no longer involved; a missing engine is a skip, not a silent pass
  solver <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote("pkgload::load_all(quiet = TRUE); quit(status = as.integer(!decision_solver_available()))"))
  )
  if (!identical(as.integer(solver), 0L)) {
    return(message("  skip (layout engine not runnable - see .tools/elk/README.md)"))
  }
  # the figures are code: regenerate every one from the engine, then fail only if a
  # *committed* figure changed - so the update is deliberate and committed, never
  # a silent drift.  Both homes are covered: the README pair in man/figures and
  # the per-chemistry decision SVGs the vignettes embed.
  figs <- c(
    list.files("man/figures", pattern = "[.](png|svg)$", full.names = TRUE),
    Sys.glob("vignettes/figures/decision-*.svg")
  )
  md5 <- function() vapply(figs, function(f) unname(tools::md5sum(f)), character(1))
  before <- md5()
  # fast, offline default: the decision trees (README-facing) + their per-scheme
  # siblings in vignettes/figures; `--render` adds the heavy panels and the demo
  # traversal
  args <- if (render) character() else c("--skip-heavy", "--skip-example")
  system2(file.path(R.home("bin"), "Rscript"), c("dev/render-all.R", args))
  after <- md5()
  changed <- figs[is.na(before) | before != after]
  if (length(changed)) {
    stop("generated figures changed - review and commit ", paste(basename(changed), collapse = ", "))
  }
  ok("generated figures are in sync with the engine")
})

# 6. install smoke ----------------------------------------------------------
step("install (smoke)", function() {
  if (quick) return(message("  skip (--quick; the install smoke is in the default/release gate)"))
  lib <- tempfile("kasp-check-lib-")
  dir.create(lib)
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  # --no-byte-compile keeps the smoke fast; `--full` still byte-compiles via
  # R CMD check, so compile problems are not lost
  out <- suppressWarnings(system2("R",
    c("CMD", "INSTALL", "--no-docs", "--no-help", "--no-byte-compile", "-l", shQuote(lib), "."),
    stdout = TRUE, stderr = TRUE
  ))
  if (!identical(attr(out, "status"), NULL) && !identical(as.integer(attr(out, "status")), 0L)) {
    cat(out, sep = "\n")
    stop("R CMD INSTALL failed - the package does not build/load")
  }
  ok("the package installs cleanly")
})

# 7. R CMD check ------------------------------------------------------------
step("R CMD check", function() {
  if (!full) return(message("  skip (pass --full to build and check)"))
  status <- system2("R", c("CMD", "check", "--no-manual", "--as-cran", "."))
  if (!identical(as.integer(status), 0L)) stop("R CMD check failed (status ", status, ")")
  ok("R CMD check clean")
})

# summary -------------------------------------------------------------------
message()
if (length(failed)) {
  message("FAILED: ", paste(failed, collapse = ", "))
  quit(status = 1L)
}
message("all checks passed")
