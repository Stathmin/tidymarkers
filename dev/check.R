# =============================================================================
#  dev/check.R  --  the internal gate: check a change before it lands
# =============================================================================
#  The local gate: CI runs `R CMD check --as-cran` too, and this script adds the
#  checks CI does not render (style, figures, artifact hashes).  Runs the checks
#  a "best package" runs, in order:
#
#    1. NAMESPACE / man/ are in sync with the roxygen comments
#    2. R/ and tests/ are styler-clean and lintr-clean
#   2b. every hand-rolled sequence primitive is registered with a reason
#   3b. each registered replacement becomes a `supersedes` edge in the graph
#    3. every criterion is falsifiable; docs/INDEX.md and the artifact manifest
#       are consistent and fresh
#    4. spelling -- only with `--spelling`: English prose (README, man/, vignettes/)
#       is checked against `inst/WORDLIST`; the Russian NEWS.md/ROADMAP.md are out
#       of scope because `spelling` has no per-file exclusion
#    5. figures -- regenerate the figures from the engine and fail if a committed
#       one (man/figures, vignettes/figures/decision-*.svg) changed, so the
#       diagrams never drift from the pipeline (see dev/render-all.R).  By
#       default only the fast, offline decision trees are refreshed; `--render`
#       also redraws the heavy panels (target panel - needs genome
#       data/OligoArrayAux).  Skip the step with `--no-render`.
#    6. install smoke -- R CMD INSTALL must succeed
#    7. R CMD check -- only with `--full`: dev/test.R --full (parallel, ~7 min)
#       runs the complete suite, then the tarball is built and checked as CI does
#       (`R CMD build` applies .Rbuildignore and prebuilds + indexes the
#       vignettes; a source-directory check does neither and fails on .tmp/)
#
#  Modes:
#    --quick  routine loop: steps 1-3 only (roxygen, styler, lintr, criterion
#             fixtures, docs/artifacts); --spelling still adds step 4
#    --no-artifacts  skip step 3 (docs/artifacts).  Use in an isolated worktree:
#             the manifest hashes all of R/, so any edit makes 20 figures stale;
#             the coordinator regenerates it once, after the merge
#    default  the gate before a change lands: steps 1-3 + figures + install smoke
#             (spelling stays opt-in with --spelling)
#    --full   release gate: default + the whole suite via dev/test.R --full, then
#             R CMD build + R CMD check --as-cran on the tarball.  The suite is
#             not re-run inside R CMD check (that serial rerun cost ~30 min for
#             the same files); CI still runs it inside `R CMD check`.
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

# 2b. hand-roll registry + docs/artifacts graph ------------------------------
step("hand-roll registry", function() {
  # a function that reimplements a sequence primitive Biostrings/TmCalculator
  # already provides must be registered with a reason (see dev/check-handrolled.R)
  status <- system2(file.path(R.home("bin"), "Rscript"), "dev/check-handrolled.R")
  if (!identical(as.integer(status), 0L)) {
    stop("dev/check-handrolled.R failed (status ", status, ")")
  }
  ok("every hand-rolled sequence primitive is registered or replaced")
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

step("hand-roll graph edges", function() {
  # injects `package::function` nodes + `supersedes` edges so the graph answers
  # "what replaced this hand-roll". Mutates the gitignored graphify-out/graph.json
  # (a developer tool, as GWASILIY does) and self-skips when the graph is absent.
  if (!file.exists("graphify-out/graph.json")) {
    return(message("  skip (no graphify-out/graph.json - CI has no graph)"))
  }
  status <- system2(file.path(R.home("bin"), "Rscript"), "dev/graphify-handrolled.R")
  if (!identical(as.integer(status), 0L)) {
    stop("dev/graphify-handrolled.R failed (status ", status, ")")
  }
  ok("hand-roll replacements are in the graph as `supersedes` edges")
})

step("docs/artifacts", function() {
  # the code graph's limits are part of the docs contract; this self-skips when
  # graphify-out/ is absent (CI has no graph, and it is a developer tool)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("dev/graphify-health.R", "--check"))
  if (!identical(as.integer(status), 0L)) {
    stop("dev/graphify-health.R failed (status ", status, ")")
  }
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
  # English prose only: `NEWS.md` and `ROADMAP.md` are Russian and
  # `spell_check_package()` has no per-file exclusion. `inst/WORDLIST` holds the
  # technical terms (chemistries, package names, gene/locus IDs).
  ignore <- readLines("inst/WORDLIST", warn = FALSE)
  files <- c(
    "README.md",
    list.files("man", pattern = "[.]Rd$", full.names = TRUE),
    list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)
  )
  words <- spelling::spell_check_files(files, ignore = ignore, lang = "en_US")
  if (nrow(words)) {
    print(words)
    stop(nrow(words), " unknown word(s)")
  }
  ok("no unknown words in the English prose")
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
  # `targets` is the freshness holder when it is installed: it reruns only the
  # targets the change invalidated (see dev/targets/_targets.R).  Without it,
  # render everything exactly as before.  The md5 compare below is the check in
  # both paths, so a new engine that skips a figure still fails the gate.
  rscript <- file.path(R.home("bin"), "Rscript")
  # the heavy panels are deliberately outside the pipeline, so `--render` keeps
  # the direct runner; the fast default uses targets when it is installed
  if (!render && requireNamespace("targets", quietly = TRUE)) {
    system2(rscript, c("dev/targets.R", "--make"))
  } else {
    args <- if (render) character() else c("--skip-heavy", "--skip-example")
    system2(rscript, c("dev/render-all.R", args))
  }
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
  # Run the suite through the parallel runner once (the SOCK files are sequenced
  # there); R CMD check's own serial rerun costs ~30 min for the same files.
  if (!identical(as.integer(system2("Rscript", c("dev/test.R", "--full"))), 0L)) {
    stop("dev/test.R --full failed")
  }
  ok("the full test suite passes (parallel)")

  work <- tempfile("tm-check-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  old <- setwd(work)
  on.exit(setwd(old), add = TRUE, after = FALSE)
  # Build first, as CI does: a check of the source directory ignores
  # .Rbuildignore (it descends into .tmp/, graphify-out/ and stale *.Rcheck
  # trees and fails on non-portable paths), and the build is what prebuilds the
  # vignettes into inst/doc -- the index R CMD check requires.
  build <- suppressWarnings(system2("R", c("CMD", "build", "--no-manual", root),
    stdout = TRUE, stderr = TRUE
  ))
  tarball <- Sys.glob("tidymarkers_*.tar.gz")
  if (!length(tarball)) {
    cat(build, sep = "\n")
    stop("R CMD build failed")
  }
  # `_R_CHECK_FORCE_SUGGESTS_=false` matches CI's *installed* suggests: the
  # heavy annotation/BLAST/Rsubread packages are optional here (see AGENTS.md).
  status <- system2("R", c(
    "CMD", "check", "--no-manual", "--as-cran", "--no-tests",
    shQuote(tarball[[1L]])
  ), env = "_R_CHECK_FORCE_SUGGESTS_=false")
  if (!identical(as.integer(status), 0L)) {
    log <- file.path(getwd(), "tidymarkers.Rcheck", "00check.log")
    if (file.exists(log)) {
      bad <- grep("ERROR|WARNING", readLines(log, warn = FALSE), value = TRUE)
      cat(bad, sep = "\n")
    }
    stop("R CMD check failed (status ", status, ")")
  }
  ok("R CMD check clean")
})

# summary -------------------------------------------------------------------
message()
if (length(failed)) {
  message("FAILED: ", paste(failed, collapse = ", "))
  quit(status = 1L)
}
message("all checks passed")
