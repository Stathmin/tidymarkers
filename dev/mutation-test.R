#!/usr/bin/env Rscript
# =============================================================================
#  dev/mutation-test.R  --  ZY14 mutation testing of three critical subsystems
# =============================================================================
#  `covr` is not installed and must not become a dependency, so the mutators and
#  the runner are part of this script: each mutant is one text-level edit of a
#  single token, applied to a throwaway copy of the package under
#  `fs::path_temp()`.  The real tree is never patched, and nothing is restored by
#  hand - the copy dies with the mutant.
#
#    Rscript dev/mutation-test.R                    # cap 8 mutants per function
#    Rscript dev/mutation-test.R --max-mutants 4    # smaller, faster pass
#    Rscript dev/mutation-test.R --dry-run          # mutant counts only
#
#  Targets (subsystem -> function -> the test file that must catch a mutant):
#    policy      R/marker_arms.R::arms_stab_policy          test-arms-mismatch-position.R
#                R/marker_decision.R::marker_decision_graph test-falsifiable.R
#    invariants  R/marker_outcomes.R::marker_record_outcome test-z8-outcomes.R
#                R/marker_outcomes.R::marker_outcomes       test-z8-outcomes.R
#                R/marker_artefact.R::artefact_id           test-trace-schema.R
#                R/marker_contracts.R::marker_contract      test-contracts.R
#    chemistry   R/marker_lamp.R::design_lamp_marker        test-lamp.R
#                R/marker_mix.R::kasp_mix_recipe            test-mix-recipe.R
#                R/marker_arms.R::design_arms_marker        test-arms.R
#
#  A mutant is killed when any expectation in its target's test file fails or
#  errors.  An all-skip run is a survivor with a note, never a kill.  Exit
#  status is 1 when any subsystem scores below --threshold (default 0.8).
# =============================================================================

targets <- list(
  policy = list(
    label = "policy (classification / routing)",
    functions = list(
      list(
        file = "R/marker_arms.R", name = "arms_stab_policy",
        tests = "tests/testthat/test-arms-mismatch-position.R"
      ),
      list(
        file = "R/marker_decision.R", name = "marker_decision_graph",
        tests = "tests/testthat/test-falsifiable.R"
      )
    )
  ),
  invariants = list(
    label = "invariants (artefact hashing / outcome binding)",
    functions = list(
      list(
        file = "R/marker_outcomes.R", name = "marker_record_outcome",
        tests = "tests/testthat/test-z8-outcomes.R"
      ),
      list(
        file = "R/marker_outcomes.R", name = "marker_outcomes",
        tests = "tests/testthat/test-z8-outcomes.R"
      ),
      list(
        file = "R/marker_artefact.R", name = "artefact_id",
        tests = "tests/testthat/test-trace-schema.R"
      ),
      list(
        file = "R/marker_contracts.R", name = "marker_contract",
        tests = "tests/testthat/test-contracts.R"
      )
    )
  ),
  chemistry = list(
    label = "chemistry (fallback cascade: LAMP / mix / ARMS)",
    functions = list(
      list(
        file = "R/marker_lamp.R", name = "design_lamp_marker",
        tests = "tests/testthat/test-lamp.R"
      ),
      list(
        file = "R/marker_mix.R", name = "kasp_mix_recipe",
        tests = "tests/testthat/test-mix-recipe.R"
      ),
      list(
        file = "R/marker_arms.R", name = "design_arms_marker",
        tests = "tests/testthat/test-arms.R"
      )
    )
  )
)

# Text-level operators, `perl = TRUE` for the lookarounds.  `TRUE`/`FALSE` are
# matched whole-word so `isTRUE` is left alone, and `>`/`<` never eat `|>` or
# `<-`.  A match that does not change the line, or does not parse, is rejected
# and never reaches the score.
mutators <- list(
  list(name = "gte_to_gt", pattern = "(?<![<>=!])>=(?!=)", replacement = ">", perl = TRUE),
  list(name = "lte_to_lt", pattern = "(?<![<>=!])<=(?!=)", replacement = "<", perl = TRUE),
  list(name = "gt_to_gte", pattern = "(?<![-=<>!|])>(?!=)", replacement = ">=", perl = TRUE),
  list(name = "lt_to_lte", pattern = "(?<![-=<>!])<(?![=<-])", replacement = "<=", perl = TRUE),
  list(name = "eq_to_neq", pattern = "(?<![=!<>])==(?!=)", replacement = "!=", perl = TRUE),
  list(name = "neq_to_eq", pattern = "(?<![=!<>])!=(?!=)", replacement = "==", perl = TRUE),
  list(name = "and_to_or", pattern = "&&", replacement = "||", perl = FALSE),
  list(name = "or_to_and", pattern = "\\|\\|", replacement = "&&", perl = TRUE),
  list(name = "true_to_false", pattern = "\\bTRUE\\b", replacement = "FALSE", perl = TRUE),
  list(name = "false_to_true", pattern = "\\bFALSE\\b", replacement = "TRUE", perl = TRUE),
  list(name = "plus_to_minus", pattern = "(?<=[0-9])\\s*\\+\\s*(?=[0-9])", replacement = " - ", perl = TRUE),
  list(name = "times_to_divide", pattern = "(?<=[0-9])\\s*\\*\\s*(?=[0-9])", replacement = " / ", perl = TRUE),
  list(name = "oneL_to_zeroL", pattern = "\\b1L\\b", replacement = "0L", perl = TRUE),
  list(name = "zeroL_to_oneL", pattern = "\\b0L\\b", replacement = "1L", perl = TRUE),
  list(name = "drop_negation", pattern = "(?<![=<>!])!(?!=)(?=\\s*[A-Za-z_(])", replacement = "", perl = TRUE)
)

# --- arguments ---------------------------------------------------------------

parse_args <- function(args) {
  opts <- list(threshold = 0.8, max_mutants = 8L, dry_run = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--dry-run")) {
      opts$dry_run <- TRUE
      i <- i + 1L
      next
    }
    if (arg %in% c("--threshold", "--max-mutants")) {
      if (i == length(args)) {
        cli::cli_abort(gettext("Missing value after {arg}."))
      }
      value <- args[[i + 1L]]
      if (identical(arg, "--threshold")) {
        opts$threshold <- suppressWarnings(as.numeric(value))
      } else {
        opts$max_mutants <- suppressWarnings(as.integer(value))
      }
      i <- i + 2L
      next
    }
    cli::cli_abort(gettextf("Unknown argument: %s", arg))
  }
  if (is.na(opts$threshold) || opts$threshold < 0 || opts$threshold > 1) {
    cli::cli_abort(gettext("`--threshold` must be a number in [0, 1]."))
  }
  if (is.na(opts$max_mutants) || opts$max_mutants < 1L) {
    cli::cli_abort(gettext("`--max-mutants` must be a positive integer."))
  }
  opts
}

# --- mutant generation -------------------------------------------------------

# Line range of a top-level `name <- function(...)` from the parse srcref, so
# the target is the function, not a hand-maintained line number.
function_range <- function(path, name) {
  exprs <- parse(path, keep.source = TRUE)
  srcrefs <- attr(exprs, "srcref")
  for (i in seq_along(exprs)) {
    expr <- exprs[[i]]
    if (!is.call(expr) || !identical(as.character(expr[[1L]]), "<-")) next
    if (!identical(as.character(expr[[2L]]), name)) next
    rhs <- expr[[3L]]
    if (!is.call(rhs) || !identical(as.character(rhs[[1L]]), "function")) next
    return(c(srcrefs[[i]][[1L]], srcrefs[[i]][[3L]]))
  }
  cli::cli_abort(gettextf("Function %s not found in %s.", name, path))
}

replace_at <- function(line, start, length, replacement) {
  paste0(
    substr(line, 1L, start - 1L),
    replacement,
    substr(line, start + length, nchar(line))
  )
}

parses <- function(text) {
  !inherits(tryCatch(parse(text = text), error = identity), "error")
}

# Every candidate mutant of one function: one token replaced per line.  Returns
# the accepted mutants plus the count rejected as a no-op or a syntax error.
mutants_for_function <- function(root, file, name, tests) {
  path <- fs::path(root, file)
  lines <- readLines(path, warn = FALSE)
  range <- function_range(path, name)
  found <- list()
  rejected <- 0L
  for (i in seq(range[[1L]], range[[2L]])) {
    line <- lines[[i]]
    if (grepl("^\\s*#", line)) next
    # Only the code before a trailing comment is mutated, so a comment cannot
    # survive as a fake mutant.
    code <- sub("#.*$", "", line)
    for (mutator in mutators) {
      hits <- gregexpr(mutator$pattern, code, perl = mutator$perl)[[1L]]
      if (hits[[1L]] == -1L) next
      lengths <- attr(hits, "match.length")
      for (k in seq_along(hits)) {
        candidate <- replace_at(line, hits[[k]], lengths[[k]], mutator$replacement)
        if (identical(candidate, line)) {
          rejected <- rejected + 1L
          next
        }
        patched <- lines
        patched[[i]] <- candidate
        text <- paste(patched, collapse = "\n")
        # A mutant that does not parse is a syntax error, not a killed mutant.
        if (!parses(text)) {
          rejected <- rejected + 1L
          next
        }
        found[[length(found) + 1L]] <- list(
          file = file, line = i, original = line, mutated = candidate,
          mutator = mutator$name, tests = tests
        )
      }
    }
  }
  list(mutants = found, rejected = rejected)
}

# Deterministic evenly spaced sample so a capped run covers the whole function
# instead of only its first lines.
cap_mutants <- function(mutants, cap) {
  n <- length(mutants)
  if (n <= cap) {
    return(mutants)
  }
  mutants[unique(round(seq(1, n, length.out = cap)))]
}

# The smallest changed token, so a survivor prints `1L -> 0L` instead of two
# whole source lines.
token_change <- function(original, mutated) {
  a <- strsplit(original, "", fixed = TRUE)[[1L]]
  b <- strsplit(mutated, "", fixed = TRUE)[[1L]]
  n <- min(length(a), length(b))
  prefix <- which(a[seq_len(n)] != b[seq_len(n)])[1L]
  prefix <- if (is.na(prefix)) n else prefix - 1L
  suffix <- which(rev(a)[seq_len(n)] != rev(b)[seq_len(n)])[1L]
  suffix <- if (is.na(suffix)) n else suffix - 1L
  widen <- function(x) {
    lo <- max(prefix, 1L)
    hi <- max(length(x) - suffix, lo)
    while (lo > 1L && grepl("[A-Za-z0-9_.]", x[[lo - 1L]])) lo <- lo - 1L
    while (hi < length(x) && grepl("[A-Za-z0-9_.]", x[[hi + 1L]])) hi <- hi + 1L
    trimws(paste(x[lo:hi], collapse = ""))
  }
  show <- function(x) if (identical(x, "")) "(none)" else x
  paste0(show(widen(a)), " -> ", show(widen(b)))
}

# --- runner ------------------------------------------------------------------

copy_package <- function(root, dir) {
  fs::dir_create(dir)
  fs::dir_copy(fs::path(root, "R"), fs::path(dir, "R"))
  fs::file_copy(fs::path(root, c("DESCRIPTION", "NAMESPACE")), dir)
  for (extra in c("inst", "data")) {
    if (fs::dir_exists(fs::path(root, extra))) {
      fs::dir_copy(fs::path(root, extra), fs::path(dir, extra))
    }
  }
  invisible(dir)
}

# One mutant -> "killed" / "survived" plus the note (skips, error class).
# `mutant$mutated = NULL` runs the unmutated code: the per-target baseline, so a
# mutant "kill" cannot be a broken copy instead of a caught change.
run_mutant <- function(root, mutant, index) {
  dir <- fs::path_temp(paste0("tidymarkers-mut-", Sys.getpid(), "-", index))
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
  copy_package(root, dir)
  if (!is.null(mutant$mutated)) {
    path <- fs::path(dir, mutant$file)
    lines <- readLines(path, warn = FALSE)
    lines[[mutant$line]] <- mutant$mutated
    writeLines(lines, path)
  }

  tryCatch(
    {
      pkgload::load_all(dir, quiet = TRUE, helpers = FALSE)
      results <- lapply(mutant$tests, function(test) {
        testthat::test_file(fs::path(root, test), reporter = "silent")
      })
      classes <- unlist(lapply(results, function(file_result) {
        unlist(lapply(file_result, function(test) {
          vapply(test$results, function(expectation) class(expectation)[[1L]], character(1))
        }))
      }))
      failed <- any(classes %in% c("expectation_failure", "expectation_error"))
      skipped <- sum(classes == "expectation_skip")
      note <- if (failed) {
        ""
      } else if (length(classes) == 0L) {
        "no expectations ran"
      } else if (skipped == length(classes)) {
        "all expectations skipped"
      } else {
        ""
      }
      list(status = if (failed) "killed" else "survived", note = note)
    },
    error = function(e) {
      list(status = "killed", note = paste0("error: ", conditionMessage(e)))
    }
  )
}

# --- main --------------------------------------------------------------------

main <- function() {
  root <- fs::path_real(".")
  if (!fs::file_exists(fs::path(root, "DESCRIPTION"))) {
    cli::cli_abort(gettext("Run this script from the package root."))
  }
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  plan <- lapply(targets, function(target) {
    lapply(target$functions, function(spec) {
      generated <- mutants_for_function(root, spec$file, spec$name, spec$tests)
      list(
        spec = spec,
        total = length(generated$mutants),
        rejected = generated$rejected,
        batch = cap_mutants(generated$mutants, opts$max_mutants)
      )
    })
  })

  if (opts$dry_run) {
    for (name in names(plan)) {
      for (entry in plan[[name]]) {
        cli::cli_inform(gettextf(
          "{name}: {entry$spec$file}::{entry$spec$name} -- {entry$total} mutants ({entry$rejected} rejected), {length(entry$batch)} after the cap"
        ))
      }
    }
    return(invisible(NULL))
  }

  cli::cli_inform(gettextf(
    "ZY14 mutation testing: threshold {opts$threshold}, cap {opts$max_mutants} mutants per function."
  ))

  rows <- list()
  survivors <- list()
  index <- 0L
  for (name in names(plan)) {
    entries <- plan[[name]]
    tests <- unique(unlist(lapply(entries, function(entry) entry$spec$tests)))
    baseline <- run_mutant(root, list(mutated = NULL, tests = tests), 0L)
    if (!identical(baseline$status, "survived")) {
      cli::cli_abort(gettextf(
        "Baseline for {name} does not pass, so a mutant kill would be meaningless: {baseline$note}"
      ))
    }
    killed <- 0L
    total <- 0L
    for (entry in entries) {
      for (mutant in entry$batch) {
        index <- index + 1L
        total <- total + 1L
        started <- Sys.time()
        result <- run_mutant(root, mutant, index)
        elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)
        if (identical(result$status, "killed")) killed <- killed + 1L
        label <- sprintf("%s:%d", mutant$file, mutant$line)
        note <- if (nzchar(result$note)) paste0(" (", result$note, ")") else ""
        cli::cli_inform(gettextf(
          "[{name}] {index}. {label} {mutant$mutator}: {result$status}{note} [{elapsed}s]"
        ))
        if (!identical(result$status, "killed")) {
          survivors[[length(survivors) + 1L]] <- c(
            subsystem = name, label = label,
            change = token_change(mutant$original, mutant$mutated),
            note = result$note
          )
        }
      }
    }
    if (total == 0L) {
      cli::cli_abort(gettextf("No mutants generated for {name}; the target list is stale."))
    }
    rows[[name]] <- c(
      subsystem = name, mutants = total, killed = killed, survived = total - killed,
      score = killed / total
    )
  }

  score_table <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  score_table$mutants <- as.integer(score_table$mutants)
  score_table$killed <- as.integer(score_table$killed)
  score_table$survived <- as.integer(score_table$survived)
  score_table$score <- as.numeric(score_table$score)

  cli::cli_inform(gettext("Per-subsystem mutation score:"))
  print(score_table, row.names = FALSE)

  if (length(survivors)) {
    cli::cli_inform(gettext("Surviving mutants:"))
    for (survivor in survivors) {
      note <- if (nzchar(survivor[["note"]])) paste0(" [", survivor[["note"]], "]") else ""
      cli::cli_inform(gettextf(
        "{survivor[['subsystem']]}: {survivor[['label']]}  {survivor[['change']]}{note}"
      ))
    }
  }

  below <- score_table$score < opts$threshold
  if (any(below)) {
    cli::cli_inform(gettextf(
      "Below threshold {opts$threshold}: {paste(score_table$subsystem[below], collapse = ', ')}."
    ))
    quit(status = 1L)
  }
  cli::cli_inform(gettext("All subsystems meet the threshold."))
  invisible(score_table)
}

if (sys.nframe() == 0L) {
  main()
}
