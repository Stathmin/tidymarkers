# =============================================================================
#  dev/single-score-baseline.R  --  can one scalar rank design quality?
# =============================================================================
#  The perspectives draft (§5) claims that any single scalar "design quality"
#  score is falsifiable: it does not separate accepted from rejected
#  constructions.  This harness makes that claim runnable instead of prose.
#
#  Method: build one plausible scalar — a weighted sum of checks the package
#  already has, evaluated from sequence alone — and score a panel whose outcome
#  was established by an operator before the package saw it.  The number
#  reported is how many constructions the scalar orders WRONGLY: every pair
#  (rejected, accepted) where the rejected construction gets a better score.
#
#  Panel: inst/extdata/sharplab2021_assays.csv — 50 KASP assays from
#  Sharplab 2021 (10.1186/s12870-021-02932-8), 46 polymorphic between 1212 and
#  B73 and 4 not polymorphic.  Both classes carry full primer sequences, so the
#  scalar is computable for every row.
#
#  The weights below are a GUESS — one plausible weighting, not "the" scalar.
#  A different weighting is a different falsification attempt; the claim only
#  needs one credible scalar to fail.  To show the failure is not a property of
#  one hand-picked weighting, four alternative weightings are also run.  They
#  were fixed before the run and every outcome is reported: this is not a fit,
#  and no weighting is selected afterwards.
#
#  The FAM/HEX tails are not a fixed additive offset: they enter the normalised
#  Tm and GC features of the whole oligo, so stripping them reorders the panel
#  (147 wrong pairs instead of 71).  That row is kept as a caveat, not a bug:
#  both weightings fail, which is the claim, and the size of the swing says a
#  scalar is sensitive to a decision the sequence alone cannot settle.
#
#  Scalar (cost, lower is better), summed over the three primers:
#    tm_dev    |NN Tm - 60|            x 1.0   per degree C
#    gc_dev    |GC - 0.5|              x 20    per fraction
#    run       longest homopolymer - 4 x 2     per base over the KASP bar
#    hairpin   max(0, -dG37)           x 1.0   per kcal/mol of hairpin
#    clamp     3' base not G/C         x 1.0
#  All five are existing checks: `primer_tm()`, `gc_frac()`, `longest_run()`
#  and `internal_hairpin_dg()`; no new criterion is introduced.
#
#  Blocker, stated in the output and not hidden: a corpus-wide run needs
#  outcomes from `marker_record_outcome()`, and the tree holds none.  Only
#  published panels have a known outcome, so the corpus row is `not_run` with
#  that reason.
#
#  Rejected panels (no primer sequences, or no rejected class):
#    soybean_worked_panel_kasp_gates.csv  all 13 assays worked
#    z8_gate_outcomes.csv                 all 36 constructions worked
#    soybean_worked_panel_assays.csv      all 36 assays worked
#    soybean_est_caps_markers.csv         only the 26 successes recorded
#    polyploid_cotton2023_kasp_screen.csv GOOD/BAD classes but no sequences
#    meng2026_soybean_per_assay.csv       4 not_attempted, 2 worked
#
#  Run:  Rscript dev/single-score-baseline.R [--check]
#        --check runs the panel twice and errors if the row count or the
#        verdict distribution moves between runs.
#
#  Writes .tmp/single-score-baseline.csv (git-ignored):
#    task | family | expected | observed | verdict | where | note
#
#  Rows: one per construction, two summary rows, one per sweep variant, and one
#  `not_run` corpus row.
# =============================================================================

# ---- plumbing ----------------------------------------------------------------

baseline_root <- function() {
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))))
  }
  wd <- getwd()
  repeat {
    if (file.exists(file.path(wd, "DESCRIPTION")) && dir.exists(file.path(wd, "R"))) {
      return(wd)
    }
    parent <- dirname(wd)
    if (identical(parent, wd)) {
      return(getwd())
    }
    wd <- parent
  }
}

baseline_extdata <- function(name) {
  shipped <- system.file("extdata", name, package = "tidymarkers")
  if (nzchar(shipped)) {
    return(shipped)
  }
  file.path(baseline_root(), "inst", "extdata", name)
}

baseline_path <- function(...) file.path(baseline_root(), ...)

baseline_rows <- function(task, family, expected, observed, verdict, where, note) {
  n <- length(task)
  recycle <- function(x) {
    if (length(x) == n) x else rep(x, length.out = n)
  }
  tibble::tibble(
    task = as.character(task),
    family = recycle(as.character(family)),
    expected = recycle(as.character(expected)),
    observed = recycle(as.character(observed)),
    verdict = recycle(as.character(verdict)),
    where = recycle(as.character(where)),
    note = recycle(as.character(note))
  )
}

baseline_not_run <- function(task, family, expected, reason) {
  baseline_rows(task, family, expected, NA_character_, "not_run", NA_character_, reason)
}

# ---- the scalar --------------------------------------------------------------

# One plausible weighting, named so the reader can re-run with another.
BASELINE_WEIGHTS <- c(tm_dev = 1.0, gc_dev = 20, run = 2, hairpin = 1.0, clamp = 1.0)

BASELINE_TARGET_TM <- 60
BASELINE_MAX_RUN <- 4

# Four other weightings, written down before the run.  Reported as-is: a second
# weighting that also fails is the point, and one that separated the panel
# would have been reported too.
BASELINE_ALT_WEIGHTS <- list(
  tm_only = c(tm_dev = 1.0, gc_dev = 0, run = 0, hairpin = 0, clamp = 0),
  tm_dominant = c(tm_dev = 3.0, gc_dev = 10, run = 1.0, hairpin = 0.5, clamp = 0.5),
  gc_dominant = c(tm_dev = 0.5, gc_dev = 60, run = 2, hairpin = 1.0, clamp = 1.0),
  gentle = c(tm_dev = 0.5, gc_dev = 10, run = 1.0, hairpin = 0.5, clamp = 0)
)

# The two 5' tails KASP puts on its allele primers.  Stripping them is a
# sensitivity variant, not a neutral control — see the header.
BASELINE_TAILS <- c("GAAGGTGACCAAGTTCATGCT", "GAAGGTCGGAGTCAACGGATT")

baseline_strip_tail <- function(seq) {
  for (tail in BASELINE_TAILS) {
    seq <- sub(paste0("^", tail), "", seq)
  }
  seq
}

# Cost contribution of one oligo.  `NA` propagates: an unmeasurable feature
# makes the construction unscored instead of silently free.
baseline_primer_cost <- function(seq, w = BASELINE_WEIGHTS) {
  seq <- as.character(seq)
  ok <- grepl("^[ACGT]+$", seq)
  tm <- primer_tm(seq)
  gc <- gc_frac(seq)
  run <- vapply(seq, longest_run, numeric(1))
  dg <- internal_hairpin_dg(seq)
  clamp <- !grepl("[GC]$", seq)
  cost <- w[["tm_dev"]] * abs(tm - BASELINE_TARGET_TM) +
    w[["gc_dev"]] * abs(gc - 0.5) +
    w[["run"]] * pmax(0, run - BASELINE_MAX_RUN) +
    w[["hairpin"]] * pmax(0, -dg) +
    w[["clamp"]] * clamp
  ifelse(ok, cost, NA_real_)
}

# Cost of a KASP assay: the three oligos carry the FAM/HEX tails, which is what
# a sequence-only score can see.  The sum, not the mean, so a longer tail set
# does not read as a better design.
baseline_assay_cost <- function(fa, fd, r, w = BASELINE_WEIGHTS) {
  n <- length(fa)
  out <- numeric(n)
  for (i in seq_len(n)) {
    parts <- c(fa[[i]], fd[[i]], r[[i]])
    out[[i]] <- sum(baseline_primer_cost(parts, w))
  }
  out
}

# ---- the panel ---------------------------------------------------------------

baseline_panel <- function() {
  path <- baseline_extdata("sharplab2021_assays.csv")
  if (!file.exists(path)) {
    return(NULL)
  }
  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("no", "marker", "primer_fa", "primer_fd", "primer_r", "assay_worked")
  if (!all(need %in% names(raw)) || nrow(raw) == 0L) {
    return(NULL)
  }
  rejected <- !as.logical(raw$assay_worked)
  cost <- baseline_assay_cost(raw$primer_fa, raw$primer_fd, raw$primer_r)
  tibble::tibble(
    no = raw$no,
    marker = raw$marker,
    rejected = rejected,
    fa = raw$primer_fa,
    fd = raw$primer_fd,
    r = raw$primer_r,
    cost = cost
  )
}

# ---- the verdict -------------------------------------------------------------

baseline_falsify <- function(panel) {
  scored <- panel[is.finite(panel$cost), , drop = FALSE]
  n_rej <- sum(scored$rejected)
  n_acc <- nrow(scored) - n_rej

  # A pair (rejected, accepted) is ordered wrongly when the rejected
  # construction gets the strictly better (lower) cost.
  inversions <- if (n_rej == 0L || n_acc == 0L) {
    NA_integer_
  } else {
    sum(vapply(
      which(scored$rejected),
      function(i) sum(scored$cost[!scored$rejected] > scored$cost[[i]]),
      integer(1)
    ))
  }
  pairs <- n_rej * n_acc

  # Accepted constructions that fall below at least one rejected one.
  worst_rej <- if (n_rej) min(scored$cost[scored$rejected]) else NA_real_
  offenders <- if (n_rej) sum(!scored$rejected & scored$cost > worst_rej) else NA_integer_

  # Where does the worst accepted construction sit, and how many rejected
  # constructions are inside the best n_rej by cost (a clean separation would
  # put all of them there)?
  ord <- order(scored$cost)
  rejected_rank <- rank(scored$cost, ties.method = "average")

  list(
    n = nrow(scored),
    n_used = nrow(panel),
    n_rej = n_rej,
    n_acc = n_acc,
    inversions = inversions,
    pairs = pairs,
    offenders = offenders,
    rejected_rank = rejected_rank[scored$rejected],
    top_rej = sum(scored$rejected[ord][seq_len(n_rej)]),
    scored = scored
  )
}

# One row per alternative weighting plus the tail-stripped variant.  Same panel,
# same metric; only the scalar or the sequences change.
baseline_variant_rows <- function(panel) {
  variants <- c(
    list(tails_stripped = baseline_assay_cost(
      baseline_strip_tail(panel$fa), baseline_strip_tail(panel$fd), baseline_strip_tail(panel$r)
    )),
    lapply(BASELINE_ALT_WEIGHTS, function(w) baseline_assay_cost(panel$fa, panel$fd, panel$r, w))
  )
  rows <- lapply(names(variants), function(nm) {
    sub <- panel
    sub$cost <- variants[[nm]]
    res <- baseline_falsify(sub)
    observed <- if (is.na(res$inversions)) {
      sprintf("not comparable: %d rejected, %d accepted", res$n_rej, res$n_acc)
    } else {
      sprintf(
        "%d of %d pairs; %d of %d accepted below a rejected one",
        res$inversions, res$pairs, res$offenders, res$n_acc
      )
    }
    baseline_rows(
      sprintf("baseline:sensitivity %s", nm),
      "panel",
      "0 of the pairs",
      observed,
      if (!is.na(res$inversions) && res$inversions == 0L) "pass" else "fail",
      if (identical(nm, "tails_stripped")) "BASELINE_TAILS" else "BASELINE_ALT_WEIGHTS",
      if (identical(nm, "tails_stripped")) {
        "tails feed the normalised Tm and GC of the whole oligo; stripping them is another scalar, not a control"
      } else {
        "a different weighting is a different falsification attempt, not a refit"
      }
    )
  })
  do.call(rbind, rows)
}

baseline_tasks <- function() {
  panel <- baseline_panel()
  if (is.null(panel)) {
    return(baseline_not_run(
      "baseline:sharplab2021 panel", "panel", "a scored construction",
      "inst/extdata/sharplab2021_assays.csv absent or has no primer columns; no construction to score"
    ))
  }
  res <- baseline_falsify(panel)
  if (res$n == 0L) {
    return(baseline_not_run(
      "baseline:sharplab2021 panel", "panel", "a scored construction",
      "no construction in the panel produced a finite scalar; a missing value is not a zero"
    ))
  }
  if (is.na(res$inversions)) {
    return(baseline_not_run(
      "baseline:sharplab2021 panel", "panel", "two outcome classes",
      sprintf(
        "panel has %d rejected and %d accepted constructions; an ordering needs both", res$n_rej, res$n_acc
      )
    ))
  }

  # Per-construction rows: a construction is `pass` when it takes part in no
  # inversion, `fail` when at least one pair is ordered the wrong way.
  worst_rej <- min(res$scored$cost[res$scored$rejected])
  per <- res$scored
  wrong <- ifelse(
    per$rejected,
    per$cost < max(per$cost[!per$rejected]),
    per$cost > worst_rej
  )
  rank_all <- rank(per$cost, ties.method = "average")
  rows <- baseline_rows(
    sprintf("sharplab2021:%s %s", per$no, per$marker),
    "panel",
    ifelse(per$rejected, "rejected", "accepted"),
    sprintf("cost %.2f rank %s", per$cost, format(rank_all, trim = TRUE)),
    ifelse(wrong, "fail", "pass"),
    "primer_fa/primer_fd/primer_r",
    ifelse(
      per$rejected,
      "the 4 assays the paper reports as not polymorphic between 1212 and B73",
      "the 46 assays the paper reports as polymorphic"
    )
  )

  summary <- baseline_rows(
    c("baseline:pairwise inversions", "baseline:rejected inside the best k"),
    "panel",
    c("0 of the pairs", sprintf("%d of %d rejected", res$n_rej, res$n_rej)),
    c(
      sprintf(
        "%d of %d pairs; %d of %d accepted below a rejected one",
        res$inversions, res$pairs, res$offenders, res$n_acc
      ),
      sprintf(
        "%d rejected in the best %d by cost (ranks %s)",
        res$top_rej, res$n_rej,
        paste(format(sort(res$rejected_rank), trim = TRUE), collapse = ", ")
      )
    ),
    c(
      if (res$inversions == 0L) "pass" else "fail",
      if (res$top_rej == res$n_rej) "pass" else "fail"
    ),
    "dev/single-score-baseline.R",
    c(
      "a rejected construction scoring above an accepted one is the falsification",
      "clean separation would put all rejected constructions in the best k"
    )
  )

  corpus <- baseline_not_run(
    "baseline:corpus", "corpus", "recorded outcomes from marker_record_outcome()",
    "the tree holds no marker_record_outcome() rows, so a corpus-wide run would be invented; published panels are the only known-outcome set"
  )

  sweep <- baseline_variant_rows(res$scored)

  rbind(rows, summary, sweep, corpus)[, c("task", "family", "expected", "observed", "verdict", "where", "note")]
}

# ---- driver ------------------------------------------------------------------

baseline_write <- function(x, path = NULL) {
  if (is.null(path)) {
    path <- baseline_path(".tmp", "single-score-baseline.csv")
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

baseline_distribution <- function(x) {
  tab <- table(factor(x$verdict, levels = c("pass", "fail", "not_run", "gap")))
  as.integer(tab)
}

main <- function() {
  check <- "--check" %in% commandArgs(trailingOnly = TRUE)
  if (!"tidymarkers" %in% loadedNamespaces()) {
    suppressMessages(pkgload::load_all(baseline_root(), quiet = TRUE))
  }
  x <- baseline_tasks()
  path <- baseline_write(x)
  n <- nrow(x)
  dist <- baseline_distribution(x)
  names(dist) <- c("pass", "fail", "not_run", "gap")
  message(sprintf(
    "single-score baseline: %d row(s); %s", n,
    paste(sprintf("%s %d", names(dist), dist), collapse = ", ")
  ))
  message("wrote ", path)
  summary <- x[grepl("^baseline:", x$task), , drop = FALSE]
  for (i in seq_len(nrow(summary))) {
    message(sprintf("  %s: %s", summary$task[[i]], summary$observed[[i]]))
  }
  if (check) {
    y <- baseline_tasks()
    if (nrow(y) != n || !identical(baseline_distribution(y), unname(dist))) {
      stop("--check: the row count or verdict distribution moved between runs", call. = FALSE)
    }
    message("--check: reproduced")
  }
  invisible(x)
}

if (sys.nframe() == 0L) {
  main()
}
