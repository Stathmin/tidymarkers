# The layered plan grammar: build with +, solve, tubes, plot, export.
# Solved against the bundled, real Arabidopsis example panel (marker_example()).

# The shipped panel needs the TAIR10 genome; use it when it is available
# (cached, so offline), otherwise skip - a synthetic genome would not exercise
# the real targets.
# The TAIR10 genome is ~35 MB compressed: read it once per file, not once per plan.
example_genome <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      value <<- tryCatch(ath_genome(), error = function(e) NULL)
    }
    value
  }
})

example_plan <- function(ids = c("ENSVATH01332998", "ENSVATH02162468"), preset = "kasp") {
  ex <- marker_example()
  ex <- ex[ex$id %in% ids, , drop = FALSE]
  genome <- example_genome()
  if (is.null(genome) || !nrow(ex)) {
    return(NULL)
  }
  marker_plan(genome) + marker_targets(data = ex) + marker_tech(preset = preset)
}

# Solving is deterministic, and every report/plot test needs the same solved
# panel.  Solve it once so the tests assert the same design without repeating
# the search (and keep its cache for the cache-reuse test).
solved_example <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      plan <- example_plan()
      if (is.null(plan)) {
        return(NULL)
      }
      value <<- suppressWarnings(
        marker_solve(plan, flank = 150L, cores = 1L, progress = FALSE)
      )
    }
    value
  }
})

# a tiny stand-in genome for the pure-logic layers (nothing is designed)
dummy_plan <- function() {
  marker_plan(Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGT")))
}

test_that("the real example panel targets carry id labels and their alt2", {
  ex <- marker_example()
  p <- dummy_plan() + marker_targets(data = ex)
  expect_equal(nrow(p$targets), nrow(ex))
  expect_identical(p$targets$label, ex$id)
  expect_true(all(c("label", "chrom", "pos", "ref", "alt", "alt2", "type") %in% names(p$targets)))
  # the multi-allelic target keeps its second ALT
  multi <- ex$id[ex$type == "SNP" & !is.na(ex$alt2)][1L]
  expect_false(is.na(p$targets$alt2[p$targets$label == multi]))
})

test_that("a plan composes in layers with + and prints without solving", {
  p <- dummy_plan() +
    marker_targets(chrom = "Chr1", pos = 4L, ref = "T", alt = "A") +
    marker_context(variants = tibble::tibble(chrom = "Chr1", pos = 6L, alt = "G"))
  expect_s3_class(p, "marker_plan")
  expect_equal(nrow(p$targets), 1L)
  expect_output(print(p), "marker_plan")
  expect_error(p + 1, "layer")
})

test_that("marker_tech orders schemes cheapest-first and honours exclude", {
  p <- dummy_plan() + marker_tech(preset = "kasp")
  s <- plan_schemes(p)
  expect_identical(s[1L], "snp")
  expect_false(is.unsorted(unname(KASP_TECH_COSTS[s])))
  p2 <- dummy_plan() + marker_tech(preset = "lab", exclude = "caps")
  expect_false("caps" %in% plan_schemes(p2))
  expect_true("hrm" %in% plan_schemes(p2))
})

test_that("marker_costs overrides merge over the defaults", {
  c1 <- marker_costs(hrm = 1.2)
  expect_equal(unname(c1[["hrm"]]), 1.2)
  expect_equal(unname(c1[["snp"]]), 1)
  expect_error(marker_costs(1), "named")
})

test_that("marker_solve designs the real example targets and reports", {
  skip_if_no_oligoarrayaux()
  solved <- solved_example()
  skip_if(is.null(solved), "TAIR10 genome not available")
  expect_s3_class(solved, "marker_solved")
  expect_equal(nrow(solved$report), 2L)
  expect_true(all(solved$report$ok))
  expect_true(all(c("label", "scheme", "cost", "ok", "tried", "marker") %in% names(solved$report)))
  expect_output(print(solved), "marker_solved")
  expect_gte(nrow(solved$attempts), 2L)
})

test_that("marker_solve reuses its cache", {
  skip_if_no_oligoarrayaux()
  s1 <- solved_example()
  skip_if(is.null(s1), "TAIR10 genome not available")
  # re-solving with the cache the first solve used must return the same answer
  s2 <- suppressWarnings(
    marker_solve(s1$plan, flank = 150L, cores = 1L, progress = FALSE, cache = s1$cache)
  )
  expect_identical(s1$report$scheme, s2$report$scheme)
  expect_identical(s1$report$marker[[1L]], s2$report$marker[[1L]])
})

test_that("marker_primers pulls the oligos out of a solved marker", {
  skip_if_no_oligoarrayaux()
  solved <- solved_example()
  skip_if(is.null(solved), "TAIR10 genome not available")
  pm <- marker_primers(solved$report$marker[[1L]])
  expect_true(all(c("AS_ref", "AS_alt", "COM") %in% names(pm)))
  expect_true(all(nzchar(pm)))
})

test_that("marker_tubes screens every pair of solved markers", {
  skip_if_no_oligoarrayaux()
  solved <- solved_example()
  skip_if(is.null(solved), "TAIR10 genome not available")
  tb <- marker_tubes(solved)
  expect_true(all(c("a", "b", "n_pairs", "worst_run", "pass") %in% names(tb)))
  expect_equal(nrow(tb), 1L) # one pair of two solved targets
})

test_that("marker_plot and marker_export work off a solved plan", {
  skip_if_no_oligoarrayaux()
  skip_if_not_installed("ggplot2")
  solved <- solved_example()
  skip_if(is.null(solved), "TAIR10 genome not available")
  expect_s3_class(marker_plot(solved), "ggplot")
  dir <- withr::local_tempdir()
  marker_export(solved, dir = dir)
  for (f in c("markers.csv", "attempts.csv", "primers.csv", "solved.rds")) {
    expect_true(fs::file_exists(fs::path(dir, f)), info = f)
  }
})
