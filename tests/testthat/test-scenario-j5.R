# J5: the lab answers back and the next design uses it.  The loop is
# artefact -> outcomes -> priors -> tech -> solve; the prior may only reorder
# the attempts, never change a verdict.

test_that("J5: a recorded outcome reorders the next solve and keeps its verdict", {
  skip_if_no_oligoarrayaux()
  scenario_timer("J5")
  genome <- read_genome_fasta(scenario_file("j1_genome.fa"))
  target <- utils::read.csv(scenario_file("j1_target.csv"), stringsAsFactors = FALSE)
  target$pos <- as.integer(target$pos)

  plan <- marker_plan(genome) +
    marker_targets(target) +
    marker_tech(schemes = "snp")
  solved <- scenario_step("j5", "design", "solved", marker_solve(plan, progress = FALSE))
  art <- scenario_step(
    "j5", "artefact", "written",
    marker_artefact(solved, dir = withr::local_tempdir(), scenario = "J5")
  )
  label <- solved$report$label[[1L]]
  marker_record_outcome(art, "works", target = label)

  read_back <- scenario_step(
    "j5", "artefact", "read",
    marker_artefact_read(art$dir)
  )
  expect_s3_class(read_back, "marker_artefact")
  expect_identical(read_back$manifest$contract$scenario, "J5")

  outcomes <- scenario_step("j5", "artefact", "outcomes", marker_outcomes(read_back))
  expect_s3_class(outcomes, "data.frame")
  expect_equal(nrow(outcomes), 1L)
  expect_identical(outcomes$result, "works")
  expect_identical(outcomes$scheme, "snp")

  priors <- scenario_step("j5", "design", "priors", marker_priors(outcomes))
  expect_s3_class(priors, "priors")
  expect_equal(priors$works, 1L)
  expect_equal(priors$fails, 0L)
  expect_false(any(priors$calibrated))

  # hrm is cheaper than snp under these costs, so an uninformed solve tries it
  # first; the laboratory success of snp moves it ahead, without dropping hrm.
  costs <- marker_costs(snp = 3, hrm = 2)
  base <- marker_plan(genome) + marker_targets(target) +
    marker_tech(schemes = c("hrm", "snp"), costs = costs)
  informed <- marker_plan(genome) + marker_targets(target) +
    marker_tech(schemes = c("hrm", "snp"), costs = costs, priors = priors)
  expect_identical(tidymarkers:::plan_schemes(base), c("hrm", "snp"))
  expect_identical(tidymarkers:::plan_schemes(informed), c("snp", "hrm"))

  base_solved <- marker_solve(base, progress = FALSE)
  next_solved <- scenario_step(
    "j5", "design", "next_solved",
    marker_solve(informed, progress = FALSE)
  )
  expect_s3_class(next_solved, "marker_solved")
  expect_identical(next_solved$report$ok, base_solved$report$ok)
  # the order decides which workable scheme wins; the verdict does not move
  expect_identical(next_solved$attempts$scheme[[1L]], "snp")
  expect_identical(next_solved$plan$priors, priors)
  expect_identical(names(next_solved$attempts)[1L], "label")
})
