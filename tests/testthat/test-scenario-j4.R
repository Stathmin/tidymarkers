# J4: a specificity scan bounds which markers are worth designing.  The
# committed fixture is small, so the memory index is the declared degradation
# when bowtie is absent; the scan contract is still asserted either way.

test_that("J4: an index scan feeds the next design decision", {
  skip_if_no_oligoarrayaux()
  genome <- read_genome_fasta(scenario_file("j1_genome.fa"))
  index <- scenario_step(
    "j4", "specificity", "index",
    genome_index(genome, strategy = "memory")
  )
  expect_s3_class(index, "genome_index")
  expect_identical(index$strategy, "memory")
  expect_identical(capabilities(index, "max_mismatch"), index$capabilities$max_mismatch)

  primers <- c(
    AS_ref = "CATGGACCTAGCATTGACGTA",
    AS_alt = "CATGGACCTAGCATTGACGTT",
    COM = "GCTTAGCCATGGTCAATGCTA"
  )
  hits <- scenario_step("j4", "specificity", "scan", specificity_of(primers, index))
  expect_s3_class(hits, "specificity_scan")
  expect_true(all(c("primer", "seed", "hits", "engine") %in% names(hits)))
  expect_true(all(hits$hits == 0L))

  # the scan is advisory here: it narrows what is designed, it does not design
  target <- utils::read.csv(scenario_file("j1_target.csv"), stringsAsFactors = FALSE)
  target$pos <- as.integer(target$pos)
  report <- scenario_step(
    "j4", "design", "report",
    design_kasp_markers(genome, target, cores = 1L)
  )
  expect_s3_class(report, "marker_batch")
  expect_identical(nrow(report), 1L)
  expect_identical(report$verdict, "pass")

  # a bowed-out primer must not hit the single 600 bp chromosome anywhere
  clean <- specificity_of(c(unique = "TTTTTTTTTTTTTTTTTTTT"), index)
  expect_equal(clean$hits, 0L)
})
