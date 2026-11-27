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

test_that("J4: three builds, six shards, and a missing shard names both", {
  skip_if_not_installed("Rbowtie")
  # committed fixture: the site is on A only, so a pooled scan would invent it
  site <- "AAGTCTGGGCAGAACCTTCACCAA"
  fa <- c(
    A = scenario_file("j4_a.fa"),
    B = scenario_file("j4_b.fa"),
    D = scenario_file("j4_d.fa")
  )
  dir <- withr::local_tempdir()

  pan <- lapply(names(fa), function(build) {
    genome_index(
      marker_genome(fa[[build]], build = build),
      "bowtie",
      chroms = c("chr1", "chr2"),
      dir = file.path(dir, build)
    )
  })
  names(pan) <- names(fa)

  expect_identical(
    sum(vapply(pan, function(x) nrow(x$manifest), integer(1))),
    6L
  )
  for (build in names(pan)) {
    expect_true(all(pan[[build]]$manifest$genome_build == build))
    expect_identical(pan[[build]]$build, build)
  }

  hits <- specificity_of(c(AS = site), pan$A, k = 15L, max_mismatch = 0L)
  expect_identical(unique(hits$build), "A")
  expect_identical(hits$hits, 1L)
  for (build in c("B", "D")) {
    other <- specificity_of(c(AS = site), pan[[build]], k = 15L, max_mismatch = 0L)
    expect_identical(unique(other$build), build)
    expect_identical(other$hits, 0L)
  }

  # a missing B shard names the build and the chromosome
  b_chr2 <- pan$B$manifest$path[pan$B$manifest$chrom == "chr2"]
  unlink(Sys.glob(paste0(b_chr2, "*")))
  err <- tryCatch(index_verify(pan$B), error = identity)
  expect_s3_class(err, "tidymarkers_index_missing_error")
  expect_match(conditionMessage(err), "build B")
  expect_match(conditionMessage(err), "chr2")
})
