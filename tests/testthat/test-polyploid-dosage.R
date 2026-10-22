# =============================================================================
#  test-polyploid-dosage
# =============================================================================
#  Roadmap Z10.2: dosage / copy-number as its own axis for polyploid KASP, so
#  that it is NOT folded into `pop_freq` (Z6.1).  Nothing under `R/` implements
#  dosage yet, and the roadmap's own blocker is that `fitPoly` is a Suggests
#  that may only be introduced with a real implementation (ZA4).  So this file
#  pins what the LANDED DATA can support, and refuses what it cannot:
#
#    * no source reports allele dosages in a number form.  Jaimes states the
#      eleven fixed proportions (0.0, 0.1, ..., 1.0) and counts how many of 76
#      markers got a call (42); it does not tabulate a dose per marker.  Cotton
#      and durum never use the word; gao2022 has no dosage column;
#
#    * the four sources run FOUR different mixture/call problems - decaploid
#      (11 classes), diploid (2), allopolyploid-with-homoeologs (unspecified),
#      and none at all (cotton, cluster separation only).  A single pooled
#      polyploid rate is therefore not constructible, and these tests forbid
#      one instead of inventing a tolerance;
#
#    * "a" / "ab" / "aab" mean what the PLODY says, not what the species says:
#      on sugarcane "a" is 1 of ~10, on diploid cotton "a" is 1 of 2, and the
#      durum array never separates the A/B homoeologs enough to say which of the
#      two copies an "a" is;
#
#    * cotton is a NEGATIVE datapoint on this axis: it does not report dosage
#      at all, so its 8-of-11 screen must not be averaged into any dosage rate;
#
#    * `marker_population_checks()` is a two-allele frequency table today.  It
#      has no dosage axis, and passing it a polyploid frequency silently prints
#      "50%" for a state the data cannot distinguish from a diploid heterozygote.
#
#  Hard-evidence note: no new `inst/extdata` file is landed by this pass.  The
#  per-marker dose vector Z10.2 would want is not in the fetched files - see the
#  `blocked` row this file cites for `data-raw/science_predictions.R`.
# =============================================================================

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

# the dosage / copy-number vocabulary of a free-text column set
DOSAGE_WORDS <- "dosage|copy.?number|dose\\b|simplex|duplex|triplex"

## --- 1. which dosage classes each source actually reports -------------------

test_that("dosage: Jaimes states eleven fixed decaploid proportions, not doses", {
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  # the fixed category set is the paper's own words; it lives in the column the
  # data-raw script writes, so that the 11 classes are traceable to a row
  pd <- unique(j$ploidy_detail)
  expect_length(pd, 1L)
  expect_true(grepl("11 dosage categories", pd, fixed = TRUE))
  expect_true(grepl("0.0-1.0", pd, fixed = TRUE))
  expect_true(grepl("decaploid", pd, fixed = TRUE))
  # 0.0..1.0 in steps of 0.1 is ELEVEN values, and the paper lists all eleven
  expect_equal(length(seq(0, 1, by = 0.1)), 11L)
  # the landed count is a call count (42 of 76), never a per-marker dose:
  # the file has no column that could carry a dose value
  expect_false(any(grepl(DOSAGE_WORDS, tolower(setdiff(
    names(j), c("outcome_definition")
  )))))
  expect_equal(j$n_worked[j$metric == "dosage assigned (n)"], 42)
  expect_equal(unique(j$n_tested[j$key == "jaimes2024"]), 76)
})

test_that("dosage: the sugarcane split is 42 assigned / 34 unassigned, both named", {
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  conv <- j[j$key == "jaimes2024", ]
  # n_worked repeats the 42 on every conversion row (worked + failed == tested)
  expect_equal(unique(conv$n_worked), 42)
  # the "others" side is a row of its own; its count is carried in `value`, the
  # metric's number, next to the 42-count row that uses the same column
  unassigned <- conv[conv$metric == "dosage not assigned (n)", ]
  expect_equal(nrow(unassigned), 1L)
  expect_equal(unassigned$value, 34)
  expect_equal(unassigned$n_failed, 34)
  # the 34 are a CLUSTER-FIT failure and the paper says so; a reader must not
  # be able to reuse this as an amplification failure
  expect_true(all(grepl("did not cluster", conv$outcome_definition, fixed = TRUE)))
  expect_true(all(grepl("all 76 SNPs amplified", conv$outcome_definition, fixed = TRUE)))
})

test_that("dosage: only one of the four sources is a dosage source", {
  files <- c(
    jaimes2024 = "polyploid_jaimes2024_sugarcane_dosage.csv",
    gao2022 = "polyploid_gao2022_sugarcane_chain.csv",
    durum2019 = "polyploid_durum2019_hexaploid_vs_tetraploid.csv",
    cotton2023 = "polyploid_cotton2023_kasp_screen.csv"
  )
  says_dosage <- vapply(names(files), function(nm) {
    d <- ext(files[[nm]])
    txt <- paste(c(
      d$outcome_definition,
      if ("ploidy_detail" %in% names(d)) d$ploidy_detail
    ), collapse = " ")
    grepl("dosage|dose", txt, ignore.case = TRUE)
  }, logical(1))
  expect_equal(unname(says_dosage), c(TRUE, FALSE, FALSE, FALSE))
})

## --- 2. no single polyploid conversion/typing rate is defensible ------------

test_that("dosage: the four sources are not one rate (distinct denominators)", {
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  g <- ext("polyploid_gao2022_sugarcane_chain.csv")
  d <- ext("polyploid_durum2019_hexaploid_vs_tetraploid.csv")
  c <- ext("polyploid_cotton2023_kasp_screen.csv")
  # SNP markers, nested SNP tiers, array assays, KASP assays: four units
  expect_equal(sort(unique(j$n_tested)), c(42, 76))
  expect_equal(g$n, c(77, 49, 24))
  expect_equal(d$n_assays, c(73647, 7940, 81587))
  expect_equal(unique(c$n_tested), 11)
  # and their success fractions are far apart, so a pooled rate would be an
  # artefact of whichever unit it was computed in
  rates <- c(
    jaimes_dosage = 42 / 76,
    gao_two_parent = 49 / 77,
    gao_on_map = 24 / 77,
    durum_tetraploid_nonfailed = (7940 - 703) / 7940,
    cotton_screen = 8 / 11
  )
  expect_gt(max(rates) - min(rates), 0.3)
  # no file may carry a cross-source pooled rate column
  for (nm in c("jaimes2024", "gao2022", "durum2019", "cotton2023")) {
    d2 <- ext(list(
      jaimes2024 = "polyploid_jaimes2024_sugarcane_dosage.csv",
      gao2022 = "polyploid_gao2022_sugarcane_chain.csv",
      durum2019 = "polyploid_durum2019_hexaploid_vs_tetraploid.csv",
      cotton2023 = "polyploid_cotton2023_kasp_screen.csv"
    )[[nm]])
    expect_false(
      any(grepl("pooled|overall|combined|across_crops", names(d2))),
      label = nm
    )
  }
})

test_that("dosage: the mixture-model class count differs by ploidy, so does the rate", {
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  d <- ext("polyploid_durum2019_hexaploid_vs_tetraploid.csv")
  # Jaimes fits a decaploid mixture: 11 classes for one biallelic SNP
  expect_true(grepl("11 dosage categories", unique(j$ploidy_detail), fixed = TRUE))
  # durum's contrast is AABBDD-in-design vs AABB-catalogue on the SAME panel;
  # the design ploidy changes the number the assay can even report
  expect_equal(d$ploidy, c("allo-hexaploid", "allo-tetraploid", "mixed"))
  # the hexaploid assay on a tetraploid panel and the tetraploid catalogue are
  # NOT the same question: only the durum catalogue is homoeolog-resolved
  expect_equal(d$homeolog_handled, c("no", "yes", "mixed"))
  # so a "polyploid conversion rate" is already two different denominators here
  hex <- d[d$ploidy == "allo-hexaploid", ]
  tet <- d[d$ploidy == "allo-tetraploid", ]
  expect_gt(
    tet$polymorphic_pct_of_nonfailed - hex$polymorphic_pct_of_nonfailed,
    15
  )
  # and the cell nobody may quote as a rate: the combined row
  expect_equal(d$ploidy[[3L]], "mixed")
})

## --- 3. "a" / "ab" / "aab" are per-ploidy states, not per-species ----------

test_that("dosage: a dose is defined against the ploidy of the panel, not the species", {
  j <- ext("polyploid_jaimes2024_sugarcane_dosage.csv")
  c <- ext("polyploid_cotton2023_kasp_screen.csv")
  d <- ext("polyploid_durum2019_hexaploid_vs_tetraploid.csv")
  # the three panels do not even share a ploidy word
  expect_equal(unique(j$ploidy), "auto-polyploid")
  expect_equal(unique(c$ploidy), "allo-tetraploid")
  expect_equal(unique(d$ploidy), c("allo-hexaploid", "allo-tetraploid", "mixed"))
  # the sugar cane detail names the level a dose is a fraction OF
  expect_true(grepl("~10x", unique(j$ploidy_detail), fixed = TRUE))
  # "ab" on an AADD panel is two homoeologs, not two homologous copies; the
  # data says nothing about which, so the file must not claim it does
  expect_true(all(c$homeolog_handled == "not stated"))
  # the only source that could name hexaploid vs tetraploid dosage is durum,
  # and its polymorphism counts are assay OUTCOMES, not doses - assert that
  # the vocabulary which would make per-ploidy dose semantics available is
  # ABSENT rather than invent a tolerance for it
  expect_false(any(grepl(DOSAGE_WORDS, names(d), ignore.case = TRUE)))
})

test_that("dosage: the hexaploid/tetraploid contrast carries no dose number", {
  d <- ext("polyploid_durum2019_hexaploid_vs_tetraploid.csv")
  # a hexaploid has 7 dose states per biallelic SNP and a tetraploid 5; the
  # durum paper reports none of them, so the test pins the absence explicitly
  expect_equal(7L - 5L, 2L)
  # every durum row defines its outcome as the array's call class, not a dose
  expect_true(all(grepl("Array assay OUTCOME", d$outcome_definition, fixed = TRUE)))
  # neither the failed nor the polymorphic count is a dosage class count
  expect_true(all(c("n_failed", "n_polymorphic") %in% names(d)))
  expect_false(any(c("n_dosage_classes", "dose_0", "dose_1") %in% names(d)))
})

## --- 4. cotton is a NEGATIVE datapoint and is never averaged in ------------

test_that("dosage: the cotton screen is logged as a negative datapoint", {
  c <- ext("polyploid_cotton2023_kasp_screen.csv")
  # its outcome is cluster separation, and the file says it is NOT a dose
  expect_true(all(grepl("cluster separation", c$outcome_definition, fixed = TRUE)))
  expect_true(all(grepl("not a SNP-to-KASP conversion", c$outcome_definition, fixed = TRUE)))
  # a screen that reports 8 usable of 11 cannot stand in for a dosage rate
  expect_equal(sum(c$usable), 8L)
  expect_equal(nrow(c), 11L)
  # and it is the only source whose usable assays include a deliberately
  # degraded class (GOOD* = best as a dominant marker), which a dosage rate
  # has no room for
  expect_equal(sum(c$outcome_class == "GOOD*"), 1L)
  expect_true(all(c$usable[c$outcome_class == "GOOD*"]))
})

## --- 5. the axis must be NEW, and pop_freq must not quietly stand in -------

test_that("dosage: marker_population_checks() today is a two-allele frequency table", {
  pkgload::load_all(test_path("../.."), quiet = TRUE)
  pop <- tibble::tibble(
    alt_freq = 0.5, min_alt_freq = 0.01,
    n_haplotypes = 2L, n_alt_haplotypes = 1L
  )
  mk <- structure(
    list(
      population = pop, snp = tibble::tibble(pos = 100L),
      haplotypes = NULL, haplotype_hits = NULL
    ),
    class = c("marker_kasp", "marker")
  )
  chk <- tidymarkers:::marker_population_checks(mk)
  expect_true(all(c("level", "name", "status", "value") %in% names(chk)))
  # TODAY: no dosage row exists.  When Z10.2 lands the separate axis this
  # expectation is the one to flip - the test is the task's acceptance gate.
  expect_false(any(grepl("dos|copy|ploid", tolower(c(chk$name, chk$value)))))
  # and the frequency row cannot tell 0.5 of 10 copies from 0.5 of 2 copies:
  # it prints "50%" either way, which is exactly why dosage is a separate axis
  expect_equal(chk$value[chk$name == "ALT freq"], "50%")
  expect_false("dose" %in% names(pop))
})

test_that("dosage: the sources' own vocabulary is pooled-safe only per source", {
  # a table that would pool the four sources must not exist yet: the roadmap
  # names Z10.2 as a NEW axis, and a pre-existing pooled dosage table would
  # mean an earlier pass averaged incomparable outcomes
  pooled <- list.files(
    test_path("../../inst/extdata"),
    pattern = "^polyploid.*dosage.*\\.csv$", full.names = FALSE
  )
  expect_true(all(grepl("jaimes2024", pooled)))
  # every other dosage-worded file would be a candidate average; there is none
  expect_length(setdiff(pooled, "polyploid_jaimes2024_sugarcane_dosage.csv"), 0L)
})
