# =============================================================================
#  test-polyploid-scope
# =============================================================================
#  Cross-source invariants for the polyploid breadth pass.  Each polyploid CSV
#  speaks its own vocabulary; these tests are the only place the family is read
#  together, and what they protect is the COMPARABILITY the brief asks for:
#    * `ploidy` and `homeolog_handled` are present and drawn from a closed
#      vocabulary in every file (so the same failure mode is not "no" in one
#      paper and "not stated" in another for the same reason);
#    * `outcome_definition` names what "worked" means, and no two sources that
#      define "worked" differently are pooled;
#    * autopolyploids say "not applicable" for homeologs, NOT "not stated" -
#      an autopolyploid has no homoeolog axis to have missed;
#    * the 42-of-76 sugarcane rate is not comparable with the 8-of-11 cotton
#      rate, because one is dosage and the other is cluster separation;
#    * the wheat-specific sources this pass did NOT overturn are still named.

PLOIDY <- c("allo-hexaploid", "allo-tetraploid", "auto-polyploid", "mixed")
HOMEOLOG <- c("yes", "no", "not stated", "not applicable", "mixed")

POLY_FILES <- c(
  jaimes2024 = "polyploid_jaimes2024_sugarcane_dosage.csv",
  gao2022 = "polyploid_gao2022_sugarcane_chain.csv",
  durum2019 = "polyploid_durum2019_hexaploid_vs_tetraploid.csv",
  cotton2023 = "polyploid_cotton2023_kasp_screen.csv"
)

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("polyploid scope: every source carries ploidy and homeolog_handled", {
  for (nm in names(POLY_FILES)) {
    d <- ext(POLY_FILES[[nm]])
    expect_true(
      all(c(
        "ploidy", "homeolog_handled", "outcome_definition",
        "evidence_tier", "source_file", "key"
      ) %in% names(d)),
      label = nm
    )
    expect_true(all(d$ploidy %in% PLOIDY), label = nm)
    expect_true(all(d$homeolog_handled %in% HOMEOLOG), label = nm)
    expect_true(all(nzchar(d$outcome_definition)), label = nm)
    expect_true(all(nzchar(d$source_file)), label = nm)
    expect_true(all(d$evidence_tier %in% c("measured", "weak")), label = nm)
  }
})

test_that("polyploid scope: an autopolyploid has no homoeolog axis to miss", {
  j <- ext(POLY_FILES[["jaimes2024"]])
  g <- ext(POLY_FILES[["gao2022"]])
  expect_true(all(j$homeolog_handled == "not applicable"))
  # Gao does not state it, and that is a DIFFERENT statement from Jaimes's
  expect_true(all(g$homeolog_handled == "not stated"))
  expect_false(any(g$homeolog_handled == "not applicable"))
})

test_that("polyploid scope: cotton is allotetraploid and makes no homoeolog claim", {
  c <- ext(POLY_FILES[["cotton2023"]])
  expect_true(all(c$ploidy == "allo-tetraploid"))
  expect_true(all(c$homeolog_handled == "not stated"))
  # the two diploid-relative assays are the negative datapoint, so they must
  # exist and be counted as usable in the landed file
  expect_equal(sum(grepl("^i[0-9]+(Gt|Gm)$", c$marker)), 2L)
  expect_true(all(c$usable[grepl("^i[0-9]+(Gt|Gm)$", c$marker)]))
})

test_that("polyploid scope: 'worked' means different things and is not pooled", {
  j <- ext(POLY_FILES[["jaimes2024"]])
  c <- ext(POLY_FILES[["cotton2023"]])
  d <- ext(POLY_FILES[["durum2019"]])
  # three different outcome axes, three different phrases
  expect_true(all(grepl("dosage", j$outcome_definition)))
  expect_true(all(grepl("cluster separation", c$outcome_definition)))
  expect_true(all(grepl("Array assay OUTCOME", d$outcome_definition)))
  # each source's own denominator: 76 SNPs / 11 assays / 81,587 assays
  expect_true(all(j$n_tested %in% c(76, 42)))
  expect_true(all(c$n_tested == 11))
  expect_equal(d$n_assays[[3L]], 81587)
  # no file carries a column named like a pooled cross-source rate
  for (nm in names(POLY_FILES)) {
    expect_false(any(grepl(
      "pooled|overall_rate|combined_rate",
      names(ext(POLY_FILES[[nm]]))
    )), label = nm)
  }
})

test_that("polyploid scope: the wheat-only claims have no non-wheat support", {
  # the family's non-wheat sources must not be usable to justify any statement
  # about subgenome-specific (homeolog-aware) DESIGN - none of them reports one
  d <- ext(POLY_FILES[["durum2019"]])
  # durum's "yes" is about SNP ORIGIN (homoeolog-resolved catalogue), not design
  expect_equal(d$homeolog_handled, c("no", "yes", "mixed"))
  expect_true(grepl("homoeolog-resolved", d$homeolog_detail[[2L]]))
  # and the two allotetraploid KASP sources both say "not stated"
  c <- ext(POLY_FILES[["cotton2023"]])
  expect_true(all(c$homeolog_handled == "not stated"))
})
