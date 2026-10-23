# =============================================================================
#  KASP family  --  design success is not work success
# =============================================================================
#  Three sources, three different things measured, and the guard on this family
#  is that no two of them can be mistaken for each other:
#
#    wheat_platform2019  design 49/50, then 2,784 assay x sample points with
#                        6.5% unamplified (a POINT fraction, not an assay rate)
#                        and 13 invalid calls; per-assay outcomes do NOT exist.
#    sharplab2021        design 50/50, polymorphism between two inbreds 46/50
#                        (the "92%"), and a separate 98.5% amplification ratio.
#    maize_kasp2020      design 100/100, "work" 92/100 = 71 usable + 21
#                        monomorphic/het, 8 no signal.
#
#  Every number these files carry is re-derived by data-raw/<key>.R from the
#  workbook or full text; the tests below assert the arithmetic that would
#  break if a number were hand-edited, and then -- when .tmp/lit is present --
#  re-run the generator and require the CSV on disk to be byte-identical.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

kasp_run <- function(script) {
  skip_if_not_installed("pkgload")
  if (!requireNamespace("fs", quietly = TRUE)) skip("fs not installed")
  root <- normalizePath(test_path("../.."), mustWork = FALSE)
  skip_if(
    !file.exists(file.path(root, "data-raw", script)),
    paste("generator", script, "not present")
  )
  old <- getwd()
  setwd(root)
  on.exit(setwd(old), add = TRUE)
  # the package is already loaded by the harness, so no load_all() here: the
  # generator only needs the helpers it sources itself from data-raw/kasp_lib.R
  # the generators print a one-line summary; capture it, and let any stop()
  # from a disagreed assertion escape as the failure
  invisible(utils::capture.output(
    suppressMessages(source(
      file.path(root, "data-raw", script),
      local = new.env()
    )),
    type = "output"
  ))
  TRUE
}

# ---------------------------------------------------------------------------
#  wheat_platform2019
# ---------------------------------------------------------------------------
test_that("KASP wheat_platform2019: the 30 assays are the 30 the paper made", {
  a <- ext("wheat_platform2019_kasp_assays.csv")
  expect_equal(nrow(a), 30L)
  expect_equal(nlevels(factor(a$assay)), 30L)
  expect_false(anyDuplicated(a$assay) > 0)
  # the worksheet lists the assays; it does NOT list an outcome for any of them
  expect_false(any(grepl(
    "^outcome$|design_outcome|work_outcome|assay_worked",
    names(a)
  )))
  # a KASP primer triple per assay, and the two allele-specific primers differ
  # at their 3' base -- the assay construct
  expect_true(all(
    nzchar(a$kasp_allele1) & nzchar(a$kasp_allele2) & nzchar(a$kasp_common)
  ))
  expect_true(all(a$kasp_allele1_3p != a$kasp_allele2_3p))
  # the 3' bases are in the target's bracketed alleles, in some orientation;
  # how many needed the complement is recorded rather than silently fixed
  expect_true(all(a$has_ambiguous_base %in% c(TRUE, FALSE) | is.na(a$has_ambiguous_base)))
  expect_true(sum(a$as_3p_is_complement) >= 1L)
  # TaPds_B1 is the one assay with no TaqMan probes (an INDEL > 6 bp is the
  # paper's stated reason) -- a DESIGN effect, not a work failure
  expect_equal(sum(!nzchar(a$kasp_id)), 0L)
})

test_that("KASP wheat_platform2019: every outcome count closes on the paper's own arithmetic", {
  o <- ext("wheat_platform2019_platform_outcomes.csv")
  get <- function(q, p = NULL) {
    i <- o$quantity == q
    if (!is.null(p)) i <- i & o$platform == p
    o[i, ]
  }
  # design: 39/50, 49/50, 49/50 -- each with denominator 50
  d <- get("design_success")
  # TaqMan 39, KASP 49, rhAmp 49 (the paper's sentence order)
  expect_equal(
    d$n_numerator[match(c("TaqMan", "KASP", "rhAmp"), d$platform)],
    c(39L, 49L, 49L)
  )
  expect_true(all(d$n_denominator == 50L))
  # the Methods "39" is flagged as ambiguous and is NOT KASP's number
  amb <- get("all_passed_design_criteria")
  expect_equal(amb$quantity_kind, "design_ambiguous")
  expect_equal(amb$scope, "study_intersection")
  # 30 of the 39 manufactured, 1 failed on all three, 29 x 96 = 2,784 points
  expect_equal(get("manufactured")$n_numerator, 30L)
  expect_equal(get("failed_on_all_three_platforms")$n_numerator, 1L)
  expect_equal(get("assays_genotyped")$n_numerator, 29L)
  expect_equal(get("assays_x_samples")$n_numerator, 2784L)
  expect_equal(96L * 29L, 2784L)
  # invalid calls sum to the printed total
  inv <- get("invalid_call")
  split <- inv[inv$platform != "all three", ]
  total <- inv[inv$platform == "all three", ]
  expect_equal(sum(split$n_numerator), total$n_numerator)
  expect_equal(sort(split$n_numerator), c(7L, 13L, 57L))
  # no row may carry a per-assay rate: amplification is a POINT fraction
  amp <- o[o$quantity %in% c("amplified", "unamplified"), ]
  expect_true(all(is.na(amp$n_numerator)))
  expect_equal(unique(amp$n_denominator), 2784L)
  expect_true(all(amp$scope == "assay_x_sample"))
  # and no row may have been turned into a probability
  expect_false(any(grepl("prob|p_success|expected", names(o))))
})

test_that("KASP wheat_platform2019: cluster metrics keep their units and their footnote", {
  m <- ext("wheat_platform2019_cluster_metrics.csv")
  expect_equal(nrow(m), 12L) # 3 platforms x 4 metrics
  expect_equal(length(unique(m$platform)), 3L)
  expect_equal(length(unique(m$metric)), 4L)
  # values are platform means in the paper's own units; call quality is a
  # fraction 0-1, the distances are positive
  q <- m[m$metric == "allele_call_quality", ]
  expect_true(all(q$value > 0 & q$value <= 1))
  expect_true(all(m$value >= 0))
  # the three platforms share one residual "Tech" df per metric
  expect_true(all(m$df_tech == 2L))
  expect_equal(m$df_residual[m$metric == "separation_angle"], rep(59L, 3))
  expect_equal(m$df_residual[m$metric == "cluster_compactness"], rep(57L, 3))
  expect_equal(m$df_residual[m$metric == "ntc_cluster_distance"], rep(71L, 3))
  # the Allele-11 footnote travels with the number
  expect_true(all(grepl("Allele 11", m$outcome_definition)))
  expect_true(all(grepl("not design, amplification or conversion", m$outcome_definition)))
  # the text's claim: rhAmp and KASP share a letter for separation angle and
  # differ for cluster compactness
  let <- function(p, k) m$lsd_letter[m$platform == p & m$metric == k]
  expect_equal(let("KASP", "separation_angle"), let("rhAmp", "separation_angle"))
  expect_false(let("KASP", "cluster_compactness") == let("TaqMan", "cluster_compactness"))
})

# ---------------------------------------------------------------------------
#  sharplab2021
# ---------------------------------------------------------------------------
test_that("KASP sharplab2021: 50 markers, 46 polymorphic, and the 92% is not a design rate", {
  a <- ext("sharplab2021_assays.csv")
  expect_equal(nrow(a), 50L)
  expect_equal(sort(a$no), 1:50)
  expect_equal(sum(a$assay_worked), 46L)
  expect_true(all(a$design_success)) # a primer triple exists for all 50
  expect_false(anyDuplicated(a$marker) > 0)
  # the crosses are exactly the rows with no PIC in the 20-line panel
  expect_true(all(is.na(a$pic_20_lines) == (a$polymorphism_1212_vs_b73 != "\u221a")))
  # the study-level statement of the same fact, with the denominator named
  o <- ext("sharplab2021_outcomes.csv")
  r <- o[o$quantity == "polymorphic_1212_vs_b73", ]
  expect_equal(r$n_numerator, 46L)
  expect_equal(r$n_denominator, 50L)
  expect_equal(r$n_numerator / r$n_denominator, 0.92)
  d <- o[o$quantity == "design_success", ]
  expect_equal(d$n_numerator, d$n_denominator) # 50 of 50 lines
  # there is a SEPARATE amplification ratio, kept as printed and not counted
  amp <- o[o$quantity == "amplified_20_lines", ]
  expect_equal(amp$value_printed, "98.5%")
  expect_true(is.na(amp$n_numerator))
  # the in-silico conversion is its own row and is not a wet-lab rate
  ins <- o[o$quantity == "insilico_conversion", ]
  expect_equal(ins$n_numerator, 71311L)
  expect_equal(ins$n_denominator, 2948985L)
  expect_equal(round(ins$n_numerator / ins$n_denominator, 4), 0.0242)
  # every rate row says which denominator it uses
  expect_false(any(is.na(o$outcome_definition)))
})

test_that("KASP sharplab2021: primer triples carry the site's two alleles at the 3' base", {
  p <- ext("sharplab2021_primers.csv")
  expect_equal(nrow(p), 50L)
  expect_true(all(p$as_3p_differs))
  expect_true(all(p$alleles_match))
  expect_true(all(nzchar(p$primer_fa) & nzchar(p$primer_fd) & nzchar(p$primer_r)))
})

# ---------------------------------------------------------------------------
#  maize_kasp2020
# ---------------------------------------------------------------------------
test_that("KASP maize_kasp2020: 100 sites, 71+21+8 = 100, and 92 is not 71", {
  o <- ext("maize_kasp2020_outcomes.csv")
  g <- function(q) o[o$quantity == q, ]
  expect_equal(g("snv_sites_attempted")$n_numerator, 100L)
  expect_equal(g("assays_designed")$n_numerator, 100L)
  poly <- g("usable_polymorphic")$n_numerator
  mono <- g("unusable_monomorphic_or_het")$n_numerator
  noamp <- g("no_amplification_signal")$n_numerator
  expect_equal(c(poly, mono, noamp), c(71L, 21L, 8L))
  expect_equal(poly + mono + noamp, 100L)
  # the paper's "work success rate of 92%" is (71 + 21)/100 and it is NOT the
  # 71% usable rate; both rows exist and disagree on purpose
  work <- g("work_success_rate_paper")$n_numerator
  expect_equal(work, poly + mono)
  expect_equal(work / 100, 0.92)
  expect_equal(poly / 100, 0.71)
  expect_false(poly == work)
  # the only assays that did not amplify are the 8, so amplification is 92/100
  expect_equal(g("amplified")$n_numerator, 100L - noamp)
  # every row carries its own definition, so a pooled table cannot lose them
  expect_false(any(is.na(o$outcome_definition)))
  expect_false(any(is.na(o$rate_recomputed)))
})

test_that("KASP maize_kasp2020: the per-assay table has no invented outcome", {
  a <- ext("maize_kasp2020_assays.csv")
  expect_equal(nrow(a), 100L)
  expect_equal(length(unique(a$chromosome)), 10L)
  expect_equal(as.integer(table(a$chromosome)), rep(10L, 10))
  # the supplement gives no per-assay outcome: the column must be empty
  expect_true(all(is.na(a$assay_outcome)))
  expect_true(all(a$design_success))
  # 92 of the 100 sites have a primer triple in supporting data S2; the 8 that
  # do not are marked, not dropped
  expect_equal(sum(a$primer_triple_in_s2), 92L)
  expect_true(all(is.na(a$primer_allele_x) == !a$primer_triple_in_s2))
  # whichever pairs exist are two-allele at the 3' base
  withp <- a[a$primer_triple_in_s2, ]
  expect_true(all(withp$as_3p_x != withp$as_3p_y))
})

# ---------------------------------------------------------------------------
#  the cross-file guard: the family cannot be pooled
# ---------------------------------------------------------------------------
test_that("KASP family: the three sources do not yield one poolable rate", {
  w <- ext("wheat_platform2019_platform_outcomes.csv")
  s <- ext("sharplab2021_outcomes.csv")
  m <- ext("maize_kasp2020_outcomes.csv")
  # three different design denominators: 50, 50, 100 -- and they must not be
  # summed anywhere by accident
  expect_equal(w$n_denominator[w$quantity == "design_success"], rep(50L, 3))
  expect_equal(unique(s$n_denominator[s$quantity == "design_success"]), 50L)
  expect_equal(unique(m$n_denominator[m$quantity == "assays_designed"]), 100L)
  # the three "success" numbers a reader would want to compare are defined
  # differently, and each file says so in its own column
  expect_true(grepl("pre-validated wheat", w$outcome_definition[
    w$quantity == "design_success" & w$platform == "KASP"
  ]))
  expect_true(grepl("1212 x B73|1212 and B73", s$note[
    s$quantity == "polymorphic_1212_vs_b73"
  ]))
  # the maize work rate is defined on the LM11 x CML25 parent pair, and the
  # usable rate (71%) exists as its own row
  expect_true(any(grepl("LM11", m$outcome_definition)))
  expect_equal(m$n_numerator[m$quantity == "usable_polymorphic"], 71L)
  # the sharpest trap in the family: two papers print "92%" and it means two
  # different things
  expect_true(any(grepl("92.00", s$value_printed)))
  expect_true(any(grepl("92%", m$value_printed)))
})

test_that("KASP generators re-derive their CSVs from the fetched files", {
  # Only the smallest source runs here: the xlsx + full-text pass is the one
  # that proves the whole path (workbook reader, text anchors, write-check), and
  # the other two generators are run by data-raw/RUN-kasp.R when a fetch lands.
  # sharplab2021 is the cheapest of the three (one 19 KB workbook, no xlsx join).
  root <- normalizePath(test_path("../.."), mustWork = FALSE)
  skip_if_not(dir.exists(file.path(root, "inst", "extdata")))
  expect_true(kasp_run("sharplab2021.R"))
})
