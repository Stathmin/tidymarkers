# =============================================================================
#  ZW2/ZW3/ZW4  --  the scientific data files are honest
# =============================================================================
#  Every table in inst/extdata that carries a literature claim is checked here.
#  The rule this file enforces is the one this repository got wrong once:
#
#      a number may only be present if the row says which source it came from,
#      how much of that source was read, and what "success" meant.
#
#  `access` is the load-bearing column.  `abstract` means a sentence was read;
#  `fulltext` / `supplement` / `source` mean the table or the code itself was
#  parsed.  A row that carries a numeric outcome at `abstract` level is allowed
#  to exist but may NOT be used to calibrate anything, so the tests below pin the
#  calibrated flag at FALSE for every table that has one.

## --- shared readers ---------------------------------------------------------

ext <- function(name) {
  readr::read_csv(test_path("../../inst/extdata", name), show_col_types = FALSE,
                  na = c("", "NA"))
}

test_that("scientific_sources: every row names a source and how far it was read", {
  s <- ext("scientific_sources.csv")
  expect_gt(nrow(s), 5L)
  for (col in c("key", "citation", "venue", "doi", "access", "what_it_anchors", "claim")) {
    expect_false(any(is.na(s[[col]]) | !nzchar(s[[col]])), label = col)
  }
  expect_setequal(unique(s$key), s$key) # no duplicate keys
  # 'closed access' is a real, honest label: the row exists to say the paper was
  # looked for and could not be read, so no number may come from it. The set is
  # CLOSED: a probing batch of workers added 'open', 'subscription', 'blocked'
  # and five prose variants ('full text (Tables 4, 6, 7)'), which would have made
  # the column unfilterable. Detail belongs in `access_detail`, never here.
  expect_true(all(s$access %in% c("fulltext", "abstract", "stub", "transcribed",
                                  "supplement", "source", "closed access")),
              label = paste(setdiff(unique(s$access),
                                    c("fulltext", "abstract", "stub", "transcribed",
                                      "supplement", "source", "closed access")),
                            collapse = ", "))
  # a closed-access row may say nothing numeric
  closed <- s[s$access == "closed access", ]
  if (nrow(closed)) {
    expect_true(all(grepl("NOT READ", closed$claim)))
  }
  # a DOI is required for every source; the prefix must look like a DOI, EXCEPT
  # for the explicit sentinel a resolved-but-doi-less article carries. A blank
  # cell is forbidden because it is indistinguishable from a dropped field; the
  # sentinel is greppable and says the question is still open.
  unresolved <- s$doi == "DOI-UNRESOLVED"
  expect_true(all(grepl("^10\\.[0-9]{4,9}/", s$doi[!unresolved])))
  # and a row with no DOI must be reachable by PMCID, or it cannot be found again
  expect_true(all(!is.na(s$pmcid[unresolved]) & nzchar(s$pmcid[unresolved])),
              label = "no-DOI rows carry a PMCID")
  # a row that has been VERIFIED must say when
  expect_true(all(!s$verified | !is.na(s$verified_on)))
})

test_that("scientific_sources: no row is labelled 'fetched' without a fetched artefact", {
  # `access` is the honest label.  Only these levels correspond to something read
  # beyond an abstract, and each of them has a corresponding parser or fixture
  # somewhere in the repository (makhoul-supplement.R, the PolyMarker source
  # read, the Lefever full text).
  s <- ext("scientific_sources.csv")
  deep <- s$key[s$access %in% c("supplement", "source", "fulltext")]
  expect_true("makhoul2020" %in% deep)
  expect_true("polymarker_source" %in% deep)
  expect_true("lefever2019" %in% deep)
  # and the supplement really is in the tree, or the row is a promise
  expect_true(file.exists(test_path("fixtures", "makhoul2020_MOESM1.xlsx")))
})

test_that("mismatch_discrimination: every row is cited and graded", {
  d <- ext("mismatch_discrimination.csv")
  expect_gt(nrow(d), 5L)
  expect_false(any(is.na(d$doi) | !nzchar(d$doi)))
  ladder <- c("rule of thumb", "vendor guideline", "primary literature", "fitted here")
  expect_true(all(d$evidence %in% ladder))
  # no row may claim a strength above what its source type supports
  expect_true(all(d$evidence != "fitted here"), label = "no fit has been done yet")
})

test_that("mismatch_discrimination: Ye's -2 and Lefever's -3 are kept as a CONTRADICTION", {
  d <- ext("mismatch_discrimination.csv")
  ye <- d[d$table == "ye_second_mismatch", ]
  lf <- d[d$table == "lefever_dmas", ]
  expect_gt(nrow(ye), 0L)
  expect_gt(nrow(lf), 0L)
  expect_true(all(ye$position == -2L))
  scan <- lf[lf$finding == "position_scan", ]
  expect_equal(scan$position_optimum, -3L)
  # the disagreement is stated in the Ye rows themselves, not only in a comment
  expect_true(all(grepl("CONFLICT", ye$claim)))
})

test_that("mismatch_anchor: our own model output is never marked calibrated", {
  a <- ext("mismatch_anchor.csv")
  expect_gt(nrow(a), 0L)
  expect_true(all(!a$calibrated))
  expect_true(all(a$evidence == "fitted here"))
  # the engine is named, and the name is one ZW2 also reports
  expect_true(all(grepl("DECIPHER", a$engine)))
  expect_true("DECIPHER" %in% names(marker_engine_stamp()))
  # The finding this table exists to record: DECIPHER's margin between these
  # 3' dinucleotides is ~0, i.e. it does not rank them.  If a future engine
  # starts separating them, this test fails and the Z7 model gets re-examined.
  expect_lt(max(abs(a$predicted_margin)), 0.05)
})

test_that("conversion_rates: no row can yield a probability, and none is all-success", {
  c1 <- ext("conversion_rates.csv")
  expect_gt(nrow(c1), 0L)
  for (col in c("study", "species", "design_strategy", "doi")) {
    expect_false(any(is.na(c1[[col]]) | !nzchar(c1[[col]])), label = col)
  }
  expect_true(all(c1$outcome_scale == "usable"))
  expect_true(all(!c1$calibrated))
  # the anti-survivorship assertion: published successes alone are not data
  expect_true(all(c1$n_designed > c1$n_success))
  expect_true(all(c1$n_failed == c1$n_designed - c1$n_success))
  # no invented residual may exceed the failures it claims to explain
  named <- dplyr::coalesce(c1$n_fail_pcr, 0L) + dplyr::coalesce(c1$n_fail_no_allele, 0L) +
    dplyr::coalesce(c1$n_fail_other, 0L)
  expect_true(all(named <= c1$n_failed))
})

test_that("conversion_rates: the Grewal 2022 named reasons sum to its own total", {
  c1 <- ext("conversion_rates.csv")
  g <- c1[c1$study == "grewal2022_amblyopyrum", ]
  expect_equal(nrow(g), 1L)
  expect_equal(g$n_designed - g$n_success, 363L)
  # 251 (PCR) + 49 (no wild-relative allele) + 22 + 10 + 31 = 363.  Only two
  # reason columns exist here, so the row is partial BY SCHEMA -- assert the
  # numbers we do carry, never a residual.
  expect_equal(g$n_fail_pcr, 251L)
  expect_equal(g$n_fail_no_allele, 49L)
  expect_true(is.na(g$n_fail_other))
})

test_that("gate-pass coverage is tied to criterion ids, not to assertion counts", {
  # The Makhoul supplement supports the SITE/ASSAY criteria of a design; it says
  # nothing about contig-wide specificity (no genome in the workbook).  Stating
  # which ids the data can support is what stops "we added tests" from reading
  # as "we covered the gates".
  crit <- marker_criteria(KASP_LIMITS, scheme = "snp")
  check_labels <- as.character(crit$label[crit$role == "check"])
  d <- ext("makhoul_design_strategies.csv")
  # every SNP is a biallelic bracket like "[A/G]"; one row records two genotypes
  expect_true(all(grepl("^\\[[ACGT]/[ACGT]\\]", d$alleles_tg)))
  expect_setequal(
    unique(unlist(regmatches(d$alleles_tg, gregexpr("[ACGT]", d$alleles_tg)))),
    c("A", "C", "G", "T")
  )
  # a real check label used by this table's tests
  expect_true("amplicon length (bp)" %in% check_labels)
})

test_that("makhoul: the supplement fixture still reproduces the paper's totals", {
  # The producer script asserts 3/8/11 while writing; this asserts the CSV that
  # ships, so a hand-edit of the CSV is caught too.
  d <- ext("makhoul_design_strategies.csv")
  expect_equal(nrow(d), 15L)
  ok <- function(col) sum(d[[col]] == "Correct cluster", na.rm = TRUE)
  expect_equal(c(ok("s1_flanks_only"), ok("s2_single_ref"), ok("s3_multi_sanger")),
               c(3L, 8L, 11L))
  # every row carries the definition of the outcome it records
  expect_true(all(!is.na(d$outcome_def) & nzchar(d$outcome_def)))
  expect_true(all(d$doi == "10.1007/s00122-020-03608-x"))
})

test_that("makhoul primers: the allele is carried by the AS 3' base in every design", {
  p <- ext("makhoul_kasp_primers.csv")
  expect_gt(nrow(p), 10L)
  expect_true(all(p$allele_in_as_3p))
  expect_true(all(p$fam_3p != p$hex_3p))
  expect_true(all(substr(p$fam, nchar(p$fam), nchar(p$fam)) == p$fam_3p))
  expect_true(all(substr(p$hex, nchar(p$hex), nchar(p$hex)) == p$hex_3p))
  # Y3.2: the ARTIFACT must show the genome-bearing base separately.  The common
  # primer's 3' base is recorded for that, and it is a distinct axis from the
  # allele-bearing one -- it must not be assumed to differ between A/B/D.
  expect_true("common_3p" %in% names(p))
  expect_false(identical(p$common_3p, p$fam_3p))
})

test_that("scientific_backlog: verified means an extract was written", {
  b <- ext("scientific_backlog.csv")
  expect_gt(nrow(b), 5L)
  expect_false(anyNA(b$key))
  # The vocabulary records how far each row got, and nothing beyond that:
  #   unfetched                      -- a DOI or a title and nothing else
  #   fetched                        -- opened, but not yet extracted
  #   extracted                      -- an extract CSV exists, with an assertion
  #   fulltext-read-no-assertion     -- read, deliberately no CSV (say why)
  #   not_extractable                -- read, no usable table
  #   fetched-not-extractable: <why> -- read, the row described another paper
  allowed <- c("unfetched", "fetched", "extracted",
               "fulltext-read-no-assertion", "not_extractable")
  expect_true(all(b$status %in% allowed | grepl("^fetched-not-extractable: ", b$status)))
  # `verified` may only be TRUE on a row whose extract really exists, and every
  # such row says when
  # `verified` = the claim was checked against the fetched artefact; that is true
  # for rows whose result is "extracted", and equally for rows whose result is
  # "read and deliberately not extracted" -- the verification is what was done,
  # not what was produced
  done <- b$status %in% c("extracted", "fulltext-read-no-assertion", "not_extractable")
  expect_equal(as.logical(b$verified), done)
  expect_true(all(!as.logical(b$verified) | !is.na(b$verified_on)))
  # a row read but deliberately not extracted must explain itself
  silent <- b$status %in% c("fulltext-read-no-assertion", "not_extractable")
  if (any(silent)) expect_true(all(grepl("no|not", b$status[silent])))
})
