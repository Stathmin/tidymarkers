# Acceptance checks: a good assay passes everything, each broken criterion fails.
# Synthetic assay objects keep this pure (no genome, no network, no DECIPHER).

good_assay <- function() {
  p <- kasp_primers()
  list(
    assay = tibble::tibble(
      amplicon = 90L, as_ref = p[["AS_ref"]], as_alt = p[["AS_alt"]], com = p[["COM"]],
      com_tm = 60, as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L,
      as_n_var = 0L, as_dropout = 0, as_foot_dropout = 0, com_dropout = 0
    ),
    thermo = tibble::tibble(
      name = names(p), tm = c(57, 57, 60), gc = c(0.5, 0.5, 0.5),
      hairpin_dg = c(-1, -1, -1)
    ),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3), b = rep(names(p), times = 3),
      dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
    ),
    amplicons = tibble::tibble(
      primer = c("AS_ref", "AS_alt"),
      n_products = c(1L, 0L), width = c(90L, NA_integer_)
    ),
    # the intended locus is one exact hit per primer; without this the
    # off-target check has no data (unknown), not a silent pass
    ideogram = list(hits = tibble::tibble(
      primer = c("AS_ref", "AS_alt", "COM")
    )),
    efficiency = tibble::tibble(discrimination = 0.9)
  )
}

test_that("a well-formed assay passes every acceptance check", {
  ch <- marker_design_checks(good_assay())
  expect_true(all(ch$pass))
  expect_true(all(ch$status == "pass"))
  expect_true(all(c("criterion", "observed", "limit", "pass", "status") %in% names(ch)))
  expect_equal(evaluation_verdict(ch), "usable")
})

test_that("check statuses keep pass/fail, and missing data is unknown not pass", {
  ch <- marker_design_checks(good_assay())
  expect_s3_class(ch$status, "factor")
  expect_setequal(levels(ch$status), c("pass", "fail", "unknown", "not_applicable"))
  # `pass` is the logical view of `status` - an unknown is not a success
  expect_equal(ch$pass, ch$status == "pass")
  # no discrimination evidence -> unknown (and pass FALSE), not a pass
  a <- good_assay()
  a$efficiency$discrimination <- NA_real_
  ch <- marker_design_checks(a)
  i <- ch$criterion == "discrimination"
  expect_equal(as.character(ch$status[i]), "unknown")
  expect_false(ch$pass[i])
  expect_equal(evaluation_verdict(ch), "unknown")
  # no ideogram -> the off-target count is unknown, never a silent pass
  a <- good_assay()
  a$ideogram <- NULL
  ch <- marker_design_checks(a)
  i <- ch$criterion == "off-target exact primer hits (n)"
  expect_equal(as.character(ch$status[i]), "unknown")
  expect_false(ch$pass[i])
})

test_that("the evaluation constructor owns the status policy", {
  out <- marker_evaluation_table(
    stage = c("alignment", "feature", "feature"),
    criterion = c("amplicon length (bp)", "target site class", "worst hairpin dG37"),
    observed = c(NA, NA, NA),
    limit = c("[50, 100]", NA, "-3"),
    pass = c(TRUE, TRUE, NA),
    detail = c(NA, "biallelic", NA),
    level = c("alignment", "feature", "feature"),
    not_applicable = c(FALSE, FALSE, TRUE),
    textual = "target site class"
  )
  # missing data beats a raw TRUE; a textual check is exempt; a family check
  # that does not apply is not_applicable rather than unknown
  expect_equal(as.character(out$status), c("unknown", "pass", "not_applicable"))
  expect_equal(out$pass, c(FALSE, TRUE, FALSE))
})

test_that("evaluation_verdict dooms on a fail and ignores not_applicable", {
  # a failed criterion dooms
  a <- good_assay()
  a$efficiency$discrimination <- 0.1
  expect_equal(evaluation_verdict(a), "doomed")
  # a not_applicable row is not a failure and does not doom (nor can it certify
  # on its own): every applicable check still passes -> usable
  ch <- marker_design_checks(good_assay())
  ch$status[ch$criterion == "discrimination"] <- marker_check_status("not_applicable")
  ch$pass <- marker_check_pass(ch$status)
  expect_equal(evaluation_verdict(ch), "usable")
  # ... but a table with nothing evaluated is unknown, not usable
  na <- ch
  na$status <- marker_check_status(rep("not_applicable", nrow(na)))
  na$pass <- marker_check_pass(na$status)
  expect_equal(evaluation_verdict(na), "unknown")
  # legacy tables without `status` still honour an NA pass as unknown
  expect_equal(evaluation_verdict(tibble::tibble(pass = NA)), "unknown")
  expect_equal(evaluation_verdict(tibble::tibble(pass = FALSE)), "doomed")
  expect_equal(evaluation_verdict(tibble::tibble(pass = TRUE)), "usable")
})

test_that("family checks are not_applicable to another marker class", {
  # a CAPS-classed object graded as LAMP: the LAMP-only checks do not apply
  p <- kasp_primers()
  x <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "G"),
      amplicon = tibble::tibble(start = 200L, end = 350L),
      primers = list(forward = p[["AS_ref"]], reverse = p[["COM"]]),
      ideogram = list(hits = tibble::tibble(
        primer = c("forward", "reverse")
      ))
    ),
    class = "marker_caps"
  )
  ch <- marker_design_checks(x, scheme = "lamp")
  fam <- ch$criterion %in% MARKER_CHECK_LAMP
  expect_true(all(ch$status[fam] == "not_applicable"))
  expect_false(any(ch$pass[fam]))
  expect_false(any(ch$status == "fail"))
  expect_equal(evaluation_verdict(ch), "usable")
})

test_that("each broken criterion flips the verdict to fail", {
  mutators <- list(
    function(a) {
      a$assay$amplicon <- 150L
      a
    },
    function(a) {
      a$assay$com_tm <- 50
      a
    },
    function(a) {
      a$assay$as_gc <- 0.2
      a
    },
    function(a) {
      a$thermo$hairpin_dg <- c(-1, -1, -6)
      a
    },
    function(a) {
      a$complementarity$dg[3] <- -15
      a
    },
    function(a) {
      a$complementarity$dg[1] <- -15
      a
    },
    function(a) {
      a$efficiency$discrimination <- 0.1
      a
    },
    function(a) {
      a$amplicons$n_products <- c(2L, 0L)
      a
    },
    function(a) {
      a$amplicons$n_products <- c(1L, 1L)
      a
    },
    function(a) {
      a$assay$com <- "ACGTACGTACGTACGT"
      a
    }
  )
  for (m in mutators) {
    expect_false(all(marker_design_checks(m(good_assay()))$pass))
  }
})

test_that("COM Tm must not fall below the AS primers (EasyKASP)", {
  a <- good_assay()
  a$thermo$tm <- c(62, 62, 60)
  a$assay$com_tm <- 60
  ch <- marker_design_checks(a)
  expect_false(ch$pass[ch$criterion == "COM Tm >= AS Tm"])
})

test_that("the amplicon window is the scheme's own, not the KASP optimum (C2/C5)", {
  amp <- function(ch) ch$criterion == "amplicon length (bp)"
  a <- good_assay()
  # a wet-lab CAPS band (120-1200 bp) passes as a caps marker, fails as a SNP
  a$assay$amplicon <- 600L
  caps <- marker_design_checks(a, scheme = "caps")
  expect_true(caps$pass[amp(caps)])
  expect_equal(caps$limit[amp(caps)], "[120, 1200]")
  snp <- marker_design_checks(a, scheme = "snp")
  expect_false(snp$pass[amp(snp)])
  # 600 bp is inside the Sanger window (400-800) but outside SSR's (100-300)
  sanger <- marker_design_checks(a, scheme = "sanger")
  expect_true(sanger$pass[amp(sanger)])
  expect_equal(sanger$limit[amp(sanger)], "[400, 800]")
  ssr <- marker_design_checks(a, scheme = "ssr")
  expect_false(ssr$pass[amp(ssr)])
})

test_that("the inverted-repeat screen fails a self-symmetric primer", {
  ch <- marker_design_checks(good_assay())
  expect_true("inverted repeat (arms>=4, loop<=6)" %in% ch$criterion)
  a <- good_assay()
  a$assay$com <- "ACGTACGTACGTACGT"
  bad <- marker_design_checks(a)
  expect_false(bad$pass[bad$criterion == "inverted repeat (arms>=4, loop<=6)"])
})
