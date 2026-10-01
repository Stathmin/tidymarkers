# Acceptance checks: a good assay passes everything, each broken criterion fails.
# Synthetic assay objects keep this pure (no genome, no network, no DECIPHER).

good_assay <- function() {
  p <- c(
    AS_ref = "CATGGACCTAGCATTGACGTA",
    AS_alt = "CATGGACCTAGCATTGACGTT",
    COM = "GCTTAGCCATGGTCAATGCTA"
  )
  list(
    assay = tibble::tibble(
      amplicon = 90L, as_ref = p[["AS_ref"]], as_alt = p[["AS_alt"]], com = p[["COM"]],
      com_tm = 60, as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L
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
    efficiency = tibble::tibble(discrimination = 0.9)
  )
}

test_that("a well-formed assay passes every acceptance check", {
  ch <- kasp_design_checks(good_assay())
  expect_true(all(ch$pass))
  expect_true(all(c("criterion", "observed", "limit", "pass") %in% names(ch)))
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
    expect_false(all(kasp_design_checks(m(good_assay()))$pass))
  }
})

test_that("COM Tm must not fall below the AS primers (EasyKASP)", {
  a <- good_assay()
  a$thermo$tm <- c(62, 62, 60)
  a$assay$com_tm <- 60
  ch <- kasp_design_checks(a)
  expect_false(ch$pass[ch$criterion == "COM Tm >= AS Tm"])
})

test_that("the inverted-repeat screen fails a self-symmetric primer", {
  ch <- kasp_design_checks(good_assay())
  expect_true("inverted repeat (arms>=4, loop<=6)" %in% ch$criterion)
  a <- good_assay()
  a$assay$com <- "ACGTACGTACGTACGT"
  bad <- kasp_design_checks(a)
  expect_false(bad$pass[bad$criterion == "inverted repeat (arms>=4, loop<=6)"])
})
