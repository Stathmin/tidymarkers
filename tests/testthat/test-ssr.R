# SSR / microsatellite as a marker scheme: one flanking pair over the repeat.

test_that("design_ssr_marker returns a clean flanking pair outside the repeat", {
  # a random 600 bp template with an explicit (AC)5 repeat inserted at 296, so
  # the designer has a guaranteed clean, located SSR
  W <- bseq(600L, 7L)
  at <- 296L
  ssr <- "ACACACACAC"
  W <- paste0(substr(W, 1L, at - 1L), ssr, substr(W, at + nchar(ssr), 600L))
  pos <- at + nchar(ssr) %/% 2L
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = "A", alt = "G")
  m <- design_ssr_marker(genome, snp, flank = 250L)
  rep_span <- m[["repeat"]]
  expect_s3_class(m, "kasp_ssr")
  # the repeat span includes the anchor, and both primers flank the repeat
  expect_lte(rep_span$start, pos)
  expect_gte(rep_span$end, pos)
  expect_lt(m$primers$start + nchar(m$primers$forward), rep_span$start)
  expect_gt(m$primers$end - nchar(m$primers$reverse), rep_span$end)
  # the size-marker amplicon window, and clean primers (no run / no SSR)
  expect_gte(m$primers$product, 100L)
  expect_lte(m$primers$product, 300L)
  expect_true(grepl("^[ACGT]+$", m$primers$forward) && grepl("^[ACGT]+$", m$primers$reverse))
  expect_false(any(has_run(c(m$primers$forward, m$primers$reverse))))
  expect_false(any(has_repeat(c(m$primers$forward, m$primers$reverse))))
  expect_output(print(m), "kasp_ssr")
})

test_that("design_ssr_marker declines when no repeat is in the window", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(600L, 11L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "G")
  # 20 tandem copies cannot occur by chance in a 201 bp window
  expect_message(
    expect_null(design_ssr_marker(genome, snp, flank = 100L, min_copies = 20L)),
    "No SSR"
  )
})

test_that("print.kasp_ssr reports the repeat and the amplicon", {
  x <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 10L, ref = "A", alt = "G"),
      primers = list(
        forward = "ACGTACGTACGTACGTAC", reverse = "TGCATGCATGCATGCATG",
        start = 1L, end = 120L, product = 120L
      ),
      amplicon = tibble::tibble(start = 1L, end = 120L)
    ),
    class = "kasp_ssr"
  )
  x[["repeat"]] <- tibble::tibble(start = 40L, end = 52L)
  expect_output(print(x), "kasp_ssr")
  expect_output(print(x), "repeat")
  expect_output(print(x), "120 bp")
})

test_that("the ssr scheme definition carries its gates, tiers and labels", {
  d <- scheme_def_ssr()
  expect_identical(d$scheme, "ssr")
  expect_identical(d$family, "ssr")
  expect_identical(d$loci, 1L)
  expect_true(all(c(
    "gate_ssr_amplicon", "gate_ssr_flanks", "gate_ssr_clean",
    "gate_ssr_tm", "gate_ssr_offtarget"
  ) %in% d$gates))
  expect_true(all(d$gates %in% names(d$gate_short)))
  expect_true(all(d$gates %in% names(d$gate_limit)))
  expect_true(all(d$relax %in% d$gates))
  expect_true(all(names(d$gate_floor_label) %in% d$relax))
  expect_true(all(names(d$gate_hard_label) %in% setdiff(d$gates, d$relax)))
  expect_setequal(names(d$tiers), c("tier1", "tier2", "tier3"))
  expect_true(is.character(d$stage_label) && nzchar(d$stage_label))
})
