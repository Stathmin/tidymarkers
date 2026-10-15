# Tetra-primer ARMS-PCR: the outer pair, the inner allele-specific primers, the bands.

test_that("an outer flanking pair sits outside the region and codes a product", {
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  p <- tidykasp:::design_flanking_pair(w, lo = 300L, hi = 300L)
  expect_true(p$start < 300L && p$end > 300L)
  expect_true(p$product >= 200L && p$product <= 800L)
})

test_that("design_arms_marker gives two resolvable, different-sized diagnostic bands", {
  skip_if_no_oligoarrayaux()
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 800L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  ref <- substr(w, 400L, 400L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 400L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_arms_marker(genome, snp, flank = 250L)
  expect_s3_class(m, "kasp_arms")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  b <- m$bands
  expect_setequal(b$band, c("control", "REF", "ALT"))
  expect_gte(abs(b$size[b$band == "REF"] - b$size[b$band == "ALT"]), 20L)
  expect_true(all(b$size > 0))
})

test_that("design_arms_marker declines non-substitution alleles", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 40L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "A", alt = "AA")
  expect_message(suppressWarnings(design_arms_marker(genome, snp, flank = 40L)), "single-base")
})

test_that("plot_arms_gel draws the AA / AB / BB readout", {
  marker <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 10L, ref = "A", alt = "G"),
      bands = tibble::tibble(
        band = c("control", "REF", "ALT"), size = c(400L, 200L, 250L)
      )
    ),
    class = "kasp_arms"
  )
  expect_s3_class(plot_arms_gel(marker), "ggplot")
})

test_that("the arms scheme is registered and yields a DAG with its own stage", {
  expect_true("arms" %in% kasp_schemes()$scheme)
  crit <- kasp_criteria(scheme = "arms")
  expect_true(all(paste0("arms:gate:", KASP_GATES_ARMS) %in% crit$id))
  s <- kasp_decision_schema("arms")
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(ord[s$edges$to] > ord[s$edges$from]))
})
