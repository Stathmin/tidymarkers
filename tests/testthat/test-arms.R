# Tetra-primer ARMS-PCR: the outer pair, the inner allele-specific primers, the bands.

test_that("an outer flanking pair sits outside the region and codes a product", {
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  p <- tidymarkers:::design_flanking_pair(w, lo = 300L, hi = 300L)
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
  expect_s3_class(m, "marker_arms")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  b <- m$bands
  expect_setequal(b$band, c("control", "REF", "ALT"))
  expect_gte(abs(b$size[b$band == "REF"] - b$size[b$band == "ALT"]), 20L)
  expect_true(all(b$size > 0))
  expect_true(inherits(plot(m), c("patchwork", "ggplot")))
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
    class = "marker_arms"
  )
  p <- plot_arms_gel(marker)
  expect_s3_class(p, "ggplot")
  expect_true(all(p$data$col %in% c(GEL_DIAG, GEL_SHARED)))
  expect_false(any(p$data$col %in% unname(KASP_STATUS[c("neutral", "broken")])))
  expect_true(all(p$data$width > 0.2 & p$data$width <= 0.4))
})

test_that("the arms scheme is registered and yields a DAG with its own stage", {
  expect_true("arms" %in% marker_schemes()$scheme)
  crit <- marker_criteria(scheme = "arms")
  expect_true(all(paste0("arms:gate:", KASP_GATES_ARMS) %in% crit$id))
  s <- marker_decision_schema("arms")
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(ord[s$edges$to] > ord[s$edges$from]))
})

test_that("the arms DAG obeys the contract", {
  expect_dag_invariants("arms")
})

test_that("arms declares its gate set, hardness, both orders and real floors", {
  crit <- marker_criteria(scheme = "arms")
  gates <- crit[crit$role == "gate", ]
  expect_setequal(gates$id, paste0("arms:gate:", KASP_GATES_ARMS))
  # application order (grave / cheap first) is the declared vector
  expect_equal(gates$id[order(gates$order)], paste0("arms:gate:", KASP_GATES_ARMS))
  # only "no allele-specific inner primer" is hard (C8)
  expect_equal(gates$id[gates$hard], "arms:gate:gate_arms_disc")
  expect_setequal(gates$id[!gates$hard], paste0("arms:gate:", KASP_RELAX_ARMS))
  # the ladder is least grave first, and the band-gap floor is *below* its bar -
  # the audit's "floor equals threshold, drawn soft, can never relax" bug
  expect_equal(KASP_RELAX_ARMS, c("gate_arms_outer", "gate_arms_bands"))
  expect_lt(KASP_GATE_FLOOR_ARMS[["gate_arms_bands"]], 0)
  expect_lt(KASP_ARMS_GAP_FLOOR, KASP_ARMS_WINDOW$band_gap)
  # limits are real values, never a `[lo, hi]` placeholder (C9)
  expect_false(any(grepl("\\[lo, ?hi\\]", gates$limit)))
})

test_that("design_arms_marker evaluates the arms gates and gives a resolvable gel", {
  skip_if_no_oligoarrayaux()
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 800L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  ref <- substr(w, 400L, 400L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 400L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_arms_marker(genome, snp, flank = 250L)
  # the gates are *evaluated*, not decorative (C7): the marker carries the outcome
  expect_setequal(names(m$selection), c("relaxed", "violated"))
  expect_true(all(unlist(m$selection) %in% KASP_GATES_ARMS))
  # the control band is the largest product (the inner primers are nested), and
  # both diagnostic bands clear the accepted separation from each other and it
  b <- m$bands
  ctrl <- b$size[b$band == "control"]
  ref_size <- b$size[b$band == "REF"]
  alt_size <- b$size[b$band == "ALT"]
  expect_gt(ctrl, max(ref_size, alt_size))
  expect_gte(min(abs(ref_size - alt_size), ctrl - max(ref_size, alt_size)), KASP_ARMS_WINDOW$band_gap)
  expect_gte(m$outer$product, KASP_ARMS_WINDOW$control[1L])
  expect_lte(m$outer$product, KASP_ARMS_WINDOW$control[2L])
})
