# A composed technology (roadmap Y1): a base chemistry + designer params +
# a specificity policy, first-class wherever a scheme is.

test_that("a composed technology is discovered with a valid DAG and its policy", {
  marker_technology("giga_snp",
    base = "snp", params = list(as_lengths = 24:30),
    specificity = "unique", label = "giga SNP KASP"
  )
  on.exit(rm("giga_snp", envir = .marker_tech_store), add = TRUE)

  expect_true("giga_snp" %in% marker_schemes()$scheme)
  expect_identical(marker_scheme("giga_snp")$family, "kasp")
  # the specificity policy is carried (Y3) and visible to the DAG contract
  expect_identical(marker_specificity("giga_snp"), "unique")
  expect_dag_invariants("giga_snp")

  t <- marker_technologies()
  expect_true("giga_snp" %in% t$scheme)
  expect_identical(t$base[t$scheme == "giga_snp"], "snp")
  expect_identical(t$specificity[t$scheme == "giga_snp"], "unique")
})

test_that("a composed technology carries its base's whole criterion set", {
  # A built-in stage is run by a *family*, and its gate/tier rows are tagged with
  # that family's schemes.  A technology that only re-parameterises a base is a
  # scheme of the same family, so it must see the same rows - not just the
  # acceptance checks.  When the membership misses it, `marker_criteria()` (and
  # every view built on it: the evaluation table, the decision text) silently
  # drops the whole selection stage.
  cols <- c("id", "role", "stage", "label", "level", "limit", "hard", "order")
  bases <- c("snp", "duplex", "arms")
  ids <- paste0("probe_", bases)
  for (i in seq_along(bases)) {
    marker_technology(ids[i], base = bases[i], overwrite = TRUE)
  }
  on.exit(rm(list = ids, envir = .marker_tech_store), add = TRUE)
  for (i in seq_along(bases)) {
    id <- ids[i]
    base <- bases[i]

    expect_true(id %in% marker_schemes()$scheme)
    # the same rows, in the same order, with the same bars as the base
    expect_equal(
      marker_criteria(scheme = id)[, cols],
      marker_criteria(scheme = base)[, cols],
      info = base
    )
    # and the same decision DAG, since the schema reads those rows
    expect_equal(
      marker_decision_schema(id)$nodes$name,
      marker_decision_schema(base)$nodes$name,
      info = base
    )
  }
})

test_that("a composed technology designs via its base designer", {
  skip_if_no_oligoarrayaux()
  marker_technology("giga_test", base = "snp", params = list())
  on.exit(rm("giga_test", envir = .marker_tech_store), add = TRUE)
  f <- snp_fixture()
  genome <- f$genome
  snp <- f$snp
  a <- suppressWarnings(marker_design_at(genome, snp, "giga_test"))
  expect_false(is.null(a))
})

test_that("marker_technology validates its arguments", {
  expect_error(marker_technology("", base = "snp"), "non-empty")
  expect_error(marker_technology("x", base = "nope"), "Unknown marker scheme")
  expect_error(marker_technology("snp"), "already exists")
  expect_error(marker_technology("y", base = "snp", specificity = "wat"), "Unknown specificity")
  expect_error(marker_technology("z", base = "snp", params = list(1, 2)), "named list")
})

test_that("marker_filter_blocks names the composable space", {
  b <- marker_filter_blocks()
  expect_true(all(c("base", "specificity", "params") %in% b$block))
  expect_setequal(b$values[b$block == "specificity"], paste(KASP_SPECIFICITY_MODES, collapse = " | "))
})

# ---- Z6: the population / haplotype layer reaches the verdict ----------------
# A well-formed KASP mock, optionally carrying a population.  `base` is the base
# each of the two haplotypes carries at the target (H1 REF, H2 usually ALT) and
# `freq` what they weigh.
pop_assay <- function(base, freq) {
  p <- kasp_primers()
  snp <- tibble::tibble(chrom = "Chr1", pos = 50L, ref = "A", alt = "G")
  marker <- list(
    scheme = "snp",
    assay = tibble::tibble(
      amplicon = 90L, as_ref = p[[1]], as_alt = p[[2]], com = p[[3]], com_tm = 60,
      as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L, as_n_var = 0L,
      as_dropout = 0, as_foot_dropout = 0, com_dropout = 0
    ),
    thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), hairpin_dg = c(-1, -1, -1)),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3), b = rep(names(p), times = 3),
      dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
    ),
    amplicons = tibble::tibble(primer = c("AS_ref", "AS_alt"), n_products = c(1L, 0L)),
    efficiency = tibble::tibble(discrimination = 0.9),
    # one exact hit per primer: the off-target check is evaluated, so the assay
    # on its own is `usable` and any change in the verdict is the population's
    ideogram = list(hits = tibble::tibble(primer = c("AS_ref", "AS_alt", "COM"))),
    selection = list(),
    snp = snp
  )
  if (is.null(base)) {
    return(marker)
  }
  haps <- tibble::tibble(
    haplotype = c("H1", "H2"), freq = freq, chrom = "Chr1", pos = 50L,
    ref = "A", alt = base
  )
  suppressWarnings(marker_attach_haplotypes(marker, haps, snp))
}

# the three biallelic populations these tests read: polymorphic, monomorphic
# (both haplotypes REF, so a *known* ALT frequency of 0) and an ALT haplotype
# whose frequency nobody supplied
poly_pop <- function() pop_assay(c("A", "G"), c(0.5, 0.3))
mono_pop <- function() pop_assay(c("A", "A"), c(0.5, 0.5))
unfreq_pop <- function() pop_assay(c("A", "G"), c(0.5, NA))

pop_rows <- function(ch) {
  ch[ch$criterion %in% c("ALT freq", "haplotypes", "ALT haplotypes", "hap covered"), ]
}

test_that("a polymorphic population reaches the verdict as pass", {
  m <- poly_pop()
  ch <- marker_design_checks(m)
  pop <- pop_rows(ch)
  # the population rows are checks, not a warning beside them
  expect_true(nrow(pop) >= 3L)
  expect_true(all(pop$status == "pass"))
  expect_true(all(pop$pass))
  # a pass never downgrades the verdict: the assay itself is usable and stays so
  expect_equal(evaluation_verdict(pop_assay(NULL)), "usable")
  expect_equal(evaluation_verdict(m), "usable")
})

test_that("a monomorphic site fails, whatever the marker's schema", {
  # both haplotypes carry REF at the target, so the ALT frequency is a *known*
  # zero: the assay cannot genotype anything and must not be recorded as pass
  m <- mono_pop()
  ch <- marker_design_checks(m)
  expect_equal(as.character(ch$status[ch$criterion == "ALT freq"]), "fail")
  expect_false(ch$pass[ch$criterion == "ALT freq"])
  expect_equal(evaluation_verdict(m), "doomed")
  # the warning is unchanged by the status (it stays additive)
  expect_warning(
    marker_attach_haplotypes(
      pop_assay(NULL),
      tibble::tibble(
        haplotype = c("H1", "H2"), freq = c(0.5, 0.5), chrom = "Chr1",
        pos = 50L, ref = "A", alt = c("A", "A")
      ),
      pop_assay(NULL)$snp
    ),
    "monomorphic"
  )
})

test_that("an absent ALT frequency is unknown, never a pass", {
  # the ALT haplotype carries no frequency: `marker_attach_haplotypes()` sums the
  # frequencies with `na.rm = TRUE`, so the stored `0` is not a measurement and
  # must neither certify nor doom the marker - the verdict drops to unknown
  m <- unfreq_pop()
  ch <- marker_design_checks(m)
  expect_equal(as.character(ch$status[ch$criterion == "ALT freq"]), "unknown")
  expect_false(ch$pass[ch$criterion == "ALT freq"])
  expect_equal(evaluation_verdict(m), "unknown")
})

test_that("the verdict reacts to a population stored after the evaluation", {
  # a marker whose `$evaluation` was stored before the population rows existed
  # must still not read as usable: the verdict folds the population in itself
  mono <- mono_pop()
  mono$evaluation <- marker_design_checks(pop_assay(NULL))
  expect_false("ALT freq" %in% mono$evaluation$criterion)
  expect_equal(evaluation_verdict(mono), "doomed")

  unk <- unfreq_pop()
  unk$evaluation <- marker_design_checks(pop_assay(NULL))
  expect_equal(evaluation_verdict(unk), "unknown")

  poly <- poly_pop()
  poly$evaluation <- marker_design_checks(pop_assay(NULL))
  expect_equal(evaluation_verdict(poly), "usable")
})

test_that("a marker with no population carries neither rows nor a doom", {
  m <- pop_assay(NULL)
  ch <- marker_design_checks(m)
  expect_length(pop_rows(ch)$criterion %||% character(0), 0L)
  expect_equal(evaluation_verdict(m), "usable")
})
