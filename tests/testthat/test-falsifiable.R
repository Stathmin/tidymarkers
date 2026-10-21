# ZW1: every criterion that the registry can publish must be able to FAIL, and
# one that cannot be driven to `fail` in this suite must say so, with a reason,
# in one place.  The ledger (and the rule that checks it) lives in
# `helper-criterion-fixtures.R`, which `dev/check-criteria.R` sources too, so a
# new criterion id fails the pull-request gate as well as this test.
#
# The checks below walk the REGISTRY, not a hand-kept list, so a new criterion
# appears here with no edit to this file: if it has no proof, it fails by name.

# The live KASP acceptance table, built the way `test-kasp-checks.R` builds its
# `good_assay()`: a clean in-memory assay, then one mutation per criterion.
zw1_good_assay <- function() {
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
    ideogram = list(hits = tibble::tibble(primer = c("AS_ref", "AS_alt", "COM"))),
    efficiency = tibble::tibble(discrimination = 0.9)
  )
}

# One driver per site / assay criterion id: a mutation of the clean table that
# makes exactly the named criterion outright `fail` (not `unknown`).
zw1_check_drivers <- function() {
  list(
    site_class = function(a) {
      a$site <- marker_site_class("A", c("AT"))
      a
    },
    site_alt_alleles = function(a) {
      a$site <- list(n_alt = 0L, class = "biallelic", disc = TRUE)
      a
    },
    site_linked_skipped = function(a) {
      a$selection <- list(linked_skipped = 3L)
      a
    },
    site_alt_covered = function(a) {
      a$site <- list(n_alt = 2L, alt_covered = 1L, class = "biallelic", disc = TRUE)
      a
    },
    site_linked_dependent = function(a) {
      a$selection <- list(linked_dependent = TRUE)
      a
    },
    amplicon = function(a) {
      a$assay$amplicon <- 150L
      a
    },
    as_tm = function(a) {
      a$thermo$tm <- c(48, 48, 60)
      a
    },
    com_tm = function(a) {
      a$assay$com_tm <- 50
      a
    },
    com_tm_ge_as_tm = function(a) {
      a$thermo$tm <- c(62, 62, 60)
      a
    },
    as_gc = function(a) {
      a$assay$as_gc <- 0.2
      a
    },
    com_gc = function(a) {
      a$assay$com_gc <- 0.2
      a
    },
    max_run = function(a) {
      a$assay$as_ref <- "AAAAAAAAAAAAAAAA"
      a
    },
    ssr = function(a) {
      a$assay$as_ref <- "ATATATATATATATAT"
      a
    },
    inverted_repeat = function(a) {
      a$assay$com <- "ACGTACGTACGTACGT"
      a
    },
    hairpin_dg = function(a) {
      a$thermo$hairpin_dg <- c(-1, -1, -6)
      a
    },
    cross_dimer_dg = function(a) {
      a$complementarity$dg[3] <- -15
      a
    },
    homodimer_dg = function(a) {
      a$complementarity$dg[1] <- -15
      a
    },
    as_n_3p = function(a) {
      a$assay$as_n_3p <- 2L
      a$assay$as_dropout <- NA_real_
      a
    },
    as_n_var = function(a) {
      a$assay$as_n_var <- 2L
      a$assay$as_foot_dropout <- NA_real_
      a
    },
    com_n_var = function(a) {
      a$assay$n_var <- 2L
      a$assay$com_dropout <- NA_real_
      a
    },
    as_dropout = function(a) {
      a$assay$as_dropout <- 0.2
      a
    },
    as_foot_dropout = function(a) {
      a$assay$as_foot_dropout <- 0.2
      a
    },
    com_dropout = function(a) {
      a$assay$com_dropout <- 0.2
      a
    },
    discrimination = function(a) {
      a$efficiency$discrimination <- 0.1
      a
    },
    in_silico_pcr = function(a) {
      a$amplicons$n_products <- c(1L, 1L)
      a
    },
    offtarget_hits = function(a) {
      a$ideogram$hits <- tibble::tibble(primer = c("AS_ref", "AS_ref"))
      a
    }
  )
}

# A mock of the class a family check is graded on, with the named fields
# replaced.  `base[names(extra)] <- extra` matters: appending would let the
# base field shadow the override, since `$` reads the first match.
zw1_family_mock <- function(class, extra = list()) {
  p <- kasp_primers()
  base <- list(
    snp = tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "G"),
    amplicon = tibble::tibble(start = 200L, end = 400L),
    primers = list(
      forward = p[[1]], reverse = p[[3]], start = 200L, end = 400L,
      as_discrimination = 1
    ),
    ideogram = list(hits = tibble::tibble(primer = c("forward", "reverse")))
  )
  base[names(extra)] <- extra
  structure(base, class = class)
}

# One driver per family check id: the scheme it is graded under and a mock whose
# defining field is wrong, so the generic evaluator reports `fail`.
zw1_family_drivers <- function() {
  list(
    family_caps_site_defined = list("caps", list()),
    family_caps_site_separates = list("caps", list(
      fragments_ref = c(100L, 101L), fragments_alt = c(100L, 101L)
    )),
    family_caps_band_gap = list("caps", list(gap = 5L)),
    family_caps_digest_spans = list("caps", list(
      fragments_ref = c(50L), fragments_alt = c(60L)
    )),
    family_caps_site_conserved = list("caps", list(
      site = list(start = 290L, end = 310L), haplotypes = tibble::tibble(pos = 295L)
    )),
    family_lamp_fip_alt = list("lamp", list(
      primers = list(F3 = "GGGG", B3 = "CCCC", FIP = "ACGTACGTAC", BIP = "TTTT"),
      tm = c(F3 = 57, B3 = 58, FIP = 60, BIP = 60)
    )),
    family_lamp_inner_tm_balance = list("lamp", list(
      primers = list(F3 = "GGGG", B3 = "CCCC", FIP = "ACGTACGTAG", BIP = "TTTT"),
      tm = c(F3 = 57, B3 = 58, FIP = 60, BIP = 70)
    )),
    family_lamp_inner_inside = list("lamp", list(
      primers = list(F3 = "GGGG", B3 = "CCCC", FIP = "ACGTACGTAG", BIP = "TTTT"),
      lamp_arms = list(list(arm5 = 100L, arm3 = 400L))
    )),
    family_lamp_fip_variant = list("lamp", list(
      primers = list(F3 = "GGGG", B3 = "CCCC", FIP = "ACGTACGTAC", BIP = "TTTT"),
      tm = c(F3 = 57, B3 = 58, FIP = 60, BIP = 60)
    )),
    family_rpa_primer_length = list("rpa", list()),
    family_rpa_thf_off = list("rpa", list(probe = list(start = 295L, thf_index = 6L))),
    family_rpa_probe_spans = list("rpa", list(probe = list(start = 250L, end = 260L))),
    family_rpa_variant_under = list("rpa", list(
      probe = list(start = 250L, end = 320L),
      primers = list(
        forward = kasp_primers()[["AS_ref"]], reverse = kasp_primers()[["COM"]],
        start = 200L, end = 400L, variant_under = c(alt = 1L, ref = 1L, other = 0L)
      )
    )),
    family_crispr_pam = list("crispr", list(gates = list(gate_crispr_pam = FALSE))),
    family_crispr_seed = list("crispr", list(gates = list(gate_crispr_seed = FALSE))),
    family_crispr_guide_spans = list("crispr", list(guide_start = 100L, guide_end = 200L)),
    family_crispr_conserved = list("crispr", list(
      guide_start = 290L, guide_end = 320L, haplotypes = tibble::tibble(pos = 310L)
    )),
    family_crispr_offtarget = list("crispr", list(off_target = 3L)),
    family_probe_spans = list("taqman", list(
      probes = tibble::tibble(start = 150L, end = 160L)
    )),
    family_probe_tm_delta = list("taqman", list(
      probes = tibble::tibble(start = 295L, end = 310L, tm = 40)
    )),
    family_probe_asp_discrimination = list("taqman", list(
      primers = list(
        forward = kasp_primers()[["AS_ref"]], reverse = kasp_primers()[["COM"]],
        start = 200L, end = 400L, as_discrimination = 0
      ),
      probes = tibble::tibble(start = 295L, end = 310L, tm = 65)
    )),
    family_probe_extension = list("taqman", list(
      probes = tibble::tibble(start = 295L, end = 310L, tm = 65),
      extension = list(start = 100L, end = 200L)
    )),
    family_sanger_read_distance = list("sanger", list(
      amplicon = tibble::tibble(start = 290L, end = 500L),
      primers = list(
        forward = kasp_primers()[["AS_ref"]], reverse = kasp_primers()[["COM"]],
        start = 290L, end = 500L
      )
    )),
    family_sanger_primers_inside = list("sanger", list(
      seq_primers = list(
        forward = "ACGT", reverse = "ACGT", forward_start = 150L, reverse_end = 400L
      )
    )),
    family_sanger_pair_defined = list("sanger", list(
      seq_primers = list(forward = "ACGT")
    )),
    family_size_repeat_unit = list("ssr", list(
      "repeat" = tibble::tibble(start = 300L, end = 301L, unit = 1L, copies = 10L)
    )),
    family_size_repeat_copies = list("ssr", list(
      "repeat" = tibble::tibble(start = 300L, end = 303L, unit = 2L, copies = 2L)
    )),
    family_size_delta = list("ssr", list(
      "repeat" = tibble::tibble(start = 300L, end = 301L, unit = 2L, copies = 1L)
    )),
    family_size_fragment = list("ssr", list(
      amplicon = tibble::tibble(start = 1L, end = 5000L)
    ))
  )
}

# The status of `criterion` (a registry label) in a checks table.
zw1_status_of <- function(ch, criterion) {
  i <- match(criterion, ch$criterion)
  if (is.na(i)) NA_character_ else as.character(ch$status[i])
}

test_that("every registered criterion is exempt with a reason or has a fail fixture", {
  ledger <- zw1_ledger()
  coverage <- zw1_ledger_coverage(ledger)
  # the registry is walked, not a hand-kept list: gates and checks both count
  expect_gt(length(ledger$registered), 0L)
  expect_equal(nrow(ledger$fixtures) + length(ledger$exempt) + length(coverage$uncovered),
    length(ledger$registered),
    info = "every registered gate/check id is exactly one of fixture-backed, exempt or uncovered"
  )
  expect_length(intersect(ledger$fixtures$id, names(ledger$exempt)), 0L)
  # an exemption is a named reason, never a bare key
  expect_true(all(nzchar(trimws(ledger$exempt))))
  # a stale row (an id the registry no longer has) is a lie, not evidence
  expect_length(coverage$stale_exemptions, 0L)
  expect_length(coverage$stale_fixtures, 0L)
  # nothing may be left without proof; report the ids and labels, never a count
  expect_equal(length(coverage$uncovered), 0L,
    info = paste(zw1_uncovered_problems(ledger), collapse = "\n")
  )
})

test_that("every cited fail fixture names a real test in this repo", {
  citations <- unique(c(zw1_fail_fixtures()$test, zw1_status_citations()$test))
  names(citations) <- citations
  exists <- vapply(citations, zw1_citation_exists, logical(1))
  expect_true(all(exists),
    info = paste("a fixture cites a test that is not there:", paste(citations[!exists], collapse = "; "))
  )
  # ... and the path is a test file under tests/testthat, not a stray string
  expect_true(all(grepl("^tests/testthat/test-.+[.]R::.+$", citations)))
})

test_that("a new criterion id without a fixture fails the gate, not just the test", {
  # the pull-request gate applies the same rule to the same table
  expect_length(zw1_ledger_violations(), 0L)
  # a criterion the registry publishes but the ledger does not cover is reported
  # by id, with its label and the exact remedy - so the gate names what to add
  # instead of a bare count
  ledger <- zw1_ledger()
  published <- c(ledger$registered, "gate_zw1_fake")
  problems <- zw1_uncovered_problems(ledger, published)
  expect_length(problems, 1L)
  expect_match(problems, "gate_zw1_fake", fixed = TRUE)
  expect_match(problems, "add a fail fixture to tests/testthat/test-falsifiable.R or an exemption with a reason",
    fixed = TRUE
  )
  expect_true(any(grepl("gate_zw1_fake", zw1_ledger_violations(ledger, published), fixed = TRUE)))
  # a fixture row for an id the registry dropped is stale and also fails
  stale <- ledger
  stale$fixtures <- rbind(
    stale$fixtures[0L, ],
    tibble::tibble(id = "gate_zw1_fake", test = ledger$fixtures$test[1L])
  )
  expect_true(any(grepl("stale fixture", zw1_ledger_problems(stale))))
  expect_match(zw1_ledger_problems(stale), "gate_zw1_fake", fixed = TRUE)
  # an exemption for an id the registry dropped is stale too
  stale_exempt <- ledger
  stale_exempt$exempt <- stats::setNames("needs a genome this suite does not carry", "gate_zw1_fake")
  expect_match(zw1_ledger_problems(stale_exempt), "stale exemption", fixed = TRUE)
  expect_match(zw1_ledger_problems(stale_exempt), "gate_zw1_fake", fixed = TRUE)
  # an exemption with an empty reason is not an exemption
  blank <- ledger
  blank$exempt <- stats::setNames("", "gate_zw1_fake")
  expect_true(any(grepl("empty reason", zw1_ledger_problems(blank))))
})

test_that("coverage is derived from the ledger, so a removed row is a violation", {
  # a cached coverage answer would let a deleted row pass silently: this is the
  # hole the gate exists to close, so both removals are checked by id
  ledger <- zw1_ledger()
  expect_length(zw1_ledger_violations(ledger), 0L)

  dropped <- ledger$fixtures$id[1L]
  without_fixture <- ledger
  without_fixture$fixtures <- ledger$fixtures[-1L, ]
  problems <- zw1_ledger_violations(without_fixture)
  expect_length(problems, 1L)
  expect_match(problems, paste0("criterion ", dropped, " ("), fixed = TRUE)
  expect_match(problems, "has no fail fixture", fixed = TRUE)

  # the same hole via an exemption: an id covered only by an exemption is clean
  # until that exemption is dropped, and then reported by name
  target <- ledger$fixtures$id[1L]
  exempted <- ledger
  exempted$fixtures <- ledger$fixtures[-1L, ]
  exempted$exempt <- stats::setNames("needs a genome this suite does not carry", target)
  expect_length(zw1_ledger_violations(exempted), 0L)

  exempted$exempt <- stats::setNames(character(0), character(0))
  problems <- zw1_ledger_violations(exempted)
  expect_length(problems, 1L)
  expect_match(problems, paste0("criterion ", target, " ("), fixed = TRUE)
})

test_that("each site / assay criterion really fails when its driver is applied", {
  drivers <- zw1_check_drivers()
  labels <- stats::setNames(
    criterion_spec(KASP_LIMITS)$label,
    criterion_spec(KASP_LIMITS)$id
  )
  # the clean table passes every check: a driver's effect is the mutation, not
  # a fixture that was already failing
  clean <- marker_design_checks(zw1_good_assay())
  expect_false(any(clean$status == "fail"))
  expect_false(any(clean$status == "unknown"))
  for (id in names(drivers)) {
    ch <- marker_design_checks(drivers[[id]](zw1_good_assay()))
    expect_equal(zw1_status_of(ch, labels[[id]]), "fail", info = id)
  }
})

test_that("each family check really fails when its defining field is wrong", {
  drivers <- zw1_family_drivers()
  labels <- stats::setNames(
    criterion_spec(KASP_LIMITS)$label,
    criterion_spec(KASP_LIMITS)$id
  )
  for (id in names(drivers)) {
    d <- drivers[[id]]
    ch <- marker_design_checks(zw1_family_mock(paste0("marker_", d[[1L]]), d[[2L]]), scheme = d[[1L]])
    expect_equal(zw1_status_of(ch, labels[[id]]), "fail", info = id)
  }
})

test_that("each gate reaches ok, relaxed and violated through its own DAG", {
  gates <- criterion_spec(KASP_LIMITS)
  gates <- gates$id[gates$role == "gate"]
  schemes <- marker_schemes()$scheme
  names_of <- lapply(schemes, function(s) {
    tibble::as_tibble(marker_decision_graph(list(scheme = s)))$name
  })
  names(names_of) <- schemes
  for (node in gates) {
    # the KASP core gates are shared by every KASP-family scheme; a scheme's own
    # gate belongs to that scheme alone.  Either way it is drawn on a real DAG.
    owners <- names(names_of)[vapply(names_of, function(n) node %in% n, logical(1))]
    expect_gt(length(owners), 0L)
    scheme <- owners[[1L]]
    seg <- sub(":gate:.*$", "", node)
    bare <- sub("^.*:gate:", "", node)
    # no record -> ok, a relaxed record -> relaxed, a violated record -> violated
    status <- function(sel) {
      g <- tibble::as_tibble(marker_decision_graph(list(scheme = scheme, selection = sel)))
      as.character(g$status[g$name == node])[1L]
    }
    expect_equal(status(list()), "ok", info = node)
    expect_equal(status(list(relaxed = bare)), "relaxed", info = node)
    expect_equal(status(list(violated = bare)), "violated", info = node)
    # a violated gate is not decorative: the DAG routes its stage somewhere,
    # either to the shared DOOMED terminal or to the ladder's refinement branch
    edges <- decision_edges(marker_decision_graph(list(scheme = scheme)))
    doomed <- any(edges$from == node & edges$to == "doomed")
    refine <- any(edges$from == paste0(seg, ":decision") & edges$kind == "fail")
    expect_true(doomed || refine, info = node)
  }
})

test_that("unknown and not_applicable are reachable, and neither is a pass", {
  cites <- zw1_status_citations()
  expect_setequal(cites$status, c("unknown", "not_applicable"))
  # unknown: missing discrimination evidence is not a success
  a <- zw1_good_assay()
  a$efficiency$discrimination <- NA_real_
  ch <- marker_design_checks(a)
  expect_equal(zw1_status_of(ch, "discrimination"), "unknown")
  expect_equal(evaluation_verdict(ch), "unknown")
  # not_applicable: a family check on another marker class neither passes nor dooms
  x <- zw1_family_mock("marker_caps")
  lm <- marker_design_checks(x, scheme = "lamp")
  fam <- lm$criterion %in% unname(names(MARKER_CHECK_LAMP))
  expect_true(all(lm$status[fam] == "not_applicable"))
  expect_false(any(lm$pass[fam]))
})
