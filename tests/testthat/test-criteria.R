# The criterion registry is the single vocabulary for selection + acceptance,
# and check membership is family-based (KASP_FAMILY_CHECKS: marker-agnostic,
# KASP-selector *site*, and allele-specific).  These tests iterate over
# marker_schemes() and read each scheme's spec, so a new chemistry is covered with
# no edit here.

crit_mock <- function(scheme = "snp") {
  p <- kasp_primers()
  # the amplicon acceptance window is per-family, so a "clean" mock must sit inside
  # the scheme's own window (a CAPS band is not 90 bp)
  amp <- as.integer(mean(marker_amplicon_window(marker_scheme(scheme)$family)))
  list(
    assay = tibble::tibble(
      amplicon = amp, as_ref = p[[1]], as_alt = p[[2]], com = p[[3]], com_tm = 60,
      as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L, as_n_var = 0L,
      as_dropout = 0, as_foot_dropout = 0, com_dropout = 0
    ),
    thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), hairpin_dg = c(-1, -1, -1)),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3), b = rep(names(p), times = 3),
      dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
    ),
    amplicons = tibble::tibble(primer = c("AS_ref", "AS_alt"), n_products = c(1L, 0L)),
    efficiency = tibble::tibble(discrimination = 0.9)
  )
}

# the check a scheme's family should carry (used by the membership tests below)
family_checks <- function(scheme) {
  fam <- marker_scheme(scheme)$family
  classes <- KASP_FAMILY_CHECKS[[fam]] %||% character(0)
  c(
    KASP_CHECK_PRIMER,
    if ("site" %in% classes) KASP_CHECK_SITE,
    if ("as" %in% classes) KASP_CHECK_AS,
    # the family constants hold criterion ids; the registry labels are their names
    unname(names(MARKER_FAMILY_CHECKS[[fam]])) %||% character(0)
  )
}

test_that("marker_criteria declares gates, tiers and checks with metadata", {
  crit <- marker_criteria()
  # `rule` is the snapshot multiplex-rule class (§U): stated guidance with no
  # continuous limit, so it is neither a gate, a tier, nor a check.
  expect_setequal(unique(crit$role), c("gate", "tier", "check", "rule"))
  expect_true(all(c("stage", "id", "label", "level", "hard", "limit", "target", "order", "schemes") %in% names(crit)))
  # ZW4: the spec carries the evidence rung (the public view does not expose it yet),
  # right after the deprecated `calibrated` column so no existing column moves
  spec <- tidymarkers:::criterion_spec(KASP_LIMITS)
  expect_true(all(c("source", "calibrated", "evidence") %in% names(spec)))
  expect_equal(
    match("evidence", names(spec)), match("calibrated", names(spec)) + 1L
  )
  expect_true(all(crit$stage %in% c(
    "site", "as", "pair", "duplex", "asq", "caps", "arms", "hrm", "assay",
    marker_schemes()$scheme
  )))
  expect_true(any(crit$hard))
  # every criterion names at least one scheme, and all of them are registered
  expect_true(all(lengths(crit$schemes) >= 1L))
  expect_true(all(unlist(crit$schemes) %in% marker_schemes()$scheme))
  # the gates are the constants, and the hard ones are exactly those off the ladder
  pgate <- crit$id[crit$role == "gate" & crit$stage == "pair"]
  expect_setequal(pgate, paste0("pair:gate:", KASP_GATES_PAIR))
  expect_setequal(
    crit$id[crit$role == "gate" & crit$hard & crit$stage %in% c("as", "pair")],
    c(
      sprintf("as:gate:%s", setdiff(KASP_GATES_AS, KASP_RELAX_AS)),
      sprintf("pair:gate:%s", setdiff(KASP_GATES_PAIR, KASP_RELAX_PAIR))
    )
  )
})

test_that("ZW4.6: a criterion names a source only where one states the bar", {
  spec <- tidymarkers:::criterion_spec(KASP_LIMITS)
  map <- tidymarkers:::criterion_sources()
  # every override points at a real criterion, with a non-empty citation
  expect_true(all(names(map) %in% spec$id))
  expect_false(anyNA(map) || any(!nzchar(map)))
  expect_identical(anyDuplicated(names(map)), 0L)
  # ZW4.6 moved these rows off the default; every other row keeps it
  expect_true(all(spec$source[match(names(map), spec$id)] == map))
  # a row off the map either keeps the default, or names a registry key
  unnamed <- unique(spec$source[!spec$id %in% names(map)])
  expect_true(all(
    unnamed == "package assumption" |
      unnamed %in% utils::read.csv(tm_extdata("scientific_sources.csv"))$key
  ))
  # the citation is a KEY of the source registry, not prose a reader cannot look up
  expect_true(all(map %in% utils::read.csv(tm_extdata("scientific_sources.csv"))$key))
})

test_that("ZW4.6: a raised rung always has a printed source", {
  spec <- tidymarkers:::criterion_spec(KASP_LIMITS)
  ev <- tidymarkers:::criterion_evidence_sources()
  # every raised rung names a real criterion and a legal ladder value
  expect_true(all(names(ev) %in% spec$id))
  expect_identical(anyDuplicated(names(ev)), 0L)
  expect_true(all(ev %in% tidymarkers:::CRITERION_EVIDENCE))
  # the rung on the row is the one the map names
  expect_true(all(as.character(spec$evidence[match(names(ev), spec$id)]) == ev))
  # a raised rung is never above the base without a citation: every raised row
  # is also in criterion_sources(), so it names a registry key
  expect_true(all(names(ev) %in% names(tidymarkers:::criterion_sources())))
  # and every raised row is off the default source
  expect_false(any(spec$source[match(names(ev), spec$id)] == "package assumption"))
})

test_that("ZW4.6: a tier splits into a printed target and a package policy", {
  spec <- tidymarkers:::criterion_spec(KASP_LIMITS)
  tiers <- spec[spec$role == "tier", ]
  expect_gt(nrow(tiers), 0L)
  # every tier has a reference point (the citable half), sourced or not
  expect_false(anyNA(tiers$target))
  expect_true(all(nzchar(tiers$target)))
  # the aggregation policy is a package decision, declared once, and no source
  # overlay ever promotes it: `evidence` on a tier row is the target's rung
  expect_identical(
    tidymarkers:::CRITERION_TIER_POLICY_EVIDENCE, "rule_of_thumb"
  )
  expect_identical(tiers$limit, rep(tidymarkers:::CRITERION_TIER_LIMIT, nrow(tiers)))
  # a promoted tier names its target's source; an unsourced tier keeps the
  # package default and the weakest rung
  promoted <- as.character(tiers$evidence) != "rule_of_thumb"
  expect_true(any(promoted))
  expect_true(all(tiers$source[promoted] != "package assumption"))
  expect_true(all(tiers$source[!promoted] == "package assumption"))
  expect_true(all(as.character(tiers$evidence[!promoted]) == "rule_of_thumb"))
  # only a tier row carries a target
  expect_true(all(is.na(spec$target[spec$role != "tier"])))
  # the public view exposes the same split
  expect_true("target" %in% names(marker_criteria()))
})

test_that("ZW4.6: every criterion carries a machine-checkable provenance statement", {
  spec <- tidymarkers:::criterion_spec(KASP_LIMITS)
  prov <- utils::read.csv(tm_extdata("criterion_provenance.csv"), stringsAsFactors = FALSE)
  expect_setequal(names(prov), c("id", "source", "reason_class", "reason"))
  # one row per criterion, and no criterion left out
  expect_identical(anyDuplicated(prov$id), 0L)
  expect_length(prov$id, length(spec$id))
  expect_setequal(prov$id, spec$id)
  # the shipped source is the spec's source, so the two cannot drift
  src <- spec$source[match(prov$id, spec$id)]
  expect_identical(prov$source, src)
  # every named source joins the registry, so offsets/ledger can key on it
  keys <- utils::read.csv(tm_extdata("scientific_sources.csv"))$key
  named <- !is.na(prov$source) & prov$source != "package assumption"
  expect_true(all(prov$source[named] %in% keys))
  allowed <- c(
    "named_source", "conditions_unstated", "simulator_threshold",
    "inherited_geometry", "scheme_constant"
  )
  expect_false(anyNA(prov$reason_class))
  expect_true(all(prov$reason_class %in% allowed))
  expect_true(all(nzchar(prov$reason)))
  # a named row says so; every default row names why no source states the bar
  unnamed <- prov$source == "package assumption"
  expect_true(all(prov$reason_class[!unnamed] == "named_source"))
  expect_true(all(prov$reason_class[unnamed] %in% setdiff(allowed, "named_source")))
  expect_true(any(unnamed))
})

test_that("the check classes partition the registered checks", {
  labels <- marker_criteria()$label[marker_criteria()$role == "check"]
  family_labels <- unique(unlist(lapply(MARKER_FAMILY_CHECKS, names), use.names = FALSE))
  expect_setequal(
    c(KASP_CHECK_PRIMER, KASP_CHECK_SITE, KASP_CHECK_AS, family_labels), labels
  )
  # 9 marker-agnostic + 5 site + 14 allele-specific + 30 family-specific = 58
  # (the 14th is `as_mismatch_type`, Z0.4: the terminal mismatch's Kwok class;
  # the 30th family check is `family_single_base`, §W)
  expect_length(KASP_CHECK_PRIMER, 9L)
  expect_length(KASP_CHECK_SITE, 5L)
  expect_length(KASP_CHECK_AS, 14L)
  expect_length(family_labels, 30L)
  expect_length(labels, 58L)
  # distinct classes: a base check is never also family-specific
  expect_length(intersect(c(KASP_CHECK_PRIMER, KASP_CHECK_SITE, KASP_CHECK_AS), family_labels), 0L)
  expect_length(intersect(KASP_CHECK_PRIMER, KASP_CHECK_SITE), 0L)
  expect_length(intersect(KASP_CHECK_PRIMER, KASP_CHECK_AS), 0L)
  expect_length(intersect(KASP_CHECK_SITE, KASP_CHECK_AS), 0L)
})

test_that("every scheme carries exactly its family's check classes", {
  for (sc in marker_schemes()$scheme) {
    chk <- marker_criteria(scheme = sc)
    chk <- chk$label[chk$role == "check"]
    # C2: a scheme carries only the checks its chemistry obeys
    expect_true(all(chk %in% family_checks(sc)), info = sc)
    # ... and every marker-agnostic check is there (it holds for any primer pair)
    expect_true(all(KASP_CHECK_PRIMER %in% chk), info = sc)
    # the AS-only and site checks track the family, never the scheme
    info <- marker_scheme(sc)
    as_family <- "as" %in% (KASP_FAMILY_CHECKS[[info$family]] %||% character(0))
    expect_equal(any(KASP_CHECK_AS %in% chk), as_family, info = sc)
  }
  # the multiplex rules are duplex-only
  expect_true(all(paste0("duplex:gate:", KASP_GATES_DUPLEX) %in% marker_criteria(scheme = "duplex")$id))
  expect_false(any(grepl("duplex", marker_criteria(scheme = "snp")$id)))
})

test_that("the AS/COM oligo checks reach only families with an AS primer (C2)", {
  # an AS/COM Tm or GC window is meaningless where those oligos do not exist, so
  # the moved labels follow the `as` class and never leak to a flanking pair.
  as_com <- c("AS Tm (C)", "COM Tm (C)", "AS GC", "COM GC")
  for (sc in marker_schemes()$scheme) {
    crit <- marker_criteria(scheme = sc)
    chk <- crit$label[crit$role == "check"]
    classes <- KASP_FAMILY_CHECKS[[marker_scheme(sc)$family]] %||% character(0)
    expect_equal(all(as_com %in% chk), "as" %in% classes, info = sc)
  }
})

test_that("marker_design_checks returns exactly the registry's checks, in order", {
  for (sc in marker_schemes()$scheme) {
    crit <- marker_criteria(scheme = sc)
    checks <- crit[crit$role == "check", ]
    ch <- marker_design_checks(crit_mock(sc), scheme = sc)
    expect_equal(ch$criterion, checks$label, info = sc)
    expect_equal(ch$limit, checks$limit, info = sc)
    # a KASP mock cannot supply CAPS/LAMP/RPA/... fields; those checks are
    # unknown/not_applicable and must never be reported as a failure.  Every
    # check that is actually evaluable must pass, and `pass` is exactly the
    # logical view of `status`.
    expect_false(any(ch$status == "fail"), info = sc)
    expect_true(all(ch$pass[ch$status == "pass"]), info = sc)
    expect_equal(ch$pass, ch$status == "pass", info = sc)
  }
})

test_that("two in-silico products fail the product-count check", {
  mock <- crit_mock("primer")
  mock$amplicons$n_products <- c(2L, 0L)
  ch <- marker_design_checks(mock, scheme = "primer")
  row <- ch[ch$criterion == "in-silico PCR products (n)", ]
  expect_equal(row$observed, 2)
  expect_equal(as.character(row$status), "fail")
  expect_equal(row$limit, "<= 1")
})

test_that("a criterion's threshold and its rendered limit cannot diverge (Z2)", {
  spec <- tidymarkers:::criterion_spec()
  crit <- marker_criteria()
  # the public view is a projection of the spec, row for row
  expect_equal(nrow(crit), nrow(spec))
  expect_equal(crit$id, spec$id)
  expect_equal(crit$limit, spec$limit)

  numeric_bar <- vapply(
    spec$bounds, function(b) is.numeric(b) && length(b) >= 2L, logical(1)
  )
  # a KASP bar is rendered from its own bound, by the one formatter
  kasp <- numeric_bar & spec$scope == "kasp"
  expect_gt(sum(kasp), 0L)
  expect_equal(
    crit$limit[kasp],
    vapply(spec$bounds[kasp], tidymarkers:::criterion_bound_text, character(1))
  )
  # a family bar keeps its published prose, but its numbers must be in that prose
  fam <- numeric_bar & spec$scope == "family"
  expect_gt(sum(fam), 0L)
  for (i in which(fam)) {
    for (v in spec$bounds[[i]][is.finite(spec$bounds[[i]])]) {
      expect_true(
        grepl(as.character(v), crit$limit[[i]], fixed = TRUE),
        info = paste(spec$id[[i]], "bound", v, "absent from", crit$limit[[i]])
      )
    }
  }
})

test_that("an open-bound family bar renders its finite end, never `NA`", {
  spec <- tidymarkers:::criterion_spec()
  open <- Filter(
    function(b) is.numeric(b) && length(b) >= 2L && any(is.infinite(b)),
    spec$bounds
  )
  expect_gt(length(open), 0L)
  txt <- vapply(open, tidymarkers:::criterion_bound_text, character(1))
  expect_false(anyNA(txt))
  expect_false(any(grepl("NA", txt, fixed = TRUE)))
})

test_that("the scheme filter moves the amplicon bound and its limit together (Z2)", {
  for (sc in c("caps", "lamp", "ssr", "taqman")) {
    win <- tidymarkers:::marker_amplicon_window(marker_scheme(sc)$family)
    lim <- KASP_LIMITS
    lim$amplicon <- win
    spec <- tidymarkers:::criterion_spec(lim)
    i <- match("amplicon", spec$id)
    # one object: the bound the selector grades against and the text the results
    # show are the same numbers
    expect_equal(spec$bounds[[i]], win, info = sc)
    crit <- marker_criteria(scheme = sc)
    expect_equal(
      crit$limit[match("amplicon length (bp)", crit$label)],
      spec$limit[[i]],
      info = sc
    )
  }
})

test_that("the criterion spec keeps its typed class and every documented column (ZA6)", {
  spec <- tidymarkers:::criterion_spec()
  expect_s3_class(spec, "marker_criteria")
  # the constructor must not drop a documented column (it did: `family`)
  expect_true(all(c("family", "bounds", "evidence", "schemes") %in% names(spec)))

  # family check rows name their family; the marker-agnostic primer/site
  # checks are also scoped `family` (they apply to every family) and leave the
  # column `NA` - absent only if the constructor dropped it
  fam <- spec$scope == "family"
  expect_gt(sum(fam), 0L)
  named <- grepl("^family_", spec$id)
  expect_true(all(named <= fam))
  expect_false(anyNA(spec$family[named]))
  expect_true(all(is.na(spec$family[!named])))

  # subsetting restores the class through vctrs, for both row and column cuts
  expect_s3_class(spec[1:3, ], "marker_criteria")
  expect_s3_class(spec[spec$role == "gate", ], "marker_criteria")
  expect_s3_class(spec[, c("id", "label")], "marker_criteria")

  # two classed specs bind to a classed spec, in either order
  expect_s3_class(dplyr::bind_rows(spec, spec), "marker_criteria")
  expect_s3_class(dplyr::bind_rows(spec[0, ], spec), "marker_criteria")
  expect_s3_class(vctrs::vec_c(spec, spec), "marker_criteria")

  # the vctrs ptype is the contract: a plain tibble on either side downgrades
  # rather than half-upgrading (dplyr, not vctrs, keeps the class when the typed
  # frame comes first, so the contract is asserted where it is defined)
  plain <- tibble::as_tibble(spec)
  expect_false(inherits(vctrs::vec_rbind(spec, plain), "marker_criteria"))
  expect_false(inherits(vctrs::vec_rbind(plain, spec), "marker_criteria"))
  expect_false(inherits(vctrs::vec_c(spec, plain), "marker_criteria"))
  expect_false(inherits(vctrs::vec_c(plain, spec), "marker_criteria"))

  # base rbind() is not vctrs: it keeps the class of its first argument, so the
  # mixed order is a named limit (the `vec_c()`/`vec_rbind()` paths above are
  # the contract; a `rbind.marker_criteria` method would add S3 surface for no
  # caller)
  expect_s3_class(rbind(spec, spec), "marker_criteria")
  expect_s3_class(rbind(spec, plain), "marker_criteria")
})

test_that("the criterion constructor refuses a table missing a documented column", {
  spec <- tidymarkers:::criterion_spec()
  bad <- spec[, setdiff(names(spec), "family")]
  expect_error(tidymarkers:::new_marker_criteria(bad), "missing column")
  expect_error(tidymarkers:::new_marker_criteria(spec[, c("id", "label")]), "missing column")
})
