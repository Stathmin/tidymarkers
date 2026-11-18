# The one-pot matrix is data: these tests assert that the optimiser reads it,
# that a rejected pair never shares a group, and that every group re-passes the
# same matrix object it was built with.

one_pot_fixture <- function() {
  data.frame(
    label = c("m1", "m2", "m3", "m4", "m5"),
    scheme = c("kasp", "kasp", "hrm", "caps", "kasp"),
    tm = c(60, 61, 60.5, 68, 60.2),
    amplicon = c(100, 140, 101, 80, 200),
    channel = c("FAM", "HEX", "FAM", "HEX", NA)
  )
}

one_pot_same_group <- function(groups, a, b) {
  ga <- groups$group[groups$label == a]
  gb <- groups$group[groups$label == b]
  length(ga) == 1L && identical(ga, gb)
}

test_that("the matrix ships exactly the three declared rules as data", {
  m <- one_pot_matrix()
  expect_s3_class(m, "tbl_df")
  expect_named(m, c("rule", "field", "comparison", "value", "source"))
  expect_equal(m$rule, c("tm_span", "amplicon_gap", "channel_unique"))
  expect_equal(m$field, c("tm", "amplicon", "channel"))
  expect_equal(m$comparison, c("range_max", "min_gap", "unique"))
  expect_equal(m$value, c(3, 10, NA_real_))
  # honest provenance: no citation is invented for a locally chosen threshold
  expect_true(all(nzchar(m$source)))
  expect_true(all(m$source == "local policy"))
})

test_that("one_pot_compatible answers per rule and per set", {
  assays <- one_pot_fixture()
  expect_true(one_pot_compatible(assays[1:2, , drop = FALSE]))
  expect_false(one_pot_compatible(assays[c(1, 3), , drop = FALSE]))
  expect_true(one_pot_compatible(assays[3:2, , drop = FALSE]))

  verdict <- one_pot_compatible(assays[c(1, 3), , drop = FALSE], require = "any")
  expect_type(verdict, "logical")
  expect_named(verdict, one_pot_matrix()$rule)
  expect_true(verdict[["tm_span"]])
  expect_false(verdict[["amplicon_gap"]]) # 1 bp apart on the same channel
})

test_that("a rejected pair never shares a group, and every group re-passes", {
  assays <- one_pot_fixture()
  m <- one_pot_matrix()
  res <- marker_one_pot(assays, matrix = m)
  pairs <- utils::combn(nrow(assays), 2L, simplify = FALSE)

  for (p in pairs) {
    pair <- assays[p, , drop = FALSE]
    if (!one_pot_compatible(pair, m)) {
      expect_false(
        one_pot_same_group(res$groups, pair$label[1L], pair$label[2L]),
        info = paste("rejected pair shares a group:", paste(pair$label, collapse = " + "))
      )
    }
  }

  # every group re-passes the same matrix object the grouping was built with
  expect_identical(res$matrix, m)
  for (g in seq_len(res$n_groups)) {
    members <- res$groups[res$groups$group == g, , drop = FALSE]
    expect_true(one_pot_compatible(members, m), info = paste("group", g))
  }
  # nothing is silently dropped: every assay is either placed or reported
  expect_equal(nrow(res$groups) + nrow(res$unassigned), nrow(assays))
})

test_that("max_markers bounds group size and the trade-off is monotone", {
  assays <- one_pot_fixture()
  m <- one_pot_matrix()
  capped <- marker_one_pot(assays, matrix = m, max_markers = 2L)
  expect_lte(max(tabulate(capped$groups$group)), 2L)

  loose <- marker_one_pot(assays, matrix = m, max_markers = 4L)
  strict <- marker_one_pot(assays, matrix = m, max_markers = 1L)
  expect_lte(loose$n_groups, strict$n_groups)
  expect_equal(strict$n_groups, nrow(strict$groups))
})

test_that("an assay no pot can hold is unassigned with a reason, not dropped", {
  assays <- data.frame(
    label = c("solid1", "solid2", "volatile"),
    scheme = c("kasp", "hrm", "kasp"),
    tm = c(60, 60.5, NA),
    amplicon = c(100, 140, 120),
    channel = c("FAM", "HEX", NA)
  )
  m <- one_pot_matrix()
  res <- marker_one_pot(assays, matrix = m)

  # it clashes with every other assay ...
  for (i in 1:2) {
    expect_false(one_pot_compatible(assays[c(i, 3L), , drop = FALSE], m))
  }
  # ... and is reported, with a reason naming the rule that failed
  expect_equal(res$unassigned$label, "volatile")
  expect_true(all(nzchar(res$unassigned$reason)))
  expect_match(res$unassigned$reason, "tm_span")
  expect_setequal(res$groups$label, c("solid1", "solid2"))
})

test_that("same tm, same channel and a close amplicon separate the assays", {
  m <- one_pot_matrix()
  # both channels unassigned: nothing but the amplicon gap can split them
  unassigned_channels <- data.frame(
    label = c("a", "b", "c"),
    scheme = "kasp",
    tm = c(60, 60, 60),
    amplicon = c(100, 105, 140),
    channel = NA_character_
  )
  res <- marker_one_pot(unassigned_channels, matrix = m)
  expect_false(one_pot_same_group(res$groups, "a", "b"))
  expect_true(one_pot_same_group(res$groups, "a", "c"))
  for (g in seq_len(res$n_groups)) {
    expect_true(one_pot_compatible(res$groups[res$groups$group == g, , drop = FALSE], m))
  }

  # identical assigned channels: only one assay per channel per pot
  one_channel <- within(unassigned_channels, channel <- "FAM")
  res_channel <- marker_one_pot(one_channel, matrix = m)
  expect_equal(res_channel$n_groups, 3L)
  for (g in seq_len(res_channel$n_groups)) {
    expect_true(one_pot_compatible(res_channel$groups[res_channel$groups$group == g, , drop = FALSE], m))
  }
})

test_that("the grouping follows the matrix values, not constants in the grouper", {
  assays <- one_pot_fixture()
  strict <- marker_one_pot(assays)
  relaxed_matrix <- one_pot_matrix()
  relaxed_matrix$value[relaxed_matrix$rule == "tm_span"] <- 100
  relaxed <- marker_one_pot(assays, matrix = relaxed_matrix)

  # m1 (tm 60) and m4 (tm 68) are 8 C apart: only the edited matrix allows them
  expect_false(one_pot_compatible(assays[c(1, 4), , drop = FALSE]))
  expect_true(one_pot_compatible(assays[c(1, 4), , drop = FALSE], relaxed_matrix))
  expect_lt(relaxed$n_groups, strict$n_groups)

  unknown <- one_pot_matrix()
  unknown$comparison[1L] <- "best_effort"
  expect_error(marker_one_pot(assays, matrix = unknown), "comparison")
})

test_that("bad input aborts and says what is wrong", {
  assays <- one_pot_fixture()
  expect_error(one_pot_compatible(assays[, setdiff(names(assays), "channel")]), "channel")
  expect_error(one_pot_compatible(assays[, setdiff(names(assays), c("tm", "amplicon"))]), "tm")
  expect_error(one_pot_compatible(within(assays, tm <- as.character(tm))), "numeric")
  expect_error(one_pot_compatible(assays, require = "most"), "all")
  expect_error(marker_one_pot(assays, max_markers = 0L), "positive")
  expect_error(marker_one_pot(assays, max_markers = 2.5), "whole number")
  expect_error(marker_one_pot(assays, matrix = one_pot_matrix()[0, ]), "no rules")
})

test_that("the result carries the experimental caveat and prints", {
  res <- marker_one_pot(one_pot_fixture())
  expect_s3_class(res, "marker_one_pot")
  expect_true(res$experimental)
  expect_named(res$groups, c("group", "label", "scheme", "tm", "amplicon", "channel"))
  expect_named(res$unassigned, c("label", "reason"))
  # called by name: dispatch needs the `S3method(print, marker_one_pot)` entry
  # that roxygen adds from the `@export` on the method
  expect_output(print.marker_one_pot(res), "<marker_one_pot>")
  expect_output(print.marker_one_pot(res), "experimental")
})
