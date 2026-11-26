# =============================================================================
#  Criterion offsets  --  where a named source and the package disagree
# =============================================================================
#  `inst/extdata/criterion_offsets.csv` is the honest comparison behind every
#  `named_source` criterion: the bar the paper states against the bar the
#  registry renders.  It exists because "the package uses the source's limit"
#  was true of some rows and false of others, and the difference was invisible.
#  These tests keep the file a measurement, not a memory: the package side is
#  read back from `criterion_spec()`, the source side from `access_detail`, and
#  the relation is recomputed from both rather than trusted.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

offsets <- function() ext("criterion_offsets.csv")

numbers_in <- function(txt) {
  regmatches(txt, gregexpr("-?[0-9]+([.][0-9]+)?", txt))[[1L]]
}

# A rendered bar read as the allowed interval it denotes.  A strict `< x` is
# read as a closed upper end: no row pairs `< 6` with `<= 5`, so the closed/open
# distinction never decides a relation in this table.
parse_bar <- function(x) {
  x <- trimws(x)
  if (grepl("^>=", x)) {
    return(c(as.numeric(sub("^>= *", "", x)), Inf))
  }
  if (grepl("^>", x)) {
    return(c(as.numeric(sub("^> *", "", x)), Inf))
  }
  if (grepl("^<=", x)) {
    return(c(-Inf, as.numeric(sub("^<= *", "", x))))
  }
  if (grepl("^<", x)) {
    return(c(-Inf, as.numeric(sub("^< *", "", x))))
  }
  if (grepl("^\\[", x)) {
    v <- as.numeric(strsplit(gsub("[][ ]", "", x), ",")[[1L]])
    return(c(v[1L], v[2L]))
  }
  v <- as.numeric(x)
  c(v, v)
}

# `a` is the source's allowed set, `b` the package's.
relation_of <- function(a, b) {
  if (identical(a, b)) {
    return("same")
  }
  if (max(a[1L], b[1L]) > min(a[2L], b[2L])) {
    return("disjoint")
  }
  if (a[1L] >= b[1L] && a[2L] <= b[2L]) {
    return("source_stricter")
  }
  if (b[1L] >= a[1L] && b[2L] <= a[2L]) {
    return("package_stricter")
  }
  "partial_overlap"
}

test_that("offsets: every row joins a live criterion to a registered source", {
  o <- offsets()
  expect_gt(nrow(o), 10L)
  expect_true(all(
    c(
      "id", "source_key", "quantity", "source_value",
      "package_value", "relation", "note"
    ) %in% names(o)
  ))
  expect_true(all(o$id %in% criterion_spec()$id))
  expect_true(all(o$source_key %in% ext("scientific_sources.csv")$key))
  expect_true(all(o$relation %in% c(
    "same", "source_stricter", "package_stricter", "disjoint",
    "partial_overlap"
  )))
  # one row per (claim, source): a second row would silently overwrite the first
  expect_false(any(duplicated(paste(o$id, o$source_key))))
  expect_false(any(is.na(o$note) | !nzchar(o$note)))
})

test_that("offsets: package_value is what the registry renders today", {
  o <- offsets()
  s <- criterion_spec()
  m <- match(o$id, s$id)
  for (i in seq_len(nrow(o))) {
    if (s$kind[m[i]] == "numeric") {
      # the bar itself: compared verbatim, so a registry edit cannot leave the
      # table asserting a limit the package no longer enforces
      expect_identical(
        o$package_value[i], criterion_bound_text(s$bounds[[m[i]]]),
        info = o$id[i]
      )
    } else {
      # boolean/categorical rows carry no `bounds`; the prose label is the limit.
      # Signs are stripped before comparing: a range written `40-60` splits into
      # `40` and `-60`, and the minus is a separator, not a sign.
      announced <- sub("^-", "", numbers_in(o$package_value[i]))
      published <- sub("^-", "", numbers_in(s$limit[m[i]]))
      expect_true(
        all(announced %in% published),
        info = paste(
          o$id[i], "announces", paste(announced, collapse = ","),
          "-- registry says", s$limit[m[i]]
        )
      )
    }
  }
})

test_that("offsets: the declared relation follows from the two bars", {
  o <- offsets()
  for (i in seq_len(nrow(o))) {
    got <- relation_of(parse_bar(o$source_value[i]), parse_bar(o$package_value[i]))
    expect_identical(got, o$relation[i], info = o$id[i])
  }
})

test_that("offsets: every source value is traceable to its access_detail", {
  o <- offsets()
  src <- ext("scientific_sources.csv")
  detail <- src$access_detail[match(o$source_key, src$key)]
  expect_false(any(is.na(detail) | !nzchar(detail)))
  for (i in seq_len(nrow(o))) {
    need <- sub("^-", "", numbers_in(o$source_value[i]))
    have <- sub("^-", "", numbers_in(detail[i]))
    missing <- setdiff(need, have)
    expect_true(
      length(missing) == 0L,
      info = paste(o$id[i], "cites", o$source_key[i], "but", paste(missing, collapse = " "), "is absent")
    )
  }
})
