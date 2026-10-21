# The pairwise primer-pool compatibility report and its tube-level wrapper.

test_that("primer_set_compatibility reports one row per unordered pair", {
  cc <- primer_set_compatibility(c(AS1 = "GAAGGTGACCAAGTTCATGCT", AS2 = "GAAGGTCGGAGTCAACGGATT"))
  expect_true(all(c("a", "b", "kind", "run", "dg", "pass") %in% names(cc)))
  expect_equal(sum(cc$kind == "self"), 2L)
  expect_equal(sum(cc$kind == "cross"), 1L)
  # a primer and its own reverse complement form a full-length cross-dimer
  p <- "ACGTACGTACGTACGTACGT"
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(p)))
  bad <- primer_set_compatibility(c(x = p, y = rc))
  expect_false(bad$pass[bad$kind == "cross"])
})

test_that("tube_compatibility pools oligos, drops duplicates and tolerates modified bases", {
  expect_null(tube_compatibility(character(0)))
  expect_null(tube_compatibility(c(only = "ACGTACGTACGT")))
  expect_null(tube_compatibility(c(a = "ACGTACGTACGT", b = "ACGTACGTACGT"))) # same oligo
  cc <- tube_compatibility(c(
    fwd = "GAAGGTGACCAAGTTCATGCT", rev = "GAAGGTCGGAGTCAACGGATT",
    probe = "FAM-ACGTXACGT" # a modified probe (non-ACGT marks) must not error
  ))
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(cc)))
  expect_equal(nrow(cc), 6L) # 3 self + 3 cross
})

test_that("primer_set_compatibility reports dg at the caller's reaction conditions (Z3)", {
  skip_if_no_oligoarrayaux()
  p <- c(AS1 = "GAAGGTGACCAAGTTCATGCT", AS2 = "GAAGGTCGGAGTCAACGGATT")
  base <- primer_set_compatibility(p)
  # no conditions = the historical 50 mM Na / 37 C, byte for byte
  expect_equal(primer_set_compatibility(p, conditions = reaction_conditions())$dg, base$dg)
  expect_equal(primer_set_compatibility(p, conditions = NULL)$dg, base$dg)
  # one row per (scheme, oligo pair): the reported dg is the scheme-temperature
  # dimer_thermo() value, and every scheme that runs warmer than 37 C differs
  # from the default (a warmer duplex is weaker, i.e. a less negative dG)
  grid <- tidyr::expand_grid(
    scheme = c("snp", "lamp", "rpa", "castpcr", "taqman", "caps"),
    row = seq_len(nrow(base))
  )
  out <- purrr::map(seq_len(nrow(grid)), function(k) {
    scheme <- grid$scheme[k]
    r <- grid$row[k]
    ref <- dimer_thermo(p, conditions = scheme_conditions(scheme))
    same <- (ref$a == base$a[r] & ref$b == base$b[r]) |
      (ref$a == base$b[r] & ref$b == base$a[r])
    got <- primer_set_compatibility(p, conditions = scheme_conditions(scheme))
    tibble::tibble(
      scheme = scheme, row = r, expected = ref$dg[same][[1L]],
      observed = got$dg[r], default = base$dg[r]
    )
  }) |>
    purrr::list_rbind()
  expect_equal(out$observed, out$expected)
  # the `pass` verdict reads the base-pair runs only: temperature never moves it
  for (s in unique(grid$scheme)) {
    expect_equal(primer_set_compatibility(p, conditions = scheme_conditions(s))$pass, base$pass)
  }
  hot <- out[out$scheme %in% c("snp", "lamp", "rpa"), ]
  expect_true(all(hot$observed != hot$default))
  # an ambient scheme (25 C) is *more* stable than the 37 C default
  amb <- out[out$scheme == "caps", ]
  expect_true(all(amb$observed != amb$default))
  expect_true(all(amb$observed < amb$default))
})

test_that("tube_compatibility forwards its reaction conditions", {
  skip_if_no_oligoarrayaux()
  set <- c(fwd = "GAAGGTGACCAAGTTCATGCT", rev = "GAAGGTCGGAGTCAACGGATT")
  hot <- tube_compatibility(set, conditions = scheme_conditions("lamp"))
  cold <- tube_compatibility(set)
  expect_equal(hot$dg, primer_set_compatibility(set, conditions = scheme_conditions("lamp"))$dg)
  expect_true(all(hot$dg > cold$dg)) # 65 C duplexes are weaker than at 37 C
  expect_equal(hot$pass, cold$pass)
})
