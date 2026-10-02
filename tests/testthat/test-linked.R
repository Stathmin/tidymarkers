# Linked variants: a co-segregating substitution is baked into the ALT primer.

test_that("a linked substitution raises ALT discrimination and leaves REF alone", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  withr::with_seed(3L, w <- paste(sample(c("A", "C", "G", "T"), 200L, TRUE), collapse = ""))
  W <- Biostrings::DNAString(w)
  v <- 100L
  ref <- substr(w, v, v)
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  # a linked SNP four bases 5' of the target, carrying the complement of the ref base
  linked <- tibble::tibble(
    offset = v - 4L,
    alt = as.character(Biostrings::complement(Biostrings::DNAString(substr(w, v - 4L, v - 4L))))
  )
  base <- design_as_primers(W, v, ref = ref, alt = alt, lengths = 20:22, stab_ats = NA, strands = "+")
  both <- design_as_primers(W, v, ref = ref, alt = alt, lengths = 20:22, stab_ats = NA, strands = "+", linked = linked)
  ob <- base[order(base$as_length), ]
  nb <- both[order(both$as_length), ]
  expect_equal(ob$ref_primer, nb$ref_primer) # REF primer untouched
  expect_true(all(ob$alt_primer != nb$alt_primer)) # ALT carries the linked base
  expect_true(all(nb$discrimination >= ob$discrimination))
  # the linked base sits at primer position (L - 4) from the 5' end, i.e. 4 nt from the 3' end
  expect_equal(substr(nb$alt_primer[nb$as_length == 22L], 22L - 4L, 22L - 4L), linked$alt)
})
