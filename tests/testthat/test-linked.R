# Linked variants: a co-segregating substitution is baked into the ALT primer.

test_that("a linked substitution raises ALT discrimination and leaves REF alone", {
  skip_if_no_oligoarrayaux()
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

test_that("a design that holds only on linkage is reported, not silently accepted", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(bseq(500, 2L))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  v <- 170L
  ref <- as.character(Biostrings::subseq(W, v, v))
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  # one linked substitution just inside the AS 3' window / footprint: the clean
  # window is bought by the linkage (the base rides with the ALT haplotype and is
  # baked into the ALT primer), not by the target sequence.
  lpos <- v + 3L
  lref <- as.character(Biostrings::subseq(W, lpos, lpos))
  lalt <- setdiff(c("A", "C", "G", "T"), lref)[1L]
  variants <- tibble::tibble(chrom = "Chr1", pos = lpos, ref = lref, alt = lalt)
  linked <- tibble::tibble(target = v, linked = lpos)
  snp <- tibble::tibble(chrom = "Chr1", pos = v, ref = ref, alt = alt)
  msgs <- character()
  a <- withCallingHandlers(
    assemble_kasp_assay(genome, snp,
      variants = variants, linked = linked,
      flank = 150L, cores = 1L
    ),
    warning = function(cond) {
      msgs <<- c(msgs, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(inherits(a, "kasp_doomed"))
  expect_true(isTRUE(a$selection$linked_dependent))
  expect_true(any(grepl("relies on linkage", msgs, fixed = TRUE)))
  # it is a first-class acceptance criterion, not only a warning
  chk <- kasp_design_checks(a)
  expect_false(chk$pass[chk$criterion == "design relies on linkage"])
})
