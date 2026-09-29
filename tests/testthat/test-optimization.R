# "Optimisation optimises": the AS grid must return candidates sorted by
# discrimination, with the top one genuinely the best and on-target beating
# off-target; and when there is no variant there is nothing to discriminate.

test_that("design_as_primers ranks candidates and the top one is the best", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))  # 80 bp, SNP at 40 = T
  d <- design_as_primers(W, 40L, ref = "T", alt = "A",
                         lengths = 20:22, stab_ats = c(NA, 2L))
  expect_gt(nrow(d), 1L)
  expect_true(all(diff(d$discrimination) <= 1e-9))        # non-increasing
  expect_equal(d$discrimination[1L], max(d$discrimination))
  expect_gt(d$discrimination[1L], 0)                      # a real allele mismatch discriminates
  expect_gt(d$ref_on[1L], d$ref_off[1L])                  # on-target beats off-target
})

test_that("with no variant (ref == alt) there is nothing to discriminate", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))
  d <- design_as_primers(W, 40L, ref = "T", alt = "T",
                         lengths = 20:21, stab_ats = NA)
  expect_true(all(abs(d$discrimination) < 0.05))
})
