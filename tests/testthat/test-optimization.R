# "Optimisation optimises": the AS grid must return candidates sorted by
# discrimination, with the top one genuinely the best and on-target beating
# off-target; and when there is no variant there is nothing to discriminate.

test_that("design_as_primers ranks candidates and the top one is the best", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = "")) # 80 bp, SNP at 40 = T
  d <- design_as_primers(W, 40L,
    ref = "T", alt = "A",
    lengths = 20:22, stab_ats = c(NA, 2L)
  )
  expect_gt(nrow(d), 1L)
  # the order is the decided position policy first (-3, then the -2 fallback,
  # then no mismatch), then discrimination within each policy class
  rank <- match(d$stab_policy, c("default", "fallback", "none", "other"))
  expect_equal(rank, sort(rank))
  nonincreasing <- function(x) all(diff(x) <= 1e-9)
  expect_true(all(vapply(split(d$discrimination, rank), nonincreasing, logical(1))))
  # the default (-3) is ranked ahead of the fallback wherever both exist; this
  # grid asks for no -3, so every mismatched row is the fallback and the policy
  # has nothing to reorder
  expect_false(any(d$stab_policy == "default"))
  expect_gt(d$discrimination[1L], 0) # a real allele mismatch discriminates
  expect_gt(d$ref_on[1L], d$ref_off[1L]) # on-target beats off-target
})

test_that("with no variant (ref == alt) there is nothing to discriminate", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))
  d <- design_as_primers(W, 40L,
    ref = "T", alt = "T",
    lengths = 20:21, stab_ats = NA
  )
  expect_true(all(abs(d$discrimination) < 0.05))
})
