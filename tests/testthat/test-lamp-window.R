# What the LAMP time window actually supports (data-raw/proposed/
# lamp_status.decision.md clauses 1 and 6). Data-only, from inst/extdata.
#
# Two mosquito AS-LAMP papers report a measured Tt for their target arms and a
# CENSORED non-target arm ("no amplification within 75 min"): badolo_kdr2012
# (Malar J 2012;11:227, PMC3407793) for the L1014F kdr-w SNP in Anopheles
# gambiae s.l., and ace1r_aslamp2015 for G119S in ace-1. The censor is an
# inequality (Tt > 75), not a measurement, so the only defensible statement is
# a window with an open upper side - not an allele gap, and not a Ct/Tt
# difference between two alleles on one template.

lw_read <- function(name) {
  utils::read.csv(testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
}

lw_data <- function() {
  list(
    badolo = lw_read("lamp_badolo_kdr2012_time_to_positive.csv"),
    ace1 = lw_read("ace1r_aslamp2015_time_to_positive.csv")
  )
}

test_that("each paper's censored arms are non-targets with no point value", {
  d <- lw_data()
  # badolo_kdr2012: 10 arms = two primer sets x (2 plasmids + 3 gDNA genotypes);
  # 4 are censored, and every censored arm is the non-target allele
  b <- d$badolo
  expect_identical(nrow(b), 10L)
  expect_identical(sum(b$is_target %in% TRUE), 6L)
  expect_identical(sum(b$is_target %in% FALSE), 4L)
  cen <- b[!(b$is_target %in% TRUE), ]
  expect_identical(
    sort(paste(cen$template, cen$primer_set)),
    sort(c(
      "gDNA +/+ kdr-w", "gDNA kdr-w/kdr-w wild-type",
      "plasmid kdr-w wild-type", "plasmid wild-type kdr-w"
    ))
  )
  expect_true(all(is.na(cen$time_to_positive_min)))
  expect_true(all(cen$observation_window_min == 75))
  expect_true(all(grepl("not amplified within 75 min", cen$censoring, fixed = TRUE)))
  expect_true(all(grepl("not measured", cen$censoring, fixed = TRUE)))

  # ace1r_aslamp2015: 4 arms = two primer sets x (homozygote target + heterozygote);
  # EVERY arm has a measured Tt, and the non-target genotype is named per arm
  a <- d$ace1
  expect_identical(nrow(a), 4L)
  expect_true(all(is.finite(a$time_to_positive_min)))
  expect_identical(a$template_genotype, c("Ace1SS", "Ace1RS", "Ace1RR", "Ace1RS"))
  expect_identical(a$non_target_genotype, c("Ace1RR", "Ace1RR", "Ace1SS", "Ace1SS"))
  expect_true(all(a$non_target_genotype != a$template_genotype))
  # the paper never saw its non-target amplified before the cut-off
  expect_true(all(a$non_target_amplified_before_cutoff %in% FALSE))
  expect_true(all(a$observation_window_min == 75))
  expect_true(all(grepl("not amplified within 75 min", a$censoring, fixed = TRUE)))
})

test_that("every measured arm is a target arm", {
  d <- lw_data()
  b <- d$badolo
  expect_identical(b$censoring == "measured", b$is_target %in% TRUE)
  expect_true(all(is.finite(b$time_to_positive_min[b$is_target %in% TRUE])))
  expect_true(all(!is.na(b$time_to_positive_min[b$is_target %in% TRUE])))
  # the census ace-1 file has no is_target column: targetness is the template
  # against the arm's named non-target
  expect_false("is_target" %in% names(d$ace1))
  a <- d$ace1
  wt_set <- a$primer_set[a$template_genotype == "Ace1SS"]
  rr_set <- a$primer_set[a$template_genotype == "Ace1RR"]
  expect_identical(wt_set, "wild-type (ace-1S) specific")
  expect_identical(rr_set, "resistance (ace-1R) specific")
})

test_that("measured target Tt is at most the 75-min window", {
  d <- lw_data()
  b <- d$badolo
  bt <- b$time_to_positive_min[b$is_target %in% TRUE]
  # exact values, as printed: plasmid 50 / 52, gDNA 60 / 75 / 65 / 60
  expect_identical(sort(bt), c(50L, 52L, 60L, 60L, 65L, 75L))
  expect_identical(range(bt), c(50L, 75L))
  expect_true(all(bt <= b$observation_window_min[b$is_target %in% TRUE]))

  a <- d$ace1
  at <- a$time_to_positive_min
  expect_identical(at, c(50L, 64L, 51L, 55L))
  expect_identical(range(at), c(50L, 64L))
  expect_true(all(at < a$observation_window_min))
  expect_true(all(a$time_is_approximate %in% TRUE))

  # The upper side is CLOSED, not open: the badolo gDNA heterozygote under the
  # wild-type set sits ON the cut-off (75 min at a 75-min window), so the
  # statement is "<= the window", never "strictly below it".
  het_wt <- b$template == "gDNA kdr-w/+" & b$primer_set == "wild-type"
  expect_identical(b$time_to_positive_min[het_wt], 75L)
  expect_identical(b$time_to_positive_min[het_wt], b$observation_window_min[het_wt])
  expect_identical(sum(bt == b$observation_window_min[b$is_target %in% TRUE]), 1L)
  # Per primer set, in row order: the plasmid (50 wt / 52 kdr-w) plus the two
  # gDNA target plates.  The heterozygote is the slowest detected arm of BOTH
  # sets (75 min wild-type, 65 min kdr-w), so the window is set by the
  # heterozygote - whose non-target arm is censored, never timed.
  wt_arm <- b$primer_set == "wild-type" & b$is_target %in% TRUE
  kdrw_arm <- b$primer_set == "kdr-w" & b$is_target %in% TRUE
  expect_identical(b$time_to_positive_min[wt_arm], c(50L, 60L, 75L))
  expect_identical(b$time_to_positive_min[kdrw_arm], c(52L, 65L, 60L))
  for (sel in list(wt_arm, kdrw_arm)) {
    het <- b$template[sel] == "gDNA kdr-w/+"
    expect_identical(b$time_to_positive_min[sel][het], max(b$time_to_positive_min[sel]))
  }
})

test_that("in both papers the heterozygote is the slowest arm of its own set", {
  a <- lw_data()$ace1
  ps <- a$primer_set
  wt <- ps == "wild-type (ace-1S) specific"
  rr <- ps == "resistance (ace-1R) specific"
  het <- a$template_genotype == "Ace1RS"
  expect_identical(a$time_to_positive_min[wt & het], 64L)
  expect_identical(a$time_to_positive_min[wt & !het], 50L)
  expect_identical(a$time_to_positive_min[rr & het], 55L)
  expect_identical(a$time_to_positive_min[rr & !het], 51L)
  # the heterozygote is the slowest arm of BOTH ace-1 sets, as in badolo ...
  expect_identical(a$time_to_positive_min[wt & het], max(a$time_to_positive_min[wt]))
  expect_identical(a$time_to_positive_min[rr & het], max(a$time_to_positive_min[rr]))
  # ... but the SIZE of the delay is a property of each paper and is never the
  # same: 14 / 4 min here against 15 / 15 min in badolo.  No model may average
  # them, and neither is an allele gap (a gap needs one template under BOTH
  # alleles, which no landed row is).
  expect_identical(
    a$time_to_positive_min[wt & het] - a$time_to_positive_min[wt & !het], 14L
  )
  expect_identical(
    a$time_to_positive_min[rr & het] - a$time_to_positive_min[rr & !het], 4L
  )
})

test_that("the censoring threshold is 75 min in both papers", {
  d <- lw_data()
  expect_identical(unique(d$badolo$observation_window_min), 75L)
  expect_identical(unique(d$ace1$observation_window_min), 75L)
  # equal thresholds are also what the censoring sentences say.  The ace-1
  # file states the window in every row (all four are measured); the badolo
  # file states it only in its censored rows, the others carry the token
  # "measured".  Extract the number only where a number is stated.
  w <- function(x) {
    x <- x[grepl("within [0-9]+ min", x)]
    unique(sub("^.*within ([0-9]+) min.*$", "\\1", x))
  }
  expect_identical(w(d$badolo$censoring), "75")
  expect_identical(w(d$ace1$censoring), "75")
  expect_identical(sum(d$badolo$censoring == "measured"), 6L)
  expect_identical(sum(d$ace1$censoring == "measured"), 0L)
})

test_that("the window is censored on one side, so no allele gap can be claimed", {
  d <- lw_data()
  b <- d$badolo
  cen <- b[!(b$is_target %in% TRUE), ]
  # the censor carries no number above the window, only the window as a bound
  expect_true(all(is.na(cen$time_to_positive_min)))
  expect_true(all(grepl("exceeds the window", cen$censoring, fixed = TRUE)))
  a <- d$ace1
  expect_identical(length(unique(a$censoring)), 1L)
  expect_true(grepl("not measured", a$censoring[1L], fixed = TRUE))

  # no row in either file is one template measured under BOTH alleles, which
  # is why clause 1 forbids an allele gap.  Every (template, primer set) pair
  # appears at most once.  All five badolo templates are run under both sets,
  # so "two rows for one template" means two PRIMER SETS on that genotype in
  # two reactions - never one template read for two alleles at once.  In
  # is_target terms the template is still the same genotype, and where the
  # other set's row carries a number it is that set's own target allele.
  key_b <- paste(b$template, b$primer_set)
  expect_identical(anyDuplicated(key_b), 0L)
  expect_identical(
    sort(names(which(table(b$template) == 2L))),
    sort(unique(b$template))
  )
  expect_identical(nrow(b), 2L * length(unique(b$template)))
  # `is_target` is a property of each REACTION, not of the template: the single
  # alleles (plasmid, homozygote) are each a target under one set and a
  # non-target under the other; only the heterozygote carries both, so both of
  # its rows are target arms - and it is the slowest of each.
  hom <- c("plasmid wild-type", "plasmid kdr-w", "gDNA +/+", "gDNA kdr-w/kdr-w")
  expect_identical(
    sort(b$is_target[b$template %in% hom]),
    c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
  )
  expect_identical(b$is_target[b$template == "gDNA kdr-w/+"], c(TRUE, TRUE))
  # and the censored rows are each set against the allele it does NOT read: the
  # wild-type set on the two kdr-w templates, the kdr-w set on the two
  # wild-type templates
  expect_identical(paste(
    b$template[!(b$is_target %in% TRUE)],
    b$primer_set[!(b$is_target %in% TRUE)]
  ), c(
    "plasmid wild-type kdr-w", "plasmid kdr-w wild-type",
    "gDNA kdr-w/kdr-w wild-type", "gDNA +/+ kdr-w"
  ))
  for (ps in c("wild-type", "kdr-w")) {
    sel <- b$primer_set == ps & b$is_target %in% TRUE
    expect_identical(sum(sel), 3L)
  }
  key_a <- paste(a$assay, a$primer_set, a$template_genotype)
  expect_identical(anyDuplicated(key_a), 0L)
  expect_true(all(a$assay == "AS-LAMP, G119S ace-1R"))
  # no file carries a per-allele Tt column, a Ct, or a predicted gap
  for (f in list(b, a)) {
    expect_false(any(grepl("ct|allele_gap|delta_tt|dpred", tolower(names(f)))))
  }
})

test_that("the censored arms have an independent anchor in the landed counts", {
  # badolo: AS-LAMP called 70 +/+ of 70 sequenced, i.e. the kdr-w set never
  # fired on a wild-type template - the same statement as its censored arms.
  # ace-1: AS-LAMP called 8 of 9 sequenced heterozygous samples; the
  # heterozygote is the one template whose non-target arms that paper censors.
  cb <- lw_read("lamp_badolo_kdr2012_counts.csv")
  lamp <- cb[cb$method == "AS-LAMP", ]
  wt <- lamp[lamp$genotype == "+/+", ]
  expect_identical(wt$n_called, 70L)
  expect_identical(wt$n_sequenced, 70L)
  b <- lw_read("lamp_badolo_kdr2012_time_to_positive.csv")
  expect_identical(nrow(b), 10L)

  ca <- lw_read("ace1r_aslamp2015_counts.csv")
  rs <- ca[ca$method == "AS-LAMP" & ca$genotype == "Ace1RS", ]
  expect_identical(rs$n_called, 8L)
  expect_identical(rs$n_sequenced, 9L)
  m <- lw_read("ace1r_aslamp2015_metrics.csv")
  rs_sens <- m[m$genotype == "Ace1RS" & m$metric == "sens" & m$method == "AS-LAMP", ]
  expect_identical(nchar(as.character(rs_sens$value)), 4L)
  expect_identical(round(rs_sens$value, 2), round(8 / 9, 2))
  expect_true(rs_sens$matches_counts %in% TRUE)
})
