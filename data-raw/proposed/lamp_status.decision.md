# Decision: what the LAMP-SNP lane may claim

Status: **decided 2026-10-22.** Closes `lamp_status.proposal.md`. **No `R/`
model change**: LAMP stays a construct designer with a time-window readout. One
data correction is applied (point 6).

## The decision

1. **No kinetics, no specificity number.** The package may say that a LAMP
   design is a structural construct whose readout is a time window. It may not
   report a predicted allele gap, a Ct/Tt difference between alleles, or a
   false-call probability, because no landed row is a per-template Tt under more
   than one allele. `readout_assumption("lamp")` is the only place the readout
   speaks, and it says `simulated = FALSE`.
2. **Do not fit a model to 8 arms over 2 SNPs** on one isothermal turbidimeter.
   That would be fitting the points the model is meant to predict.
3. **The 2.7x ratio is not "LAMP specificity".** `lamp_comp2018` is endpoint
   fluorescence between synthetic oligos at a fixed duration; it is never quoted
   beside a Tt-based window.
4. **One placement rule is evidenced: FIP/F2 3' terminal.**
   `mismatch_rtlamp2019`: a 3'-terminal FIP mismatch delays Tt by 13.5-17.3 min
   against 1.1-3.5 min at F3 and 6.9-13.4 min at F1c. This matches what
   `R/marker_lamp.R` already implements.
5. **A second deliberate mismatch is not calibrated.** The two sources that use
   one disagree on its offset (penultimate vs 3' end in `badolo_kdr2012`; B2
   positions 1 and 3 in `snp_lamp_lfd2025`, whose per-base figure is an image),
   and neither isolates the second mismatch's effect. The package implements no
   second-mismatch offset.
6. **Discrimination is a window with censoring, not a binary call.** Non-target
   arms in both mosquito papers are "no amplification within 75 min", i.e.
   `Tt > 75`; a REF-only template and a failed reaction stay indistinguishable.
   This is already stated in `R/marker_lamp.R` L44-48 and the landed censored
   arms support it.
7. **The negative result is stated, not hidden.** No systematic scan of LAMP
   primer mismatch tolerance by position, for allele discrimination, was found
   (searches s8-s13; s9 and s12 return 0). The nearest artefact
   (`mismatch_rtlamp2019`) is a tolerance scan under 0.15 U Q5 proofreading
   polymerase per 25 uL — a condition the package's designs do not assume. The
   package advertises no mismatch-position guidance beyond point 4.

## The `badolo2012` key correction (data)

`inst/extdata/scientific_backlog.csv` row 5 keys `badolo2012` to
`10.1186/1475-2875-11-232`, a resistance-monitoring survey with no primer
content, and stamps it `extracted` with a note since corrected. The LAMP paper
the row's anchor text describes is `10.1186/1475-2875-11-227` (PMC3407793,
PMID 22770418), landed under its own key `badolo_kdr2012`. The backlog row is
already right; the **source** record was not: `scientific_sources.csv` carried
`10.1186/1475-2875-11-231`, which resolves to neither article, with no PMID.
Verified against Europe PMC (PMC3489511 -> PMID 22799568,
doi ...-11-232) and corrected in place. This is a key/DOI correction, not a
reversal of the verdict.

## Resulting `R/` change

None. `R/marker_lamp.R` already states the window and the failed-reaction
indistinguishability; the decision freezes that as the claim boundary.

## Test that asserts it

`tests/testthat/test-lamp-status.R`: `readout_assumption("lamp")` is not
simulated; no LAMP evaluation row carries a predicted allele gap or a Ct/Tt
difference; the FIP 3' base is the one the designer places the variant at; the
2.7x separation row is never labelled a specificity figure; and every source
record's DOI, PMID and PMCID resolve to the same article on Europe PMC
(`badolo2012` = ...-11-232 / PMID 22799568 / PMC3489511, `badolo_kdr2012` =
...-11-227 / PMID 22770418 / PMC3407793; the two are different articles).
