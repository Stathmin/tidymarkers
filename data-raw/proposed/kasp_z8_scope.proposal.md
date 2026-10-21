# KASP conversion: which rates may be decomposed into design vs work

> **DECIDED 2026-10-21** — see `kasp_z8_scope.decision.md`: no `R/` change; the guard lives in `tests/testthat/test-z8-outcomes.R`. The status
> line below is the state this file was written in; it is kept as the evidence
> record and is not edited after the decision.

Status: proposal from the sci-kasp worker. **Nothing in this file has been
applied** to `inst/extdata/conversion_rates.csv` or `R/`. Per the brief, that
table is not edited; this note states what the three new extracts do and do not
license.

## 1. The family, in one paragraph

A KASP "conversion rate" is at least four different measurements, and papers
print whichever one they have:

| stage | what it counts | who reports it |
|---|---|---|
| **design** | a primer/probe triple exists for the SNV | Makhoul (3/8/11 of 15), Sharplab (50/50), Jagtap (100/100), Ayalew (49/50 KASP) |
| **manufacture** | the triple was made and put on a plate | Ayalew only (30 of 39) |
| **work** | on real DNA the assay produced a call | Ayalew (per sample point), Jagtap (92/100) |
| **usable / polymorphic** | the call was two-allele and informative on the material | Jagtap (71/100), Sharplab (46/50 between two inbreds), Ayalew (implicitly, 13 invalid) |

The new extracts make this explicit with four columns that every KASP row in
this family now carries: `design_attempted`, `design_success`, `assay_worked`,
`outcome_definition`. `design_success` is never a wet-lab result, and
`assay_worked` is never a design result.

## 2. The three new sources, precisely

| key | design | work | usable | what the denominator is |
|---|---|---|---|---|
| `wheat_platform2019` | KASP 49/50 (TaqMan 39/50, rhAmp 49/50) | 6.5% of **2,784 assay x sample points** unamplified | 13 invalid calls out of those points; **no per-assay outcome is printed** | 50 pre-validated wheat SNVs; then 96 x 29 = 2,784 points |
| `sharplab2021` | 50/50 (a synthesised triple per marker) | 98.5% amplification ratio, **no count printed** | 46/50 polymorphic between **the two inbreds 1212 and B73** | 50 selected markers |
| `maize_kasp2020` | 100/100 | 92/100 = (71 + 21)/100 | 71/100 polymorphic on **LM11 x CML25** | 100 RNA-seq SNVs |

Two traps the extracts are built to avoid:

* `wheat_platform2019`'s **6.5%** is a fraction of assay x sample points, not
  of assays. It must never be written as "KASP failed 6.5% of its assays".
* `sharplab2021`'s **92%** and `maize_kasp2020`'s **92%** are numerically equal
  and mean different things (polymorphism between one inbred pair; calls
  produced on a parent pair). Any table that puts them in one column is wrong.

## 3. Which existing rows become decomposable? None.

Row by row, against the six rows of `conversion_rates.csv`:

| row | decomposable into design vs work after this pass? | why |
|---|---|---|
| `grewal2020` | **No** | 1000 -> 730 is a conversion of a *validated* SNP to a usable chromosome-specific assay in a wild-relative cross. It contains no separate design stage: primer design is not reported as a funnel of its own. The 270 failures are not split by the paper into design / amplification / polymorphism. |
| `ramirezgonzalez2014` | **No** | 35 tested / 28 polymorphic is a polymorphism count in a biparental cross with no design stage reported at all. |
| `grewal2022_amblyopyrum` | **No** | 698 -> 335 already names five failure classes (251 PCR, 49 no Am. muticum allele, 22 wheat-internal, 10 monomorphic, 31 heterozygote). Those are **work** classes only; design is not reported as a separate denominator. It is the closest existing row to a decomposition and still cannot be split. |
| `makhoul2020_bp_flanks` | **No** | 3/15 is a **design-and-cluster-inspection** rate: "Correct cluster" is a judgement made at the design stage on simulated clusters. There is no separate amplification stage and no sample panel, so `assay_worked` is undefined for it. |
| `makhoul2020_single_reference` | **No** | same, 8/15. |
| `makhoul2020_multi_assembly_sanger` | **No** | same, 11/15. This is the only row with `fully_decomposed = TRUE`, and that flag refers to the *per-SNP outcome reasons* (from `makhoul_design_strategies.csv`), not to a design/work split. The two must not be confused: `fully_decomposed` here means "we know why each of the 15 failed", not "we know design and work separately". |

**Conclusion:** no row of `conversion_rates.csv` can be split into design and
work with the files we now hold, because none of the six reports both stages on
the same denominator. The new sources are **additive**, not corrective: they are
new rows, and they must **not** be pooled with the CAPS/raw-band rows already in
`conversion_rates.csv` (Liu 2016, soybean, Snap2Caps), which measure a band on a
gel, nor with the six KASP rows above, whose rate is a different stage.

## 4. What must stay flagged as not comparable

1. **Any rate whose denominator is assay x sample points** (`wheat_platform2019`
   amplification and invalid-call fractions). Not comparable to an assay-level
   conversion rate at all; keep it in a separate table with `scope =
   "assay_x_sample"`, which the extract now does.
2. **Any rate defined on a single inbred pair** (`sharplab2021` 46/50,
   `maize_kasp2020` 71/100). These are polymorphism rates *on that pair* and
   move with the pair. Two different pairs giving 92% proves nothing about
   either.
3. **In-silico yields** (`sharplab2021` 71,311 / 2,948,985 = 2.42%). A sequence
   criterion, not an assay.
4. **Cluster-quality metrics** (`wheat_platform2019` separation angle,
   compactness, NTC distance, call quality). These are read off the homozygous
   cluster only (the paper's Allele-11 footnote); they are not success rates.
5. **"Design ambiguous" numbers.** In `wheat_platform2019`, "Thirty-nine out of
   the 50 ... passed assay design criteria across the three platforms" (Methods,
   an intersection) reuses the same 39 as TaqMan's per-platform 39/50
   (Results). The extract records both and flags the Methods one
   `design_ambiguous`. Do not cite it as KASP's design rate.

## 5. What a future pooled table would need

If a design/work decomposition is ever placed in `conversion_rates.csv`, each
row needs: the design denominator (SNVs attempted), the design numerator, the
work denominator (assays that produced a call), the work numerator, the usable
numerator (informative calls), and `conversion_definition`. Only
`maize_kasp2020` and `wheat_platform2019` supply a design numerator at all, and
only `maize_kasp2020` supplies a work numerator on the same denominator. That is
one source; a pooled rate built on one source is not a pooled rate.

Recommended guard for `test-z8-outcomes.R` when this is merged: any future
KASP row in `conversion_rates.csv` must carry a non-NA `conversion_definition`
and must not share a value with a row whose definition differs.

## 6. Not found

The 648-locus wheat set (326/648 = 50.3%, from ~791,000 SNPs) was **not
located**, and no row is proposed for it. Searches run are listed verbatim in
`.tmp/sdd/report-kasp.md` section 5, with the hits they returned. The closest
thing EPMC returned for "648" + KASP + Triticum was a set of papers that merely
mention 648; none reports that conversion. Two candidates were inspected at
title level only and rejected (a cowpea fingerprinting paper and a watermelon
powdery-mildew paper), because neither is wheat.

Also not found in `sharplab2021`: the "cross-crop conversion survey" the brief
expected. The paper compares its average PIC (0.13) with sesame, radish,
grapevine and sweet orange PIC values -- a PIC survey, not a conversion survey.
Nothing cross-crop was landed.
