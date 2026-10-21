# LAMP-SNP — what this package may and may not claim

Status: **proposal only.** Nothing in `R/` has been edited. The coordinator decides.
Every clause cites a landed CSV row, a line of `R/`, or a recorded negative search.

## 1. The package has no LAMP kinetics or specificity model — and must not gain one from these sources

`R/marker_readout.R` returns `simulated = FALSE` for LAMP with the note
`"isothermal amplification kinetics not modelled"`. `readout_assumption("lamp")`
is the only place `marker_readout()` can speak about a LAMP design.

The four sources landed in this pass **cannot calibrate a kinetics model**:

| source | file | what the readout actually is |
|---|---|---|
| `lamp_comp2018` | `lamp_comp2018_separation.csv` | **endpoint fluorescence**, 2.7× SNP/non-SNP, at a fixed 30–75 min duration (`lamp_comp2018_durations.csv`). No Tt, no Ct. |
| `snp_lamp_lfd2025` | `lamp_snp_lamp_lfd2025_detection.csv` | **visual dipstick**, 60 min + dipstick read; a 2×2 call table (6/0 and 0/4), no magnitude. `evidence_tier = weak`. |
| `badolo_kdr2012` | `lamp_badolo_kdr2012_time_to_positive.csv` | **turbidity Tt**, 4 measured arms (50/52/60/65–75 min), 4 non-target arms **censored** at 75 min. |
| `ace1r_aslamp2015` | `ace1r_aslamp2015_time_to_positive.csv` | **turbidity Tt**, 4 arms, 75-min window, same censoring. |

**Clause 1.** The package may say LAMP designs are structural constructs with a
time-window readout. It may **not** report a predicted allele gap, a Ct
difference, or a probability of a false call, because no landed row is a
per-template Tt under more than one allele.

**Clause 2.** Two of the four sources (badolo, ace1r) each supply **four**
time-bearing arms on **one** SNP, in **one** species, on an isothermal
turbidimeter. That is 8 arms total, over 2 SNPs. A model fitted to that would be
fitted to the two points it is meant to predict. Do not fit one.

**Clause 3.** Only `lamp_comp2018` reports a signal *ratio* (2.7×;
`lamp_comp2018_separation.csv`), and that ratio is between **synthetic DNA
oligos** in an endpoint fluorescence readout — not a Tt, and not genomic DNA.
It must not be quoted as a "LAMP specificity" figure next to a Tt-based window.

## 2. What the mismatch-placement evidence actually supports

**Clause 4 — FIP/F2 3′ end is evidenced as the discriminating position.**
`mismatch_rtlamp2019_position_rules.csv`: a 3′-terminal mismatch at FIP's F2 arm
delays Tt by **13.5–17.3 min** under the conventional method, versus **1.1–3.5
min** at F3 and **6.9–13.4 min** at the F1c 5′ end. This is the only source with
a measured position ranking, and it agrees with `R/marker_lamp.R`'s choice of
the FIP 3′ base (lines 45–48, `lamp_fip_alt` in `R/marker_criteria.R`).

**Clause 5 — BIP/B2 3′ end is evidenced, with a position conflict.**
`badolo_kdr2012`: the two BIP sets differ only at the last two bases of B2,
`…CGA gT` vs `…CGA gA` (`lamp_badolo_kdr2012_primers.csv`). The paper calls the
extra mismatch **penultimate** (Methods, Results Q17) and **3′ end** (Discussion
Q28) — the contradiction is recorded in
`lamp_badolo_kdr2012_contradictions.csv`, not resolved.
`snp_lamp_lfd2025` puts modifications at **B2 positions 1 and 3 from the 3′
end** (`lamp_snp_lamp_lfd2025_primers.csv`), but its Figure 1 (the per-base
diagram) is an image, so the per-base identity is **not** landed, and its own
Methods vs Discussion disagree on what those positions hold
(`lamp_snp_lamp_lfd2025_contradictions.csv`).

**Clause 6.** The package may claim **one** placement rule it already
implements (FIP/F2 3′). It may **not** claim a second-mismatch rule as
calibrated: the two papers that use one disagree on its offset, and neither
reports the effect of the second mismatch alone.

## 3. The negative result (priority c) — state it, with the searches

**Clause 7.** No systematic scan of LAMP primer mismatch tolerance **by
position, for allele discrimination**, was found. The nearest artefact,
`mismatch_rtlamp2019`, is a mismatch-**tolerance** scan (how much one mismatch
slows one primer set) and its "novel" arm requires **0.15 U Q5 proofreading
polymerase per 25 µL** on top of Bst (`mismatch_rtlamp2019`, Methods) — a
condition the package's designs do not assume. Its own numbers are a delay,
never a discrimination ratio: in the conventional arm the mutants are 1.1–22.6
min slower, but the slowest (44.3 min) is still positive.

Searches run this pass against Europe PMC (`.tmp/probe/s1..s13.json`; hit counts
are the API's, queries verbatim in the corresponding JSON `nextPageUrl`):

| id | query | hits |
|---|---|---|
| s8 | `TITLE:"loop-mediated" AND ABSTRACT:"mismatch" AND ABSTRACT:"allele"` | 1 |
| s9 | `ABSTRACT:"LAMP" AND ABSTRACT:"single nucleotide" AND ABSTRACT:"mismatch" AND ABSTRACT:"position"` | 0 |
| s10 | `ABSTRACT:"loop-mediated isothermal" AND ABSTRACT:"allele-specific" AND ABSTRACT:"specificity"` | 17 |
| s11 | `ABSTRACT:"isothermal" AND ABSTRACT:"SNP" AND ABSTRACT:"primer" AND ABSTRACT:"mismatch" AND ABSTRACT:"discrimination"` | 1 |
| s12 | `TITLE:"loop-mediated isothermal amplification" AND ABSTRACT:"mismatch tolerance"` | 0 |
| s13 | `ABSTRACT:"LAMP" AND ABSTRACT:"primer" AND ABSTRACT:"penultimate"` | 1 (ace1r) |

The 17 s10 hits were screened by title; none reports Tt per mismatch position
across ≥2 positions with a stated enzyme, except `mismatch_rtlamp2019` itself.
**This absence is the finding.** The package must not advertise mismatch-
position guidance beyond Clause 6.

## 4. The `badolo2012` correction — the prior verdict was right for the DOI, wrong for the row

`inst/extdata/scientific_backlog.csv` row 5 keys `badolo2012` to
`10.1186/1475-2875-11-232` (a resistance-**monitoring** survey, no primer
content) and stamps it `not_extractable` with the note "the row describes a
different paper". That is correct: `badolo2012` is **not** a LAMP paper.

But the row's `what_it_would_anchor` ("mutation at the 3′ end of B2 plus a
second deliberate mismatch for specificity (kdr-w)") and its `identifier`
("Badolo A, … Three years of insecticide resistance monitoring …") are two
different articles. The LAMP paper the anchor text describes is
**`10.1186/1475-2875-11-227`, PMC3407793**, and its evidence is landed here
under the new key `badolo_kdr2012` (`lamp_badolo_kdr2012_*`).

**Proposal:** keep row 5 as the correction record with its status unchanged, add
a new backlog row for `badolo_kdr2012` pointing at `…11-227`, and record in
row 5's note that the key was mispaired with a monitoring paper rather than
that the LAMP evidence does not exist. This is a key/DOI correction, not a
reversal of the earlier verdict.

## 5. Summary of what the package may say today

- LAMP is offered as a **construct designer** (F3/B3/FIP/BIP with the variant at
  the FIP 3′ terminus), with explicit non-simulated readout. **(Clause 1)**
- Allele discrimination in LAMP is a **time window with censoring**, not a
  binary call: non-target arms in both mosquito papers are "no amplification
  within 75 min", i.e. `Tt > 75`, never a measured Tt. **(Clause 1, tables)**
- A REF-only template and a failed reaction are indistinguishable — already
  stated in `R/marker_lamp.R` lines 44–48; the landed data supports it
  (every censored arm).
- Nothing about LAMP **specificity numbers**, **mismatch-position tolerance**,
  or **second-mismatch offsets** beyond Clause 6.
