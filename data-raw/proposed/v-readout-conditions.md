# ROADMAP §V — readout simulation: which condition → outcome sets now exist

Agent A6. Date 2026-10-29. No edit to `R/marker_readout.R`: every chemistry below
except hydrolysis probes still returns `simulated = FALSE`, and this document says
what would have to change (and reproduce what) before that flips.

The rule from §V is kept literal: a model may be switched on **only** when a table
exists whose rows the model must reproduce, with every condition a model reads
printed in the same source. A paper that only *draws* the relationship is not a set.

## Summary verdict

| chemistry | new set written here | existing sets | sufficient for `simulated = TRUE`? | blocker |
|---|---|---|---|---|
| HRM melt | `readout_hrm_diatom2025.csv` (28 rows) | `hrm_jsbbs2018_*`, `hrm_ijms2022_*`, `hrm_gundry2008_*` | **partly**: enough for a Tm-prediction model, **not** for an HRM discrimination model | no Mg / salt printed (proprietary master mix) |
| dPCR Poisson | none — no third independent raw-count set found | `dpcr_dong_poisson_residual.csv` (60 rows, 4 platforms), `dpcr_mps_ecoli_dilution.csv` (102 rows) | **yes, the existing two suffice** | none for the Poisson term; mutant-fraction estimation untested |
| CRISPR-Dx collateral | none | `communbiol2024_crrna_ratio.csv` (320 rows, ratio only) | **no** | every kinetics/concentration number lives in figure panels or a proof-of-work-gated PDF supplement |
| virtual gel | none | `hrm_liew2004_restated_bands.csv` (citation restatement, not migration) | **no** | no open paper prints bp → migration distance (mm) with gel %, voltage and time |
| LAMP/RPA kinetics | `readout_lamp_erysiphe2026.csv` (30 rows) | `lamp_badolo_kdr2012_time_to_positive.csv` (10), `ace1r_aslamp2015_time_to_positive.csv` (4), `lamp_mismatch_rtlamp2019_tt.csv` (26) | **no (Q1 not closed)** | the dilution series is a *mass* dilution of an unquantified extract: no copies per reaction |
| hydrolysis probe | — | existing | **yes (already `simulated = TRUE`)** | — |

## HRM melt

**Set.** `inst/extdata/readout_hrm_diatom2025.csv`, 28 rows (Zhou 2025,
`10.1371/journal.pone.0321885`, `PMC12002481`).

One HRM run melts 28 amplicons: 5 diatom species × 4 screening primer pairs
(Table 6, 20 cells) plus psaA-2 (5) and Dia18S (3) whose *measured* Tm is printed
in Table 9. For every amplicon the paper also prints the Sanger product sequence
(Table 7), its GC%, bp, summed ΔH° and ΔS° (Table 8). Conditions are constant and
stated for the whole run: EvaGreen, 35 cycles, melt ramp 60 → 95 °C at 0.1 °C/s.

**Is it sufficient?** For a **Tm-prediction** model, yes — this is the first set in
the tree where the *input sequence, its nearest-neighbour parameters and the
measured Tm* are all printed in one place, so a model can be scored per amplicon
(`tm_measured_c` vs a prediction from `dh_kcal_mol`/`ds_cal_k_mol`, or from the
sequences themselves). For an **HRM discrimination** model (what `marker_readout()`
would be claiming), no: there is no `dTm`, and no salt term can be fitted because
Mg²⁺ is inside a proprietary 2× master mix and the paper prints no concentration.
`mg_mm` is therefore `NA` in every row, deliberately.

**What a model would have to reproduce.** (1) Each of the 28 Tm values within a
tolerance the paper itself sets: its nearest-neighbour method averages
0.67 °C error over the eight Table 9 cells (mean of |measured − predicted|).
(2) The ordering within a primer pair, e.g. Primer1 ranks Chaetoceros > Odontella >
Nitzschia > Navicula > Skeletonema. (3) Any salt correction must be invisible to
this set — if a model needs Mg to reproduce these numbers it will be fitting noise.

**Paper defects kept as columns, not repaired.** Table 8 GC% disagrees with its own
Table 7 sequence in 3 of 25 sequenced amplicons (`gc_printed_vs_recomputed_slip`);
Table 7 sequences are +1 nt relative to Table 8 bp in 23 of 25
(`seq_len_vs_table_bp`); Dia18S prints no sequence; Table 8 prints 59 and 154 bp
while the prose claims "60 to 152". Table 4 (Tm vs DNA input and cycle number) is
**excluded** because it contradicts Table 6 for psaA-2 Navicula (82.38 vs 83.40),
and its PLOS supplementary S1 Table could not be read unambiguously.

## dPCR Poisson

**Set.** No new file. The two existing sets are sufficient and are already the test
fixture (`tests/testthat/test-readout-physics.R` asserts
`|lambda_poisson − printed| ≤ 0.005` on Dong):
`dpcr_dong_poisson_residual.csv` prints raw positive/negative/total counts *and* the
instrument's own λ for 4 platforms × 15 replicates, so the Poisson inversion is
checked against the instrument at λ ≈ 1.4–1.6; `dpcr_mps_ecoli_dilution.csv` covers
the same inversion across a 6.6–7900 copies/µL dynamic range (102 rows).

**Is it sufficient?** Yes for the counting model. A third independent set was
searched for and **not found** in open access: papers printing "positive partitions"
almost always print a *concentration* plus a partition acceptance threshold, not the
per-well counts (checked `PMC13188893`, `PMC12622311`, `PMC12644700`, `PMC12361875`,
`PMC11472245`). What is *not* constrained is mutant-fraction estimation from two
channels (`jak2v617f_ddpcr_*` gives limits, not counts).

**What a model would have to reproduce.** λ = −ln(N_neg/N_total) from the printed
counts within 0.005 of the printed λ in all 60 Dong rows, and linear recovery of the
gravimetric dilution in the MPS rows.

## CRISPR-Dx collateral cleavage

**Set.** None written. `communbiol2024_crrna_ratio.csv` (320 rows) gives
mutant/WT fluorescence *ratios* per crRNA, not a concentration or time course.

**Is it sufficient?** No. Every source read puts the kinetics in figure panels or
supplements: `PMC12604663` (NAR 2025, Cas12a trans-cleavage kinetics — the best
candidate found) has all per-RNP k_act, k_ss and LOD numbers in Supplementary
Tables S2–S17, whose PDF is behind a proof-of-work interstitial on PMC and 403 on
OUP; the search for "target concentration → endpoint fluorescence" in
`PMC12304900` and similar returned figures with no printed values.

**What a model would have to reproduce** (so the next attempt knows the target): a
table of target concentration → initial rate, or time → fraction of reporter
cleaved, at a stated enzyme/reporter concentration and temperature. Without that,
`simulated = TRUE` for Cas12a would be a claim about curves that are not numbers.

## Virtual gel (CAPS / ARMS / indel / SSR products)

**Set.** None written.

**Is it sufficient?** No. The candidate papers read
(`PMC8888377` Raspberry-Pi imaging, `PMC6418129` indel genotyping on 4% agarose,
`PMC8427238` hydrodynamic chromatography, `PMC12053679` GelGenie) all report gel
images, "resolved / not resolved" statements or chromatographic mobility — none
prints band position (mm or px) against bp together with agarose %, voltage and
run time. `hrm_liew2004_restated_bands.csv` is a *citation* of Liew's dTm bands and
is not a migration table.

**What a model would have to reproduce.** A ladder calibration: ≥8 points of
(bp, migration distance) in one lane, with agarose %, buffer, voltage and duration
printed. Sequence length alone is not a gel; without distance the "band sizes
reported" note stays an arithmetic report, not a simulation.

## LAMP / RPA isothermal kinetics

**Set.** `inst/extdata/readout_lamp_erysiphe2026.csv`, 30 rows (Barone 2026,
`10.3390/jof12010079`, `PMC12842720`). Two grids from one SYTO-dye LAMP run
(60 min, 2 replicates, 60 °C unless noted):
- `experiment = temperature`: 3 templates × 60/63/65 °C, all nine cells detected;
- `experiment = dilution`: 3 templates × ten-fold **mass** dilution (20 ng/µL down to
  10⁻⁶), with the 10⁻⁶ cells printed "ND" and stored as censored at the 60-min run.

**Is it sufficient?** **No — Q1 is not closed.** The set improves on what was in the
tree (it is the first LAMP table with a dilution axis *and* a temperature axis in one
run, and it is the only one where the censoring is a stated run length rather than a
"not determined"), but `copies_per_reaction` is `NA` for every row: the paper dilutes
an extract whose concentration it never prints. A time-to-positive model needs
molecules per reaction; against relative mass it cannot be scored, only ordered.

**What a model would have to reproduce.** Tp rising monotonically with dilution
(asserted on the real rows) and slowing at 60 → 63 °C. A buildable set would need a
paper printing **copies (or genomes) per reaction and Tp together**, with the
non-detections given as ">60 min". Existing `lamp_badolo_kdr2012_time_to_positive.csv`
has exactly that censoring form but only two templates at one copy level.

RPA: nothing new. Existing extracts are call-agreement and mismatch-position tables,
not kinetics.

## What stayed unbuildable

CRISPR-Dx (no printed concentration series), virtual gel (no printed migration
distance), and LAMP/RPA **as a function of copies** — the last is the only one where
the qualitative relationship *is* tabulated but the driving variable is not.
