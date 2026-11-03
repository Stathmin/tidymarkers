# ROADMAP §V — readout models: third search pass (2026-11-03, Agent B)

No edit to `R/marker_readout.R`, no criteria change, no flip of `simulated`. This pass
re-reads only open-access full text (Europe PMC `fullTextXML`) and confirms what the two
prior sweeps (`v-readout-conditions.md`, `v-readout-sources.md`) already established, with
one genuinely new source. Rule from §V kept literal: a model may switch on only when a
table exists whose rows it must reproduce, with every condition it reads printed in the
same source.

## Summary verdict

| chemistry | new set found? | file written | conditions printed | sufficient for `simulated = TRUE`? | remaining blocker |
|---|---|---|---|---|---|
| HRM melt | no | — | — | no | no amplicon-level table prints Mg²⁺ / salt next to a measured Tm; the 28 diatom amplicons use a proprietary master mix (`mg_mm` stays `NA`) |
| dPCR Poisson | n/a (already sufficient) | — | — | yes (unchanged) | none for the Poisson term |
| CRISPR-Dx collateral | **yes (partial)** | `readout_crispr_nalefski2024_pass2.conditions.csv` (109 rows) | kcat/KM, kcat, KM per trans-substrate × RNP × target × 37/25 °C (Table 1); kact(max), Kd, kon(app) full/limiting, koff(net), n per RNP × DNA target (Table 2) | **no** | a kinetic-constants table, not a printed `target concentration → fraction cleaved at time t` series; that time course is still in figures/supplement |
| virtual gel | no | — | — | no | no open paper prints bp → migration distance (mm) with gel %, voltage and time; Southern 1979 / Elder 1983 are paywalled |
| LAMP/RPA kinetics | no | — | — | no (Q1 open) | no open-access table prints absolute copies *and* time-to-positive in one series; closest are regression-only (fitted) or endpoint (no time) |

## CRISPR-Dx collateral cleavage — new source (Nalefski 2024)

**Found.** Nalefski EA, Kooistra RM, Parikh I, Hedley D, Rajaraman S, Madan D, *et al.*
"Determinants of CRISPR Cas12a nuclease activation by DNA and RNA targets." *Nucleic
Acids Res* 2024;52(8):4502–4522. `10.1093/nar/gkae152`, `PMC11077072`. Read from the
Europe PMC `fullTextXML` on 2026-11-03.

This is the first CRISPR-Dx kinetics source read in this repo that prints the
trans-cleavage and activation numbers **in the body as tables**, not as figure panels or
a gated supplement. It is *not* the iScience 2021 paper (`PMC8411246`, already read) and
*not* the NAR 2025 paper whose kinetics are all in a proof-of-work-gated supplement.

- **Table 1** ("Steady-state kinetic constants for C10 trans-cleavage by target-activated
  crRNA-Cas12a RNP") — 19 rows: `kcat/KM`, `kcat`, `KM` per trans-substrate (FB-C10,
  FQ-C10) × RNP (F4, F2, A2, A19, E27, A1) × target (dsDNA, ssDNA, RNA-t) × temperature
  (37 °C, 25 °C, both printed). Example verbatim: FB-C10, 37 °C, F4, dsDNA → kcat/KM
  0.97 (±0.05) ×10⁸ M⁻¹s⁻¹, kcat 2.77 (±0.05) s⁻¹, KM 29 (±2) nM.
- **Table 2** ("Activation parameters of crRNA-Cas12a RNP reacted with DNA targets") —
  11 rows: `kact(max)` (min⁻¹), `Kd` (nM), `kon(app)` full-RNP and limiting-RNP range
  (×10⁷ M⁻¹s⁻¹), `koff(net)` (min⁻¹), and `n`, per RNP (A2, A19, F2, F3, F4, F23, F25) ×
  DNA target (ds, ss, ss-t). Section B (A19/ds) prints only `kon(app)` 0.0051 (±0.0001)
  and `koff(net)` 0.043 (±0.001).
- **Buffer (Methods):** 10 mM Tris–Cl pH 7.5, **10 mM MgCl₂**, 1.0 mM TCEP, 0.01 % Igepal
  CA-630, 40 µg/mL BSA. Table 1 prints temperature; Table 2 does not, so its
  `temperature_c` is `NA` (not guessed).

**Why it is still not sufficient.** The §V target for Cas12a remains a printed
**target concentration → fraction cleaved at time t** (or target concentration → initial
rate) series at a stated enzyme/reporter concentration and temperature. Table 1/2 give
the steady-state and activation constants, but the concentration × time course itself is
still only in figures (e.g. Fig. 4 time courses, Supplementary Table S1). A model could
be *parameterised* from these rows, but it would not be *reproducing a printed
concentration→fraction-cleaved table*, so `simulated = TRUE` stays off.

## HRM melt — not found

Repeated the amplicon-level Mg/salt search (Europe PMC): no open paper prints a
multi-row `MgCl₂/salt → measured Tm` table for amplicons. The best prior candidates
(Oliveira 2020 `PMC8163305`, fixed salt; Huguet 2017 `PMC5728412`, salt→energy single
molecule) remain the only salt-bearing sources and were already rejected as not
amplicon-level HRM. The blocker is unchanged: the 28 diatom amplicons run in a
proprietary 2× master mix, so `mg_mm` is `NA` by design.

## Virtual gel — not found

No open-access paper prints a `bp → migration distance (mm)` ladder with agarose %,
voltage and run time. Ziraldo 2019 (`PMC6895257`) still supplies only the *formula*
(Eq. 2) and the printed conditions (1 % agarose, 3.0 V cm⁻¹, 16 h / 2 h), no numeric
bp→mm table. The original calibrations — Southern 1979 `10.1016/0003-2697(79)90235-5`
and Elder & Southern 1983 `10.1016/0003-2697(83)90369-x` — are **paywalled, not open
access**; they were not read and remain "not in open access", not "does not exist".

## LAMP / RPA kinetics — not found (Q1 still open)

No open-access table prints **absolute copies (or genomes) per reaction and
time-to-positive together** with non-detections as `> run length`.

- `PMC12941333` (LAMP vs ddPCR, *Y. pseudotuberculosis*) prints absolute DNA copies per
  reaction for the LAMP detection limit (Table 2: blue dye 1.43×10², violet 1.43×10²;
  Table 5: 1.25×10³ / 1.48×10²), but the readout is **endpoint colorimetric** — no
  time-to-positive column. Copies without time does not close Q1.
- `PMC9349938` (quantitative RT-qLAMP, SARS-CoV-2) prints the copy-number → TT (threshold
  time) relationship only as **regression equations** (Table 1: slope, intercept, R² per
  variant), not the underlying per-dilution (copies, TT) rows; the raw pairs are in
  Fig. 1B/2B. A fitted line is not a row table.
- `PMC12093978` (*M. ulcerans* P-LAMP) prints only the Tp↔Ct correlation (r = 0.9818) and
  "4 vs 2 genome copies" LoD in prose; the Tp values are in Fig. 1B, not tabulated.
- The classic quantitative LAMP source — Cai et al. 2008, *J Clin Virol*,
  `10.1016/j.jcv.2007.11.025`, "real-time LAMP for hepatitis B virus DNA quantification"
  — is **paywalled (not open access)**: Europe PMC reports `isOpenAccess=N`,
  `inEPMC=N`, so its copies→Tt table could not be read.

RPA: nothing new; the `PMC13377618` (PEG crowding) fulltext has **0 table-wraps**.

## Files written

- `data-raw/proposed/readout_crispr_nalefski2024_pass2.conditions.csv` (109 rows).
- `data-raw/proposed/v-readout-pass2.md` (this file).

No `R/`, `tests/`, `inst/extdata/`, `ROADMAP.md`, `NEWS.md` or existing CSV touched. No
R run.
