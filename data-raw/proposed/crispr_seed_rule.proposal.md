# Proposed correction to the Cas12a seed rule — from measured data

> **DECIDED 2026-10-21** — see `crispr_seed_rule.decision.md`: the gate stays a POSITION rule; the measured magnitude is `crispr_seed_evidence()`. The status
> line below is the state this file was written in; it is kept as the evidence
> record and is not edited after the decision.

Status: **proposal only.** Nothing in `R/` has been edited. The coordinator decides.
Every clause below cites a row (or set of rows) in a landed extract.

## What the package says now

`R/marker_crispr.R`: `CRISPR_SEED_LEN <- 8L`, and `crispr_seed_gate()` requires
`v_off < 8` for Cas12a — i.e. a variant in the PAM-proximal 8 nt is treated as
discriminated. The docstring calls the seed "the mismatch-critical region".

The rule is a **position** rule with no magnitude. It says a variant in the seed
is detected; it never says how big the signal gap is.

## What was measured (source: Commun Biol 2024, PMC11541961)

`inst/extdata/communbiol2024_crrna_ratio.csv`, 320 rows = 40 activators × 8 crRNAs
(`doi:10.1038/s42003-024-07173-7`, OA full text + source-data workbook fetched).
Outcome: mutant endpoint trans-cleavage fluorescence ÷ WT fluorescence within 1 h
("inverse multiples of specificity", Fig. 1b legend). 1.0 = no discrimination.

### Clause 1 — the unmodified Cas12a guide does not discriminate at all

`crrna = "WT crRNA"`, all 20 single-mismatch rows (`n_mismatch = 1`):
the ratio stays in **0.767 – 1.394** across every spacer position 1–20.
The lowest is `mut1-C1A`…`mut15-T14C` etc.; concretely `T13C` = 0.767 and
`A4G` = 1.207 (`crrna = "WT crRNA"`).

*Proposal:* a bare position rule over an unmodified crRNA must not be described
as conferring specificity. Position alone explains almost nothing here.

Rows: `ratio` rows with `crrna == "WT crRNA"` and `n_mismatch == 1`.

### Clause 2 — even inside the seed, a single mismatch leaks most of the signal

Recommended crRNA `M39`, positions 1–8 (`n_mismatch = 1`), ratio range
**0.373 – 0.968**; the *worst* seed single mismatch (`A6G`) still transmits
**0.373** of WT signal, and the best (`A2G`) transmits **0.968** — i.e. a single
seed mismatch at position 2 is essentially undiscriminated.

*Proposal:* replace "seed < 8 nt is intolerant" with a statement of magnitude —
a single seed mismatch reduces, but does not abolish, collateral signal; the
observed floor is ~0.37 of WT on the paper's own optimised crRNA.

Rows: `ratio` where `crrna == "M39"`, `n_mismatch == 1`,
`mismatch_position` in 1–8.

### Clause 3 — two adjacent mismatches are usually stronger, and this is a trend, not a rule

On `M39`, **17 of 19** adjacent doubles fall below 0.5, versus **8 of 20**
single mismatches. The strongest adjacent double in the whole set is
`A16G&C17A` on `M39` = **0.0158**.

But the ordering **fails at the spacer ends**: `A19G&C20A` is below both of its
constituent singles on only **4 of 8** crRNAs, and on `M65` it reads **1.191** —
worse than either single (`A19G` 0.704, `C20A` 0.897). `T18C&A19G` also holds on
only 6 of 8.

*Proposal:* do **not** encode "two adjacent seed mismatches are always better"
as a hard filter. It holds robustly for pairs starting at positions 7–13
(all 8 crRNAs), weakly at 1–6, and not at all for the last two positions.

Rows: `ratio` where `n_mismatch = 2`, `adjacency = "adjacent"`, joined to the
constituent single rows at the same positions.

### Clause 4 — the registry cannot express a seed rule anyway

`data-raw/proposed/crisprbase_registry.csv` (regenerated from the crisprBase
1.13.2 source tarball; all 8 rows): the `CrisprNuclease` class has slots
`nucleaseName, targetType, motifs, weights, metadata, pam_side, spacer_gap,
spacer_length`. **There is no PFS slot and no seed slot.** "PFS" occurs in
`R/CrisprNuclease-class.R` only as `show()` label text. Both RNA nucleases
(`CasRx` = Cas13d, `Csm`) store the catch-all motif `N`, which constrains nothing.

*Proposal:* a seed model must live in the package, not be delegated to
`crisprBase`. Any `seed_model` column the registry grows is ours to justify.

## Not yet answerable from these data

* Cas13a position weights — Gootenberg 2017 (`10.1126/science.aam9321`) per-position
  numbers are in the non-OA supplement; not landed.
* Cas12a PAM-proximal **PFS** effects — no fetched source gives per-position PFS
  trans-cleavage signal.
* Whether the ratio magnitude predicts LOD — the paper reports LOD separately and
  this extract does not carry it.

## Caveats the coordinator should weigh

1. One target gene (CPSIT_0429), one RPA-Cas12a format, endpoint fluorescence at 1 h.
   Not a kinetic measurement.
2. The ratios are **computed here** from the sheet's raw triplicates, because the
   paper draws the heatmap but never tabulates the ratio. The triplicates reproduce
   the sheet's own `mean` column exactly (max relative difference 0).
3. The crRNA column mixes scaffold/loop/stem engineering (`M37`, `M39`, `M42+M43`,
   `M45`, `M48`, `M63`) with a 3′ spacer truncation (`M65`, 18 nt). Only the 20 nt
   spacers are strictly comparable for position rules; the extract flags spacer
   length but does not silo it.
