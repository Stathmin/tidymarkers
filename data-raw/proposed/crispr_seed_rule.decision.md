# Decision: the Cas12a seed rule

Status: **decided 2026-10-21.** Closes `crispr_seed_rule.proposal.md`.

## What the code said, and what is wrong with saying only that

`CRISPR_SEED_LEN <- 8L` and `crispr_seed_gate()` require `v_off < 8` for Cas12a,
with the docstring calling the seed "the mismatch-critical region". That is a
**position** rule with no magnitude: it says a variant in the seed is
discriminated, and never how much signal separates the alleles. The measured data
(Commun Biol 2024, `communbiol2024_crrna_ratio.csv`, 320 rows = 40 activators ×
8 crRNAs) says position alone explains almost nothing:

1. On an **unmodified** Cas12a crRNA, all 20 single-mismatch rows keep the
   mutant/WT endpoint fluorescence ratio in **0.77–1.39** — no position
   discriminates.
2. On the paper's optimised crRNA (`M39`), seed positions 1–8 range **0.373–0.968**:
   a single seed mismatch reduces, but does not abolish, collateral signal; the
   *worst* seed mismatch (`A6G`) still transmits 0.373 of WT.
3. Two adjacent mismatches are usually stronger (17 of 19 below 0.5 on `M39` vs
   8 of 20 singles), but the ordering **fails at the spacer ends** — `A19G&C20A`
   is better than both singles on only 4 of 8 crRNAs, and on `M65` it is worse
   than either (1.191 against 0.704 and 0.897).

## The decision

1. **The seed stays a POSITION gate.** It is a coarse filter on where a variant
   must sit; it is not a claim that a seed mismatch abolishes Cas12a activity.
   The docstring, the gate label and the report text now say so.
2. **The magnitude is carried as evidence, not as a second gate.** A single seed
   mismatch retains a substantial fraction of signal (0.37–0.97 on the optimised
   crRNA). Since this package cannot predict the ratio for a given target, the
   honest report is `unknown` on the magnitude — never `pass` because the
   position is inside the seed.
3. **"Two adjacent mismatches are always better" is NOT a filter.** It holds for
   pairs starting at positions 7–13, weakly at 1–6, and not at all for the last
   two positions. Encoding it would delete designs the evidence cannot condemn.
4. **A seed model lives in this package, not in `crisprBase`.** The crisprBase
   1.13.2 registry (8 rows, regenerated from the source tarball) has no PFS slot
   and no seed slot; "PFS" occurs only as `show()` label text, and both RNA
   nucleases store the catch-all motif `N`.
5. **No PFS effect is claimed for Cas12a.** No fetched source gives per-position
   PFS trans-cleavage signal, so the PFS handling stays Cas13-only (Cas13a's
   H preference, Cas13b's D-PFS), which is `Z0.1`.

## Contradiction named

The package's own field data (`rpa_asrpa_*` is RPA, not CRISPR; but the Cas12a
literature it draws on) presents the seed as the discriminating window, while the
2024 trans-cleavage measurement on the same enzyme shows singles across the whole
spacer within 1.4× of WT. Both are kept: the position rule as a **design
constraint**, the ratio as **evidence about magnitude**, and the report says which
is which.

## Why we are lenient here

One target gene (`CPSIT_0429`), one RPA-Cas12a format, endpoint fluorescence at
1 h — not a kinetic measurement, and the crRNA column mixes scaffold/loop/stem
engineering with a 3′ spacer truncation (`M65`, 18 nt). A hard magnitude gate
built on that would be a guess dressed as a rule. The gate therefore stays
positional and the magnitude is reported as unknown, which keeps real designs
alive without claiming evidence we do not have.

## Resulting `R/` change

* `R/marker_crispr.R`: `crispr_seed_evidence()` — the measured seed-mismatch
  signal window and the adjacent-pair trend as a small evidence object
  (`ratio_range`, `n_seed_singles`, `adjacent_better`, `source`), and the
  `crispr_seed_gate()` docstring restated as position-only.
* The gate itself is **unchanged** in behaviour; only its claim changes.

## Test that asserts it

`tests/testthat/test-crispr-seed.R`: `crispr_seed_gate()` is position-only
(inside/outside the seed, both effects); `crispr_seed_evidence()` matches the
extract (seed single-mismatch ratios for `M39` positions 1–8 span 0.373–0.968,
the unmodified-crRNA singles span 0.767–1.394, 17 of 19 adjacent doubles below
0.5); and the DAG's seed gate carries no magnitude claim.
