# HRM resolvability criterion — proposal (no code changes made)

> **DECIDED 2026-10-21** — see `hrm_resolvability.decision.md`: implemented as `HRM_DTM_RESOLVABLE` / `hrm_resolvability()` and the `family_hrm_resolvable` criterion row. The status
> line below is the state this file was written in; it is kept as the evidence
> record and is not edited after the decision.

Status: **proposal only**. `R/` is untouched. Every clause below cites a CSV row
that exists in `inst/extdata/` and a generator that re-derives it.

## What the evidence supports

**Clause 1 — the default resolvability test is a *difference of melting peaks*,
not a class label.** Use `dTm_obs` between the two allele homozygotes as the
quantity of interest.
Rows: `hrm_jsbbs2018_dtm_exp1.csv` (38 markers, 0.198–1.118 °C),
`hrm_jsbbs2018_dtm_exp3.csv` (18 markers, mean 0.74 °C),
`hrm_ijms2022_dtm_single_snp_controls.csv` (0–1.7 °C).

**Clause 2 — require `dTm >= 0.5 °C` as the *pass* threshold on the Roche
LightCycler 96.** This is the only numeric threshold any fetched source states
as its own empirical rule ("a ΔTm value >0.5 °C guaranteed accurate genotyping"),
and it is corroborated by the observed minimum: the smallest dTm the jsbbs2018
authors still called polymorphic was 0.198 °C *in a screened set*, while the
co-dominant Exp. 3 set (pre-filtered on |ΔΔG°3bp| > 1.43) has a 0.74 °C mean.
Rows: `hrm_jsbbs2018_selection_rule.csv` (amplicon-length row, basis column),
`hrm_jsbbs2018_summary.csv` (`dtm_exp1_min_c`).
*Caveat to ship with the threshold:* the source itself calls this
instrument-specific ("using the Roche LightCycler 96 system") and adds that some
class III/IV markers still genotyped correctly at 0.2–0.3 °C. So 0.5 °C is a
**conservative default for that platform**, not a physical constant.

**Clause 3 — the pre-design filter is a free-energy proxy, and it is a filter,
never a guarantee.** Use `|ΔΔG°3bp| > 1.0 kcal/mol` with amplicon < 80 bp
(the paper: this "would generally guarantee ΔTm of >0.5 °C"). Do **not** use
1.43 kcal/mol as a general threshold — that value was chosen for one experiment
and only 18/24 of those markers segregated co-dominantly.
Rows: `hrm_jsbbs2018_selection_rule.csv` (both `|ΔG3bp|` rows).

**Clause 4 — 3-allele / multi-locus panels must be reported as a *partition*,
not as a per-SNP flag.** With 3 C>T loci, HRM alone left 14/36 controls in 7
indistinguishable pairs (notably `CTC/TCT` vs `CCT/TTC`); all 22 others were
unique. A per-SNP resolvability flag cannot express this, so any multi-locus
input must yield the partition.
Rows: `hrm_ijms2022_partition_2snp.csv`, `hrm_ijms2022_unresolvable_3snp.csv`,
`hrm_ijms2022_discrimination.csv`.

**Clause 5 — nothing here is a general SNP-resolvability table.** All ijms2022
SNPs are class I (C>T) and all jsbbs2018 markers are soy/rice class I–IV
markers under one instrument. The widely cited cross-transversion HRM table
(Liew 2004: A::T↔G::C exchanges 0.8–1.4 °C, others <0.4 °C, ~4% unresolvable)
could **not** be fetched (closed access), so **no clause cites it** and the
package must not encode its numbers.
Block: `data-raw/proposed/readout.rows.csv` key `liew2004`.

## What this proposal deliberately does NOT add

- No per-SNP "resolvable: yes/no" column, because Clause 4 shows the outcome is
  a property of the panel, not of the SNP.
- No transversion penalty table (no fetched evidence).
- No dPCR or Sanger term in the HRM gate.
