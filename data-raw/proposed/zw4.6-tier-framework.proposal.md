# ZW4.6 — the `*:tier:*` block: analytical framework, and the 55 non-tier rows

Status: proposal. No code or data changed. Date: 2026-11-05.

Scope: all 54 `criterion_spec()` rows with `role == "tier"` (every one carries
`limit = CRITERION_TIER_LIMIT`, "maximise (per-tier bottleneck min)"), plus the
55 non-tier rows still at `evidence == "rule_of_thumb"`. The question is how to
score a cascade of multi-parameter optimisations where some parts are hard gates
and others are weightable, and how those rows can honestly leave `rule_of_thumb`.

---

## 1. What the selector does today (verified in source)

Four functions decide everything; there is no hidden weighting layer.

| where                                          | behaviour                                                                                                                                                             |
| ---------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `R/marker_select.R:397` `order_by_margins()`   | `dplyr::arrange()` tier columns descending, then tie columns ascending. This is a **strict lexicographic order**.                                                     |
| `R/marker_select.R:404` `rank_candidates()`    | thin wrapper over `order_by_margins()`.                                                                                                                               |
| `R/marker_select.R:435` `select_with_ladder()` | gates first (all TRUE, AND); survivors ranked by `order_by_margins()`; if none survive, soft gates lowered least-grave-first as far as needed but never past a floor. |
| each `grade_*()`                               | builds margins (`m_* = observed - threshold`), then sets `tier_k` to a margin, sometimes `pmin()` of several.                                                         |

So the current objective is exactly **lexicographic max-min**:

    rank(c) = ( min_{j in T1} m_j(c),  min_{j in T2} m_j(c),  ..., ties )

- _within_ a tier, the aggregate is the bottleneck `min` (Rawlsian: the assay is
  as good as its weakest dimension);
- _across_ tiers, a strict priority order, so no gain in a later tier can pay
  for a loss in an earlier one;
- tie-breakers (`KASP_TIES_*`) only order otherwise-equal candidates.

Gate handling is an **ε-constraint** scheme in disguise: the gates are the hard
ε-bounds, and the relaxation ladder lowers one ε at a time in a declared order of
"least grave", stopping at a per-gate floor. The selector comment at
`R/marker_select.R:4-15` already states the intent ("some flaws are fatal, others
are just costly"). The implementation is coherent; what is missing is a declared,
testable _policy_ behind the tier order and the within-tier aggregate.

### Row shape

`R/marker_criteria.R:1135-1150` (`tier_rows()`) gives every tier row
`kind = "categorical"`, `bounds = NULL`, `values = CRITERION_TIER_LIMIT`,
`evidence = "rule_of_thumb"`. `R/marker_criteria.R:132-136`
(`CRITERION_TIER_EVIDENCE`) maps `moderate/weak/none` **down** the ladder and
never up. There is therefore no path by which any paper can raise a tier row:
the row's stated "limit" is the aggregation policy, and no paper prints
"maximise the per-tier bottleneck minimum".

**This is structural, not a sourcing failure.** No amount of literature reading
makes the _current_ tier row citable, because the row conflates two different
quantities (see §2). This is why the floor of 54 was reported earlier, and why a
different fix is required.

---

## 2. The key modelling error: one row, two quantities

A tier currently mixes:

- **(a) the reference point / target** — e.g. "amplicon length near the ideal
  50-100 bp", "probe Tm above primer Tm", "AS Tm near 63 °C", "hairpin ΔG ≥ -3
  kcal/mol", "≥ 10 bp band gap". These are quantities a paper or a vendor manual
  _does_ print, so they are sourceable.
- **(b) the aggregation policy** — the bottleneck `min`, the tier membership,
  and the priority order. These are package decisions no paper prints.

Because (a) and (b) share one `limit` string, the row is pinned at
`rule_of_thumb` even when (a) is fully sourced. The fix is to **split them**:
keep the printed target as a first-class, citable field, and declare the
aggregation policy once, in code, as an auditable design decision.

This split is what converts the user's observation — "the papers rank or predict,
they do not measure this package's threshold" — into a rule we can live with: the
_target_ can be `primary_literature`/`vendor_guideline`; the _policy_ stays
`rule_of_thumb` **by construction** and is documented, not disguised as a
threshold.

---

## 3. State-of-the-art for this exact shape of problem

The problem shape is: hard feasibility constraints + a non-compensatory ranking
over several normalised dimensions, where some dimensions are clearly more
important than others. This is textbook multi-criteria decision making (MCDM).
The relevant SOTA building blocks:

| method                                           | what it gives us                                                                                                                                      | canonical reference                                                                                                                                                      | adopt?                                                                     |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------- |
| Constraint satisfaction / constraint programming | the hard-gate block: prune the candidate grid before scoring, propagate feasibility                                                                   | standard CP; Rossi, van Beek & Walsh, _Handbook of Constraint Programming_ (2006)                                                                                        | **yes**, as Layer 0                                                        |
| ε-constraint method                              | formalises the relaxation ladder: optimise the ranking subject to `gate ≤ ε`, lower ε stepwise                                                        | Haimes, Lasdon & Wismer, _IEEE Trans. SMC_ SMC-1(3):296-297 (1971), DOI 10.1109/TSMC.1971.4308298                                                                        | **yes**, it is already what the ladder does                                |
| Lexicographic optimisation                       | strict priority between bands; no cross-band compensation                                                                                             | Fishburn, _Mathematics of Decision Theory_ (1972); Ehrgott, _Multicriteria Optimization_ (2005)                                                                          | **yes**, as Layer 2                                                        |
| Pareto dominance / non-dominated sorting         | keep the genuine trade-off front instead of hiding it; stop a single winner from implying false certainty                                             | Deb, Pratap, Agarwal & Meyarivan, _IEEE Trans. Evol. Comput._ 6(2):182-197 (2002), DOI 10.1109/4235.996017                                                               | **yes**, as Layer 4 (report top-k)                                         |
| Weighted Chebyshev / Tchebycheff scalarisation   | non-compensatory aggregation with weights: minimise the worst weighted regret; the bottleneck still dominates but dimensions can differ in importance | Steuer & Choo, _Math. Oper. Res._ 8(2):252-262 (1983), DOI 10.1287/moor.8.2.252; Miettinen, _Nonlinear Multiobjective Optimization_ (1999)                               | **yes**, as Layer 3                                                        |
| Weighted sum                                     | the usual default; fully compensatory and scale-sensitive                                                                                             | Miettinen (1999), above                                                                                                                                                  | **no** — it can let a great lower-priority dimension mask a near-fatal one |
| AHP (analytic hierarchy process)                 | derives defensible priority weights from structured pairwise comparisons                                                                              | Saaty, _Eur. J. Oper. Res._ 45(1):9-23 (1990), DOI 10.1016/0377-2217(90)90057-I                                                                                          | **yes**, for eliciting band priorities                                     |
| PROMETHEE / ELECTRE (outranking)                 | handles indifference/preference thresholds and mixed hard/soft criteria without forcing commensurability                                              | Brans & Vincke, _Manage. Sci._ 31(6):647-656 (1985), DOI 10.1287/mnsc.31.6.647; Roy (1968)                                                                               | alternative, only if indifference thresholds are wanted                    |
| Max-min fairness / price of fairness             | formalises "as good as its weakest part" and quantifies what the bottleneck costs                                                                     | Bertsimas, Farias & Trichakis, _Oper. Res._ 59(1):17-31 (2011), DOI 10.1287/opre.1100.0865                                                                               | **yes**, justifies Layer 3                                                 |
| Robust optimisation                              | makes the winner stable to margin estimation error (Tm, ΔG, dropout are model outputs)                                                                | Ben-Tal & Nemirovski, _SIAM J. Optim._ 13(3):769-805 (2002), DOI 10.1137/S1052623499363220; Bertsimas & Sim, _Oper. Res._ 52(1):35-53 (2004), DOI 10.1287/opre.1030.0065 | **yes**, as a robustness margin                                            |
| Successive halving / Hyperband                   | cheap pruning of a large candidate grid before the expensive score                                                                                    | Li, Jamieson, DeSalvo, Rostamizadeh & Talwalkar, _Evol. Comput._ 26(1):129-149 (2018), DOI 10.1162/evco_a_00201                                                          | only if the grid explodes; not needed for correctness                      |

The combination that matches this package's stated intent — _gate the grave,
score the minor, never compensate a fatal flaw_ — is:

> **constraint filtering → lexicographic priority bands → weighted Chebyshev
> inside a band → Pareto front reporting → calibratable weights.**

We call it **banded non-compensatory compromise (BNCC)**.

---

## 4. Proposed framework: BNCC, layer by layer

Each layer maps to an existing touchpoint, so the change is local.

### Layer 0 — Feasibility (unchanged semantics, optional speed-up)

Hard gates stay an AND over boolean predicates. Formally this is the feasible
set `F = {c : g_i(c) = TRUE ∀ hard i}`. Optional: prune the candidate grid by
constraint propagation before computing any margin (performance only; no change
to which candidate wins).

The relaxation ladder stays, restated as ε-constraint: when `F = ∅`, lower one
soft gate `g_k` from its declared ε, in the declared least-grave order, only as
far as needed, never past `floor_k`. This is already `select_with_ladder()`.

### Layer 1 — Normalise margins to regret

Margins have mixed units (kcal/mol, °C, bp, counts). Convert every tier
dimension `j` to a **regret** in `[0, 1]` on the current candidate set:

    r_j(c) = (z*_j - z_j(c)) / (z*_j - z^n_j)

with `z*_j` the ideal (the printed target) and `z^n_j` the nadir on `F`. This
makes dimensions commensurable without discarding their physical meaning. Where
no nadir is available (tiny candidate set), use the declared window width as the
denominator, which is already what the margins encode.

### Layer 2 — Priority bands (lexicographic across bands)

Group tiers into a small number of **bands** with a strict order. The principle:

- **Band A — validity**: would the assay work at all if this is bad?
  structural specificity, discrimination, Tm-vs-target agreement.
- **Band B — robustness**: will it work reliably across the sample?
  Tm balance, amplicon window, probe/Tm ordering.
- **Band C — operational**: cost and convenience.
  footprint load, 3'-end weakness, variant centrality, off-target count.

Candidates are compared lexicographically **by band score**, so a candidate that
wins Band A always beats one that loses it. This preserves the current
non-compensation guarantee but makes the guarantee _explicit and auditable_,
instead of implicit in an arbitrary 5-deep tier list.

### Layer 3 — Within a band: weighted Chebyshev

Inside a band, score a candidate by

    S_B(c) = max_{j in B}  w_j · r_j(c)      (minimise)

The Chebyshev (min-max) metric is the canonical non-compensatory scalarisation:
it minimises the worst weighted regret. The _bottleneck still dominates_ (so the
"as good as its weakest part" property survives), but weights `w_j` now express
that, e.g., a 1 °C Tm imbalance matters more than 10 bp of amplicon length.

Backward compatibility: if each band contains exactly one tier, all `w_j = 1`
and all ranges are equal, then `argmin max_j r_j = argmax min_j m_j`. The
proposal is therefore a **strict generalisation** of the current lexicographic
max-min, not a replacement: the current behaviour is the diagonal case.

### Layer 4 — Report the Pareto front

Compute the non-dominated set of `F` over the band scores. The winner is the
BNCC argmin; but if `k > 1` alternatives lie on the front within a declared
tolerance, report them (top-k), ordered by the existing tie-breakers. This turns
"we picked one" into "we picked one and here is the trade-off", which is exactly
what the selector header promises ("we never hide it"). It also gives the
calibration stage (§6) its natural response variable.

### Layer 5 — Calibration hook

The only things left to fit are the **band boundaries** and the **weights** `w_j`
(and, if wanted, the robustness radius). They are the parameters that make a row
`fitted_here`. Everything above is deterministic given them.

---

## 5. Recommended priority order per stage

The current per-scheme tier orders mostly already follow the band principle.
The audit below flags the exceptions. "Order OK" means the declared order is
already validity → robustness → operational.

| stage        | current tiers (order)                                             | band reading                                                  | flag                                                              |
| ------------ | ----------------------------------------------------------------- | ------------------------------------------------------------- | ----------------------------------------------------------------- |
| `as`         | discrimination, Tm, 3' weakness, footprint                        | A: discrimination, Tm; C: 3' weakness, footprint              | OK                                                                |
| `pair`       | structure, Tm, amplicon, footprint, COM 3' weak                   | A: structure, Tm; B: amplicon; C: footprint, 3' weak          | OK                                                                |
| `duplex`     | self-complementarity, amplicon budget                             | A: selfcomp; B: amplicon                                      | OK                                                                |
| `asq`        | AS Tm toward 63 °C optimum                                        | A: Tm                                                         | OK (single tier)                                                  |
| `caps`       | band gap, Tm balance, amplicon                                    | A: band gap; B: Tm balance, amplicon                          | OK                                                                |
| `arms`       | inner discrimination, band separation, tube Tm                    | A: discrimination; B: band sep, tube Tm                       | OK                                                                |
| `hrm`        | variant centrality, Tm balance                                    | A: Tm balance; B: centrality                                  | **flip** — Tm balance is validity, centrality is operational      |
| `castpcr`    | 3' specificity, Tm balance, amplicon                              | A: specificity; B: Tm balance, amplicon                       | OK                                                                |
| `crispr`     | guide quality, seed, specificity, amplicon                        | A: guide, seed, specificity; B: amplicon                      | OK                                                                |
| `dpcr`       | probe Tm margin, probe self-structure, amplicon                   | A: probe Tm, self-structure; B: amplicon                      | OK                                                                |
| `indel_size` | Tm balance, amplicon, 3' end                                      | A: Tm balance; B: amplicon; C: 3' end                         | OK                                                                |
| `lamp`       | discrimination, inner Tm set, structure, amplicon                 | A: discrimination, structure; B: inner Tm, amplicon           | **minor**: structure before inner Tm is defensible but declare it |
| `lax`        | resolvable ALT haplotype freq                                     | A: discrimination                                             | OK (single tier)                                                  |
| `rpa`        | discrimination placement, primer length, probe geometry, amplicon | A: discrimination, probe geometry; B: primer length, amplicon | OK                                                                |
| `sanger`     | variant centrality, Tm balance, amplicon                          | A: Tm balance; B: centrality, amplicon                        | **flip** — Tm balance before centrality                           |
| `snapshot`   | Tm centrality, self-structure, amplicon                           | A: Tm centrality, self-structure; B: amplicon                 | OK                                                                |
| `ssr`        | Tm balance, amplicon, 3' end                                      | A: Tm balance; B: amplicon; C: 3' end                         | OK                                                                |
| `taqman`     | probe discrimination, Tm ordering, amplicon                       | A: discrimination, Tm ordering; B: amplicon                   | OK                                                                |

The two flips (`hrm`, `sanger`) are the only order changes the framework
recommends; everything else is an explicit declaration of what is already there.

---

## 6. What makes a tier row citable, and what makes it `fitted_here`

### Citable (target split, §2)

Once the target is a separate field, a tier row's target can carry the rung of
the source that prints the target. This does **not** require any new literature
for most tiers: the registry already holds the anchors (§7). The aggregation
policy is recorded once, in code and in a design note, as a package decision —
honest, auditable, and not pretending to be a threshold.

### `fitted_here` (calibration protocol)

The weights and band boundaries are fit by:

1. **Factors**: the tier margins (`m_*`) of candidate assays actually designed.
2. **Response**: a wet-lab outcome — binary call success, cluster separation, or
   ΔCt between alleles, sampled per candidate.
3. **Design**: a D-optimal or response-surface design over the feasible margin
   space (Box & Wilson 1951 lineage), with replicates; never a convenience
   sample of only successes.
4. **Model**: monotone ordinal/logistic regression of the outcome on the
   normalised regrets; the coefficients are the weights `w_j`, the intercepts the
   thresholds. Bootstrap CIs.
5. **Artefact**: store the dataset hash, the design, the fit and the CIs; add a
   `science_ledger.csv` row with `use="report"` and a `tests/testthat` pin. Only
   then does `evidence` move to `fitted_here`.

Until step 5 exists, no tier row can honestly be `fitted_here`; the framework
above is what makes step 5 well-posed.

---

## 7. Tier targets that are already anchorable from the collected registry

**Status 2026-11-06 — phases A and B executed.** The split is in
`R/marker_criteria.R`: every tier row carries a `target` (the printed reference
point) and keeps `limit` as the aggregation policy; `CRITERION_TIER_POLICY_EVIDENCE`
declares that policy `rule_of_thumb` once. 36 of 54 tiers were raised after the
quote was re-grepped from the cached full texts (18 remain `rule_of_thumb`).
`rule_of_thumb` 109 -> **73**, `primary_literature` 91 -> **126**,
`vendor_guideline` 12 -> **13**. One new source was needed,
`kalendar2022_competitive` (DOI 10.3389/fmolb.2022.773956), for the 63 C ASQ
optimum. The table below is the original work list; the rows not raised are
listed under "Not anchorable this way".

These are **candidate** anchors: each must be re-grepped in `/tmp/lit` and the
target sentence quoted before any row is raised. Rows are marked by confidence:
**A** = a gate already cites this key for the same quantity; **B** = the source
prints the reference window/target; **C** = likely, needs a read.

| tier rows                                                                                                                                  | target (the `z*`)                                  | candidate key                        | DOI                                                    | conf |
| ------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------- | ------------------------------------ | ------------------------------------------------------ | ---- |
| `pair:tier:tier1`; `dpcr:tier:tier2`; `snapshot:tier:tier2`                                                                                | hairpin/homodimer ΔG ≥ −3 kcal/mol                 | `degenprime2024`                     | 10.1093/bioadv/vbae044                                 | A    |
| `pair:tier:tier1`                                                                                                                          | dimer/cross-dimer ΔG                               | `shimizu2025_kasp`                   | 10.1270/jsbbs.24052                                    | A    |
| `pair:tier:tier2`                                                                                                                          | COM primer Tm +2 °C above AS                       | `easykasp2025`                       | 10.1186/s12859-025-06322-x                             | A    |
| `pair:tier:tier3`; `caps:tier:tier3`; `castpcr:tier:tier3`                                                                                 | amplicon 50-100 bp (KASP); 120-1200 bp (CAPS)      | `easykasp2025` / `liu2016_caps`      | 10.1186/s12859-025-06322-x / 10.1270/jsbbs.66.244      | A    |
| `as:tier:tier2`; `caps:tier:tier2`; `castpcr:tier:tier2`; `hrm:tier:tier2`; `indel_size:tier:tier1`; `ssr:tier:tier1`; `sanger:tier:tier2` | primer Tm window / Tm balance                      | `vavilov2024` / `easykasp2025`       | 10.18699/vjgb-24-40 / 10.1186/s12859-025-06322-x       | B    |
| `duplex:tier:tier1`; `duplex:tier:tier2`                                                                                                   | inter-pair complementarity < 6 bp; selfcomp ≤ 4 bp | `wang2025_multikasp`                 | 10.1186/s12870-025-07005-8                             | A    |
| `asq:tier:tier1`                                                                                                                           | AS-primer Tm optimum (63 °C)                       | `fret_asqpcr2021`                    | 10.3389/fpls.2021.747886                               | B    |
| `caps:tier:tier1`; `arms:tier:tier2`                                                                                                       | gel band gap ≥ 10 bp                               | `ssr_gel2026`                        | 10.3389/fpls.2026.1873607                              | A    |
| `hrm:tier:tier1`                                                                                                                           | variant centrality / melt-domain placement         | `jsbbs2018`                          | 10.1270/jsbbs.18048                                    | C    |
| `hrm:tier:tier2`                                                                                                                           | measured dTm                                       | `gundry2008`                         | 10.1093/nar/gkn204                                     | B    |
| `crispr:tier:guide`, `:seed`, `:specificity`, `:amplicon`                                                                                  | PAM/PFS, spacer length, seed window                | `cas12a_beech2026` / `cas13_pfs2026` | 10.1038/s42003-026-10805-9 / 10.1093/nar/gkag616       | A/B  |
| `dpcr:tier:tier1`                                                                                                                          | probe Tm above primers (68-70 °C)                  | `cov_a4_qpcr2026`                    | 10.3390/v18101120                                      | A    |
| `dpcr:tier:tier3`                                                                                                                          | amplicon 70-150 bp                                 | `ppr_ddpcr2026`                      | 10.3390/vetsci13080841                                 | A    |
| `indel_size:tier:tier2`; `:tier:tier3`                                                                                                     | amplicon 200-320 bp; length mutation               | `indel_melon2025`                    | 10.3390/plants14243782                                 | A    |
| `lamp:tier:tier2`; `:tier:tier3`; `:tier:tier4`                                                                                            | F2/B2 60-65 °C; no structure; target ≤ 300 nt      | `nemeth2025`                         | 10.3389/fpls.2025.1568657                              | A    |
| `lamp:tier:tier1`                                                                                                                          | mismatch position / Tt delay                       | `mismatch_rtlamp2019`                | 10.3389/fmicb.2019.01056                               | B    |
| `rpa:tier:tier1`                                                                                                                           | 3'-anchor mismatch class                           | `rpa_mismatch_class`                 | 10.1016/j.jmoldx.2022.08.005                           | A    |
| `rpa:tier:tier2`                                                                                                                           | minimum primer length 30 nt                        | `riordan2022`                        | 10.1186/s13007-022-00923-w                             | A    |
| `snapshot:tier:tier1`                                                                                                                      | extension-primer Tm window 53-60 °C                | `zhao2026_snapshot`                  | 10.1186/s12864-026-13173-0                             | A    |
| `snapshot:tier:tier1`                                                                                                                      | 3' end one base upstream of variant                | `tso2022_snapshot`                   | 10.1186/s13007-022-00923-w                             | A    |
| `ssr:tier:tier1`; `:tier:tier2`; `:tier:tier3`                                                                                             | amplicon/product sizes; repeat lengths             | `ssr_gel2026` / `moso_bamboo2026`    | 10.3389/fpls.2026.1873607 / 10.1186/s12870-026-09292-1 | A    |
| `taqman:tier:tier1`                                                                                                                        | probe GC/design; mismatch ΔG                       | `vavilov2024`                        | 10.18699/vjgb-24-40                                    | A    |
| `taqman:tier:tier2`                                                                                                                        | probe Tm 68-70 °C, above primers                   | `cov_a4_qpcr2026`                    | 10.3390/v18101120                                      | A    |
| `taqman:tier:tier3`                                                                                                                        | amplicon 70-150 bp                                 | `ppr_ddpcr2026`                      | 10.3390/vetsci13080841                                 | A    |

**Not anchorable this way (keep at `rule_of_thumb`, by design):**

- `as:tier:tier1` (discrimination) — the tier maximises a _computed_ ratio; no
  paper prints a minimum. Gate already covers the floor.
- `as:tier:tier3`; `arms:tier:tier1`; `castpcr:tier:tier1` — Kwok-style
  mismatch _ranking_ (user decision: ranking ≠ threshold).
- `pair:tier:tier4`; `as:tier:tier4` — footprint load is a package simulation.
- `pair:tier:tier5`; `indel_size:tier:tier3`; `ssr:tier:tier3` — 3'-end
  weakness is the README/manual qualitative rule.
- `lax:tier:tier1` — resolvable ALT haplotype frequency is a package metric.
- `sanger:tier:tier1`; `sanger:tier:tier3` (centrality, read length) — needs a
  dedicated Sanger read-length source.
- `arms:tier:tier3` (tube Tm balance) — vendor/manual range, needs reading.

**Net effect on the metric.** About 40 of the 54 tier rows have an existing
anchor for their target. With the §2 split, those targets move off
`rule_of_thumb`; only the ~14 rows above stay, and they stay _honestly_. This is
the concrete route from 109 toward the low 70s the user asked for — and it does
not require a single new source for most of it, only reading what is already in
the registry.

---

## 8. The 55 non-tier rows: plan by group

Counts are exact from `criterion_spec()` (`role != "tier"`, `evidence ==
"rule_of_thumb"`). Actions assume the no-fabrication rules already in force.

### G1 — Structural / bookkeeping (10 rows)

`site_class`, `site_alt_alleles`, `site_linked_skipped`, `site_alt_covered`,
`site_linked_dependent`, `as_n_3p`, `as_n_var`, `com_n_var`, `family_single_base`,
`family_sanger_pair_defined`.

These are counts or definitions, not external thresholds ("exactly one ALT",
"zero variants under a primer"). Two honest options:

1. keep `rule_of_thumb` and record in `criterion_provenance.csv` that the bar is
   **structural** (a design definition), or
2. add a distinct provenance kind (`structural`) rather than abusing an evidence
   rung.
   Recommendation: option 1 now (no schema churn), option 2 next major version.

### G2 — Simulator / computed outputs (10 rows)

`as_dropout`, `as_foot_dropout`, `com_dropout`, `discrimination`,
`in_silico_pcr`, `family_lamp_inner_tm_balance`, `family_crispr_offtarget`,
`family_size_fragment`, `family_probe_asp_discrimination`,
`family_probe_extension`.

The _model_ is citation-backed (`santalucia1998`, OligoArrayAux, DECIPHER), but
the cutoff is a package choice. Plan: cite the model (not the cutoff), document
the cutoff as a package assumption, and target `fitted_here` via the §6 protocol
once outcome data exist. `santalucia1998` (10.1073/pnas.95.4.1460) anchors the
thermodynamics; it does not print these cutoffs.

### G3 — Kwok ranking (6 rows)

`as:gate:gate_wobble`, `arms:gate:gate_arms_disc`, `as_mismatch_type`,
`castpcr:gate:gate_castpcr_as3`, `castpcr:gate:gate_castpcr_mismatch`,
`castpcr:gate:gate_castpcr_blocker`.

Already decided: keep `rule_of_thumb`; the source ranks, it does not measure this
package's bar. Do not spend effort here. (`kwok1990`, DOI
10.1093/nar/18.4.999; `newton1989`, DOI 10.1093/nar/17.7.2503; `liu2012_arms`,
DOI 10.1186/1746-4811-8-34.)

### G4 — Chemistry geometry / family rules (11 rows)

`family_caps_digest_spans`, `family_caps_site_conserved`, `family_rpa_thf_off`,
`family_rpa_probe_spans`, `family_rpa_variant_under`, `family_crispr_seed`,
`family_crispr_guide_spans`, `family_crispr_conserved`,
`family_sanger_primers_inside`, `family_sbe_primer_ratio_reported`,
`family_probe_extension`.

These state structural requirements that the defining assay paper or the
registry already prints ("probe spans the variant", "guide/PAM conserved", "THF
off the variant", "extension 3' base at the variant"). Plan: one read per row and
raise where the defining paper states the requirement. Likely anchors:
`tso2022_snapshot`/`atanasovska2012`/`fanis2014` (SBE, DOIs
10.1186/s13007-022-00923-w / 10.1371/journal.pone.0048167 /
10.1186/1471-2164-15-108), `cas12a_beech2026` (CRISPR), `snap2caps` (CAPS
digest, DOI 10.1093/nar/gnh006), `lamp_comp2018b`/`nemeth2025` (LAMP).

### G5 — Numeric bars with no anchor yet (18 rows)

`as:gate:gate_discrimination`, `as:gate:gate_as_footprint`,
`as:gate:gate_as_window`, `pair:gate:gate_com_footprint`,
`arms:gate:gate_arms_outer`, `crispr:gate:gate_crispr_pam`,
`crispr:gate:gate_crispr_seed`, `crispr:gate:gate_crispr_specificity`,
`dpcr:gate:gate_dpcr_probe_specific`, `lamp:gate:gate_lamp_outer_tm`,
`lamp:gate:gate_lamp_cross_dimer`, `lax:gate:gate_lax_coverage`,
`lax:gate:gate_lax_alt_specific`, `lax:gate:gate_lax_cross_priming`,
`sanger:gate:gate_sanger_amplicon`, `sanger:gate:gate_sanger_unique`,
`snapshot:gate:gate_snapshot_unique`, `snapshot:gate:gate_snapshot_amplicon`,
`taqman:gate:gate_taqman_discrimination`.

Each needs a targeted literature pass with the existing tooling, in this order:

1. **Vendor/manual windows** first (`lamp_outer_tm`, `sanger_amplicon`,
   `snapshot_amplicon`, `arms_outer`) — most likely `vendor_guideline`.
2. **Primary paper structural bars** next (`crispr_pam`, `crispr_seed`,
   `lax_*`, `dpcr_probe_specific`).
3. **Physical bounds** last (`sanger_unique`, `snapshot_unique`,
   `lamp_cross_dimer`, `as_footprint`, `as_window`, `com_footprint`): if only a
   physical/mechanical constraint exists, record it as such with a rationale, not
   as a fabricated citation.
   `family_sbe_primer_ratio_reported` (a `rule` row) is in G4 by content; it has a
   printed ratio in the SNaPshot sources and should be raised there.

**Realistic ceiling for the 55.** G4 (~11) and the vendor-first half of G5 (~8)
are the plausible raises; G1/G2/G3 (~26) should stay or move only via
`fitted_here`. So, counting the whole `rule_of_thumb` block (109 = 54 tiers +
55 non-tier):

- tier targets raised (§7): −40 → 14 tiers remain;
- non-tier raised (G4 + vendor-first G5): −19 → 36 remain (26 of them
  structurally stuck in G1/G2/G3).

Honest floor after Phases A–B+E: **14 + 36 = 50 rows**; if only the tier split
lands first, 69. The gap from 50 to 70 must come from `fitted_here`, not from
more citation hunting.

---

## 9. Recommended execution order

1. **Phase A — schema split (§2).** Add the target/reference field and a single
   declared aggregation policy. No rung moves yet; add tests for the split.
   This is the prerequisite for every raise.
2. **Phase B — source the tier targets (§7).** Re-grep each candidate anchor in
   `/tmp/lit`, quote the printed target, raise only where the quote holds.
   Expected: ~40 rows off `rule_of_thumb`.
3. **Phase C — swap to BNCC (§4) behind a flag.** Keep the lexicographic path as
   the default until the band definitions and weights are pinned by tests;
   add the Pareto top-k report. Order-equivalence test for the diagonal case is
   the acceptance gate.
4. **Phase D — calibration (§6).** Once wet-lab outcomes exist, fit weights and
   band boundaries; move the first rows to `fitted_here` with a ledger row.
5. **Phase E — G4/G5 literature pass** for the non-tier rows, vendor-first.

Phases A and B are independent of C/D and deliver the metric movement the user
asked for. Phase C is the correctness improvement; D is the only route to
`fitted_here`.

---

## 10. Risks and what this does not do

- **The split changes the public criterion schema.** `marker_criteria()`'s
  `limit` column is a rendering of `bounds`; adding a target field touches the
  registry, the plotter and the text export. That is a real change, not a
  cosmetic one, and needs its own tests.
- **BNCC is order-equivalent to the current selector only in the diagonal case.**
  With weights ≠ 1 it _will_ reorder candidates. That is the point, but it must
  be behind a flag and pinned by tests before it becomes the default.
- **No tier row becomes `primary_literature`.** The aggregation policy never
  will; only the target does. Stating otherwise would be the fabrication the
  project forbids.
- **~26 non-tier rows are structurally stuck.** More searching will not move
  them; only `fitted_here` will. The 70 target is reachable only with calibration
  data, not with citations alone.
- This document changes nothing. It is a plan; the code is untouched.
