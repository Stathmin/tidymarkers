# Decision audit: where the numbers come from, and where the invariant does not hold

Two questions a user asks before trusting a verdict: what is each threshold based
on, and which symmetries was it allowed not to have. The first is the section
“The evidence behind every threshold”, the second is the section “Where the
invariant deliberately does not hold”.

## The evidence behind every threshold

Every criterion row carries two different things: where the number comes from
(`source`) and how strongly the number is fixed (`evidence`, the ladder
`CRITERION_EVIDENCE`). Counted through `criterion_spec()`, never by `grep` (the
string is recycled by `criterion_tbl()`):

| column             | value                                         | rows of 212 |
| ------------------ | --------------------------------------------- | ----------- |
| `source`           | `package assumption` — no printed source      | 63          |
| `source`           | a named key of the source registry (132 rows) | 149         |
| `evidence`         | `rule_of_thumb` (weakest rung)                | 71          |
| `evidence`         | `vendor_guideline`                            | 13          |
| `evidence`         | `primary_literature`                          | 128         |
| `evidence`         | `fitted_here`                                 | 0           |
| `asymmetry_reason` | non-`NA` (a declared exception, below)        | 9           |

### Reading the flags

- **`source = package assumption`** — the bar is what this package chose. It is
  not wrong, but it is not cited, and no paper should be implied. 63 of 212 rows
  are in this class.
- **A `*:tier:*` row carries two quantities, and only one is citable.** The
  `target` column is the reference point the tier ranks toward (an amplicon
  window, a dG floor, a printed design rule); the `limit` column is the
  aggregation policy (bottleneck `min`, tier membership, order), a package
  decision no paper prints. A tier's `evidence` is the rung of its target; the
  policy stays `rule_of_thumb` by construction
  (`CRITERION_TIER_POLICY_EVIDENCE`). 36 of 54 tiers now name a source for their
  target; the other 18 rank a computed quantity (Kwok ranking, a simulator, a
  population metric) and stay `rule_of_thumb`.
- **The aggregation policy is declared once, and there are two of them.** The
  default is the lexicographic tier order; an opt-in BNCC policy
  (`R/marker_bncc.R`) normalises margins to a regret, compares declared priority
  bands lexicographically, scores a band by its weighted Chebyshev (min-max)
  regret and reports the Pareto front. Both are package policies: neither is a
  printed threshold, so neither raises a rung, and `BNCC_CALIBRATION` stays
  `fitted = FALSE` until wet-lab outcomes exist (proposal §6).
- **`source = kwok1990`** — the 3'-terminal mismatch _ranking_ is Kwok's (Kwok et
  al., NAR 1990, doi 10.1093/nar/18.4.999). It does not state the
  `discrimination >= 0.5` or the 0.1 dropout number; it ranks mismatch classes by
  yield, so these rows stay `rule_of_thumb`.
- **The ΔG numbers still use SantaLucia 1998's nearest-neighbour parameters**
  (doi 10.1073/pnas.95.4.1460); those are a **model**, not a measurement. Round 4
  found papers that print the _threshold_ the model output has to clear —
  DeGenPrime (a combined hairpin/dimer dG below -3 kcal/mol), Vavilov (primers
  and probes at ΔG >= -9 kcal/mol) and Shimizu (self-/cross-dimer below -10
  kcal/mol triggers redesign) — so the threshold rows rise while the parameter
  source does not.
- **Naming a source does not raise the rung.** A paper that only ranks or
  predicts is not enough; the same-quantity rule (below) requires a printed
  threshold. `fitted_here` is still 0, so nothing in this package is calibrated
  to its own measurements. Check the bars on your own instrument before a
  genotype call.
- The one legacy string `primary literature` on `as_mismatch_type` is not a rung
  name; that row names `kwok1990` in the source key map.

### What a figure proves

The decision graphs (`decision-{scheme}.svg` / `.txt`) are the **scheme**, built
by `plot_decision_tree(scheme = <id>)` without an `assay`. Every node is
`neutral`: the figure shows the order and hardness of the gates, and no candidate
walk with statuses. `marker_figure()` panels show the construct and the checks,
not clusters, droplets, electropherograms or ΔCt. "The graph passed" and "the
figure is green" are not claims these figures support.

### The same statement per run

`marker_decision_text(evidence = TRUE)` appends each check's rung, and the
weakest rung the scheme relies on, to the text graph. That is the per-assay
version of the table above: most rows still print `rule_of_thumb`, and the
SNaPshot rows print the rung their printed source fixes.

## Where the invariant deliberately does not hold

The metamorphic tests ask whether a verdict survives a change that cannot change
the answer: M1 reverse-complement the locus, M2 swap REF/ALT, M3 shift the
coordinate frame, M4 rename the chromosome. For a few criteria the answer is
_no_, and the exception is structural, not a bug.

### Declared in the registry (ZW0.1)

`tidymarkers:::criterion_spec()` carries one column, `asymmetry_reason`: `NA` for
a symmetric row, a one-sentence reason for an asymmetric one. Nine rows carry it.
`tests/testthat/test-invariants.R` ("ZW0.1") enforces that the set is exactly
these nine, that each names a real criterion id once, that no reason is empty,
and that the column is appended last.

<!-- prettier-ignore -->
| criterion id                   | what breaks                                               | why                                                                                                                                                                                                                      |
| ------------------------------ | --------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `as:gate:gate_discrimination`  | a REF/ALT swap moves an allele to another dye             | KASP: `AS_ref` carries the FAM tail and `AS_alt` the HEX tail, at the same AS concentration. The gates grade the chosen allele-specific primer, so the swapped probe set is a different set, not the same one relabelled |
| `as:gate:gate_wobble`          | the wobble gate is not mirrorable across REF/ALT          | it grades the chosen AS primer once, while that primer's allele is pinned to FAM or HEX by its tail                                                                                                                      |
| `as:gate:gate_as_footprint`    | a swap re-assigns the graded primer to the other dye slot | the footprint is judged on `AS_ref` / `AS_alt` as placed on the FAM/HEX tails                                                                                                                                            |
| `as:gate:gate_as_window`       | the AS window reports the primer in its fixed dye slot    | same reason: the REF/ALT image is a different probe set                                                                                                                                                                  |
| `family_lamp_fip_alt`          | only FIP reads the variant                                | the ALT allele sits at the 3' end of FIP's displacement arm; BIP does not read it                                                                                                                                        |
| `family_lamp_inner_tm_balance` | FIP and BIP Tm patterns are constrained separately        | each inner primer's 5' loop arm must be hotter than its 3' arm, and BIP is placed downstream of FIP                                                                                                                      |
| `family_lamp_fip_variant`      | same as `family_lamp_fip_alt`                             | the variant-reading role belongs to FIP alone                                                                                                                                                                            |
| `lamp:gate:gate_lamp_variant`  | the variant gate grades FIP only                          | same reason                                                                                                                                                                                                              |
| `lamp:gate:gate_lamp_inner_tm` | the inner Tm gate is not symmetric in FIP and BIP         | same reason                                                                                                                                                                                                              |

### Structural, outside the registry

Each row is asserted by a named test, not inferred.

<!-- prettier-ignore -->
| where                     | what breaks                                                                                                                        | why                                                                                                                                                                       | test                                                                                        |
| ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| the locus model           | a reverse-complemented locus is not representable                                                                                  | `locus_at()` has `chrom` / `pos` / `ref` / `alt` and no strand column, so a revcomp locus is only another REF/ALT pair; M1 is not testable until the model grows `strand` | `test-invariants.R` asserts `!"strand" %in% names(snp)`                                     |
| two contig-wide criteria  | "in-silico PCR (REF 1 / ALT 0)" and "off-target exact primer hits (n)" differ under a coordinate shift                             | they read the whole contig, and a shifted frame is a different contig; M3 excludes them **by label**, not by weakening the assertion                                      | `test-invariants.R`, M3                                                                     |
| REF/ALT swap (M2)         | the AS primers swap _roles_, not sequences, and their 3' bases exchange between FAM and HEX                                        | the discriminating base must follow the allele, not the primer slot                                                                                                       | `test-invariants.R`, M2                                                                     |
| the `-` strand AS primer  | its 3' base is the **complement** of the allele                                                                                    | the primer binds the plus strand, so using the allele base itself would invert discrimination                                                                             | `test-as-primer.R`                                                                          |
| the CAPS restriction scan | every catalogue site must be its own reverse complement                                                                            | the scanner reads the plus strand only, so an asymmetric cutter would be modelled wrongly                                                                                 | `test-caps.R`                                                                               |
| a LAMP design             | it carries a paired REF set (`$paired`, `$ref_fip`, `$discriminating_base`); an indel gives `status = "not_applicable"` and a note | a REF-only template and a failed reaction both give no amplification, so an ALT-only set cannot tell them apart; an indel has no single 3' base to flip                   | `test-lamp.R` (Z9.1)                                                                        |
| status, every scheme      | `NA` and `not_applicable` never become `pass`                                                                                      | "no evidence" is not "evidence in favour": an unknown row stamps `unknown`, and a non-applicable one is neutral, not ok                                                   | `test-decision.R`; user-facing table in [reading results](../vignettes/reading-results.Rmd) |

### Known gap, not an exception

The design search is **not** frame-invariant: candidate ranking differs across
frames, and the evaluation row must be the chosen candidate. That is a defect,
kept as a precondition-guarded test (`test-invariants.R`, "ZW0 KNOWN GAP"), not
an accepted asymmetry. Do not add it to `asymmetry_reason`.
