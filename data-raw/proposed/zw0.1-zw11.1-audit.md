# ZW0.1 + Z11.1 — deliberate invariant exceptions and the user-facing decision audit

> **Superseded (2026-11-25).** The frame-invariance "known gap" recorded below
> is closed: see `docs/decision-audit.md` §«Frame invariance: closed by ZW0 M5»
> and `ROADMAP.md` №29 / ZW0 M5. This file stays as a provenance record; do not
> re-derive live numbers from it.

Proposal, 2026-11-03. Lands in `docs/decision-audit.md` (ZW0.1 §"Where the
invariant deliberately does not hold", Z11.1 §"Decision audit") plus one line in
the superseded inventory §W. Nothing here changes `R/`; the registry half
of ZW0.1 (`criterion_spec()$asymmetry_reason`) already exists.

## Where every number came from

Counted through `tidymarkers:::criterion_spec()` after `pkgload::load_all(".")`,
never by `grep` (the string is recycled by `criterion_tbl()`):

| column             | value                                                     | rows |
| ------------------ | --------------------------------------------------------- | ---- |
| (total)            | —                                                         | 193  |
| `source`           | `package assumption`                                      | 177  |
| `source`           | `kwok1990`                                                | 7    |
| `source`           | `santalucia1998`                                          | 9    |
| `evidence`         | `rule_of_thumb`                                           | 193  |
| `evidence`         | `vendor_guideline` / `primary_literature` / `fitted_here` | 0    |
| `asymmetry_reason` | non-`NA`                                                  | 9    |

The "176 of 193" printed in four vignettes and "191 rows" printed in ROADMAP
ZW1.3 are both stale. Naming a source raises `source`, not `evidence`: no row on
`primary_literature` appears in `source` because that column carries a key.

## ZW0.1 — the exceptions, and which test asserts each

Registry-declared (9, the rows `criterion_asymmetry()` names and
`test-invariants.R:190` enforces): `as:gate:gate_discrimination`,
`as:gate:gate_wobble`, `as:gate:gate_as_footprint`, `as:gate:gate_as_window`
(KASP FAM/HEX), `family_lamp_fip_alt`, `family_lamp_inner_tm_balance`,
`family_lamp_fip_variant`, `lamp:gate:gate_lamp_variant`,
`lamp:gate:gate_lamp_inner_tm` (LAMP FIP/BIP).

Outside the registry, each grounded in a named test:

| #   | where                   | what breaks                                                    | why                                                                                                               | test                                |
| --- | ----------------------- | -------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | ----------------------------------- |
| 1   | `test-invariants.R:181` | a revcomp locus is not representable                           | the locus model is strandless (`chrom`/`pos`/`ref`/`alt`)                                                         | asserts `!"strand" %in% names(snp)` |
| 2   | `test-invariants.R:124` | two criteria differ under a coordinate shift                   | they read the whole contig; excluded by label `in-silico PCR (REF 1 / ALT 0)`, `off-target exact primer hits (n)` | M3                                  |
| 3   | `test-invariants.R:33`  | REF/ALT swap changes the probe set                             | the AS primers swap _roles_, not sequences; 3' bases exchange FAM↔HEX                                             | M2                                  |
| 4   | `test-as-primer.R:1`    | the `-` strand AS 3' base is the complement                    | the primer binds the plus strand                                                                                  | makes the complement explicit       |
| 5   | `test-caps.R:214`       | every CAPS site must be its own reverse complement             | the scanner reads the plus strand only                                                                            | asserts `all(rc(site) == site)`     |
| 6   | `test-lamp.R:238`       | LAMP emits a paired REF set; an indel pair is `not_applicable` | REF template and failed reaction look alike; an indel has no single 3' base                                       | Z9.1                                |
| 7   | `test-decision.R:161`   | `NA` / `not_applicable` are never `pass`                       | "no evidence" is not "evidence in favour"; both stamp `unknown`/neutral                                           | asserts `unknown` and `neutral`     |

> **Superseded 2026-11-25:** closed by ZW0 M5 — see
> `docs/decision-audit.md` §«Frame invariance: closed by ZW0 M5» and ROADMAP №29.

Known gap, **not** an exception: `test-invariants.R:165`, the design search is
not frame-invariant (candidate ranking differs across frames; the evaluation row
must be the chosen candidate). It is a defect with a precondition guard.

## Z11.1 — what the four vignettes must say

Replace each `### Decision audit` body with a link to `docs/decision-audit.md`
and keep the gate list. The audit page carries:

- the graph is the scheme (`plot_decision_tree(scheme =)` without `assay`),
  every node `neutral`: "the graph passed" is not a claim the figure supports;
- the ladder `CRITERION_EVIDENCE` verbatim, and that all 193 rows sit on its
  weakest rung;
- 177/193 rows name no source at all; 7 `kwok1990` (3'-mismatch rank only, not
  the 0.5 discrimination bar) and 9 `santalucia1998` (nearest-neighbour ΔG) name
  one without moving the rung (ZW4.6);
- per scheme: rows whose only named source is `santalucia1998` are the
  structure checks, not the Tm/GC bars;
- `marker_decision_text(evidence = TRUE)` prints the rung and the weakest rung
  the scheme relies on — the per-run version of the same statement;
- fix the duplicated paragraph at `markers-crispr-sequencing.Rmd:203-211`.
