# Handoff — tidymarkers science worktree

## Where things stand

| Path | Branch | Role |
|---|---|---|
| `/root/development/tidymarkers` | `wip/roadmap-base` | main clone. **Uncommitted**: `ROADMAP.md`, `NEWS.md`, `.gitignore` (adds `.tmp/`) |
| `/root/development/tidymarkers-agent-science` | `feat/science-invariants` @ `200fd54` | science work (this) |
| `/root/development/tidymarkers-sci-{arms,hrm,tools}` | — | merged into `feat/science-invariants`; safe to remove |

Science worktree branches are fully merged. Nothing is on `master`.

## What landed (`feat/science-invariants`)

- **Merges** of `feat/sci-tools`, `feat/sci-hrm`, `feat/sci-arms` (one at a time).
  The first two merges of `scientific_backlog.csv` were **wrong** — one left a
  conflict marker and dropped the resolved rows. Fixed in `3dab306` by a row-wise
  3-way merge (`.tmp/merge-backlog.py` is the throwaway that did it): 23 rows, no
  markers, every advanced status preserved.
- **`data-raw/polymarker-source.R`** → `inst/extdata/polymarker_design_rules.csv`.
  Ten rules read from the PolyMarker source at VERSION 1.3.3 / `d4a4578`, each
  with method + line, file digests pinned. Settles Y3.2: allele on the AS
  3′ (`seq_original[pr.snp_pos]`), subgenome on the **common** primer
  (`right_pos => pos` → `SEQUENCE_FORCE_RIGHT_END`). The construct has no
  "homoeologous primer" field at all.
- **`data-raw/z7-priors.R`** → `z7_priors.csv` (8 predictions, contradictions
  kept), `z7_engine_probe.csv` (match vs **off**-allele; the blocker measured:
  terminal floor < 0.7%, all discrimination in the off allele).
- **`data-raw/science_predictions.R`** → `science_ledger.csv` (source → extract →
  test, 23 rows, 1 deliberate `blocked`).
- **Tests**: `test-z7-priors.R`, `test-z8-outcomes.R`, `test-y32-homeolog.R`,
  `test-science-ledger.R`. Plus vocabulary fixes in `test-scientific-data.R`.
  Science-only run (~1 min) is green: 8 files, 2 pre-existing documented skips.

## Corrections worth knowing

- **PolyMarker** is **`10.1111/pbi.12281`**, RNA-Seq BSA, PBJ 2015 (PMID
  25382230) — *not* the Bioinformatics DOI some notes cite. The
  `@misc{polymarker_source}` key keeps the source URL as its howpublished.
- **`mismatch_panel.csv`** cannot be regenerated: two different generator
  constructions in the tree, and neither reproduces the committed numbers. Not
  used by any test; `z7_engine_probe.csv` replaces it for Z7. Regenerating
  `data-raw/mismatch-discrimination.R` therefore rewrites that file — do not run
  it "to check".
- **`crisprBase`**: verified from the repo, not the docs. Cas12a ×3
  (`AsCas12a`, `enAsCas12a`, `MAD7`), Cas13 **only** Cas13d (`CasRx`). No PFS
  field. Do not promise it replaces the `cas` registry.
- **NIST STRBase** is now a client-rendered app; no static stutter table. Not
  verified.
- Do not trust the user's or the review text's paper identifiers either —
  `desilva2009` and `ace1r_aslamp2015` were both wrong in the backlog.

## Next

1. **Z7 itself** — the parameter. Must be built on `margin = match − off`
   (`z7_engine_probe.csv`), never on a single efficiency.
2. **`homeolog_diagnostics` / `allele_base` / `subgenome_anchor`** fields in the
   candidate; the rule table is already the source of truth.
3. Wire `engine_stamp` into `Evaluation` and the artifact manifest (stamp exists,
   wiring does not).
4. `ZW0` revcomp invariant is still untestable (strandless locus model) — that is
   the blocker to state, not to skip.
5. Review the ROADMAP/NEWS diff in the main clone and commit it there.

## Pitfalls (unchanged + new)

- `c()` **drops duplicate names** — that is how `has()`'s `quoted` column vanished.
- `utils::readLines` does not exist; `chartr("ACGT","TGCA", …)` maps A→T *and*
  T→A into T, giving `U` downstream. Use an explicit complement map.
- `unquoted commas in the backlog CSV` silently create rows; always check
  `anyNA(b$key)`.
- R regex is TRE: use `[[:space:]]+`, not `\s`.
- `DECIPHER::CalculateEfficiencyPCR` needs primer and target **the same length**
  and the target is **forward** (open with the complement of the primer's 3′ base).
