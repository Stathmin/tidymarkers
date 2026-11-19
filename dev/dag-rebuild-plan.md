# DAG rebuild plan — 20 marker chemistries, 10 agents

This is the orchestrator's brief. It exists so ten agents can each fix one chunk of the
decision DAGs (the trees behind `vignettes/figures/decision-<scheme>.svg`) without
colliding: how they **search**, how they **report**, how they **edit**, and how they
**test only their own chunk**.

Source of the "scientific requirements" used throughout: `README.md` §3 (_Краткий
мануал по подбору KASP_ — LGC/KBioscience + EasyKASP 2025), the lecture deck
`dev/lecture/kasp-slides.qmd`, and each scheme file's own `# Design rules (researched, not
invented)` header comment. Agents must treat those as the truth and cite them.

---

## 1. What a "DAG" is here (read this first)

One decision DAG per marker chemistry. It is _generated_, never hand-drawn:

```
marker_criteria(scheme)      R/marker_criteria.R   → 5 site checks + 21 assay checks + the stage's gates/tiers
marker_decision_stages(scheme) R/marker_decision.R → the scheme's selection stages
      │
      ▼
marker_decision_schema(scheme)  R/marker_decision.R  → nodes + edges (order, x/y, hard)
      ├─ marker_decision_graph()  → tidygraph, stamps status from assay$selection / $evaluation
      ├─ plot_decision_tree()   → the PNG (dev/render-decision.R)
      └─ marker_decision_text()   → the greppable text export
```

The scheme's _stage_ (its gates / relaxation ladder / tiers) comes from either the
built-in constants in `R/marker_select.R` (`KASP_GATES_*`, `KASP_RELAX_*`, `KASP_TIERS_*`)
or from a discovered `scheme_def_<scheme>()` in `R/marker_<scheme>.R`. The **checks** are
the fixed lists in `R/marker_criteria.R`, tagged with schemes by `check_schemes()` in
`R/marker_schemes.R`. Acceptance is run by `marker_design_checks()` in `R/marker_checks.R`.

Inspecting any DAG locally (needs a working R env, see §10):

```r
pkgload::load_all(".")
cat(marker_decision_text("lamp"), sep = "\n")   # the whole DAG as text
marker_decision_schema("lamp")$nodes            # the node table
```

---

## 2. The 20 schemes and their files

`marker_schemes()` has 9 built-in rows + 11 discovered (`scheme_def_*`). Every scheme has
one test file. "Kind" matters for the write scopes in §5.

| scheme       | family     | kind       | source of the stage                               | test                                  |
| ------------ | ---------- | ---------- | ------------------------------------------------- | ------------------------------------- |
| `primer`     | primer     | built-in   | none (marker-agnostic)                            | `test-decision.R`, `test-criteria.R`  |
| `snp`        | kasp       | built-in   | `KASP_*_AS` + `_PAIR`                             | `test-selection.R`, `test-decision.R` |
| `indel`      | kasp       | built-in   | same as `snp`                                     | `test-indel.R`                        |
| `multi`      | kasp       | built-in   | same as `snp`                                     | `test-mnp.R`, `test-alleles.R`        |
| `duplex`     | multiplex  | built-in   | `KASP_*_DUPLEX`                                   | `test-duplex.R`                       |
| `asq`        | asq        | built-in   | `KASP_*_ASQ`                                      | `test-decision.R`                     |
| `caps`       | caps       | built-in   | `KASP_*_CAPS`                                     | `test-caps.R`                         |
| `arms`       | arms       | built-in   | `KASP_*_ARMS`                                     | `test-arms.R`                         |
| `hrm`        | hrm        | built-in   | `KASP_*_HRM`                                      | `test-hrm.R`                          |
| `castpcr`    | castpcr    | discovered | `scheme_def_castpcr` · `R/marker_castpcr.R`       | `test-castpcr.R`                      |
| `crispr`     | crispr     | discovered | `scheme_def_crispr` · `R/marker_crispr.R`         | `test-crispr.R`                       |
| `dpcr`       | dpcr       | discovered | `scheme_def_dpcr` · `R/marker_dpcr.R`             | `test-dpcr.R`                         |
| `indel_size` | indel_size | discovered | `scheme_def_indel_size` · `R/marker_indel_size.R` | `test-indel-size.R`                   |
| `lamp`       | lamp       | discovered | `scheme_def_lamp` · `R/marker_lamp.R`             | `test-lamp.R`                         |
| `lax`        | kasp       | discovered | `scheme_def_lax` · `R/marker_lax.R`               | `test-lax.R`                          |
| `rpa`        | rpa        | discovered | `scheme_def_rpa` · `R/marker_rpa.R`               | `test-rpa.R`                          |
| `sanger`     | sanger     | discovered | `scheme_def_sanger` · `R/marker_sanger.R`         | `test-sanger.R`                       |
| `snapshot`   | snapshot   | discovered | `scheme_def_snapshot` · `R/marker_snapshot.R`     | `test-snapshot.R`                     |
| `ssr`        | ssr        | discovered | `scheme_def_ssr` · `R/marker_ssr.R`               | `test-ssr.R`                          |
| `taqman`     | taqman     | discovered | `scheme_def_taqman` · `R/marker_taqman.R`         | `test-taqman.R`                       |

---

## 3. The contract — invariants every DAG must satisfy

These are the acceptance criteria for the whole effort. Agent 1 (foundation A) turns
them into a reusable, machine-checked helper (§7) so each scheme agent tests them for
its own scheme only. **No agent may relax an invariant to make its scheme pass.**

- **C1 — one criterion, one node.** A physical criterion is either a selection _gate_
  or an acceptance _check_, never both. Today gates restate checks (e.g. `gate_discrimination`
  ↔ `check:discrimination`; `gate_tm` ↔ `check:com_tm_as_tm`; `gate_amplicon` ↔
  `check:amplicon_length_bp`), sometimes one gate ↔ two checks (footprints). Kill the
  duplication: the graph shows the criterion once.
- **C2 — check membership is family-based.** A scheme carries only the checks that its
  chemistry obeys. Non-KASP chemistries must not carry AS/COM/discrimination/3'-window
  checks; the 5 _site_ checks (`target site class`, `ALT alleles (n)`, … `design relies
on linkage`) are KASP-selector notions and apply only where an ALT primer/linkage exists.
- **C3 — graph == evaluation.** `marker_design_checks(assay)` evaluates exactly the checks
  the scheme's DAG displays (today it ignores `scheme` entirely, `R/marker_checks.R:113`).
- **C4 — hard failure dooms.** A failed _hard_ gate must have a direct path to `DOOMED`.
  A mixed hard/soft stage must not draw a hard gate as merely "relaxable" (today only
  all-hard stages emit `decision --fail--> doomed`).
- **C5 — a relaxed win still passes its check.** If the ladder relaxes a gate past the
  acceptance limit, the winner fails the corresponding check (e.g. amplicon floor `-20`
  admits 120 bp, but the check demands `[50,100]`). Floors and `KASP_LIMITS` must be
  reconciled so "relaxed" and "accepted" never contradict.
- **C6 — two declared orders.** Gate order = the order the method applies them (grave /
  cheap first); relaxation order = _least grave first_. Both are declared per scheme and
  labelled. They need not be equal.
- **C7 — no decorative gates.** Every gate the DAG shows is actually evaluated by the
  scheme's designer/selector. No declaration-only gates (today all `gate_caps_*`,
  `gate_arms_*`, `gate_asq_*`, `gate_duplex_*`, `gate_crispr_*`, … are never emitted).
- **C8 — hardness is justified.** `hard = TRUE` only when the method says the assay is
  _meaningless_ without it (and give the `gate_hard_label`); otherwise soft, with a
  bounded, reachable floor.
- **C9 — no placeholders.** Limits are real, parameterised values (`KASP_LIMITS`),
  never a literal `[lo, hi]` (see `gate_caps_tm`).
- **C10 — ungated list is honest.** `KASP_UNGATED_CHECKS` must not call a criterion
  ungated while a scheme hard-gates it (off-target is hard in `ssr`/`indel_size`).
- **C11 — acyclic, strictly downward, reachable.** Every edge endpoints exists; `ord[to] >
ord[from]`; `assay:verdict` present (already tested — keep it).
- **C12 — the runtime graph matches the scheme.** `assemble_kasp_assay()` stores its
  scheme and builds _that_ scheme's `$decision` (today always `snp`), and
  `plot_decision_tree(assay)` honours it.
- **C13 — reachable or explicitly documentation-only.** A scheme is reachable through
  `marker_design_at()` / the ladder, or it is marked documentation-only with a stated
  reason (`asq`, `primer`). Same for status stamping (gate statuses are stamped only for
  `as`/`pair` today).
- **C14 — tiers are non-vacuous.** A tier must not merely restate a hard gate already
  enforced (e.g. `sanger` tier1 "unique amplicon" while uniqueness is a hard gate).

---

## 4. Phase 0 — foundation (Agents 1–2, serial, before anything else)

Everything below depends on these two, and they must land first so the scheme agents
have a stable, shared contract and non-overlapping files.

### Agent 1 — foundation A (taxonomy, checks, harness)

Files: `R/marker_schemes.R`, `R/marker_criteria.R`, `R/marker_checks.R`,
`R/marker_scheme_registry.R`, `tests/testthat/test-criteria.R`,
`tests/testthat/test-decision.R` (only the shared/generic assertions),
new `tests/testthat/helper-dag.R`, new `dev/dump-dags.R`, `NEWS.md`, `ROADMAP.md`.

1. Define the family → check-set taxonomy (C2): which of the 5 site + 21 assay checks
   each `family` serves. Rewrite `check_schemes()` so a discovered scheme gets its
   family's set, not "everything".
2. Make `marker_design_checks()` scheme-aware (C3): read the checker's `scheme` and filter
   `marker_criteria(scheme = )`, computing only the rows that scheme owns. Keep the
   REF/ALT/footprint _values_ only where the chemistry has them.
3. Introduce the DAG-invariant helper `expect_dag_invariants(scheme)` in `helper-dag.R`
   covering C1–C11 as far as statically checkable, so scheme agents call it in one line.
4. `dev/dump-dags.R`: write `marker_decision_text()` for **all** schemes to a file, so an
   agent can diff its own scheme before/after.
5. Rewrite the shared tests to iterate over `marker_schemes()` reading each scheme's spec
   (scheme-agnostic), so scheme agents never edit them.

### Agent 2 — foundation B (DAG generator, graph, floors)

Files: `R/marker_decision.R`, `R/marker_select.R` (generic selector + shared floors only),
`R/marker_assemble.R` (only the `$decision` wiring), `R/marker_scheme_registry.R`,
`tests/testthat/test-selection.R`, `tests/testthat/test-decision.R` (generic),
plus the mechanical extraction in step 3.

1. Fix the generator: (a) remove the gate↔check duplicate nodes (C1); (b) route a failed
   hard gate straight to `DOOMED` even in a mixed stage (C4); (c) parameterise the
   `[lo, hi]` placeholder (C9); (d) make the site/gate/report order explicit (C6).
2. Reconcile the AS/PAIR floors with `KASP_LIMITS` (C5) so a relaxed winner still passes
   its check (or the check reads the relaxed bar) — today `gate_amplicon` floor `-20`
   and `gate_structure` floor `-3` both admit candidates the checks reject.
3. Make all 20 schemes uniform so the write scopes are disjoint: move the built-in stage
   constants (`KASP_*_DUPLEX/ASQ/CAPS/ARMS/HRM`) out of `R/marker_select.R` into
   `scheme_def_duplex/asq/caps/arms/hrm()` in their own files, leaving `marker_select.R`
   with only the generic selector + the shared AS/PAIR constants.
4. Scheme-aware graph (C12/C13): store `scheme` on the assay object and build/stamp its
   `$decision` graph; make `plot_decision_tree(assay)` honour it.
5. Update `KASP_UNGATED_CHECKS` to respect scheme-level off-target gates (C10).

**Gate:** after Phase 0, `Rscript dev/check.R --no-render` must pass and every scheme's
`marker_decision_text()` must be free of duplicate-criterion pairs, `[lo, hi]`, and
hard-gate-reachable-only-via-refine. Freeze §3 as the contract reference.

### 4.5 Seams and interfaces (verified from the code graph, 2026-11-06)

Read this before assigning work. It lists the boundaries Phase-1 agents must _not_ move,
where each boundary lives today, and which briefs in §4/§8 are already stale.

**Method.** `graphify-out/` was rebuilt at HEAD `fd0780c` (7809 nodes, 10406 edges, 411
communities) and the health gate is green: `Rscript dev/graphify-health.R --check`.
One caution: in this graph the `calls` relation is **intra-file only** — 0 of 2586 `calls`
edges cross a file — so it cannot answer "who calls this in another file". Cross-file
consumers below were read from parsed source (`parse()` + `all.names()`), which is the same
static view `codetools` uses. The graph contributes what it is good at: ownership
communities (each seam file is its own community, so file-scoped write scopes are real)
and call fan-in (god nodes).

#### The nine seams

A seam is a function one owner changes and many others call. Its signature and return
shape are frozen by Phase 0; a scheme agent may _use_ it, never change it.

| id  | seam              | signature and return                                                                                                                                                                                                         | defined (owner)                                                                                                                                                  | consumers                                                                                                                                                                                                                                                                 | contract                                                                                                           |
| --- | ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| S1  | criteria registry | `marker_criteria(limits = KASP_LIMITS, scheme = NULL)` → `marker_criteria` tibble, one row per criterion (`id`, `label`, `role`, `limit`, …)                                                                                 | `R/marker_criteria.R:1444` (A1)                                                                                                                                  | `marker_design_checks` (`R/marker_checks.R:277`), `marker_generic_checks` (`:731`), `KASP_UNGATED_CHECKS` (`R/marker_criteria.R:39`), `marker_decision_schema` (`R/marker_decision.R:84`), `marker_decision_text` (`:448`)                                                | one criterion id → one row; `role %in% c("gate","check")`; `scheme=` returns only that family's checks (C1, C2)    |
| S2  | stage spec        | `marker_decision_stages(scheme = "snp")` → named list of `list(id, label, gates, relax, tiers)`                                                                                                                              | `R/marker_decision.R:39`; the built-in table is `R/marker_stages.R:37` (A2)                                                                                      | `marker_decision_schema`, `marker_decision_text`, `tests/testthat/helper-dag.R:26`                                                                                                                                                                                        | `gates`/`relax`/`tiers` are gate _ids_ resolvable through `scheme_label()`; `relax` is least-grave-first (C6)      |
| S3  | DAG schema        | `marker_decision_schema(scheme = "snp")` → `list(nodes = tibble(name, stage, kind, label, order, x, y, hard, detail), edges = tibble(from, to, kind))`                                                                       | `R/marker_decision.R:84` (A2)                                                                                                                                    | `marker_decision_text`, `marker_decision_graph`, `marker_scheme_adapter` (`R/marker_scheme_registry.R:220`)                                                                                                                                                               | endpoints exist, `ord[to] > ord[from]`, one `assay:verdict` (C11)                                                  |
| S4  | graph object      | `marker_decision_graph(assay = NULL, scheme = NULL)` → the plottable graph                                                                                                                                                   | `R/marker_decision.R:289` (A2)                                                                                                                                   | `assemble_kasp_assay` (`R/marker_assemble.R:1024`, `:1410`), `assemble_primer_assay`, `plot_decision_tree` (`R/marker_decision.R:401`), `marker_scheme_adapter`                                                                                                           | the graph is built for `assay$scheme`, never a default (C12)                                                       |
| S5  | check evaluation  | `marker_design_checks(assay, limits = KASP_LIMITS, scheme = NULL)` → check tibble (`criterion`, `status`, `value`, …)                                                                                                        | `R/marker_checks.R:277` (A1)                                                                                                                                     | widest fan-in: `assemble_kasp_assay`, `assemble_primer_assay`, `evaluation_verdict`, `marker_decision_graph`, `marker_evaluate`, `passes_all`, `design_kasp_markers`, `plot_checks`, `marker_risk`, `design_search_frames`, `R/marker_duplex.R`, `R/marker_search_anim.R` | evaluated check set == `marker_decision_schema(scheme)` check labels (C3)                                          |
| S6  | scheme definition | `scheme_def_<scheme>()` → `list(scheme, family, label, discrimination, readout, zygosity, amplicon, loci, verdict, gates, relax, tiers, gate_short, gate_limit, gate_floor_label, gate_hard_label, analytical, stage_label)` | one per scheme file (the scheme agent); discovered by `extra_scheme_defs()` (`R/marker_scheme_registry.R:17`), designer resolved by `scheme_designer()` (`:172`) | `extra_scheme_rows`, `extra_stage_rows`, `extra_decision_stages`, `marker_strategy`, `marker_target`, `marker_performance`                                                                                                                                                | scheme ids unique; every gate id has an entry in each `gate_*` map; `gates ⊆` that scheme's criteria ids           |
| S7  | selection ladder  | `select_with_ladder(data, gates, relax_order, tiers, ties, …)` → ranked selection with the relax trail                                                                                                                       | `R/marker_select.R:441` (A2 generic, A3 AS/PAIR)                                                                                                                 | `assemble_kasp_assay`, `R/marker_dpcr.R`, `R/marker_duplex.R`, `R/marker_lax.R`, `R/marker_sanger.R`                                                                                                                                                                      | a relaxed winner must still pass the check it traded against (C5); `warn_relaxed()` (`:499`) is the reporting side |
| S8  | status stamping   | `marker_check_status()`, `marker_evaluation_table()`, `marker_with_population_checks()`                                                                                                                                      | `R/marker_checks.R:30`, `:65`, `:586` (A1)                                                                                                                       | `R/marker_duplex.R:290`, `R/marker_hrm.R`, `assemble_kasp_assay`, `evaluation_verdict`                                                                                                                                                                                    | every stage with gates is stamped, not only `as`/`pair` (C13)                                                      |
| S9  | test harness      | `expect_dag_invariants(scheme)` → invisible `scheme`, asserts the statically visible invariants                                                                                                                              | `tests/testthat/helper-dag.R:56` (A1)                                                                                                                            | one call in each scheme's test file (20 today)                                                                                                                                                                                                                            | must stay genome-free and assay-free, so it runs in the inner loop                                                 |

**S5 is the load-bearing seam.** It is where the graph and the evaluator must agree, and
it has the largest fan-in. Freeze its return shape first; everything else can follow.

#### Write-scope collisions that must be resolved before Phase 1

The write scopes in §4/§5 are file-disjoint for Agents 3–10, but three overlaps remain _inside_
Phase 0:

1. `R/marker_scheme_registry.R` is listed under **both** Agent 1 (§4 Agent 1, line 134) and
   Agent 2 (§4 Agent 2, line 153). Give it to Agent 2 (it owns `extra_decision_stages()`),
   or split the file; do not let both edit it.
2. `tests/testthat/test-decision.R` is listed under both as well ("shared/generic assertions"
   in A1, "generic" in A2). Split by ownership in the file header, or hand the whole file to
   A1 and have A2 add cases only after A1 lands.
3. Section splits within one file: `R/marker_select.R` (A2 generic + floors, A3 AS/PAIR) and
   `R/marker_assemble.R` (A2 `$decision` wiring, A3 the file). Serial order already prevents a
   write race, but A2 must leave a banner-delimited section boundary so A3 can rebase on it.

#### Status delta — §4 and §8 have been overtaken by the tree

Measured on the current tree (not from the brief): all 20 schemes pass
`expect_dag_invariants()`, no DAG has a duplicate criterion label or a `[lo, hi]`
placeholder, and `caps`/`arms`/`hrm` emit 4/3/2 gate nodes (so their gates are no longer
declarative-only). Already done relative to this document:

| brief item          | evidence in the tree                                                                                                                                                                                                                                           |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| C3, Agent 1 task 2  | `marker_design_checks(assay, limits, scheme)` — `R/marker_checks.R:277`                                                                                                                                                                                        |
| C10, Agent 2 task 5 | `KASP_UNGATED_CHECKS(scheme)` is scheme-aware — `R/marker_criteria.R:39`                                                                                                                                                                                       |
| C12, Agent 2 task 4 | `marker_decision_graph(..., scheme = scheme)` — `R/marker_assemble.R:1024`, `:1410`                                                                                                                                                                            |
| Agent 2 task 3      | per-scheme stage constants already live in their own files (`R/marker_caps.R:42`, `R/marker_arms.R:38`, `R/marker_hrm.R:50`, `R/marker_duplex.R:49`, `R/marker_sanger.R:93`, `R/marker_asq.R:139`), and `R/marker_stages.R` is the single built-in stage table |
| C1, C7, C9          | no duplicate labels, every scheme's gates are emitted, no placeholder limits                                                                                                                                                                                   |

The static helper did **not** cover C5 (relaxed winner still passes), C8 (hardness
justified), C13 (reachability and status stamping beyond `as`/`pair`) or C14 (non-vacuous
tiers) when this was written. All four are enforced now: commit `4be0b7c` widened
`tests/testthat/helper-dag.R` (C8 at `:223`, C14 at `:289`) and added the C5 row to
`test-selection.R` and the C13 loop to `test-decision.R:244`; see “Status — §8 closed”
below. Do not hand §8's defect list to a new agent as if it were current.

---

## 5. Phase 1 — per-scheme agents (Agents 3–10, parallel, disjoint files)

Each agent owns exactly its scheme files + its test file(s). **Hard rule:** if an agent
believes it must edit a file outside its scope (especially `marker_select.R`,
`marker_criteria.R`, `marker_decision.R`, `marker_checks.R`, `marker_primer_pairs.R`), it must
stop and request a serial slot from the orchestrator, not edit it inline and collide.

| Agent | Schemes                           | Owned files (write scope)                                                                                                                                                                             |
| ----- | --------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 3     | `snp`, `indel`, `multi`, `primer` | `R/marker_assemble.R`, `R/marker_select.R` (AS/PAIR section), `R/marker_design.R`, `R/marker_target.R` (core path), `tests/testthat/test-selection.R`, `test-indel.R`, `test-mnp.R`, `test-alleles.R` |
| 4     | `duplex`, `asq`                   | `R/marker_duplex.R`, new `R/marker_asq.R`, `tests/testthat/test-duplex.R`, new `test-asq.R`                                                                                                           |
| 5     | `lax`                             | `R/marker_lax.R`, `tests/testthat/test-lax.R`                                                                                                                                                         |
| 6     | `caps`, `arms`                    | `R/marker_caps.R`, `R/marker_arms.R`, `tests/testthat/test-caps.R`, `test-arms.R`                                                                                                                     |
| 7     | `hrm`, `castpcr`                  | `R/marker_hrm.R`, `R/marker_castpcr.R`, `tests/testthat/test-hrm.R`, `test-castpcr.R`                                                                                                                 |
| 8     | `crispr`, `rpa`, `lamp`           | `R/marker_crispr.R`, `R/marker_rpa.R`, `R/marker_lamp.R`, `tests/testthat/test-crispr.R`, `test-rpa.R`, `test-lamp.R`                                                                                 |
| 9     | `dpcr`, `taqman`, `snapshot`      | `R/marker_dpcr.R`, `R/marker_taqman.R`, `R/marker_snapshot.R`, `tests/testthat/test-dpcr.R`, `test-taqman.R`, `test-snapshot.R`                                                                       |
| 10    | `sanger`, `ssr`, `indel_size`     | `R/marker_sanger.R`, `R/marker_ssr.R`, `R/marker_indel_size.R`, `tests/testthat/test-sanger.R`, `test-ssr.R`, `test-indel-size.R`                                                                     |

Agent 3 also updates the committed README figures (`man/figures/decision-tree.svg` for
`snp`, `decision-multiplex.svg` for `duplex` via `Rscript dev/render-decision.R`). The
per-chemistry SVGs and their `.txt` twins are committed under `vignettes/figures` and
embedded by the vignettes; only the PNG rasters there are git-ignored.

---

## 6. What each Phase-1 agent must produce

### 6a. Search (evidence first)

1. Read the scheme's own header (`# Design rules (researched, not invented)` block) in
   `R/marker_<scheme>.R` — it cites the method literature; that is the authority.
2. Read `README.md` §3 and, where relevant, `dev/lecture/kasp-slides.qmd`.
3. Read the scheme's `scheme_def_<scheme>()` (gates/relax/tiers/labels) and its designer.
4. Dump the current DAG: `cat(marker_decision_text("<scheme>"), sep = "\n")` (or
   `dev/dump-dags.R`).
5. `grep` for the scheme's gate ids across `R/` to confirm which are actually emitted by
   the engine (C7) — e.g. the audit found `gate_caps_*`/`gate_arms_*` are never emitted.

### 6b. Report (fill this template _before_ editing — this is the deliverable)

```
## <scheme>
Method basis: <where the requirements come from; cite file:line + external refs>
DAG now: <ordered node list; flag duplicates, decorative gates, placeholder limits>
Defects: <numbered; each with file:line, the invariant it violates (C…), and the fix>
Proposed DAG: gate order (grave/cheap first) | hard/soft per gate | relax order + floors
              | tiers | site checks | assay checks (family-correct set)
Biological risk if unfixed: <what unviable result the current DAG returns>
Tests to add/extend: <file + what they assert; must call expect_dag_invariants()>
```

### 6c. Edit

- Fix the designer and/or `scheme_def_<scheme>()` so the _declared_ DAG is real (C7) and
  method-true. Implement the proposed gates/tiers/hard-soft/order/floors.
- Add only the checks your chemistry obeys (via the family taxonomy from Phase 0).
- Do **not** silence a check or widen a limit just to pass — fix the cause.
- Keep tidyverse style: `|>`, `.data$`, `snake_case`; comments explain _why_.
- Add a `NEWS.md` bullet and a `ROADMAP.md` line for each scheme you touch.

### 6d. Test **only your chunk**

- Inner loop (fast, your files only): the suite is one file per scheme, so
    ```sh
    Rscript -e 'testthat::test_local(filter = "crispr|rpa|lamp")'
    ```
    runs only `test-crispr.R`, `test-rpa.R`, `test-lamp.R`.
- Add to each of your test files a call to the shared invariant helper:
    ```r
    test_that("<scheme> DAG obeys the contract", expect_dag_invariants("<scheme>"))
    ```
    plus explicit assertions for _your_ ordering, hardness and gate set (Phase 0 made the
    shared tests scheme-agnostic, so you never edit them).
- Before declaring done, also run the two cheap shared DAG tests (they read every
  scheme; do not edit them):
    ```sh
    Rscript -e 'testthat::test_local(filter = "decision|criteria|selection")'
    ```
- Full gate is the orchestrator's, not yours: `Rscript dev/check.R --no-render`
  (adds styler, lintr, roxygen sync, install smoke, figure sync). Expect to fix the
  figure step if you own `snp`/`duplex`.

---

## 7. Test harness (built once in Phase 0, used by everyone)

`tests/testthat/helper-dag.R` exposes:

```r
expect_dag_invariants(scheme)  # C1 dup-criteria, C2 check membership, C4 hard→doomed,
                               # C6 orders declared, C7 gates emitted, C9 no placeholder,
                               # C10 ungated honesty, C11 acyclic/downward/verdict
```

`dev/dump-dags.R` writes every scheme's `marker_decision_text()` to `results/dags.txt`
(git-ignored), so an agent can `diff` its scheme's DAG before/after and paste the delta
into its report.

---

## 8. Per-scheme defect briefs (grounded; each agent expands its own)

These were found by static audit (read-only; the code was not changed). They are a
starting point, not gospel: **run each scheme and confirm** before editing (`cat(marker_decision_text("<scheme>"), sep = "\n")`, and `dev/dump-dags.R`). Line numbers are
from the current tree — re-verify.

**Core (Agent 3).** `snp`: amplicon floor `-20` admits 120 bp but the check demands
`[50,100]`; structure floor `-3` admits hairpin `-6` the check rejects (C5); relax order
puts `wobble` before footprint though the manual calls wobble the _worst_ mismatch (C6);
tier3 rewards _more_ stable 3' ends, against the KASP destabilising rationale; no COM
3'-stability criterion. `indel`: DAG is byte-identical to `snp` although it declares
"allele-specific amplicon length" — **there is no REF/ALT size-resolvability gate**
(cf. `gate_caps_resolvable`, `gate_arms_bands`); `gate_wobble` is vacuous (`wobble = NA`).
`multi`: `wobble`/`disc` computed on `alts[1]` only, ignoring the degenerate IUPAC 3'
base; no check that the 3' code is the consensus of _all_ ALTs. `primer`: unreachable
(`marker_design_at()` has no branch; absent from the ladder) — mark documentation-only (C13)
or wire it; missing an in-silico-product-count check.

**`duplex`/`asq` (Agent 4).** `duplex`: relax node advertised but there is no
`KASP_GATE_FLOOR_DUPLEX` and the designer never calls `select_with_ladder`; a failed
_hard_ `gate_duplex_complement` still returns `"warn"` not DOOMED (C4); amplicon gate
promises ≤150 while the check demands `[50,100]` (C5). `asq`: no engine at all
(`R/marker_asq.R` missing, not in `per_locus_schemes()`); both gates
(`gate_asq_cassette`, `gate_asq_readout`) are panel/format properties, not
candidate-computable, yet hard (C7/C8); "Kalendar 2022 / 4-plex cassette" citation
unverified.

**`lax` (Agent 5).** `family = "kasp"`, so it carries the _whole_ as+pair machinery plus
its own stage; no hard gate and the coverage floor literally says "leave a haplotype
uncovered (warn)" — defeats the scheme's promise (C8); no LD/pop-freq guard that baked-in
linked variants are ALT-specific (a common REF variant ⇒ dropout); the pooled tube is
screened for dimers only, not cross-priming.

**`caps`/`arms` (Agent 6).** Both stages are **never evaluated** — the designers
(`design_caps_marker`, `design_arms_marker`) never call `grade_*`/`select_with_ladder`
and never emit `gate_caps_*`/`gate_arms_*` (C7). `gate_caps_tm` limit is the literal
`[lo, hi]` (C9); `gate_caps_amplicon` `[120,1200]` vs `check:amplicon_length_bp` `[50,100]`
(C5); relaxing `gate_caps_resolvable` risks an unresolvable digest despite being the whole
point. `arms`: `gate_arms_bands` floor equals its threshold (can never relax, drawn soft);
`gate_arms_outer` has no engine support. Both silently _lose_ the site + 9 allele-specific
checks (the mirror of C2) and their AS/COM-labelled checks don't fit a flanking pair.

**`hrm`/`castpcr` (Agent 7).** `hrm`: gate set matches practice but hard/soft is
debatable (amplicon, HRM's core requirement, is the relaxable one); the "variant under
neither primer" rule is structural, not a gate. `castpcr`: **targets the wrong allele** —
the designer picks `best$ref_primer` (REF 3' base) although castPCR (and its own header:
Roma 2013, Barbano 2015) detects the low-abundance ALT; `zygosity="co-dominant"`
unsupported by a single AS primer; probe is locus-specific, not allele-specific; the MGB
blocker underpinning the 0.1–1% LOD is unmodelled; `gate_castpcr_probe_tm` hard/soft
contradicts the physically identical `gate_dpcr_probe_tm`.

**`crispr`/`rpa`/`lamp` (Agent 8).** `crispr`: the single always-present **hard**
`gate_crispr_pam` is Cas12a-specific, yet PAM-free Cas13 is supported — a Cas13 marker is
misrepresented (C8); no guide-specificity gate; pre-amp screened with KASP PCR Tm/GC;
`discrimination`/`REF1/ALT0` checks are actively wrong (pre-amp amplifies both alleles).
`rpa`: the probe gate is unconditional though the exo probe is optional; the _governing_
primer length is soft while the designer hard-requires 30–35 nt (contradiction); the
no-probe path accepts the variant anywhere under the primer, not at the 3' end. `lamp`:
the designer forces the variant into an unprimed loop (inner arms stop short), so no
primer 3' end reads it — the marker cannot discriminate REF/ALT, contradicting the
header's "ends of the inner primers"; loop primers LF/LB absent despite "4–6 primers"; no
cross-dimer gate.

**`dpcr`/`taqman`/`snapshot` (Agent 9).** `dpcr`: the probe is real but no probe
Tm/GC/structure/off-target is in the acceptance set; inherits AS-primer checks that don't
apply; `readout` claims Poisson/partitioning never computed. `taqman`: `gate_taqman_probe_tm`
is declared **hard** while the header explicitly says the probe-Tm rule is _not_ a hard
gate (tier2) — internal contradiction; the gate is never evaluated (only used as a sort
key); the real probe-discrimination (mismatch ΔG) is neither gate nor check. `snapshot`:
`gate_snapshot_tm` says `[50,55]` but the designer relaxes to `[48,60]`; no uniqueness gate
for the extension primer; relax order not "least grave first"; tier1 restates a hard gate.

**`sanger`/`ssr`/`indel_size` (Agent 10).** `sanger`: no structure/GC/`>2` identical-3'-base
gate despite the header; designs no sequencing primer; tier1 "unique amplicon" is vacuous
(uniqueness is hard-gated); `gate_sanger_unique` duplicates `check:in_silico_pcr`. `ssr`:
`gate_ssr_flanks` only asserts "outside the repeat", not flank _uniqueness_; no
polymorphism/allele-size gate; `gate_ssr_clean` is vacuous (guaranteed upstream);
`gate_ssr_offtarget` is soft yet README calls a full-3'-end off-target fatal, and it is
also in `KASP_UNGATED_CHECKS` (C10). `indel_size`: `gate_indel_size_gap` is soft but the
designer hard-rejects `< min_gap` (unreachable rung); both hard gates (`structure`,
`unique`) are never enforced; `gate_indel_size_unique` conflicts with `KASP_UNGATED_CHECKS`
(C10).

### Status — §8 closed (2026-11-19)

The briefs above are the static audit's *pre-rebuild* view. They are kept as the
record of what was wrong, not as a task list: every row is now either fixed in
`R/` or pinned by a test that fails if it regresses. The end-to-end guard is
`tests/testthat/test-decision.R:238`, which runs `expect_dag_invariants()` over
`marker_schemes()` (all 20) for C1–C11 and C14, and `:244`, which checks C13
(`scheme_designer()` resolves) over `per_locus_schemes()`. C12 is pinned by the
constructor signature row in the status delta above. Each scheme also carries
its own explicit call, so a failure names the scheme:

| Scheme | Own assertion |
| --- | --- |
| `snp` | `test-selection.R:478` |
| `indel` | `test-indel.R:74` |
| `multi` | `test-mnp.R:70` |
| `primer` | `test-decision.R:100` |
| `duplex` | `test-duplex.R:182` |
| `asq` | `test-asq.R:10` |
| `lax` | `test-lax.R:248` |
| `caps` | `test-caps.R:171` |
| `arms` | `test-arms.R:130` |
| `hrm` | `test-hrm.R:91` |
| `castpcr` | `test-castpcr.R:152` |
| `crispr` | `test-crispr.R:213` |
| `rpa` | `test-rpa.R:147` |
| `lamp` | `test-lamp.R:236` |
| `dpcr` | `test-dpcr.R:145` |
| `taqman` | `test-taqman.R:147` |
| `snapshot` | `test-snapshot.R:129` |
| `sanger` | `test-sanger.R:226` |
| `ssr` | `test-ssr.R:130` |
| `indel_size` | `test-indel-size.R:113` |

C5, C8, C13 and C14 — the four the delta table called “the real remaining work” —
landed as test-only changes in `4be0b7c`, with no `R/` edit, because the tree had
already converged on the contract; the commit only made the silent half
enforced. Anything still listed as a “blocker” in a brief above is therefore a
property of the host or of the underlying physics, not an open code defect.

---

## 9. Sequencing & integration

```mermaid
graph TD
    A1[Agent 1 foundation A: taxonomy + checks + harness] --> A2[Agent 2 foundation B: DAG generator + floors + extraction]
    A2 --> A3[Agent 3 snp/indel/multi/primer]
    A2 --> A4[Agent 4 duplex/asq]
    A2 --> A5[Agent 5 lax]
    A2 --> A6[Agent 6 caps/arms]
    A2 --> A7[Agent 7 hrm/castpcr]
    A2 --> A8[Agent 8 crispr/rpa/lamp]
    A2 --> A9[Agent 9 dpcr/taqman/snapshot]
    A2 --> A10[Agent 10 sanger/ssr/indel_size]
    A3 --> I[Integrator: full gate + figures + NEWS/ROADMAP]
    A4 --> I
    A5 --> I
    A6 --> I
    A7 --> I
    A8 --> I
    A9 --> I
    A10 --> I
```

- Agents 3–10 run in parallel only after Agent 2 reports done.
- Each returns its §6b report + a green `test_local(filter=...)` for its chunk.
- The **integrator** (the user, or one more serial pass) runs the full gate
  `Rscript dev/check.R`, regenerates committed figures, and reconciles `NEWS.md`/`ROADMAP.md`.
- If two agents need the same shared file, serialize those two (the orchestrator decides).

---

## 10. How to run anything (the env is verified working)

The deps are installed in a **renv cache library** and the external tools live in the
project `.tools/`. They are **not** in the default R library, so every command must be
prefixed with the env below or `pkgload::load_all()` fails at import and nothing runs.
This is already recorded in `AGENTS.md`.

```sh
export R_LIBS=/home/stathmin/.cache/R/renv/library/kasp-example-cee55197/linux-debian-forky/R-4.6/x86_64-pc-linux-gnu
export PATH=/home/stathmin/Documents/devel/tidymarkers/.tools/bin:$PATH
export UNAFOLDDAT=/home/stathmin/Documents/devel/tidymarkers/.tools/share/oligoarrayaux

Rscript -e 'pkgload::load_all(".")'                       # load from source
Rscript -e 'testthat::test_local(filter = "<scheme>")'   # one scheme's tests
Rscript dev/check.R                                        # the full merge gate
```

`Rscript dev/check.R` is green on a clean checkout today, so agents can and must run
tests. `Rscript dev/dump-dags.R` (from Agent 1) works offline — use it as the source of
truth for a scheme's current DAG.
