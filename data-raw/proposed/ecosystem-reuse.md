# Z4.2 — ecosystem reuse, not reimplementation: a concrete audit

**Verdict: keep own today; wrap nothing; adopt nothing new.** Every ecosystem
piece Z4.2 names was checked against the code that would consume it, and each
one either does not cover the need (different input object, different question)
or has no consumer in the code yet — which is Z4.2's own blocker ("a second
dependency with no call site is dead code"). No `DESCRIPTION` change follows
from this note. Read 2026-11-03 against the working tree; `crisprBase` rows are
pinned by `tests/testthat/test-crispr-crisprbase-registry.R`.

Three labels are used per row: **covers** (the ecosystem function answers the
same question on the same input), **partial** (same question, different input
object or a narrower answer), **no** (different object or different question).

---

## 1. What tidymarkers reimplements, and what sits next to it

Current state, from the code and `DESCRIPTION`, not from the roadmap's plan:

| need                          | owned code                                                                                          | ecosystem candidate                                             | already a dependency?          |
| ----------------------------- | --------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- | ------------------------------ |
| PAM / PFS scan                | `R/marker_crispr.R` (`^TTT[ACG]$` Cas12a; Cas13a 3′ H, Cas13b 5′ D)                                 | `crisprBase` nuclease objects                                   | no (`crisprVerse` in Suggests) |
| seed / position gate          | `crispr_seed_gate()`, position only                                                                 | `crisprScore` on-target models                                  | no                             |
| CRISPR specificity            | `crispr_spacer_hits()` — exact match, both strands, minus 1                                         | `crisprBowtie` / `crisprBwa`                                    | no                             |
| primer design grid            | `design_as_primers()`, `design_as_primers_indel()` + KASP filters                                   | `DECIPHER::DesignPrimers` / `DesignProbes` / `DesignSignatures` | `DECIPHER` in Imports          |
| degenerate / alignment panels | `R/marker_alignment.R` (`DECIPHER::AlignSeqs`), `R/marker_rprimer.R` (design + `checkMatch` report) | `rprimer`                                                       | `rprimer` in Suggests          |
| Tm                            | `primer_tm()` via `TmCalculator::tm_calculate()`                                                    | `TmCalculator` (27 NN, salt/DMSO)                               | Imports, used                  |
| allelic discrimination        | `DECIPHER::CalculateEfficiencyPCR()`                                                                | same                                                            | Imports, used                  |
| oligo ΔG (hairpin/dimer)      | `oligo_dg()` via OligoArrayAux binary                                                               | `DECIPHER` (same engine)                                        | external binary                |
| wheat homeolog rules          | `inst/extdata/polymarker_design_rules.csv` (method + line)                                          | PolyMarker                                                      | reference only, by design      |
| alternative Tm oracle         | optional `oligotm` binary (Primer3)                                                                 | Primer3 / MFEprimer / NUPACK / ViennaRNA                        | external binary, optional      |

## 2. `crisprBase` — nucleases, PAM, spacer

**Covers: partial.** `crisprBase` has eight nuclease objects, and the registry
table `data-raw/proposed/crisprbase_registry.csv` (crisprBase 1.13.2) pins what
they carry: `motifs`, `pam_side`, `spacer_length`, but **no seed slot** (`seed_model ==
"none (class has no seed slot)"` on every row) and no PFS field — a PFS appears
only for the RNA-targeting classes as a catch-all `N`.

Against `R/marker_crispr.R`: PAM motif and spacer length would be covered; the
PFS the package needs for Cas13a (3′ H, G rejected) and Cas13b (5′ D) is **not**
a slot, and Cas13 is Cas13d (`CasRx`) only — there is no Cas13a and no Cas13b
object. The seed gate is positional and would still be ours.

**Recommendation: wrap as a data cross-check only** (a test, not an adapter).
Reason: the PAM/spacer half can be verified against the registry; the half Z4.2
would need it for (PFS, seed, Cas13a/b) it does not carry, and
`tests/testthat/test-crispr-crisprbase-registry.R` already asserts the table this
would check — a wrapper would restate that test with a dependency attached.

## 3. `crisprScore` and `crisprDesign` — seed, activity, annotation

**Covers: no.** `crisprScore` predicts **on-target cleavage activity** from
sequence and chromatin context (and, for the annotation functions, puts guides on
variants). The package's question is different: does one reference spacer
carry a _second_ exact site, and where does the variant sit in the spacer. Its
trans-cleavage readout is not the on-target activity these models are fitted to,
so an activity score is not a predictor of the collateral signal — Z4.2 says the
same ("not as a trans-cleavage activity predictor").

`crisprDesign` adds an annotation layer with a human-genome bias; the package
works on plant genomes, and `Y3.5` (wheat, ~15 Gb) would want a genome-agnostic
off-target search, not a variant annotator.

**Recommendation: keep own.** Reason: adopting them would add two Suggests with
no call site, and the one number they supply (on-target activity) is not the one
the design gates on.

## 4. `crisprBowtie` / `crisprBwa` — off-target search

**Covers: partial.** Both do the thing the package does by hand — enumerate
spacers with a PAM and count near-matches — but as alignment against a human
index. The package's `crispr_spacer_hits()` counts **exact** matches on
`Biostrings`, deliberately (see the decision audit in
`vignettes/markers-crispr-sequencing.Rmd`): a conservative count, no mismatch
search.

**Recommendation: keep own until a real 15 Gb consumer exists (`Y3.5` / `Z5.1`).**
Reason: `Rsubread` (FM index) and `rBLAST` already cover "is this oligo unique in
a big genome" without a new dependency; an aligner backend earns its place only
when its index memory is _measured_ on the 15 Gb genome first.

## 5. `DECIPHER` `DesignPrimers` / `DesignProbes` / `DesignSignatures`

**Covers: no.** Already evaluated in the superseded inventory §N2 and
rejected there: `DesignPrimers` and `DesignSignatures` work on `tiles` (a
`TileSeqs` database of sequence groups) and aim at primer pairs covering groups,
not at a reverse primer sitting opposite a single SNP in one reference window.
`DesignProbes` has the same shape — probes chosen from a set of sequences for
group discrimination.

**Recommendation: keep own.** Reason: the package already takes the physics from
`DECIPHER` (`CalculateEfficiencyPCR`, and the same OligoArrayAux engine for ΔG);
what the grid adds is the candidate enumeration and the KASP filters, which the
tile-based API does not express.

## 6. `rprimer` — degenerate panels and an independent match report

**Covers: partial.** `R/marker_rprimer.R` already wraps it: design on an
alignment and a `checkMatch()` report are available. It cannot re-run
`checkMatch()` on the package's own primers — the S4 objects are locked, so a
foreign oligo cannot be handed back to the checker. That limit is documented as a
scope limit, not a bug.

**Recommendation: keep the existing wrapper (already Suggests).** Reason: the
cross-check it does supply is worth having and is already wired; there is nothing
left to adopt.

## 7. `TmCalculator` and the third Tm oracle

**Covers: yes, and it is used.** `primer_tm()` calls `TmCalculator::tm_calculate()`
with `tm_nn` (SantaLucia 2004 family). `VALIDATION.md` records the three-engine
comparison: `oligotm -tp 1 -sc 1 -mv 50 -dv 0 -n 0` agrees to 4.8e-7 °C, the
Primer3 defaults differ by 0.40 / 0.98 °C, and the package caps at 1.0 °C.

**Recommendation: already adopted; no second use.** Reason: the independent-check
role belongs to `oligotm`, which is an optional binary and stays out of
`Imports` — a second R Tm package would not add an independent oracle, only a
third copy of the same NN family.

## 8. `PolyMarker` — the wheat baseline

**Covers: no, and that is the point.** The rules were transcribed with method and
line into `inst/extdata/polymarker_design_rules.csv` and are exercised by
`test-y32-homeolog.R`. `VALIDATION.md` is explicit that this is a **rule
adoption, not a benchmark**.

**Recommendation: benchmark and rule source, never reimplementation.** Reason: it
is the only published reference for wheat genome-specific KASP; its capability is
the target the package is measured against, not code to vend into it.

## 9. `openPrimeR` and the external tools

`openPrimeR` is not a dependency, and ZA4 forbids adding it without a
maintenance check and a second adapter. It demands MAFFT / ViennaRNA / MELTING
and has no KASP-specific layer; `openPrimeR`/`primeR` were already "сняты".

Primer3, MFEprimer, NUPACK and ViennaRNA stay behind the capability interface —
an optional binary that reports `NA` and a reason when absent, never an `Imports`
entry.

**Recommendation: do not adopt `openPrimeR`; keep the external tools
capability-only.** Reason: neither answers the allele-specific 3′ question, and a
heavy dependency without that answer costs the install time YP is trying to cut.

## 10. What stays ours, and when to revisit

Stays ours after this audit:

- the PFS and seed gates (no ecosystem slot carries them);
- the conservative exact-match CRISPR specificity count (an aligner is a
  different trade, not a superset);
- the KASP design grid and its thresholds;
- the transcribed PolyMarker rules.

Revisit when, and only when, one of these becomes true:

1. a **15 Gb index memory is measured** for `Rsubread` vs a bowtie/BWA backend
   (`Y3.5`, `Y5.1`) — then adopt that aligner as the specificity backend and keep
   the 3′-aware post-filter;
2. `crisprBase` gains a **PFS/seed slot** or a Cas13a/b object — then extend the
   registry cross-check test, still without an adapter;
3. `Z8` publishes negative examples — then metric/resampling packages (`yardstick`,
   `rsample`) are unlocked under `ZA6` step 7, unrelated to this note's rows;
4. a consumer for a CRISPR activity score appears — the score is then recorded as
   `calibration_evidence`, not as a gate.

Until then, adding any of the nine to `DESCRIPTION` is a dependency without a
call site.
