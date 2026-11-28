# 0019 — The marker crate: GWAS loci and a dossier become one reproducible crate

- Status: proposed
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0009](0009-panel-orchestrator-contract.md) (the panel
  orchestrator), [ADR 0010](0010-index-strategy-at-scale.md) (indexed
  specificity), [ADR 0015](0015-a-record-resolves-only-its-own-directory.md)
  (a record resolves only its own directory)

## Context

`marker_design_panel()` starts at a target table. A GWAS programme starts one
step earlier: a locus table, a locus dossier with candidate genes and evidence,
and a population reference that decides which variant at the locus is worth
assaying. The panel can design a KASP assay and write an artefact, but it
cannot say which GWAS run produced the locus, which dossier row supported it,
which linked candidate was rejected, or whether the chosen primers are specific
in the reference build. Those four facts are the difference between a designed
assay and a crate another laboratory can reproduce.

The package already owns every brick the crate needs. The new contract is one
orchestrator over them, not a second panel pipeline.

## Decision

### 1. `marker_locus_crate()` is the only crate entry point

```r
marker_locus_crate(gwas, dossier, genome, variants = NULL, haplotypes = NULL,
                   chemistries = "kasp", index = NULL, maf_min = NULL,
                   ld_min = NULL, dir = tempfile("marker-crate"),
                   created = NULL, ...)
```

It returns the crate directory, a character vector of length one. The directory
is an RO-Crate; `marker-crate.json` is the crate-level record. No new S3 class
is introduced: the crate is a directory of files, and the existing
`marker_artefact` and `marker_export_rocrate()` carry the typed part.

`gwas` and `dossier` are a path or a data frame. `genome` is a `DNAStringSet`
or a `marker_genome`. `variants` and `haplotypes` use the shapes the panel
already accepts; `variants` gains the crate key `locus`. `index` is a
`genome_index` or a named list of indexes (a pangenome); `NULL` builds
`genome_index(genome, strategy = "memory")`, which is the honest default only
for a fixture-sized genome. A production caller passes the indexed genome.

`chemistries` is a scheme id vector or a `KASP_TECH_PRESETS` name; the default
`"kasp"` expands to the preset's schemes. `maf_min` and `ld_min` are optional
inclusive lower bounds. No threshold is invented: `NULL` records
`filter = "not_applied"`.

### 2. The input contract is checked before any design runs

**GWAS locus table** — one lead row per locus:

| column | required | meaning |
| --- | --- | --- |
| `locus` | yes | locus id, unique in the table |
| `chrom`, `pos`, `ref`, `alt` | yes | lead variant |
| `population` | yes | the population the LD and MAF refer to |
| `genome_build` | yes | build of `chrom`/`pos` |
| `gwas_run` | yes | run id or accession of the summary-statistics run |
| `p`, `beta`, `se`, `maf` | no | signal statistics, carried through |
| `lead` | no | logical; when present it must be `TRUE` on the one row per locus that is the lead |

A repeated locus, an empty locus id, a non-positive `pos`, a REF equal to ALT,
a missing build or run aborts and names the row. The table is hashed and copied
to `gwas-loci.csv`; the hash and run id enter both the artefact inputs and
`marker-crate.json`.

**Locus dossier** — one or more rows per GWAS locus:

| column | required | meaning |
| --- | --- | --- |
| `locus` | yes | joins the GWAS table |
| `gene` | yes | candidate gene id |
| `source` | yes | database or publication that supplied the row |
| `source_version` | yes | release or version read |
| `retrieved_at` | yes | ISO date the row was read |
| `confidence` | yes | one finite number in `[0, 1]` |
| `pathway`, `citation`, `doi` | no | enrichment and evidence, carried through |

Every GWAS locus needs at least one dossier row; a dossier locus absent from
the GWAS table aborts. Confidence is copied, never recomputed and never
thresholded. The whole table is copied to `dossier.csv`; the artefact records
its hash and row count, and `marker-crate.json` records the distinct sources,
versions, retrieval dates and confidence range. The dossier is evidence, not a
filter: no row is dropped.

**Candidate variants** — optional, keyed by locus:

`locus`, `chrom`, `pos`, `ref`, `alt` are required; `maf`, `ld_r2` and
`population` are optional. The GWAS lead is always a candidate even when the
table omits it. A candidate whose `population` is present and differs from the
locus population is dropped and counted in `selection.csv` with
`status = "population"`. A candidate that duplicates the lead is collapsed to
one row.

**Haplotypes** — the ADR 0006 shape (`haplotype`, `freq`, `chrom`, `pos`,
`ref`, `alt`). They are the dosage source for LD and the population model for
MAF. Without them, LD and MAF can only come from the candidate columns; an
unmeasured value stays `NA` and is never read as zero.

### 3. One candidate per locus, selected deterministically

For each locus, candidates are ranked by:

1. `ld_r2` descending against the lead (the lead itself is `1`);
2. `maf` descending;
3. the lead before a non-lead candidate;
4. `pos` ascending, then `alt` ascending.

`ld_r2` is the candidate column when present, otherwise the squared Pearson
correlation of per-haplotype carrier indicators at the lead and candidate
positions — the same dosage rule as `panel_ld()`, factored into one helper. A
haplotype absent at a position is a reference call, as in the panel.
`maf = min(p, 1 - p)`, where `p` is the candidate `maf`/`alt_freq` column,
otherwise the haplotype `alt_freq` from `marker_informativeness()`, otherwise
the GWAS `maf` for the lead. The source of each value is written to
`selection.csv` (`ld_source`, `maf_source`).

`ld_min` drops a candidate whose measured r² is below the bound. `maf_min`
does the same for MAF. A bound applied to an unmeasured value aborts and names
the locus and candidate: "cannot filter" is not "passes". The best remaining
candidate is selected; the selection table keeps every candidate with its rank,
values, sources and status (`selected`, `ld`, `maf`, `population`, `duplicate`).
If no candidate survives, the crate aborts naming the locus and the filter that
emptied it.

### 4. Design, specificity and multiplex reuse the existing stages

The selected rows go to `marker_design_panel()` with the same `genome`,
`variants`, `haplotypes`, `chemistries`, `created` and `...`. The panel owns the
haplotype contract, chemistry choice, per-target cause chain, one-pot grouping
and artefact. The crate does not design an assay or group a mix itself.

After a successful design, each designed assay's primers come from
`marker_primers()`, are prefixed with the locus id, and go to
`specificity_of(primers, index)`. The scan is advisory: `specificity.csv` keeps
one row per primer with `locus`, `target`, `primer`, `seed`, `hits`, `engine`
and `build`. The crate records the index build and strategy; it does not turn a
hit count into a pass/fail threshold, because the package has no validated
specificity threshold. The one-pot groups are the panel's `one_pot` record.

### 5. Provenance reaches the GWAS run and the dossier

The crate attaches `plan$provenance` records:

- `role = "gwas"`: source `gwas:<run>`, path, SHA-256, bytes, locus count and
  genome build;
- `role = "dossier"`: source, path, SHA-256, bytes, row count and source count;
- `role = "candidates"` when a candidate table was read.

`artefact_inputs()` appends `plan$provenance`, so the hashes enter
`artefact.json`, the artefact id and the RO-Crate PROV entities. The artefact is
then rewritten in the same directory with the same `created`, so the manifest
and the crate agree. `marker-crate.json` additionally names the run,
population, build, dossier sources/versions/retrieval dates/confidence range,
selection rows, panel id, one-pot groups, specificity summary and limitations.

### 6. The crate directory is self-contained and deterministic

```
<dir>/
  ro-crate-metadata.json
  artefact.json
  decisions.csv checks.csv why.csv signals.csv targets.csv replay.csv
  gwas-loci.csv dossier.csv selection.csv specificity.csv marker-crate.json
  artefact/            # the panel's artefact directory, including solved.rds
```

`solved.rds` stays in `artefact/` and is never copied: it is byte-unstable and
executing it is not part of reading a crate. `marker_export_rocrate()` gains a
small `extra =` argument: the crate sidecars are copied beside the manifest
files and added to the root `hasPart` and the run `result`, while
`marker_import_rocrate()` still verifies only declared file hashes and the root
identifier. With a fixed `created`, every text file is byte-identical across
runs; the crate id hashes the contract, timestamp, input hashes, selection,
panel id, one-pot groups and specificity table.

## Alternatives rejected

- **A second panel pipeline inside the crate.** It would fork the haplotype
  contract, chemistry choice and one-pot rule that ADR 0009 just unified.
- **Dropping the dossier rows after computing confidence.** The crate exists to
  show the evidence behind the locus, including rows that did not change the
  choice.
- **A specificity pass/fail threshold.** The package has no validated threshold;
  a number invented here would look like one.
- **Selecting several markers per locus by default.** One lead-linked marker is
  the reproducible minimum. A future `max_markers` is a new contract, not a
  hidden default.
- **Copying `solved.rds` into the crate.** It is excluded by name for the same
  reason the artefact excludes it.

## Consequences

- A crate names its GWAS run, dossier evidence, candidate selection, chosen
  chemistry, primers, specificity scan and multiplex grouping in one directory.
- The panel's typed artefact is still the design record; the crate adds the
  programme context around it.
- LD and MAF are only as good as the supplied haplotypes or candidate columns;
  the sources are written, and a missing value is `NA`, never a silent zero.
- The default memory index is for fixtures; a production genome passes an
  index, and the crate records which one.

## Verification

- `tests/testthat/test-crate.R`: a committed synthetic soybean-named locus
  fixture (the reference sequence is not the real Wm82 assembly and says so)
  produces a crate twice with the same `created`; every text file is compared
  byte for byte, the selected candidate, LD/MAF sources, GWAS run, dossier
  provenance and RO-Crate root id are asserted.
- Error tests cover a missing column, a repeated locus, a dossier locus absent
  from the GWAS table, a bound on an unmeasured value and an empty selection.
