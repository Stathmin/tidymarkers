# tidymarkers — documentation index

Canonical entry point for every document in the repository. `dev/check-docs.R`
verifies that each document is linked from this page and that every relative
link resolves; `dev/check-artifacts.R` verifies the artifact graph.

## Project

- [README](../README.md) — what the package is, install, quick start, decision trees.
- [AGENTS](../AGENTS.md) — agent instructions and the environment to run in.
- [CONTRIBUTING](../CONTRIBUTING.md) — conventions, scope and the merge gate.
- [CREDITS](../CREDITS.md) — contributor record and how work is ranked here.
- [ROADMAP](../ROADMAP.md) — the open queue (ZA8 is this docs/artifacts layer).
- [NEWS](../NEWS.md) — changelog.
- [VALIDATION](../VALIDATION.md) — what is checked against external evidence,
  what is not, and where the gaps are.

## Archive

- [archived inventory](archive/inventory-ru.md) — the pre-0.4.1 working
  inventory (Russian); kept as history, not as current documentation.

## Vignettes

Marker chemistries:

- [KASP and allele-specific designs](../vignettes/markers-kasp.Rmd) — `snp`,
  `indel`, `multi`, `duplex`, `lax`, `asq`, `primer`.
- [Gel and melt designs](../vignettes/markers-gel.Rmd) — `caps`, `arms`, `hrm`,
  `ssr`, `indel_size`.
- [Probe designs](../vignettes/markers-probe.Rmd) — `taqman`, `dpcr`,
  `snapshot`, `castpcr`.
- [Isothermal designs](../vignettes/markers-isothermal.Rmd) — `lamp`, `rpa`.
- [CRISPR and Sanger sequencing](../vignettes/markers-crispr-sequencing.Rmd) —
  `crispr`, `sanger`.

Other vignettes:

- [KASP introduction](../vignettes/kasp-intro.Rmd)
- [Custom genome](../vignettes/custom-genome.Rmd)
- [Reading results](../vignettes/reading-results.Rmd)

## Developer docs

- [dev/README](../dev/README.md) — the internal checks, renderers and benchmarks.
- [dev/dag-rebuild-plan](../dev/dag-rebuild-plan.md) — the decision-DAG rebuild:
  plan, invariants (C1–C14), agent scopes.

## Docs and artifacts layer

- [docs/INDEX](INDEX.md) — this page.
- [docs/artifacts.yml](artifacts.yml) — the artifact manifest (source of truth).
  Committed artifacts only: the logo, the two README decision trees, the 20
  per-chemistry decision SVGs in `vignettes/figures` and the 20 vignette gallery
  SVGs. The `.txt` twins of the decision SVGs are committed too; only the
  `vignettes/figures/*.png` rasters are generated locally and ignored.
- [docs/generated/artifacts.md](generated/artifacts.md) — the generated artifact
  graph; one heading per artifact with links to the artifact, its generator,
  its code and the docs that consume it.

Refresh and check the manifest after changing a figure, its generator or the
engine it reads:

```sh
Rscript dev/artifacts.R --write    # hashes + docs/generated/artifacts.md
Rscript dev/check-artifacts.R      # stale/orphan gate (exit 1)
Rscript dev/check-docs.R           # orphans, links, include_graphics, backlinks
```
