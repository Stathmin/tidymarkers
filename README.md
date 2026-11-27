<img src="man/figures/logo.png" align="right" width="150" alt="tidymarkers logo" />

# tidymarkers

<!-- badges: start -->

[![R-CMD-check](https://github.com/Stathmin/tidymarkers/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Stathmin/tidymarkers/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
<!-- badges: end -->

Design and quality control of genetic markers in pure R.

## What problem this solves

Choosing a marker is a search over sequence, chemistry, and genomic
specificity, and most tools answer only one part of it: a primer designer
returns oligos, a specificity script returns hit counts, and the decision
between them stays undocumented. `tidymarkers` takes a reference (FASTA, an
indexed `marker_genome()`, a `DNAStringSet` or a Bioconductor `BSgenome`) and
variation (a VCF or the tidy variant table), designs allele-specific and common
primers on either strand, scores them against hard gates, and returns the
candidate that survived together with the reason every other candidate did not.
It covers KASP, CAPS/dCAPS, tetra-primer ARMS, HRM, Sanger, SSR, InDel-size,
TaqMan/dPCR, SNaPshot, LAMP/RPA and CRISPR-Dx under one grammar, with
thermodynamics from `TmCalculator` and OligoArrayAux, Taq-aware discrimination
from `DECIPHER`, and genome-wide specificity from `Biostrings` (optionally
through `rBLAST` or an `Rsubread` FM-index). tidymarkers is three things at
once: a marker design tool for a single researcher (Level 1), a unified grammar
for marker design across chemistries (Level 2), and infrastructure for
community-wide method refinement through aggregated design traces (Level 3).
Each level stands on its own; together they define a new way to work.

## Key validated results

- **3′-seed specificity: 8/8 equal to exhaustive enumeration.** `Rscript dev/specificity-oracle.R` counts every ≤1-substitution seed hit on TAIR10 by brute force and compares it with `specificity_seed()`; `tests/testthat/test-specificity.R` holds the equality on fixtures.
- **Exact in-silico PCR: 3/3 equal.** The same script enumerates exact primer pairs with `vmatchPattern()` and compares them with `in_silico_pcr(max_mismatch = 0)`; `results/specificity-oracle.csv` records both columns.
- **Measured Tm ranking: r = 0.958 on 4096 19-mers.** `Rscript dev/test.R tm-mm2020` reproduces the Oliveira et al. table and its correlation; the measurement is committed as `inst/extdata/tm_mm2020_s1_measured.csv`.
- **134 registered literature sources.** `Rscript dev/check-docs.R` re-derives the count from `inst/extdata/scientific_sources.csv`, where every row carries its DOI and access level.
- **223 ledger `test` links across 377 ledger rows.** `Rscript dev/test.R science-ledger` requires each link to name an existing test that actually reads the named extract; `inst/extdata/science_ledger.csv` is the join.
- **213 criterion rows, 0 fitted here.** `Rscript dev/check-criteria.R` walks the registry, and `docs/decision-audit.md` prints the rung counts: 128 `primary_literature`, 13 `vendor_guideline`, 72 `rule_of_thumb`.
- **54 tier rows, 36 with a printed target.** The other 18 rank a quantity the package computes or fixes as a scheme default and stay `rule_of_thumb`; `docs/decision-audit.md` names the split.

The thermodynamic engines are compared, not blended. `thermo_crosscheck()`
reports SantaLucia 2004, Allawi 1998, Sugimoto 1996 and Breslauer 1986
separately, Primer3's `oligotm` and `rprimer` join as external oracles when
installed, and no disagreement is averaged into a correction factor. Against
the 4096 measured 19-mer Tm values the ranking is intact while the level sits
low; three recorded steps account for it (−16.9 °C from the default salt term,
−9.6 °C after the NN-consistent `SantaLucia1998-2` correction, ≈ −5.4 °C after
matching the set's own 1 µM strand), and the offset is recorded, not corrected
away; the full comparison is in
[`VALIDATION.md`](VALIDATION.md#cross-tool-comparison).

## Install

`tidymarkers` is not on CRAN; install it from GitHub. `pak` also resolves the
system libraries (BLAST+, cairo):

```r
install.packages("pak")
pak::pak("Stathmin/tidymarkers", sysreqs = TRUE)
```

`DECIPHER` needs [OligoArrayAux](https://www.unafold.org/Dinamelt/software/oligoarrayaux.php);
one call builds it into your user cache, so nothing has to be on `PATH`:

```r
tidymarkers::marker_install_tools()   # once; hybrid-min / hybrid-ss-min
```

**R >= 4.5** is required. The hard `Imports` are Bioconductor plus the
tidyverse; the heavy optional pieces (a `TxDb`/`BSgenome`, `Rsubread`,
`ensembldb`, `rBLAST`, `Rsamtools`, `systemfonts`) live in `Suggests` and are
never installed implicitly: each is reached only through a guarded call, so a
missing one degrades a single feature instead of breaking the load.
`marker_cache_dir()` reports where downloaded genomes, variation, indexes and
the external tools are kept (`tools::R_user_dir()`), so nothing is written
inside the package or the project.

## Quick start

One shipped FASTA, one assay, one verdict. The snippet is the executable
`quickstart` chunk of `vignette("kasp-intro")`, line for line, so it runs on
every build where the Tm engine is present rather than being prose that merely
looks right (`dev/check-docs.R` compares the two):

```r
fa <- system.file("extdata", "sat2_transcripts.fa", package = "tidymarkers")
genome <- read_genome_fasta(fa) # Rsamtools index
snp <- tibble::tibble(chrom = names(genome)[[1L]], pos = 600L, ref = "G", alt = "A")
assay <- assemble_kasp_assay(genome, snp, flank = 150L, cores = 1L)
assay # prints a summary
marker_design_checks(assay) # row by row, what and why
passes_all(assay) # TRUE when all criteria pass
```

Any other reference is passed in directly: a FASTA path through
`marker_genome()` / `read_genome_fasta()`, a `DNAStringSet`, or a Bioconductor
`BSgenome` (see `vignette("custom-genome")`); `ath_genome()` is the cached TAIR10
demo. Two project-level knobs are arguments rather than constants: `strands`
(search the AS primer on `"+"`, `"-"` or both) and `tm_range` (the permitted
melting window, °C).

## Pipeline in six steps

1. **Input.** `read_genome_fasta()` / `marker_genome()` / `read_variants_vcf()`
   read a reference and variation; `ensembl_gene()` / `ensembl_variants()` pull
   a gene model and SNPs from Ensembl Plants (1001 Genomes), cached and retried.
2. **Candidates.** `pick_target_snp()` chooses a designable substitution
   (clean 3′ AS window first, then missense > splice > synonymous > UTR);
   `design_as_primers()` and `design_common_primers()` scan strand × length ×
   mismatch position × base.
3. **Specificity.** `specificity_seed()` counts 3′-seed hits,
   `in_silico_pcr()` checks product uniqueness and size, and `blast_primers()`
   adds a BLAST baseline when an index is available.
4. **Assembly.** `assemble_kasp_assay()` builds one assay for a SNP, MNP or
   InDel; `design_kasp_markers()` runs a batch and returns a
   `pass`/`warn`/`fail` verdict per locus.
5. **Decision.** Candidate pairs face hard gates first (known variant under the
   primer, an allele that cannot be discriminated, `Tm_COM < Tm_AS`, an amplicon
   outside the window, impossible secondary structure). Survivors are ranked
   lexicographically by tier; within a tier the bottleneck `min` decides. If
   nothing clears the gates they are lifted one at a time, least grave first,
   and the least-bad candidate is returned with a warning naming every violated
   gate. `marker_design_checks()` reports each criterion, `marker_explain()`
   orders the reasons, and `marker_decision_text()` prints the graph.
6. **Artifact.** `marker_export()` writes `markers.csv`, `attempts.csv`,
   `primers.csv` and `solved.rds`; `marker_artefact()` adds `artefact.json`
   (hashed inputs and provenance), `decisions.csv`, `checks.csv` and `why.csv`;
   `marker_replay()` re-runs a saved artefact, and `marker_record_outcome()`
   appends a lab result bound to its id.

<img src="man/figures/decision-tree.svg" width="470" alt="KASP selection decision tree" />
<img src="man/figures/decision-multiplex.svg" width="470" alt="duplex KASP decision tree" />

The decision graph is data, not prose: `marker_decision_schema()` reads the same
constants as the engine, so it cannot drift from the code. Each tier criterion
carries a citable `target` (an amplicon window, a ΔG floor, a vendor's Tm range)
and a package `limit` (the bottleneck policy, tier membership and order); the
policy is declared once in `CRITERION_TIER_POLICY_EVIDENCE` and stays
`rule_of_thumb` however well sourced the target is.

## For the manuscript

| Manuscript element | Command | Artifact |
| --- | --- | --- |
| 3′-seed and exact-PCR specificity | `Rscript dev/specificity-oracle.R` | `results/specificity-oracle.csv`, `tests/testthat/test-specificity.R` |
| Primer3 2.6.1 head-to-head | `Rscript dev/primer3-crosscheck.R` | `results/primer3-thermo.csv`, `results/primer3-design.csv` |
| Allele-specific mismatch panel | `Rscript dev/as-tm-crosscheck.R` | `results/as-tm-panel.csv`, `results/as-tm-crosscheck.csv` |
| Measured Tm benchmark | `Rscript dev/test.R tm-mm2020` | `inst/extdata/tm_mm2020_s1_measured.csv`, `tests/testthat/test-tm-mm2020.R` |
| Threshold provenance and tiers | `Rscript dev/check-criteria.R` | `docs/decision-audit.md`, `inst/extdata/criterion_provenance.csv` |
| Evidence registry and ledger | `Rscript dev/check-docs.R` | `VALIDATION.md`, `inst/extdata/scientific_sources.csv`, `inst/extdata/science_ledger.csv` |
| Decision figures | `Rscript dev/render-all.R` | `man/figures/decision-tree.svg`, `man/figures/decision-multiplex.svg`, `vignettes/figures/decision-*.svg` |
| Artifact graph | `Rscript dev/artifacts.R --write` | `docs/generated/artifacts.md` |
| Validation benchmark | `Rscript dev/validation-benchmark.R --check` | `.tmp/validation-benchmark.csv`, `VALIDATION.md` |
| Performance baseline | `Rscript dev/bench.R` | `dev/bench-baseline.csv` |

**Environment.** R 4.6.1 and Bioconductor 3.23, pinned in `renv.lock` (the
package itself requires R >= 4.5); no `.Rprofile` or `renv/` directory is
committed, so restore explicitly. The specificity panel and the Primer3 oligo
set are seeded (`set.seed(20261023L)`, `set.seed(20040301L)`); the mm2020 table
is a fixed 4096-row measurement and the design pipeline carries no stochastic
step. Generated `results/` and `.tmp/` files are git-ignored, and every command
above re-creates them.

**Data availability.** No genome ships with the package. TAIR10 and the 1001
Genomes variation are downloaded from Ensembl Plants into `marker_cache_dir()`;
the literature extracts that back the rules are committed under
`inst/extdata/` with DOIs and access levels in `scientific_sources.csv`. The
rule-to-source mapping and every known gap are in
[`VALIDATION.md`](VALIDATION.md).

## Container and MCP

A pinned image (`rocker/r-ver:4.6.1`) builds the package with its `Imports` and
starts a stdio MCP server; see [`docker/README.md`](docker/README.md). The
container opens no port (JSON-RPC on stdin/stdout), mounts data **read-only**
(`-v /path/to/data:/data:ro`), and needs no network at run time. The artefact's
decision log carries `outcome`, `verdict`, `relaxed`, `violated` and `risk` as
typed columns, so an agent reads a field instead of parsing a sentence.

The registry (`inst/mcp/tools.R`, 38 tools) follows the vignette workflow: read
a genome/VCF, screen variants, run thermodynamics and QC, check specificity,
design one marker (`design_scheme`) or a batch (`design_markers`), and read the
verdict. The design tools take `genome_path`, `variants_path`, `blast_db` and
`oligoarrayaux_dir` as mounts under `/data`. Plotting calls, the `marker_plan()`
/ `marker_solve()` objects and the writers stay out: they return structures
JSON cannot carry, or they write outside the mount.

## Reproducibility

`renv.lock` pins the exact versions used in development. It is deliberately not
activated by an `.Rprofile` or a `renv/` directory: activation makes
`R CMD build` / `R CMD check` on a fresh clone resolve against an empty
library. Use `renv::init()` then `renv::restore()` for the pinned snapshot;
development installs use `Rscript dev/install-deps.R`.

External data never ships with the package: TAIR10, the 1001 Genomes variation
and the thermodynamic binaries are fetched or built once inside
`marker_cache_dir()`, from an FTP path on `ftp.ebi.ac.uk` rather than this
repository. Genome indexes are **not** cached: `marker_bowtie_index()` and
`marker_fm_index()` default to a basename under `tempdir()`, so a large index
dies with the session. That default is cheap for TAIR10 and wrong for wheat,
where the index is about 1.4× the FASTA (≈21.5 GB for 16 Gb); build once into
`marker_cache_dir()` and delete both sides when done. A BLAST database is
written beside its FASTA (or into `tempdir()`); delete the `.nhr`/`.nin`/`.nsq`
files to force a rebuild.

## Limitations

- `DECIPHER` amplification efficiency is a model (hybridisation plus Taq
  extension), a proxy for allelic discrimination rather than a wet-lab result;
  it depends on buffer and temperature.
- The ALT allele is scored in the reference context (there is one genome), which
  is not the real ALT template; frequencies and MAF come from 1001 Genomes as
  published, so production use needs your own QC (depth, MAF, Mendelian errors).
- Hairpins and dimers are real nearest-neighbour thermodynamics at 37 °C
  (OligoArrayAux); matrix secondary structure is not modelled.
- The 3′ cleanliness of an AS primer is checked against the 1001 Genomes
  catalogue (`n_3p`), which is not exhaustive: rare variants in your own samples
  are invisible to it.
- No `tidymarkers` design has been validated at the bench by the authors. The
  repository supplies reproducibility evidence, not experimental evidence.

## How to cite

`tidymarkers` is research software by Daniil Ulyanov at the All-Russian
Research Institute of Agricultural Biotechnology ([VNIISB](http://www.vniisb.ru/ru/news/)),
Moscow. If it contributed to a publication, please cite the package and the
version you ran:

```r
citation("tidymarkers")
```

which at version 0.4.0 prints:

```bibtex
@Manual{,
  title = {tidymarkers: Tidy Design and Quality Control of Genetic Markers},
  author = {Daniil Ulyanov},
  organization = {All-Russian Research Institute of Agricultural Biotechnology (VNIISB)},
  year = {2026},
  note = {R package version 0.4.0},
  url = {https://github.com/Stathmin/tidymarkers},
}
```

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff) (GitHub and
Zenodo), so a release can be cited by DOI once archived; the author's ORCID is
[0000-0002-5880-5931](https://orcid.org/0000-0002-5880-5931). Beyond the package,
cite the primary source for any rule that shaped your design, not the package
alone; the rule-to-source mapping is in [`VALIDATION.md`](VALIDATION.md).

## License

MIT, copyright 2026 Daniil Ulyanov. The full text is in [`LICENSE.md`](LICENSE.md);
[`LICENSE`](LICENSE) holds the `YEAR`/`COPYRIGHT HOLDER` stub that R's
`DESCRIPTION` (`MIT + file LICENSE`) requires.
