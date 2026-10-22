# Validation and evidence

What in this package is checked against something outside itself, how that
check is enforced, and what has not been validated at all.

## Summary

- **The maintainers have produced no wet-lab data.** Every measured number in
  this repository was measured by the authors of the cited paper, not here.
  There is no bench validation of a `tidymarkers` design.
- **There is no head-to-head design benchmark against Primer3 or PolyMarker.**
  Both are used as *sources*, not as competitors: Primer3's `oligotm`
  nearest-neighbour family is the Tm engine lineage, and PolyMarker's selection
  rules are read from its source and encoded as rules.
- **What does exist** is a literature-anchored evidence layer: 80 registered
  sources, 130 extracted tables in `inst/extdata/`, and a ledger that binds 174
  extracts to the 37 test files that consume them. Each link from a source, to
  the exact table read from it, to the test that asserts a property of it is
  machine-checked.
- **Agreement between models is not evidence about a real duplex.** That
  sentence is the rule the thermodynamic tests are written to respect.

## How the evidence is organised

Three tables carry the provenance, and three tests keep them honest.

| File | Rows | What it holds |
|---|---|---|
| `inst/extdata/scientific_sources.csv` | 80 | One row per paper or software source: citation, DOI/PMIM, **how much was read** (`access`), `what_it_anchors`, and the verbatim `claim` |
| `inst/extdata/science_ledger.csv` | 284 | The join `source_key → extract → test_file`, plus `verdict_moved`: which design or ranking decision the extract actually changed |
| `inst/extdata/scientific_backlog.csv` | 45 | Search results and their status: 35 `extracted`, 6 `not_extractable`, 4 `unfetched` |

`access` is a closed vocabulary and it is load-bearing: 54 sources were read in
full text, 11 through a supplement, 3 as an abstract only, 6 are closed access,
2 are source code, 2 are stubs, 1 was transcribed, 1 was not found. 77 of the 80
rows are marked `verified`, each with the date it was checked. A source read only
as an abstract may not calibrate anything, and the tests enforce that.

Ledger rows by role: 174 `test`, 78 `report` (a claim the paper makes, kept for
comparison), 29 `blocked` (read, extracted, consumed by nothing yet — each says
why), 2 `rule`, 1 `artefact`.

The gates:

| Test | What it refuses |
|---|---|
| `tests/testthat/test-science-ledger.R` | A ledger row naming a test that does not exist or does not actually read the named extract; an extract whose own DOI belongs to a different paper; a source with no ledger entry |
| `tests/testthat/test-scientific-data.R` | A table with a numeric outcome and no source, no `access` level, or a `calibrated = TRUE` flag on something read only as an abstract |
| `tests/testthat/test-falsifiable.R` | A published acceptance criterion that cannot be driven to `fail` by the suite — every criterion in the registry needs a mutation that trips it, or a recorded reason why not |

## Cross-tool comparison

### Thermodynamics: same lineage, four parameter sets, three engines

`primer_tm()` is the Primer3/`oligotm` nearest-neighbour family (TmCalculator,
SantaLucia 2004). That lineage is shared rather than compared, so the package
tests the thing that is actually uncertain: how much the *window* wobbles
between published parameter sets.

```r
thermo_crosscheck(c("GAAGGTGACCAAGTTCATGCT", "TTTTGGCGCGCCGATATATAT"))
```

runs the same oligos through SantaLucia 2004, Allawi 1998, Sugimoto 1996 and
Breslauer 1986 and reports `tm_min`, `tm_max`, `spread`. `tests/testthat/test-tm-oracle.R`
adds Primer3's `oligotm` binary as a third engine when it happens to be
installed, with an explicit tolerance per model pair and structural anchors that
carry no tolerance (Tm rises with GC at fixed length, with length at fixed
pattern, and is symmetric under reverse complement). The three-way disagreement
is reported, never averaged, and never asserted against a single constant.
`dev/compare-thermo.R` prints the same matrix for inspection.

### rprimer

`tests/testthat/test-rprimer.R` designs with the independent `rprimer` package on
the same alignment and checks the match table it returns (perfect, one-mismatch
and off-target columns). This is an external cross-check and runs in the full
suite only.

### PolyMarker

PolyMarker 1.3.3 (commit `d4a4578`) is read in `data-raw/polymarker-source.R`;
every extracted rule names the file, method and line it came from
(`inst/extdata/polymarker_design_rules.csv`). The rule that matters is the
two-axis model: the allele sits on the allele-specific primer's 3′ base, and
subgenome specificity sits on the common primer's 3′ base — they are separate
columns, not one "3′ base" field. `tests/testthat/test-y32-homeolog.R` asserts
the axes stay separate. This is rule adoption, not a benchmark.

### What is missing

No comparison of *designs* against Primer3 or PolyMarker on a common locus
panel: neither tool produces a KASP assay with the same acceptance criteria, so
a like-for-like comparison needs a shared panel and a scoring definition that
does not exist here yet. This is tracked in `ROADMAP.md`.

## Measured versus predicted

| Anchor | Data | What is asserted | What is not |
|---|---|---|---|
| Oliveira et al., *Chem Sci* 2020 (`10.1039/d0sc01700k`) | 4096 measured 19-mer duplex Tm (`tm_mm2020_s1_measured.csv`) | Ordering and discrimination: model and measurement rank the 64 canonical duplexes alike (r = 0.958) and the model separates a perfect duplex from a triple mismatch | Absolute agreement. `primer_tm()` sits ~16.9 °C below the measurements under that buffer (sd 0.64 °C); the offset is recorded, not corrected away |
| *PLoS ONE* 2025 (`10.1371/journal.pone.0321885`) | 5 amplicon Tm, psaA-2 | The table's printed arithmetic; the authors' ordering claim | A test of `primer_tm()` on a primer duplex — this is an amplicon measurement |
| Readout physics (HRM, ddPCR, Sanger) | `hrm_*`, `dpcr_*`, `sanger_*` tables | Detection floors, partition arithmetic and class resolution reproduce the printed tables | That a design will reach those floors in your lab |

## Published outcomes: design → working assay

`inst/extdata/conversion_rates.csv` records six studies that report both a
denominator and its failures, for example 730/1000 (Grewal 2020) and 335/698
(Grewal 2022). Every row keeps `n_fail_pcr`, `n_fail_no_allele` and
`n_fail_other` separate, every row is `calibrated = FALSE`, and no two studies
are pooled into one rate. "Success" means a different thing in each paper, so
the definition travels with the count.

`tests/testthat/test-z8-outcomes.R` enforces that rule. It exists because of a
review point: published markers are successes, so a pass/fail dataset built from
them is survivorship-biased. The response is narrow — keep the denominators,
keep the failures, never present the result as a probability.

Per-family outcome tables do the same for wheat KASP platforms
(`test-kasp-platforms.R`), rice KASP validation, tetra-primer ARMS in wheat and
sorghum, HRM, LAMP, RPA, TaqMan, dPCR, SNaPshot and Sanger.

## Evidence by marker family

Counts are ledger rows whose consumer is a test.

| Area | Rows | Representative anchors | Tests |
|---|---|---|---|
| LAMP / AS-LAMP | 33 | time-to-positive windows, terminal-mismatch rules, B2 dinucleotide | `test-lamp-*` |
| Readout physics | 18 | HRM partition counts, ddPCR LoB/LoD, Sanger VAF floor | `test-readout-physics.R` |
| Tetra-primer ARMS | 16 | wheat 2023 funnel (39 − 7 = 32), sorghum 2021 mismatch positions verified from sequence, blocked barley record, cycle-decay negative result | `test-arms2-*.R` |
| KASP gates, platforms and priors | 13 + 8 + 2 | wheat 2019, Sharplab 2021, maize 2020 design→work outcomes; Kwok 1990 mismatch severity; PolyMarker two axes | `test-kasp-platforms.R`, `test-z7-priors.R`, `test-y32-homeolog.R` |
| Design → working assay (published outcomes) | 12 | Makhoul 2020 strategies, Grewal 2020/2022 and Ramirez-Gonzalez 2014/2015 conversion with denominators, De Silva 2009 HRM, CAPS/LAMP/dPCR funnels | `test-z8-outcomes.R` |
| HRM | 12 + 4 | jsbbs2018 6 × 16 context grid, IJMS 2022 partition arithmetic, Liew 2004 bands | `test-hrmplant-hrm.R`, `test-hrm-jsbbs2018.R` |
| dPCR | 8 | Poisson residuals (Dong), E. coli MPS dilution series | `test-dpcr-*.R` |
| TaqMan / probes | 7 + 6 | wheat 2019 concordance and composition, rice 2018/2024 rates, SARS 2022, C. elegans 2018 | `test-taqman-*`, `test-probes-*` |
| SNaPshot / typing | 10 | 14-locus × 5-input grid, stutter model | `test-typing-*.R` |
| RPA / AS-RPA | 10 | mismatch-position and terminal-mismatch claims | `test-rpa-sources.R` |
| CRISPR-Dx | 4 | Cas12a single and adjacent-double mismatch ratios | `test-crispr-communbiol.R` |
| Sanger | 4 | measured VAF floor for a named enrichment-plus-caller pipeline | `test-sanger-bda-braf.R` |
| Polyploid dosage | 4 | cotton negative control, durum, sugarcane | `test-polyploid-*.R` |
| Tm oracles | 2 | 4096 measured duplex Tm; HRM amplicon Tm | `test-tm-mm2020.R`, `test-tm-hm2025.R` |

The rows above account for 173 of the 174 `test` rows; the last one is
`test-scientific-data.R`, a data-integrity gate rather than a marker family.

## Reproducing

The fast suite (`Rscript dev/test.R`) runs the core contracts and skips the
external cross-checks. The full suite runs them:

```sh
Rscript dev/test.R --full          # or TIDYMARKERS_FULL_TESTS=1
```

`R CMD check --as-cran` sets `_R_CHECK_PACKAGE_NAME_`, so both CI
(`.github/workflows/R-CMD-check.yaml`) and `Rscript dev/check.R --full` exercise
the full suite. Engines that are not installed (BLAST+, `rprimer`, `oligotm`,
annotation databases) skip with a named reason rather than pass silently.

The evidence tables themselves are regenerated from the fetched papers by
`data-raw/*.R`; each generator `stop()`s if a printed number disagrees with the
source, and the tests re-assert the arithmetic on the installed CSV so a
hand-edited table fails without the papers on disk.

## Known gaps

- No bench validation, and none claimed in the README.
- No head-to-head design benchmark (see above).
- 29 ledger rows are `blocked`: read, extracted, consumed by nothing.
- 6 sources are closed access; they exist so that the number they would have
  anchored is visibly absent rather than silently invented.
- `DECIPHER`'s amplification efficiency is a model of hybridisation plus Taq
  extension. It is a proxy for allelic discrimination, not a replacement for an
  experiment.
