# Validation and evidence

What in this package is checked against something outside itself, how that
check is enforced, and what has not been validated at all.

## Summary

- **The maintainers have produced no wet-lab data.** Every measured number in
  this repository was measured by the authors of the cited paper, not here.
  There is no bench validation of a `tidymarkers` design.
- **There is now a head-to-head run against Primer3 2.6.1**, on the layer the
  two tools share: the same oligos through both Tm and dG engines, and the same
  20 SNP loci through both designers with Primer3 constrained to the package's
  own amplicon window. The numbers are in *Primer3, end to end* below, and are
  reproducible with `Rscript dev/primer3-crosscheck.R`.
- **Specificity is scored against ground truth, not against PolyMarker.** On
  TAIR10 the correct answer to "how many places does this primer bind" is
  computable, so `dev/specificity-oracle.R` enumerates it and the fast screens
  are checked against the enumeration, with BLAST+ as a speed/recall baseline.
- **What does exist** is a literature-anchored evidence layer: 80 registered
  sources, 130 extracted tables in `inst/extdata/`, and a ledger that binds 196
  extracts to the 44 test files that consume them. Each link from a source, to
  the exact table read from it, to the test that asserts a property of it is
  machine-checked.
- **Agreement between models is not evidence about a real duplex.** That
  sentence is the rule the thermodynamic tests are written to respect.

## How the evidence is organised

Three tables carry the provenance, and three tests keep them honest.

| File | Rows | What it holds |
|---|---|---|
| `inst/extdata/scientific_sources.csv` | 80 | One row per paper or software source: citation, DOI/PMIM, **how much was read** (`access`), `what_it_anchors`, and the verbatim `claim` |
| `inst/extdata/science_ledger.csv` | 286 | The join `source_key → extract → test_file`, plus `verdict_moved`: which design or ranking decision the extract actually changed |
| `inst/extdata/scientific_backlog.csv` | 45 | Search results and their status: 35 `extracted`, 6 `not_extractable`, 4 `unfetched` |

`access` is a closed vocabulary and it is load-bearing: 54 sources were read in
full text, 11 through a supplement, 3 as an abstract only, 6 are closed access,
2 are source code, 2 are stubs, 1 was transcribed, 1 was not found. 77 of the 80
rows are marked `verified`, each with the date it was checked. A source read only
as an abstract may not calibrate anything, and the tests enforce that.

Ledger rows by role: 196 `test`, 78 `report` (a claim the paper makes, kept for
comparison), 9 `blocked` (read, extracted, consumed by nothing yet — each says
why), 2 `rule`, 1 `artefact`. Every table in `inst/extdata/` that carries a
literature claim now has at least one test reading it; the nine remaining
`blocked` rows name no extract (`extract = none`) and say which paper could not
be landed.

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
installed. `oligotm` is a separate implementation rather than another parameter
set, so it has its own tolerance table (`tm_primer3_tolerances()`): both sides
are pinned to the same physical model — 50 mM monovalent, no Mg²⁺, no dNTPs and
the SantaLucia (1998) salt correction (`oligotm -tp 1 -sc 1 -mv 50 -dv 0 -n 0`)
— and only arithmetic remains. On the 32-oligo oracle set, TmCalculator's
Allawi (1998) table and `-tp 1` agree to the printed precision (max
|ΔTm| 4.8 × 10⁻⁷ °C), while the package default (SantaLucia 2004) differs by
0.40 °C on the 20-mers and 0.98 °C on the degenerate poly-A probe; that probe
sets the 1.0 °C cap. Pinning the flags matters because `oligotm`'s `-n` is the
dNTP concentration, not the salt concentration: a run that leaves them at their
defaults folds in a divalent term this model does not have. Structural anchors
carry no tolerance at all (Tm rises with GC at fixed length, with length at
fixed pattern, and is symmetric under reverse complement). The three-way
disagreement is reported, never averaged, and never asserted against a single
constant. `dev/compare-thermo.R` prints the same matrix for inspection.

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

### Primer3, end to end

`dev/primer3-crosscheck.R` runs Primer3 2.6.1 beside this package on one machine
and prints, rather than asserts, how far apart they are. `tests/testthat/test-primer3-crosscheck.R`
asserts only the bands that run supports.

**Thermodynamics, same oligos, same salt and temperature** (34 oligos: 24 seeded
20-mers, 8 GC/homopolymer edge probes, the 2 universal KASP tails):

| Quantity | Oligotm / ntthal − this package |
|---|---|
| Tm, matched salt (Na 50 mM, Mg 1.5 mM, dNTPs 0.6 mM, 37 °C, 50 nM) | mean **+3.92 °C**, sd 0.25, max 4.66; r = 1.000, ρ = 0.999 |
| Tm, package called at Mg = 0 | mean **+12.48 °C** — the package's divalent term is worth ≈8.5 °C here |
| Hairpin dG (OligoArrayAux vs `ntthal -a HAIRPIN`) | mean **−0.84 kcal/mol**, sd 0.76, range [−2.86, 0.27], ρ = 0.67 |
| Heterodimer dG (`ntthal -a ANY`) | n = 4 pairs, mean **−1.75 kcal/mol**, sd 1.16, ρ = 1.000 |

The ~4 °C Tm residual is the nearest-neighbour table (SantaLucia 1998 on
Primer3's side, 2004 in `primer_tm()`) plus the salt formula (Primer3 defaults
to SantaLucia 1998, the package to Schildkraut 2010). The crosscheck reports the
total and does not apportion it. The oracle test above does, because it runs
both engines on one salt model: at 50 mM monovalent with the same correction
the two tables differ by at most 1.0 °C (0.40 °C on the 20-mers), so at most
that much of the residual is the table and the rest is the salt model. The
offset is recorded, not corrected. `ntthal` reports no hairpin at all on 16 of
the 34 oligos; the package reports dG = 0 on the 2 that genuinely cannot fold,
which is what the same statement means.

**Design, same 20 loci.** The panel is the first 20 TAIR10-matching SNVs of the
cached Chr4 variation window, thinned to ≥1 kb apart (the cache region carries
one SNV per ~10 bp, and `design_kasp_markers()` passes its own `loci` table as
the `variants` table for footprint accounting). Primer3 is run twice per locus:
with its usual 60–150 bp product range, and constrained to the package's own
50–100 bp window. The head-to-head uses the constrained run, so a length
disagreement cannot pose as a quality disagreement.

| Result | Count |
|---|---|
| Primer3 returned ≥1 pair in the 50–100 bp window | 17/20 |
| …and that pair cleared all seven package gates (GC, Tm, amplicon, homopolymer, 3′-placement, hairpin, dimer) | 8/20 |
| …and this package also designed an assay | **5/8** (the agreement floor) |
| Primer3 returned no pair at all | 3/20 |
| Primer3 pair rejected by a package gate | 9/20 |
| Primer3 would have used a >100 bp product on its own 60–150 bp range | 12/20 |
| Amplicon size, both solved (n = 7) | median difference **+1 bp**, range [−5, +25], ρ = 0.847 |
| Runtime per locus, single core | Primer3 0.06 s, this package **9.5 s** |

The three loci where Primer3 had a usable pair and this package returned none
fail on their own stated rules (`No AS primer survived even the relaxed
constraints`, `No common primer passed the filters`, `No AS primer passed the
GC/Tm filter`), not on an unstated preference. What the package adds after that
layer — allele-specific discrimination, 3′-mismatch class, off-target counting —
is exactly what Primer3 does not do, so this table is the *shared* layer and not
a ranking of the two tools.

### Specificity against ground truth

`dev/specificity-oracle.R` hits TAIR10 (119.7 Mb) three ways. PolyMarker is not
an oracle: it is another implementation with its own heuristics, so the right
answer is enumerated instead.

| Check | Ground truth | Found | Time |
|---|---|---|---|
| 3′-seed screen (15 nt, ≤1 substitution, both strands), 8 primers | 312, 12, 43, 37, 73, 91, 9, 18 hits | **8/8 equal** | 4.7 s vs 169 s to enumerate |
| Exact in-silico PCR (`max_mismatch = 0`), 3 genomic primer pairs | 1, 2, 1 products | **3/3 equal** | 10.6 s |
| `blastn -task blastn-short` (word 7) as a baseline | 8 full-length exact hits | recall **1.000** | 0.47 s |

BLAST is faster than the seed screen on this panel, but it answers a narrower
question (exact full-length hits) and returns no per-primer mismatch profile.
The seed screen's value is the 3′-seed definition itself: it is the quantity the
specificity gates are written in, and it is now pinned to enumeration.

The third table above is not a benchmark of `blast_primers()`; that path needs
`rBLAST` and an index, and is covered by its own test.

### Two bugs the cross-checks found

Both were found by running an independent answer, and both are fixed with a
test that fails without the fix:

1. **UNAFold's "no structure" sentinel** was passed through as a dG: a poly-A
   primer carried **+999.999 kcal/mol** into any mean or sum over primers.
   `oligo_dg()` now maps the ~+1000 sentinel to 0 (no stabilisation relative to
   the unfolded state), and `tests/testthat/test-oligo-dg.R` asserts it.
2. **`in_silico_pcr()` paired hits from different chromosomes.** The whole
   genome was scanned inside a per-chromosome loop and the hits were never
   filtered by chromosome, so one real product came back as twelve phantom
   ones — at 7× the scan cost. Found by enumerating the same answer. Fixed, with
   `tests/testthat/test-insilico-pcr.R` pinning the cross-chromosome case.

### What is still missing

Wet-lab validation of a `tidymarkers` design: none exists, and nothing in the
tables above substitutes for it. A *predictor* of published per-assay success is
also not possible from this corpus yet: the tables that carry published primer
triples do not carry both outcomes. The maize 2020 table, for instance, has 100
rows with primer triples and `design_success = TRUE` on all 100, with no
`assay_outcome`; a confusion matrix needs the failures too. The denominator-level
outcome tables (Grewal 2020/2022, Ramirez-Gonzalez 2014/2015, the wheat and
rice KASP conversion studies) are consumed by `tests/testthat/test-z8-outcomes.R`
and are the right place to start once per-assay outcomes are landed.

## Measured versus predicted

| Anchor | Data | What is asserted | What is not |
|---|---|---|---|
| Oliveira et al., *Chem Sci* 2020 (`10.1039/d0sc01700k`) | 4096 measured 19-mer duplex Tm (`tm_mm2020_s1_measured.csv`) | Ordering and discrimination: model and measurement rank the 64 canonical duplexes alike (r = 0.958) and the model separates a perfect duplex from a triple mismatch; that the artifact itself records the stance (`calibration_stance` on every row) and forbids a calibration | Absolute agreement. `primer_tm()` sits ~16.9 °C below the measurements under that buffer (sd 0.64 °C); the offset is recorded, not corrected away. The candidate explanations (salt term, measurement standard; Mg ruled out) and their arithmetic are documented in `data-raw/mm2020.R`, and the ~9.3 °C left after the more complete salt model stays unexplained |
| *PLoS ONE* 2025 (`10.1371/journal.pone.0321885`) | 5 amplicon Tm, psaA-2 | The table's printed arithmetic; the authors' ordering claim | A test of `primer_tm()` on a primer duplex — this is an amplicon measurement |
| Soybean worked panel (local design document, no DOI) | 85 ordered oligos with the Tm printed by the source design pipeline (`soybean_worked_panel_primers.csv`) | That `primer_tm()` reproduces a real design table: SantaLucia ↔ Allawi max \|ΔTm\| 0.079 °C, and the `SantaLucia1998-2` salt correction reproduces the printed Tm to 0.12 °C. Also that the 30 printed exact-genome-hit counts are recomputed from GCF_000004515.6 (exact full-length body matches, both strands, 20 chromosomes) and agree, each hit reported with its chromosome, `chrom:start-end` coordinates and strand (`genome_hit_coord`, `genome_hit_strand`); that each of the 23 Sanger pairs amplifies exactly one exact product at the printed coordinates (confirmed with `in_silico_pcr()`), since both primers are genome-wide unique; and that the off-target check found one ALT primer whose only match is off-target (chr 1, not chr 18) | That the Tm is thermodynamically true — this is an external *computation*, not a measured duplex; the hit and product checks are exact full-length only (near-matches are a separate question); the source text (`.tmp/`) and the 946 MB genome (`marker_cache_dir()`) are not in git, so both reproduction tests skip in CI without them |
| Soybean worked panel — KASP quality gates (same local document) | 13 KASP sets × the package's 27 `snp` acceptance checks (`soybean_worked_panel_kasp_gates.csv`) | That each set is graded by `marker_design_checks()`, not by the document's prose: the amplicon span, in-silico PCR and off-target counts are recomputed from GCF_000004515.6; the AS on/off windows give the discrimination; Tm uses the `SantaLucia1998-2` model the printed Tm pins. Off-target and in-silico PCR pass all 13; the measured failures are recorded (11 `COM Tm >= AS Tm`, 5 discrimination < 0.5, 4 amplicons > 100 bp, 4 AS GC, 4 homopolymer, 3 SSR, 3 inverted repeat). The six footprint/dropout checks and the two linkage checks are `unknown` because the panel carries no population variant table | That a design gate guarantees a working assay; the footprint/dropout and linkage checks cannot be evaluated without population variants; the genome is not in git, so the reproduction test skips in CI |
| Readout physics (HRM, ddPCR, Sanger) | `hrm_*`, `dpcr_*`, `sanger_*` tables | Detection floors, partition arithmetic and class resolution reproduce the printed tables | That a design will reach those floors in your lab |

## Published outcomes: design → working assay

`inst/extdata/conversion_rates.csv` records six studies that report both a
denominator and its failures, for example 730/1000 (Grewal 2020) and 335/698
(Grewal 2022). Every row keeps `n_fail_pcr`, `n_fail_no_allele` and
`n_fail_other` separate, every row is `calibrated = FALSE`, and no two studies
are pooled into one rate. "Success" means a different thing in each paper, so
the definition travels with the count: `conversion_definition` names which of
`design` / `manufacture` / `work` / `usable` the numerator counts, and `scope`
keeps the assay-level count apart from the assay×sample rate (which lives only
in `wheat_platform2019_platform_outcomes.csv`).

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
| LAMP / AS-LAMP | 37 | time-to-positive windows, terminal-mismatch rules, B2 dinucleotide, the ace1r contradiction/rule/metric extracts | `test-lamp-*` |
| Readout physics | 21 | HRM partition counts, ddPCR LoB/LoD, Sanger VAF floor, JAK2 CV reproduction, HRM minimum input | `test-readout-physics.R` |
| Tetra-primer ARMS | 17 | wheat 2023 funnel (39 − 7 = 32), sorghum 2021 mismatch positions verified from sequence, blocked barley record, cycle-decay negative result, Liu 2012's 1686 published ASP primers | `test-arms2-*.R`, `test-arms-mismatch-position.R` |
| KASP gates, platforms and priors | 25 | wheat 2019, Sharplab 2021, maize 2020 design→work outcomes; Kwok 1990 mismatch severity; PolyMarker two axes; LefeVER's three rules and the 4-cell mismatch panel | `test-kasp-platforms.R`, `test-z7-priors.R`, `test-y32-homeolog.R`, `test-kasp-mismatch-type.R` |
| Design → working assay (published outcomes) | 15 | Makhoul 2020 strategies, Grewal 2020/2022 and Ramirez-Gonzalez 2014/2015 conversion with denominators, De Silva 2009 HRM amplicons/SNPs/strategies, CAPS/LAMP/dPCR funnels | `test-z8-outcomes.R` |
| HRM | 19 | jsbbs2018 6 × 16 context grid and its summary/recomputed split, IJMS 2022 dTm controls and partition arithmetic, Liew 2004 bands | `test-hrm-curve.R`, `test-hrm-jsbbs2018.R` |
| dPCR | 8 | Poisson residuals (Dong), E. coli MPS dilution series | `test-dpcr-*.R` |
| TaqMan / probes | 13 | wheat 2019 concordance and composition, rice 2018/2024 rates, SARS 2022, C. elegans 2018 | `test-taqman-*`, `test-probes-*` |
| SNaPshot / typing | 10 | 14-locus × 5-input grid, stutter model | `test-typing-*.R` |
| RPA / AS-RPA | 10 | mismatch-position and terminal-mismatch claims | `test-rpa-sources.R` |
| CRISPR-Dx | 4 | Cas12a single and adjacent-double mismatch ratios | `test-crispr-communbiol.R` |
| Sanger | 4 | measured VAF floor for a named enrichment-plus-caller pipeline | `test-sanger-bda-braf.R` |
| Polyploid dosage | 4 | cotton negative control, durum, sugarcane | `test-polyploid-*.R` |
| Tm oracles | 2 | 4096 measured duplex Tm; HRM amplicon Tm | `test-tm-mm2020.R`, `test-tm-hm2025.R` |
| CAPS | 4 | 1236 palindromic enzymes with their cut offsets, 26 soybean primer pairs, SNP2CAPS predicted-vs-observed patterns, Liu 2016's within-family SNP composition | `test-caps.R` |
| Registry and IO | 3 | the 80-row source registry, the 45-entry backlog's three statuses, the example locus table's transition chemistry | `test-scientific-data.R`, `test-lamp-status.R`, `test-io.R` |

These are all 196 `test` rows: each row is assigned to the first pattern in the
last column that matches its `test_file`, so the counts sum to the ledger total
and no row is counted twice. `test-science-ledger.R` and `test-scientific-data.R`
themselves are the data-integrity gates rather than marker-family tests, and are
counted under *Registry and IO*.

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

The two cross-checks in this document are separate scripts, because they print
disagreement rather than assert it:

```sh
Rscript dev/primer3-crosscheck.R   # Primer3 -> results/primer3-{thermo,design}.csv
Rscript dev/specificity-oracle.R   # TAIR10 ground truth -> results/specificity-oracle.csv
```

Both need binaries that are not package dependencies (`oligotm`, `ntthal`,
`primer3_core`, `blastn`, `makeblastdb`) and the OligoArrayAux build; the
matching tests skip with a named reason when one is absent.

## Known gaps

- No bench validation, and none claimed in the README.
- The Primer3 comparison covers the *shared* layer only: Primer3 does not design
  allele-specific primers, so the allelic-discrimination layer has no external
  comparator at all.
- No predictor of published per-assay success is fitted: the corpus has
  denominators at study level, but the tables that carry primer triples do not
  carry both outcomes (maize 2020 has 100/100 `design_success = TRUE` and no
  `assay_outcome`).
- 9 ledger rows are `blocked`: papers that could not be landed, each naming
  what it would have anchored.
- 6 sources are closed access; they exist so that the number they would have
  anchored is visibly absent rather than silently invented.
- `DECIPHER`'s amplification efficiency is a model of hybridisation plus Taq
  extension. It is a proxy for allelic discrimination, not a replacement for an
  experiment.
