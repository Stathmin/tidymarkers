# Validation and evidence

What in this package is checked against something outside itself, how that
check is enforced, and what has not been validated at all.

## Summary

- **The maintainers have produced no wet-lab data.** No number that comes from
  a paper was measured here. What the repository measures itself — engine
  agreement, the specificity oracle, the validation benchmark — is marked
  `measured` below and is reproducibility evidence, not experimental evidence.
  There is no bench validation of a `tidymarkers` design.
- **There is now a head-to-head run against Primer3 2.6.1**, on the layer the
  two tools share: the same oligos through both Tm and dG engines, and the same
  20 SNP loci through both designers with Primer3 constrained to the package's
  own amplicon window. The numbers are in _Primer3, end to end_ below, and are
  reproducible with `Rscript dev/primer3-crosscheck.R`.
- **The allele-specific layer has one external comparator.** Primer3 does not
  design allele-specific primers, but Biopython's `Tm_NN` scores a 3′-terminal
  mismatch through its own penalty table. On the same eight panel rows both
  engines put the mutant duplex lower, and within every primer the penultimate
  mismatch costs more than the terminal one; the two metrics rank the rows
  differently (ρ = 0.238), so only direction is asserted. See _Allele-specific
  discrimination_ below.
- **Specificity is scored against ground truth, not against PolyMarker.** On
  TAIR10 the correct answer to "how many places does this primer bind" is
  computable, so `dev/specificity-oracle.R` enumerates it and the fast screens
  are checked against the enumeration, with BLAST+ as a speed/recall baseline.
- **What does exist** is a literature-anchored evidence layer: 134 registered
  sources, 149 extracted tables in `inst/extdata/`, and a ledger that binds 209
  `test` rows to the 46 test files that consume them. Each link from a source,
  to the exact table read from it, to the test that asserts a property of it is
  machine-checked.
- **Agreement between models is not evidence about a real duplex.** That
  sentence is the rule the thermodynamic tests are written to respect.

## How the evidence is organised

Three tables carry the provenance, and three tests keep them honest.

| File                                  | Rows | What it holds                                                                                                                           |
| ------------------------------------- | ---- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `inst/extdata/scientific_sources.csv` | 134  | One row per paper or software source: citation, DOI/PMIM, **how much was read** (`access`), `what_it_anchors`, and the verbatim `claim` |
| `inst/extdata/science_ledger.csv`     | 369  | The join `source_key → extract → test_file`, plus `verdict_moved`: which design or ranking decision the extract actually changed        |
| `inst/extdata/scientific_backlog.csv` | 45   | Search results and their status: 35 `extracted`, 6 `not_extractable`, 4 `unfetched`                                                     |

`access` is a closed vocabulary and it is load-bearing: 104 sources were read in
full text, 11 through a supplement, 4 as an abstract only, 6 are closed access,
3 are source code, 2 are stubs, 1 was transcribed, 3 were not found. 128 of the
134 rows are marked `verified`, each with the date it was checked. A source read
only as an abstract may not calibrate anything, and the tests enforce that.

Ledger rows by use: 209 `test`, 132 `report` (a claim the paper makes, kept for
comparison), 19 `blocked` (read, extracted, consumed by nothing yet — each says
why), 3 `rule`, 4 `artefact`, 2 `config`. Across the ledger 146 extract tables are named — `example_region_targets.csv`
is the shipped slice's target table, `steele2025_rice_kasp_calls.csv` the
per-marker matrix above, and the 149th CSV in `inst/extdata/` is the ledger
itself — and `test` rows account for 130 of them.
Of the 19 `blocked` rows, 9 name no extract (`extract = none`) and say which
paper could not be landed, and 10 name a real extract that carries no verdict of
its own yet.

The gates:

| Test                                    | What it refuses                                                                                                                                                                  |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tests/testthat/test-science-ledger.R`  | A ledger row naming a test that does not exist or does not actually read the named extract; an extract whose own DOI belongs to a different paper; a source with no ledger entry |
| `tests/testthat/test-scientific-data.R` | A table with a numeric outcome and no source, no `access` level, or a `calibrated = TRUE` flag on something read only as an abstract                                             |
| `tests/testthat/test-falsifiable.R`     | A published acceptance criterion that cannot be driven to `fail` by the suite — every criterion in the registry needs a mutation that trips it, or a recorded reason why not     |

## Cross-tool comparison

### Thermodynamics: same lineage, four parameter sets, three engines

`primer_tm()` is the Primer3/`oligotm` nearest-neighbour family (TmCalculator,
SantaLucia 2004). That lineage is shared rather than compared, so the package
tests the thing that is actually uncertain: how much the _window_ wobbles
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

| Quantity                                                           | Oligotm / ntthal − this package                                        |
| ------------------------------------------------------------------ | ---------------------------------------------------------------------- |
| Tm, matched salt (Na 50 mM, Mg 1.5 mM, dNTPs 0.6 mM, 37 °C, 50 nM) | mean **+3.92 °C**, sd 0.25, max 4.66; r = 1.000, ρ = 0.999             |
| Tm, package called at Mg = 0                                       | mean **+12.48 °C** — the package's divalent term is worth ≈8.5 °C here |
| Hairpin dG (OligoArrayAux vs `ntthal -a HAIRPIN`)                  | mean **−0.84 kcal/mol**, sd 0.76, range [−2.86, 0.27], ρ = 0.67        |
| Heterodimer dG (`ntthal -a ANY`)                                   | n = 4 pairs, mean **−1.75 kcal/mol**, sd 1.16, ρ = 1.000               |

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

| Result                                                                                                       | Count                                                   |
| ------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------- |
| Primer3 returned ≥1 pair in the 50–100 bp window                                                             | 17/20                                                   |
| …and that pair cleared all seven package gates (GC, Tm, amplicon, homopolymer, 3′-placement, hairpin, dimer) | 8/20                                                    |
| …and this package also designed an assay                                                                     | **5/8** (the agreement floor)                           |
| Primer3 returned no pair at all                                                                              | 3/20                                                    |
| Primer3 pair rejected by a package gate                                                                      | 9/20                                                    |
| Primer3 would have used a >100 bp product on its own 60–150 bp range                                         | 12/20                                                   |
| Amplicon size, both solved (n = 7)                                                                           | median difference **+1 bp**, range [−5, +17], ρ = 0.883 |
| Runtime per locus, single core                                                                               | Primer3 0.05 s, this package **6.1 s**                  |

The amplicon-size row's `n = 7` is the set of loci where _both_ tools produced
any product — Primer3 on its own 60–150 bp range and this package on its
50–100 bp gate. That set overlaps the eight gate-clearing loci above; it is not
a subset of them.

Both figures are wall-clock and move with the host and with designer changes, so
they are recorded with the tree they were measured on: the first crosscheck run
read ρ = 0.847 over [−5, +25] and 9.5 s per locus; the current values come from a
run after `263cb50` (hairpin dG in R), `0383fe3`/`7529c4a` (screening) and
`98502f5` (canonical locus row). Exactly one package amplicon moved across those
commits — 81 → 73 bp at `ENSVATH10569559` (Chr4:1513124) — which turns its
difference to Primer3 from +25 to +17 and leaves the Primer3 side unchanged.
That single row is the whole of the drift; every acceptance band still holds.

The three loci where Primer3 had a usable pair and this package returned none
fail on their own stated rules (`No AS primer survived even the relaxed
constraints`, `No common primer passed the filters`, `No AS primer passed the
GC/Tm filter`), not on an unstated preference. What the package adds after that
layer — allele-specific discrimination, 3′-mismatch class, off-target counting —
is exactly what Primer3 does not do, so this table is the _shared_ layer and not
a ranking of the two tools.

### Allele-specific discrimination: a terminal-mismatch Tm cross-check

Primer3 does not design allele-specific primers, so the Primer3 run above cannot
score the discrimination layer. The one external engine that can is Biopython's
`Bio.SeqUtils.MeltingTemp.Tm_NN` (1.85), which applies a terminal-mismatch
penalty table (`DNA_TMM1`) on top of its nearest-neighbour set. `ntthal` cannot:
it truncates the primer at the mismatch instead of scoring it, and three
different 3′ substitutions returned the identical dG = −4333.7 with it.

`dev/as-tm-crosscheck.R` builds an 8-row panel — four 6-nt allele-specific
regions, each scored with the mismatch at the 3′ terminal base and again at the
penultimate base — computes the package's `discrimination` and the on/off
efficiency gap `dgap`, calls `dev/as-tm-crosscheck.py` (`Tm_NN` at
dnac1 = dnac2 = 25 nM, Na 50 mM, Mg 1.5 mM, dNTPs 0.6 mM, saltcorr 5), and merges
the two.

| primer region | mismatch at | `discrimination` | `dgap` | `tm_ref` | `tm_alt` | `dtm` |
| ------------- | ----------- | ---------------- | ------ | -------- | -------- | ----- |
| CATGCT        | 3′ terminal | 0.0296           | 0.1154 | 58.49    | 57.55    | 0.95  |
| CAGGAT        | 3′ terminal | 0.1618           | 0.8740 | 59.82    | 59.38    | 0.45  |
| CGATGA        | 3′ terminal | 0.4303           | 0.8644 | 58.74    | 57.83    | 0.91  |
| CCGTAA        | 3′ terminal | 0.3243           | 0.8695 | 58.55    | 57.97    | 0.58  |
| CATGCT        | penultimate | 0.3542           | 0.6966 | 58.49    | 49.69    | 8.80  |
| CAGGAT        | penultimate | 0.2992           | 0.6035 | 59.82    | 54.50    | 5.33  |
| CGATGA        | penultimate | 0.3175           | 0.6324 | 58.74    | 51.89    | 6.85  |
| CCGTAA        | penultimate | 0.3029           | 0.6021 | 58.55    | 52.75    | 5.80  |

All eight rows point the same way in both engines: the mismatch destabilises the
duplex (`dtm > 0`) and lowers amplification (`dgap > 0`), and within every primer
the penultimate mismatch costs more than the terminal one (median `dtm` 6.33 °C
against 0.74 °C). Beyond that the two metrics do not agree. Spearman ρ(`dtm`,
`discrimination`) is 0.238, and `discrimination` inverts the terminal/penultimate
ordering for two of the four primers. The test therefore asserts the direction
and the within-primer ordering, reports ρ, and asserts no magnitude: a Tm drop is
not an efficiency drop, and the second depends on where the mismatch sits with
respect to the polymerase, which a two-state nearest-neighbour model does not
contain.

Routes that were tried and do not work, recorded so a later reader does not
retry them: the Bioconductor package `melting` cannot be installed here (`curl`
to bioconductor.org times out, exit 28); `nupack` has no PyPI distribution;
`pydna` models perfect duplexes only; and `primer3-py` is the same engine as the
cross-check above rather than an independent one.

### Specificity against ground truth

`dev/specificity-oracle.R` hits TAIR10 (119.7 Mb) three ways. PolyMarker is not
an oracle: it is another implementation with its own heuristics, so the right
answer is enumerated instead.

| Check                                                            | Ground truth                        | Found            | Time                        |
| ---------------------------------------------------------------- | ----------------------------------- | ---------------- | --------------------------- |
| 3′-seed screen (15 nt, ≤1 substitution, both strands), 8 primers | 312, 12, 43, 37, 73, 91, 9, 18 hits | **8/8 equal**    | 4.7 s vs 169 s to enumerate |
| Exact in-silico PCR (`max_mismatch = 0`), 3 genomic primer pairs | 1, 2, 1 products                    | **3/3 equal**    | 10.6 s                      |
| `blastn -task blastn-short` (word 7) as a baseline               | 8 full-length exact hits            | recall **1.000** | 0.47 s                      |

BLAST is faster than the seed screen on this panel, but it answers a narrower
question (exact full-length hits) and returns no per-primer mismatch profile.
The seed screen's value is the 3′-seed definition itself: it is the quantity the
specificity gates are written in, and it is now pinned to enumeration.

The third table above is not a benchmark of `blast_primers()`; that path needs
`rBLAST` and an index, and is covered by its own test.

### One row per external tool

`dev/cross-tool-sweep.R` answers ZY20's question — one row per tool, numbers or an
explicit `not_run` — and writes `.tmp/cross-tool-sweep.csv`. It recomputes the
Primer3, UNAFold, BLAST and Biopython rows from the CSVs the scripts above produced, and
probes the host for every other binary in the same run, so each cause is a fact
observed then, not a note copied from an earlier session. The Primer3 block is a
regression anchor: it re-checks the bands this file records (Tm r ≥ 0.999, mean
+3.92 °C, 17/20 loci with a pair, 5/8 agreement) and exits non-zero if one slips,
because a slip means the harness changed, not the world.

| Tool              | Kind                     | Status                     | Result on this host                                                                                                                                                                |
| ----------------- | ------------------------ | -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Primer3 2.6.1     | designer + thermo oracle | measured                   | 17/20 loci get a pair in the 50–100 bp window; amplicon median +1 bp (n = 7), ρ = 0.883; 3/20 no pair; 9/20 rejected by a package gate                                             |
| UNAFold / mfold   | thermo oracle (dG)       | measured (not independent) | hairpin dG mean −0.84 kcal/mol (n = 34); heterodimer −1.75 (n = 4); `hybrid-ss-min` is the engine `oligo_dg()` calls, so this checks the wrapper, not the method                   |
| BLAST+ (`blastn`) | specificity baseline     | measured                   | seed 8/8 equal, exact PCR 3/3 equal, 8 exact full-length hits, 0.47 s; answers a narrower question than the seed screen and is not a designer                                      |
| Biopython `Tm_NN` | allele-specific oracle   | measured                   | 8/8 panel rows point the same way as the package's efficiency gap; median `dtm` 6.33 °C penultimate vs 0.74 °C terminal; ρ(`dtm`, `discrimination`) = 0.238, reported not asserted |
| Primer-BLAST      | designer                 | not_run                    | web service only — needs network and manual submission; no local binary to run unattended                                                                                          |
| BatchPrimer3      | designer                 | not_run                    | same: web form only                                                                                                                                                                |
| MFEprimer-3.0     | designer + specificity   | not_run                    | `Sys.which("mfeprimer") == ""` — no binary in `PATH`                                                                                                                               |
| isPcr (UCSC)      | specificity              | not_run                    | `Sys.which("ispcr") == ""` — `blastn` is not a substitute for its product model                                                                                                    |
| EasyKASP          | designer (KASP)          | not_run                    | `Sys.which("easykasp") == ""`                                                                                                                                                      |
| NUPACK            | thermo oracle            | not_run                    | no binary; licence-gated                                                                                                                                                           |
| MELTING           | thermo oracle            | not_run                    | `Sys.which("melting") == ""`                                                                                                                                                       |
| VCF2CAPS          | designer                 | not_run                    | publication code, not distributed as an installable tool here; source is `vcf2caps2021`                                                                                            |
| KASP-IEva         | designer                 | not_run                    | publication code, not installable here; no source row yet                                                                                                                          |
| LAMPrimers iQ     | designer                 | not_run                    | publication code, not installable here; no source row yet                                                                                                                          |
| ThermoPlex        | designer                 | not_run                    | publication code, not installable here; no source row yet                                                                                                                          |
| ARMSprimer3       | designer                 | not_run                    | publication code, not installable here; no source row yet                                                                                                                          |
| PrimerWeaver      | designer                 | not_run                    | publication code, not installable here; no source row yet                                                                                                                          |
| PyPCRtool         | designer                 | not_run                    | publication code, not installable here; no source row yet                                                                                                                          |
| PRISM             | designer                 | not_run                    | publication code, not installable here; no source row yet                                                                                                                          |

Four of nineteen run here. The fifteen that do not are named with the cause the
run measured — a missing binary, a web-only form, a licence, or code that is not
distributed — and none of them is reported as a comparison. The return gate is
unchanged: install a named tool, the sweep reports it as newly present, then its
block is written against the same anchors.

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
tables above substitutes for it. A _predictor_ of published per-assay success is
also not possible from this corpus yet: the tables that carry published primer
triples do not carry both outcomes. The maize 2020 table, for instance, has 100
rows with primer triples and `design_success = TRUE` on all 100, with no
`assay_outcome`; a confusion matrix needs the failures too. The denominator-level
outcome tables (Grewal 2020/2022, Ramirez-Gonzalez 2014/2015, the wheat and
rice KASP conversion studies) are consumed by `tests/testthat/test-z8-outcomes.R`
and are the right place to start once per-assay outcomes are landed.

## Measured versus predicted

| Anchor                                                          | Data                                                                                                     | What is asserted                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | What is not                                                                                                                                                                                                                                                                                                                                                         |
| --------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Oliveira et al., _Chem Sci_ 2020 (`10.1039/d0sc01700k`)         | 4096 measured 19-mer duplex Tm (`tm_mm2020_s1_measured.csv`)                                             | Ordering and discrimination: model and measurement rank the 64 canonical duplexes alike (r = 0.958) and the model separates a perfect duplex from a triple mismatch; that the artifact itself records the stance (`calibration_stance` on every row) and forbids a calibration                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | Absolute agreement. `primer_tm()` sits ~16.9 °C below the measurements under that buffer (sd 0.64 °C); the offset is recorded, not corrected away. The candidate explanations (salt term, measurement standard; Mg ruled out) and their arithmetic are documented in `data-raw/mm2020.R`, and the ~9.3 °C left after the more complete salt model stays unexplained |
| _PLoS ONE_ 2025 (`10.1371/journal.pone.0321885`)                | 5 amplicon Tm, psaA-2                                                                                    | The table's printed arithmetic; the authors' ordering claim                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | A test of `primer_tm()` on a primer duplex — this is an amplicon measurement                                                                                                                                                                                                                                                                                        |
| Soybean worked panel (local design document, no DOI)            | 85 ordered oligos with the Tm printed by the source design pipeline (`soybean_worked_panel_primers.csv`) | That `primer_tm()` reproduces a real design table: SantaLucia ↔ Allawi max \|ΔTm\| 0.079 °C, and the `SantaLucia1998-2` salt correction reproduces the printed Tm to 0.12 °C. Also that the 30 printed exact-genome-hit counts are recomputed from GCF_000004515.6 (exact full-length body matches, both strands, 20 chromosomes) and agree, each hit reported with its chromosome, `chrom:start-end` coordinates and strand (`genome_hit_coord`, `genome_hit_strand`); that each of the 23 Sanger pairs amplifies exactly one exact product at the printed coordinates (confirmed with `in_silico_pcr()`), since both primers are genome-wide unique; and that the off-target check found one ALT primer whose only match is off-target (chr 1, not chr 18) | That the Tm is thermodynamically true — this is an external _computation_, not a measured duplex; the hit and product checks are exact full-length only (near-matches are a separate question); the source text (`.tmp/`) and the 946 MB genome (`marker_cache_dir()`) are not in git, so both reproduction tests skip in CI without them                           |
| Soybean worked panel — KASP quality gates (same local document) | 13 KASP sets × the package's 27 `snp` acceptance checks (`soybean_worked_panel_kasp_gates.csv`)          | That each set is graded by `marker_design_checks()`, not by the document's prose: the amplicon span, in-silico PCR and off-target counts are recomputed from GCF_000004515.6; the AS on/off windows give the discrimination; Tm uses the `SantaLucia1998-2` model the printed Tm pins. Off-target and in-silico PCR pass all 13; the measured failures are recorded (11 `COM Tm >= AS Tm`, 5 discrimination < 0.5, 4 amplicons > 100 bp, 4 AS GC, 4 homopolymer, 3 SSR, 3 inverted repeat). The six footprint/dropout checks and the two linkage checks are `unknown` because the panel carries no population variant table                                                                                                                                  | That a design gate guarantees a working assay; the footprint/dropout and linkage checks cannot be evaluated without population variants; the genome is not in git, so the reproduction test skips in CI                                                                                                                                                             |
| Readout physics (HRM, ddPCR, Sanger)                            | `hrm_*`, `dpcr_*`, `sanger_*` tables                                                                     | Detection floors, partition arithmetic and class resolution reproduce the printed tables                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | That a design will reach those floors in your lab                                                                                                                                                                                                                                                                                                                   |

## Published outcomes: design → working assay

`inst/extdata/conversion_rates.csv` records sixteen rows from twelve publications that report both a
denominator and its failures, for example 730/1000 (Grewal 2020) and 335/698
(Grewal 2022). Every row keeps `n_fail_pcr`, `n_fail_no_allele`,
`n_fail_other` and (for Grewal 2022) the three further named classes
`n_fail_wheat_polymorphic`, `n_fail_monomorphic`, `n_fail_heterozygote`
separate; five rows reconstruct their failure total exactly, every row is
`calibrated = FALSE`, and no two studies
are pooled into one rate. "Success" means a different thing in each paper, so
the definition travels with the count: `conversion_definition` names which of
`design` / `manufacture` / `work` / `usable` the numerator counts, and `scope`
keeps the assay-level count apart from the assay×sample rate (which lives only
in `wheat_platform2019_platform_outcomes.csv`).

`tests/testthat/test-z8-outcomes.R` enforces that rule. It exists because of a
review point: published markers are successes, so a pass/fail dataset built from
them is survivorship-biased. The response is narrow — keep the denominators,
keep the failures, never present the result as a probability.

One file moves from study-level counts to per-marker calls:
`inst/extdata/steele2025_rice_kasp_calls.csv` is the per-marker call matrix of
Steele 2025, 3,366 KASP designs × 178 rice varieties parsed from the paper's
Supplementary File S2 by `data-raw/steele2025_rice_s2.R`. Each row is one marker
partitioned into ref / alt / het / missing / other calls that sum to 178. `?` is
the only no-call (9,810 cells, 1.64 %); `-:-` is a real homozygous deletion, not
a miss; `Uncallable` and `Bad` are annotations rather than calls.
`tests/testthat/test-steele-s2.R` enforces the per-marker split, the measured
panel totals and three hand-checked markers. The table is a published,
measured outcome, not a simulation. It adds no ledger row of its own, so it does
not change the family counts below; the paper's prose (>90 % called for all
3,366 designs, 635 designs failed) does not match the shipped matrix (91 markers
at ≤90 % called, one with no calls), and that mismatch is recorded rather than
adjudicated.

Per-family outcome tables do the same for wheat KASP platforms
(`test-kasp-platforms.R`), rice KASP validation, tetra-primer ARMS in wheat and
sorghum, HRM, LAMP, RPA, TaqMan, dPCR, SNaPshot and Sanger.

## Evidence by marker family

Counts are ledger rows whose consumer is a test.

| Area                                        | Rows | Representative anchors                                                                                                                                                             | Tests                                                                                           |
| ------------------------------------------- | ---- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| LAMP / AS-LAMP                              | 37   | time-to-positive windows, terminal-mismatch rules, B2 dinucleotide, the ace1r contradiction/rule/metric extracts                                                                   | `test-lamp-*`                                                                                   |
| Readout physics                             | 21   | HRM partition counts, ddPCR LoB/LoD, Sanger VAF floor, JAK2 CV reproduction, HRM minimum input                                                                                     | `test-readout-physics.R`                                                                        |
| Tetra-primer ARMS                           | 17   | wheat 2023 funnel (39 − 7 = 32), sorghum 2021 mismatch positions verified from sequence, blocked barley record, cycle-decay negative result, Liu 2012's 1686 published ASP primers | `test-arms2-*.R`, `test-arms-mismatch-position.R`                                               |
| KASP gates, platforms and priors            | 25   | wheat 2019, Sharplab 2021, maize 2020 design→work outcomes; Kwok 1990 mismatch severity; PolyMarker two axes; LefeVER's three rules and the 4-cell mismatch panel                  | `test-kasp-platforms.R`, `test-z7-priors.R`, `test-y32-homeolog.R`, `test-kasp-mismatch-type.R` |
| Design → working assay (published outcomes) | 15   | Makhoul 2020 strategies, Grewal 2020/2022 and Ramirez-Gonzalez 2014/2015 conversion with denominators, De Silva 2009 HRM amplicons/SNPs/strategies, CAPS/LAMP/dPCR funnels         | `test-z8-outcomes.R`                                                                            |
| HRM                                         | 19   | jsbbs2018 6 × 16 context grid and its summary/recomputed split, IJMS 2022 dTm controls and partition arithmetic, Liew 2004 bands                                                   | `test-hrm-curve.R`, `test-hrm-jsbbs2018.R`                                                      |
| dPCR                                        | 8    | Poisson residuals (Dong), E. coli MPS dilution series                                                                                                                              | `test-dpcr-*.R`                                                                                 |
| TaqMan / probes                             | 13   | wheat 2019 concordance and composition, rice 2018/2024 rates, SARS 2022, C. elegans 2018                                                                                           | `test-taqman-*`, `test-probes-*`                                                                |
| SNaPshot / typing                           | 10   | 14-locus × 5-input grid, stutter model                                                                                                                                             | `test-typing-*.R`                                                                               |
| RPA / AS-RPA                                | 10   | mismatch-position and terminal-mismatch claims                                                                                                                                     | `test-rpa-sources.R`                                                                            |
| CRISPR-Dx                                   | 4    | Cas12a single and adjacent-double mismatch ratios                                                                                                                                  | `test-crispr-communbiol.R`                                                                      |
| Sanger                                      | 4    | measured VAF floor for a named enrichment-plus-caller pipeline                                                                                                                     | `test-sanger-bda-braf.R`                                                                        |
| Polyploid dosage                            | 4    | cotton negative control, durum, sugarcane                                                                                                                                          | `test-polyploid-*.R`                                                                            |
| Tm oracles                                  | 2    | 4096 measured duplex Tm; HRM amplicon Tm                                                                                                                                           | `test-tm-mm2020.R`, `test-tm-hm2025.R`                                                          |
| CAPS                                        | 4    | 1236 palindromic enzymes with their cut offsets, 26 soybean primer pairs, SNP2CAPS predicted-vs-observed patterns, Liu 2016's within-family SNP composition                        | `test-caps.R`                                                                                   |
| Registry and IO                             | 3    | the 134-row source registry, the 45-entry backlog's three statuses, the example locus table's transition chemistry                                                                 | `test-scientific-data.R`, `test-lamp-status.R`, `test-io.R`                                     |

The table sums to 196 `test` rows, the state when it was last recomputed: each
row is assigned to the first pattern in the last column that matches its
`test_file`, so no row is counted twice. The ledger held 208 `test` rows on
2026-11-23, so 12 rows postdate this snapshot and the per-area split is not
recomputed here. `test-science-ledger.R` and `test-scientific-data.R` themselves
are the data-integrity gates rather than marker-family tests, and are counted
under _Registry and IO_.

## Claim → test

Every claim below is stated somewhere in this package's own text (README, this
document, a vignette or a roster). The `Held by` column names the test that
actually asserts it, checked by reading the test body rather than by inferring
from its filename. A dash in that column means no test exists; the reason is in
the same cell. `Kind` is one of:

- `code` — a machine check that fails if the claim stops holding;
- `measured` — a number this repository measured on this host, with the script
  or test that measured it (the test holds the band, the script prints the
  value);
- `literature` — the number comes from an external source in `inst/extdata/`
  and is held by the ledger tests;
- `not-checkable` — no machine check can settle it;
- `gap` — a machine check would need data the repository does not have.

| Claim                                                                                                                                                                                                                                           | Held by                                                                                                                                                                                                                                          | Kind          |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------- |
| `primer_tm()` accepts all four parameter sets and its default is the pinned table                                                                                                                                                               | `tests/testthat/test-tm-oracle.R`                                                                                                                                                                                                                | code          |
| SantaLucia and Allawi agree within the pinned tolerance, including at GC edges and homopolymers; Sugimoto sits above SantaLucia; Breslauer's gap appears only where it exists                                                                   | `tests/testthat/test-tm-oracle.R`                                                                                                                                                                                                                | code          |
| `oligotm` agrees with both tables at matched salt (max \|ΔTm\| 4.8e-7) and with the default table within 1.0 °C                                                                                                                                 | `tests/testthat/test-tm-oracle.R` (skippable without `oligotm`)                                                                                                                                                                                  | code          |
| Tm rises with GC, rises with length and is strand-invariant                                                                                                                                                                                     | `tests/testthat/test-tm-oracle.R`                                                                                                                                                                                                                | code          |
| Four parameter sets and three engines are reported as a matrix, never averaged                                                                                                                                                                  | `tests/testthat/test-tm-oracle.R` (`oligotm` column skippable without the binary)                                                                                                                                                                | code          |
| The model ranks the 64 canonical mm2020 duplexes (r > 0.9, Spearman > 0.85) and separates perfect from triple mismatch                                                                                                                          | `tests/testthat/test-tm-mm2020.R`                                                                                                                                                                                                                | code          |
| The Mm2020 model reads about 16.9 °C below the measured Tm (sd 0.64 < 2); three steps, −16.9 → −9.6 → ≈ −5.4 °C, are recorded, not corrected away                                                                                               | `tests/testthat/test-tm-mm2020.R`                                                                                                                                                                                                                | measured      |
| Every mm2020 row carries a `calibration_stance`, and calibration is forbidden                                                                                                                                                                   | `tests/testthat/test-tm-mm2020.R`                                                                                                                                                                                                                | code          |
| Primer3 2.6.1: Tm r 1.000, ρ 0.999, mean +3.92 °C, sd 0.25, max 4.66; hairpin dG −0.84 kcal/mol (ρ 0.67); heterodimer dG −1.75 kcal/mol (n = 4, ρ 1.000)                                                                                        | `dev/primer3-crosscheck.R` prints it (`results/primer3-crosscheck.log`); `tests/testthat/test-primer3-crosscheck.R` holds the bands                                                                                                              | measured      |
| Primer3 design solves 17/20, 8/20 and 5/8 across the three loci, median +1 bp, 0.06 s against 9.5 s                                                                                                                                             | `dev/primer3-crosscheck.R` prints it; `tests/testthat/test-primer3-crosscheck.R` (“design_kasp_markers and Primer3 both solve the same three loci”) holds the shape on three loci (full suite plus binaries)                                     | measured      |
| `specificity_bowtie()` is a two-strand `countPattern` equal to a manual count, and an absent primer scores zero                                                                                                                                 | `tests/testthat/test-bowtie.R` (skippable without `Rbowtie`/full suite)                                                                                                                                                                          | code          |
| Ground truth: 3′-seed 8/8, exact PCR 3/3, BLAST recall 1.000                                                                                                                                                                                    | `dev/specificity-oracle.R` prints it (`results/specificity-oracle.csv`); `tests/testthat/test-specificity.R` holds the equalities on fixtures                                                                                                    | measured      |
| A missing genome or engine gives `unknown`, never a silent pass; the BLAST path needs a database                                                                                                                                                | `tests/testthat/test-specificity.R` (“a missing genome or engine gives unknown, never a silent pass”) and `tests/testthat/test-scan.R` (“specificity_of routes an indexed genome to BLAST, and needs a db”; skippable without BLAST+/full suite) | code          |
| A genomic `N` never satisfies a primer, and a mismatch radius does not buy a hit across it                                                                                                                                                      | `tests/testthat/test-specificity.R`                                                                                                                                                                                                              | code          |
| `in_silico_pcr()` never pairs primers across chromosomes                                                                                                                                                                                        | `tests/testthat/test-insilico-pcr.R`                                                                                                                                                                                                             | code          |
| UNAFold's ~+1000 sentinel maps to dG 0, and poly-A cannot fold                                                                                                                                                                                  | `tests/testthat/test-oligo-dg.R` (skippable without OligoArrayAux)                                                                                                                                                                               | code          |
| LAMP and RPA report the additives no engine models: `scheme_conditions()` carries `unmodelled` = `"betaine"` for LAMP and `"crowding"` for RPA/CRISPR, so the gap is named rather than silently dropped                                         | `tests/testthat/test-oligo-dg.R` ("scheme_conditions() names the additives no engine models")                                                                                                                                                    | code          |
| `rprimer` designs independently and reports a match table                                                                                                                                                                                       | `tests/testthat/test-rprimer.R` (“rprimer_assays designs and reports a match table”; full suite only)                                                                                                                                            | code          |
| PolyMarker's two axes stay separate, so one 3′ column cannot pass                                                                                                                                                                               | `tests/testthat/test-y32-homeolog.R`                                                                                                                                                                                                             | code          |
| Conversion rows keep their denominator and named failures, stay ordinal and uncalibrated, and are never pooled                                                                                                                                  | `tests/testthat/test-z8-outcomes.R`                                                                                                                                                                                                              | code          |
| The mm2020 table reproduces the paper's totals (4096 duplexes; 64/576/1728/1728) and the canonical minimum 60.7 °C                                                                                                                              | `tests/testthat/test-tm-mm2020.R`                                                                                                                                                                                                                | literature    |
| Every criterion is fixture-backed or exempt with a reason, and a cited fixture names a real test                                                                                                                                                | `tests/testthat/test-falsifiable.R`                                                                                                                                                                                                              | code          |
| The ledger joins source → extract → a test file that exists and reads it; `use` is a closed vocabulary                                                                                                                                          | `tests/testthat/test-science-ledger.R`                                                                                                                                                                                                           | code          |
| An extract's DOI belongs to the row's source                                                                                                                                                                                                    | `tests/testthat/test-science-ledger.R`                                                                                                                                                                                                           | code          |
| No measured source lacks an entry, every `extracted` backlog row is covered, and the three review axes plus the PolyMarker rule are covered                                                                                                     | `tests/testthat/test-science-ledger.R`                                                                                                                                                                                                           | code          |
| Conversion rows reconstruct the named failure totals (Grewal 2022 sums to its own total)                                                                                                                                                        | `tests/testthat/test-scientific-data.R` (“conversion_rates: the Grewal 2022 named reasons sum to its own total”)                                                                                                                                 | literature    |
| A trace or artefact carries every field, the manifest names the trace version, an unknown schema aborts, and no absolute path or hostname is written                                                                                            | `tests/testthat/test-trace-schema.R`                                                                                                                                                                                                             | code          |
| A criterion has a source only where one states the bar, its tier is the printed target plus policy, and the threshold matches the rendered limit                                                                                                | `tests/testthat/test-criteria.R`                                                                                                                                                                                                                 | code          |
| Preconditions fail at their own boundary                                                                                                                                                                                                        | `tests/testthat/test-contracts.R`                                                                                                                                                                                                                | code          |
| The engine stamp names every engine with a non-empty version, and a missing engine stays `NA` rather than being dropped                                                                                                                         | `tests/testthat/test-engine-stamp.R`                                                                                                                                                                                                             | code          |
| A source row names its citation, DOI, access depth and what it anchors; no row claims more strength than its source type, and no conversion row yields a probability or is `calibrated`                                                         | `tests/testthat/test-scientific-data.R`                                                                                                                                                                                                          | code          |
| No bench validation is claimed                                                                                                                                                                                                                  | — no wet-lab data exists in this repository; the statement is true by construction                                                                                                                                                               | not-checkable |
| No fitted success predictor is claimed                                                                                                                                                                                                          | — maize 2020 carries 100/100 `design_success` and no `assay_outcome`, so there are no pairs to fit                                                                                                                                               | gap           |
| A population- or subgenome-specificity threshold works on A/B/D homeologs                                                                                                                                                                       | — roadmap `Y3`: the ledger has no primer × subgenome → outcome pairs, and footprint/dropout plus linkage remain `unknown`                                                                                                                        | gap           |
| The plot's Tm status uses a declared per-scheme window (a static proxy for the realised-primer-relative probe gate): dPCR 44-50 (probe floor 58), TaqMan 58-63 / 56-72, LAMP 55-65 / 55-62, SNaPshot 50-60 / 58-63, castPCR 55-65 (probe 63-80) | `tests/testthat/test-plot-methods.R` (“Tm guidelines are method-specific, not KASP gates”)                                                                                                                                                       | code          |
| A 16 Gb wheat genome can be scanned                                                                                                                                                                                                             | — roadmap `Y3`: not measured; the largest genome exercised here is far smaller                                                                                                                                                                   | gap           |
| Every commons-licensed file is CI-checked                                                                                                                                                                                                       | — roadmap `ZY18`: no such CI step exists; this is a defect, not a promise kept                                                                                                                                                                   | gap           |

Two limits on what the table says. Engine agreement is not evidence about a real
oligo duplex: every thermodynamic and specificity row above compares
implementations with each other or with a fixture, never with a measured
hybridisation. And a `measured` row states what this host produced on
2026-11-19; the test holds a band wide enough to survive a tool upgrade, so a
pass means the number moved no further than the band, not that it was reproduced
to the last digit.

## Validation benchmark

`dev/validation-benchmark.R` runs a fixed set of 25 tasks whose verdicts were
established _before_ the package ran them. Each row compares that `expected`
verdict with what the package reports (`observed`) and writes
`task | family | expected | observed | verdict | where | note` to
`.tmp/validation-benchmark.csv` (git-ignored). `Rscript dev/validation-benchmark.R
--check` runs the set twice and errors if the task count or the verdict
distribution moves.

Verdict vocabulary: `pass`, `fail`, `not_run` (with its reason), `gap` (no
independently established verdict exists). A task is never invented, and a
missing verdict is never reported as a pass.

As of 2026-11-24 the run is **25 tasks: pass 12, fail 13, not_run 0, gap 0**.
The thirteen `fail` rows are the soybean panel, and they are intended: the
operator recorded all 13 assays as `worked`, while every one carries at least
one failed gate in `soybean_worked_panel_kasp_gates.csv`. The benchmark records
that the package's current rules reject assays the field accepted; it does not
mask that as a pass.

| Task                                                            | Family   | Expected                                                                 | Observed                                                      | Verdict | Where it diverged                                                                                                                                                              |
| --------------------------------------------------------------- | -------- | ------------------------------------------------------------------------ | ------------------------------------------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `fixture:Chr1:300 A>C`                                          | fixture  | designed                                                                 | designed                                                      | `pass`  | snp:doomed — one flanking pair must survive on a clean template                                                                                                                |
| `fixture:Chr2:150 A>C`                                          | fixture  | rejected with a cause from every chemistry                               | rejected with a cause                                         | `pass`  | snp:doomed; hrm:NA; caps:NA — no chemistry may design on the homopolymer block                                                                                                 |
| `soybean:KASP_AGP10_branch_chr6_13.303311Mb`                    | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | AS GC x1; COM Tm >= AS Tm x1; longest homopolymer (bp) x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                                                      |
| `soybean:KASP_AGP10_height_chr18_58.02Mb`                       | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | amplicon length (bp) x1; COM Tm >= AS Tm x1; longest homopolymer (bp) x1 — 9 unknown gate rows; 28 gate rows; outcome recorded as worked                                       |
| `soybean:KASP_AGP9_habit_chr19_44.43Mb`                         | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | COM Tm >= AS Tm x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                                                                                             |
| `soybean:KASP_AGP9_height_chr19_46.10Mb`                        | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | AS GC x1; discrimination x1; longest homopolymer (bp) x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                                                       |
| `soybean:KASP_AGPO2_height_chr19_49.74Mb`                       | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | COM Tm >= AS Tm x1; discrimination x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                                                                          |
| `soybean:KASP_AGPO2_height2_chr19_51.07Mb`                      | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | AS GC x1; COM Tm >= AS Tm x1 — 9 unknown gate rows; 28 gate rows; outcome recorded as worked                                                                                   |
| `soybean:KASP_AGPO2_nodes_chr19_49.65Mb`                        | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | amplicon length (bp) x1; COM Tm >= AS Tm x1; discrimination x1; SSR (2-4 nt, >=3 copies) x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                    |
| `soybean:KASP_AGPO2_seed_chr19_48.03Mb`                         | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | COM Tm >= AS Tm x1; discrimination x1; inverted repeat (arms>=4, loop<=6) x1; SSR (2-4 nt, >=3 copies) x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked      |
| `soybean:KASP_AGPO8_protein_chr15_10.93Mb`                      | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | amplicon length (bp) x1; discrimination x1; inverted repeat (arms>=4, loop<=6) x1; longest homopolymer (bp) x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked |
| `soybean:KASP_DT2_flank_high_chr18_56.61Mb`                     | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | COM Tm >= AS Tm x1; SSR (2-4 nt, >=3 copies) x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                                                                |
| `soybean:KASP_DT2_flank_low_chr18_55.26Mb`                      | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | AS GC x1; COM Tm >= AS Tm x1; inverted repeat (arms>=4, loop<=6) x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                                            |
| `soybean:KASP_DT2_near_chr18_56.14Mb`                           | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | amplicon length (bp) x1; COM Tm >= AS Tm x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                                                                    |
| `soybean:KASP_JOINT_protein_chr15_11.28Mb`                      | soybean  | worked                                                                   | gates_reject                                                  | `fail`  | COM Tm >= AS Tm x1 — 8 unknown gate rows; 28 gate rows; outcome recorded as worked                                                                                             |
| `measured:mm2020 ranking`                                       | measured | r > 0.9, spearman > 0.85, -25 < mean(tm - measured) < -10, sd < 2        | n 64, r 0.958, rho 0.941, mean -16.86, sd 0.64                | `pass`  | test-tm-mm2020.R:88-113 — 64 canonical duplexes; the level bias is the salt term, not corrected                                                                                |
| `measured:primer3 Tm offset`                                    | measured | n 34, mean +3.92, sd 0.25, max 4.66 (primer3 - pkg)                      | n 34, mean +3.92, sd 0.25, max 4.66                           | `pass`  | results/primer3-thermo.csv — cross-tool offset recorded from a local primer3 run; additive, not modelled                                                                       |
| `measured:primer3 design (one locus)`                           | measured | 8/20 variants have a primer3 pair passing every gate                     | 8/20                                                          | `pass`  | results/primer3-design.csv — only one of the three loci in the VALIDATION.md claim is recorded here; the other two need a live primer3_core and are not invented               |
| `measured:specificity oracle`                                   | measured | 3'-seed 8/8, exact PCR 3/3, BLAST recall 1.000                           | seed 8/8, pcr 3/3, blast recall 1.000                         | `pass`  | results/specificity-oracle.csv — BLAST rows carry equal = NA by design; recall is compared instead                                                                             |
| `contract:crispr_spacer_hits counts raw occurrences`            | contract | 1 occurrence -> 1; 2 occurrences -> 2                                    | once 1, twice 2                                               | `pass`  | R/marker_crispr.R:268 — counts forward and reverse-complement occurrences; the design subtracts 1                                                                              |
| `contract:off_target == hits - 1`                               | contract | design_crispr_marker()$off_target equals crispr_spacer_hits(guide) - 1   | off_target 0, hits 1                                          | `pass`  | R/marker_crispr.R:471 — the on-target occurrence is not an off-target                                                                                                          |
| `contract:crispr flow uses a preamp pair, no RPA gate`          | contract | no gate_rpa_structure in scheme_def_crispr(); preamp defaults to "rpa"   | rpa gate present: FALSE; preamp default: "rpa"                | `pass`  | R/marker_crispr.R:600 — the vignette prose already matches this                                                                                                                |
| `contract:design_duplex warns on undeclared homology`           | contract | warning 'undeclared', gate_duplex_targets not_applicable, verdict fail   | warned TRUE, gate_duplex_targets not_applicable, verdict fail | `pass`  | R/marker_duplex.R:400-405 — an undeclared topology cannot pass, so the pairing fails                                                                                           |
| `contract:decision graph stamps gate status from the selection` | contract | one gate: no record ok, relaxed record relaxed, violated record violated | as:gate:gate_discrimination -> ok/relaxed/violated            | `pass`  | tests/testthat/test-falsifiable.R:405-440 — the graph reads the selection records; it does not restate the rules                                                               |
| `contract:select_with_ladder refuses an NA margin`              | contract | NA COM Tm -> best NULL; 60 -> best non-NULL                              | NA -> NULL; 60 -> non-NULL                                    | `pass`  | R/marker_select.R:458-465 — an unmeasured margin is never a pass                                                                                                               |

## Single-score baseline

The perspectives draft (§5) argues that any single scalar "design quality" score
is falsifiable: it does not separate accepted from rejected constructions. This
section makes that claim runnable. `dev/single-score-baseline.R` builds one
plausible scalar from checks the package already has and scores a panel whose
outcome was established before the package saw it.

_Panel:_ `inst/extdata/sharplab2021_assays.csv` — 50 KASP assays from Sharplab
2021 (doi:10.1186/s12870-021-02932-8), of which the paper reports 46 as
polymorphic between 1212 and B73 and 4 as not. Both classes carry full primer
sequences, so the scalar is computable for every row. Six other candidate panels
were rejected for having no rejected class or no sequences; they are listed in
the script header.

_Scalar_ (cost, lower is better), summed over the three primers:

| Term      | Feature                             | Weight           |
| --------- | ----------------------------------- | ---------------- |
| `tm_dev`  | `abs(primer_tm(seq) - 60)`          | 1.0 per °C       |
| `gc_dev`  | `abs(gc_frac(seq) - 0.5)`           | 20 per fraction  |
| `run`     | `max(0, longest_run(seq) - 4)`      | 2 per base       |
| `hairpin` | `max(0, -internal_hairpin_dg(seq))` | 1.0 per kcal/mol |
| `clamp`   | 3' base is not G/C                  | 1.0              |

_Metric:_ every pair (rejected, accepted) ordered the wrong way is an
inversion; `pairs = n_rejected × n_accepted`.

As of 2026-11-24 the run is **58 rows: pass 5, fail 52, not_run 1, gap 0**. The
falsification succeeds:

| Row                                   | Observed                                                    | Verdict   |
| ------------------------------------- | ----------------------------------------------------------- | --------- |
| `baseline:pairwise inversions`        | 71 of 184 pairs; 41 of 46 accepted below a rejected one     | `fail`    |
| `baseline:rejected inside the best k` | 0 of 4 rejected in the best 4 by cost (ranks 6, 33, 41, 43) | `fail`    |
| `baseline:sensitivity tm_only`        | 47 of 184 pairs; 27 of 46 accepted below a rejected one     | `fail`    |
| `baseline:sensitivity tm_dominant`    | 50 of 184 pairs; 29 of 46 accepted below a rejected one     | `fail`    |
| `baseline:sensitivity gc_dominant`    | 106 of 184 pairs; 46 of 46 accepted below a rejected one    | `fail`    |
| `baseline:sensitivity gentle`         | 66 of 184 pairs; 38 of 46 accepted below a rejected one     | `fail`    |
| `baseline:sensitivity tails_stripped` | 147 of 184 pairs; 40 of 46 accepted below a rejected one    | `fail`    |
| `baseline:corpus`                     | —                                                           | `not_run` |

The four alternative weightings were fixed before the run and all outcomes are
reported: no weighting separates the panel, and none was selected afterwards.
The tail-stripped row is a caveat rather than a control. The FAM/HEX tails are
not an additive offset — they enter the normalised Tm and GC of the whole
oligo — so stripping them moves the score from 71 to 147 wrong pairs. Both
still fail, which is the claim, but the size of the swing says a single scalar
is sensitive to a decision the sequence alone cannot settle.

The corpus row is `not_run` with its reason, not `0`: a corpus-wide run needs
recorded outcomes from `marker_record_outcome()`, and the tree holds none. Only
published panels have a known outcome, so the falsification runs on fixtures
and the corpus run waits for outcomes to accumulate.

_Consequence for the multi-axis verdict:_ the scalar did not survive, so the
multi-axis verdict is not shown to be over-complicated by this run; the
falsification is ready to cite in the perspectives draft.

## Success model

The single-score section above asks whether one number predicts anything. This
section asks the same question of the checks themselves, which is the direction
`dev/success-model.R` opens: a recorded outcome is a label, a check from
`checks.csv` is a predictor, and the corpus from `marker_outcome_corpus()` is
the join.

The join is by artefact and target. An outcome names the artefact it came from
and, when the design had one, the target it was built for; the script reads that
artefact's `checks.csv` and, per criterion, scores the row `1` when every
recorded instance passed, `0` when any instance failed, and leaves it `NA` when
the criterion was not evaluated. A criterion that never appears in the corpus is
not a predictor there.

Discrimination is reported as leave-one-out AUC per check, with ties counted as
one half. The threshold is twelve labelled rows with at least four in each
class; below it a predictor reads `not_run` with the reason spelled out, not a
number. The `success-model:sum` row is the complete-case sum over checks, so a
row is only compared with another row on the checks both were evaluated on.

```sh
Rscript dev/success-model.R            # writes .tmp/success-model.csv
Rscript dev/success-model.R --check    # recomputes twice, writes nothing
```

The table below is generated; `--check` reproduces it byte for byte and fails if
the file has drifted.

<!-- success-model:begin -->

_The corpus holds no recorded outcome, so no predictor can be scored. This run reports `not_run` with its reason rather than a number._

| Predictor              | n   | works | fails | AUC | LOO mean | LOO min-max | Verdict   | Reason                                                                      |
| ---------------------- | --- | ----- | ----- | --- | -------- | ----------- | --------- | --------------------------------------------------------------------------- |
| `success-model:corpus` | 0   | 0     | 0     | —   | —        | —           | `not_run` | no outcomes.csv under .tmp/corpus; record some with marker_record_outcome() |

The harness is `dev/success-model.R`; its header names the join rule, the tie convention and the thresholds above.
<!-- success-model:end -->

## Formal verification of the ranking

`ZY13`'s invariant is the lexicographic ordering behind `order_by_margins()`
(`R/marker_select.R:405`). `dev/formal/lexicographic_order.py` restates it as a
comparator over rank keys `d * v` (`d = +1` for a `desc` column, `-1` for an
`asc` one), and `dev/formal/README.md` argues the R-to-model mapping. Six
properties are confirmed by CrossHair (`crosshair-tool 0.0.114`, Python 3.13.5):
tie reflexivity, asymmetry, transitivity, totality, agreement with "the first
differing coordinate decides", and invariance when a value and its direction are
mirrored together.

The proof covers a bounded domain — key tuples and direction lists of at most
three entries — and it covers the model, not the R code; the README states both
limits and carries the command lines. Mutating the comparator's `>` to `>=`
produces a counterexample, so the checks are not vacuous. `dev/check.R` step 3d
runs only the model's dependency-free selftest, because the six CrossHair runs
need a virtualenv and about five minutes.

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
Rscript dev/as-tm-crosscheck.R     # Biopython Tm_NN -> results/as-tm-crosscheck.csv
Rscript dev/specificity-oracle.R   # TAIR10 ground truth -> results/specificity-oracle.csv
```

These need binaries that are not package dependencies (`oligotm`, `ntthal`,
`primer3_core`, `blastn`, `makeblastdb`), the OligoArrayAux build and, for the
allele-specific pair, a `python3` with Biopython; the matching tests skip with a
named reason when one is absent.

## Known gaps

- **No bench (wet-lab) validation, and none claimed.** The gap is wet-lab, not
  computational: `README.md` §Validation says in as many words that the short
  answer is published measurements, not a bench run by the maintainers, and that
  no wet-lab data was produced here. What would close it is a single operator
  experiment — the eight-row mismatch panel already used for the allele-specific
  cross-check, run through a real assay, scored against the package's predicted
  sign and within-primer ordering — and that is an operator action, not agent
  work. Separately, the wall-clock and RSS figures quoted in `ROADMAP.md` are
  host measurements from a single machine, not a benchmark suite:
  `capabilities("profmem")` is `FALSE`, so `bench::mark()$mem_alloc` is `NA`, and
  `dev/proc-mem.R` reads `VmHWM` from `/proc` instead.
- The allele-specific layer has one external comparator — Biopython's
  terminal-mismatch nearest-neighbour Tm — and it is a proxy: it agrees with the
  package's efficiency gap on the direction of all eight panel rows but ranks
  them differently (ρ = 0.238), so only the sign and the within-primer ordering
  are asserted. No external engine scores allele-specific amplification
  efficiency itself, and no wet-lab discrimination has been measured.
- No predictor of published per-assay success is fitted: the corpus has
  denominators at study level, but the tables that carry primer triples do not
  carry both outcomes (maize 2020 has 100/100 `design_success = TRUE` and no
  `assay_outcome`).
- 20 ledger rows are `blocked`: 10 name no extract at all, 10 name an extract
  that carries no verdict of its own yet, and each says why the paper could not
  be landed.
- 6 sources are closed access; they exist so that the number they would have
  anchored is visibly absent rather than silently invented. Every route tried for
  the closed and not-found rows on 2026-11-25 failed at the publisher: IA Scholar
  served its bot challenge, CSHL, OUP and Canadian Science Publishing answered
  403, Nature returned only its landing page, Springer returned its paywall
  notice, `api.fatcat.wiki` and `web.archive.org` timed out, and keyless web
  search was challenged. Each row's `access_detail` records the Unpaywall and
  OpenAlex measurement behind its `access` level, so a later reader can retry
  without re-walking a route that was already tried. The workspace's own network
  allowlist was measured on 2026-11-25: `www.ebi.ac.uk` answered 200, while
  `api.core.ac.uk` and `ouci.dntb.gov.ua` did not connect (000) and
  `academic.oup.com` returned 403 — so the six PDFs are an operator action
  (institutional access or a supplied copy), not an untried route.
- `DECIPHER`'s amplification efficiency is a model of hybridisation plus Taq
  extension. It is a proxy for allelic discrimination, not a replacement for an
  experiment.
- The formal model bounds key and direction lengths at three entries and covers
  the comparator only: `bncc_rank()`'s Pareto step and the gate cascade are
  checked by examples and mutation tests, not proven.
