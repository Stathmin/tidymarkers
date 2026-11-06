# ZW4.6 — naming the source of each criterion bar

State as of 2026-10-29, agent A3. Recount and change are in `R/marker_criteria.R`
(`criterion_sources()`, applied once in `criterion_spec()` right after
`bind_rows()`); the test is the `ZW4.6` block in
`tests/testthat/test-criteria.R`.

## Before / after

Counted through `criterion_spec()`, never by `grep` (the string is recycled by
`criterion_tbl()`, so `grep -c '"package assumption"' R/marker_criteria.R`
counts constructors, not rows).

|                                                                   | rows                                       |
| ----------------------------------------------------------------- | ------------------------------------------ |
| before: `package assumption`                                      | 191                                        |
| before: `primary literature`                                      | 1 (`as_mismatch_type`)                     |
| after: `package assumption`                                       | 176                                        |
| after: `kwok1990`                                                 | 7                                          |
| after: `santalucia1998`                                           | 9                                          |
| after: `primary literature` (unchanged)                           | 1                                          |
| **rows moved off the default**                                    | **16** (15 new + `as_mismatch_type`, whose |
| row already said `primary literature` and now also names its key) |

`evidence` is untouched: every row is still `rule_of_thumb`. A named source
does not raise the rung by itself — ZW4's ladder is about how the number was
fixed, not about whether a paper can be cited next to it.

## Rows that got a source

`kwok1990` — Kwok S, Kellogg DE, McKinney N, Spasic D, Goda L, Levenson C,
Sninsky JJ. _Effects of primer-template mismatches on the polymerase chain
reaction: HIV-1 model studies._ NAR 1990;18(4):999–1005,
doi 10.1093/nar/18.4.999, PMID 2179874, PMC330356.

Abstract read at `https://pmc.ncbi.nlm.nih.gov/articles/PMC330356/` (PMC bodytext;
the Europe PMC `fullTextXML` endpoint returns HTTP 500 and the PDF shell is behind
a 403). Verbatim: _"A:G, G:A, and C:C mismatches reduced overall PCR product
yield about 100-fold, A:A mismatches about 20-fold … mismatches of T with either
G, C, or T had a minimal effect on PCR product yield."_ The row
`kwok_3prime_severity` in `inst/extdata/mismatch_discrimination.csv` already
carries that sentence.

- `as:gate:gate_wobble` — "no weak G:T / A:C 3' mismatch". Same basis as
  `as_mismatch_type`: a 3'-terminal T mismatch is the minimal-effect class.
- `as:gate:gate_discrimination`, `arms:gate:gate_arms_disc`,
  `castpcr:gate:gate_castpcr_as3`, `castpcr:gate:gate_castpcr_mismatch`,
  `castpcr:gate:gate_castpcr_blocker` — each requires the allele-specific 3'
  base to sit on a weakened mismatch; that ranking is Kwok's.
- `as_mismatch_type` — unchanged in value; the map now also owns its key.

`kwok` does **not** state the numbers these rows threshold on (`discrimination >=
0.5`, the 0.1 dropout tolerance): it ranks mismatch classes by yield. The rows
say "whose rule decides the 3' base", which is what `as_mismatch_type` has said
since Z0.4.

`santalucia1998` — SantaLucia J Jr. _A unified view of polymer, dumbbell, and
oligonucleotide DNA nearest-neighbor thermodynamics._ PNAS 1998;95(4):1460–1465,
doi 10.1073/pnas.95.4.1460.

Not read directly (no PMC copy was fetched this session). It is cited as
Primer3 cites it, and that text **was** read:

- Primer3 manual, `https://primer3.org/manual.html`: `PRIMER_TM_FORMULA`
  "A value of 1 (_RECOMMENDED_) directs Primer3 to use the table of
  thermodynamic values and the method for melting temperature calculation
  suggested in the paper [SantaLucia JR (1998) … Proc Natl Acad Sci 95:1460-65
  http://dx.doi.org/10.1073/pnas.95.4.1460]", with `PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT`
  = 1 by default governing hairpin and oligo-oligo interactions.
- Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG.
  _Primer3 — new capabilities and interfaces._ NAR 2012;40(15):e115,
  doi 10.1093/nar/gks596, PMID 22730293, PMC3424584 (full text read). Reference
  19 is SantaLucia 1998; the body states "All thermodynamic models now use
  up-to-date nearest-neighbor parameters and salt correction formulae … Primer3
  filters out primers likely to form hairpins" and formulates a primer-pair
  quality score less than the sum of the individual primer qualities.

- `hairpin_dg`, `cross_dimer_dg`, `homodimer_dg`, `pair:gate:gate_structure`,
  `indel_size:gate:gate_indel_size_structure`,
  `dpcr:gate:gate_dpcr_probe_structure`, `taqman:gate:gate_taqman_probe_structure`,
  `snapshot:gate:gate_snapshot_structure`, `rpa:gate:gate_rpa_structure` — these
  gate on a ΔG that only the NN parameter set can produce, and the package
  computes it from that set (`salt_method = "SantaLucia1998-2"`, see ROADMAP
  ZW3.2/ZW3.8). The _method_ is the source; the bars (−3 hairpin, −9 dimer,
  the pairing rule "sum of qualities < 4") remain package assumptions.

## Why every other row stays `package assumption`

Every candidate was opened before it was rejected. None of them states the
package's number.

- **Primer length / GC.** Primer3 manual defaults are `PRIMER_MIN_SIZE 18`,
  `PRIMER_MAX_SIZE 27`, `PRIMER_MIN_GC 20.0`, `PRIMER_MAX_GC 80.0`. The package
  uses 18–25 and 40–60 %. Our GC window is _inside_ Primer3's default and ours
  is narrower — a default that does not equal the bar is not its source.
  Untergasser 2012 states no hard number at all ("Primer3Plus … reasonable
  default values"). Koressaar & Remm 2007 (Bioinformatics 23:1289–1291,
  doi 10.1093/bioinformatics/btm091, PMID 17379693) has no PMC copy and its
  abstract is only about a Tm and salt-correction formula, not length/GC.
- **Homopolymer `< 4`.** Primer3's `PRIMER_MAX_POLY_X` default is **5**, not 4.
- **Tm windows** (`as_tm` 55–65, `com_tm` 58–63, per-scheme primer windows).
  Primer3's own defaults are 57–63; no read source prints our window. SantaLucia
  and Koressaar give the _calculation_, not a chosen window.
- **Amplicon windows** (50–100, 60–150, 120–1200 bp, …) — design choices, no
  printed source.
- **Discrimination 0.5 and footprint dropout 0.1** — thresholds on our own
  simulator's output; no external measurement defines them.
- **`offtarget_hits <= 1`** — deliberate: a marker on a duplicated locus has to
  be _seen_, so the bar is not "zero perfect hits" and no paper states "at most
  one". Primer3's `PRIMER_MAX_NS_ACCEPTED = 0` is the opposite rule, and citing
  it here would misdescribe our gate.
- **Repeat / inverted-repeat / TaqMan probe checks, band gaps, RPA and CRISPR
  geometry, Sanger read distance, SSR copy count, the family checks** — scheme
  constants in `R/marker_stages.R` / the scheme files; no sentence stating the
  number was located, so they stay unnamed rather than borrow a neighbour's
  citation.

A vendor-guideline default that happens to sit inside our bar (e.g. a
`PRIMER_MAX_POLY_X = 4` in some other tool) was **not** used: the bar predates
any such reading, and naming it would claim provenance the code does not have.
That is a literature gap, not a code one — the same blocker ROADMAP records
for `Z8.1`.

## Unknowns a reviewer should know

- `santalucia1998` is cited from Primer3's own attribution, not from a page of
  the paper read in this session. Its `access` is `source` in
  `data-raw/proposed/a3.sources.csv`, with the route spelled out in
  `access_detail`.
- The two keys (`kwok1990` exists in `inst/extdata/scientific_sources.csv`;
  `santalucia1998` does not yet) must be merged there by whoever closes ZW4.6,
  and `tests/testthat/test-criteria.R` carries a one-line exclusion for
  `santalucia1998` that must be removed with that merge.

---

## Second pass, 2026-10-29

Read **only** through the Merovingen read-only MCP lane (collection `4RSE6KC2`,
ADR-0037). Nothing below is quoted from memory.

### Method (and why the first pass was short)

The first pass read `zotero_get_item_fulltext`, which for a PDF item returns at
most the first 10 pages and stops silently mid-sentence. 23 of the 30
full-text items in this collection are PDFs, so the _Discussion_, the
tables after page 10 and every appendix were out of reach, and the sweep also
drew in 13 items that belong to a different literature review (LLM abstract
screening, RAG retrieval, research-software FAIR) — they are in the
collection but carry no marker-design sentence.

This pass read all pages with `zotero_read_pdf_pages`
(`zotero_get_pdf_outline` is unavailable here: the sidecar lacks PyMuPDF).
28 of the 30 items were re-read end to end and cached under
`/tmp/read/pdf/<KEY>.md`; `677TI8JD` was complete in ten pages, and
`2W9JCWDM` has only a 5 kB abstract.

### Accepted

`fret_asqpcr2021` — Kalendar R, Shustov A, Akhmetollayev I, Kairov U.
_Designing allele-specific qPCR (ASQ) primers for C/T and A/G SNPs._
Front Plant Sci 2022;12:747886, doi 10.3389/fpls.2021.747886 (item
`6D9UWLMW`; 16-page PDF, read in full).

| row                                              | statement                                                                                                        | where                       |
| ------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------- | --------------------------- |
| `amplicon` (50-100 bp)                           | "the common reverse primer … is positioned such that an **amplicon of 50-100 bp** (**Table 2A**)"                | p. 7, _The First Component_ |
| `asq:tier:tier1` ("toward the 63 C ASQ optimum") | "the two forward primers … are designed to have a **similar range of melting temperature (Tm) between 55-56 C**" | p. 7, _Optimization_        |

Table 2A prints the same window for the worked HvSAP16 assay. The paper's
Prerequisites paragraph names ASQ, KASP and Amplifluor together, which is the
reason the row is KASP's `amplicon` and not a family-specific copy.

The second sentence is a **window, not an optimum**: 55-56 C there is the ASP
Tm _in the absence of Mg2+_, computed for the untagged ASP at 100 nM primer
with the Table 2A footnote spelled out (50 mM K+, 3 mM Mg2+, without Mg for the
untagged ASP). It is narrower and lower than the package's `as_tm` 55-65 and
`com_tm` 58-63, so it supports the **direction and floor** of the ASQ tier, not
the 63 C number the tier's label prints. If the tier label is ever re-worded to
carry its real bar, this is its citation.

### Rejected, with the sentence that was actually there

| candidate row                                                                                   | what the full text states                                                                                                     | why not our bar                                                                                                                    |
| ----------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `max_run` (< 4)                                                                                 | "PSQ doesn't allow genotyping SNPs in **homopolymeric stretches**" (`I3EXS422`, p. 3, ref. 24, no length)                     | no number anywhere in the 29 texts                                                                                                 |
| `ssr` (unit 2-4 nt, >= 3 copies)                                                                | SSR markers are _used_, never _thresholded_ (`MGT8CSL3` 53 mentions, `J2D3QY8K` 16, `KYA69U56` 24, `I3EXS422` 54)             | design practice, not a stated minimum                                                                                              |
| `inverted_repeat`                                                                               | "homologous or homopolymeric regions due to the limited assay design space" (`JHLSJPAH`)                                      | no arm/length/loop number                                                                                                          |
| `family_crispr_pam` (TTTV)                                                                      | "the PAM sequence **(TTTV)** of Cas12a" (`TKHWSS46`, p. 8)                                                                    | the bar is right; the package already cites `TKHWSS46` on the assay rows, so this is an existing attribution, **not** a new source |
| `family_sanger_read_distance` (>= 50 bp)                                                        | amplicons 85 / 101 bp (`WGHRARTK`), 400 bp / < 400 bp (`PSS8SYZT`)                                                            | amplicon sizes, never the variant-to-primer distance                                                                               |
| `as_gc`, `com_gc` (0.40-0.60)                                                                   | "GC content 40-60%" (`MGT8CSL3`, p. 5, its own CAPS criteria)                                                                 | a different tool's criteria for a different marker; first-pass rule still stands                                                   |
| `family_caps_band_gap` (>= 20 bp)                                                               | "the resolution ratio of the 1% agarose gel electrophoresis (**100 bp**) … >100 bp away from each [other]" (`MGT8CSL3`, p. 5) | the paper's own 100 bp on a 1 % gel is _not_ our 20 bp on capillary                                                                |
| `family_probe_tm_delta` (>= 5 C)                                                                | "probes were designed to have different Tm (63 C and 54 C)" (`PSS8SYZT`, p. 4)                                                | separates two probes in one tube; says nothing about probe versus primer                                                           |
| `as_tm` / `com_tm`                                                                              | 55-56 C (`6D9UWLMW`, accepted above); 63 / 54 C and 75-80 / 90-95 C amplicon melts (`PSS8SYZT`)                               | no source prints 55-65 or 58-63                                                                                                    |
| every tier row that reads "Tm balance" (`caps:tier:tier2`, `arms:tier:tier3`, `hrm:tier:tier2`) | no sentence pairs a tolerance in C with a balance rule                                                                        | the 1 C / 5 C windows are package constants                                                                                        |

That is **2 rows**: `amplicon` and `asq:tier:tier1`.

**Both stay unwired, and the reason is structural, not laziness.**
`amplicon` has `scope = "family"` and is inherited by **all 20 schemes**
(`hrm`, `sanger`, `lamp`, `rpa`, …), where 50-100 bp is not the bar. Citing an
ASQ-specific window on that row would misattribute the other 19; a correct wire
needs a scheme-specific row first. `asq:tier:tier1`'s label prints "toward the
63 C ASQ optimum", while the source states 55-56 C — the citation only becomes
true after the label is reworded to its real floor.

### What remains unsourceable, and the class of reason

Reconciled 2026-10-29: `criterion_spec()` has **193 rows**, of which **177**
carry `source = "package assumption"` (the ROADMAP's 192/176 predates the §W
`family_single_base` row). Of those 177, the 2 above move to
`fret_asqpcr2021` when wired; **175** stay `package assumption`, in four classes:

1. **Proprietary or unstated reaction conditions** — the Tm windows
   (`as_tm`, `com_tm`, every tier that ranks on Tm) are the same bars Primer3
   ships with different numbers. No paper states _this_ window under _this_
   buffer.
2. **Our simulator's output thresholds** — `discrimination` 0.5, the three
   dropout tolerances 0.1, `offtarget_hits <= 1`. These threshold a number the
   package itself computes; no wet-lab paper can state them.
3. **Geometry that the schemes inherit, not choose** — PAM `TTTV`,
   `family_lamp_fip_alt`, `family_rpa_thf_off`, `family_probe_extension`
   "adjacent". Each is the definition of the chemistry it cites; a citation
   would describe the method, not the bar, and the method's own key
   (`cas12a_*`, `lamp_*`, `rpa_*`, `fret_asqpcr2021`) is already on the assay row.
4. **Below-the-radar scheme constants** — band gap 20 bp, Sanger 50 bp read
   distance, SSR 2-4 nt x >= 3 copies, RPA 30-35 nt, LAMP inner-arm balance
   5 C, size delta 3 bp. None of the 29 texts states one of these numbers, and
   the ones that come closest (MGT8CSL3's 100 bp gel, WGHRARTK's 85 bp
   amplicon) state a different quantity.

### Costs and caveats a reviewer should know

- Six items in the collection have no usable attachment in this lane:
  Ayyadevara 2000 (`5LLSU4HN`), Oliveira 2020 (`MMARQH5C`), Stadhouders 2010
  (`8CSIT9FS`), Liew 2004 (`6W4GLHNY`), plus a 5 kB abstract-only record
  (`2W9JCWDM`) and a DOI-only stub (`GVJ3E2MH`). Each is already cited,
  excluded or marked closed-access in `inst/extdata/scientific_sources.csv`,
  and none of them is load-bearing for a row here.
- `HJ77998N` (ZNA) prints "Tm raised to 67 C": a Zn2+ additive, no bearing
  on any row.
- `amplicon` is shared by every scheme, so a single string change reaches
  many `scope`s; the row is one row in `criterion_spec()` but several gates in
  the schemes. Whichever merge path is used, the source string is set once.

## Machine-readable reason classes (2026-11-03)

The four prose classes above now ship as a per-row table,
`inst/extdata/criterion_provenance.csv` (columns `id, source, reason_class,
reason`, one row per `criterion_spec()` id, sorted by id). It is generated by
`data-raw/proposed/gen-criterion-provenance.R` from the live spec, so ids and
sources cannot be typed out of date, and `tests/testthat/test-criteria.R`
(the ZW4.6 block) fails if the two drift.

Measured totals on 2026-11-03 (`criterion_spec()`, 212 rows): **171** rows
carry `source = "package assumption"` and **41** rows carry a named source
(`kwok1990` 7, `santalucia1998` 9, plus 25 rows citing a source inline as a
`doi: … PMCID: …` snapshot rule). The 193/177/175 counts in "What remains
unsourceable" above are stale and kept only as history.

Classes:

| `reason_class`        | rows | meaning                                                                                                                     |
| --------------------- | ---- | --------------------------------------------------------------------------------------------------------------------------- |
| `named_source`        | 41   | the row's `source` states the bar                                                                                           |
| `conditions_unstated` | 119  | the Tm/GC/composition windows and length bands: no cited source states this bar under these reaction conditions (class 1)   |
| `simulator_threshold` | 18   | discrimination, dropout, off-target counts: thresholds on a number the package simulates (class 2)                          |
| `inherited_geometry`  | 17   | PAM, LAMP FIP, RPA THF, probe adjacency: the geometry of the chemistry, whose own key is already on the assay row (class 3) |
| `scheme_constant`     | 17   | band gap, Sanger read distance, SSR, LAMP arm balance, size delta: scheme defaults no text prints (class 4)                 |

No threshold, bound or `evidence` rung changed to produce this table; it
records provenance only.

## Round 4, 2026-11-05 — same-quantity sourcing

The user authorised a second sourcing pass and an explicit acceptance rule:

> Raise a row when a source prints a design rule for the **same quantity**, even
> if the printed range is wider or offset from the package bar. Record the exact
> printed number and the offset in `access_detail`. **Never** change `bounds` /
> `limit`; the rung rises, the threshold does not. Every raised row must name a
> registry key whose `claim` is a quotation re-runnable from the cached full text.

Result (`criterion_spec()`, 212 rows): `rule_of_thumb` 146 → **109**,
`primary_literature` 54 → **91**, `vendor_guideline` 12, `fitted_here` 0;
`source = package assumption` 123 → **100**. 37 rows moved.

### New registry keys (all DOI verified via Crossref, quotes grepped from the full text)

| key | DOI | anchors |
| --- | --- | --- |
| `wang2025_multikasp` | 10.1186/s12870-025-07005-8 | `duplex:gate:gate_duplex_complement`, `duplex:gate:gate_duplex_selfcomp` (exact: "<6 bp inter-pair" / "<=5 self-complementary") |
| `shimizu2025_kasp` | 10.1270/jsbbs.24052 | (registered for the -10 kcal/mol self/cross-dimer rule) |
| `degenprime2024` | 10.1093/bioadv/vbae044 | `hairpin_dg` and every `*_structure` gate (exact -3 kcal/mol combined hairpin/dimer filter) |
| `zhao2026_snapshot` | 10.1186/s12864-026-13173-0 | `snapshot:gate:gate_snapshot_tm` (printed 53-60 C vs [50,60]) |
| `tso2022_snapshot` | 10.1186/s13007-022-00923-w | `snapshot:gate:gate_snapshot_placement` (exact: "one base upstream of the SNP") |
| `li2026_indel` | 10.1186/s12863-026-01486-8 | `indel_size:gate:gate_indel_size_unique`, `offtarget_hits` (non-specific primers removed) |
| `zheng2025_ssr` | 10.3390/plants14233589 | `ssr:gate:gate_ssr_offtarget` (multiple in-silico amplicons excluded), `ssr:gate:gate_ssr_flanks` (unique flanking sequence) |

### Rows raised onto existing keys

- `caps:gate:gate_caps_amplicon` and `caps:gate:gate_caps_site`, plus the family
  `caps_site_defined` / `caps_site_separates` — `liu2016_caps`: cut site in the
  middle of the amplicon, PCR products 200-1000 bp (bar [120,1200]).
- `cross_dimer_dg`, `homodimer_dg` — `vavilov2024`: ΔG >= -9 kcal/mol for primers
  and probes (exact).
- `taqman:gate:gate_taqman_centrality`, `taqman:gate:gate_taqman_probe_specific`,
  `family_probe_spans` — `vavilov2024`: polymorphic site in the middle third of
  an 18-30 nt probe; a probe must bind only one region.
- `rpa:gate:gate_rpa_amplicon`, `rpa:gate:gate_rpa_probe`,
  `rpa:gate:gate_rpa_primer_len`, `family_rpa_primer_length` — `lobato2018`:
  amplicons 100-200 bp (exact), nfo probe 46-52 nt with >=30/>=15 nt THF arms
  (exact), primers 30-35 nt (the review then softens it).
- `amplicon` — `fret_asqpcr2021`: "an amplicon of 50-100 bp" (exact).
- `snapshot:gate:gate_snapshot_comp` — `dekeyzer2026`: probe GC >30, homopolymer
  <5, no repeat motif >6 (offset).
- `crispr:gate:gate_crispr_amplicon` — `lobato2018`: RPA 100-200 bp under the
  80-250 bp bar.
- `sanger:gate:gate_sanger_primer` and the generic `ssr` — `faasa2026`: avoid
  homopolymers / hairpins / dimers; avoid repetitive sequences.
- `family_lamp_fip_alt`, `family_lamp_fip_variant` — `nemeth2025`: the SNP sits
  at the FIP 3' end (or BIP 5' end); variants at FIP/BIP ends inhibit LAMP.
- `family_lamp_inner_inside` — `nemeth2025`: FIP/BIP 3' ends lie in regions
  flanked by F3/B3.
- `family_size_repeat_unit` — `ssr_gel2026`: motifs typically 1-6 bp (bar 2-4 nt).

### Printed bars still not found (stayed `rule_of_thumb`)

- **All 54 `*:tier:*` rows** — the limit is `maximise (per-tier bottleneck min)`,
  an optimisation objective, not a threshold any paper prints. Hard floor.
- `arms:gate:gate_arms_outer` [200,800] — every tetra-primer ARMS paper read
  prints only concrete band sizes (191/494, 382/667, ...), no window.
- `lamp:gate:gate_lamp_outer_tm` [55,62], `lamp:gate:gate_lamp_cross_dimer`
  (run >=8), `family_lamp_inner_tm_balance` (<=5 C) — the LAMP reviews read say
  only "somewhat lower" / "slightly higher"; no number.
- `sanger:gate:gate_sanger_amplicon` [400,800], `snapshot:gate:gate_snapshot_amplicon`
  [100,300] — no printed amplicon window (Sanger protocols print read length,
  ~700 bp-1 kb, a different quantity).
- `taqman:gate:gate_taqman_discrimination` (dG >=1.5) — Vavilov gives a direction
  (bigger ΔTm = better) but no numeric threshold.
- `crispr:gate:gate_crispr_pam` / `_seed` — the combined row needs Cas12a + Cas13a
  + Cas13b tokens; no single source prints all three.
- `snapshot:gate:gate_snapshot_unique` / `sanger:gate:gate_sanger_unique` /
  `dpcr:gate:gate_dpcr_probe_specific` — only generic "check specificity", no
  printed site count.
- The rest are the package's own bookkeeping and simulator outputs:
  `site_*`, `*_dropout`, `discrimination`, `in_silico_pcr`, `as_n_*`,
  `com_n_var`, `lax:*`, `family_crispr_*`, `family_rpa_thf_off`,
  `family_probe_extension`, `family_single_base`, `family_sanger_pair_defined`,
  `castpcr_as3/mismatch/blocker` (Kwok ranks, does not measure), and the
  chemistry geometry rows.

### Why 70 was not reached

`rule_of_thumb` stands at 109. The floor under the same-quantity rule is 54 (the
tiers) plus the simulator/geometry rows above, which have no external printed
threshold at all. Reaching 70 would require either re-labelling simulator
thresholds as "literature" (false citation) or dropping/loosening the tier rows —
neither is a sourcing change. The gap is reported, not papered over.

**Zotero.** The user authorised importing the new DOIs into the library
(collection `tidymarkers — marker design sources`, tag `merovingen:tidymarkers`).
Writes are the operator's action per ADR-0037; this round only adds the records
the registry cites.
