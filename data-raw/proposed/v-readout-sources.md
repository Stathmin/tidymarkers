# ROADMAP §V — readout models: which "conditions → outcome" sources exist (2026-10-29)

Source-gathering only. No criterion changed, no code changed, no `inst/extdata/`
file written or edited. Prior state and verdicts: `data-raw/proposed/v-readout-conditions.md`
(read, not rewritten here) and `ROADMAP.md` §V.

Rule kept literal: a **set** is a table a model must reproduce, with every condition
the model reads printed in the same source. Prose that states a measured effect is
cited as evidence but marked **not a set**. "Absent" below means absent from the
executed search, not proof that it does not exist.

## Search performed

- **Zotero, read-only MCP lane** (`zotero_get_collection_items` on `4RSE6KC2`, then
  `zotero_get_item_fulltext` / `zotero_get_item_metadata` per item). Collection
  listed live: **72 items**. Fulltext returned for 29 of them (44 short/no
  attachment); the rest title/abstract only.
- Queries: `trans-cleav|collateral`, `migration distance|relative mobility`,
  `agarose`, `recombinase|RPA`, `Mg2|MgCl`, `hybridiz`, plus a keyword scan of all
  128 cached `/tmp/read/*.md`.
- `zotero_semantic_search` is **not installed** in this Merovingen venv (returns
  "Install it with pip install 'zotero-mcp-server[semantic]'"), so no vector search
  was possible; `zotero_search_by_tag` only returns the same 72 items.
- Open access outside the library via the Europe PMC REST API (search, then
  `fullTextXML`). Every extra citation below was read in full text, not from memory.

## Verdict per readout

| readout              | status now                                                      | sources read | measured conditions→outcome table?                                                                                                                                              |
| -------------------- | --------------------------------------------------------------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| HRM melt             | **partial** (improved: a salt/Mg axis now exists)               | 6            | **yes, two** — Oliveira 2020 Table 1 (ΔTm per mismatch, fixed salt); Huguet 2017 Table 1 (Δg and salt-correction `m` per NN motif vs [Mg²⁺]). Neither is on a **PCR amplicon**. |
| dPCR Poisson         | **sufficient** (unchanged)                                      | 3            | **yes** — already in tree (Dong 2015 Suppl. S5–S8; Tere-Peña 2024 Suppl. S1–S7). No third raw-count set found.                                                                  |
| CRISPR-Dx collateral | **absent** as a kinetics set, but two factor-level tables found | 3            | **yes, two** — Nalefski 2021 Table 1 + Method table; Huang 2022 Table 2 (+ KCl × signal). Both are **endpoint/rate constants**, not a time course.                              |
| virtual gel          | **absent** for bp→mm; one calibration method found              | 3            | **yes** for the _method_ (Ziraldo 2019 Eq. 2 log₁₀ M → y, with gel %, V cm⁻¹, duration); **no numeric table**.                                                                  |
| LAMP/RPA kinetics    | **partial** (Q1 still open)                                     | 4            | **no** — no source prints copies per reaction **and** time-to-positive **and** Mg/conc. in one table.                                                                           |

## HRM melt

Set in tree: `inst/extdata/readout_hrm_diatom2025.csv` (Zhou 2025,
`10.1371/journal.pone.0321885`, `PMC12002481`, Tables 6–9). Read, unchanged.

What was missing in the ROADMAP verdict was **not** a Tm table (Diatom has one) but
any measured salt axis with the salt concentration printed. Two now exist:

1. **Oliveira LM, Long AS, Brown T, Fox KR, Weber G.** _Chem Sci_ 2020;11:8273–8287.
   `10.1039/d0sc01700k`, `PMC8163305`. Read in full.
    - **Table 1** (p. 8279): mean melting temperature ± SD per central base pair
      over 4096 sequence contexts — GC 66.1 ± 1.9 °C, AT 64.1 ± 1.8 °C, GT 59.9 ± 2.2 °C,
      CC 54.7 ± 2.6 °C. ESI Table S1 has all 4096 measured Tm values.
    - Conditions (Materials and methods, "Experimental data used"): **50 mM NaCl,
      10 mM sodium phosphate pH 7.4, 1.0 µM total strand**, 9-bp duplex, LightCycler,
      Tm from the first-derivative maximum.
    - **Measured table: yes.** Salt is **fixed**, so it constrains a mismatch-ΔTm
      term but cannot fit a salt term.
2. **Huguet JM, Ribezzi-Crivellari M, Bizarro CV, Ritort F.** _Nucleic Acids Res_
   2017;45:12921–12931. `10.1093/nar/gkx1161`, `PMC5728412`. Read in full.
    - **Table 1**: measured Δg° and the salt-correction prefactor `m` in
      **Na⁺ and Mg²⁺** for all ten NN motifs (e.g. CC/GG `m_Mg` = 0.032 kcal/mol,
      AA/TT 0.086). **Table 2**: Δh°, Δs°, `m` per motif used for Tm prediction; the
      fit uses 92 oligos of 10–30 bp at **[Mg²⁺] = 0.5, 1.5, 3, 10 mM** (Owczarzy 2008
      data). Reported prediction error χ² = 3.37 (°C)² vs UO 9.1.
    - Conditions: 10 mM Tris pH 7.5, 0.01 % NaN₃, MgCl₂ varied; unzipping at 25 °C.
    - **Measured table: yes**, and it is a _salt → energy_ grid, which is exactly the
      axis HRM lacked. Caveat: single-molecule unzipping, not a PCR amplicon.

Still absent for HRM: an **amplicon-level** table with Mg²⁺/salt printed next to
measured Tm, i.e. something that ties the Diatom 28 amplicons to these parameters.
The in-tree Diatom rows keep `mg_mm = NA` for that reason.

Prose only (cited, not a set): **Kalendar et al.** _Front Plant Sci_ 2022,
`10.3389/fpls.2021.747886` — "magnesium concentration within optimal limits
(1.5–3 mM) makes no difference"; ASPs/UP Tm 55–56 °C; annealing 55 °C then 62 °C.
Library item `6D9UWLMW`.

## dPCR Poisson

Set in tree: `dpcr_dong_poisson_residual.csv` (Dong 2015, `10.1038/srep13174`,
`PMC4548226`, Suppl. Tables S5–S8) and `dpcr_mps_ecoli_dilution.csv` (Tere-Peña 2024,
`10.3390/mps7060094`, Suppl. S1–S7). Read, unchanged.

No third independent raw-partition-count set was found: papers found in the lane
print a _concentration_ and a partition-acceptance threshold rather than per-well
positive/negative counts. Verdict stays **sufficient**; nothing to add.

## CRISPR-Dx collateral cleavage

No set in tree (`communbiol2024_crrna_ratio.csv` is ratio-only). Two measured
factor-level tables were found and read:

1. **Nalefski EA, … Madan D.** _iScience_ 2021;24:102996. `10.1016/j.isci.2021.102996`,
   `PMC8411246`. Read in full.
    - **Table 1** (Results): Cas12a LbCas12a trans-cleavage kcat/KM per substrate —
      C10 2.77 s⁻¹ / 28 nM, C20 2.26 / 28, DNaseAlert 1.30 / 26, dsDNA 0.044 / 13,
      U20 0.118 / 133 nM; Cas13a LbuCas13a U5/U10/U20/RNaseAlert rows.
    - **Method table** (STAR Methods): per-figure **Cas (nM), crRNA (nM), target,
      trans-substrate (nM), time (h)** for every figure, e.g. 5D = 1.0 nM Cas,
      1.0 nM crRNA, 10 pM target, 2 h.
    - Buffer: 10 mM Tris-HCl pH 7.5, **10 mM MgCl₂**, 1 mM TCEP, 50 µg/mL BSA,
      0.01 % Igepal; 37 °C unless noted; 20 mM HEPES pH 6.8, 50 mM KCl, 5 mM MgCl₂
      for Cas13a.
    - **Measured table: yes** (enzyme constants + a condition grid). **Not** a
      time → fluorescence course, and **not** a target-concentration → rate series.
    - Note: this is the _iScience_ kinetics paper, not the NAR 2025 paper
      (`nalefski2025cas12a`) that the ledger rejected for having everything in its
      supplement. This one prints the numbers in the body.
2. **Huang M, … Wang T.** _Front Microbiol_ 2022;13:754995. `10.3389/fmicb.2022.754995`,
   `PMC8822122`. Read in full.
    - **Table 2**: primer pair → positive-judgment time (Primer1 16:00, Primer2 16:30,
      Primer3 16:30, Primer4 13:30, **Primer5 13:00**, Primer6 15:30, blank negative).
    - **KCl × signal** (Results, Fig. 2E): 1000 copies/µL tested at **0.5, 1, 2 µM KCl**,
      2 µM chosen — a printed factor-level outcome (values are plotted, not tabulated).
    - Reaction: RPA 50 µL, 2.5 µL each 10 µM primer, 39 °C 20 min; Cas12a 1 µM +
      crRNA 1 µM 25 °C 20 min; 20 min at 37 °C, read every 30 s; ssDNA FQ/FB probes.
    - **Measured table: yes** for primer→time. **No** concentration → endpoint value
      table (that is in Fig. 3C/3D, image panels).

Still absent: a printed **target concentration → fraction cleaved at time t** table
at a stated enzyme/reporter concentration and temperature. That is the set a
`simulated = TRUE` Cas12a model would have to reproduce.

Related, library item `TKHWSS46` = **Mao X, … Wang F.** _Commun Biol_ 2024;7:1454,
`10.1038/s42003-024-07173-7`, `PMC11541961` (repo `communbiol2024_crrna`). Read.
Ceiling values used elsewhere in the tree are confirmed here: fold-change in
fluorescence **within 1 h** (Fig. 1b, 2a legends), LOD **10² copies per reaction**
(Results, Supplementary Fig. 18b,c). **Prose + heatmaps, not a set.**

## Virtual gel

No set in tree; `hrm_liew2004_restated_bands.csv` is a citation restatement. The
library scan (128 full texts, 63 with "gel") returned **no** printed band-distance
calibration; `DZAUXPRX` (Badolo 2012, `10.1186/1475-2875-11-227`) says only
"2 % agarose gel electrophoresis at 100 V" (Fig. 2B/3D), qualitatively.

The sizing _method_ exists and is printed with all four conditions:

- **Ziraldo R, Shoura MJ, Fire AZ, Levene SD.** _Nucleic Acids Res_ 2019;47:e92.
  `10.1093/nar/gkz534`, `PMC6895257`. Read in full.
    - **Methods, "Library preparation and gel-electrophoresis-based size selection"**:
      **1 % agarose**, TBE (50 mM Tris-borate, 1 mM Na₂EDTA), **3.0 V cm⁻¹**, ethidium
      bromide 0.5 µg/mL; ladders HiLo, NEB 100-bp, λ-AciI (206 discrete sizes, 2–1086 bp).
    - **Fig. 2A/B vs 2C/D**: same field strength, **16 h** (long format) vs **2 h**
      (mini gel) — a printed duration contrast at fixed voltage.
    - **Eq. (2)**: piecewise-linear log₁₀(M) → migration distance y (mm), with
      M(l) = 607.4·l + 157.9; band width σ interpolated from the ladder.
    - **Measured table: no numeric bp→mm table.** The relation is a formula fitted to
      a ladder whose coordinates are only in the figures. Supplementary Table S1/S2
      hold _average fragment size_, not distances.
- Cited by that paper as the original calibration: **Southern EM.** _Anal Biochem_
  1979;100:319–323, `10.1016/0003-2697(79)90235-5`; **Elder J, Southern E.**
  _Anal Biochem_ 1983;128:227–231, `10.1016/0003-2697(83)90369-x`. **Not read**
  (not open access here); cited only as the method's origin, not as data.

A buildable virtual-gel set would be a ladder's (bp, mm) pairs at the printed
1 % / 3 V cm⁻¹ / 16 h. Neither the library nor the Europe PMC query
(`"migration distance" AND agarose AND "fragment length"`, 104 hits, top 10 checked)
printed one.

## LAMP / RPA isothermal kinetics

Set in tree: `inst/extdata/readout_lamp_erysiphe2026.csv` (Barone 2026,
`10.3390/jof12010079`, `PMC12842720`, Tables 2–3). Read, unchanged. Q1 stays open:
`copies_per_reaction` is `NA`.

Searched this round and found **no** new buildable set:

- **Ma Z, … Shen J.** _Front Microbiol_ 2026;17:1934648. `10.3389/fmicb.2026.1934648`,
  `PMC13597815` (the previous handoff's open lead). Read in full.
    - **Table 2** is a _recipe_, not an outcome table: 50 µL = 29.4 µL A buffer,
      2 µL of 2 µM F, 2 µL of 2 µM R, 5 µL template, 9.1 µL water, **2.5 µL of 280 mM
      magnesium acetate**; 37 °C, 30 min. **Table 5**: qPCR LOD (10⁴ 22.96 … 10⁰ 36.98,
      negative). LOD of RT-RPA-LFNA = **10¹ copies/µL**, T/C grey ratio threshold >0.1.
    - **Measured table: no time course.** Single Mg₂(OAc)₂ concentration, single time,
      one LOD endpoint. Not a set.
- In-tree LAMP tables were re-checked for a copy axis and none carries one:
  `lamp_snp_lamp_lfd2025_conditions.csv` (8 rows, all `evidence_tier = weak`,
  outcome_definition says "not a readout outcome"),
  `lamp_comp2018_*` (design-space and fixed-duration rows),
  `lamp_mismatch_rtlamp2019_tt.csv` (Tt at **3 × 10⁴ RNA copies** for every row —
  one copy level, so no copy axis),
  `ace1r_aslamp2015_time_to_positive.csv` (4 rows, unquantified template, 75-min window).

RPA: nothing new. The library item is **Natoli ME et al.** _Anal Chem_ 2021;93:4832–4840,
`10.1021/acs.analchem.0c04191` (`MR64AR9J`), but its fulltext is not retrievable
through the lane (644 bytes), so it is cited for its title and DOI only.

What a buildable LAMP set would need: one table with **copies (or genomes) per
reaction, time-to-positive, primer set and Mg²⁺/KCl**, with non-detections written as
"> run length". No source read this round has that shape.

## 2026-11-03 — second search pass (Europe PMC open access, read the tables)

Nothing promoted: every chemistry still returns `simulated = FALSE` except hydrolysis
probes. `R/`, `inst/extdata/`, `ROADMAP.md` and `NEWS.md` untouched.

| readout              | rows found                                       | verdict                                      | missing column                                                |
| -------------------- | ------------------------------------------------ | -------------------------------------------- | ------------------------------------------------------------- |
| HRM melt             | 28 (in tree) + 20 (Oliveira/Huguet salt rows)    | **partial**                                  | `mg_mm`, `salt_mm` on an amplicon-level Tm table              |
| dPCR Poisson         | 60 + 102 (in tree)                               | **sufficient** (unchanged, not re-litigated) | —                                                             |
| CRISPR-Dx collateral | 0                                                | **none**                                     | a printed `target_conc × time → fraction_cleaved` table       |
| virtual gel          | 0                                                | **none**                                     | a printed `bp → migration_mm` ladder at stated gel %, V, time |
| LAMP                 | **24** (`readout-lamp-zhang2026.conditions.csv`) | **partial** (Q1 still open)                  | `time_to_positive` and `mg_mm`/`salt_mm`/`primer_conc_nm`     |

### LAMP — new set (Zhang W et al. 2026)

`10.1093/nar/gkaf1378`, `PMC12774662`, _Nucleic Acids Res_ 2026, "A universal
kinetic framework for quantitative isothermal amplification". Tables 1 and 2 read in
the body text (both are printed, not plotted).

- One row per cell of 3 genes (N 145 bp, N1 175 bp, S 151 bp) × {63, 65, 67} °C ×
  {[E], [E]/2} × {clean, SVNA, dSVNA} = 24 rows.
- Outcomes: doubling time `T` (min, Table 1) and `T_θ`, `N₀` (copies/µL), `K_A`, `ν`,
  MSE (Table 2 — printed for the 63–65 °C blocks only, hence 12 non-`NA` rows).
- **Why it is still partial.** It is the first LAMP table read that puts copies/µL
  _and_ a kinetics parameter _and_ temperature in one grid, but the outcome is a
  **doubling time**, not a time-to-positive, and `N₀` is a _fitted_ quantity, so the
  copy axis is a model output rather than an independent condition. Mg²⁺, KCl and
  primer concentrations are not printed: the methods defer to Supplementary Table S2
  (Bst 2.0 WarmStart, NEB standard protocol), which was not read. Table 1 and Table 2
  are the same curves fitted twice, so their rows are not independent observations.
- Pre-existing LAMP blockers are unchanged (`copies_per_reaction = NA` in
  `readout_lamp_erysiphe2026.csv`).

### CRISPR-Dx — still none

Checked this pass and rejected (figure panels only, no printed values): `PMC13696524`
(Adv Sci 2026, phosphate reporters, 0 `<table-wrap>`), `PMC12802953` (NAR 2026, steric
trans-cleavage regulation), `PMC8667580`, `PMC10376970`, `PMC10452489` (its one table
is a literature LOD comparison), `PMC13129545` (NAR 2026, sequence-background — its
Table 1 is the gRNA list, not kinetics), `PMC13320335` (iScience 2026, hCMV — the
buffer × temperature trans-cleavage rate grid is **Table S3**, a supplementary file
reached only through the DOI, not read; the body prints only K_m = 1760 nM and
V_max = 2.16 s⁻¹ for one optimized condition). `PMC13320335` is the closest to a set
and is the next thing to open if the supplement link is unblocked.

### Virtual gel — still none

Queries run this pass: `TITLE:"agarose gel" AND ("migration distance" OR "migrated
distance")` (0 hits); `"migration distance" AND "agarose gel" AND fragments` (290
hits, top 10 checked, none a ladder calibration); `("band migration" OR "migration
distance") AND ("agarose concentration" OR "gel concentration") AND ladder` (5 hits,
all SDS-PAGE or immunoblotting). Ziraldo 2019 (`PMC6895257`) still supplies the
_formula_ with gel %, voltage and duration but no numeric bp→mm table.

## Skipped / risk

- **2026-11-03 pass:** `PMC13320335` Table S3, `PMC12774662` Supplementary Table S2
  and every other proof-of-work-gated supplement were **not read**; the LAMP recipe
  columns and the CRISPR buffer grid are empty for that reason, not because the
  papers lack them.
- `PMC12774662` Table 2 was read from the JATS XML body, where the last (67 °C) block
  is cut off; its `T_θ`/`N₀` cells were left `NA` rather than guessed.
- No `dev/check.R` run: no R code, no `inst/extdata/` file and no `.csv` read by the
  package was touched, so the gate has nothing new to check.

- `zotero_semantic_search` unavailable: the library sweep is keyword-based, so a
  relevant item whose abstract uses different wording could be missed.
- Kalendar 2022, Oliveira 2020 and Huguet 2017 were read through Europe PMC, not the
  Zotero attachment; the library records for Oliveira (`MMARQH5C`) and Huguet
  (`TZQ33APF`) return no fulltext through the lane.
- Southern 1979 and Elder 1983 are cited second-hand from Ziraldo 2019 and were not read.
- The Cas12a numbers in Nalefski 2021 are _in vitro_ RNP trans-cleavage; the numbers in
  Huang 2022 and Mao 2024 are in-tube assay endpoints. A model mixing both would be
  fitting two different systems.
