# ZW3.4 — can a published measured Tm set anchor the package's Tm scale?

**Verdict: no.** A measured, buffer-stated oligo-Tm set exists and is already in
the tree (`tm_mm2020_s1_measured.csv`), and `primer_tm()` ranks it almost
perfectly — but it sits **-9.33 °C** low under the package's own best salt
model with the residual a _level_, not a slope. A level offset cannot be
absorbed by a calibration factor chosen from the same set, and no second
measured set on the designer's length range exists in the operator's library
to break the tie. ZW3.4 therefore stays recorded-and-not-corrected.

Read 2026-10-29 through the Merovingen read-only MCP lane (collection
`4RSE6KC2`, ADR-0037); every page of the mm2020 and hm2025 PDFs was read with
`zotero_read_pdf_pages`. Nothing here is from memory.

---

## 1. What an anchor has to satisfy

An "anchor" is a set of **measured** duplex melting temperatures that can fix
the _level_ of `primer_tm()`, i.e. turn it from a ranker into an absolute scale.
Four conditions, all necessary:

1. **Measured, not computed.** A UV-melt / first-derivative maximum from the
   paper, not a design tool's output. Calibrating against another model's Tm
   reproduces that model's bias instead of removing ours.
2. **Numeric buffer on the paper's own scale.** `[Na+]` (or K+) and `[Mg2+]` in
   mM, plus dNTP and DMSO if used, so the measurement can be put under the same
   `reaction_conditions()` the model reads.
3. **The designer's length range.** `design_as_primers()` works in 20–27 nt and
   the schemes add 17–25 nt RPA/Sanger primers; a 19-mer-only set measures a
   different duplex length, and an amplicon set measures a completely different
   object (ds product melt, not oligo duplex).
4. **A residual that is a slope, not a level.** Only a slope can be absorbed
   by a salt term; a level says the _reference state_ differs (what "1 M" the
   paper meant, how the strand concentration was set), and that is a
   correspondence with the authors, not a search.

## 2. Inventory — every measured-Tm candidate in the library

| source                                                                                                                                                                                                                               | what is measured                                                                                                                 | buffer printed                                                                                                     | length         | verdict                                                                             |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ | -------------- | ----------------------------------------------------------------------------------- |
| **mm2020** — Oliveira LM, Long AS, Brown T, Fox KR, Weber G, _Melting temperature measurement and mesoscopic evaluation of single, double and triple DNA mismatches_, Chem Sci 2020;11:11419, doi 10.1039/d0sc01700k, **PMC8163305** | 4096 measured duplex Tm, all central-trimer contexts; ESI **Table S1**                                                           | yes: **50 mM NaCl + 10 mM Na-phosphate pH 7.4, Mg2+ 0**, 1.0 µM total strand, LightCycler first-derivative maximum | **19 nt** only | **the only true anchor candidate — and it fails condition 3 and 4**                 |
| hm2025 — _Establishment of a simple prediction method for DNA melting temperature_, PLoS ONE, doi 10.1371/journal.pone.0321885, **PMC12002481**, **Table 9**                                                                         | 5 HRM amplicon Tms (psaA-2 pair, 5 diatom species)                                                                               | no: "2x MonAmp Taq mix + Eva Green", Na/K/Mg/dNTP not stated numerically                                           | ds amplicon    | direction only (condition 1 **and** 2 fail)                                         |
| dnamate — dnaMATE, doi 10.1093/nar/gki379, **IMULDMZ7**                                                                                                                                                                              | nothing: the paper prints no Tm table; its "348 measured oligo Tm, 16–30 nt" are _borrowed_ from a 2013 Bioinformatics benchmark | —                                                                                                                  | 16–30 nt       | the one set that satisfies condition 3 — **not read, no attachment in the library** |
| vonahsen2001 — doi 10.1093/clinchem/47.11.1956, **K3ZWDD8K**                                                                                                                                                                         | closed access                                                                                                                    | —                                                                                                                  | —              | registry calls it "the external Tm oracle ZW3 wants"; cannot be read                |
| kutyavin2000 — doi 10.1093/nar/28.2.655, **PMC102528**                                                                                                                                                                               | MGB-modified duplex Tms (a 15-mer MGB melts at 60 °C vs 42 °C unmodified)                                                        | partly                                                                                                             | 12–15 nt, MGB  | MGB chemistry the package does not model; short-mer                                 |
| huguet2017 — doi 10.1093/nar/gkx1161, **PMC5728412**                                                                                                                                                                                 | 58 validation Tms, a raster in Fig. S10 with no text layer                                                                       | yes (Mg series)                                                                                                    | —              | numbers were never transcribed; not usable as an anchor                             |
| jrc2025_h5nx — doi 10.2807/1560-7917.es.2025.30.33.2500183, **PMC12372895**, Table 1                                                                                                                                                 | 6 dPCR primer/probe Tms                                                                                                          | no                                                                                                                 | 22–25 nt       | **computed** (Primer3Plus), condition 1 fails; see `z0.2-dpcr-tm.md` §3             |

## 3. Why mm2020 fails conditions 3 and 4 — arithmetic, recomputed

`tm_mm2020_s1_measured.csv`, 64 canonical 19-mers, under
`reaction_conditions(Na = 50, Mg = 0, dNTPs = 0, DMSO = 0)`:

| quantity                                                   | value            |
| ---------------------------------------------------------- | ---------------- |
| measured Tm range                                          | 60.7 – 69.3 °C   |
| `primer_tm()` range (Schildkraut default)                  | 43.42 – 52.21 °C |
| `bias = mean(model - measured)`                            | **-16.8643 °C**  |
| sd of that residual (a **level**, not scatter)             | 0.6414 °C        |
| Pearson r / Spearman rho                                   | 0.9579 / 0.9413  |
| salt term `16.6·log10(0.05)`                               | -21.5971 °C      |
| **`bias - salt`** (what the 1 M reference state is off by) | **+4.7328 °C**   |
| `bias` with `salt_method = "Owczarzy2008"`                 | **-9.3348 °C**   |
| Owczarzy residual minus its own salt term                  | +12.2623 °C      |

Reproduce exactly as in `z0.2-dpcr-tm.md` §1 (same script, same file). The
rank is excellent — the model is a good _ordinal_ instrument on this set — and
the residual is ±0.64 °C about a fixed -9.33 or -16.86 °C. That is the
signature of a reference-state disagreement:

- Under the **default** salt model the salt term overshoots by +4.73 °C.
- Under **Owczarzy** — the more complete model — **-9.33 °C survives**, and its
  own `-salt` remainder is +12.26 °C. Mg is zero in this set, so Mg is excluded
  by construction (ROADMAP records the same rejection: 3 mM Mg would close
  ~12 °C and overshoot the mismatch classes).

Two measured sources that disagree on the reference state by ~9 °C cannot both
be right, and the paper does not state which "1 M" convention its comparison
uses. Anchoring means _adopting_ one of them; the honest action is to record
it, which is exactly what `calibration_stance = "recorded, not corrected: do
not calibrate into the minus"` does on every row.

## 4. What would turn this into an anchor

In order of cheapness:

1. **A second measured oligo-Tm set on 18–25 nt, buffer printed**, ≥ 8 duplexes
   spanning 50–70 °C. Two independent measured sets that _agree_ fix the
   reference state; that is a search, and it has not succeeded here.
2. **The mm2020 authors' stated measurement standard** — which 1 M reference
   and which strand-concentration convention they used. One correspondence,
   not a search, and it would convert the -9.33 °C from debt into a parameter.
3. **A measured, home-buffered dPCR/primer set** — needs a lab, and the vendor
   mixes are proprietary (see `z0.2-dpcr-tm.md` §2, §5).

Until (1) or (2): the scale is **ranked and uncalibrated**. Every downstream
window (the dPCR 44–50 °C window in `Z0.2`, the KASP 55–65/58–63 windows) is a
window on _this_ scale, not on a vendor's, and must be read that way.

## 5. Cross-references

- `Z0.2` — the dPCR Tm window is calibration debt on this scale:
  `data-raw/proposed/z0.2-dpcr-tm.md`.
- `ZW4.6` — the same literature pass that rejected a per-row source for the Tm
  windows: `data-raw/proposed/zw4.6-sources.md`, "Second pass, 2026-10-29".
- `VALIDATION.md` row "Oliveira et al., Chem Sci 2020" records the
  ordering-only framing; this note is the longer prose behind it.
