# Tm oracle status — what an external Tm test can honestly assert today

ZW3 anchors `primer_tm()`. The rule is: a number lands only if it was read from a
downloaded file **and** re-derived by a `data-raw/` script that stops on disagreement.
A DOI is provenance, not evidence.

## Landed

### `mm2020` — the only measured, tabular, oligo-level set found
- File: `inst/extdata/tm_mm2020_s1_measured.csv` (4096 rows). Generator `data-raw/mm2020.R`.
- 19-mer duplexes, all sharing `5'-CGACGTGCN1N3N5ATGTGCTG-3'`; only the central trimer varies.
- Buffer **stated**: 50 mM NaCl + 10 mM Na-phosphate pH 7.4, 1.0 µM total strand, **no Mg²⁺**.
- Measured (LightCycler first-derivative maximum), not computed.
- Asserted: 4096 rows, `i` a permutation of 1..4096, classes exactly
  64/576/1728/1728 (matches the prose), the 64 canonical rows are exactly the
  complement-of-top rows and are written head-to-head, rebuilt duplexes are 19 nt,
  measured Tm non-increasing in `i`.
- **What the model does**: r = 0.958, Spearman = 0.941, but mean residual −16.86 °C
  (sd 0.64) on the canonical 64. This is a **level offset**, not a ranking failure.
- **Mg is not the explanation.** Raising Mg from 0 → 3 mM closes only ~12 °C
  (bias −16.86 → −5.04) and then makes the mismatch classes *worse*
  (single −9.88 → +3.9). The offset is real and is **not corrected anywhere**.
- Therefore the tests assert **ranking, spread, and cross-class separation**, and
  explicitly assert that the bias stays in (−25, −10) °C so a silent correction fails.

### `hm2025` — amplicon-level, directional only
- File: `inst/extdata/tm_hm2025_amplicon_psaa2.csv` (5 rows). Generator `data-raw/hm2025.R`.
- psaA-2 primer pair, 5 diatom species; measured HRM amplicon Tm vs the
  **authors' own** predicted Tm. All 5 printed differences recompute exactly.
- **Cannot test `primer_tm()`**: these are ~100 bp amplicons, not primer duplexes.
- Dia18S half excluded — its primers contain `N` (degenerate).
- Buffer gap: 2× Taq mix + Eva Green; Na⁺/K⁺/Mg²⁺/dNTP/DMSO **not stated numerically**.
- A 30-point Mg²⁺ series (S2 Table, 1.5–5 mM) exists and its 30 values were confirmed
  present in the DOCX, but it was **not landed**: Mg²⁺ is the only varied ion and the
  rest of the buffer is unquantified, so it cannot be re-derived as a clean condition set.

## Blocked / pointer only — nothing landed

### `vonahsen2001` (Clin Chem 2001, `10.1093/clinchem/47.11.1956`)
- **Closed access** (Crossref `oa=FALSE`, EPMC `isOpenAccess=N`, no PMCID). No bypass attempted.
- The abstract was read and confirms: the 475 experimental Tms were *"obtained in our
  laboratories or were compiled from the literature"* — **not printed as a table**.
  Data are available on request.
- **Can assert nothing numerically.** Cannot be cited as a validated reference set.

### `huguet2017` (NAR 2017, `gkx1161`, PMC5728412)
- OA; full text and supplement fetched. The 58 validation Tms are in **Figure S10,
  a raster image** — no text layer, so no number was transcribed.
- Tables S2–S6 are **NN parameter sets, not Tm sets**.
- It is a *pointer* to Owczarzy 2008 (`10.1021/bi702363u`), which has no PMCID and is
  not OA. **Nothing landed.**
- **Mg²⁺ consequence**: Huguet's thesis is that Mg²⁺ changes the parameters, but its
  supporting data are unreadable here, so the package's Mg²⁺ handling remains
  **unvalidated against an external measured set**.

### `dnamate` (NAR 2005, `10.1093/nar/gki379`, PMC1160140)
- OA full text read: **0 Tm tables**. The 348 measured Tms belong to its reference 1
  (*Bioinformatics* 2013, `10.1093/bioinformatics/bdt066`), not fetched.
- **Nothing landed.**

## The rule for future work
`primer_tm()` may now be described as *rank-validated* against a measured set at
sodium-only conditions, and **must not** be described as reproducing measured Tm in
absolute terms, nor as validated at 1.5 mM Mg²⁺. The −17 °C offset is an open defect
to investigate, not something to paper over with a constant.
