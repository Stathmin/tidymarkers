# ZW3.4 second pass — a measured duplex Tm set with printed salt, and why it is still not an anchor

Literature pass, 2026-11-03. Read through the Zotero lane (collection `4RSE6KC2`,
ADR-0037, PDFs read end to end) and Europe PMC `fullTextXML`. Prior state:
`data-raw/proposed/zw3.4-anchor.md` (mm2020 is the only measured 50 mM-NaCl set;
it fails the 18–25 nt and level-vs-slope conditions).

## 1. The two closest sets, and the exact condition each fails

| set | Tm: measured? | buffer printed? | oligo length | verdict |
| --- | --- | --- | --- | --- |
| mm2020 (PMC8163305, ESI Table S1) | yes, UV melt | 50 mM NaCl + 10 mM Na-phosphate pH 7.4, Mg 0 | 19 nt | length OK; residual is a slope |
| **Zhou 2025 (PMC12002481, item `F96SK3PY`)** | **yes**, HRM first-derivative peak | **Mg 1.5–2 mM recommended; Mg varied 1.5–3.5 mM (S2 Table)** | 60–152 bp **amplicons** | **length fails** — cannot test a primer model |
| Specht 2020 (PMC7260945) | yes, reconfiguration rate | 25 mM KCl, 10 mM MgCl2 | full-length | not an oligo Tm |

**New finding (Zhou 2025).** The paper measures Tm by HRM and prints what no other
paper in this collection prints: *the dependence of the measured Tm on [Mg2+]*.
It states the averaged model error is **1.77 °C** once Mg exceeds 2 mM and that
Mg must be held in **1.5–2 mM** for its formula to hold. It also prints the
absolute error of its own prediction (average **0.67 °C**, Table 9) and states
the residual grows with Mg.

That is the **first printed estimate of how much of a Tm-model residual is a
Mg term** in this literature. It is a bound, not a calibration: the objects are
60–152 bp amplicons, so it says nothing about an 18–25 nt duplex, and [Na+]/[K+]
are never printed.

## 2. Proposed rows

`data-raw/proposed/zw3.4-anchor.md` keeps its 4 admission conditions. No row is
proposed that passes them. What is proposable, against the **bias explanation**
(not the bias value):

| # | source | doi / pmcid | claim as stated | number | tier | location |
| - | ------ | ----------- | --------------- | ------ | ---- | -------- |
| A | Zhou 2025 | 10.1371/journal.pone.0321885 / PMC12002481 | "a significant augmentation in the average error is observed when the Mg2+ concentration surpasses 2mM … (with the average error reaching 1.77)" | 1.77 °C at >2 mM Mg | measured (HRM), amplicon scale | p. 13, Discussion |
| B | Zhou 2025 | same | "it is imperative to regulate the Mg2+ concentration in the PCR system to within the range of 1.5 mM to 2 mM" | Mg 1.5–2 mM | stated protocol | p. 13, Discussion |
| C | Kalendar 2022 | 10.3389/fpls.2021.747886 / PMC8784781 | Table 2A footnote: "50 mM K+, 3 mM Mg2+" (the buffer for the 55–56 °C ASP Tm already accepted in ZW4.6) | 50 mM K / 3 mM Mg | stated buffer for a computed Tm | p. 7 |
| D | Digital MCA 2025 | 10.3390/bios15010036 / PMC11764372 | "an average maximum difference of 1.1 °C" (ddPCR melt vs HRM on 6 amplicons) | 1.1 °C | measured paired comparison | Sec. 3.2, Table 1 |

D is the one that speaks to Z0.2 rather than ZW3.4: it shows the **readout**
contributes ~1 °C between platforms, i.e. two orders of magnitude below the
−16.86 °C bias.

## 3. Verdict — unchanged, with new negative evidence

The bias cannot be re-anchored from this literature, and the reason is now
demonstrated rather than inferred:

1. Every **measured oligo** set that prints a buffer is 19 nt or full-length;
   the one set in the dPCR/HRM field that prints *both* a measured Tm and a Mg
   series measures **amplicons** (Zhou 2025), which a nearest-neighbour primer
   model does not predict.
2. The only printed *quantitative* statement about the size of a Mg-driven model
   error is that it becomes **1.77 °C** above 2 mM Mg (Zhou 2025). The known
   bias is −16.86 °C at Mg = 0. A Mg term therefore **cannot** be the missing
   ~17 °C, in either sign or magnitude.
3. Set C is the useful by-product: it is the first **numeric buffer** attached to
   a sentence the package already cites (the 55–56 °C ASP Tm), and it shows the
   field's working K/Mg values (50 mM K, 3 mM Mg) are nothing like the 1 M Na+
   reference the Schildkraut salt term assumes.

**Recommendation.** Keep `flanking tm_range = c(44, 50)` and keep ZW3.4 open.
Cite A–C in the ZW3.4 note as the reason the residual is *not* a buffer-
correction bug; do not add any of them as a Tm-value anchor.

## 4. Access notes

`F96SK3PY`, `6D9UWLMW` (Kalendar 2022) read from Zotero PDFs, all pages.
`PMC11764372` read from Europe PMC fullTextXML. `PMC7269394` full text was
unavailable (HTTP 500); its numbers are cited from the abstract only and are
**not** used in the rows above.
