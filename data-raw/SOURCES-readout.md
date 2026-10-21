# SOURCES readout — scientific extraction, branch `feat/sci-hrm`

Per paper: what was actually fetched, the load-bearing assertion each
`data-raw/*-supplement.R` passes, what is **not** modelled, and the outcome
definition quoted verbatim from the source.

**Identifier resolution is provenance, not evidence.** A DOI proves a paper
exists and let us fetch it; it licenses no number. Every value in
`inst/extdata/` comes from a file that was downloaded and read, and is
re-derived by the script that wrote it.

Test discipline: each script was verified only with
`Rscript -e 'source("data-raw/<key>-supplement.R")'`. The full suite was not run.

---

## 1. `jak2v617f_ddpcr` — JAK2 V617F ddPCR

- **Citation / DOI**: Liu et al., *Global Medical Genetics* 2024;11:132,
  `10.1055/s-0044-1785537`.
- **Fetched**: `.tmp/lit/jak2v617f_ddpcr/fulltext.xml`. The supplementary
  archive (`supplements.zip`) contains **figures only** (`.gif`/`.jpg`) — there
  is no numeric supplement, so nothing here depends on one.
- **Written**: `jak2v617f_ddpcr_precision.csv` (16 cells = 4 target VAF ×
  4 CLSI EP05-A3 conditions), `jak2v617f_ddpcr_limits.csv` (LoB, LoD, LoQ as
  three distinct definitions).
- **Assertion**: `CV = 100·SD/Mean` reproduces within 0.005 p.p. in
  **10 of the 12** cells with VAF ≥ 1%; the four *total-precision* cells
  reproduce the prose's SD range **0.03–0.36%** and total-precision CV range
  **0.25–20.17%**; and the four 0.1%-VAF cells must be **off** by more than the
  tolerance (checked, so the gap is recorded, not silently smoothed).
- **Not modelled** (all stated in the script header):
  - Table 2 is **not** a partition-model record and does not validate the
    Poisson model. Poisson droplet statistics predict a different CV; that gap
    is the reason to keep the table.
  - The four 0.1%-VAF cells **do not satisfy** `CV = 100·SD/Mean` (printed
    6.5 / 23.916 / 16.79 / 20.174; recomputed 6.748 / 24.000 / 16.667 /
    19.863). At a mean of ~0.15 on a 2-d.p. SD the identity cannot survive
    rounding, so we cannot tell whether a different estimator was used. Printed
    values are kept unchanged.
  - Two 1% cells miss by 0.01–0.03 p.p. (5.046 vs 5.044; 3.911 vs 3.925) and
    are also kept as printed.
  - The per-level LoD data and the WHO 7-level accuracy set exist only as prose
    counts and a **figure** (Fig 3); the per-level read-outs are not tabulated,
    so no values were copied.
  - The 39-patient ddPCR-vs-qPCR comparison is aggregate only
    (r = 0.988; bias +4.12%, LoA −5.76 to 14).
  - No CV is converted into a per-call probability. **The 76% CV at the LoD is
    carried verbatim**, not tidied.
- **Verbatim outcome definition** (precision): the paper reports
  "CV ranging from 0.25 to 20.17% (total precision in Table 2)" and
  "SD ranging from 0.03 to 0.36% (total precision in Table 2)".
  For the limits: "LoB value of 0.002% VAF, 0.043 copies/μL, with a minimum
  requirement of more than one positive droplet"; "an LoD of 0.005% VAF with a
  CV of approximately 76%"; "LoQ … as low as 0.01% VAF, with a concentration of
  0.193 copies/μL in a 20-μL system".

---

## 2. `snap2caps` — SNP2CAPS

- **Citation / DOI**: Thiel, Kota, Grosse, Stein, Graner, *Nucleic Acids
  Research* 2004;32:W205, `10.1093/nar/gnh006` (PMC373308).
- **Fetched**: Europe PMC `fullTextXML` for PMC373308 returns **HTTP 500**
  (kept as `fulltext.epmc-500.json`); OUP is 403. The working artifact is the
  PMC **HTML** page (`fulltext.html`) and its cells as `fulltext.txt`.
- **Written**: `snap2caps_markers.csv` (14 rows), `snap2caps_accession_groups.csv`
  (29 rows).
- **Assertion**: exactly **14** marker blocks / 14 distinct marker × enzyme
  pairs; and for **every** accession group the enumerated confirmed fragments
  **sum to that group's PCR product size** (e.g. 204 + 131 = 335; 157 + 106 + 88
  = 351). This is the paper's own claim, so a mis-transcribed fragment number
  fails the run.
- **Not modelled**:
  - The outcome is **prediction confirmation**, not conversion: the paper says
    "all 14 EST–restriction enzyme pairs … revealed the predicted restriction
    pattern". Table 2 deliberately includes accessions predicted **not** to be
    cut; those uncut lanes are part of the confirmation, **not failures**.
  - The 14 were randomly chosen *after* pre-screening with common enzymes, so
    this is not a conversion rate over arbitrary SNPs.
  - SNP2CAPS itself does not judge gel resolvability — the paper warns that
    some predicted sites are too close to a fragment end or to another site.
  - Which accessions are cut is recorded but **unasserted** (no independent
    second count exists in the text).
  - Distinct from the 2021 VCF2CAPS paper — different program, never pooled.
- **Verbatim outcome definition**: "Upon digestion of the corresponding PCR
  fragments, all 14 EST–restriction enzyme pairs that were tested revealed the
  predicted restriction pattern." Table 2 footnote: "In all cases the
  predictions of SNP2CAPS could be confirmed."

---

## 3. `vcf2caps2021` — VCF2CAPS

- **Citation / DOI**: Wesołowski et al., *PLoS Comput Biol*
  2021;17(5):e1008980, `10.1371/journal.pcbi.1008980`.
- **Fetched**: `.tmp/lit/vcf2caps2021/fulltext.xml`.
- **Correction to the registry row.** The backlog says this paper "reports no
  wet-lab outcomes". **That is false.** Its Results contains "Experimental
  verification of the identified CAPS markers": 10 single-cut markers assayed on
  30 plants each (Fig 8), with a per-marker verdict in prose.
- **Written**: `vcf2caps2021_verification.csv` — 10 rows, per-marker verbatim
  verdict, `n_plants_per_marker = 30`, and **`extractable = FALSE` with the
  reason**. No count was written.
- **Assertion**: the ledger's marker keys equal **exactly** the 10 abbreviated
  ids of Table 1 ("CAPS markers selected for experimental verification").
- **Not modelled — and why no rate exists**: the per-marker result is a
  qualitative prose verdict ("exactly matched genotypic predictions", "the only
  discrepancy was … a single supposed heterozygote", "five supposed
  non-cleavable homozygotes were in fact heterozygotes"). There is **no
  per-plant call table** and no supplementary data table; Fig 8 is a gel image.
  Any "x/30 correct" would be invented, so none was computed.
- **Verbatim outcome definition**: "Each of these markers was tested on 30
  plants from population 398 – 10 non-cleavable homozygotes, 10 heterozygotes
  and 10 cleavable (by the restriction enzyme devised for a given marker)
  homozygotes (Fig 8). … so were the sizes of the derived restriction
  fragments."

---

## 4. `soybean_est_caps` — soybean EST-SNP → CAPS

- **Citation / DOI**: 束永俊 et al., *作物学报 (Acta Agronomica Sinica)*
  2010;36(4):574–579, `10.3724/SP.J.1006.2010.00574` (printed on p. 1).
- **Fetched**: `.tmp/lit/soybean_est_caps/paper.pdf` (6 pp) → `paper.txt`
  (`pdftotext -layout`). The paper has **no supplementary files**; the paper is
  the dataset.
- **Written**: `soybean_est_caps_markers.csv` (26 primers + enzyme),
  `soybean_est_caps_funnel.csv` (537 → 48 → 44 → 36 → 26).
- **Assertion**: Table 1 yields exactly **26** primer rows with ids contiguous
  `ESNP01…ESNP26`; the funnel is strictly decreasing and nested
  (26 ≤ 36 ≤ 44 ≤ 48 ≤ 537); every enzyme is one of the four named in Methods
  (HindIII, SacI, VspI, PstI).
- **Not modelled**:
  - **26/48 is not a success rate.** The stages discard different things:
    48→44 is amplification failure, 44→36 is SNP-calling failure on Sanger
    traces, 36→26 is **digest-scoring** failure. The losses are not
    interchangeable and are never pooled with another paper.
  - Only the **survivors** are named. The paper does not say *which* 22 of 48
    failed, nor which 8 of 44 gave no SNP, nor which 10 digests were ambiguous —
    per-stage "reasons" would be invented.
  - Table 1 has primers + enzyme only: no amplicon size, no digest fragment
    sizes, no per-variety genotype calls, so **no** predicted-vs-observed
    cross-check is possible here (unlike SNP2CAPS).
- **Verbatim outcome definition**: "获得的 26 个 CAPS 标记" / "the products of
  26 pair primers with **unequivocal restriction pattern** were identified as
  CAPS markers"; abstract: "酶切分析发现 26 个 PCR 产物具有酶切多态性".
  (The English abstract's "36 discordant ones" refers to the Sanger traces
  showing the expected EST-SNP — a different stage.)

---

## 5. `liu2016_caps` — watermelon CAPS

- **Citation / DOI**: Liu et al., *Breeding Science* 2016;66(2):244–259,
  `10.1270/jsbbs.66.244` (PMC4785002).
- **Fetched**: `.tmp/lit/liu2016_caps/fulltext.xml`. It has **no**
  `<supplementary-material>` node and no "Table S"/"suppl" mention — the
  article's tables are the whole record.
- **Written**: `liu2016_caps_funnel.csv` (14,246 → 3,192 → 532 → 516 → 296 →
  271), `liu2016_caps_snp_types.csv` (Table 2).
- **Assertion** (three linked checks):
  1. Table 2's three family totals re-summed from their components:
     **72,285 / 153,408 / 19,268**;
  2. those sum to **244,961**, i.e. exactly Table 1's genome-wide Total
     (225,693 SNPs + 19,268 indels), tying the two tables together;
  3. the funnel is monotone, and the paper's own "polymorphic rate of 50.94%"
     reproduces from 271/532.
- **Not modelled**:
  - **There is no per-primer table.** Tables 1–7 are SNP density, SNP types,
    linkage map, phenotype means, QTLs, candidate-gene annotation and trait
    correlations — none lists the 532 primer pairs or their outcomes. Only
    stage totals exist, so no per-marker CSV was written.
  - The paper's "approximately 42.64% … amplified additional fragments or
    failed in enzyme digestion" **has no stated denominator** and does not
    complement 50.94% over 532; the two percentages cannot both be right
    against the same denominator. We record the counts and flag the conflict
    rather than resolving it.
  - The 301 mapped markers / 1836.51 cM are linkage-map facts, not conversion
    outcomes, and are excluded.
- **Verbatim outcome definition**: "271 pairs of primers revealed expected
  fragment lengths and polymorphisms, with a polymorphic rate of 50.94%,
  according to the results of the restriction enzyme reaction" — i.e. a
  distinct band of the expected length **and** polymorphic across P1 / P2 / F1.

---

## Not obtained / nothing written

None of the five is missing. Note the two **corrections to earlier assumptions**:

| Registry claim | Actual |
|---|---|
| VCF2CAPS "reports NO wet-lab outcomes" | It does report a 10-marker wet-lab verification — in prose only, hence `extractable = FALSE` |
| SNP2CAPS "14 markers experimentally validated" | The paper's claim is confirmation of **predicted restriction patterns** (14/14), which is a different outcome |
