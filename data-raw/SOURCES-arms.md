# Sources for the allele-discrimination "sci-arms" sweep

Per paper: what was fetched, what the load-bearing assertion is, what is NOT
modelled, and any recorded contradiction. Everything written from these sources
is **uncalibrated**. Contradictions between sources are recorded, never averaged.

Scratch downloads live under `.tmp/lit/<key>/` and are not committed.

---

## desilva2009 — Garritano et al., BMC Genet 2009 — EXTRACTED

**Provenance correction (important).** The backlog row carried
`doi:10.1186/1471-2156-10-70`, which resolves to a different paper
(Hallander & Waldmann, tree-breeding selection contribution). The TP53 HRM paper
is **`10.1186/1471-2156-10-5`** (PMID 19173742, PMC2648999), and **De Silva D is
a middle co-author, not first author**. Corrected in `scientific_backlog.csv`
and in the script header.

- **Fetched:** `fulltext.xml` (Europe PMC REST, JATS) and
  `supplements/1471-2156-10-5-S1.doc` (Additional file 1). The BMC and PMC PDF
  endpoints serve an HTML interstitial in place of the PDF, so no PDF was used
  and none is needed. The `.doc` is legacy Word and is read with `strings`, not
  a Word library; the text streams are plain ASCII in table order.
- **Files written:** `desilva2009_amplicons.csv` (22), `desilva2009_snps.csv`
  (25 primer rows, 4 multi-SNP amplicons), `desilva2009_outcomes.csv` (2 pooled
  rows), `desilva2009_strategies.csv` (5).
- **Assertion:** Table 1's 22 sizes, split at the paper's 400 bp, must reproduce
  its own printed group statistics: **9 amplicons <400 bp with mean exactly
  286 bp**, and **13 at >400 bp**, i.e. 9+13 = 22. Plus: every printed fraction
  must repeat the counts we read from it (`(38 TP/(38 TP + 0 FN))`), and the
  recomputed percentages must match the printed ones at the printed precision.
- **Boundary finding:** Table 1 contains one amplicon at *exactly* 400 bp, which
  belongs to neither of the paper's groups as stated. Reading it as <400 gives
  9/13 and mean 286 — the printed mean — so it is assigned there, and the
  assignment is **inferred, not observed**. Under the other reading the paper's
  own 286 is unreproducible.
- **Recorded contradictions (unresolved):**
  1. `mean(>400 bp)` = 545 from Table 1 vs **544 printed** (0.2%).
  2. Printed percentages vs counts: >400 bp sensitivity 105/128 = 0.82 vs
     printed 0.81; specificity 339/408 = 0.83 vs printed 0.84.
- **NOT modelled:** per-amplicon TP/FP/TN/FN (the paper prints only two pooled
  rows; no per-amplicon contingency exists, so none is invented); the 45
  amplicons of the 67 outside Table 1; the 30 SNP × 270 sample genotyping arm
  (qualitative concordance only). Table 1 and Additional file 1 are
  **independent assay sets** — they share one forward primer used with a
  different reverse, and no amplicon pair in common — so neither corroborates
  the other's sizes. The script asserts that independence.

## ace1r_aslamp2015 — Badolo et al., Malar J 2015 — EXTRACTED

**Provenance.** The backlog had no DOI. `10.1186/1475-2875-15-66` does not
exist (Crossref 404). The real DOI, read from the fetched article's own
`<article-id>`, is **`10.1186/s12936-015-0968-9`** (PMID 26620269, PMC4665897).

- **Fetched:** `fulltext.xml` (Europe PMC). Supplements are figures only.
- **Files written:** `ace1r_aslamp2015_counts.csv` (15),
  `ace1r_aslamp2015_metrics.csv` (18),
  `ace1r_aslamp2015_time_to_positive.csv` (4),
  `ace1r_aslamp2015_contradictions.csv` (3),
  `ace1r_aslamp2015_design_rules.csv` (5).
- **Assertion:** the sequencing row re-derives the paper's 104 samples
  (1 + 9 + 94); **every method column also sums to 104**; the AS-LAMP metric
  cells must equal the six values printed in Table 1 (1 / 0.89 / 1 for
  sensitivity, 1 / 1 / 0.90 for specificity); and Table 2's matrix sums to the
  24 samples the text says were randomly drawn. The script also checks its own
  DOI against the fetched XML.
- **What this anchors for the package:** discrimination is a **time window**,
  not a binary call. Time-to-positive is stored with its censoring
  (non-target allele not amplified within 75 min → Tp > 75 min, not measured).
- **Recorded contradictions (unresolved):**
  1. Results text says PCR-RFLP found **3 RR / 9 RS / 95 SS** (total 107);
     Table 1's PCR-RFLP row reads **1 / 9 / 94**. Both are kept.
  2. Table 1's PCR-RFLP SS specificity (0.90) contradicts its own counts
     (94/94 = 1.00) and duplicates the neighbouring methods' 0.90 — most likely
     a source-table duplication, recorded not corrected.
  3. Only 14 of 18 printed metric cells follow from the counts under a
     one-vs-rest reading (the four exceptions: PCR-RFLP SS sensitivity 0.99 and
     SS specificity 0.90, RT-PCR RS sensitivity 0.75 and RS specificity 1.00).
     Those four are carried verbatim and flagged
     `matches_one_vs_rest = FALSE`; the rule is not assumed to hold for them.
  4. Paper-internal wording: the Discussion says the penultimate mismatch sits
     at the "3′ (BIP) or 5′ end (FIP)" while the Results and Fig. 1 legend put
     it at the penultimate position. The Results/legend wording is the one
     recorded; the Discussion sentence is inconsistent with them.
- **NOT modelled:** the AS-LAMP **primer sequences**. They exist only in
  Figure 1, a raster image, so nothing is transcribed from it. The paper's
  design statements (mutation on the 5′ end of FIP/BIP; **an additional
  mismatch at the penultimate position**; outer primers identical between
  sets; a 3′-BIP design with one extra mismatch *failed*) are recorded verbatim
  as text, tagged uncalibrated, without sequences.

## ye2001 — Ye et al., Nucleic Acids Res 2001 — READ, NOT EXTRACTED

- **Fetched:** the PMC55900 HTML article page (Europe PMC `fullTextXML`
  returns HTTP 500 for this PMCID and the OA package endpoint 404s). Table 1,
  Table 2 and the full body text are readable from that page.
- **Why no new CSVs:** the full text was read and its rules recorded, but no
  **algebraic** assertion could be constructed. The only numbers are (a)
  author-reported Tm, annealing temperature, Mg²⁺ and amplicon sizes — a
  transcription already held as `ye2001_primer_set`, with no independent total
  to re-derive; and (b) the inner/outer amplicon sizes, which cannot be
  checked without the reference sequence and the outer-primer distances, which
  the paper does not give. Writing a CSV here would mean asserting a rule with
  nothing to verify it against, so nothing is written.
- **Rule, verbatim:** "a second deliberate mismatch at **position −2 from the
  3′-terminus**", from ref. 8: strong (G/A, C/T) 3′ mismatch needs a weak
  (C/A, G/T) second mismatch and vice versa; medium (A/A, C/C, G/G, T/T) needs
  medium; inner primers 26 nt or longer; inner:outer primer ratio 10:1.
- **CONTRADICTION (recorded, NOT resolved, NOT averaged):** Ye's −2 position
  versus Liu 2012's "3rd base from the 3′ end" is best versus Lefever's
  position 3 (= 4th nt). Three different sources, three different optima.

## badolo2012 — Badolo et al., Malar J 2012 — NOT EXTRACTABLE

- **Fetched:** `fulltext.xml` and both supplements (S1 `.doc`, S2 `.rtf`).
- **Finding:** the paper is an insecticide-resistance monitoring report. Its
  full text contains **no LAMP content at all** — the only "B2" matches are
  inside Table 4's ace-1R allele frequencies. Tables 1–5 are bioassays and
  kdr/ace-1R frequencies; S1 and S2 carry no primer or LAMP material.
- **Conclusion:** the backlog description ("LAMP-SNP, mutation at the 3′ end of
  B2 plus a second mismatch, kdr-w") does **not** describe this paper. Nothing
  was written. The LAMP second-mismatch rule we actually hold comes from
  `ace1r_aslamp2015`, which cites Badolo for it.
- Set `unfetched` → `fetched-not-extractable` in the backlog: fetched, opened,
  and found to contain nothing extractable.
