# Tm oracle sources — what was fetched, and why no CSV was written

Status: **no experimental Tm table was obtained.** Nothing in this file is a
number copied from a paper we did not open; every claim below is tied to a URL
we actually requested and a response we actually read. No file was written to
`inst/extdata/` for any of the three sources.

The purpose was an **independent Tm oracle** for `primer_tm()` / the TmCalculator
engine. The best available anchor remains `tests/testthat/test-tm-oracle.R`,
which states its tolerance **per nearest-neighbour table** (same-lineage
agreement, cross-lineage offsets, monotonicity) — not a blanket |ΔTm| ≤ 2 °C,
which the measured offsets (Sugimoto +1.5…+7.3 °C, Breslauer +8…+18 °C on 32
Makhoul primers) make false.

---

## 1. `vonahsen2001` — resolved, but the table does not exist as a table

| | |
|---|---|
| DOI | `10.1093/clinchem/47.11.1956` |
| PMID / PMCID | 11673362 / **none** |
| Journal | Clin Chem 2001;47(11):1956–1961 |
| Fetched | Crossref `works/10.1093/clinchem/47.11.1956` → HTTP 200 |
| Fetched | Unpaywall `v2/10.1093/clinchem/47.11.1956` → HTTP 200 |
| Fetched | `Rscript data-raw/fetch-sources.R --only vonahsen2001,dnamate` → `crossref oa=FALSE`, assets **0/2** |

Confirmed from the fetched Crossref record (title, authors von Ahsen / Wittwer /
Schütz, volume/issue/pages above). The fetched abstract states, verbatim:

> "The experimental Tms of 475 matched or mismatched target/probe duplexes were
> obtained in our laboratories **or were compiled from the literature** based on
> studies using the same real-time PCR platform."
> SE 1.76 °C for perfect matches (n = 217); 2.06 °C for all duplexes (n = 475).

Two reasons there is no CSV:

1. **The 475 are not a published per-row table.** They are a compiled set of
   duplexes; the abstract describes an accuracy study (standard errors), not a
   sequence × Tm listing. Transcribing 475 Tms would require the full text
   *and* the authors' compilation.
2. **The full text could not be read.** Unpaywall on the fetched record:
   `is_oa: false`, `oa_status: closed`, `has_repository_copy: false`,
   `oa_locations: 0`. Every OUP and PMC PDF endpoint returned **HTTP 403 with a
   Cloudflare "Just a moment..." interstitial**. That is bot-blocking, not
   simply a paywall — but either way no document reached us.

**Trap recorded, because it nearly shipped:** `fetch-sources.R` left
`.tmp/lit/vonahsen2001/fulltext.pdf`, 5737 bytes. It is **not a PDF** — `file(1)`
reports `HTML document`, and `pdftotext` fails on it. It is the Cloudflare
challenge page saved under a `.pdf` name. Any future agent must check the magic
bytes before treating an asset as a paper; a non-empty file is not a fetched
paper.

## 2. `dnamate` — resolved, but it is not the paper that holds the 348 Tms

The backlog identifier was the placeholder `"dnaMATE benchmark"`, with no DOI.
Resolution:

| | |
|---|---|
| Identifier | **dnaMATE** server paper: Panjkovich A, Norambuena T, Melo F. NAR 2005;33(suppl_2):W570–W572 |
| DOI / PMID / PMCID | `10.1093/nar/gki379` / 15980538 / **PMC1160140** (OA) |
| Fetched | Europe PMC full text → `.tmp/lit/dnamate/fulltext.xml`, 28389 bytes, HTTP 200 |
| Fetched | `supplements.zip`, 39809 bytes, extracted to `.tmp/lit/dnamate/supplements/` |

What the fetched full text contains:

* `table-wrap` elements: **0**. Legacy `<table>` elements: **0**.
  `<supplementary-material>` references: **0**.
* The only extracted supplement files are **`gki379f1.gif` / `gki379f1.jpg`** —
  Figure 1. There is no data table of any kind.
* The string "348" occurs exactly once, inside a run-together digit sequence
  (`2003314834851252005710`, a volume/year/page blob) — **not** a result.

So dnaMATE yields **no Tm table**. The 348 Tms belong to dnaMATE's ref. 1:

| | |
|---|---|
| True home | Panjkovich A, Melo F. "Comparison of different melting temperature calculation methods for short DNA sequences." Bioinformatics 2005;21(6):711–722 |
| DOI / PMID / PMCID | `10.1093/bioinformatics/bti066` / 15501913 / **none** |
| Fetched | Europe PMC `EXT_ID:15501913` (core) → HTTP 200 |
| Fetched | Unpaywall → `is_oa: true`, `oa_status: bronze`, `has_repository_copy: false`, 1 publisher location |
| **Not fetched** | the paper itself (Europe PMC `isOpenAccess: N`, `inEPMC: N`, `hasSuppl: N`; OUP link is Cloudflare-gated) |

The fetched abstract, verbatim:

> "Using a total of 348 DNA sequences in the size range between 16mer and 30mer,
> for which the experimental Tm data are available, we demonstrated that the
> consensus Tm is a robust and accurate measure."

Read carefully: the 348 sequences are ones **"for which the experimental Tm data
are available"** — i.e. **compiled by the authors from elsewhere**, not measured
by them. So this is a benchmark, not a primary measurement set, and we have not
opened the paper that lists it.

## 3. `huguet2017` — correct DOI is `gkx1161`, not `gkx447`

| | |
|---|---|
| DOI | **`10.1093/nar/gkx1161`** (matches `source_resolution.csv`) |
| Paper | Huguet, Ribezzi-Crivellari, Bizarro, Ritort. "Derivation of nearest-neighbor DNA parameters in magnesium from single molecule experiments." NAR 2017;45(22):12921–12931; PMC5728412 |
| Fetched | `fulltext.xml` (copied read-only from the main worktree, 155056 bytes) |
| Fetched | `supplements.zip` re-downloaded from Europe PMC; `supplements/gkx1161_supp.pdf` |

**Do not use `10.1093/nar/gkx447`.** That DOI is a different paper (SAGA
complex biology), not Huguet. The backlog and `source_resolution.csv` both say
`gkx1161`; that is the one used here.

What is in it: the main text carries **nearest-neighbour parameter tables**
(ΔH, ΔS per motif), not Tm measurements. The only Tm data are in the
supplement, §"independent validation", verbatim:

> "we have challenged our newly derived set of energy values with an independent
> set of 22 oligos at various salt conditions and a total of 58 melting
> temperatures."

and Figure S10's caption, verbatim:

> "The first 12 sequences were obtained from Table 4 from Ref. [16] at a total
> oligo concentration of CT ~ 2 µM. The rest 10 sequences are given in Table S2
> of Ref. [16] …"

Ref. [16] is **Owczarzy et al. 2008**. Those 58 Tms exist only as a **raster
table rendered inside a figure image**. Recovering them would mean OCR of a
figure, or fetching Owczarzy 2008 — both outside this task. OCR output is not an
acceptable source for `inst/extdata/`.

---

## Why no `vonahsen2001-supplement.R` / `dnamate-supplement.R`

The style of this directory is: a `*-supplement.R` script **asserts** that the
file it parsed reproduces a number printed in the paper, and only then writes a
CSV. Here there is no parsed table to assert against:

* von Ahsen: no document reached us (403 Cloudflare / closed access), and the
  475 Tms are a compilation, not a printed per-row table.
* dnaMATE: the paper has zero tables and a figure-only supplement; the 348 Tms
  live in a different, unopened paper and are themselves compiled.
* Huguet: 58 Tms exist only as pixels in Figure S10.

Writing any CSV from these would mean typing numbers from memory or from a
raster — precisely the invented-citation failure this file exists to prevent.
Per the task rule, absence is recorded here instead.

## What would unblock it

1. Institutional access to Clin Chem 2001;47(11) (von Ahsen) **and** the
   authors' supplementary compilation of the 475 duplexes — only useful if such
   a compilation was actually published.
2. The full text or supplement of Bioinformatics 2005;21(6):711–722
   (`bti066`), which lists the 348 sequences **and** their source references.
3. Owczarzy et al. 2008 Table 4 and Table S2, which would give Huguet's 58
   validation Tms from their primary source.

Any of these, once fetched and parsed, should arrive as a `*-supplement.R` with
a load-bearing assertion — and the resulting test must keep a **per-table**
tolerance, never a blanket one.
