# ZW4.6 round 3 — one new named source (Lobato & O'Sullivan 2018), RPA rows

Literature pass, 2026-11-03. Round 1+2 (see `zw4.6-sources.md`) moved 16 rows to
`kwok1990` / `santalucia1998` and left 175 at `package assumption` across four
documented reason classes. This round searched the **Europe PMC open-access
full text** of the one RPA review the code comments already cite but the
registry never wired, plus a set of new term searches over the local full-text
scratch and the Zotero lane.

## Result up front

**One new named source, four rows.** Lobato IM & O'Sullivan CK,
_Recombinase polymerase amplification: Basics, applications and recent advances_,
Trends in Analytical Chemistry 2018;98:19–35, doi `10.1016/j.trac.2017.10.015`,
**PMC7112910** (Europe PMC `fullTextXML`, read in full this session). It states,
verbatim, the RPA primer-length, amplicon and probe-geometry numbers the package
uses. The source key `lobato2018` is **not yet** in
`inst/extdata/scientific_sources.csv` and must be merged there when these rows
are wired (the same merge caveat `zw4.6-sources.md` already records for
`santalucia1998`).

Every other of the 175 rows stays `package assumption`. Nothing new was found
for the Tm windows, the simulator thresholds, the inherited geometry, or the
CAPS / Sanger / SSR / LAMP / InDel scheme constants after the searches below.

## Rows that got a source

| criterion id | current source | proposed key | verbatim quote | DOI / PMCID | section | evidence tier |
| --- | --- | --- | --- | --- | --- | --- |
| `rpa:gate:gate_rpa_primer_len` | `package assumption` | `lobato2018` | "primers of 30–35 bases in length were necessary for RPA" | 10.1016/j.trac.2017.10.015 / PMC7112910 | §2.2.1 _Primer design_ | review; stated, nuanced |
| `family_rpa_primer_length` | `package assumption` | `lobato2018` | "primers of 30–35 bases in length were necessary for RPA" | 10.1016/j.trac.2017.10.015 / PMC7112910 | §2.2.1 _Primer design_ | review; stated, nuanced |
| `rpa:gate:gate_rpa_amplicon` | `package assumption` | `lobato2018` | "RPA can amplify targets up to 1,5 kb but is better suited to amplicons between 100 and 200bp" | 10.1016/j.trac.2017.10.015 / PMC7112910 | §2.2.1 _Primer design_ | review; stated |
| `rpa:gate:gate_rpa_probe` | `package assumption` | `lobato2018` | "The probe is recommended to be a 46–52 oligonucleotide … The abasic nucleotide … is placed at least 30 nucleotides from the 5′ end and 15 nucleotides from the 3′end" | 10.1016/j.trac.2017.10.015 / PMC7112910 | §3.1.1 _Lateral flow_ | review; stated (nfo probe) |

### Caveats a reviewer must see

1. **Primer length is nuanced.** The full sentence is "Whilst it was initially
   believed that specifically designed primers of 30–35 bases in length were
   necessary for RPA, there are several reports demonstrating that normal PCR
   primers can be used and efficient amplification achieved". The number 30–35
   is stated, but as a *historical* belief the review immediately softens. The
   package's `marker_rpa.R` comment already carries this ("Wikipedia allows to
   38; we use the review's 30–35 and note longer oligos are tolerated"). Citing
   the number is honest only with this caveat attached.
2. **The probe geometry is the nfo probe, not the exo probe.** The 46–52 nt
   length and the ≥30 nt 5′ / ≥15 nt 3′ abasic-arm rule are stated for the
   lateral-flow **nfo** probe (§3.1.1). The package applies the same geometry to
   the **exo** probe; `marker_rpa.R` marks that extension an ASSUMPTION. The
   fluorophore/quencher gap "separated by a maximum of 2–4 bases" is the exo
   probe (§3.2) but is a constant (`RPA_PROBE_FQ_GAP`), not a gate-limit bar.
3. **The categorical RPA rows are not covered.** `family_rpa_thf_off`,
   `family_rpa_probe_spans`, `family_rpa_variant_under`, and
   `rpa:gate:gate_rpa_variant` are placement properties, not numbers; Lobato
   2018 states none of them as a bar, so they stay `package assumption`.

## Why the other 171 rows still stay `package assumption`

The four reason classes from `zw4.6-sources.md` §"What remains unsourceable"
hold. New searches this round, none of which produced a number matching a bar:

- Zotero lane (`zotero_quote.py --grep`): `homopolymer`, `band gap`,
  `read distance`, `30-35`, `30 to 35`, `trinucleotide repeat`,
  `simple sequence repeat`, `2-4 nt`, `inner arm`, `5 °C`, `5 C`,
  `size delta`, `run of`. Hits were the already-known rejections
  (`I3EXS422` / `JHLSJPAH` "homopolymeric" with no length; SSR papers that
  *use* repeats without thresholding them).
- Local full-text scratch `/tmp/tidymarkers-scratch/lit/` (38 PMC XML) greps for
  `30-35`, `read distance`, `band gap`, `homopolymer`, `size delta`, `3 bp`,
  SSR copy thresholds — no marker-design bar stated.
- Lobato 2018 full text (the one genuinely new full text read this round) states
  only the RPA numbers above; it prints no CAPS band gap, Sanger read distance,
  SSR copy floor, LAMP inner-arm balance, or size delta.

Tm-window rows (`as_tm` 55–65, `com_tm` 58–63) and the simulator thresholds
(`discrimination` 0.5, dropout 0.1, `offtarget_hits <= 1`) remain unanchored for
the reasons already documented: no paper states *this* window under *this*
buffer, and no wet-lab paper can state a threshold on the package's own
simulator output.

## Access notes

`PMC7112910` read from Europe PMC `fullTextXML` (open access, Elsevier COVID
collection). No Zotero write, no record added. The key `lobato2018` needs a row
in `inst/extdata/scientific_sources.csv` and the four criterion rows need a
`criterion_sources()` entry before any merge; this document proposes the source,
it does not wire it.
