# TaqMan / 5'-nuclease family — status proposal

**Recommendation: `R/marker_taqman.R`'s mismatch claims are UNSUPPORTED BY DATA.**
The family is no longer empty of extracts (four open papers are landed, see
`.tmp/sdd/report-taqman.md`), but **no landed source measures mismatch
discrimination with a TaqMan/5'-nuclease readout**. Nothing that would
calibrate `TAQMAN_MISMATCH_DG`, or justify the centrality rule, was found in
the open literature.

## Claims with no data behind them

| location | claim | evidence found |
|---|---|---|
| `R/marker_taqman.R` L49-53 | a single mismatch costing <1.5 kcal/mol (dG37) "is not discriminated by the exonuclease step" | **none** — no open paper reports a dG threshold, and none reports a ΔCt between matched and mismatched probe |
| L58-63, L110 (`taqman_probe_candidates`, `design_taqman_marker`) | the probe is centred "so the mismatch is maximally destabilising" | **none, and contradicted as a universal rule**: in `citrus2020` (8 TaqMan-MGB assays, 384/384 calls made) the allele base is never central — centrality 0.80, 0.86, 0.95, 0.80, 0.93, 0.73, 0.80, NA — even though the probe pairs differ by one base in six of eight families |
| `scheme_def_taqman()` tier2 | probe Tm must exceed primer Tm | **rule only** (`znavsmgb2014`: "the probes should be longer than primers"), evidence_tier weak, no measurement |
| probe-length default 18:25 | — | 5 of 8 `citrus2020` FAM probes are 13-18 nt and still call all 48 cultivars; `znavsmgb2014` uses 17/19-nt probes |
| `apol1_2019` (landed) | the modelled failure mode | a real, observed failure: a 6-bp deletion allele removes probe bases, so 1.4% of 1,110 samples do not amplify and 5.9% are mis-typed — a failure mode the current design does not model |

The one open source that states the trade-off directly (`lna2004`,
`yao2006`, `stadhouders2010`) or measures it (`yao2006`'s five-mismatch
tolerance, `stadhouders2010`'s position dependence) is closed:

* `yao2006` doi:10.1016/j.mcp.2006.03.003, PMID 16704921 — EPMC core:
  `inEPMC:N, inPMC:N, hasPDF:N, isOpenAccess:N`; the only `fullTextUrlList`
  entry is `Subscription required` (doi.org). No PMCID.
* `stadhouders2010` doi:10.2353/jmoldx.2010.090035, PMID 19952061,
  PMC2797725 — PMC OA service → 404 HTML; EPMC `fullTextXML` → HTTP 500;
  `textMinedTerms`/`unicode` → 404; PMC article page → reCAPTCHA; eutils
  `efetch db=pmc` → publisher-restriction stub; `jmdjournal.org` → 403.
* `lna2004` doi:10.1093/nar/gnh046, PMID 15047860 — closed; same metadata
  pattern as yao2006, no OA full text.

## Searches run (Europe PMC REST, `resultType=core`)

```
TITLE:"TaqMan"
TITLE:"5 nuclease"
TITLE:"TaqMan" AND (TITLE:SNP OR TITLE:genotyping) AND HAS_FT:Y
(TITLE:wheat OR rice OR maize OR soybean OR cotton) AND ABSTRACT:TaqMan AND ABSTRACT:SNP
ABSTRACT:"minor groove binder" AND (ABSTRACT:TaqMan OR ABSTRACT:probe)
ABSTRACT:TaqMan AND ABSTRACT:mismatch AND HAS_FT:Y
ABSTRACT:"delta Ct" AND ABSTRACT:probe AND ABSTRACT:mismatch        (0 hits)
TITLE:"LNA and MGB probes"
TITLE:"TaqMan-MGB"
```

Open hits were opened and read; the full-text survivors are the four landed
sources plus `maize2017` (PMC5370143), `pharm2017` (PMC5403119),
`snp_tm2017` (PMC5325530), `cyp2d7_2015` (PMC4709848), `cyp2d6_2015`
(PMC4365394) — all OA but none reporting a ΔCt, an endpoint ratio or an
MGB/non-MGB arm — and `alleletag2024` (PMC11611941), whose "100% call
rate/accuracy" numbers exist only as supplementary **images**. Restricted
XML stubs: `kutyavin2000`, `calling2018`, `snp5n2012`, `mg2016`.

## What would settle it

A source (open, or a licensed copy obtained by the user — do **not** bypass a
paywall) reporting at least one of:

1. ΔCt (or endpoint signal ratio) between a perfectly matched and a
   single-mismatched TaqMan probe, for a set of mismatch positions;
2. the same probe chemistry with and without an MGB/Tm-stabilising conjugate;
3. the same comparison at two probe lengths.

Until then, treat `TAQMAN_MISMATCH_DG` and the centrality gate as
**engineering defaults**, not as calibrated science, and do not cite the
landed wheat2019/citrus2020 rows as support for either.
