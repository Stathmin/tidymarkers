# ZW3.4 anchor, pass 4 — still no external measured oligo-Tm anchor with printed buffer

Literature pass, 2026-11-03. Passes 1–3 are
`zw3.4-anchor.md` and `zw3.4-anchor-pass3.md`. Pass 4 searched **new** material
only: the Europe PMC full text of Lobato & O'Sullivan 2018
(`PMC7112910`, `10.1016/j.trac.2017.10.015`) and a fresh set of term greps over
the local full-text scratch and the Zotero lane. Nothing already read was
re-read.

## Answer up front

**Not found.** No source printing an oligo Tm **and** its own numeric buffer/salt
was located that would satisfy the four anchor conditions
(measured-not-computed; numeric `[Na+]`/`[Mg2+]`; 18–25 nt; a residual that is a
slope not a level). The only true anchor candidate remains mm2020
(`PMC8163305`, 19-mers, 50 mM NaCl, Mg 0), which fails length and level-vs-slope
exactly as pass 3 records. Lobato 2018 states no measured oligo Tm at all — it is
an RPA *methods* review with no duplex-Tm table and no numeric salt.

## What was searched this pass (new, not a re-sweep)

- Lobato & O'Sullivan 2018 full text (`PMC7112910`): read in full; RPA primer
  length / amplicon / probe geometry numbers only, no measured Tm, no `[Na+]` /
  `[Mg2+]`.
- Local scratch `/tmp/tidymarkers-scratch/lit/` (38 PMC XML): greps for salt and
  Tm co-occurrence (`mM NaCl`, `mM Mg`, `melting temperature`, `Tm`) and for the
  remaining scheme constants. No measured-oligo-Tm-plus-buffer hit.
- Zotero lane greps for `30-35`, `5 °C`, `band gap`, `read distance`,
  `homopolymer`, `size delta`: no new Tm-anchor candidate.

## Verdict — unchanged

The −16.86 °C (Schildkraut) / −9.33 °C (Owczarzy) bias is still not re-anchorable
from the searched literature, for the structural reason already stated: applied
design papers use a tool and a proprietary master mix and print neither a
measured Tm nor a numeric buffer; the one measured set (mm2020) is 19-mers only
and its residual is a level. No window change is proposed; `ZW3.4` stays open.
