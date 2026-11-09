# Language purity — prose tics and verbosity

English, matching `dev/dag-rebuild-plan.md`; `dev/README.md` is the Russian
index and links here. This file records what `dev/prose-tics.R` measures, how
each tier was triaged in 2026-11, and the ordered work it leaves. The script is
the source of truth for the numbers; re-measure with:

```sh
Rscript dev/prose-tics.R            # counts + 8 samples per tier
Rscript dev/prose-tics.R --list     # every hit
Rscript dev/prose-tics.R --check    # gate the tiers that must stay empty
```

`filler`, `comment`, `verbosity`, `prose_verbose` and `encoding` are gated at
zero. The other tiers are report-only: each hit was read in 2026-11 and either
edited or kept because it carries technical information. A named source does
not make a sentence a tic, and a tic does not make a sentence wrong; the meter
keeps the counts visible so a regression is noticed.

## Measured tiers

| tier            | count | files | verdict                                               |
| --------------- | ----: | ----: | ----------------------------------------------------- |
| filler          |     0 |     0 | gated; keep at zero                                   |
| comment         |     0 |     0 | gated; keep at zero                                   |
| negation        |   397 |   110 | reviewed; kept — each hit states a technical contrast |
| meta            |     0 |     0 | narrowed to genuine signposts, none survive           |
| hedge           |     4 |     4 | reviewed; four technical qualifiers kept              |
| verbosity       |     0 |     0 | gated; keep at zero                                   |
| punctuation     |     0 |     0 | markup excluded + prose edited; `!!` never fired      |
| encoding        |     0 |     0 | gated; non-Unicode signal scanner                     |
| comment_verbose |     6 |     5 | reviewed; all six are constraint or reference blocks  |
| prose_verbose   |     0 |     0 | gated; long audit paragraphs split                    |

The counts move only when the regexes change or a file is edited, so a diff in
this table is a signal to re-run the audit, not a bug.

## Per-tier triage (2026-11)

### negation — 397

The patterns are the "not X, it's Y" family in every shape: `is, not`, `not
but`, `not just`, `rather than`, `instead of`, `unlike`, `far from`. Every
sampled hit carries technical information: "data, not prose", "the gate is not
Tm", "computed, not copied". This is the repository's contrastive definition
style, so it is reported and never gated. **Kept after review.** Do not
mass-rewrite: a bulk diff flattens the voice and can lose the distinction the
sentence draws. Edit one only when a specific sentence is genuinely poor.

### meta — 0 (was 13, raw 22)

The regex now matches only genuine signposts: `the upshot|takeaway|bottom
line|lesson`, `the key/main point`, `the point is that`, `here's the thing`,
`at its heart`, `serves as`, and the `it is important to` family. Ordinary
connectives (`in practice`, `precisely because`, `which is why`, `the reason
is`) are technical English and are excluded. No hit survives the narrowed
patterns, so the tier is zero. Gate it if it stays zero through the next
audit.

### hedge — 4 (was 80)

The regex keeps only genuine softeners: `simply|basically|obviously|clearly|
trivially|essentially|literally`, `quite|somewhat|pretty|extremely|truly`,
`seems to|appears to|arguably|perhaps|maybe|presumably`, `of course`, `in fact`.
Technical qualifiers that distinguish a nominal from an effective value
(`actually`, `very`, `really`, `particularly`, `especially`, `highly`) are
excluded. The four survivors are legitimate:

- `vignettes/markers-gel.Rmd:26` "clearly distinguishable" — a requirement.
- `R/marker_decision_render.R:299` "not to be pretty" — a constraint.
- `R/marker_dpcr.R:43` "somewhat sticky" — measured probe behaviour.
- `R/marker_select.R:31` "_somewhat_ relaxed" — deliberate nuance.

### punctuation — 0 (was 104)

The two patterns are em dash (`—`) and `!!`. In 2026-11 the scanner learned to
skip markup lines (table rows `|`, list items `- * +`, numbered items `1.`),
which removed the false positives where a dash separated a term from its
definition. The prose that remained was edited case by case: every em dash in
`README.md` and the six vignettes became a colon, comma, semicolon or pair of
parentheses. `!!` never fired.

### encoding — 0 (gated)

`encoding_hits()` scans `DESCRIPTION`, `README.md`, `man/*.Rd`,
`vignettes/*.Rmd`, `R/*.R`, `tests/` and `inst/{csv,txt,md,yml}` for invalid
UTF-8, U+FFFD (replacement character), mojibake (`Ã`, `Â`, `â€`), zero-width
characters (U+200B/C/D, U+FEFF) and ASCII control characters. It reports zero
and is gated, so a future non-Unicode paste fails `--check`.

### comment_verbose — 6

Ordinary-`#` blocks of 15+ lines whose file banner is excluded. All six are
`marker_caps.R:51`, `marker_core.R:32`, `marker_criteria.R:53`,
`marker_criteria.R:90`, `marker_lax.R:52`, `marker_specificity.R:1289`. Each
explains a non-obvious constraint (source-file ordering, why a band gap is 20
bp not 50, the criterion schema, the cache signature). **Reviewed 2026-11: all
six stay.** The tier is a guard against new rambling blocks, not a queue.

### prose_verbose — 0

A paragraph over 120 words with no link, `` `code` `` or call, a paragraph over
150 words, or a `##` output block over 15 lines. The four long audit
paragraphs (`markers-gel.Rmd` 173/240/328, `markers-kasp.Rmd` 380) were split
at their existing verdict and scheme boundaries, facts intact. Gated.

## Actionable roadmap

1. **Done — regex false positives.** `not (a|an)` was dropped from negation;
   meta was narrowed to genuine signposts (now 0); hedge dropped technical
   qualifiers (80 → 4).
2. **Done — the verbosity meter.** `comment_verbose` and `prose_verbose`; a
   URL-only comment line is not "wide".
3. **Done — punctuation.** Markup lines are skipped, so only running prose is
   measured; the remaining prose dashes were edited.
4. **Done — encoding gate.** A non-Unicode signal fails `--check`.
5. **Guard, not queue.** Keep `filler`, `comment`, `verbosity`, `prose_verbose`
   and `encoding` at zero with `--check`; keep the six comment hits as reviewed
   exceptions.
6. **Do not** bulk-rewrite negation. It is the deliberate contrastive voice and
   each hit carries technical information; edit one sentence at a time only
   when it reads poorly.
7. **Em dash rule.** Prose may use `—` only when a colon or parenthesis would
   not read better; `#` comments and roxygen `#'` under `R/` should not.
8. **Re-audit triggers.** A change to a regex, a new vignette, or a `NEWS.md`
   entry that adds prose: re-run the meter and update the table above.
