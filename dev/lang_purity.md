# Language purity — prose tics and verbosity

English, matching `dev/dag-rebuild-plan.md`; `dev/README.md` is the Russian
index and links here. This file records what `dev/prose-tics.R` measures, how
each tier was triaged, and the ordered work it leaves. The script is the
source of truth for the numbers; re-measure with:

```sh
Rscript dev/prose-tics.R            # counts + 8 samples per tier
Rscript dev/prose-tics.R --list     # every hit
Rscript dev/prose-tics.R --check    # gate the tiers that must stay empty
```

`filler`, `comment` and `verbosity` are gated at zero. Every other tier is
report-only: its hits are the repository's teaching voice or dense technical
prose, not noise. A named source does not make a sentence a tic, and a tic does
not make a sentence wrong — the point of the meter is to keep the counts
visible so a regression is noticed, not to drive a mass rewrite.

## Measured tiers

| tier | count | files | verdict |
|---|---:|---:|---|
| filler | 0 | 0 | gated; keep at zero |
| comment | 0 | 0 | gated; keep at zero |
| negation | 399 | 110 | house voice — keep, do not rewrite |
| meta | 13 | 11 | mostly regex false positives — keep report-only |
| hedge | 80 | 46 | deliberate technical emphasis — keep |
| verbosity | 0 | 0 | gated; keep at zero |
| punctuation | 104 | 8 | em-dash as term—definition separator — documented style |
| comment_verbose | 6 | 5 | long comment blocks, all justified — review only |
| prose_verbose | 4 | 2 | long audit paragraphs — split when next edited |

The counts move only when the regexes change, so a diff in this table is a
signal to re-run the audit, not a bug.

## Per-tier triage

### negation — 399

The patterns and their raw hit counts (before two fixes):

| pattern | hits | what it matches |
|---|---:|---|
| `is|are|was|were … , not` | 178 | “That policy is **data, not prose**” |
| `not … (but|rather than|instead)` | 18 | “not by haplotype, but …” |
| `not (just|only|merely|simply)` | 39 | “computed, not merely scored” |
| `not (a|an)` | 247 | **removed** — plain English, not the tell |
| `rather than|instead of|as opposed to` | 174 | “ranks rather than prints” |
| `is not a … it’s` | 1 | |
| `unlike` | 5 | “unlike a KASP amplicon” |
| `more than just`, `not because`, `far from` | 0 | |

The `not (a|an)` alternation was counting every ordinary “not a numeric
threshold” / “not a replacement”. It is not the “not X, it’s Y” rhetorical
tell, so it was dropped; that is the 541 → 399 drop. What remains is the
repository’s deliberate contrastive definition style: a sentence says exactly
what a thing is by naming what it is not. Every sampled hit carries technical
information. **Keep. Do not mass-rewrite**; a rewrite flattens the voice and
risks losing the distinction the sentence was drawn for.

### meta — 13

Raw hits were 22; the generic `the (point|key|…)` pattern caught identifiers
and cross-references (“the `key` of the row”, “the key is the sequence’s
content signature”, “the point of roadmap Y5.3”, “in practice: tune…”, “which
is why”). The pattern was narrowed to genuine signposts — `the upshot`,
`the key/main point`, `the point is that` — leaving 13, and every survivor is
still one of those technical uses. **No genuine signposting remains.** Leave
report-only; gate it if a future audit reaches zero, matching the house rule
“gate only categories at zero”.

### hedge — 80

Sampled: `actually` 41, all declared-vs-effective contrasts (“the temperature
the chemistry actually runs at”); `very` = “this very namespace”, “very large
files”; `really` = “really is monomorphic”; `simply` = “simply drops the
panel”. `clearly distinguishable`, `somewhat sticky`, `*somewhat* relaxed` are
technical requirements or deliberate nuance. **Keep.**

### punctuation — 104

All 104 are em dashes (`—`), none in an R comment; the tier never fires on
code. They separate a term from its definition in docs, vignettes and
`README.md` (`markers.csv — the per-locus report`). This is a consistent
typographic style, not a tell. **Documented style:** em dash in prose, never in
code comments. The `!!` pattern is at zero.

### comment_verbose — 6

Ordinary-`#` blocks of 15+ lines whose file banner is excluded. All six
(`marker_caps.R:51`, `marker_core.R:32`, `marker_criteria.R:53`,
`marker_criteria.R:90`, `marker_lax.R:52`, `marker_specificity.R:1289`) explain
a non-obvious constraint: source-file ordering, why a band gap is 20 bp not 50,
the criterion schema, the cache signature. This is what a comment is for.
Wide comment lines: max 117 columns, one line, and it is a URL that cannot
wrap. **Keep.** The tier is a guard against new rambling blocks, not a queue.

### prose_verbose — 4

Thresholds: a paragraph over 120 words with **no** link, `` `code` ``, call or
DOI (measured: **0** — the vignettes already cite), a paragraph over 150 words
(measured: 4), or a `##` output block over 15 lines (measured: **0**). The four
long paragraphs are audit prose — `markers-gel.Rmd` 173/240/328 and
`markers-kasp.Rmd` 380 — each lead by a bold verdict and each carrying code
references. **Split only when the paragraph is next edited**; the content is
dense, not padded, so a rewrite for length alone would lose facts.

## Actionable roadmap

Ordered, smallest first. Each item is either a guard already in place or a
deliberate non-action with its reason.

1. **Done — regex false positives.** Dropped `not (a|an)` from negation and
   narrowed the meta `key/point` pattern. These were the only counts that moved
   for a reason other than a real edit.
2. **Done — the verbosity meter.** `comment_verbose` and `prose_verbose` above;
   URL-only comment lines are not “wide”.
3. **Guard, not queue.** Keep every remaining gated tier at zero with
   `--check`; keep the four prose and six comment hits as reviewed exceptions.
4. **When the worst paragraph is touched**, split `markers-gel.Rmd:328`
   (271 words) at its existing sentence boundaries; do not shorten the facts.
5. **Comment-block review is per file, not global.** A new 15-line block is
   acceptable only if it states a constraint the code cannot; a narrative block
   should be cut.
6. **Em dash rule** (already true): prose tier may use `—`; `#` comments and
   roxygen `#'` under `R/` should not. Enforce by eye at review, not by gate.
7. **Do not** mass-edit negation, hedge or punctuation. Prior audits show they
   are the deliberate voice; a bulk diff would be churn with no reader benefit.
8. **Re-audit triggers**: a change to a regex, a new vignette, or a `NEWS.md`
   entry that adds prose — re-run the meter and update the table above if a
   count changes.
