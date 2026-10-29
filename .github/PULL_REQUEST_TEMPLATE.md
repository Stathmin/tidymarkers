<!--
Keep this short. The merge gate and the four kinds of change are in
CONTRIBUTING.md; this template only asks what a reviewer cannot infer.
-->

## What and why

<!-- One paragraph: the change, and the problem it solves. -->

**Kind of change:** <!-- rule / fixture / wiring / prose -->

## Evidence

<!--
For a rule: the source, the values it reports, and the data-raw script that
re-derives them. For the rest: the test or call site that proves it.
-->

## Checklist

- [ ] `Rscript dev/check.R` passes (or I name the pre-existing failure I did not touch)
- [ ] New behaviour has a test in `tests/testthat/`
- [ ] `styler` and `lintr` are clean on the files I changed
- [ ] `roxygen` is regenerated if a doc comment changed
- [ ] `NEWS.md` has an entry for user-visible change
- [ ] `ROADMAP.md` is updated if this closes or changes a queued item
- [ ] No fabricated source, number, or quotation; gaps are labelled as gaps
- [ ] No secrets, private sequences, or patient data in the diff

## How I verified it

<!--
The exact commands you ran and what they printed. One test file at a time,
not the full suite: `Rscript dev/test.R '<stem>'`. If this is docs-only, say so.
-->

```sh
# commands
```

<!-- Output summary, or "docs only". -->

## Risks and follow-ups

<!--
What could break for a caller? Migration needed? Attributes that should become
first-class fields? Anything you deliberately left out of scope.
-->
