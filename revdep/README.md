# Reverse dependencies

This directory is where `revdepcheck::revdep_check()` writes its report at
release time. It is currently empty on purpose: `tidymarkers` has no CRAN
release yet, so there are no reverse dependencies to check.

CRAN asks for the reverse-dependency check only once the package is on CRAN and
a release would break a downstream package. Before the first release this file
is the statement that the check does not apply, not a missing report.

To regenerate the report after the first CRAN release:

```r
revdepcheck::revdep_check(num_workers = 4)
```

The report lands here and is summarised in `cran-comments.md`.
