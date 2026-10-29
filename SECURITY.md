# Security policy

## Supported versions

`tidymarkers` is experimental. Only the latest release and the current `main`
receive security fixes.

| Version | Supported |
|---|---|
| 0.4.x   | yes       |
| < 0.4   | no        |

## Reporting a vulnerability

Do **not** open a public issue. Email <uldas1508@gmail.com> with:

- what is affected (function, assay scheme, or dependency);
- how to reproduce it (a minimal reprex, see below);
- the impact you think it has.

You will get an acknowledgement within 7 days and an assessment within 30 days.
Please give a reasonable window for a fix before public disclosure.

## What counts as a security issue here

The package is pure R and reads local files, but these are in scope:

- **Unsafe file handling** — a path from a VCF, FASTA header, or `marker_genome()`
  argument that escapes an intended directory or writes where you did not ask.
- **Command injection** — `system()`/`system2()` calls to BLAST+,
  OligoArrayAux, or `Rsubread` where an argument is interpolated from input.
- **Dependency compromise** — a malicious or typosquatted package in
  `renv.lock` or `DESCRIPTION`.
- **Data exfiltration** — any code path that sends sequence or variant data off
  the machine. The package is offline by design; there is no telemetry.

Not in scope: a wrong primer, a thermodynamic model disagreeing with a bench
result, or a source that turns out to be misquoted. Those are scientific bugs —
please file them with the bug template and the evidence.

## A minimal reprex

Never include patient data, unpublished sequences, or credentials. Reduce the
input to a small synthetic genome:

```r
tidymarkers::marker_genome(...)  # or a short FASTA you can share
```

Run it, paste the code and the exact error.
