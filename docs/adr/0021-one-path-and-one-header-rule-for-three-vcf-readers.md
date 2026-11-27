# 0021 — One path and one header rule for three VCF readers

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)
  (one allele and `POS` rule for all three readers),
  [ADR 0016](0016-one-vcf-admission-rule.md) (one admission rule),
  [ADR 0019](0019-an-absent-allele-is-not-a-base.md) (an absent allele is
  absent before it is normalised)

## Context

ADR 0018 put the three VCF readers — `read_variants_vcf()`,
`read_variant_genotypes()` and `read_variant_depth()` — under one allele and
`POS` rule. The record-level rule was applied to all three. Two rules that are
not about a record were not:

1. **The path.** `reader_path()` names the argument, refuses `NULL`, an empty
   vector and a directory before any file is opened. The two record readers use
   it; `read_variant_depth()` kept a bare `fs::file_exists()`:
   `read_variant_depth(NULL)` failed in `enc2utf8()` and a directory failed in
   `file()` with a warning, while its siblings returned `rlang_error`
   (`"`path` must be one VCF file path."`, `"names a directory"`).

2. **The header.** `read_variant_depth()` scanned at most 1000 lines for the
   `#CHROM` line (`readr::read_lines(path, n_max = 1000L)`), then looked for a
   line starting with `#CHROM`. A legal GATK-style header carries more than a
   thousand `##INFO`/`##FORMAT` lines; on such a file the depth reader aborted
   with `"has no '#CHROM' header line"` while `read_variants_vcf()` read the
   same file. The bound was an implementation detail leaking into the contract:
   a header is terminated by `#CHROM`, not by a line count.

Both defects are silent in the happy path and only show on an input the reader
claims to accept. The same class was found twice before (ZV17, ZV20): a rule
that holds for the family is written in one member and not the others.

## Decision

A rule shared by the VCF family lives in one helper, and every reader uses it.

- **Path:** every reader starts with `reader_path(path, "VCF")`. A reader does
  not call `fs::file_exists()` itself.
- **Header:** `vcf_header_line(path)` scans line by line until `#CHROM` and
  refuses a file whose first non-`##` line is not `#CHROM`. There is no line
  cap. `vcf_sample_names()` and `read_variant_depth()` both take the line from
  it, so sample columns and record columns come from one header.

The normalisation of `chrom`/`pos`/`ref`/`alt` stays where ADR 0018 put it; this
decision only covers the read before the first record.

## Consequences

- `read_variant_depth()`, `read_variants_vcf()` and `read_variant_genotypes()`
  now fail the same way on the same bad path and read the same legal header.
- The header scan costs one line per `##`, no memory for the whole header, so a
  large VCF is still streamed.
- A VCF whose first non-`##` line is not `#CHROM` is refused by name instead of
  looking like an empty file.
- Tests: `test-depth.R` ("the depth reader names a path that is not one existing
  file") and `test-io.R` ("a VCF header past 1000 lines is still found", 1200
  `##INFO` lines) fail before this change on both counts.
