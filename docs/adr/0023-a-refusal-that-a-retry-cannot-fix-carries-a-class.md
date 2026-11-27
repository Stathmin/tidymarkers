# 0023 — A refusal that a retry cannot fix carries a class

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0022](0022-a-fasta-base-outside-the-alphabet-is-refused-by-every-reader.md)
  (a bad FASTA base is refused by name)

## Context

`ath_genome()` treats an unreadable cache as a transient problem: the file may
be truncated or half-written, so it deletes the cache, downloads once and reads
again. That is the right response to a broken transfer.

ADR 0022 then made a cached FASTA with a letter outside the DNA/IUPAC alphabet
a named refusal. The refusal arrives through the same channel as a truncated
gzip — an error from the read — so `ath_genome()` would delete a good download,
fetch the same letters again and report the failure as a read or download
problem. Measured before the fix: a cache containing `ACGTXCGT` produced the
download path, not the alphabet sentence.

The distinction was previously made by matching the message text. That is
fragile: `abort_msg()` runs user-facing messages through the translation layer,
so the sentence is not a stable key.

## Decision

When a refusal cannot be fixed by retrying or refetching, the condition carries
a stable class, and a caller that retries branches on the class, not on the
message.

- `read_fasta_strict()` raises its alphabet refusal with class
  `"tidymarkers_fasta_alphabet"` in addition to `rlang_error`.
- `ath_genome()` re-signals that condition unchanged: it does not delete the
  cache and does not call `tair10_download()`.
- Any other read failure keeps the old behaviour: delete once, refetch once,
  then report with the underlying condition as `parent`.
- Tests use the class indirectly through `expect_error(..., "sentence")` and
  assert that the mocked download is never called.

## Consequences

- A corrupt download and a wrong alphabet now have different outcomes: the
  first is repaired by refetching, the second names the file.
- The rule is about retry semantics: a class is added only where a caller must
  distinguish "retry may help" from "retry cannot help", not to every error.
- Message text stays user-facing and translatable; no code depends on it.
- Test: `test-ath-genome-cache.R` — "a cache with a base outside the alphabet
  is refused, not refetched".
