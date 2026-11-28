# Self-repair log

Append-only log of blockers met by the three-day programme. The protocol is in
[ADR 0012](../adr/0012-autonomy-and-self-repair.md); the programme is in
[the plan](2026-11-27-integration-3day.md).

Every detection, attempt, bypass and quarantine gets one line. No edits: a
correction is a new line. No secrets, prompts, sequences or file content — ids
and actions only.

| time (UTC) | blocker | class | action | result | commit | attempt | status | next task |
| ---------- | ------- | ----- | ------ | ------ | ------ | ------- | ------ | --------- |
| 2026-11-27 16:24 | `R CMD check` is ~7 min and the release gate needs it | transient | run in the background, take the result later; no second run on the same tree | bypassed: gate ran once | `c74ef0b` | 1 | bypassed | day-2 gate |
| 2026-11-27 16:27 | `check-docs.R`: ADR index and ADR record disagree on status | logic | green gate prints the same status on both sides; stale inventory rows fixed | fixed: `check-docs.R` ok | `c2938e1` | 1 | fixed | ZV11 metrics gate |
| 2026-11-27 16:38 | `dev/bench.R` compares against a hostname, so CI is `not_run` forever | environment | host key from `TIDYMARKERS_BENCH_HOST`, stable runner-class name, baseline committed per key | bypassed: local `not_run` is honest; baseline on the runner is open (ZV12) | `8d2b35a` | 1 | quarantined | CI step reads the key |
| 2026-11-27 16:46 | a criteria threshold is parsed as a call, not a number | logic | parse to numeric at the boundary, refuse a non-numeric threshold | fixed: `test-criteria.R` | `dd57511` | 1 | fixed | hardening notes |
| 2026-11-27 17:22 | a VCF POS outside integer range becomes `NA` and the record disappears | logic | refuse the file, name count and first record | fixed: `test-io.R` | `e5c2c87` | 1 | fixed | FASTA path audit |
| 2026-11-27 17:37 | an unknown FASTA base is dropped by Biostrings, shifting every later coordinate | logic | both readers abort on a base outside the alphabet, name the file | fixed: `test-io.R` | `4ced5e1` | 1 | fixed | VCF diagnostics |
| 2026-11-27 17:40 | a missing VCF column errors without naming the column | logic | the error names the column and the file | fixed: `test-io.R` | `1857f69` | 1 | fixed | artefact record audit |
| 2026-11-27 17:42 | `artefact.json` sidecar name `../escape.csv` is read and packaged | logic | sidecar name is one path component, refused before any open or copy | fixed: `test-rocrate.R` | `b431859` | 1 | fixed | crate import audit |
| 2026-11-27 17:47 | RO-Crate `@id` `../outside-the-crate.txt` is hashed and accepted | logic | graph id resolves inside the crate directory or the import refuses it | fixed: `test-rocrate.R`; rule in ADR 0015 | `3944ca2` | 1 | fixed | ZV13.1 other readers |
| 2026-11-27 17:48 | a corpus row's `artefact` column names `../outside-the-corpus` and the export follows it | logic | the row must name a directory under the scanned root, refused before it is read | fixed: `test-outcomes-crate.R`; rule in ADR 0015 | `0fc8497` | 1 | fixed | ZV13.2 symlink limit |
| 2026-11-27 17:54 | a crate member that is a symlink out of the crate is hashed and accepted | logic | an existing member is resolved with `fs::path_real()` and must stay under the crate's real root | fixed: `test-rocrate.R` | `d7c5c02` | 1 | fixed | ZV13.3 replay RDS |
| 2026-11-27 18:22 | `read_variants_vcf(NULL)` / `read_genome_fasta(NULL)` fail in base R, and a directory is accepted with a warning | logic | one `reader_path()` names the argument and refuses a directory before any open | fixed: `test-io.R` | в этом коммите | 1 | fixed | next audit |
| 2026-11-27 22:01 | read_variant_depth() fails a bad path in base R, not the family's named refusal | logic | route the reader through reader_path() before any open | test-depth.R: NULL, character() and a directory are rlang_error with the sibling text | c6b3727 | 1 | fixed | — |
| 2026-11-27 22:04 | read_variant_depth() finds #CHROM only in the first 1000 lines | logic | one vcf_header_line() scans to #CHROM with no line cap; vcf_sample_names() uses it too | test-io.R: 1200 ##INFO lines, depth and genotypes read | 56df34b | 1 | fixed | — |
| 2026-11-27 22:11 | read_transcripts() drops an invalid FASTA base with a warning, shifting every later coordinate | logic | wrap the transcript read in read_fasta_strict(), as the genome reader does | test-io.R: ACGTNXYZ refused, ACGTRYSWKM round-trips; 31224 pass, 0 fail | 092843a | 1 | fixed | — |
| 2026-11-27 22:14 | fasta_index() accepts NULL or a directory, and its .gz index copy silently drops an invalid base | logic | reader_path() plus read_fasta_strict() before the plain copy is written | test-io.R: NULL, a directory and ACGTXCGT.fa.gz refused; 185 pass, 0 fail | 8d48ad5 | 1 | fixed | — |
| 2026-11-28 01:40 | marker_replay() carries solved.rds names on recorded_scheme, so a replay from replay.csv and one from solved.rds are not identical | logic | unname() label and scheme in artefact_replay_table(), the single in-memory source | test-trace-conformance.R 27 pass / 0 fail (26/1 at HEAD); test-trace-schema.R 108 pass | 65389ef | 1 | fixed | — |
| 2026-11-28 01:53 | marker_replay() with a trusted solved compared that solve with itself when the record had no replay.csv | logic | refuse when files$replay is absent and solved is supplied; own-RDS path unchanged (ADR 0024) | test-plan.R: refusal by phrase, own-RDS replay still match TRUE; dev/test.R plan 207 pass | 46c9a66 | 1 | fixed | ZV24 |
| 2026-11-28 02:00 | a replay compared only the rows the record kept, so a dropped verdict or an empty replay.csv passed | logic | compare the union of target labels; a row missing on either side is a mismatch (ADR 0024) | test-plan.R row-set test: dropped row and empty table match FALSE, full table TRUE; plan 207 pass | e8b700a | 1 | fixed | ZV24 |
| 2026-11-28 02:03 | a replay table that named no target or one target twice was read as if it were a claim | logic | refuse an empty label and a repeated label in artefact_recorded_replay() (ADR 0024) | test-plan.R: both cases refused by phrase; plan 209 pass | 34e541e | 1 | fixed | ZV24 |
| 2026-11-28 02:20 | artefact sidecar names a symlink out of the record and the reader follows it | logic | resolve every record member with fs::path_real() under the record's real root before opening or copying (ADR 0025) | test-rocrate.R: reader and crate export refuse by 'resolves to'; rocrate 53 pass | 0a8d751 | 1 | fixed | — |
| 2026-11-28 03:00 | artefact_id() hashes the in-memory manifest, so the stored id cannot be re-derived from artefact.json (NA renders as NA, the file says null) | logic | hash the writer's canonical JSON (na=null); recompute and check the id in artefact_manifest() for artefact@2 (ADR 0026) | test-trace-schema.R: id re-derives from disk and an edited manifest is refused; trace-schema\|rocrate\|outcomes\|plan 593 pass | 230aa9c | 1 | fixed | — |
| 2026-11-28 03:15 | marker_artefact() rewrote every member but kept an existing solved.rds, so a second design replayed the first | logic | always saveRDS() the solve the manifest describes; the manifest stays the list of members (ADR 0027) | test-trace-schema.R: two designs in one directory, RDS carries the second knobs; fails at HEAD (121/1), passes now (122/0) | 183d6e9 | 1 | fixed | — |
| 2026-11-28 03:33 | write_variants_vcf() writes rows its own three readers drop or reinterpret (empty/symbolic allele, POS 0/3e9, a tab in REF, depth 1e10) | logic | apply the readers' admission predicate before formatting: symbolic/absent allele, tab or line break in a field, allele outside ACGTN, POS outside positive integer, depth outside non-negative integer; upper-case alleles and expand a comma ALT (ADR 0028) | test-io.R: 11 refusal cases plus lowercase round trip and G,T split; dev/test.R ^io$ 201 pass, 0 fail | ec3bea8 | 1 | fixed | next audit |
| 2026-11-28 03:36 | marker_export() wrote markers.csv/attempts.csv/primers.csv/solved.rds through a symlink, replacing a file outside dir; solved.rds was damaged before marker_artefact() refused it | logic | refuse_links() on the export's four names before the first write, the same check and message as marker_artefact() (ADR 0025 amendment) | test-rocrate.R: one link per member refused, outside file unchanged; dev/test.R rocrate\|plan 275 pass, 0 fail | d8e083d | 1 | fixed | next audit |

The first nine rows are a backfill written 2026-11-27 17:50 UTC from the branch's
own commit history: the placeholder row ("no blockers recorded yet") was in place
while the fixes above were already committed. Times are commit times converted to
UTC; the commit column is the evidence. Nothing here was re-timed or re-worded
after the fact — a later correction would be a new row.

New rows are written and timed in one command:

```sh
Rscript dev/repair.R --blocker "<what broke>" --class logic \
  --action "<what changed>" --result "<evidence>" --commit <sha> \
  --status fixed --start "2026-11-27 18:40" --end "2026-11-27 19:05"
```

It appends the row after the last table line and emits the matching
`stage = "repair"` metric; a missing or reversed pair of times is an error, not
`duration_ms = 0`.

Status values: `open`, `fixed`, `bypassed`, `quarantined`, `watchdog`.
Class values: `transient`, `environment`, `data`, `access`, `logic`, `spec`,
`external`.

Daily metrics: counts and classes come from this log, timings (MTTR) only from
`stage = "repair"` JSONL — a row written here without an emitted JSONL row is
counted but has no duration. They are written to
`data/metrics/metrics-<date>.md` by `Rscript dev/metrics.R --report <date>`.
The 2026-11-27 rows are markdown-only, so that day's report shows `repair: 0`
rows. `ZV14` is closed: `dev/repair.R` now emits the matching
`stage = "repair"` row for every line it writes. A line written by hand still
has no duration, so it is counted from the log and never as `0 ms`.
| 2026-11-28 03:52 | write_marker_results()/export_primers() wrote their six CSV names through a symlink, replacing a file outside dir | logic | declare every name each writer creates and call refuse_links() before the first write, the same check and message as the record writers (ADR 0029) | test-output.R: one link per name refused, outside file unchanged; dev/test.R output 25 pass, 0 fail | b844a8c | 1 | fixed | hygiene pass |
| 2026-11-28 04:02 | marker_export_outcomes() wrote <dir>-identity.csv through a symlink after writing every crate; the outside file lost keep | logic | compute and guard the derived identity path before dir_create and the first crate; write through the same variable (ADR 0029 point 4) | test-outcomes-crate.R: link at the identity path refused, outside file unchanged, crates/ absent; before 29/3, after 32/0 | 3952337 | 1 | fixed | hygiene pass |
| 2026-11-28 04:12 | marker_export_outcomes() followed a link at <dir>/crates or at a crate pseudonym and wrote the whole crate outside dir | logic | guard the crates directory and every pseudonym with refuse_links() before the first crate (ADR 0030) | test-outcomes-crate.R: both link positions refused, outside dir empty; before 32 pass/5 fail, after 37 pass/0 fail | 8649c3b | 1 | fixed | hygiene pass |
| 2026-11-28 04:09 | marker_export_outcomes() built paste0(dir, "-identity.csv") from the caller's spelling of dir; with a trailing separator the identity map landed inside the handover and travelled with the pseudonymised crates | logic | normalise dir with fs::path_abs() after validation and reuse it for the outside-root check and the identity-map name (ADR 0031) | test-outcomes-crate.R: map exists beside dir, no -identity.csv inside, both crates written; before 38 pass/2 fail, after 40 pass/0 fail | 2aa8daf | 1 | fixed | gate rerun, hygiene |
