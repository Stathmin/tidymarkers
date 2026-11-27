# 0001 — A 16 Gb fixpur is a disk job on per-chromosome shards

- Date: 2026-11-27
- Status: accepted (measured proxy; the 16 Gb reference is not on this host)
- Closes: ROADMAP ZT T1 (ADR 0005 gate 2), ZT T2, ZX3's memory question

## Context

`ZX3` and `ZT T2` were blocked on a 16 Gb wheat reference (`fixpur`) that this
host does not have. The queue asked two different questions of the same missing
fixture:

1. **Memory:** can a host with 16 GB of RAM build the wheat index at all?
2. **Time:** what does one build and one query cost on 16 Gb?

The first question does not need the reference. It needs the peak resident set
of the path the code would take, measured as a function of chromosome size, and
a comparison against the 16 GB budget. The second question does need the
reference and stays an operator action.

## Measurements

Host: 32 cores, 29 GB RAM, 226 GB free disk, R 4.6.1. Peak is `VmHWM` from
`/proc/self/status`; reproduce with `Rscript dev/bench-index-memory.R`.

| case | input | peak RSS |
| --- | --- | --- |
| R + package baseline | no genome | 455.4 MB |
| whole TAIR10 FASTA in memory | 116 MB, 7 chromosomes | 569.7 MB |
| per-chromosome shard, TAIR10 Chr1 | 30.4 Mb, bowtie build + query | 563.2 MB |
| per-chromosome shard, synthetic | 239 MB FASTA, 250 Mb one chromosome | 885.5 MB |

The shard path holds one chromosome, not the genome. Above the 455 MB baseline
the 250 Mb run costs 430 MB, i.e. **1.7 bytes per base of peak RSS** for a
one-chromosome build plus query. For a conservative budget use **5 bytes per
base**: the largest wheat chromosome is about 0.8 Gb (16 Gb over 21
chromosomes), so

```
455 MB + 0.8e9 * 5 B/base = 4.5 GB  <  16 GB
```

A monolithic 16 Gb build is a different object: the index alone is 21.5 GB at
1.34 bytes/base, which exceeds a 16 GB RAM host before the query and the genome
are counted. Disk is not the limit: 16 GB FASTA + 21.5 GB index + ~1 GB
temporary chromosome FASTA is about 38.5 GB against 226 GB free here.

## Decision

- **Per-chromosome sharding is the supported 16 Gb path.** A 16 GB RAM host can
  build and query it; the measured proxy is 0.9 GB at 250 Mb and 4.5 GB
  conservative at 0.8 Gb.
- **A monolithic 16 Gb index is not supported on a 16 GB host.** The 21.5 GB
  index exceeds RAM; `BOWTIE_MONOLITH_LIMIT` already routes this case to shards.
- **The actual `fixpur` build-and-query run stays an operator action**, because
  no 16 Gb reference exists on this host. The blocker is data availability, not
  code or memory. When the reference is mounted, the command is
  `genome_index(marker_genome(fixpur_fa), "bowtie", dir = <cache>)` followed by
  `specificity_of(primers, index)`; record build seconds, query median and
  peak RSS in `docs/index-perf.md`.

## Consequence for the queue

`ZT T1` and `ZT T2` close on the measured shard path and the TAIR10 gold set. A
future 16 Gb run can reopen them only if it contradicts the extrapolation; the
memory budget is no longer an open question.
