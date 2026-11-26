# Index performance on one locus

One locus, two strategies, same verdict. Reproduce with `Rscript dev/bench-index.R`
(R 4.6.1, 2026-11-26). The numbers are a single host and a synthetic genome, not
TAIR10: the point is the order of magnitude and what is *not* in memory, not a
wheat-scale claim.

Fixture: 3,000,000 bp in three 1 Mb chromosomes, one planted 24-mer locus and its
reverse complement (`FWD`/`REV`), FASTA 2.9 MB. The memory strategy scans the
3'-seed (`k = 15`, `max_mismatch = 0`); the bowtie strategy scans the whole
primer on per-chromosome shards (`max_mismatch = 0`). Both find one hit per
primer, i.e. the same verdict (`unique`).

| strategy | shards | build (s) | query median (s) | hits | index object (MB) |
| --- | --- | --- | --- | --- | --- |
| memory (3'-seed) | 1 | - | 0.0800 | 1 | 2.9 |
| bowtie, per-chromosome (whole-primer) | 3 | 1.576 | 0.1180 | 1 | 0.0 |

- **The index object holds paths, not sequence.** `object.size()` of the
  `genome_index` is under 0.05 MB for three shards; the sequence lives in 27.1 MB
  of bowtie files on disk. A query opens a shard, maps the primers and closes it;
  nothing is cached beyond the per-genome chromosome cache, which a shard build
  deliberately does not populate.
- **The build is one-off and can be larger than the first query.** 1.58 s to
  build three shards against 0.08 s for the in-memory scan of 3 Mb. On TAIR10 the
  measured ratio stays in this shape (index about 1.34 bytes/base, build once,
  query flat in primer count); the index wins when the genome no longer fits in
  RAM or when the primer count grows.
- **Per-chromosome shards buy peak memory, not query speed.** Three shards are
  three bowtie invocations, so the query pays a little more than one monolithic
  index. That is the trade the 21.5 GB monolith limit makes: a wheat or A/B/D
  pan-genome run needs shards that fit the host.
- **The two counts answer different questions.** `scan_question()` reports
  `"3'-seed"` for memory and `"whole-primer"` for bowtie; the verification
  compares the verdict and the failing criterion ids, not raw counts.

Not measured here: a 16 Gb wheat build (no such reference on this host), a
pan-genome A/B/D run, and the 2x per-locus budget in ADR 0005 gate 2 (that gate
needs the TAIR10 gold set).
