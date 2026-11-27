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

### TAIR10 gold set and the per-locus budget (ADR 0005 gate 2)

Three real TAIR10 SNP targets (`marker_example("SNP")` rows: Chr1:18,519,408,
Chr3:6,475,771, Chr3:3,229,835), the KASP designer at `flank = 150`, R 4.6.1.
The budget is the measured ratio, not a wall-clock number from another host:
the index path must stay within 2x of the in-memory path per locus.

| path | per-locus (s) | median (s) | max / median |
| --- | --- | --- | --- |
| memory (`DNAStringSet`) | 20.4, 18.6, 18.5 | 18.64 | 1.09 |
| bowtie (Chr1/Chr3 shards) | 21.7, 18.3, 17.2 | 18.32 | 1.18 |

Index / memory median is **0.98x**; every locus is inside 2x. The golden index
set itself is `tests/testthat/fixtures/tair10/` (three 40 kb TAIR10 chromosomes
and `expected.csv`): memory and bowtie return the same hit count and the same
verdict for a unique primer on each chromosome and for a primer that occurs
twice.

### Peak memory of the two paths (16 Gb decision)

`Rscript dev/bench-index-memory.R` reports `VmHWM` on this host (29 GB RAM):

| case | peak RSS |
| --- | --- |
| R + package baseline | 455.4 MB |
| whole TAIR10 FASTA in memory (116 MB) | 569.7 MB |
| one TAIR10 chromosome shard (30.4 Mb) | 563.2 MB |
| one synthetic chromosome shard (250 Mb, 239 MB FASTA) | 885.5 MB |

The shard path holds one chromosome, not the genome: 1.7 bytes per base of peak
above the baseline, 5 bytes/base conservative. The largest wheat chromosome is
about 0.8 Gb, so the measured shard path needs about 4.5 GB, which fits a 16 GB
host; a monolithic 16 Gb index (21.5 GB) does not. The decision and the disk
arithmetic are in `docs/decisions/0001-16gb-fixpur.md`. The actual 16 Gb run is
an operator action because no such reference is on this host.

Not measured here: a 16 Gb wheat build (no such reference on this host), a
pan-genome A/B/D run.
