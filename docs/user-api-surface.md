# Пользовательская поверхность API и граф соединимости

Дата срез: **2026-11-27**. Пакет `tidymarkers` 0.4.0, ветка
`docs/user-call-graph-20261127`.

Это инвентарь того, что пользователь может вызвать, и карта того, как выход
одной функции стыкуется со входом другой **по классу объекта**, а не по
названию. Сценарии «J1–J7» здесь намеренно не строятся: сначала нужен честный
перечень интерфейсов и дыр в нём.

Источники: `NAMESPACE`, формы функций в `R/`, roxygen-блоки и `man/`.
Машинная версия — [`user-api-surface.csv`](user-api-surface.csv), граф —
[`user-api-graph.dot`](user-api-graph.dot). Проверка —
`Rscript dev/check-user-api.R`.

Отрисовка графа (209 узлов, 949 рёбер): `dot -Tsvg docs/user-api-graph.dot
-o user-api-graph.svg`; для плотной раскладки читаемее `sfdp -Tsvg
docs/user-api-graph.dot -o user-api-graph.svg`. Файл `.dot` — производный,
руками не правится: единственный генератор — `Rscript dev/check-user-api.R
--write`.

## 0. Что посчитано

| Что | Число |
|---|---|
| Экспортов в `NAMESPACE` | **208** |
| S3-методов (включая `+.marker_plan`) | **52** |
| Всего записей в CSV | **260** |
| Пользовательских (`user_facing = TRUE`) | **138** |
| Служебных (`user_facing = FALSE`) | **122** |
| Типизированных рёбер в графе | **949** |
| Точек входа (нет типизированных входов) | **58** |
| Листьев (нет типизированных потребителей) | **116** |
| Изолятов (ни входа, ни выхода по классам) | **29** |

Типизированное ребро `A → B` означает: `output_class(A)` совпадает с одним из
`inputs(B)`; `assay` считается подклассом `marker`. Всё, что имеет класс
`character`, `numeric`, `logical`, `tibble`, `path`, `any` или `none`, рёбер не
даёт — это и есть «дыры», перечисленные в §8.

## 1. Инвентарь

Полные аргументы — в CSV; таблица ниже обрезает длинные списки после восьми
позиций. Колонка «Польз.» означает «предназначена пользователю», а не
«экспортирована»: экспортированы все 208.

### 1.1 Ingest — чтение и приведение данных

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `alignment_from()` | — | seqs | — | `sequences` → `alignment` | Multiple alignment from a set of sequences (for rprimer) |
| `annotate_variants()` | да | variants | txdb, seqsource, gene | `variants + txdb` → `variants` | Annotate variants with a TxDb transcript model |
| `ath_genome()` | да | — | cache, release | `—` → `genome` | Download (once) and load the TAIR10 genome as a DNAStringSet |
| `ensembl_gene()` | да | gene_id | species, cache | `character` → `gene` | Fetch a gene record from Ensembl Plants REST (cached on disk) |
| `ensembl_variants()` | да | chrom, start, end | species, cache | `region` → `variants` | Natural variants overlapping a region, from Ensembl Plants (1001 Genomes) |
| `fasta_index()` | — | path | — | `path` → `path` | Ensure a FASTA has a .fai index, decompressing .gz once if needed |
| `genome_seq()` | да | genome, chrom | start, end | `genome + region` → `sequence` | Read a region as a DNAString, from an indexed source or a DNAStringSet |
| `haplotypes_from_alignment()` | да | seqs | reference, chrom, freq, maf | `sequences` → `haplotypes` | Align transcripts and call haplotypes against a reference |
| `marker_as_granges()` | — | x | extra | `variants` → `granges` | Variants tibble <-> GRanges |
| `marker_bsgenome()` | да | — | package | `—` → `genome` | Default Arabidopsis reference genome as a BSgenome |
| `marker_ensdb()` | да | — | gtf, cache, release | `path` → `txdb` | Build (or load) an EnsDb gene model from the Ensembl Plants GTF |
| `marker_example()` | да | — | type | `character` → `variants` | Bundled example target variants (one per marker type) |
| `marker_features()` | да | path | chrom, feature, chrom_map | `path` → `features` | Read a genome feature annotation for the plot feature track |
| `marker_from_granges()` | — | x | — | `granges` → `variants` | GRanges → variants (пара marker_as_granges) |
| `marker_genome()` | да | — | path, chrom_map, blast_db, bsgenome, fm_index, bowtie_index | `path` → `genome` | Indexed genome source backed by an Rsamtools FASTA index |
| `marker_ingest_approve()` | — | queue, record_id | registry, reviewer | `queue + character` → `registry` | Approve one queued record into the registry |
| `marker_ingest_literature()` | — | path | min_evidence, min_confidence, provenance_required | `path` → `queue` | Ingest an external literature corpus into a review queue |
| `marker_ingest_registry()` | — | — | registry | `—` → `registry` | The approved literature records |
| `marker_txdb()` | да | — | package | `—` → `txdb` | Default Arabidopsis transcript annotation (TxDb) |
| `read_genbank()` | да | path | region | `path` → `features` | Read a GenBank flat file through a cached feature index |
| `read_genome_fasta()` | да | path | chrom_map, seqids, use_index | `path` → `genome` | Read an arbitrary genome FASTA, optionally normalising sequence names |
| `read_transcripts()` | да | path | merge_accessions | `path` → `haplotypes` | Read a transcript FASTA into a haplotype set |
| `read_variant_depth()` | да | path | sample | `path` → `depth` | Read per-sample depth (DP/AD/GT) from a VCF |
| `read_variants_vcf()` | да | path | chrom, start, end, backend | `path` → `variants` | Read variants from a VCF/BCF file (streaming, memory-mapped) |
| `snp_template()` | — | genome, chrom, pos | flank | `genome + region` → `sequence` | Extract a DNA window centred on a SNP as a DNAString |
| `variant_confidence()` | — | depth | min_depth, min_weight | `depth` → `numeric` | Turn DP/AD/GT into a confidence weight |
| `variant_type()` | да | ref, alt | — | `character + character` → `character` | Classify a variant by its REF/ALT strings |

### 1.2 Design — поиск и дизайн кандидатов

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `allele_templates()` | — | ref_tmpl, snp_i, ref, alt | — | `sequence + variant` → `sequences` | Plus-strand REF and ALT allele templates for an indel |
| `design_arms_marker()` | да | genome, snp | flank, min_band_gap | `genome + variant` → `marker` | Design a tetra-primer ARMS-PCR marker at a variant |
| `design_as_primers()` | — | window, v_start, ref, alt | strands, lengths, stab_ats, temp, P, ions, linked, gc_range … +4 | `sequence + numeric` → `primers` | Grid-search allele-specific primers maximising Taq-based discrimination |
| `design_as_primers_indel()` | — | ref_tmpl, snp_i, ref, alt | strands, lengths, temp, P, ions, gc_range, tm_range, BPPARAM … +1 | `sequences` → `primers` | Allele-specific primer candidates for an insertion/deletion (junction scan) |
| `design_asq_marker()` | да | genome, snp | variants, haplotypes, flank, BPPARAM | `genome + variant + limits` → `assay` | Design an ASQ / PACE marker (variant at the penultimate 3' base) |
| `design_auto_marker()` | да | genome, snp | strategy, schemes, context | `genome + variant + context` → `auto_result` | Design a marker with the first chemistry that works |
| `design_caps_marker()` | да | genome, snp | flank, enzymes, min_gap, haplotypes | `genome + variant + enzyme_table` → `marker` | Design a CAPS / dCAPS marker at a variant |
| `design_castpcr_marker()` | да | genome, snp | flank, as_lengths, mate_lengths, probe_lengths, blocker_lengths, probe_margin, probe_tm_max, min_product … +3 | `genome + variant` → `marker` | Design a castPCR / allele-specific qPCR marker at a variant |
| `design_common_primers()` | — | window, snp_i | strand, lengths, gaps, target_tm, min_gc, max_gc, tm_range, n … +1 | `sequence` → `primers` | Design candidate common primers |
| `design_crispr_marker()` | да | genome, snp | flank, cas, strand, spacer_lengths, gc_range, specificity_max, preamp, min_product … +3 | `genome + variant` → `marker` | Design a CRISPR-Dx (SHERLOCK / DETECTR) marker at a variant |
| `design_dpcr_marker()` | да | genome, snp | flank, min_product, max_product, lengths, tm_range, gc_range, probe_lengths, probe_gc … +2 | `genome + variant` → `marker` | Design a dPCR hydrolysis-probe marker at a variant |
| `design_duplex()` | да | genome, loci | readout, homologous | `genome + variants` → `marker` | Design a duplex assay: two loci in one tube |
| `design_hrm_marker()` | да | genome, snp | flank, min_product, max_product, lengths, tm_range, gc_range | `genome + variant` → `marker` | Design an HRM marker at a variant |
| `design_indel_size_marker()` | да | genome, snp | flank, min_gap, min_product, max_product, lengths, tm_range, gc_range, engines | `genome + variant` → `marker` | Design an InDel size-polymorphism marker |
| `as_marker_solved()` | да | report | — | `report + panel + solved` → `marker_solved` | Continue a batch report or a panel as the typed solved object |
| `design_kasp_markers()` | да | genome, loci | variants, BPPARAM | `genome + variants` → `report` | Design KASP markers for a table of loci, with a per-locus verdict |
| `design_lamp_marker()` | да | genome, snp | flank, arm_lengths, outer_lengths, outer_tm, inner_tm, gc_range, min_product, max_product … +2 | `genome + variant` → `marker` | Design a LAMP marker at a variant |
| `design_lax_marker()` | да | genome, snp | haplotypes | `genome + variant + haplotypes` → `marker` | Design a lax KASP: one ALT primer per haplotype |
| `design_marker()` | да | genome, snp | — | `genome + variant` → `auto_result` | Устаревший алиас design_auto_marker |
| `design_marker_batch()` | да | genome, loci | design, BPPARAM | `genome + variants` → `report` | Run a per-locus design (or any single-locus call) over many loci, in parallel |
| `design_markers()` | да | genome, loci | scheme, context | `genome + variants` → `report` | Design markers for a table of loci with one chemistry (or the rescue ladder) |
| `design_rpa_marker()` | да | genome, snp | flank, primer_len, min_product, max_product, design_probe, probe_len, gc_range, tm_range … +2 | `genome + variant` → `marker` | Design an RPA (recombinase polymerase amplification) marker at a variant |
| `design_sanger_marker()` | да | genome, snp | flank, min_dist, read_length, min_product, max_product, lengths, tm_range, gc_range … +3 | `genome + variant` → `marker` | Design a Sanger sequencing marker at a variant |
| `design_snapshot_marker()` | да | genome, snp | flank, amp_lengths, min_amplicon, max_amplicon, amp_tm, amp_gc, ext_lengths, ext_tm … +3 | `genome + variant` → `marker` | Design a SNaPshot / single-base-extension marker at a variant |
| `design_ssr_marker()` | да | genome, snp | flank, min_product, max_product, lengths, units, min_copies, tm_range, gc_range | `genome + variant` → `marker` | Design an SSR / microsatellite marker at a repeat locus |
| `design_taqman_marker()` | да | genome, snp | flank, min_product, max_product, primer_lengths, probe_lengths, primer_tm, probe_tm, probe_modification … +3 | `genome + variant` → `marker` | Design a TaqMan 5'-nuclease allelic-discrimination marker at a variant |
| `find_caps_enzymes()` | да | genome, snp | flank, enzymes, min_gap | `genome + variant + enzyme_table` → `enzyme_table` | Enzymes whose digest distinguishes the REF and ALT amplicons (CAPS scan) |
| `find_dcaps_enzymes()` | да | genome, snp | flank, enzymes, min_gap | `genome + variant + enzyme_table` → `enzyme_table` | dCAPS: an enzyme the primer can turn on for one allele |
| `haplotype_cover()` | да | haplotypes, snp, design | window | `haplotypes + variant` → `primers` | Design one ALT primer per primer-window group, merging by footprint |
| `haplotype_groups()` | да | haplotypes, snp | window | `haplotypes + variant` → `groups` | Group haplotypes by the alleles within a window around the target |
| `marker_attach_haplotypes()` | да | marker, haplotypes | snp, min_alt_freq | `marker + haplotypes` → `marker` | Attach known population haplotypes to any designed marker |
| `marker_context()` | да | — | variants, haplotypes, depth, min_depth, min_weight | `variants + haplotypes + depth` → `context` | Add variant / haplotype context to a plan |
| `marker_design_panel()` | да | genome, chemistries | variants, haplotypes, targets, region, ld, costs, created, dir | `genome + variants + haplotypes + targets + costs + character` → `panel` | Design a marker panel from variants, haplotypes and chemistries |
| `marker_run_context()` | да | — | genome, variants, population, haplotypes, conditions, engines, options | `genome + variants + haplotypes` → `context` | A marker run context |
| `marker_targets()` | да | — | data, chrom, pos, ref, alt, alt2, type, label … +1 | `variants` → `targets` | Add target loci to a plan |
| `pick_target_snp()` | да | variants, genome, gene | variant_id, clean_3p, specificity_candidates, BPPARAM, scheme, max_tries | `variants + genome + gene` → `variant` | Choose a designable target variant from a variant table |
| `rprimer_assays()` | — | alignment | check | `alignment` → `tibble` | Independent oligo design + match report on an alignment (rprimer) |
| `screen_snp_targets()` | да | variants, genome | lengths, strands, clean_3p, hotzone, linked, limits | `variants + genome` → `screen` | Cheap, vectorised screen of candidate target SNPs |

### 1.3 Specificity — проверка геномной уникальности

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `assay_amplicons()` | да | assay, genome | max_mismatch, BPPARAM, products | `assay + genome` → `hits` | Products predicted for each allele-specific primer against the common one |
| `blast_primers()` | да | primers, db_path | — | `primers + index` → `hits` | BLAST every primer against a local genome database (rBLAST, optional) |
| `build_blast_db()` | — | genome | path | `genome` → `index` | Build (once) a BLAST nucleotide database from a genome |
| `exact_hit_table()` | да | sequences, genome | max_mismatch, BPPARAM, scheme, seed_len | `sequences + genome` → `hits` | Exact genome hits with coordinates |
| `genome_index()` | — | genome | strategy, chroms, dir, threads | `genome` → `genome` | Attach a declared genome-access strategy to a genome |
| `in_silico_pcr()` | да | f_primer, r_primer, genome | max_mismatch, BPPARAM, scheme, pam | `primers + genome` → `hits` | In-silico PCR: every product a primer pair would amplify |
| `marker_blast_index()` | — | genome | path | `genome` → `index` | Build (once) a BLAST database for a genome |
| `marker_bowtie_index()` | — | genome | basename, threads, chroms | `genome` → `index` | Build a bowtie1 index for a genome (Rbowtie::bowtie_build) |
| `marker_fm_index()` | — | genome | basename | `genome` → `index` | Build an FM-index for a genome (Rsubread::buildindex) |
| `offtarget_profile()` | да | hits, genome | primers, scheme, seed_len | `hits + genome` → `tibble` | Grade an off-target locus by its 3'-end mismatches |
| `specificity_blast()` | — | primers, db_path | k, BPPARAM | `primers + index` → `hits` | Genome-wide primer hit counts via a BLAST database (the indexed path) |
| `specificity_bowtie()` | — | primers, index | max_mismatches, nthreads | `primers + index` → `hits` | Count genome hits of whole primers with a bowtie1 index (Rbowtie) |
| `specificity_fm()` | — | primers, index | max_mismatches, nthreads | `primers + index` → `hits` | Count genome hits of whole primers with an FM-index (Rsubread) |
| `specificity_of()` | да | primers, genome | k, max_mismatch, BPPARAM | `primers + genome` → `hits` | Specificity by the right path for the genome at hand |
| `specificity_seed()` | да | primers, genome | k, max_mismatch, BPPARAM | `primers + genome` → `hits` | Count genome-wide near-matches of each primer's 3'-terminal k-mer (prefilter) |
| `specificity_verdict()` | да | hits | scheme, genome, pam, spacer_len | `hits + genome` → `verdict` | Grade in-silico PCR products against a scheme's mismatch policy |

### 1.4 Assembly — сборка ассея и мультиплекса

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `assemble_kasp_assay()` | да | genome, snp | variants, flank, strands, min_gc, max_gc, tm_range, clean_3p, max_footprint_variants … +19 | `genome + variant + variants + haplotypes` → `assay` | Assemble and score the best KASP assay for a target SNP |
| `assemble_primer_assay()` | да | genome, snp, primers | variants, flank, scheme, salt_method, cores, BPPARAM, diagnostic_positions, engines | `genome + variant + primers` → `assay` | Assemble an assay from fixed (already given) primers |
| `marker_one_pot()` | да | assays | matrix, max_markers | `marker + matrix` → `one_pot` | Group assays into one-pot multiplexes |
| `one_pot_compatible()` | да | assays | matrix, require | `marker + matrix` → `logical` | Can these assays share one reaction mix? |
| `one_pot_matrix()` | — | — | — | `—` → `matrix` | One-pot compatibility matrix |

### 1.5 Criteria — проверки, пороги и QC

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `amplicon_fold_check()` | да | seq | kind, conditions, min_single_stranded, max_length | `sequence` → `fold` | Amplicon/target folding check (whole amplicon, not the primer) |
| `asq_limits()` | — | — | amplicon, as_tm, com_tm | `—` → `limits` | ASQ acceptance thresholds (the ASQ twin of KASP_LIMITS) |
| `caps_fragments()` | — | seq | enzymes | `sequence` → `fragments` | Digest fragment lengths of one sequence, per enzyme |
| `dimer_thermo()` | да | primers | BPPARAM, conditions | `primers + conditions` → `thermo` | Duplex (dimer) free energy for every primer pair (OligoArrayAux, DNA) |
| `duplex_dg_thermo()` | — | seq1, seq2 | conditions, BPPARAM | `sequence + sequence` → `numeric` | Duplex (oligo pair) dG at a scheme's reaction temperature |
| `find_palindromes()` | — | seq | min_armlength, max_looplength, max_mismatch | `sequence` → `palindromes` | Inverted repeats / palindromes in a primer (library screen) |
| `find_repeats()` | — | seq | units, min_copies | `sequence` → `repeats` | Tandem repeats (simple sequence repeats) in a sequence |
| `fold_amplicon()` | да | seq | kind, conditions | `sequence` → `fold` | Fold a whole amplicon/target, not just a primer |
| `hairpin_dg_thermo()` | да | seq | conditions, BPPARAM | `sequence + conditions` → `numeric` | Hairpin (self-folding) dG at a scheme's reaction temperature |
| `has_palindrome()` | — | seq | min_armlength, max_looplength, max_mismatch | `sequence` → `logical` | Does a primer contain an inverted repeat? (see find_palindromes()) |
| `has_repeat()` | — | seq | units, min_copies | `sequence` → `logical` | Does a sequence contain an SSR? (see find_repeats()) |
| `homeolog_axes()` | — | assay | — | `assay` → `tibble` | The two 3' axes of a KASP candidate, made explicit |
| `internal_duplex_dg()` | — | seq1, seq2 | na, temp | `sequence + sequence` → `numeric` | Nearest-neighbour duplex dG without OligoArrayAux (pure R) |
| `internal_hairpin_dg()` | — | seq1 | na, temp | `sequence` → `numeric` | Hairpin dG of a single strand without OligoArrayAux (pure R screen) |
| `marker_costs()` | — | — | — | `—` → `costs` | A cost schema for the chemistries |
| `marker_criteria()` | да | — | limits, scheme | `—` → `criteria` | The registry of every candidate evaluation |
| `marker_design_checks()` | да | assay | limits, scheme, engines | `marker + limits` → `checks` | Check an assembled assay against the KASP acceptance criteria |
| `marker_explain()` | да | x | scheme | `marker + checks` → `explain` | Why this marker: the design checks as ordered, readable reasons |
| `marker_informativeness()` | да | targets | genotypes, haplotypes, min_lines | `targets + haplotypes` → `informativeness` | How informative a marker is on an available panel |
| `marker_priors()` | — | — | — | `—` → `tibble` | Provenance of every KASP threshold |
| `marker_readout()` | да | marker | — | `marker` → `readout` | The readout a marker's scheme models (or says it does not) |
| `marker_risk()` | да | assay | limits | `marker + limits` → `risk` | An experimenter-facing risk verdict for an assay |
| `marker_site_class()` | — | ref, alts | — | `character + character` → `character` | Classify a site's allele set for KASP design |
| `marker_tubes()` | да | solved | min_cross_run, max_selfcomp | `solved` → `tube_pairs` | Which solved markers can share one tube without interfering |
| `oligo_dg()` | да | seq1 | seq2, na, temp | `sequence` → `numeric` | OligoArrayAux duplex / hairpin free energy |
| `passes_all()` | да | assay | limits | `marker + limits` → `logical` | Do all acceptance checks pass? |
| `pcr_efficiency()` | да | primers, targets | temp, P, ions, cores | `primers + targets` → `numeric` | PCR amplification efficiency of primer/target pairs (DECIPHER) |
| `primer_set_compatibility()` | да | primers | min_cross_run, max_selfcomp, conditions | `primers + conditions` → `compatibility` | Pairwise compatibility of a primer pool |
| `primer_thermo()` | да | primer | BPPARAM | `sequence + conditions` → `thermo` | Full structural profile of one primer |
| `primer_tm()` | да | seq | DMSO, nn_table, conditions, salt_method, dnac | `sequence + conditions` → `numeric` | Nearest-neighbour Tm of oligos (perfect duplex), salt-corrected (TmCalculator) |
| `reaction_conditions()` | — | — | Na, Mg, dNTPs, DMSO, temperature | `numeric` → `conditions` | Reaction conditions for the thermodynamics engines |
| `restriction_enzymes()` | — | — | — | `—` → `enzyme_table` | Bundled restriction-enzyme catalogue (name, site, cut) |
| `scheme_conditions()` | — | — | scheme | `character` → `conditions` | Reaction conditions for one marker scheme |
| `specificity_policy()` | да | — | scheme, length, spacer_len | `character` → `policy` | Mismatch tolerance for a specificity scheme |
| `terminal_dg()` | да | primer | n, BPPARAM | `sequence` → `numeric` | Duplex free energy of a primer's 3'-terminal n-mer (3'-end stability) |
| `thermo_crosscheck()` | — | seq | nn_tables | `sequence` → `tibble` | Independent thermodynamic cross-check of a window's Tm across NN tables |
| `wobble_mismatch()` | — | on_allele, off_allele | — | `character + character` → `logical` | Is the 3'-terminal allelic mismatch a weak "wobble" pair (G:T / A:C)? |

### 1.6 Decision — выбор химии и плана

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `marker_decision_graph()` | да | — | assay, scheme | `marker` → `decision` | The decision tree as a tidygraph, annotated with an assay's outcome |
| `marker_decision_text()` | да | — | scheme, evidence | `character` → `character` | The decision tree of a scheme as text |
| `marker_filter_blocks()` | — | — | — | `—` → `character` | The blocks a technology is composed from |
| `marker_pareto()` | да | x | quality, linkage_risk, anchor, recombination, genotypes, haplotypes, min_lines | `solved + haplotypes` → `pareto` | A Pareto front over marker candidates |
| `marker_plan()` | да | genome | targets, region | `genome + targets` → `plan` | Start a marker-design plan |
| `marker_schemes()` | да | — | — | `—` → `tibble` | The marker chemistries this package can describe |
| `marker_solve()` | да | plan | cache, all, progress | `plan` → `solved` | Solve a marker-design plan |
| `marker_tech()` | да | — | schemes, costs, preset, exclude | `character + costs` → `tech` | Declare which chemistries a plan may use, and their costs |
| `marker_technologies()` | — | — | — | `—` → `tibble` | List the composed technologies of this session |
| `marker_technology()` | — | name | base, params, blocks, specificity, label, verdict, overwrite | `spec` → `technology` | Compose a technology from a base chemistry |
| `rescue_strategy()` | да | — | schemes | `—` → `tibble` | The rescue policy: which chemistries to try, in order |

### 1.7 Artefact — артефакт, trace, исходы, replay

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `marker_artefact()` | да | solved | dir, created | `solved` → `artefact` | Write a self-documenting artefact for a solved plan |
| `marker_artefact_read()` | да | path | — | `path` → `artefact` | Read a persisted artefact back |
| `marker_engine_stamp()` | — | — | — | `—` → `stamp` | A version stamp for an Evaluation or an artifact manifest |
| `marker_evaluation_stamp()` | — | x | — | `any` → `stamp` | Engines that produced an evaluation |
| `marker_failure_registry()` | — | — | dir, recursive | `path` → `registry` | Every rejection the artefacts under a directory recorded |
| `marker_failure_summary()` | — | registry | — | `registry` → `tibble` | How often each rule rejected a candidate |
| `marker_import_rocrate()` | да | path | — | `path` → `artefact` | Read an artefact back from an RO-Crate |
| `marker_outcome_corpus()` | да | dir | — | `path` → `corpus` | The outcomes recorded across every artefact under a directory |
| `marker_outcomes()` | да | artefact | — | `artefact` → `outcomes` | The outcomes recorded for an artefact |
| `marker_record_outcome()` | да | artefact | result, target, notes | `artefact` → `outcomes` | Record what a lab got out of a designed assay |
| `marker_replay()` | да | artefact | cache, progress, tolerance | `artefact` → `replay` | Re-run the design an artefact recorded and compare |
| `marker_reprex()` | да | x | path, preserve | `marker` → `files` | Anonymised, self-contained reproducer for a bug report |
| `marker_trace_from_artefact()` | да | x | — | `artefact` → `trace` | Build a trace from an artefact or a solved plan |
| `marker_trace_migrate()` | — | x | — | `trace` → `trace` | Migrate a record to the current trace schema |
| `marker_trace_payload()` | — | artefact | — | `artefact` → `payload` | Payload a second writer needs to reproduce an artefact byte for byte |
| `marker_trace_read()` | да | path | — | `path` → `trace` | Read a trace from an RDS file |
| `marker_trace_write()` | да | trace, path | — | `trace` → `files` | Write a trace, always at the current schema |

### 1.8 Export — файлы, отчёты, патент

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `decision_render()` | — | graph | style, out | `decision + style` → `files` | Draw a decision graph to SVG |
| `export_primers()` | да | assay, gene | dir | `assay + gene` → `files` | Write IDT, LGC and tidy order sheets for the assay primers |
| `kasp_mix_recipe()` | да | — | n_reactions, mix_uL, stock_uM, as_uM, com_uM, reaction_uL, master_mix_conc, dna_uL … +1 | `numeric` → `recipe` | Tidy KASP primer-mix + reaction recipe |
| `marker_defensive_publication()` | да | assay_or_solved | — | `marker + solved` → `files` | Draft a defensive publication from the designed markers |
| `marker_export()` | да | solved | dir, artefact | `solved` → `files` | Write a solved plan to disk |
| `marker_export_outcomes()` | да | corpus, dir | trustee | `corpus` → `files` | Export outcomes for a data trustee without a network |
| `marker_export_rocrate()` | да | artefact | dir | `artefact` → `files` | Write an artefact as an RO-Crate |
| `marker_patent_export()` | да | assay_or_solved | format, output | `marker + solved` → `files` | Export markers as a WIPO ST.26 sequence listing or a draft application |
| `marker_patent_validate()` | — | path | schema | `path` → `logical` | Validate a WIPO ST.26 sequence listing against the XSD |
| `marker_plot()` | да | solved | — | `solved` → `plot` | Plot a solved plan: which chemistry takes each target |
| `marker_primer_table()` | да | assay, gene | — | `assay + gene` → `tibble` | Tidy table of the assay oligos (with KASP tails) |
| `marker_primers()` | да | marker | — | `marker` → `primers` | The oligos a designed marker would put in a tube |
| `plot_alignment()` | да | hp | type, features | `alignment + features` → `plot` | Plot an alignment of the haplotypes (a base x sequence matrix) |
| `plot_checks()` | да | assay | limits | `assay` → `plot` | Acceptance checks as a thin list of tumblers (feature vs alignment) |
| `plot_decision_tree()` | да | — | x, scheme, out, style | `assay` → `files` | Draw the selection decision tree |
| `plot_design_search_step()` | — | frames, step | — | `frames + integer` → `plot` | One frame of the design-search animation |
| `plot_oligo_structure()` | да | assay | min_armlength, max_looplength | `assay` → `plot` | Hairpins, inverted repeats and SSRs within the assay oligos |
| `plot_primers()` | да | assay, variants, snp | genome, features | `assay + variants + variant + genome + features` → `plot` | Local assay map: styled primers, Tm, and the acceptance checks |
| `plot_self_complementarity()` | — | seq | — | `sequence` → `plot` | Self-complementarity dot-plot of one oligo |
| `plot_snp_screen()` | да | screen | gene, threshold | `screen + gene` → `plot` | Plot a target screen along a gene |
| `plot_variant_landscape()` | да | gene, variants, snp, assay | — | `gene + variants + variant + assay` → `plot` | Gene model with the variant track and a zoom on the design window |
| `render_design_search()` | — | assay, variants, snp | path, fps, hold, cores | `assay + variants + variant` → `files` | Render the design-search animation as a GIF |
| `write_kasp_results()` | да | — | — | `assay + gene` → `files` | Устаревший алиас write_marker_results |
| `write_marker_results()` | да | assay, gene | dir | `assay + gene` → `files` | Write the assay results as a few consolidated files |
| `write_variants_vcf()` | да | x | path | `variants` → `files` | Write a variants tibble as a VCF |

### 1.9 Service — служебные примитивы и движки

| Функция | Польз. | Обязательные | Опциональные | Вход → выход | Назначение |
|---|---|---|---|---|---|
| `as_3p_interval()` | — | snp_pos, strand | window | `numeric` → `interval` | Plus-strand interval of the AS primer's N 3'-terminal bases, SNP excluded |
| `as_anchor()` | — | v_start, ref_len, strand | variant_offset | `numeric` → `integer` | Plus-strand coordinate of the AS primer's 3'-terminal base for a variant |
| `capabilities()` | — | engine | capability | `character` → `capabilities` | Capabilities of a marker engine |
| `chunk_alignment()` | — | chunk | — | `chunk` → `chunk_part` | Часть чанка: выравнивание |
| `chunk_genome()` | — | chunk | — | `chunk` → `chunk_part` | Часть чанка: геном |
| `chunk_haplotypes()` | — | chunk | — | `chunk` → `chunk_part` | Часть чанка: гаплотипы |
| `chunk_reference()` | — | chunk | — | `chunk` → `chunk_part` | A chunk's parts |
| `chunk_variants()` | — | chunk | — | `chunk` → `chunk_part` | Часть чанка: варианты |
| `chunk_weight()` | — | chunk, maf | — | `chunk + numeric` → `chunk` | Weight a chunk's variants by minor-allele frequencies |
| `count_variants_in()` | — | interval, vpos | exclude | `interval + variants` → `integer` | Count variant positions inside an interval (IRanges) |
| `decision_palette()` | — | — | — | `—` → `palette` | Palette of a decision figure: fill by kind, border by status, line by edge |
| `decision_style()` | — | — | palette, font_size, line_height, wrap, pad_x, pad_y, margin, segment_penalty … +8 | `—` → `style` | The look of a decision figure |
| `design_search_frames()` | — | assay, variants, snp | — | `assay + variants` → `frames` | Build the frame table of the design search |
| `genome_chroms()` | да | genome | — | `genome` → `character` | Sequence names of a genome (indexed source or DNAStringSet) |
| `genome_length()` | да | genome, chrom | — | `genome + character` → `integer` | Length of one sequence of a genome |
| `has_blast()` | — | — | — | `—` → `logical` | Is NCBI BLAST+ available (for the indexed genome scan)? |
| `has_bowtie()` | — | — | — | `—` → `logical` | Is the bowtie1 scan available (Rbowtie)? |
| `has_oligoarrayaux()` | — | — | — | `—` → `logical` | Is OligoArrayAux (hybrid-min) available? |
| `is_marker_assay()` | — | x | — | `any` → `logical` | Is x a KASP assay object? |
| `marker_cache_dir()` | — | — | — | `—` → `path` | Per-user cache directory for downloads and indices |
| `marker_chunk()` | — | x | — | `haplotypes` → `chunk` | A local chunk of a haplotype set, in every representation |
| `marker_engine_stamp_field()` | — | — | fmt, stamp | `stamp` → `character` | A one-column form of the engine stamp |
| `marker_engine_versions()` | — | — | — | `—` → `tibble` | Versions of every engine an Evaluation depends on |
| `marker_engines()` | — | — | thermo | `—` → `tibble` | Marker engines and their capabilities |
| `marker_install_tools()` | — | — | dest, version, ask, force | `—` → `path` | Build OligoArrayAux into the user cache |
| `marker_trace_schema()` | — | — | — | `—` → `character` | The current trace schema identifier |
| `mask_positions()` | — | seq, positions | — | `sequence + integer` → `sequence` | Mask positions in a DNA sequence by replacing them with N |
| `normalise_allele()` | — | x | — | `character` → `character` | Normalise an allele string to plain ACGTN |
| `primer_interval()` | — | snp_pos, strand, length | gap, role | `numeric` → `interval` | Plus-strand interval occupied by a primer |
| `scan_question()` | — | x | — | `screen` → `character` | The question a specificity table answers |
| `strand_gc_hint()` | — | genome, chrom, pos | flank | `genome + region` → `numeric` | GC fraction of the flanks on each side of the SNP (EasyKASP direction hint) |
| `with_unafloddat()` | — | path, code | — | `path + code` → `any` | Point OligoArrayAux at its data tables for the duration of code |


### 1.10 S3-методы

S3-методы не входят в граф как отдельные узлы, кроме `+.marker_plan`: это
композиция слоёв плана, без неё `marker_plan()` не собирается. Остальные —
методы `print`/`plot` и внутренние vctrs-методы класса `marker_criteria`.

| Метод | Generic | Роль | Назначение |
|---|---|---|---|
| `+.marker_plan` | `+` | decision | пользовательский |
| `capabilities.default` | `capabilities` | service | служебный |
| `capabilities.genome_index` | `capabilities` | service | служебный |
| `capabilities.marker_engine` | `capabilities` | service | служебный |
| `capabilities.marker_engines` | `capabilities` | service | служебный |
| `marker_as_footprint.default` | `marker_as_footprint` | service | служебный |
| `marker_as_footprint.marker_assay` | `marker_as_footprint` | service | служебный |
| `marker_chunk.DNAMultipleAlignment` | `marker_chunk` | service | служебный |
| `marker_chunk.data.frame` | `marker_chunk` | service | служебный |
| `marker_chunk.list` | `marker_chunk` | service | служебный |
| `marker_chunk.marker_chunk` | `marker_chunk` | service | служебный |
| `marker_map_spec.default` | `marker_map_spec` | service | служебный |
| `marker_map_spec.marker` | `marker_map_spec` | service | служебный |
| `marker_primers.default` | `marker_primers` | service | служебный |
| `marker_primers.marker_arms` | `marker_primers` | service | служебный |
| `marker_primers.marker_assay` | `marker_primers` | service | служебный |
| `marker_primers.marker_castpcr` | `marker_primers` | service | служебный |
| `marker_primers.marker_crispr` | `marker_primers` | service | служебный |
| `marker_primers.marker_dpcr` | `marker_primers` | service | служебный |
| `marker_primers.marker_duplex` | `marker_primers` | service | служебный |
| `marker_primers.marker_lamp` | `marker_primers` | service | служебный |
| `marker_primers.marker_lax` | `marker_primers` | service | служебный |
| `marker_primers.marker_rpa` | `marker_primers` | service | служебный |
| `marker_primers.marker_sanger` | `marker_primers` | service | служебный |
| `marker_primers.marker_snapshot` | `marker_primers` | service | служебный |
| `marker_primers.marker_taqman` | `marker_primers` | service | служебный |
| `plot.marker` | `plot` | export | пользовательский |
| `print.marker` | `print` | export | пользовательский |
| `print.marker_artefact` | `print` | export | пользовательский |
| `print.marker_chunk` | `print` | export | пользовательский |
| `print.marker_crispr_seed_evidence` | `print` | export | пользовательский |
| `print.marker_genome` | `print` | export | пользовательский |
| `print.marker_hrm_resolvability` | `print` | export | пользовательский |
| `print.marker_one_pot` | `print` | export | пользовательский |
| `print.marker_panel` | `print` | export | пользовательский |
| `print.marker_plan` | `print` | export | пользовательский |
| `print.marker_replay` | `print` | export | пользовательский |
| `print.marker_solved` | `print` | export | пользовательский |
| `print.specificity_policy` | `print` | export | пользовательский |
| `print.specificity_verdict` | `print` | export | пользовательский |
| `vctrs::vec_cast.marker_criteria.marker_criteria` | `vctrs::vec_cast` | service | служебный |
| `vctrs::vec_cast.marker_criteria.tbl_df` | `vctrs::vec_cast` | service | служебный |
| `vctrs::vec_cast.tbl_df.marker_criteria` | `vctrs::vec_cast` | service | служебный |
| `vctrs::vec_ptype2.marker_criteria.marker_criteria` | `vctrs::vec_ptype2` | service | служебный |
| `vctrs::vec_ptype2.marker_criteria.tbl_df` | `vctrs::vec_ptype2` | service | служебный |
| `vctrs::vec_ptype2.tbl_df.marker_criteria` | `vctrs::vec_ptype2` | service | служебный |
| `vctrs::vec_ptype_abbr.marker_criteria` | `vctrs::vec_ptype_abbr` | service | служебный |
| `vctrs::vec_restore.marker_criteria` | `vctrs::vec_restore` | service | служебный |

## 2. Классы данных

Класс — это контракт, а не тип R. `untyped` означает, что функция возвращает
голый `tibble`/`character`/`numeric`, и граф не может проверить стык.

| Класс | Реальный объект | Производителей | Потребителей |
|---|---|---|---|
| `genome` | `marker_genome` (list с индексом) или `DNAStringSet` | 5 | 48 |
| `variant` | одна строка tibble: chrom, pos, ref, alt | 1 | 27 |
| `variants` | tibble вариантов (chrom, pos, ref, alt) | 5 | 19 |
| `gene` | list записи Ensembl (id, chrom, start, end, strand) | 1 | 7 |
| `sequence` | `DNAString` / character длины 1 | 3 | 21 |
| `sequences` | именованный `DNAStringSet` | 1 | 4 |
| `alignment` | `DNAMultipleAlignment` | 1 | 2 |
| `primers` | tibble кандидатов ИЛИ именованный character | 5 | 11 |
| `marker` | list схемы (`marker_kasp`, `marker_hrm`, …) | 16 | 13 |
| `assay` | `marker_assay`: праймеры + ампликон + `$decision` | 3 | 26 |
| `plan` | `marker_plan`: genome + layers | 2 | 2 |
| `solved` | `marker_solved`: report + attempts + signals | 1 | 7 |
| `panel` | `marker_panel`: targets + one_pot + artefact + solved | 1 | 0 |
| `targets` | tibble целей (label, chrom, pos, ref, alt, type) | 1 | 5 |
| `context` | list вариантов/гаплотипов/глубины | 2 | 2 |
| `tech` | list химий, costs и blocks | 1 | 1 |
| `haplotypes` | tibble гаплотипов по окну | 2 | 11 |
| `depth` | tibble DP/AD/GT по образцу | 1 | 2 |
| `hits` | tibble геномных попаданий | 9 | 2 |
| `decision` | `tbl_graph` решающего дерева | 1 | 1 |
| `artefact` | `marker_artefact`: json + 4 CSV | 3 | 6 |
| `trace` | trace-объект RO-Crate | 3 | 2 |
| `corpus` | каталог исходов | 1 | 1 |
| `outcomes` | tibble исходов, привязанный к id | 2 | 0 |
| `report` | tibble по локусу (design_kasp_markers и батчи) | 3 | 0 |
| `files` | пути/каталог на диске | 14 | 0 |

Отдельно: 47 записей имеют класс `none` — это S3-методы без собственного
возврата; среди экспортов таких нет. `character` — 9, `numeric` — 10,
`logical` — 10, `tibble` — 12. Это не ошибка само по себе, но каждый такой
выход требует ручной распаковки, и граф её не видит.

## 3. Правила композиции

1. Ребро существует, если класс выхода производителя входит в множество
   входных классов потребителя. Необязательные типизированные аргументы тоже
   считаются: `assemble_kasp_assay(genome, snp, variants = ...)` принимает
   `variants`.
2. `assay` — подкласс `marker`: собранный KASP-ассе несёт `$scheme` и
   проходит те же проверки, что scheme-объект. Поэтому `assay → marker_*`
   рисуется, а `marker → assay` — нет.
3. Композиция плана идёт через `+`: `marker_targets() + marker_context() +
   marker_tech()` даёт `marker_plan`, который принимает `marker_solve()`.
   Отдельного конструктора «собери план из трёх слоёв» нет.
4. `marker_design_panel()` — оркестратор: внутри он делает plan → solve →
   one-pot → artefact и возвращает `marker_panel`. Его `$solved` — обычный
   `marker_solved`, но других функций, принимающих `marker_panel`, нет.
5. Выход `report` (батч-таблицы) и `panel` не имеет типизированного
   потребителя: это терминальные отчёты, а не промежуточные объекты.

## 4. Точки входа, листья, изоляты

### 4.1 Пользовательские точки входа (24)

Функции, которым нужны только данные (путь, строка, число), без объекта
пакета:

`ath_genome()`, `ensembl_gene()`, `kasp_mix_recipe()`, `marker_artefact_read()`, `marker_bsgenome()`, `marker_criteria()`, `marker_decision_text()`, `marker_ensdb()`, `marker_example()`, `marker_features()`, `marker_genome()`, `marker_import_rocrate()`, `marker_outcome_corpus()`, `marker_schemes()`, `marker_trace_read()`, `marker_txdb()`, `read_genbank()`, `read_genome_fasta()`, `read_transcripts()`, `read_variant_depth()`, `read_variants_vcf()`, `rescue_strategy()`, `specificity_policy()`, `variant_type()`.

Часть из них — повторный вход в готовый артефакт: `marker_artefact_read()`,
`marker_import_rocrate()`, `marker_trace_read()`, `marker_outcome_corpus()`.

### 4.2 Служебные точки входа (34)

`as_3p_interval()`, `as_anchor()`, `asq_limits()`, `capabilities()`, `decision_palette()`, `decision_style()`, `fasta_index()`, `has_blast()`, `has_bowtie()`, `has_oligoarrayaux()`, `is_marker_assay()`, `marker_cache_dir()`, `marker_costs()`, `marker_engine_stamp()`, `marker_engine_versions()`, `marker_engines()`, `marker_evaluation_stamp()`, `marker_failure_registry()`, `marker_filter_blocks()`, `marker_ingest_literature()`, `marker_ingest_registry()`, `marker_install_tools()`, `marker_patent_validate()`, `marker_priors()`, `marker_site_class()`, `marker_technologies()`, `marker_trace_schema()`, `normalise_allele()`, `one_pot_matrix()`, `primer_interval()`, `reaction_conditions()`, `restriction_enzymes()`, `scheme_conditions()`, `wobble_mismatch()`.

Это предикаты движков, пути кэша, стили и пороги; они не начинают
пользовательский сценарий.

### 4.3 Листья (116)

Пользовательских листьев 64: это файлы, графики, отчёты,
вердикты и артефакты. Полный список — в CSV (`leaf = TRUE`). Типичные:
`write_marker_results()`, `export_primers()`, `marker_export()`,
`marker_artefact()`, `marker_patent_export()`,
`marker_replay()`, `marker_outcomes()`, `plot_checks()`.

### 4.4 Изоляты (29)

Ни один типизированный выход не входит и ни один типизированный вход не
выходит:

`as_anchor()`, `capabilities()`, `decision_palette()`, `fasta_index()`, `has_blast()`, `has_bowtie()`, `has_oligoarrayaux()`, `is_marker_assay()`, `kasp_mix_recipe()`, `marker_cache_dir()`, `marker_criteria()`, `marker_decision_text()`, `marker_engine_versions()`, `marker_engines()`, `marker_filter_blocks()`, `marker_install_tools()`, `marker_patent_validate()`, `marker_priors()`, `marker_schemes()`, `marker_site_class()`, `marker_technologies()`, `marker_technology()`, `marker_trace_schema()`, `normalise_allele()`, `rescue_strategy()`, `specificity_policy()`, `variant_type()`, `with_unafloddat()`, `wobble_mismatch()`.

Большинство — намеренные примитивы (`has_*`, `marker_cache_dir`,
`with_unafloddat`, `as_anchor`). Реальные дефекты среди них:

- `marker_schemes()`, `rescue_strategy()`, `marker_technologies()`,
  `marker_priors()` возвращают `tibble`; чтобы передать их в дизайн, надо
  вручную вынуть `$scheme`. Типизированного моста нет.
- `marker_criteria()` возвращает таблицу критериев, но `marker_design_checks()`
  принимает `limits`, а не этот объект.
- `marker_technology()` создаёт объект технологии, но ни `marker_tech()`, ни
  `marker_solve()` его не принимают.
- `marker_decision_text()` возвращает текст, а не объект `decision`; это
  терминальный лист, но класс в CSV — `character`, поэтому ребра нет.
- `restriction_enzymes()` теперь связан с `find_caps_enzymes()` через
  `enzyme_table`, но `design_caps_marker()` принимает только вектор имён.

## 5. Дубликаты и перекрытия

| Что | Где | Статус |
|---|---|---|
| `design_marker()` → `design_auto_marker()` | `R/marker_target.R:668`, `.Deprecated` | алиас, оставлен для совместимости |
| `write_kasp_results()` → `write_marker_results()` | `R/marker_plots_output.R:92`, `.Deprecated` | алиас |
| `design_kasp_markers()` / `design_markers()` / `design_marker_batch()` | `R/marker_batch.R`, `R/marker_strategy.R` | три батч-входа с разной степенью сохранения объектов |
| `plot_decision_tree()` vs `marker_decision_graph()` + `decision_render()` | `R/marker_decision.R:407` | обёртка; первый рисует файл, второй даёт граф |
| `marker_primer_table()` / `marker_primers()` / `export_primers()` | `R/marker_plots_output.R`, `R/marker_tubes.R` | таблица / имена / файлы — три представления одних праймеров |
| `marker_export()` vs `write_marker_results()` | `R/marker_export.R`, `R/marker_plots_output.R` | bundle `solved` против CSV одного ассея |
| `marker_patent_export()` vs `marker_defensive_publication()` | `R/marker_patent.R` | один вход `assay_or_solved`, разные выходы |
| `find_caps_enzymes()` vs `find_dcaps_enzymes()` | `R/marker_caps.R` | одна форма, разная логика ферментов |

## 6. Скрытые зависимости

Эти функции читают или пишут состояние, которого нет в аргументах:

- **Сессионный реестр технологий:** `marker_technology()`,
  `marker_technologies()`, `marker_filter_blocks()`.
- **Диск/кэш:** `marker_cache_dir()`, `ath_genome()`, `ensembl_gene()`,
  `ensembl_variants()`, `marker_ensdb()`, `marker_txdb()`, `marker_bsgenome()`.
- **Переменная окружения:** `with_unafloddat()` (`UNAFOLDDAT`).
- **Опции:** `tidymarkers.cores`, `tidymarkers.thermo_engine`.
- **Файлы по умолчанию:** `marker_record_outcome()` и `marker_replay()` пишут
  рядом с артефактом; `marker_failure_registry(dir = ".")`.
- **Реестр литературы:** `marker_ingest_*()` использует путь по умолчанию.

## 7. Подтверждённые рёбра

Рёбра ниже проверены примером, виньеткой или тестом. Остальные
906 рёбер графа — теоретическая совместимость классов.

<!-- confirmed-edges:start -->
- `read_genome_fasta` -> `assemble_kasp_assay` — README quickstart; vignettes/kasp-intro.Rmd
- `marker_example` -> `design_kasp_markers` — README batch example
- `ensembl_variants` -> `pick_target_snp` — README web quickstart
- `pick_target_snp` -> `assemble_kasp_assay` — README web quickstart
- `assemble_kasp_assay` -> `marker_design_checks` — README quickstart; test-decision
- `assemble_kasp_assay` -> `passes_all` — README quickstart
- `assemble_kasp_assay` -> `marker_risk` — roxygen example в R/marker_risk.R
- `assemble_kasp_assay` -> `write_marker_results` — roxygen; vignettes/reading-results.Rmd
- `assemble_kasp_assay` -> `export_primers` — roxygen; test-export
- `marker_targets` -> `marker_plan` — test-plan
- `marker_tech` -> `+.marker_plan` — test-plan
- `marker_context` -> `+.marker_plan` — test-plan
- `+.marker_plan` -> `marker_solve` — test-plan
- `marker_plan` -> `marker_solve` — test-plan
- `marker_solve` -> `marker_pareto` — test-pareto
- `marker_solve` -> `marker_artefact` — test-plan; test-outcomes-crate
- `marker_solve` -> `marker_export` — test-plan
- `marker_solve` -> `marker_tubes` — test-plan
- `marker_solve` -> `marker_patent_export` — test-patent
- `marker_artefact` -> `marker_record_outcome` — test-plan; test-z8-outcomes
- `marker_artefact` -> `marker_replay` — test-plan
- `marker_artefact` -> `marker_trace_from_artefact` — test-rocrate
- `marker_outcome_corpus` -> `marker_export_outcomes` — test-outcomes-crate
- `read_variants_vcf` -> `design_kasp_markers` — roxygen; vignettes/custom-genome.Rmd
- `specificity_seed` -> `specificity_verdict` — test-specificity
- `in_silico_pcr` -> `specificity_verdict` — test-specificity
- `marker_design_checks` -> `marker_explain` — test-decision
- `specificity_of` -> `specificity_verdict` — test-specificity
<!-- confirmed-edges:end -->

## 8. Дефекты поверхности

1. **`panel` — тупик.** `marker_design_panel()` возвращает `marker_panel`, но
   ни одна экспортированная функция его не принимает; работать дальше можно
   только через `panel$solved` или `panel$artefact` — это не типизированный
   стык.
2. **`outcomes` — тупик.** `marker_record_outcome()` и `marker_outcomes()`
   дают таблицу исходов, но ни `marker_priors()`, ни `marker_plan()`, ни
   `marker_solve()` её не читают. Петля «результат → следующий дизайн»
   отсутствует.
3. **`report` — тупик.** Батч-функции возвращают `report`, но из него нельзя
   типизированно достать отдельный `marker`/`assay` для проверки или экспорта.
4. **`primers` — два разных объекта.** `design_as_primers()` и
   `design_common_primers()` возвращают tibble кандидатов; `marker_primers()`
   возвращает именованный `character`. Один класс на два контракта.
5. **Справочники без моста.** `marker_schemes()`, `rescue_strategy()`,
   `marker_technologies()`, `marker_priors()` возвращают `tibble`; дизайн
   принимает `character`/`costs`, поэтому пользователь вручную извлекает
   колонку.
6. **Один `variant` — только из `pick_target_snp()`.** `read_variants_vcf()`
   даёт `variants`, но выбрать строку для `assemble_kasp_assay()` можно только
   через `pick_target_snp()` (нужен `gene`) или вручную.
7. **Нет читателя генотипов.** `marker_informativeness()` принимает
   `genotypes`, но `read_variants_vcf()` читает только CHROM/POS/REF/ALT;
   `read_variant_depth()` даёт глубину, не генотипы.
8. **Устаревшие алиасы и три батч-входа** (§5) — не ошибка, но каждый новый
   пользователь выбирает наугад.

## 9. Как поддерживать

Граф и CSV — производные от `NAMESPACE` и контрактов. После добавления
экспорта:

```sh
Rscript dev/check-user-api.R --write   # пересобрать docs/user-api-graph.dot
Rscript dev/check-user-api.R           # гейт: запись есть, классы стыкуются
Rscript dev/check-docs.R               # ссылка из docs/INDEX.md
```

Гейт падает, если экспорт есть в `NAMESPACE`, но нет строки в CSV, или если
подтверждённое ребро не стыкуется по классам. Ручная правка DOT запрещена:
`--write` — единственный генератор.
