# tidykasp

Демонстрационный пайплайн подбора KASP-маркеров для *Arabidopsis thaliana*
на чистом R. Здесь — краткая история значимых изменений; подробный инвентарь и
очередь работ — в `ROADMAP.md`.

## 0.3.0

### Гигиена дерева
- Внешние данные и индексы (геном, варианты, BLAST/FM) перенесены в пользовательский
  кэш `kasp_cache_dir()` (`tools::R_user_dir`), а не в проект; `data/` из проекта убран.
- Лекционный дек Quarto (`lecture/`) с `sh lecture/render.sh`; всё генерируемое — в
  `_site/`, в `lecture/` только исходники.
- `plot_discrimination()` переписан: выбранное размещение подсвечено, порог приёмки
  линией, метки значений — и исправлено наложение кандидатов с одинаковой подписью.

### Отбор AS-праймера
- `assemble_kasp_assay()` добавил `seed_hits` (уникальность 3'-сида, обе цепи)
  тай-брейкером после дискриминации/футпринта (`specificity_candidates = 10`):
  при прочих равных предпочитается более «тихая» цепь; скрин только лидеров.

### Анимация поиска
- `render_design_search()` — GIF: перебор AS-размещений в порядке алгоритма (и затем
  COM), каждое подписано причиной отказа, сходятся к принятому тесту; статичные
  дорожка вариантов и мишень. `design_search_frames()` (кадры) и
  `plot_design_search_step()` — библиотечные ggplots; сборка через ImageMagick
  `convert` или `gifski`. Встроено в демо, лекцию и отчёт.
- `render_design_search()` рисует только различные кадры (повторы `hold`
  переиспользуют файл) и может распараллелить отрисовку на Unix
  (`options(tidykasp.cores = ...)`) — ~2× на многоядерной машине; по умолчанию один
  поток (как требует CRAN).

### Пакет и UX
- Проект собирается как настоящий пакет: `NAMESPACE` и `man/` из roxygen
  (все аргументы документированы, `@family` для «See also», `@examples`).
- Объект теста получил класс `kasp_assay` и `print()`-метод; хелперы
  `is_kasp_assay()`, `passes_all()`, `kasp_example()`.
- Виньетка `kasp-intro`; `NEWS.md`.
- Внешний доступ к FASTA по индексу (`Rsamtools::indexFa`).

### Дизайн
- Движок InDel: аллельные шаблоны и скан стыка (`design_as_primers_indel()`,
  `allele_templates()`); размещение выбирает физика (`DECIPHER`), не ручное правило.
- `pick_target_snp()` учитывает специфичность 3'-сида (`seed_hits`) и пропускает
  комплементарные пары (A/T, C/G); тяжесть последствия — по VEP `impact` или SO-термам.
- Substitution-блоки / MNP: вариант — строка REF/ALT, единая геометрия через
  `as_anchor()` / `block_in_seg()`.
- `design_as_primers_indel()` и SNP-путь унифицированы в `assemble_kasp_assay()`
  (пер-аллельный 3'-футпринт, ампликон — худший из двух аллелей).

### Специфичность
- `specificity_seed()`: точный путь `PDict`/`vcountPDict`, мисматч-скан — параллель
  по хромосомам (`BPPARAM=`); `in_silico_pcr`/`assay_amplicons`/`assemble_kasp_assay`
  принимают `BPPARAM`.
- Батч: `design_kasp_markers()` с вердиктом на локус и параллелью по локусам.
- VCF-регион стримится через `VariantAnnotation::ScanVcfParam(which=)`.

### Термодинамика / QC
- `terminal_dg()` — 3'-концевая стабильность (GC-clamp, `OligoArrayAux`).
- `wobble_mismatch()` — флаг слабого G:T / A:C неспаривания (симметрично по REF/ALT).
- `oligo_dg()` идёт через `processx` с явной диагностикой вместо тихого `NA`.

### Ввод / вывод
- `read_variants_vcf(backend = )`: потоковый `vroom` (дефолт) или канонический
  `VariantAnnotation::readVcf` для сложных VCF (INFO/GT, phased, breakends).
- `normalise_allele()` — IUPAC-коды сворачиваются к `ACGTN`, а не отбрасываются.
- `kasp_mix_recipe()` — tidy-рецепт смеси (AS:AS:COM = 2:2:5) и реакции.

### Масштаб / большие геномы
- `kasp_genome()` — индексированный источник (`Rsamtools`): регионы читаются по
  запросу, без всей ДНК в памяти; `genome_seq()`/`genome_chroms()`/`genome_length()`
  совместимы с `DNAStringSet`. `snp_template()`, `strand_gc_hint()`, `in_silico_pcr()`
  работают с обоими представлениями.
- `kasp_blast_index()` + `specificity_blast()` — индексный (BLAST/FM-подобный) скан
  генома; `specificity_of()` выбирает путь сам. `dimer_thermo(BPPARAM=)` — параллель.
- `kasp_genome()` берёт имена/длины из индекса (`seqinfo()`), без ручного разбора `.fai`.
- `kasp_genome(bsgenome = )` — источник из Bioconductor-объекта (`BSgenome::getSeq()`),
  без FASTA на диске.

### Аннотация и экосистема (2026)
- `annotate_variants()` — настоящая аннотация из TxDb (`VariantAnnotation::locateVariants()`,
  `predictCoding()` при `seqsource`); словарь `so_rank()` перестал быть единственным
  источником. `kasp_txdb()`/`kasp_bsgenome()` грузят пакеты по умолчанию (в `Suggests`).
- `find_palindromes()`/`has_palindrome()` — библиотечный скрин инвертированных повторов
  (`Biostrings::findPalindromes()`), 14-й критерий приёмки.
- **FM-индекс** (`Rsubread`): `kasp_fm_index()`/`specificity_fm()` — индексный
  (без BLAST) сквозной скан целых праймеров; `kasp_genome(fm_index = )`,
  `assemble_kasp_assay(run_blast = TRUE)` использует индекс при наличии.
- `kasp_ensdb()` — EnsDb из GTF Ensembl Plants (актуальная модель); consequence-аннотация
  остаётся на TxDb (у ensembldb нет `locateVariants`, мост к TxDb теряет CDS).
- Экосистема проверена и задокументирована (см. ROADMAP P): `oligo` — микрочипы, не
  подходит; `rprimer`/`pwalign`/`targets`/`biomaRt` — оценены.

### Инструменты разработчика
- `dev/check.R` — локальный «CI» (доки ↔ roxygen, styler, lintr, spelling;
  `--full` → `R CMD check`); `dev/bench.R` — время и память по функциям
  (`bench`/`Rprofmem`); `dev/profile-gif.R` — профиль отрисовки GIF.
- Пакет стал lint-clean (`.lintr` в формате, понятном lintr 3.4); тесты можно
  гонять параллельно (`Config/testthat/parallel`).

### Локализация и внешние инструменты
- Все пользовательские сообщения (`stop`/`warning`/`cli`/прогресс) переведены через
  GNU gettext (`po/`, русский каталог `R-ru.po`); выбор по `LANGUAGE`, справка `.Rd`
  остаётся английской (обоснование — ROADMAP Q).
- OligoArrayAux собирается из исходников скриптом `data-raw/install-oligoarrayaux.R`
  (base R, без conda) в `.tools/`; `lecture/render.sh` и тесты
  (`helper-oligoarrayaux.R`) подхватывают сборку, поэтому термодинамические тесты
  больше не скипаются.

## 0.2.0
- `kasp_design_checks()` и `KASP_LIMITS`; пороги LGC/EasyKASP.
- `in_silico_pcr()` (`matchProbePair`), SSR-скрин, `strand_gc_hint()`.
- Пакетный ввод VCF (`vroom`) и FASTA с нормализацией имён.

## 0.1.0
- Каркас `renv`, дизайн AS/COM через `DECIPHER`/`TmCalculator`/`OligoArrayAux`,
  библиотечная графика (`ggmsa`, `gggenes`, `ggrepel`, `patchwork`), экспорт IDT/LGC.
