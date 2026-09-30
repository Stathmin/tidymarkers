# kasp-example — инвентарь и roadmap

Живой список. Первая часть — **ретроспективный инвентарь** по подтипам (что уже есть,
со ссылкой на коммит), вторая — **очередь**. Метод **не упрощаем**: где у нас пакетная
термодинамика (`TmCalculator`, `DECIPHER`, `OligoArrayAux`) — оставляем её, даже если в
сторонних инструментах есть кустарная формула.

Обозначения: `[x]` сделано и покрыто тестами · `[~]` частично · `[ ]` очередь.
Ссылки — короткие хеши коммитов (`git log --oneline -- kasp-example`).

## A. Инфраструктура и воспроизводимость
- [x] `renv`: `DESCRIPTION` + пиннутый `renv.lock`, `renv`-activated `.Rprofile` (`a234b0f`)
- [x] Вендоринг внешних тулов: OligoArrayAux в `.tools/`, запуск через `.tools/run.sh`
      с `PATH` и `UNAFOLDDAT` (`a234b0f`)
- [x] Раннер тестов без сборки пакета: `tests/testthat.R` (сорсит `R/*.R`) (`59c3f63`)
- [x] `.gitignore` (кэш, results, .tools) (`a234b0f`, `7bb69b0`)
- [x] Скрипты остаются скриптами, но с roxygen-заголовками — путь к реальному пакету (см. M)

## B. Вход и интероперабельность (`R/kasp_io.R`)
- [x] Ensembl Plants REST: ген-модель и варианты 1001 Genomes, кэш на диск + retry (`0eda4f2`)
- [x] `read_variants_vcf()`: потоково через `vroom` (memory-mapped), `.vcf.gz` без распаковки (`aa06d66`)
- [x] Мультиаллельный ALT → строки; символические аллели (`<DEL>`, `*`, `.`) отброшены (`aa06d66`, `ca75edb`)
- [x] IUPAC/неоднозначные основания → `N` (`normalise_allele()`), регистр нормализуется (`c59957e`)
- [x] `read_genome_fasta()`: нормализация имён + фильтр по `seqids` (`aa06d66`)
- [x] Загрузка и кэш TAIR10 (`ath_genome()`), нормализация хромосом (`ath_chrom()`) (`0eda4f2`)
- [x] Индекс FASTA (`fasta_index()`, `Rsamtools::indexFa`/`scanFa`) — случайный доступ;
      `read_genome_fasta(use_index=)` читает только нужные `seqids`
- [x] Маскировка «прочих» вариантов: учёт в скоринге (`n_var`/`n_3p`) и жёсткий
      предел `max_footprint_variants`; `mask_positions()` для буквальных `N`

## C. Геометрия (единый источник истины)
- [x] `primer_interval()` / `as_3p_interval()` / `count_variants_in()` на `IRanges` (`d8e3e27`)
- [x] `as_anchor()` / `block_in_seg()` — SNP и MNP идут одним путём кода (`c59957e`)
- [x] Геометрия переиспользуется дизайном, сборкой теста и графикой; спан хранится в тесте (`d8e3e27`)

## D. Движок AS-праймеров
- [x] `design_as_primers()`: сетка цепь × длина 20–27 × позиция мисматча × основание (`0eda4f2`)
- [x] `strands = c("+","-")` — поиск AS-праймера на обеих цепях как аргумент (`0eda4f2`)
- [x] Дестабилизирующий (вторичный) мисматч (`0eda4f2`)
- [x] Variant-aware отбор: минимизация вариантов под футпринтом (`n_var`), чистая 3'-зона (`n_3p`) (`0eda4f2`, `ca75edb`)
- [x] Substitution-блоки / MNP — замена блока, дискриминация с любой цепи (`c59957e`)
- [x] `design_as_primers_indel()`: аллельные шаблоны (`allele_templates()`) + скан стыка, размещение выбирает DECIPHER (`2ed509a`)
- [x] Тай-брейк по wobble и 3'-стабильности (`fba27ad`)

## E. Общий (common) праймер
- [x] `design_common_primers()`: скан кандидатов, метрики из пакетов (`0eda4f2`)
- [x] Учёт самоструктуры COM (гомодимер/шпилька), `COM Tm >= AS Tm`, ампликон ≤ 100 bp (`ee63dfa`)

## F. Термодинамика и структурный QC
- [x] Tm через `TmCalculator` (NN, SantaLucia 2004, солевая/DMSO-коррекция) (`0eda4f2`)
- [x] Аллельная дискриминация через `DECIPHER::CalculateEfficiencyPCR` (гибридизация + Taq) (`0eda4f2`)
- [x] `oligo_dg()` / `dimer_thermo()` / шпильки через OligoArrayAux (DNA) (`0eda4f2`, `7db2005`)
- [x] `has_run()` / `find_repeats()` — гомополимеры и SSR (`b93aab8`)
- [x] `terminal_dg()` — 3'-концевая стабильность (GC-clamp) (`fba27ad`)
- [x] `wobble_mismatch()` — слабое G:T / A:C неспаривание (симметрично по REF/ALT) (`fba27ad`)

## G. Специфичность
- [x] `specificity_seed()` — 3'-сид по обеим цепям (`d8e3e27`)
- [x] Быстрый точный путь `PDict`/`vcountPDict`; мисматч-скан + `BiocParallel`-хук (`bb5cea5`)
- [x] `in_silico_pcr()` через `matchProbePair`: REF → 1 продукт, ALT → 0 (`b93aab8`)
- [x] Опциональный BLAST (`rBLAST`, `blastn-short`) (`0eda4f2`)
- [x] `strand_gc_hint()` — направление по GC флангов (EasyKASP) (`b93aab8`)

## H. Отбор мишени и приёмка
- [x] `pick_target_snp()`: референс из генома; приоритет — чистая 3'-зона → тяжесть
      последствия (VEP `impact` или канонические SO-термы) → транзиция (`0eda4f2`,
      `ca75edb`, `c59957e`, `d52a253`)
- [x] Специфичность 3'-сида учтена при выборе мишени: `seed_hits` (лучшая из двух
      цепей) как критерий среди топ-`specificity_candidates`; комплементарные пары
      (A/T, C/G) пропускаются — их 3'-основание не различает (`d52a253`+)
- [x] `KASP_LIMITS` + `kasp_design_checks()` — tidy `criterion/observed/limit/pass`, пороги LGC/EasyKASP (`ee63dfa`)
- [x] Отказ от чужой эмпирической формулы Tm (EasyKASP) в пользу NN-термодинамики (`b93aab8`)

## I. Батч-дизайн
- [x] `design_kasp_markers()` — строка на локус, вердикт `pass`/`warn`/`fail` + причина (`b0ff61d`)
- [x] `variant_type()` — SNP / MNP / insertion / deletion (`b0ff61d`)
- [x] `BPPARAM` — параллель по локусам (`80ccc76`)

## J. Визуализация (`R/kasp_plots.R`, пакетные библиотеки)
- [x] `plot_discrimination()` — бар-чарт кандидатов (`ca75edb`)
- [x] `plot_primers()` — стрелки `gggenes` + локальная дорожка вариантов (`0420d26`, `2ed509a`)
- [x] `plot_variant_landscape()` — модель гена + врезка региона (`0420d26`)
- [x] `plot_alignment()` — `ggmsa` (ДНК-MSA + логотип), REF/ALT и AS (`0420d26`)
- [x] `plot_self_complementarity()` — праймер против своего revcomp (`ggmsa`) (`0420d26`)
- [x] `plot_dimer_qc()` — матрица dG дуплексов + шпильки (`0420d26`)

## K. Вывод и заказ
- [x] `export_primers()` — tidy + sheets для IDT (bulk) и LGC/KASP (`0420d26`)
- [x] `kasp_primer_table()` — join-таблица олиго с хвостами (`ca75edb`)
- [x] `kasp_mix_recipe()` — рецепт смеси (AS:AS:COM = 2:2:5) и реакции (`681f38f`)
- [x] `report/kasp_report.Rmd` — нарративный отчёт (квейзи-виньетка) (`7bb69b0`, `9a34a4c`)

## L. Гигиена
- [x] Идиоматический tidyverse-стайл: `list_rbind`/`crossing`/`.data` вместо `*_dfr`/`expand.grid` (`ca75edb`)
- [x] Тесты на всё новое; 170 проверок; ретро-чистка пустых ассертов (`011acd1`)

---

## M. Пакет и UX — «душевность» в духе tidyverse  [сделано]

Цель: превратить набор скриптов в настоящий R-пакет с приятным опытом — виньетки,
описания функций, хелперы, примеры, подсказки. Основа уже положена (roxygen-заголовки,
`.data`-стайл, общий визуальный слой), ниже — что довести.

### M1. Каркас пакета
- [x] `DESCRIPTION` с `Imports`/`Suggests` и `Config/testthat/edition: 3`
- [x] roxygen-заголовки: `@export`, `@param` (все аргументы), `@return` на публичных
- [x] Собрать настоящий пакет: `NAMESPACE` + `man/` из roxygen (`roxygen2::roxygenise()`),
      `R CMD check` — **Status: OK** (0 ERROR / 0 WARNING / 0 NOTE); `library(kaspexample)`
- [x] `@family`-теги по подтипам (geometry / input / design / thermodynamics /
      specificity / QC / checks / output / plots / utilities) — собираются в «See also»
- [x] `@examples` на публичных функциях (сетевые/геномные — в `\dontrun{}`)

### M2. Виньетки
- [x] Нарративный отчёт (`report/kasp_report.Rmd`) как прообраз виньетки
- [x] `vignettes/kasp-intro.Rmd` — быстрый старт, карта библиотек, тяжёлые шаги как код
- [x] `vignettes/custom-genome.Rmd` — свой FASTA/VCF, не Arabidopsis
- [x] `vignettes/reading-results.Rmd` — как читать `kasp_design_checks()` и файлы/графики

### M3. Хелперы и подсказки
- [x] `print()`-метод для объекта теста (класс `kasp_assay`) — человекочитаемая сводка
- [x] Предикаты-хелперы: `is_kasp_assay()`, `passes_all()`
- [x] Функция-«быстрый старт»: `kasp_example()` — пример таблицы мишеней из коробки
- [x] `cli`-сообщения: загрузка генома, распаковка для индекса, прогресс батча
      (`cli_progress_bar` в `design_kasp_markers()`)
- [x] `@section`-подсказки в докстрингах (кэш генома, шаги дизайна, `UNAFOLDDAT`)
- [x] `with_unafloddat()` — временно указать таблицы OligoArrayAux (`withr::with_envvar`)

### M4. Опыт пользователя
- [x] `NEWS.md` — changelog по версиям (tidyverse-конвенция)
- [x] `README` с бейджами и снеппетом `remotes::install_github()`
- [x] pkgdown-сайт (`_pkgdown.yml`), сгруппированный по `@family` (`has_concept`)
- [x] Линт/стиль: конфиг `.lintr` + шаблон CI (`.github/workflows/R-CMD-check.yaml`)
- [x] `CONTRIBUTING.md` — цикл разработки, стиль, документация

### M5. Уроки для лекции
- [x] Скелет-слайды «какая библиотека за что отвечает» (`inst/slides/library-map.md`)
- [x] Живой пример: `inst/scripts/quickstart.R` (один SNP → тест → вердикт, на кэше)

## N. Честный аудит: библиотеки vs самописное  [местами очередь]

Где мы **делегируем** науку пакетам (хорошо):

| Задача | Библиотека |
|--------|------------|
| Tm (nearest-neighbour) | `TmCalculator` |
| Аллельная дискриминация (гибридизация + Taq) | `DECIPHER::CalculateEfficiencyPCR` |
| Шпильки/димеры (dG) | OligoArrayAux (движок, что и у DECIPHER) |
| Специфичность | `Biostrings` (`PDict`, `vcountPattern`, `matchProbePair`), `rBLAST` |
| Последовательности/интервалы | `Biostrings`, `IRanges`, `GenomicRanges` |
| Ввод | `vroom` (VCF), `Rsamtools` (индекс FASTA), `httr2`/`jsonlite` (REST) |
| Графика | `ggmsa`, `gggenes`, `ggrepel`, `patchwork`, `ggplot2` |
| Инфра | `renv`, `testthat`, `BiocParallel`, `cli`, `styler`/`lintr` |

Где **не хватило хороших либ** — есть библиотечный канон, который мы не взяли:

- [x] *(закрыто, N1)* **VCF-парсинг**: теперь есть второй бэкенд —
      `read_variants_vcf(backend = "VariantAnnotation")` через `readVcf` (INFO/GT,
      phased, breakends, missing ALT); `vroom` остаётся дефолтом для больших файлов.
- [x] *(оценено, N2)* **COM-праймер**: `DECIPHER::DesignPrimers()` изучен — требует
      `tiles` (наборы групп последовательностей) и целится в пары праймеров
      под покрытие, а не в reverse-праймер напротив SNP; для KASP не подходит.
      Остаётся своя сетка + библиотечный скоринг `DECIPHER`. См. N2.
- [x] *(закрыто, N3)* **dG олиго**: `oligo_dg()` теперь идёт через `processx`
      с явной проверкой статуса и ошибкой вместо тихого `NA` (см. N3).
- [x] *(закрыто, N4)* **Consequence**: `so_rank()` теперь принимает Ensembl VEP
      `impact` (лестница HIGH/MODERATE/LOW/MODIFIER) либо мапит канонические SO-термы
      по полному токену (не подстрокой), включая comma-списки. См. N4.

Где самописное — **оправдано** (KASP-специфика или тонкая клейка, либы нет):

- [x] `design_as_primers()`/`design_as_primers_indel()` — сетка + вовлечение физики
      `DECIPHER`; `openPrimeR`/`primeR` сняты, готового KASP-пакета нет — это наш слой.
- [x] `KASP_LIMITS`/`kasp_design_checks()` — пороги LGC/EasyKASP, специфика KASP.
- [x] `primer_interval()`/`as_anchor()` — арифметика интервалов на `IRanges`.
- [x] `find_repeats()`/`has_run()` — SSR/гомополимеры регулярками (в нашем стеке
      библиотечного аналога нет).
- [x] `wobble_mismatch()`, `terminal_dg()`, `kasp_mix_recipe()`, `write_kasp_results()` —
      мелкая логика/арифметика/укладка таблиц.

### N1. Опциональный бэкенд `VariantAnnotation::readVcf`  [сделано]
- [x] `read_variants_vcf(backend = c("vroom", "VariantAnnotation"))`: `vroom` по
      умолчанию, `readVcf` для сложных VCF; оба возвращают одинаковый tidy-tibble.
- [x] `VariantAnnotation`/`SummarizedExperiment` добавлены в `Suggests` и в `renv.lock`
      (`snapshot.type = "all"`); тест на эквивалентность двух бэкендов на маленьком VCF.

### N2. `DECIPHER::DesignPrimers` как библиотечная альтернатива COM  [оценено, не принято]
- [x] Изучен API и назначение: работает с `tiles` (`TileSeqs` по базе
      последовательностей) и проектирует пары праймеров по региону для дискриминации
      групп; неприменим к окну одного референса и к reverse-праймеру напротив SNP.
- [x] Решение: **не адаптируем**. Физику COM уже даёт `DECIPHER`
      (`CalculateEfficiencyPCR` для дискриминации) и `TmCalculator` (Tm); наша
      переборная сетка — только перечисление кандидатов и KASP-специфичные фильтры.

### N3. Надёжный dG вместо парсинга stdout  [сделано]
- [x] `oligo_dg()` идёт через `processx::run(error_on_status = FALSE)`: явная ошибка
      при ненулевом статусе и при непарсимом выводе (вместо тихого `NA`).
- [x] `processx` в `Imports`; задокументировано поведение (warning + `NA`).

### N4. Структурированный consequence  [сделано]
- [x] `so_rank(consequence, impact = NULL)`: если есть колонка `impact` (Ensembl VEP),
      берётся стандартная лестница `IMPACT_RANK`; иначе consequence режется по запятым
      и мапится по **полному** каноническому SO-терму (`SO_TERM_RANK`), бepётся самый
      тяжёлый. `pick_target_snp()` использует `impact`, если он есть в таблице.
- [x] Тесты: канонические термы, comma-списки, fallback на неизвестном терме,
      приоритет impact и выбор более тяжёлого варианта в `pick_target_snp()`.

## Ссылки
- LGC/KBioscience KASP manual — методичка производителя (пороги, рецепт смеси).
- Zhang, Yang & Wen, *BMC Bioinformatics* 2025, **26**:292
  (doi:10.1186/s12859-025-06322-x) — EasyKASP: SSR-скрин, направление по GC флангов.
- Wickham & Bryan, *R Packages* (2-е изд.) — каркас пакета, roxygen, виньетки.
- `tidyverse`/`usethis` пакетная идеология: маленькие экспортируемые хелперы, примеры
  на каждой функции, «See also» через `@family`, пакетные виньетки как основной учебник.
