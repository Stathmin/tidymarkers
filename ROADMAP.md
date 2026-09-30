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
- [x] `pick_target_snp()`: референс из генома, missense > splice > synonymous > UTR, транзиции, чистая 3'-зона (`0eda4f2`, `ca75edb`, `c59957e`)
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
- [x] Тесты на всё новое; 144 проверки; ретро-чистка пустых ассертов (`011acd1`)

---

## M. Пакет и UX — «душевность» в духе tidyverse  [очередь]

Цель: превратить набор скриптов в настоящий R-пакет с приятным опытом — виньетки,
описания функций, хелперы, примеры, подсказки. Основа уже положена (roxygen-заголовки,
`.data`-стайл, общий визуальный слой), ниже — что довести.

### M1. Каркас пакета
- [x] `DESCRIPTION` с `Imports`/`Suggests` и `Config/testthat/edition: 3`
- [x] roxygen-заголовки: `@export`, `@param`, `@return` на ключевых функциях
- [ ] Собрать реальный пакет: `NAMESPACE` + `man/` из roxygen (`roxygen2::roxygenise()`),
      `R CMD check` без ошибок — тогда `library(kaspexample)` вместо `source()`
- [~] `@family`-теги по подтипам (design / io / specificity / checks / plots / order)
      — проставлены на всех экспортируемых функциях (KASP geometry/input/design/
      thermodynamics/specificity/QC/checks/output/plots/utilities)
- [~] `@examples` на публичных функциях (там, где нужен геном/сеть — `\dontrun{}`) —
      есть у точек входа (`assemble_kasp_assay`, `design_as_primers`, `design_kasp_markers`,
      `kasp_design_checks`, `kasp_mix_recipe`); остальные — в очереди

### M2. Виньетки
- [x] Нарративный отчёт (`report/kasp_report.Rmd`) как прообраз виньетки
- [ ] `vignettes/kasp-intro.Rmd` — за 15 минут: что такое KASP, зачем декомпозиция
      дизайна, где чья термодинамика
- [ ] `vignettes/custom-genome.Rmd` — свой FASTA/VCF, не Arabidopsis
- [ ] `vignettes/reading-results.Rmd` — как читать `kasp_design_checks()` и графики

### M3. Хелперы и подсказки
- [ ] `@section`-подсказки в докстрингах (типа «не забудьте про `UNAFOLDDAT`»)
- [ ] `cli`-сообщения вместо голого `rlang::warn/inform`: шаги пайплайна, прогресс батча
- [ ] `print()`-метод для объекта теста (сейчас это сырой список) — человекочитаемая сводка
- [ ] Предикаты-хелперы: `is_kasp_assay()`, `passes_all()` — для пайпов и юзер-кода
- [ ] Функция-«быстрый старт»: `kasp_example()` — прогон примера из коробки
- [ ] `withr`-хелперы для временных настроек (например, выключить сеть в примерах)

### M4. Опыт пользователя
- [x] `NEWS.md` — changelog по версиям (tidyverse-конвенция)
- [ ] `README` с бейджами и `devtools::install_github()`-снипетом
- [ ] pkgdown-сайт (`_pkgdown.yml`), сгруппированный по `@family`
- [ ] Стиль-гайд: `styler`/`lintr` в CI
- [ ] `usethis`-флоу для контрибьютора (`CONTRIBUTING.md`, шаблон issue)

### M5. Уроки для лекции (DX из пакета)
- [ ] Скелет-слайды «какая библиотека за что отвечает» из раздела F/handout
- [ ] Живой пример: один SNP → тест → графики за < 2 минут (кэш генома)

## Ссылки
- LGC/KBioscience KASP manual — методичка производителя (пороги, рецепт смеси).
- Zhang, Yang & Wen, *BMC Bioinformatics* 2025, **26**:292
  (doi:10.1186/s12859-025-06322-x) — EasyKASP: SSR-скрин, направление по GC флангов.
- Wickham & Bryan, *R Packages* (2-е изд.) — каркас пакета, roxygen, виньетки.
- `tidyverse`/`usethis` пакетная идеология: маленькие экспортируемые хелперы, примеры
  на каждой функции, «See also» через `@family`, пакетные виньетки как основной учебник.
