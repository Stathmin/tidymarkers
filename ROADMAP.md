# tidykasp — инвентарь и roadmap

Живой список. Первая часть — **ретроспективный инвентарь** по подтипам (что уже есть),
вторая — **очередь**. Метод **не упрощаем**: где у нас пакетная
термодинамика (`TmCalculator`, `DECIPHER`, `OligoArrayAux`) — оставляем её, даже если в
сторонних инструментах есть кустарная формула.

Обозначения: `[x]` сделано и покрыто тестами · `[~]` частично · `[ ]` очередь.

A–P — реализовано или оценено; **Q** — локализация (i18n); **R–T** — аудит 2026:
недоиспользованные возможности, переносимость (Windows/macOS), соответствие CRAN/tidyverse.

## A. Инфраструктура и воспроизводимость
- [x] `renv`: `DESCRIPTION` + пиннутый `renv.lock` (живая активация в репозиторий
      намеренно не входит)
- [x] Внешний тул OligoArrayAux собирается из исходников скриптом
      `data-raw/install-oligoarrayaux.R` (base R, без conda) в `.tools/` как
      `bin/` + `share/oligoarrayaux/`; `lecture/render.sh` выставляет
      `PATH`/`UNAFOLDDAT`, а тесты подхватывают сборку через
      `tests/testthat/helper-oligoarrayaux.R`
- [x] Раннер тестов: `tests/testthat.R` (`test_check()` по установленному пакету)
- [x] Внутренние проверки и замеры: `dev/check.R` (доки, styler, lintr, spelling,
      `--full` → `R CMD check`), `dev/bench.R` (время + память по функциям),
      `dev/profile-gif.R` (профиль отрисовки GIF)
- [x] `.gitignore` (кэш, results, .tools)
- [x] Скрипты остаются скриптами, но с roxygen-заголовками — путь к реальному пакету (см. M)

## B. Вход и интероперабельность (`R/kasp_io.R`)
- [x] Ensembl Plants REST: ген-модель и варианты 1001 Genomes, кэш на диск + retry
- [x] `read_variants_vcf()`: потоково через `vroom` (memory-mapped), `.vcf.gz` без распаковки
- [x] Мультиаллельный ALT → строки; символические аллели (`<DEL>`, `*`, `.`) отброшены
- [x] IUPAC/неоднозначные основания → `N` (`normalise_allele()`), регистр нормализуется
- [x] `read_genome_fasta()`: нормализация имён + фильтр по `seqids`
- [x] Загрузка и кэш TAIR10 (`ath_genome()`), нормализация хромосом (`ath_chrom()`)
- [x] Индекс FASTA (`fasta_index()`, `Rsamtools::indexFa`/`scanFa`) — случайный доступ;
      `read_genome_fasta(use_index=)` читает только нужные `seqids`
- [x] Маскировка «прочих» вариантов: учёт в скоринге (`n_var`/`n_3p`) и жёсткий
      предел `max_footprint_variants`; `mask_positions()` для буквальных `N`

## C. Геометрия (единый источник истины)
- [x] `primer_interval()` / `as_3p_interval()` / `count_variants_in()` на `IRanges`
- [x] `as_anchor()` / `block_in_seg()` — SNP и MNP идут одним путём кода
- [x] Геометрия переиспользуется дизайном, сборкой теста и графикой; спан хранится в тесте

## D. Движок AS-праймеров
- [x] `design_as_primers()`: сетка цепь × длина 20–27 × позиция мисматча × основание
- [x] `strands = c("+","-")` — поиск AS-праймера на обеих цепях как аргумент
- [x] Дестабилизирующий (вторичный) мисматч
- [x] Variant-aware отбор: минимизация вариантов под футпринтом (`n_var`), чистая 3'-зона (`n_3p`)
- [x] Substitution-блоки / MNP — замена блока, дискриминация с любой цепи
- [x] `design_as_primers_indel()`: аллельные шаблоны (`allele_templates()`) + скан стыка, размещение выбирает DECIPHER
- [x] Тай-брейк по wobble и 3'-стабильности

## E. Общий (common) праймер
- [x] `design_common_primers()`: скан кандидатов, метрики из пакетов
- [x] Учёт самоструктуры COM (гомодимер/шпилька), `COM Tm >= AS Tm`, ампликон ≤ 100 bp

## F. Термодинамика и структурный QC
- [x] Tm через `TmCalculator` (NN, SantaLucia 2004, солевая/DMSO-коррекция)
- [x] Аллельная дискриминация через `DECIPHER::CalculateEfficiencyPCR` (гибридизация + Taq)
- [x] `oligo_dg()` / `dimer_thermo()` / шпильки через OligoArrayAux (DNA)
- [x] `has_run()` / `find_repeats()` — гомополимеры и SSR
- [x] `terminal_dg()` — 3'-концевая стабильность (GC-clamp)
- [x] `wobble_mismatch()` — слабое G:T / A:C неспаривание (симметрично по REF/ALT)

## G. Специфичность
- [x] `specificity_seed()` — 3'-сид по обеим цепям
- [x] Быстрый точный путь `PDict`/`vcountPDict`; мисматч-скан — параллель по
      хромосомам (`BiocParallel`); `in_silico_pcr` тоже параллелен (см. O)
- [x] `in_silico_pcr()` через `matchProbePair`: REF → 1 продукт, ALT → 0
- [x] Опциональный BLAST (`rBLAST`, `blastn-short`)
- [x] Индексный скан генома: `kasp_genome()` + `specificity_blast()` + диспетчер
      `specificity_of()` (BLAST-индекс); регионный доступ без ДНК в RAM (см. O)
- [x] `strand_gc_hint()` — направление по GC флангов (EasyKASP)

## H. Отбор мишени и приёмка
- [x] `pick_target_snp()`: референс из генома; приоритет — чистая 3'-зона → тяжесть
      последствия (VEP `impact` или канонические SO-термы) → транзиция
- [x] Специфичность 3'-сида учтена при выборе мишени: `seed_hits` (лучшая из двух
      цепей) как критерий среди топ-`specificity_candidates`; комплементарные пары
      (A/T, C/G) пропускаются — их 3'-основание не различает
- [x] `KASP_LIMITS` + `kasp_design_checks()` — tidy `criterion/observed/limit/pass`, пороги LGC/EasyKASP
- [x] Отказ от чужой эмпирической формулы Tm (EasyKASP) в пользу NN-термодинамики
- [x] Отбор пар — «гейт серьёзного, скор минорного»: жёсткие гейты (дискриминация,
      структура, `Tm_COM >= Tm_AS`, ампликон, вариант под COM-футпринтом, вариант в
      3'-окне AS), лексикографический ранг по марже тира (узкое место, `min`), «лестница»
      релаксации с предупреждением на каждый нарушенный гейт (`assay$selection`)
- [x] Взвешенная нагрузка вариантов: необязательный `pop_freq` в таблице вариантов
      утяжеляет футпринт COM и масштабирует маркеры на трэке
- [x] Ожидаемая доля «выпадения» под футпринтом (`DECIPHER::CalculateEfficiencyPCR`):
      вклад каждого варианта — потеря эффективности (гибридизация + 3'-концевая
      дотация Taq), взвешенная частотой; гейт и ранг по COM-футпринту и AS 3'-окну,
      порог `KASP_LIMITS$footprint_dropout`, тумблеры `AS/COM drop`
- [x] Сцепление вариантов (co-segregation): вариант в 3'-окне AS (или футпринте COM)
      *вместе* с мишенью — не «грех», а союзник. Сделано: `screen_snp_targets(linked=)`
      прощает сцепленный вариант в 3'-окне; `assemble_kasp_assay(linked=)` /
      `design_as_primers(linked=)` «впечатывают» сцепленную замену в ALT-гаплотип и
      ALT-праймер (двойной мисматч к REF), REF не трогается; хранится в
      `assay$selection$linked`. Осталось: помечать дизайн, держащийся *только* на
      сцеплении. Гаплотипы/LD — 1001 Genomes per-accession
- [x] Дешёвый векторный скрин мишени `screen_snp_targets()`: по SNP (или целому гену)
      без термодинамики/BLAST/геномных сканов — GC, run/SSR, хотзона вариантов, вариант
      в 3'-окне; вердикт «promising/doomed» и счёт `k/16` (strand x length). ~3 ms/SNP,
      фильтрует гены до `pick_target_snp()`/`assemble_kasp_assay()`/rBLAST; `plot_snp_screen()`
- [x] Мультиаллельные сайты и вырожденные праймеры: `kasp_site_class()` формализует
      критерий «сайт с >2 гаплотипами, всё ещё годными для одного KASP»: одиночный
      ALT, не комплементарный REF → `biallelic`; ≥2 SNP-ALT → `degenerate` (один
      AS-праймер с 3'-концевым IUPAC-кодом объединения ALT, `iupac_code()` из
      `Biostrings::IUPAC_CODE_MAP`, код исключает REF) покрывает все ALT одним KASP;
      одиночный комплементарный ALT (A/T, C/G) → `complementary`, 3'-основанием
      неразделим; не-SNP ALT → `indel` (свой путь). `screen_snp_targets()` получил
      колонку `class`, поэтому мультиаллельный сайт с комплементарным первым ALT
      больше не помечается `doomed`. Фикстура разных типов (SNP, MNP, инсерция/делеция,
      два близких ALT) — `inst/extdata/example_targets.csv`, `kasp_example(type=)`.
      Осталось: сам дизайн вырожденного AS-праймера в движке (IUPAC-ALT в
      `design_as_primers()`/`pick_target_snp()`) и сверка с литературой по multi-allelic KASP
      Гаплотипы/LD — 1001 Genomes per-accession

## I. Батч-дизайн
- [x] `design_kasp_markers()` — строка на локус, вердикт `pass`/`warn`/`fail` + причина
- [x] `variant_type()` — SNP / MNP / insertion / deletion
- [x] `BPPARAM` — параллель по локусам

## J. Визуализация (`R/kasp_plots.R`, `R/kasp_search_anim.R`, пакетные библиотеки)
- [x] `plot_primers()` — локальная карта ассая: буквы последовательностей (нейтрально
      чёрные), **жирные** — 3' GC-clamp, **красные** — основания SSR/шпильки; кольцо =
      мисматч/дискриминирующее 3'-основание; Tm числом + ΔTm-тумблер; дорожка вариантов
- [x] `plot_checks()` — тонкая колонка именованных тумблеров по критерию, фасет
      feature/alignment (правая панель `plot_primers()`)
- [ ] **Контекстный вид** (крупный масштаб): геном/фича + врезка в регион — через
      `Gviz`/`karyoploteR` или `ggforce::facet_zoom`; сейчас `plot_variant_landscape()`
      — `gggenes`-приближение
- [ ] **Выравнивание на референс** с фичами (не общий MSA): референс-бэкбон, ALT —
      только отличия, праймеры на его координатах, метки критериев рядом
- [ ] План: как показывать **self-complementarity, шпильки, SSR** — dot-plot self×self,
      дуги (`ggforce::geom_arc0`) по `findPalindromes`, дорожка повторов; вместо
      текущего primer×revcomp MSA
- [ ] `plot_discrimination()` — удалён (шумный бар-чарт)
- [x] `plot_variant_landscape()` — модель гена + врезка региона
- [x] `plot_alignment()` — `ggmsa` (ДНК-MSA + логотип), REF/ALT и AS
- [x] `plot_self_complementarity()` — праймер против своего revcomp (`ggmsa`)
- [x] `plot_dimer_qc()` — матрица dG дуплексов + шпильки

- [x] **Анимация поиска** (`design_search_frames()`, `plot_design_search_step()`,
      `render_design_search()`) — GIF: перебор AS-размещений в порядке алгоритма и
      затем COM, каждое с причиной отказа, до принятого теста; статичная дорожка
      вариантов. Кадры — обычные ggplots, сборка GIF через ImageMagick `convert`
      (gganimate требует `transformr`→`sf`/GDAL). Рисуются только различные кадры
      (`hold` переиспользует файл); `cores > 1` (или `options(tidykasp.cores=)`)
      распараллеливает отрисовку на Unix — ~2× (замер: 9→5 с, `dev/profile-gif.R`).

## K. Вывод и заказ
- [x] `export_primers()` — tidy + sheets для IDT (bulk) и LGC/KASP
- [x] `kasp_primer_table()` — join-таблица олиго с хвостами
- [x] `kasp_mix_recipe()` — рецепт смеси (AS:AS:COM = 2:2:5) и реакции
- [x] `report/kasp_report.Rmd` — нарративный отчёт (квейзи-виньетка)

## L. Гигиена
- [x] Идиоматический tidyverse-стайл: `list_rbind`/`crossing`/`.data` вместо `*_dfr`/`expand.grid`
- [x] Тесты на всё новое; 170 проверок; ретро-чистка пустых ассертов

---

## M. Пакет и UX — «душевность» в духе tidyverse  [сделано]

Цель: превратить набор скриптов в настоящий R-пакет с приятным опытом — виньетки,
описания функций, хелперы, примеры, подсказки. Основа уже положена (roxygen-заголовки,
`.data`-стайл, общий визуальный слой), ниже — что довести.

### M1. Каркас пакета
- [x] `DESCRIPTION` с `Imports`/`Suggests` и `Config/testthat/edition: 3`
- [x] roxygen-заголовки: `@export`, `@param` (все аргументы), `@return` на публичных
- [x] Собрать настоящий пакет: `NAMESPACE` + `man/` из roxygen (`roxygen2::roxygenise()`),
      `R CMD check` — **Status: OK** (0 ERROR / 0 WARNING / 0 NOTE); `library(tidykasp)`
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
- [x] Линт/стиль: конфиг `.lintr` (CI в репозиторий не входит — см. `CONTRIBUTING` и T)
- [x] `CONTRIBUTING.md` — цикл разработки, стиль, документация

### M5. Уроки для лекции
- [x] Лекционный дек Quarto (`lecture/kasp-slides.qmd` + `sh lecture/render.sh`) — живой
      tidyverse-прогон; заменил прежний скелет `inst/slides/library-map.md` (убран)
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
- [x] **Настоящая аннотация** (`R/kasp_annotate.R`): `annotate_variants()` берёт
      consequence из транскрипционной модели (`VariantAnnotation::locateVariants()`,
      плюс `predictCoding()` с `seqsource`), а не из строки входного файла;
      словари `LOCATION_TO_SO`/`VA_TO_SO` сводят оба словаря к тем же SO-термам,
      которые ранжирует `so_rank()`. Словарь больше не «единственный источник».

## O. Масштаб на большие геномы (индексы, параллель, память)

Ориентир: не «забить/округлить зону», а правильно использовать ресурсы.

- [x] **Индекс точного поиска**: `PDict`/`vcountPDict` — один вызов на все сиды.
- [x] **Приблизительный (мисматч) поиск**: проверено, что `PDict(max.mismatch=)` даёт
      `MTB_PDict`, но `vcountPDict` его **не поддерживает** → вместо индекса
      параллель **по хромосомам** (`BiocParallel::bplapply`); `BPPARAM` прокинут в
      `specificity_seed`, `in_silico_pcr`, `assay_amplicons`, `assemble_kasp_assay`,
      `pick_target_snp` (тесты serial == Multicore).
- [x] **Сид-скан выбора мишени**: один вызов `specificity_seed` на все сиды обеих цепей
      всего шортлиста (вместо 2 × N сканов).
- [x] **VCF**: region-стриминг `ScanVcfParam(which=)` для `VariantAnnotation`
      (контиги заголовка сопоставляются через `ath_chrom()`).
- [x] **FASTA**: `fasta_index()` — доступ к регионам без загрузки всего генома.
- [x] **Геном как `FaFile`-источник**: `kasp_genome()` + `genome_seq()`/
      `genome_chroms()`/`genome_length()`. Функции дизайна (`snp_template()`,
      `strand_gc_hint()`, поиск референса в `pick_target_snp()`) читают только нужный
      регион; `in_silico_pcr()` идёт по хромосомам по одной. Оба представления
      (индекс и `DNAStringSet`) взаимозаменяемы — через `genome_seq()`.
- [x] **Индексный скан генома**: `kasp_blast_index()` (NCBI BLAST+ строит свой индекс)
      + `specificity_blast()`; диспетчер `specificity_of()` сам выбирает путь
      (в памяти → `specificity_seed`, индекс → BLAST). `Rsubread`/`gmapR` — опция.
- [x] **Батч/параллель `OligoArrayAux`**: `dimer_thermo(BPPARAM=)` параллелит
      n(n+1)/2 пар (тесты serial == Multicore).

## P. Экосистема R для KASP — что ещё взять  [почти закрыто]

Уже взято: `Biostrings`, `IRanges`/`GenomicRanges`/`GenomeInfoDb`, `Rsamtools`,
`VariantAnnotation`, `GenomicFeatures` + `TxDb.Athaliana.BioMart.plantsmart28`,
`BSgenome.Athaliana.TAIR.TAIR9`, `ensembldb`, `TmCalculator`, `DECIPHER`, OligoArrayAux,
`rBLAST`, `Rsubread`, `BiocParallel`, `vroom`, `httr2`/`jsonlite`,
`ggmsa`/`gggenes`/`ggrepel`/`patchwork`, `renv`, `testthat`, `cli`/`withr`.

Оценки (с вердиктом):

- [x] **Структурированные последствия**: `GenomicFeatures`/`VariantAnnotation` +
      `TxDb.Athaliana.BioMart.plantsmart28` — `annotate_variants()` (см. N4) заменяет
      словарь `so_rank()` настоящей аннотацией. Принято.
- [x] **`BSgenome`-аннотации** (`BSgenome.Athaliana.TAIR.TAIR9`) — геном как
      Bioconductor-объект: `kasp_genome(bsgenome = ...)` читает регионы через
      `BSgenome::getSeq()`, без FASTA на диске. Принято.
- [x] **`Biostrings::findPalindromes()`** — инвертированные повторы/шпильки как
      библиотечный QC-скрин (`find_palindromes()`/`has_palindrome()`, 14-й критерий
      приёмки). Принято.
- [~] **`oligo`** (Bioc) — **не подходит**: это пакет препроцессинга микрочипов Affymetrix;
      в его namespace нет ни одной функции про melting/hybrid/duplex/dimer ΔG.
      Для ΔG дуплексов/шпилек остаётся OligoArrayAux (через DECIPHER). Не принято.
- [~] **`rprimer`** (Bioc, qPCR) — дизайн/штрафы олиго **по выравниванию**
      (`RprimerProfile`/`designOligos`): ищет консервативные регионы, а не
      аллель-специфичный 3'-мисматч. Может служить независимым **валидатором**
      (`RprimerOligo`/`checkMatch`), но требует сборки выравнивания вокруг локуса.
      Оставлено как опция для лекции, в конвейер не встроено.
- [~] **`pwalign`** (`pairwiseAlignment`) — настоящее выравнивание пригодилось бы
      для визуализации инделов (`plot_alignment()`), но `allele_templates()` —
      точная ручная склейка кодирующего блока, для SNP/MNP/InDel достаточная.
      Оставлено на случай парных выравниваний аллелей.
- [~] **`targets`** — воспроизводимый инкрементальный пайплайн вместо скрипта
      (тысячи локусов через `design_kasp_markers()`). Установлен; обёртка-шаблон
      `_targets.R` — задел для лекции, в пакет не входит (не часть API).
- [~] **`biomaRt`** — канонический R-интерфейс к Ensembl, но variation-overlap
      (`overlap/region`) он не умеет → наш `httr2`/`jsonlite` остаётся оправданным.
      Не принято (функционально ничего не добавляет).
- [x] **FM-индекс** (`Rsubread`) — `kasp_fm_index()` строит индекс (`buildindex`),
      `specificity_fm()` считает локъ-хиты картой целого праймера (`align`,
      многопточно, геном не грузится в память). `kasp_genome(fm_index = )` —
      индекс прикрепляется, `assemble_kasp_assay(run_blast = TRUE)` использует его
      вместо BLAST. Принято. `gmapR` не берём (нужны внешние Gmap/Gsnap).
- [~] **`ensembldb`** — **не как бэкенд аннотации**: пресобранного EnsDb для
      Arabidopsis нет; в `ensembldb` нет `locateVariants()`/`predictCoding()`, и
      `EnsDb` не приводится к `TxDb`; мост через `txdbmaker::makeTxDbFromGRanges()`
      **теряет CDS** (проверено: транскрипты есть, CDS — нет). Значит
      consequence-аннотация остаётся на `TxDb` (N4). Наш вклад:
      `kasp_ensdb()` — собранный/кэшированный `EnsDb` из GTF Ensembl Plants
      (актуальная модель как источник для запросов генов/транскриптов).

Оценено и **не подходит**: `DECIPHER::DesignPrimers` (нужны `tiles`/группы, см. N2);
`oligo` (микрочипы, см. выше). Готового R-пакета под KASP нет (`primeR` снят,
`openPrimeR` требует MAFFT/ViennaRNA/MELTING) — значит наш слой дизайна и порогов
оправдан, а брать надо вокруг: аннотации, индексы, оркестровка.

## Q. Локализация (i18n) — сообщения и документация  [сообщения — сделано; справка — очередь]

Разбор + реализация: все пользовательские сообщения обёрнуты в `gettext()` и переведены
(`po/R-ru.po`, 64 строки; проверено end-to-end — `gettext` + `potools` → `LANGUAGE=ru`
отдаёт русский, `en` — английский). Справка (`.Rd`) намеренно не переводится. «Русские инструкции» — это два независимых вопроса: текст,
который функция **печатает** (ошибки/warning/`cli`/прогресс), и текст **справки**
(roxygen/`.Rd`, виньетки, README). Ниже — 5 пунктов по best practices.

- [x] **1. Сообщения во время выполнения — только gettext, не хардкод.** Исходники
      пишем на одном языке (английский), переводы — в `po/R-ru.po` (домен
      `R-tidykasp`), компиляция в `.mo` при установке; R выбирает по локали/`LANGUAGE`.
      Строки в `stop()`/`warning()`/`message()`/`packageStartupMessage()` переводятся
      автоматически; остальное помечаем `gettext()`/`gettextf()`/`ngettext()`;
      `cat()` **не** переводится. «Ветвиться по `Sys.getlocale()`» — неверный путь.
      Стиль сообщений (Writing R Extensions §1.7): полные предложения, первая буква
      строчная, без точки; не склеивать из фрагментов `paste()`; числа — `ngettext()`
      (в русском **3** формы, не 2). Нужны NLS-сборка R (по умолчанию) и gettext-tools
      (`msgfmt`, `msgmerge`, …).
- [x] **2. Инструмент — `potools`, а не `tools::`-вручную.** Цикл `po_extract()` →
      `po_create("ru")` → перевод → `po_compile()` (или интерактивный
      `translate_package(languages = "ru")`); полезны диагностики
      `check_cracked_messages()` / `check_untranslated_cat()` / `check_untranslated_src()`.
      Стиль извлечения фиксируем в `DESCRIPTION` (`Config/potools/style: base|explicit`);
      для `explicit` заводим свой `tr_()`-враппер. Проверка — `Sys.setenv(LANGUAGE = "ru")`
      + перезагрузка пакета. Проверено: `R CMD INSTALL` сам `po/` **не** компилирует —
      `.mo` собирает `po_compile()` в `inst/po/<lang>/LC_MESSAGES/`, откуда он и ставится.
- [x] **3. Острый угол: `cli`/`rlang` не участвуют в gettext.** У `cli` нет ни функции
      перевода, ни интеграции (нет в NAMESPACE/NEWS) — `cli_alert_*()`, `cli_abort()`,
      `cli_progress_*` сами не переводятся. `po_extract()` их **видит** как сообщения (в
      прототипе — 25 из кода), но в каталог попадают лишь обёрнутые в `gettext()`; значит
      `rlang::warn()/abort()` и `cli_*` надо оборачивать вручную, а структурные
      `cli_abort(c("msg", i = …, x = …))` — по элементу. `{?ед/мн}` и inline-markup англоцентричны. Варианты:
      (а) оборачивать шаблон целиком `gettext()`/`tr_()` и извлекать вручную; (б) переводить
      строку `gettextf()`/`ngettext()` и отдавать её `cli_text()` на форматирование;
      (в) держать сообщения только на английском (так делают tidyverse-пакеты — переводов
      не шлют). **Выбрано (а):** все сообщения обёрнуты в `gettext()` (см. `po/README.md`).
- [ ] **4. Документация — выбор одного языка, не рантайм-переключение.** Механизма
      «перевести `.Rd`» у пакета нет (каталоги покрывают R/C-сообщения, а не текст справки);
      roxygen2 одноязычен. Язык справки объявляется полем `Language` в `DESCRIPTION`
      (IETF-тег, напр. `Language: ru`; первый язык — основной для HTML-мануала). Виньетки/
      README — отдельным файлом (напр. `README.ru.md`), но это контент, а не слой перевода.
      CRAN ждёт Title/Description на английском. Текущий сплит (англ. код/roxygen/man, рус.
      README/ROADMAP/NEWS/лекция) — норма; русские `.Rd` — сознательный отказ от англ. справки.
- [x] **5. Кодировки (для русского больно).** Источники — UTF-8 + `Encoding: UTF-8` (есть);
      в **коде** держим ASCII (идентификаторы/комментарии — non-ASCII в R флагает `R CMD check`),
      кириллицу — только в строках (а лучше — в `.po`). В `Rd` — `\enc{}{}` для ASCII-фолбэка;
      кириллица в LaTeX/PDF хрупкая (`_R_CYRILLIC_TEX_`, `\usepackage[T2A]{fontenc}`).
      **Не ассертить точный текст сообщений** в тестах: их переводят, а `R CMD check` и так
      гоняет примеры/тесты с `LANGUAGE=en` и коллацией `C`. Вернуть английский в сессии —
      `Sys.setenv(LANGUAGE = "en")`.

## R. Недоиспользованные возможности (dev/QA и библиотеки)
- [x] `renv.lock`: запинены инструменты разработки (`bench`, `profmem`, `lobstr`,
      `spelling`, `styler`, `lintr`, `pkgdown`, `potools` + зависимости) — `dev/README.md`
      теперь воспроизводим
- [x] Шаг орфографии вынесен за флаг `--spelling` (`dev/check.R`): проза русская,
      англ. словарь помечает каждое слово; `inst/WORDLIST` при этом не нужен
- [x] `object_usage_linter` возвращён в `.lintr` (`codetools`: неопределённые
      глобалы/опечатки и мёртвые присваивания) — 0 линтов
- [x] `primer_tm()`: `DMSO` прокинут в `TmCalculator::tm_calculate()` (аргумент +
      `PCR_COND$DMSO`) — F больше не расходится с кодом
- [x] `dev/bench.R` переписан на `bench` (медиана + `mem_alloc`) и `lobstr` (размеры
      объектов) вместо самописного `Rprofmem`/`system.time`
- [x] `dev/bench.R`: базовая таблица `dev/bench-baseline.csv` (`--update`), вывод — с
      множителем времени к базовой и флагом замедления
- [x] Тесты — на `withr`-идиомах: `with_seed`/`local_seed`, `local_tempfile`/
      `local_tempdir`, `local_envvar(LANGUAGE = "en")` вместо `set.seed`/`tempfile`/
      ручных `on.exit(unlink)`
- [x] testthat 3e: `expect_snapshot()` для сообщений (`test-messages.R`, `_snaps/messages.md`)
- [x] `httr2`: общий `ensembl_perform()` с `req_throttle(15)`/`req_retry`/
      `req_user_agent` вместо копипасты (дисковый кэш остаётся свой; `req_cache()` — опция)
- [ ] `cli`: `cli_abort()`/`cli_warn()`/`cli_inform()` с bullet-ами (`i`/`x`) и
      плюрализацией для пользовательских ошибок (сейчас `rlang`-строки + `gettext`)
- [x] `BiocParallel`: общий `bp_lapply()` ставит `progressbar` через `BPOPTIONS`;
      `RNGseed` — на BPPARAM (`MulticoreParam(RNGseed=)`; `bpoptions` его не принимает),
      а ядра дизайна детерминированы, сеять нечего
- [x] Мёртвый/латентный вес `renv.lock`: closure-анализ → удалены 8 листьев
      (`oligo`, `oligoClasses`, `affxparser`, `affyio`, `preprocessCore`, `rprimer`,
      `shinycssloaders`, `shinyFeedback`); `R4RNA`/`ff`/`targets`/`pwalign`/`biomaRt`
      оказались достижимы (оставлены)
- [x] `Biostrings`/`Rsamtools`/`ggmsa`/`gggenes`/`ggrepel`/`patchwork` используются
      хорошо (PDict/`vcountPDict`, `matchProbePair`, `findPalindromes`, `scanFa`)

## S. Переносимость (Windows/macOS)
- [x] **GIF-сборка:** `image_tool()` берёт `magick convert …` (ImageMagick 7) либо IM6
      `convert`, а `System32\convert.exe` отвергает; тот же разбор — в `data-raw/logo.R`
- [x] **OligoArrayAux, `.exe`:** `helper-oligoarrayaux.R` видит и `hybrid-min.exe`,
      так что термотесты на Windows не скипаются молча
- [ ] **OligoArrayAux, сборка:** `./configure`/`make` под Windows требуют Rtools
      (MSYS2/UCRT64) и запуска через `bash -lc` — задокументировать
- [x] **sha256 и пробелы:** `sha256_of()` пробует `sha256sum`/`shasum`/`openssl`;
      `UNAFOLDDAT`/`--prefix` кавычатся; компилятор ищется как `cc`/`gcc`/`clang`
- [x] **`.gitattributes`:** `* text=auto eol=lf`, `*.mo`/`*.rds`/`*.png`/`*.gif` — бинарные
- [ ] **Инструкции/шелл:** README и `lecture/render.sh` — варианты cairo для
      Windows/macOS, `sh` только из Git-Bash/Rtools; `apt-get`-хинты (напр.
      `R/kasp_scan.R`) дополнить кроссплатформенными
- [ ] **Форк:** `parallel::mclapply` в `render_design_search()` (Unix-only) корректен,
      но ради macOS/CRAN рассмотреть `BiocParallel::SnowParam`
- [ ] Проверить на живой Windows/R-hub: сборку `hybrid-min.exe` и паритет `-q`-вывода,
      `PATHEXT`/processx, cairo/`ggmsa`, CRLF у `vroom`

## T. Соответствие CRAN/tidyverse
- [x] `.Rbuildignore`: `^\.i18n-lib$` (каталог попадал в тарбол → NOTE и лишняя
      зависимость `R (>= 3.5.0)`); убран дубль `^.*\.Rproj$`
- [x] `ragg` объявлен в `Suggests` (`ragg::agg_png` давал WARNING «'::' import not
      declared from: 'ragg'»)
- [x] **`@examples`** добавлены всем экспортам: короткие исполнимые у чистых хелперов,
      `\dontrun{}` у сетевых/термо/BLAST/файловых (`ath_genome()`, `ensembl_*()`,
      OligoArrayAux, `specificity_*()`, `write_kasp_results()`…)
- [x] **`@return`** дописан всем экспортам (описание типа/структуры возврата)
- [x] **Запись вне `tempdir()`:** дефолты `write_kasp_results()`/`export_primers()`/
      `render_design_search()` → `tempdir()`; `here` ушёл из `Imports` в `Suggests`
      (нужен только `inst/scripts`)
- [ ] **Стиль сообщений** противоречит Q1: строки с заглавной и точкой
      (`kasp_functions.R`, `kasp_io.R`, `kasp_scan.R`) — привести к «строчная, без точки»
- [x] **i18n / `Language:`** — решено **не** задавать: справочник (`.Rd`) английский,
      и `Language: ru` объявил бы его русским; русская проза — отдельный слой (Q4)
- [x] **Тесты не ассертят переводной текст** напрямую: в блоках с англ. подстроками
      добавлен `withr::local_envvar(LANGUAGE = "en")`
- [ ] **Бейджи/CI:** бейдж `lifecycle: experimental` без зависимости `lifecycle`;
      бейдж R-CMD-check статичный, воркфлоу нет
- [ ] **Мелочи:** `COPYRIGHT HOLDER` в `LICENSE` ≠ автор; нет `CITATION`;
      `Config/roxygen2/version` (8.1.0) ≠ установленному (8.0.0); лишний
      `utils::globalVariables(".data")` при `@importFrom rlang .data`
- [x] Уже соответствует: Title/Description (англ., title case), `Authors@R`,
      License+LICENSE, URL/BugReports, `Encoding`, `|>`/`.data$`, `rlang`/`cli`,
      `@family`/`@noRd`, тесты edition 3 без сети/снапшотов, поставка `po/` + `.mo`

## Ссылки
- LGC/KBioscience KASP manual — методичка производителя (пороги, рецепт смеси).
- Zhang, Yang & Wen, *BMC Bioinformatics* 2025, **26**:292
  (doi:10.1186/s12859-025-06322-x) — EasyKASP: SSR-скрин, направление по GC флангов.
- Wickham & Bryan, *R Packages* (2-е изд.) — каркас пакета, roxygen, виньетки.
- `tidyverse`/`usethis` пакетная идеология: маленькие экспортируемые хелперы, примеры
  на каждой функции, «See also» через `@family`, пакетные виньетки как основной учебник.
- *Writing R Extensions* §1.7 (диагностические сообщения), §1.8 (i18n: C/R-сообщения,
  подготовка переводов), §1.1.1 (`Language`), §2.14 (кодировки).
- `potools` (CRAN; Chirico & Wickham) — виньетки developers/translators: `po_extract()`,
  `translate_package()`, диагностики сообщений.
