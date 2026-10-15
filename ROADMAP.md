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
- [x] Внутренние проверки и замеры: `dev/check.R` (доки, styler, lintr, фигуры,
      spelling, `--full` → `R CMD check`), `dev/bench.R` (время + память по функциям),
      `dev/profile-gif.R` (профиль отрисовки GIF)
- [x] **Фигуры как код, без дрейфа:** `dev/render-all.R` — один вход для всех
      генерируемых картинок (деревья решений → `man/figures` + `results/plot-smoke`,
      панель мишеней, plot-smoke, screen, FRI). Шаг «фигуры» в `dev/check.R`
      вызывает его **на каждой приёмке** и падает, если перегенерированная
      закоммиченная картинка (`man/figures`) отличается, — авто-синхрон c просьбой
      закоммитить, без тихого дрейфа. По умолчанию — только деревья решений (быстро,
      оффлайн); `--render` добавляет тяжёлые панели, `--no-render` пропускает шаг.
      Тяжёлые части (панель, FRI) — сеть/`OligoArrayAux`: на голой машине скип
      сообщением, не падение. `dev/render-decision.R` умеет `TIDYKASP_SKIP_EXAMPLE=1`
- [x] `.gitignore` (кэш, results, .tools)
- [x] Скрипты остаются скриптами, но с roxygen-заголовками — путь к реальному пакету (см. M)

### A1. Внешние тулы у пользователя: «поставил одной строкой и работает»  [частично]

Сейчас `hybrid-min`/`hybrid-ss-min` (OligoArrayAux — движок `DECIPHER`) и
`blastn`/`makeblastdb` (NCBI BLAST+) ищутся **только** через `Sys.which()`, то есть по
`PATH`. В репозитории это прячется за `data-raw/install-oligoarrayaux.R` + `.tools/` +
`tests/testthat/helper-oligoarrayaux.R`, но у пользователя, поставившего пакет
`remotes::install_github()`, тулов нет: термо/QC молча скипаются, BLAST-функции падают.
Нужен тот же путь, но без ручного `export PATH/UNAFOLDDAT`.

- [x] **Слоистый резолвер** (`R/kasp_tools.R`: `oligoarrayaux_paths()`, `blast_bin()`):
      порядок — `options(tidykasp.*)` → проектный `.tools/` → пользовательский кэш
      `kasp_cache_dir()/tools` → системный `PATH`. `has_oligoarrayaux()`/`has_blast()`
      и `oligo_dg()` ходят через него; DECIPHER/rBLAST оборачиваются в
      `with_oligoarrayaux()`/`with_blast()` (`withr`), так что `UNAFOLDDAT`/`PATH`
      действуют только на время вызова. Проверено: `pcr_efficiency()`/`oligo_dg()`
      работают без ручного `export`.
- [x] **Рантайм-установщик** `kasp_install_tools()`: качает релизный тарбол OligoArrayAux,
      сверяет sha256, собирает (`make`), дымовой тест, кладёт в **`R_user_dir`** (не в
      каталог пакета — он read-only), идемпотентно; согласие `ask = interactive()`;
      Windows — через `bash -lc`/Rtools. Проверено end-to-end (сборка в temp + ветка
      `options(tidykasp.oligoarrayaux=)`). Прецеденты: `tinytex::install_tinytex()`,
      `webshot2`/`chromote` (Chrome в рантайме), `rnaturalearth`/`elevatr` (данные).
      BLAST+ есть в apt/conda (`ncbi-blast+`) — **не бандлим**, а советуем `pak` (ниже).
- [x] **Одна строка для пользователя**: README (§0/§2) даёт `pak::pak("stathmin/tidykasp",
      sysreqs = TRUE)` (системные `ncbi-blast+`, cairo; `renv::sysreqs()` тоже умеет) →
      `tidykasp::kasp_install_tools()` → всё на месте. Первый вызов без тула не молчит:
      `assemble_kasp_assay()` подсказывает `kasp_install_tools()` (сейчас `rlang::inform`,
      `cli` — в §R), а не отдаёт тихий `NA`
- [x] **CRAN-рамки**: `SystemRequirements: OligoArrayAux (optional), NCBI BLAST+ (optional)`
      в `DESCRIPTION`; сеть и запись — только в рантайме, только в `R_user_dir`, с
      согласием; `configure` не качает; отсутствие тула — корректный skip/warning.
      Осталось закрепить словами в `CONTRIBUTING`.
- [x] Единый тест-хелпер `skip_if_no_oligoarrayaux()` вместо ad-hoc
      `skip_if_not(nzchar(Sys.which("hybrid-min")))` (24 замены в 12 файлах);
      `helper-oligoarrayaux.R` ищет и `.tools/`, и пользовательский кэш.

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
      цепей) как критерий среди топ-`specificity_candidates`. **Исправлено:** A/T и C/G
      больше не пропускаются — само-пары T·T/A·A — сильнейший 3'-мисматч; слабый — G:T/A:C
      wobble (транзиции), он ловится мягким `gate_wobble`, а не классом сайта
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
- [x] Футпринт AS — по dropout, а не по числу: `gate_as_footprint` требует
      `as_foot_dropout <= tol` (позиция+идентичность, `DECIPHER`, с весами),
      как COM-футпринт и AS 3'-окно; счёт `n_var == 0` остался лишь фолбэком
      без DECIPHER. Сцепленные варианты исключены (`as_exclude`). Проверено:
      редкий (0.01) SNP/индел внутри праймера **виден** (`AS footprint variants
      (n)` > 0, `as_foot_dropout` ~ 0) и **не думит**; стоимость растёт с частотой
      (`dev/haplotype-panel.R`, `test-frequencies.R`)
- [x] **MAF: у Ensembl их нет.** Проверено для Arabidopsis `variation/<id>?pops=1`:
      `MAF`/`minor_allele` = `null`, `populations` пуст (127 SNP) — 1001 Genomes
      отдаётся как источник без частот. Значит `pop_freq` — только явный/
      пользовательский; реальные MAF дал бы только сырой VCF 1001 Genomes (крупно)
- [x] **Модель гаплотипов** (`R/kasp_haplotypes.R`): вход
      `assemble_kasp_assay(haplotypes=)` — таблица (`haplotype`, `freq`, `chrom`,
      `pos`, `ref`, `alt`). Из неё выводятся (а) `pop_freq` на позицию (сумма
      частот носителей) → веса футпринта; (б) сцепленные позиции (едут с ALT-аллелем
      мишени) → `linked=`; (в) `haplotype_hits` — какой AS-праймер какой гаплотип
      ловит (`target_allele`, `primer`, `footprint_variants`, `binds`). Трек SNP
      рисуется **слоями** (`haplotype_layers()`): строка на гаплотип, `H 0.30 ->
      AS_alt`, свои варианты подписаны. Демо/трек: `dev/haplotype-panel.R`
- [x] Модель гаплотипов: множественные ALT-гаплотипы. `linked` теперь по
      гаплотипам (`target`, `haplotype`, `linked`, `shared`): каждому ALT-гаплотипу
      — свои сцепленные позиции, `shared` = едущие со *всеми* ALT (их и впекаем в
      один ALT-праймер, см. `assemble_kasp_assay()`); приватные остаются вариантами
      футпринта, что видно в `haplotype_hits` (`n_variants`, `footprint_variants`,
      `binds`). **Сцепленные инделы**: `design_as_primers(linked=)` пересобирает
      ALT-шаблон через `replaceAt` с учётом ширины (`ref`), а `v_start` сдвигается
      на суммарную длину правок выше блока — так делеция/вставка впекается без
      потери координат
- [x] Сцепление вариантов (co-segregation): вариант в 3'-окне AS (или футпринте COM)
      *вместе* с мишенью — не «грех», а союзник. Сделано: `screen_snp_targets(linked=)`
      прощает сцепленный вариант в 3'-окне; `assemble_kasp_assay(linked=)` /
      `design_as_primers(linked=)` «впечатывают» сцепленную замену в ALT-гаплотип и
      ALT-праймер (двойной мисматч к REF), REF не трогается; хранится в
      `assay$selection$linked`. Гаплотипы/LD — 1001 Genomes per-accession
- [x] Сцепление, остаток: (а) дизайн, держащийся *только* на сцеплении, теперь **помечается**,
      а не проходит молча: если сцеплённая позиция попадает в 3'-окно или футпринт AS,
      `assemble_kasp_assay()` предупреждает (`relies on linkage`) и кладёт
      `assay$selection$linked_dependent = TRUE` (покрыто `test-linked.R`); (б) сцепленный
      **indel** вправляется пересбором `window_alt` с шириной (разная длина `ref`/`alt`
      сдвигает координаты ниже блока на `delta`), а не хрупкой однобуквенной заменой;
      (в) множественные сцеплённые позиции обрабатываются одним вектором
- [x] Дешёвый векторный скрин мишени `screen_snp_targets()`: по SNP (или целому гену)
      без термодинамики/BLAST/геномных сканов — GC, run/SSR, хотзона вариантов, вариант
      в 3'-окне; вердикт «promising/doomed» и счёт `k/16` (strand x length). ~3 ms/SNP,
      фильтрует гены до `pick_target_snp()`/`assemble_kasp_assay()`/rBLAST; `plot_snp_screen()`
- [x] Мультиаллельные сайты и вырожденные праймеры: `kasp_site_class()` формализует
      критерий «сайт с >2 гаплотипами, всё ещё годными для одного KASP»: одиночный
      SNP-ALT (включая A/T, C/G — см. исправление выше) → `biallelic`; ≥2 SNP-ALT → `degenerate`
      (один AS-праймер с 3'-концевым IUPAC-кодом объединения ALT, `iupac_code()` из
      `Biostrings::IUPAC_CODE_MAP`, код исключает REF) покрывает все ALT одним KASP;
      не-SNP ALT → `indel` (свой путь). `screen_snp_targets()` получил колонку `class`,
      поэтому мультиаллельный сайт больше не помечается `doomed`. Фикстура разных типов
      (SNP, MNP, инсерция/делеция, два близких ALT) — `inst/extdata/example_targets.csv`,
      `kasp_example(type=)`: по одной мишени на тип маркера, разбросанных по генам и
      хромосомам, отобраны по *наименьшей
      местной нагрузке вариантов* (а не свалены из одного гипер-вариабельного региона);
      колонка `gene` — AGI-локус для карты фичи. SNP-строки проверены на дизайнабельность.
      Гаплотипы/LD — 1001 Genomes per-accession
- [x] **Мультиаллельный движок (дизайн):** `design_as_primers()` принимает вектор ALT:
      конфигурация оценивается по каждому ALT (DECIPHER), а сводится консервативно —
      on-эффективность минус максимум off-эффективности по ALT-гаплотипам, REF-праймер
      видит свой худший ALT-шаблон. Возвращаются двое: `alt_primer` — конкретный
      худший ALT-член (реальные ACGT, для термодинамики/специфичности/QC) и
      `alt_degen` — IUPAC-консенсус по всем ALT (физический вырожденный олиго).
      `assemble_kasp_assay()` для all-SNP мультиаллельного сайта гонит все ALT в движок,
      кладёт `as_alt` (конкретный, внутренний) и `as_alt_degen` (заказной); экспорт
      (`kasp_primer_table()`/`export_primers()`), `print()` и подписи на графике берут
      вырожденный; специфичность/термо/`find_palindromes` остаются на конкретном, поэтому
      IUPAC в них не протекает. Проверено: `test-as-primer.R` (движок + сборка, код S/M/K)
- [x] **Уродливые QC гейтят отбор:** приёмка ловит inverted repeat, второй продукт
      in-silico PCR и точные off-target хиты — раньше *выбор* их не отбрасывал (возвращался
      «ok»-ассей с `inv. repeat=2`, `PCR=2`, `off-target=4`). Теперь `assemble_kasp_assay()`
      после лестницы **экранирует** ведущие `pair_screen` (дефолт 3) пары по этим трём
      проверкам и берёт первую чистую; если ни одна не чиста — оставляет рангового лидера
      и предупреждает (как раньше). Хиты общих AS-праймеров сканируются один раз, только
      COM — по кандидату; ПЦР-пассы по паре. `pair_screen = 0` выключает. Покрыто
      `test-selection.R` (дубль сайта COM на второй хромосоме → свап на чистую сестру).
      Проверки остаются в `kasp_design_checks()`/`KASP_UNGATED_CHECKS` и в `kasp_risk()`
- [x] **RDS рядом с картинкой:** `dev/render-targets.R` сохраняет `target_<id>.rds`
      (полный ассей: праймеры, тумблеры, off-target, граф) рядом с PNG — разбор без пересчёта
- [x] **Проверки уровня сайта:** `kasp_design_checks()` получил уровень *site*
      (`target site class`, `ALT alleles (n)`, `linked variants skipped (n)`) —
      вырожденные/мультиаллельные сайты и несшитые indel'ы видны в тумблерах;
      `plot_checks()` фасетит site/feature/alignment; `assemble_kasp_assay()` кладёт
      в объект `assay$site` (`kasp_site_class`). Плюс два новых site-критерия:
      **«ALT primer covers all ALTs»** (вырожденный движок покрывает все ALT;
      `site$alt_covered`) и **«design relies on linkage»** (дизайн держится только на
      сцеплении) — видны и в тумблерах, и в дереве/тексте решения автоматически
      (реестр один)
- [x] **Дерево решения как данные:** `kasp_criteria()` (R/kasp_criteria.R) — единый реестр
      *всех* критериев (site-проверки, гейты и тиры отбора, приёмочные проверки). Из него
      генерируются узлы/рёбра (`kasp_decision_schema()`), граф (`kasp_decision_graph()`,
      `tidygraph`), текст (`kasp_decision_text()`) и проверки (`kasp_design_checks()` —
      `role=="check"`), поэтому новый критерий появляется и в дереве, и в результатах сразу.
      DAG сверху вниз: site → AS → ASxCOM → приёмочные проверки → verdict; жёсткие лимиты
      рисуются **жирной рамкой**, мягкие спускаются по лестнице релаксации и возвращаются в
      главный путь, есть терминал `DOOMED`. `assay$decision` — граф в объекте,
      `assay$evaluation` — вычисленные проверки (без пересчёта при отрисовке).
- [x] **Схемы по химии маркера:** `kasp_schemes()` (R/kasp_schemes.R) — реестр химий
      (family, дискриминация, ридаут, зиготность, число локусов). `kasp_criteria()` помечает
      каждую проверку своими схемами, поэтому один реестр даёт деревья для SNP/indel/
      мультиаллельного KASP, **дуплексных** химий (dKASP Jiang 2023, mxKASP Meilink 2025,
      Multi-KASP Wang 2025 — два локуса в одной пробирке + стадия совместимости: Tm-окно
      всех праймеров, запрет ≥6 bp inter-pair комплементарности, ≤5 self-comp, ампликон
      ≤150 bp) и **ASQ/PACE** (Kalendar 2022), плюс голый subset **«любая пара праймеров»**
      (маркер-агностичные проверки: состав/структура/повторы/off-target — не KASP).
      `kasp_criteria(scheme=)`, `kasp_decision_schema/graph/plot_decision_tree/decision_text`
      принимают scheme; стадия со всеми жёсткими гейтами ведёт провал прямо в `DOOMED`.
      `dev/render-decision.R` рисует по картинке на химию
- [~] **Движок под дуплексные химии:** `R/kasp_duplex.R` — `design_duplex(genome, loci,
      readout, homologous)` собирает каждый локус сам ([assemble_kasp_assay()]) и считает
      **стадию совместимости** ([duplex_compatibility()]) по всем праймерам обоих локусов:
      Tm-окно (59-65 °C), межпарная комплементарность (`complementary_run()` < 6 bp — правило
      мультиплекс-дизайна; внутри локуса уже ловит dimer-скрин), само-комплементарность
      (≤ 5 bp) и бюджет ампликона (≤ 150 bp); те же gate-id, что в реестре (`KASP_GATES_DUPLEX`).
      Ридаут co-dominant (dKASP: кассета на локус, `homologous = TRUE`) vs presence/absence
      (Multi-KASP: доминантный вызов из комбинации, `homologous = FALSE`), гейт `gate_duplex_targets`
      связывает одно с другим; класс `kasp_duplex` + печать. **Осталось:** (в) ASQ/PACE как
      отдельная cassette-химия
- [x] **Релаксация с границами, а не снятие гейта:** `select_with_ladder()` опускает порог
      мягкого гейта только до лучшего достижимого значения (не до `-Inf`) и не ниже floor;
      жёсткие гейты (дискриминация, вариант в AS 3'-окне) не опускаются вообще. Сцепленный
      вариант в AS 3'-окне не считается нарушением (исключён из счётчиков AS-зоны).
- [x] **Doomed-aware отбор мишени:** `pick_target_snp(scheme=, max_tries=)` перебирает
      ранжированных кандидатов и дизайнит их движком схемы (`kasp_design_at()` →
      KASP/CAPS/ARMS), возвращая первого *дизайнабельного* — «doomed» зависит от типа
      маркера, так что KASP-провал может пройти CAPS/ARMS. Дизайн едет в атрибуте `marker`;
      когда ничего не дизайнится — предупреждение и лучший по рангу (`31259d6`)
- [x] **Реальные indel/MNP:** окно длин AS-индел-праймера расширено `20:27` → `18:32`
      (FLC-инсерция 0 → 14 GC/Tm-годных, делеция → 92). Остаток (поли-A стык) честно
      недизайнабелен как KASP-индел → CAPS/dCAPS/HRM (§U). Сообщение `DOOMED` разнесено:
      GC/Tm и гомополимер/SSR — разные тексты

## I. Батч-дизайн
- [x] `design_kasp_markers()` — строка на локус, вердикт `pass`/`warn`/`fail` + причина
- [x] `variant_type()` — SNP / MNP / insertion / deletion
- [x] `BPPARAM` — параллель по локусам

## J. Визуализация (`R/kasp_plots.R`, `R/kasp_search_anim.R`, пакетные библиотеки)
- [x] `plot_primers()` — локальная карта ассая: буквы последовательностей (нейтрально
      чёрные), **жирные** — 3' GC-clamp, **красные** — основания SSR/шпильки; кольцо =
      мисматч/дискриминирующее 3'-основание; Tm числом + ΔTm-тумблер; дорожка вариантов
- [x] `plot_primers()` — при сцепленном гаплотипе трек вариантов расщепляется на два:
      **REF- и ALT-гаплотипы** (основания мишени и сцепленных позиций подписаны), так
      что сцепление видно глазом; `assay$selection$linked` несёт `ref` и `alt`
- [x] `plot_checks()` — тонкая колонка именованных тумблеров по критерию, фасет
      site/feature/alignment (правая панель `plot_primers()`)
- [x] `plot_decision_tree()` — граф политики отбора (`ggraph`/`tidygraph`),
      зелёный/янтарный/красный по статусу узла; данные графа — в объекте (см. §H)
- [x] **Контекстный вид** (крупный масштаб): `plot_variant_landscape()` переписан на
      `ggforce::facet_zoom`: модель гена (`gggenes`) + трек вариантов в полном виде и
      **зум-панель** на окно дизайна (рамка зоны на полном виде) — вместо ручной врезки
      `patchwork`. `ggforce` добавлен в `Imports`; `facet_zoom` строит без предупреждений
      (тест билдит объект и проверяет класс `FacetZoom`). `Gviz`/`karyoploteR` **сознательно
      не ставим**: `karyoploteR` тянет тяжёлое дерево Bioconductor (`biovizBase`/`regioneR`) и
      собирается из исходников ~20+ мин (бинарников для `linux-debian-forky` нет) — контекст уже
      даёт `facet_zoom`, а общегеномный идеограмм — `assay_ideogram()`
- [x] **Выравнивание на референс** с фичами (не общий MSA): `plot_primers()` — референс-бэкбон
      (буквы последовательности), ALT — только отличия (вариантный трек / слои REF- и
      ALT-гаплотипов), праймеры на его координатах, метки критериев рядом (`plot_checks()`)
- [x] Показ **self-complementarity, шпилек, SSR**: `plot_oligo_structure()` — три олиго рядом
      (буквы, SSR красной подчёркой, шпильки/инвертированные повторы дугами постоянной
      высоты над последовательностью; высота не растёт со стеблем, чтобы не доминировать —
      длина стебля читается по горизонтальному спану). `plot_self_complementarity()` —
      dot-plot self×self (i спаривается с j, когда base i — комплемент base j; шпилька =
      внедиагональный блок). Библиотечный ландшафт проверен: `ggforce::geom_arc` — это
      центр/радиус/углы (полуокружность высотой ∝ стеблю — плохо для компактной полосы),
      `geom_arc0` не существует; `R4RNA` рисует base-R arc-диаграммы (не ggplot-слой).
      Рисуется `geom_path()` по вычисленным точкам (нативный ggplot2)
- [x] `plot_discrimination()` — удалён (шумный бар-чарт)
- [x] `plot_variant_landscape()` — модель гена + врезка региона
- [x] `plot_alignment()`, `plot_self_complementarity()`, `plot_dimer_qc()` — **удалены**
      (в текущем виде не нужны; `ggmsa` убран из `Imports`). Замена выравнивания/
      самокомплементарности — отдельная задача ниже (dot-plot/дуги), пока не делаем.

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

### O1. Замер 2026: один ассей и батч (пул вместо повторов)  [сделано]

Профиль одного локуса (FRI `Chr4:269666 G>A`, 16 ядер): 43.4 → 29.9 с серийно и
→ 14.3 с (`options(tidykasp.cores = 8)`); батч 16 локусов 231 → 58 с
(`design_kasp_markers(BPPARAM = MulticoreParam(8))`). Вердикты и числа
идентичны серийным.

- [x] **Единый knob** `cores = getOption("tidykasp.cores", 1L)` → `kasp_bpparam()`
      (`NULL` при 1 ядре); `BPPARAM` задаёт бэкенд явно и главнее. Внутри внешнего
      locus-пула `cores` форсится в 1 — вложенности нет.
- [x] **Параллель подпроцессов**: карты OligoArrayAux (`oligo_dg_vec`,
      `terminal_dg`, `primer_thermo`) и батч DECIPHER (`pcr_efficiency(cores=)`,
      `design_as_primers*`).
- [x] **Пул dropout-пар**: `zone_dropouts()` — один вызов `pcr_efficiency` на все
      зоны ASS/COM вместо ~125; откат на построчный путь при ошибке.
- [x] **Один скан генома на ассей**: `assay_pcr_products()` переиспользуется
      счётчиками ампликонов и идеограммой; `primer_hits()` параллелен по хромосомам.
- [ ] **Пул между локусами** (следующий шаг, с оговоркой): профиль 2026 показал, что
      мисматч-скан сидов (`SEED_MISMATCH = 1`) **не** ускоряется библиотечным индексом
      (`PDict(max.mismatch=)` → `MTB_PDict`: `countPDict` считает верно, но в ~О(45)
      медленнее пер-сидового `countPattern`, 1.4 с против 0.03 с на 0.5 Мб), а
      остальные геномные сканы дёшевы (идеограмма/ПЦР ~2 с, `specificity` ~2 с в
      одном ассе). Для индексных геномов пул уже есть: `specificity_blast()` берёт
      вектор праймеров одним запуском. Поэтому полноценный пул оправдан только для
      батч-API, возвращающего объекты ассеев (сотни локусов → один проход генома на
      всё): сиды/праймеры всех кандидатов → один `PDict`/`matchProbePair`-проход.
      Оценка выигрыша без потери качества — ~10-20% для уже параллельного батча.

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

- [x] **Перевод-проход по `po`:** завершён — **167 сообщений переведено, 0 fuzzy,
      0 untranslated, 0 obsolete** (`msgfmt --statistics`). Сообщения новых схем (§U:
      sanger/dpcr/taqman/castpcr/snapshot/lamp/rpa/crispr/ssr/indel_size) и служебные
      переведены, `#, fuzzy` снят, `potools::po_extract()` → `po_compile()`
      (`msgfmt -c`, 167 сообщений), проверено `LANGUAGE=ru` end-to-end (en — английский).

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
- [x] **Сообщения: `cli` + стиль + i18n — одним проходом.** Перевести пользовательские
      ошибки/предупреждения с `rlang`-строк на `cli_abort()`/`cli_warn()`/`cli_inform()`
      с bullet-ами (`i`/`x`) и плюрализацией (`cli::qty()`/`{?ед/мн}`), заодно приведя
      стиль к Q1 «строчная, без точки» (§T). Смена msgid рвёт переводы, поэтому делать
      **вместе с `po/R-ru.po`**: `potools::po_extract()` → обновить `.po` (элементы
      `cli`-векторов — по элементу) → `potools::po_compile()`, затем прогнать i18n
      end-to-end (`LANGUAGE=ru`). Диагностики `check_cracked_messages()`,
      `check_untranslated_cat()`, `check_untranslated_src()`. Начать с самого
      user-facing: валидация аргументов и подсказки установки тулов (`kasp_install_tools()`).
      **Статус:** каталог полон и синхронизирован. **Найдено:** миграция на `cli`
      **msgid-сохраняющая** — достаточно заменить `rlang::abort/warn/inform` на
      `cli::cli_abort/cli_warn/cli_inform` и убрать обёртку `glue::glue()` (cli сам
      интерполирует `{}`), **не трогая текст внутри `gettext()`**: `po` при этом не меняется
      (проверено — `po_extract`, msgid идентичны; `po/R-ru.po` не изменился). **Сделано:**
      тулы/установка + валидация ввода (`R/kasp_tools.R`, `R/kasp_io.R` — 11 сообщений).
      **Осталось:** ничего — миграция завершена тем же механическим приёмом (все
      `rlang::abort/warn/inform` → `cli::cli_abort/cli_warn/cli_inform`); `po` не трогался.
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
- [ ] **Windows: сборка OligoArrayAux + R-hub.** Задокументировать сборку на Windows
      (Rtools/MSYS2/UCRT64, запуск через `bash -lc`; см. установщик A1) и проверить на
      живом Windows/R-hub: `kasp_install_tools()` под Windows, паритет `-q`-вывода
      `hybrid-min.exe`, `PATHEXT`/processx, cairo/`ggmsa`, CRLF у `vroom`. **Нужен
      Windows-раннер** — в текущей (Linux) среде не воспроизводится.
- [x] **sha256 и пробелы:** `sha256_of()` пробует `sha256sum`/`shasum`/`openssl`;
      `UNAFOLDDAT`/`--prefix` кавычатся; компилятор ищется как `cc`/`gcc`/`clang`
- [x] **`.gitattributes`:** `* text=auto eol=lf`, `*.mo`/`*.rds`/`*.png`/`*.gif` — бинарные
- [x] **Инструкции/шелл:** README и `lecture/render.sh` — варианты cairo для
      Windows/macOS (pak `sysreqs` / brew / Rtools), запуск `sh` на Windows из
      Git-Bash/Rtools; `apt-get`-хинты в сообщениях (`R/kasp_scan.R`,
      `R/kasp_search_anim.R`) дополнены macOS/Windows
- [x] **Форк:** `render_design_search()` больше не зовёт `parallel::mclapply` (Unix-only,
      и `parallel` не был в `Imports`); теперь `BiocParallel` — `MulticoreParam` (форк) на
      Unix, `SnowParam` иначе, через общий `bp_lapply()`. Проверено на Linux (кадры в
      параллель, GIF пишется).

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
- [~] **Стиль сообщений** против Q1 (заглавная/точка, напр. `kasp_io.R`,
      `kasp_scan.R`): механизм доставлен — все сообщения идут через `cli` (§R), — но
      перерегистр ~70 текстов нашёл бы ~70 переводов *fuzzy* ради косметики, поэтому тексты
      намеренно оставлены в текущем виде (полные предложения с заглавной и точкой). base R
      пишет строчными без точки; расхождение осознанное — сохранённые переводы важнее
- [x] **i18n / `Language:`** — решено **не** задавать: справочник (`.Rd`) английский,
      и `Language: ru` объявил бы его русским; русская проза — отдельный слой (Q4)
- [x] **Тесты не ассертят переводной текст** напрямую: в блоках с англ. подстроками
      добавлен `withr::local_envvar(LANGUAGE = "en")`
- [ ] **Бейджи/CI:** бейдж `lifecycle: experimental` без зависимости `lifecycle`;
      бейдж R-CMD-check статичный, воркфлоу нет (CI намеренно не в репозиторий, см. CONTRIBUTING)
- [x] **Мелочи:** `COPYRIGHT HOLDER` в `LICENSE` = автор; добавлен `inst/CITATION`;
      `Config/roxygen2/version` совпадает с установленным (8.1.0); убран лишний
      `utils::globalVariables(".data")` при `@importFrom rlang .data`
- [x] Уже соответствует: Title/Description (англ., title case), `Authors@R`,
      License+LICENSE, URL/BugReports, `Encoding`, `|>`/`.data$`, `rlang`/`cli`,
      `@family`/`@noRd`, тесты edition 3 без сети/снапшотов, поставка `po/` + `.mo`

## U. Другие PCR-маркеры и внешние приоры (2025–2026)

Три «капризных» локуса из фикстуры (вариант в AS 3'-окне; поли-А стык; `disc = 0.497`)
— не провал, а сигнал: KASP — не единственная PCR-химия. У KASP три врождённых
ограничения: (1) аллель различается **3'-основанием** праймера — ампликон < 100 bp,
`Tm_COM > Tm_AS`, ~45% GC у варианта, праймер в 20–30 bp от варианта (EasyKASP 2025);
(2) различение — это **разница эффективностей** (порог `min_disc`), не бинарный признак;
(3) вариант якорится в **праймере**, а не в ампликоне.

Карта «почему локус провалился → чем его взять»:

| Провал (наш локус) | Причина | Чем взять |
|---|---|---|
| вариант в 3'-окне (MNP GA>AA) | нужен чистый 3'-торец | **CAPS/dCAPS/PIRA** (вариант внутри ампликона, различение — рестрикция), **HRM**, ампликонное секвенирование |
| AT-rich / поли-А стык (A>AA) | нет 20–27mer с `Tm ≥ 56` / `GC ≥ 0.35` | **CAPS/dCAPS** (праймеры в GC-богатых флангах), фрагментный анализ (SSR/indel), секвенирование; KASP — только с расширенным окном Tm/GC |
| `disc < 0.5` (TT>T) | слабая разница эффективностей | **CAPS/dCAPS** (бинарный разрез), **HRM** (короткий ампликон), **tetra-primer ARMS** (полосы по размеру), секвенирование |

Общий вывод: **CAPS/dCAPS (и PIRA-PCR)** — рабочая лошадка там, где вариант — плохой
якорь для праймера: он якорится в **ампликоне**, а различение — бинарная рестрикция.
**HRM** терпит соседей и слабую дискриминацию на коротком ампликоне (слаб в AT-богатом).
**Tetra-primer ARMS** даёт гель-ридаут (две аллельные полосы + контроль), но разделяет
те же 3'-ограничения KASP. **Ампликонное/Sanger-секвенирование** — универсальный фолбэк
(gold standard в сравнении 2025 г.).

- [~] **CAPS/dCAPS как вторая схема.** Прототип есть: `R/kasp_caps.R` —
      `restriction_enzymes()` (**данные из REBASE *bionet*, ~1200 палиндромных Type II**,
      не руками; регенерация `dev/fetch-rebase.R`), `caps_fragments()` (IUPAC-aware сайты →
      длины фрагментов),
      `find_caps_enzymes()` (REF vs ALT-ампликон, «gel resolvability» по разрыву полос,
      дефолт 50 bp) и `design_caps_marker()` (выбор фермента + фланкирующая пара с Tm/GC-скрином,
      возвращает класс `kasp_caps` с рецептом рестрикции). Схема `caps` зарегистрирована в
      `kasp_schemes()`/`kasp_criteria()`/`kasp_decision_stages()` (своя стадия `gate_caps_*`,
      маркер-агностичные праймерные проверки, без KASP 3'-машинерии); дерево рисуется.
      **Сделано:** (а) **dCAPS** — сайт, *создаваемый* праймером: `find_dcaps_enzymes()`
      ищет окна вокруг варианта, которые на один (не-вариативный) мисматч короче сайта для
      одного аллеля, применяет базу праймера к обоим ампликонам и режет;
      `design_dcaps_marker()` — фолбэк `design_caps_marker()`. **Сделано также:**
      (б-гейт) `site_change` — CAPS-вызов требует смены *сайта*, а не тривиальной разницы
      длин индела (иначе 1-bp индел «различался» бы всеми ферментами с gap=1);
      (в) виртуальный гель `plot_caps_gel()`; (б) каталог расширен: REBASE bionet, ~1200
      палиндромных Type II (`dev/fetch-rebase.R`; непалиндромные и Type IIS исключены —
      фрагментная модель берёт один топ-страндовый разрез). **Осталось:** (г) `@examples`/виньетка.
      FOSS-ориентиры: dCAPS Finder 2.0 (Neff 2002), indCAPS (Kieber lab), SNP2CAPS,
      primerblast-oss. **Эмпирика:** на 3 «капризных» локусах (MNP GA>AA, insertion A>AA,
      deletion TT>T) природного CAPS-маркера в каталоге нет (0 сайтов; у длинных делеций/
      инсерций gap=1) — им нужен dCAPS или больший каталог ферментов
- [~] **Tetra-primer ARMS / AS-PCR.** `scheme = "arms"` реализован (`R/kasp_arms.R`):
      `design_arms_marker()` — внешняя контрольная пара (`design_flanking_pair()`) + две
      внутренние аллель-специфичные на противоположных цепях (`design_as_primers()`, тот же
      3'-концевой якорь), дающие полосы REF/ALT *разных* размеров (AA/AB/BB); валидация
      разрыва полос, класс `kasp_arms`, `plot_arms_gel()`. Схема зарегистрирована (стадия
      `gate_arms_*`). **Осталось:** усиление мисматчей -2/-3 (как ARMSprimer3), учёт
      off-target/специфичности; FOSS: Primer1 (Ye 2001/2012), WASP (2007), ARMSprimer3 (2025),
      primerblast-oss
- [x] **HRM** (`scheme = "hrm"`, `R/kasp_hrm.R`): `design_hrm_marker()` — одна фланкирующая
      пара поверх варианта (вариант в домене плавления, под праймером не сидит), ампликон
      **короткий** (40–150 bp; короткий домен — резче кривую, Premier Biosoft / Thermo Fisher),
      пара Tm-сбалансирована; соседи терпимы (нет гейта 3'-окна) — поэтому HRM берёт
      локусы, которые KASP думит (вариант в AS 3'-окне, поли-A стык). Схема
      зарегистрирована (`kasp_schemes`/`kasp_criteria`/`kasp_decision_stages`, стадия `hrm`),
      класс `kasp_hrm`, диспетчер `kasp_design_at("hrm")`. **Честно:** кривая плавления *не*
      моделируется (нужны uMELT/PolyMelt-модели) — возвращается корректный ампликон, не профиль.
      FOSS-ориентиры: uMELT/PolyMelt, SNPgen (2024)
- [x] **«Ладдер» химий — KASP не единственный путь:** `design_marker(genome, snp, schemes)`
      перебирает химии по очереди (KASP → CAPS → ARMS → HRM) и возвращает первую годную
      (диспатч `kasp_design_at()` + `is_target_designable()`), отдавая tidy-вердикт по
      каждой (`tried`). `dev/render-targets.R` теперь гоняет весь ладдер на каждой мишени
      панели: где KASP думит, рисуется альтернатива (`plot_caps_gel`/`plot_arms_gel`/
      `plot_hrm_amplicon`). **На курированном наборе (9 локусов):** 7 дизайнятся KASP;
      MNP `ENSVATH15405655` спасает HRM; поли-A инсерция `ENSVATH09061610` и делеция
      `ENSVATH15371499` — **ни одной химией** (CAPS/ARMS только для замен, а в AT-богатом
      контексте нет HRM-пары) → честно `NO marker` (нужны Sanger/dPCR, «дальние методы»)
- [x] **PIRA-PCR** — частный случай CAPS (интро-праймер создаёт сайт): покрыт
      `find_dcaps_enzymes()` (тот же приём — сайт вносит праймер)
- [ ] **Дальние методы** — держать в уме под задачи, а не в ядро: SNaPshot / single-base
      extension, ампликонное/Sanger-секвенирование (фолбэк + аттестация), digital PCR
      (низкое обилие), LAMP/RPA + CRISPR (SHIELD/SHERLOCK/DETECTR) — изотермические полевые
- [ ] **Внешние FOSS для параметров/методов/приоров:**
      - **Primer3 / primer3-py** — NN-термодинамика (SantaLucia), соль-коррекция,
        hairpin/dimer; независимая сверка наших ΔG (у нас TmCalculator + OligoArrayAux).
        `Primer3_masker` — маскирование
      - **primerblast-oss (2026, MIT)** — офлайн-аналог Primer-BLAST: 3'-end-aware модель
        прайминга (**3'-окно 5 bp, ≤1 мисматч в нём, ≤4 всего, терминальное основание
        совпадает**), термодинамический гейт (**anneal Tm ≥ 40, 3'-ΔG ≥ -5**), dimer-приоры
        (**ΔG ≥ -8 / Tm ≥ 45**), **gel resolvability** (разрыв ≥ 50 bp), CAPS/dCAPS +
        AS-PCR/tetra-ARMS + **multiplex-set design** (один совместимый праймер на мишень),
        risk-скор. Наши пороги (homodimer/duplex ≥ -9, `clean_3p = 5`) — рядом; выровнять
      - **PrimerServer2** — strand-aware pairing BLAST-хитов, `--use-3-end`, мультипоточность
      - **polyMarker (2015)** — KASP для полиплоидов (локальные выравнивания)
      - **WASP (2007)** — batch AS-PCR с намеренными мисматчами (у нас `stab_ats`)
      - **FastPCR / PrimerDigital** — suite (KASP/PACE/ASQ, CAPS, STR, SNaPshot); платная,
        но пороги читаемы; **KASPioneer (2025, PyPI)** — свежий KASP-дизайнер, сверить приоры
      - **KASPspoon (2019)** — in-silico PCR для KASP
      - Обзорные сравнения: *Anal. Chim. Acta* 2025 (rs9939609 T>A: ARMS vs PIRA vs TaqMan
        vs CADMA+HRM vs HRM-snapback, Sanger gold) — карта «какой метод на какой SNP»;
        **KASP vs T-ARMS** (MDPI 2026) — флуоресцентный против гель-ридаута
- [x] **Приоры в код:** `kasp_priors()` — единый провенанс порогов (значение → источник:
      LGC/EasyKASP 2025, primerblast-oss, TmCalculator), сверка термодинамики —
      `thermo_crosscheck()`. **`min_disc` градирован, а не бинарен:** `discrimination_band()`
      (R/kasp_risk.R) отображает наблюдаемый margin в low/medium/high, привязанный к
      `min_disc` + `discrimination_margin` (0.1) — без подогнанной вероятности, для которой
      нет данных о калибровке; `kasp_risk()` возвращает `discrimination` и `disc_risk`.
      Гейт отбора остаётся жёстким («no allele discrimination = no assay»), градируется
      *вердикт*, а не приём
- [~] **Провенанс порогов + термодинамическая сверка (без новых зависимостей).**
      `kasp_priors()` — таблица «параметр → значение → источник» (LGC/EasyKASP 2025 /
      primerblast-oss), аудит для `kasp_design_checks()`. `thermo_crosscheck()` прогоняет
      олиго через несколько публичных NN-наборов TmCalculator (SantaLucia 2004 — наш
      дефолт — vs Allawi 1998 / Sugimoto 1996 / Breslauer 1986) и отдаёт разброс Tm:
      подтверждает, что окно `as_tm` (10 °C) шире межмодельного разброса. `primer_tm(nn_table=)`.
      `dev/compare-thermo.R` пишет `results/thermo-reconcile.csv`. **Замечание:** отдельного
      бинаря/пакета `primer3` в среде нет (README:63); наша Tm — это и есть NN-модель
      семейства Primer3/oligotm, так что «independent ΔG cross-check» — это NN-таблицы,
      а ΔG остаётся за OligoArrayAux. **Сделано:** `rprimer` подключён (Bioc, Suggests) —
      `R/kasp_rprimer.R`: `alignment_from()` строит выравнивание, `rprimer_assays()` -
      независимый дизайнер + `checkMatch()` (on/off-target пропорции). Честное ограничение:
      S4-классы rprimer заперты, *свои* праймеры в `checkMatch` не вставить (валидируется
      локус/выравнивание, не готовый ассей). **Сделано также:** `kasp_risk()` —
      experimenter-facing риск `low/medium/high` (grave-провал → high) с причинами
      (R/kasp_risk.R). **Сделано:** `min_disc` градирован (`discrimination_band()` +
      `disc_risk` в `kasp_risk()`), а не бинарен — см. «Приоры в код» выше

## V. Модели считывания (readout simulation)  [очередь]

Дизайнеры химий (§U) возвращают **конструкты** — праймеры, пробы, гайды, ампликоны — но
**не предсказанный сигнал**. Это осознанное ограничение, но оно дорого: лабораторный
(«мокрый») вывод остаётся за человеком, а обещание пакета — «tidy-дизайн и QC» — без модели
считывания неполно. Добавлено в очередь, по химии:
- [ ] **HRM** — модель плавления (uMELT/PolyMelt-класс): ΔTm гомо-/гетеродуплекса, форма
      кривой; только она делает «короткий ампликон» осмысленным критерием
- [ ] **dPCR** — пуассоновское разбиение капель: цифровой подсчёт, доверительный интервал,
      предел обнаружения из числа капель и целевой концентрации
- [ ] **TaqMan / castPCR** — гибридизация + 5'-нуклеазная (Taq-дотяжка) эффективность по
      аллелям: у нас уже есть `DECIPHER::CalculateEfficiencyPCR` для праймеров — расширить на пробы
- [ ] **CRISPR-Dx** — коллатеральное расщепление Cas12/13: качественная модель «есть/нет»
      сигнал-репортера (не энергия, а порог)
- [ ] **ARMS / CAPS / indel_size / SSR** — виртуальный гель: разрыв полос по размеру
      (для CAPS/ARMS уже есть `plot_caps_gel`/`plot_arms_gel` — обобщить на indel/SSR)
- [ ] **LAMP / RPA** — изотермическая кинетика (качественно: порог/время до порога)
- [ ] Явно помечать в объекте/дереве, что readout **не** смоделирован (флаг-допущение),
      чтобы `design_*` не выглядел обещанием кинетики

## W. Повторяющийся аудит (сдвиг «сухих»/«мокрых» принципов, доков, стиля)  [очередь]

Наблюдение: почти все задачи закрыты — и это тревожно: значит, аудитов было мало, и дрейф
копится **между** проходами (новые модули, закрытые пункты, изменённое поведение). Нужен
регулярный аудит (после каждой крупной группы правок), а не разовый:
- [ ] **DRY / «сухой» код против библиотек:** нет дублей-копипаст (новые модули `kasp_*.R`
      не должны копировать CAPS/ARMS/`primer_candidates()`); где есть библиотека — берём её,
      а не самопис; сверить, что discovery/registry не обошёл библиотечные примитивы
- [ ] **Дрейф документации (ROADMAP/README/NEWS/виньетки ↔ код):** каждый закрытый пункт
      перечитать против кода; `@return`/`@examples`/`@family` у **всех** экспортов (включая
      новые); винетки не обещают снятое/изменённое
- [ ] **Стиль/линт/repo-требования:** styler, lintr, `R CMD check --as-cran` (0 WARN/NOTE),
      NAMESPACE ↔ man/ sync, `renv.lock` без мёртвого веса, `Rbuildignore`
- [ ] **«дерево == код»:** гейты/тиры/пороги в дереве равны константам селектора для **каждой**
      схемы (проверено скриптом), плюс полнота `gate_labels` (нет `NA`)
- [x] **i18n-паритет:** 0 fuzzy/untranslated/obsolete (см. §Q) — 167 переведённых
- [ ] **«мокрые» принципы:** биофиз. допущения каждой химии задокументированы и не выдаются
      за симуляцию (см. §V); пороги сверить с литературой 2025–2026
- [ ] Встроить аудит в `dev/check.R` как регулярный шаг (рутина, а не разовый заход)

### Аудит 2026-10 (4 агента: QC-дыры, API, tidy/CRAN, биология схем)  [исправлено + очередь]

**Исправлено:**
- **Совместная проверка пула.** `tube_compatibility()` (`R/kasp_pool.R`) — один хелпер
  над `primer_set_compatibility()`; теперь финальный набор олиго (праймеры **и** пробы)
  проверяется в `lamp`/`arms`/`snapshot`/`taqman`/`castpcr`/`dpcr`/`rpa` (раньше — только
  `lax`/`duplex`). Модифицированные основания (THF-сайт пробы) трактуются как `N`.
- **dPCR-проба сделана аллель-специфичной** (`R/kasp_dpcr.R`): две пробы (REF/ALT,
  FAM/HEX) — вариант-основание меняется во второй пробе (strand-aware); раньше одна
  общая spanning-проба давала лишь присутствие/отсутствие, а схема заявляла co-dominant.
- **`po`-каталог вычищен** до 168 переведённых, 0 fuzzy/untranslated/obsolete (§Q).
- **Ладдер:** `per_locus_schemes()` ограничивает дефолтный `rescue_strategy()` схемами
  с одно-локусным дизайнером — `duplex` (`loci = 2`) и `asq` (без движка) больше не висят
  в ладдере как «непроходимые»; дрейф доков `design_marker`/`design_markers` снят.
- `@param strategy`; битая ссылка `[kasp_scheme_registry()]`; сталые «wired by the caller».
- `options(tidykasp.tool_hint_shown=)` → `.kasp_state` (`R/zzz.R`, namespace-env) — пакет
  больше не пишет глобальную опцию; `<<-` в CRISPR-Dx заменён на императивный probing-loop.
- Тесты: новый `test-pool.R` (хелперы) + проверка `$compatibility` в схемах.

**Осталось (очередь):**
- **probe-Tm порядок TaqMan/castPCR** — вендорское «+8..10 °C» **недостижимо** на нашей
  сжатой NN-шкале (замер: 30-мер пробы даёт лишь ≈ +4 °C над праймерами); сейчас порядок —
  *ранг* (tier2), не гейт. Нужна вендор-шкала Tm или принудительно длинные/горячие пробы.
- **CRISPR/Cas13** — гайд выбирается без ограничения по транскрибируемой (смысловой) РНК;
  Cas13 режет ssRNA, поэтому гайд вне T7-цепи инертен.
- **LAMP** — не проверяется порядок Tm F1c/B1c > F2/B2 и длина FIP/BIP (45–49 nt);
  внутренние плечи по умолчанию 18:22 дают 36–44 nt.
- **indel_size** — `min_gap = 5` bp агарозой не разрешается (≥10 или PAGE/капилляр);
  заголовок это признаёт, дефолт — нет.
- **Дубли:** probe-candidates (`taqman` ≡ `castpcr` ≡ `rpa`), homopolymer/SSR-regex (4 места),
  guard «single-base only» (~10 дизайнеров; **меняет msgid** → один проход с §Q/§T).
- **`assemble_kasp_assay()` 609 строк** — вынести этапы в приватные хелперы.
- `render_design_search()` экспортирован, без теста, течёт temp-PNG; `Sys.setenv` в
  `helper-oligoarrayaux.R` — заменить на `withr::local_envvar`.

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
- Другие PCR-маркеры (см. §U): dCAPS Finder 2.0 (Neff 2002); indCAPS (Kieber lab);
  SNP2CAPS (Thiel 2004); Primer1 (Ye 2001; Collins 2012) и ARMSprimer3 (2025) для tetra-primer
  ARMS; WASP (Wangkumhang 2007); polyMarker (Ramirez-Gonzalez 2015); KASPspoon (Alsamman 2019);
  EasyKASP (Zhang 2025); KASPioneer (2025); FastPCR / PrimerDigital.
- FOSS-ориентиры по параметрам: **Primer3** (Untergasser 2012) и primer3-py; **PrimerServer2**;
  **primerblast-oss** (2026, MIT) — 3'-end-aware прайминг (окно 5 bp, ≤1 мисматч, ≤4 всего),
  термогейт (anneal Tm ≥ 40, 3'-ΔG ≥ -5), dimer (ΔG ≥ -8 / Tm ≥ 45), gel resolvability ≥ 50 bp.
- Сравнения методов: *Anal. Chim. Acta* 2025, doi:10.1016/j.aca.2025.343994 (rs9939609 T>A:
  ARMS / PIRA / TaqMan / CADMA+HRM / HRM-snapback, Sanger gold); KASP vs T-ARMS (MDPI 2026);
  SNPgen (2024) — веб-портал ARMS + HRM.
