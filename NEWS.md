# tidykasp

Демонстрационный пайплайн подбора KASP-маркеров для *Arabidopsis thaliana*
на чистом R. Здесь — краткая история значимых изменений; подробный инвентарь и
очередь работ — в `ROADMAP.md`.

## 0.3.0

### Производительность: пул вместо повторов, параллель по ядрам
- Один knob `cores = getOption("tidykasp.cores", 1L)` (форк на Unix, SOCK иначе);
  `BPPARAM` по-прежнему задаёт бэкенд явно. `assemble_kasp_assay(cores=)`
  распараллеливает карты OligoArrayAux (`terminal_dg`, `primer_thermo`) и батч
  DECIPHER; `design_kasp_markers(BPPARAM=)` форсит `cores = 1` во внешнем пуле —
  без вложенного форка.
- `zone_dropouts()` собирает (праймер, мишень)-пары всех кандидатов зон AS/COM в
  **один** вызов `pcr_efficiency` (было ~125 вызовов DECIPHER на ассей).
- Виртуальная ПЦР считается один раз (`assay_pcr_products()`) и переиспользуется
  счётчиками ампликонов и идеограммой; `primer_hits()` параллелится по хромосомам.
- Замер (FRI `Chr4:269666 G>A`): 43.4 → 29.9 с серийно, → 14.3 с при `cores = 8`;
  батч 16 локусов 231 → 58 с (`MulticoreParam(8)`); числа идентичны серийным.

### Отбор пар: гейт серьёзного, скор минорного
- `assemble_kasp_assay()` выбирает пары не взвешенной суммой, а жёсткими гейтами
  (серьёзные грехи: дискриминация, структура, `Tm_COM >= Tm_AS`, ампликон, вариант под
  COM-футпринтом, вариант в 3'-окне AS) и лексикографическим рангом по марже тира
  (узкое место, `min`). Пустое множество — «лестница» релаксации от мягкого к серьёзному
  с предупреждением на каждый нарушенный гейт; результат пишется в `assay$selection`.
- Дизайн, держащийся **только на сцеплении**, больше не проходит молча: если сцеплённый
  вариант попадает в 3'-окно или футпринт AS, ассей предупреждает (`relies on linkage`) и
  кладёт `selection$linked_dependent = TRUE` — валиден только там, где есть этот гаплотип.
- Сообщения `DOOMED` для AS-праймера разнесены: GC/Tm и гомополимер/SSR — разные тексты
  (раньше путались); каталог `ru` синхронизирован (129 сообщений).
- `kasp_design_checks()`: новые критерии «AS 3' window variants (n)» и
  «COM footprint variants (n)»; тумблеры в `plot_checks()` подписаны соответственно.
- Необязательный `pop_freq` в таблице вариантов: взвешивает нагрузку футпринта и
  масштабирует размер маркеров на трэке вариантов (иначе вес 1).
- Ожидаемая доля «выпадения» под футпринтом праймера (`DECIPHER`): не голый счёт
  вариантов, а частотно-взвешенная потеря эффективности, где положение и тип
  мисматча учитывает 3'-концевая дотация Taq. Порог `KASP_LIMITS$footprint_dropout`;
  гейты COM-футпринта и AS 3'-окна; тумблеры `AS drop` / `COM drop`.

### Дешёвый скрин мишени
- `screen_snp_targets()` — векторный, без термодинамики/BLAST/геномных сканов: по SNP
  (или целому гену) считает, сколько канонических посадок AS (strand x length)
  проходят GC/run/SSR/хотзону и чистое 3'-окно, и выдаёт `promising`/`doomed` со счётом
  `k/16` (~3 ms/SNP). Фильтрует гены до `pick_target_snp()`/`assemble_kasp_assay()`/rBLAST.
- `plot_snp_screen()` — по-вариантный профиль счёта по гену.
- `linked=` в скрине прощает сцепленный (co-segregating) вариант в 3'-окне.
- `kasp_site_class()` формализует мультиаллельные сайты: одиночный SNP-ALT → `biallelic`
  (любой, включая A/T и C/G — см. исправление ниже), ≥2 SNP-ALT → `degenerate` (один
  AS-праймер с 3'-концевым IUPAC-кодом объединения ALT, `iupac_code()`, исключает REF,
  покрывает все ALT одним KASP); не-SNP ALT → `indel`.
  Скрин получил колонку `class`, поэтому мультиаллельный сайт больше не помечается `doomed`.

### Исправления ядра
- **A/T и C/G больше не считаются неразделимыми.** `kasp_site_class()` помечал одиночный
  комплементарный ALT (A/T, C/G) как `complementary` с `disc = FALSE`, и `pick_target_snp()`
  пропускал такие сайты — это ~1/3 SNPs, и научно наоборот. На самом деле транзверсии A/T,
  C/G дают *само-пару* (T·T, A·A) на чужой матрице — сильнейший 3'-мисматч; слабый случай —
  G:T/A:C *wobble* (транзиции A/G, C/T), и его ловит мягкий `gate_wobble`. Класс
  `complementary` убран, `pick_target_snp()`/`screen_snp_targets()` больше не пропускают A/T/C/G.
- **Один сайт — один вариант.** `count_variants_in()`/`variant_load_in()` теперь дедуплицируют
  позиции, а `zone_dropout()` — пары (позиция, ALT): мультиаллельная VCF-строка (одна строка на
  ALT) больше не считается несколько раз в `n_var`/`n_3p`/нагрузке футпринта.
- **Окно Tm AS-праймера согласовано с приёмкой:** `assemble_kasp_assay(tm_range=)` по
  умолчанию берёт `KASP_LIMITS$as_tm` (было жёстко `c(56,65)`, уже приёмочного `[55,65]` —
  годные кандидаты с Tm 55–56 отбраковывались на дизайне).
- **Мультиаллельный сайт наконец виден:** `assemble_kasp_assay()` разворачивает list-колонку
  `ref_genome` (`unlist` вместо `as.character` — последний строкифицировал список, и сайт с
  двумя ALT помечался `biallelic`); теперь класс `degenerate`, `ALT n=2`, а `plot_primers()`
  подписывает мишень как `A>T/C` и в треке SNP, и в заголовке. ALT-праймер при этом строится
  по первому ALT, и об этом выдаётся **warning** (вырожденный праймер — в roadmap).
- **RDS рядом с картинкой:** `dev/render-targets.R` сохраняет `target_<id>.rds` (полный ассей)
  рядом с PNG.

### Другие PCR-маркеры и провенанс порогов
- **CAPS/dCAPS-схема** (`R/kasp_caps.R`): `restriction_enzymes()` (данные из **REBASE
  bionet**, ~1200 палиндромных Type II, регенерация `dev/fetch-rebase.R` — не руками),
  `caps_fragments()` (IUPAC-aware сайты → длины
  фрагментов), `find_caps_enzymes()` (REF vs ALT-ампликон, gel-resolvability по разрыву полос),
  `find_dcaps_enzymes()` (**dCAPS** — сайт, который доводит праймер: окно на один мисматч
  короче сайта, база вносится праймером, вариант внутри сайта) и `design_caps_marker()`
  (природный фермент или dCAPS-фолбэк + фланкирующая пара, класс `kasp_caps` с рецептом
  рестрикции). Схема `caps` в `kasp_schemes()`/`kasp_criteria()`/`kasp_decision_stages()`:
  своя стадия `gate_caps_*`, маркер-агностичные праймерные проверки, без KASP 3'-машинерии.
- **Независимый валидатор rprimer** (`R/kasp_rprimer.R`, Bioc, `Suggests`): `alignment_from()`
  строит выравнивание (REF/ALT, гаплотипы), `rprimer_assays()` запускает независимый дизайнер
  + `checkMatch()` (on/off-target пропорции) — второй алгоритм о том же локусе.
- **Прозрачность «уродливых» проверок:** `assemble_kasp_assay()` теперь **предупреждает**
  именем проваленной структурно/специфичной проверки (inverted repeat / in-silico PCR /
  off-target — `KASP_UNGATED_CHECKS`), а `kasp_risk()` сворачивает все проверки в
  `low`/`medium`/`high` с причинами (grave-провал → high).
- **CAPS:** `site_change` (CAPS — это смена *сайта*, а не тривиальная разница длин индела)
  и виртуальный гель `plot_caps_gel()`.
- **Tetra-primer ARMS** (`scheme = "arms"`, `R/kasp_arms.R`): `design_arms_marker()` —
  внешняя контрольная пара + две внутренние аллель-специфичные на разных цепях, полосы
  REF/ALT разного размера (AA/AB/BB); `plot_arms_gel()`. Схема в реестре (стадия `gate_arms_*`).
- **Провенанс порогов и термодинамическая сверка** (`R/kasp_priors.R`): `kasp_priors()` —
  таблица «параметр → значение → источник» (LGC/EasyKASP 2025 / primerblast-oss);
  `thermo_crosscheck()` — Tm одного набора олиго в нескольких NN-таблицах TmCalculator
  (SantaLucia 2004 vs Allawi/Sugimoto/Breslauer) с разбросом; `primer_tm(nn_table=)`.
  `dev/compare-thermo.R` пишет `results/thermo-reconcile.csv`.

### Гигиена дерева
- Внешние данные и индексы (геном, варианты, BLAST/FM) перенесены в пользовательский
  кэш `kasp_cache_dir()` (`tools::R_user_dir`), а не в проект; `data/` из проекта убран.
- Лекционный дек Quarto (`lecture/`) с `sh lecture/render.sh`; всё генерируемое — в
  `_site/`, в `lecture/` только исходники.

### Внешние инструменты (OligoArrayAux, BLAST+)
- Резолвер `R/kasp_tools.R`: тулы ищутся в `options(tidykasp.*)` → `.tools/` → кэш
  пользователя → `PATH`; `UNAFOLDDAT`/`PATH` подставляются на время вызова
  (`with_oligoarrayaux()`/`with_blast()`), ручной `export` больше не нужен.
- `kasp_install_tools()` — сборка OligoArrayAux в пользовательский кэш
  (согласие, sha256, дымовой тест, идемпотентно). `SystemRequirements` в `DESCRIPTION`;
  тесты — на `skip_if_no_oligoarrayaux()`.


### Отбор AS-праймера
- `assemble_kasp_assay()` добавил `seed_hits` (уникальность 3'-сида, обе цепи)
  тай-брейкером после дискриминации/футпринта (`specificity_candidates = 10`):
  при прочих равных предпочитается более «тихая» цепь; скрин только лидеров.

### Решение как граф и проверки уровня сайта
- **Единый реестр критериев** `kasp_criteria()` (R/kasp_criteria.R): site-проверки, гейты
  и тиры отбора, приёмочные проверки. Из него генерируются дерево
  (`kasp_decision_schema()`/`kasp_decision_graph()`/`plot_decision_tree()`), текст
  (`kasp_decision_text()`) и `kasp_design_checks()` — новый критерий появляется везде
  сразу. В дереве теперь *все* проверки (состав, структура, гетеродуплексы, SSR,
  off-target, in-silico PCR, site), а не только гейты; жёсткие лимиты — **жирной рамкой**.
- **Релаксация с границами:** мягкий гейт опускается только до лучшего достижимого
  значения и не ниже floor; жёсткие (дискриминация, вариант в AS 3'-окне) — не опускаются
  вовсе; при провале — явный `DOOMED` (один на всё дерево). Сцепленный вариант в AS 3'-окне
  не считается нарушением. `assay$decision` (граф) и `assay$evaluation` (результаты проверок) —
  в объекте, без пересчёта при отрисовке.
- **Меньше дублей в дереве:** пара гейтов, которые пара-стадия лишь повторяла за AS
  (дискриминация и AS 3'-окно — они там константы), убраны; лестница релаксации — один узел
  на стадию (а не узел на ступень); один общий `DOOMED`. 49 узлов вместо 57.
- **Схемы по химии маркера** `kasp_schemes()` (R/kasp_schemes.R): реестр химий с их
  дискриминацией, ридаутом, зиготностью и числом локусов. Тот же единый реестр критериев
  помечает каждую проверку своими схемами, поэтому из одного источника получаются деревья
  для SNP/indel/мультиаллельного KASP, **дуплексных** химий (dKASP, mxKASP, Multi-KASP:
  два локуса в одной пробирке + стадия совместимости) и **ASQ/PACE**, плюс голый
  **«любая пара праймеров»** (маркер-агностичный subset: состав, структура, повторы,
  off-target — не KASP). `kasp_criteria(scheme=)`, `kasp_decision_schema(scheme=)`,
  `kasp_decision_graph(scheme=)`, `plot_decision_tree(scheme=)`, `kasp_decision_text(scheme=)`.
  Стадия, у которой все гейты жёсткие, ведёт провал прямо в `DOOMED` (без лестницы).
  `dev/render-decision.R` рисует по картинке на химию; в README — SNP и дуплексная.
- **Off-target только по геному:** специфичность считается по референсной
  последовательности (не по списку ALT-аллелей), так что фича с несколькими вариантами
  — один локус, а не повтор; закреплено тестом.
- `kasp_design_checks()` получил уровень **site**: `target site class`, `ALT alleles (n)`,
  `linked variants skipped (n)` — вырожденные/мультиаллельные сайты и несшитые indel'ы
  видны в тумблерах; `plot_checks()` фасетит site/feature/alignment; в объекте — `assay$site`.
- `plot_primers()` при сцепленном гаплотипе показывает трек вариантов как два —
  **REF- и ALT-гаплотипы** (основания мишени и сцепленных позиций подписаны).
- Зависимости: `tidygraph` → `Imports`, `ggraph` → `Suggests`.

### Анимация поиска
- `render_design_search()` — GIF: перебор AS-размещений в порядке алгоритма (и затем
  COM), каждое подписано причиной отказа, сходятся к принятому тесту; статичные
  дорожка вариантов и мишень. `design_search_frames()` (кадры) и
  `plot_design_search_step()` — библиотечные ggplots; сборка через ImageMagick
  `convert` или `gifski`. Встроено в демо, лекцию и отчёт.
- `render_design_search()` рисует только различные кадры (повторы `hold`
  переиспользуют файл) и может распараллелить отрисовку через `BiocParallel`
  (`options(tidykasp.cores = ...)`; форк на Unix, SOCK иначе) — ~2× на многоядерной
  машине; по умолчанию один поток (как требует CRAN).

### Графика (переработка)
- `plot_primers()` — локальная карта ассая: буквы последовательностей (нейтрально чёрные),
  **жирные** — 3' GC-clamp, **красные** — основания SSR/шпильки, кольцо — мисматч/
  3'-основание, Tm числом и ΔTm-тумблером; справа — тонкая колонка тумблеров
  `plot_checks()` с фасетом site/feature/alignment/specificity; в тумблерах значение
  ≤ 4 символов (site class — код `bial`/`degn`/`comp`/`indl`).
- **Удалены** `plot_alignment()`, `plot_self_complementarity()`, `plot_dimer_qc()`
  (и `plot_discrimination()` ранее) — в текущем виде не нужны; `ggmsa` ушёл из `Imports`.
- `off-target` теперь считается по **точным** совпадениям праймера в геноме (из идеограммы),
  а не по 3'-сиду с 1 мисматчем (тот в 100+ Мб геноме шумит); порог футпринта
  `footprint_dropout` 5% → 10%.
- `dev/render-targets.R` рисует карту для каждой мишени из `inst/extdata/example_targets.csv`.

### Пакет и UX
- Проект собирается как настоящий пакет: `NAMESPACE` и `man/` из roxygen
  (все аргументы документированы, `@family` для «See also», `@examples`).
- Объект теста получил класс `kasp_assay` и `print()`-метод; хелперы
  `is_kasp_assay()`, `passes_all()`, `kasp_example()`.
- Метаданные пакета: `inst/CITATION`, `COPYRIGHT HOLDER` = автор в `LICENSE`,
  убран лишний `utils::globalVariables(".data")` (`@importFrom rlang .data`).
- Виньетка `kasp-intro`; `NEWS.md`.
- Внешний доступ к FASTA по индексу (`Rsamtools::indexFa`).

### Дизайн
- Движок InDel: аллельные шаблоны и скан стыка (`design_as_primers_indel()`,
  `allele_templates()`); размещение выбирает физика (`DECIPHER`), не ручное правило.
- `pick_target_snp()` учитывает специфичность 3'-сида (`seed_hits`); тяжесть последствия —
  по VEP `impact` или SO-термам. **Исправлено:** больше не пропускает A/T и C/G —
  транзверсии A/T, C/G дают *само-пару* (T·T, A·A) на чужой матрице — сильнейший мисматч,
  дискриминируются лучше всего; слабый 3'-мисматч — это G:T/A:C *wobble* (транзиции
  A/G, C/T), который ловится мягким `gate_wobble`, а не классом сайта.
- Substitution-блоки / MNP: вариант — строка REF/ALT, единая геометрия через
  `as_anchor()` / `block_in_seg()`.
- `design_as_primers_indel()` и SNP-путь унифицированы в `assemble_kasp_assay()`
  (пер-аллельный 3'-футпринт, ампликон — худший из двух аллелей).
- Сцепленный вариант (`linked=`) «впечатывается» в ALT-гаплотип и ALT-праймер
  (`Biostrings::replaceAt`), давая двойной мисматч к REF; REF-праймер не трогается,
  сцепление хранится в `assay$selection$linked`.
- `kasp_example(type=)` — типизированная фикстура мишеней (`inst/extdata/
  example_targets.csv`): по одной мишени на тип маркера, разбросанных по генам и
  хромосомам (SNP transition/transversion, MNP, инсерция, делеция, мультиаллельный
  сайт с двумя ALT). Локусы отобраны по *наименьшей местной нагрузке вариантов* в
  реальном гене, а не свалены из одного гипер-вариабельного региона; SNP-строки
  проверены на дизайнабельность (`dev/render-targets.R`). Колонка `gene` несёт
  AGI-локус для карты фичи; `variant_type()` стал векторным.

### Специфичность
- `specificity_seed()`: точный путь `PDict`/`vcountPDict`, мисматч-скан — параллель
  по хромосомам (`BPPARAM=`); `in_silico_pcr`/`assay_amplicons`/`assemble_kasp_assay`
  принимают `BPPARAM`.
- Батч: `design_kasp_markers()` с вердиктом на локус и параллелью по локусам.
- VCF-регион стримится через `VariantAnnotation::ScanVcfParam(which=)`.

### Термодинамика / QC
- `primer_tm()` принимает `DMSO` (%), прокинутый в `TmCalculator::tm_calculate()`;
  Ensembl REST-запросы идут через общий `ensembl_perform()` с `req_throttle(15)` и retry.
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
- `dev/check.R` — локальный «CI» (доки ↔ roxygen, styler, lintr; `--full` →
  `R CMD check`, `--spelling` → орфография); `dev/bench.R` — время и память по
  функциям (`bench`/`Rprofmem`); `dev/profile-gif.R` — профиль отрисовки GIF.
- Пакет стал lint-clean (`.lintr` в формате, понятном lintr 3.4); тесты можно
  гонять параллельно (`Config/testthat/parallel`).
- `renv.lock` запинил инструменты разработки (`styler`, `lintr`, `spelling`, `bench`,
  `profmem`, `lobstr`, `pkgdown`, `potools`); `.i18n-lib/` убран из сборки
  (`.Rbuildignore`), `ragg`/`pkgdown` объявлены в `Suggests`.
- Аудит (ROADMAP R–T): недоиспользованные возможности, переносимость Windows/macOS и
  соответствие CRAN/tidyverse.
- Тесты переведены на `withr`-идиомы (`with_seed`/`local_seed`, `local_tempfile`,
  `local_envvar(LANGUAGE = "en")`); `dev/bench.R` — на `bench`/`lobstr`;
  `object_usage_linter` снова включён.
- Переносимость: `image_tool()` различает ImageMagick 6/7 и отвергает
  `System32\convert.exe`; добавлен `.gitattributes` (`eol=lf`); `sha256_of()` умеет
  `sha256sum`/`shasum`/`openssl`.
- `@return` и `@examples` — у всех экспортов (сеть/термо/BLAST — в `\dontrun{}`);
  `expect_snapshot()` на пользовательские сообщения; `dev/bench.R` — с базовой
  таблицей `dev/bench-baseline.csv`.
- Параллель: общий `bp_lapply()` добавляет прогресс через `BPOPTIONS`; дефолты записи
  файлов — `tempdir()` (CRAN), `here` переведён в `Suggests`; из `renv.lock` убран
  мёртвый вес (8 неиспользуемых листьев).

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
