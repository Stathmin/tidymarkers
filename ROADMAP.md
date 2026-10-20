# tidymarkers — roadmap (открытая работа)

Только незакрытое. Полный ретроспективный инвентарь (что уже есть, с обоснованиями и
ссылками) — в `INVENTORY.md`; история изменений — в `NEWS.md`.

Обозначения: `[ ]` очередь · `[~]` частично · `[x]` сделано.

## Открытая работа — сводка (детали — в `INVENTORY.md` по номеру раздела)

### §U. Другие PCR-маркеры и внешние приоры
- [~] CAPS/dCAPS и ARMS/AS-PCR — довести до конца.
- [ ] «Дальние» методы под задачи (SNaPshot / single-base extension); внешние FOSS-приоры.
- [~] Провенанс порогов + термодинамическая сверка.

### §V. Модели считывания (readout simulation)
- [ ] HRM — модель плавления; dPCR — пуассоновское разбиение; CRISPR-Dx — коллатеральное
  расщепление; ARMS/CAPS/indel_size/SSR — виртуальный гель; LAMP/RPA — изотермическая кинетика.

### §W. Повторяющийся аудит
- [~] DRY / «сухой» код; дрейф документации; styler/lintr/`R CMD check --as-cran`;
  «дерево == код»; «мокрые» принципы; встроить аудит в `dev/check.R`; guard «single-base only»;
  тест `render_design_search()`.
  Сделано: LAMP B3/FIP ориентация и F1c/B1c > F2/B2, FIP 3′ ALT,
  containment зондов TaqMan/dPCR/SNaPshot/RPA, жёсткая специфичность CRISPR,
  классы AS/COM-проверок, общие тестовые фикстуры, `in_silico_pcr()` с IUPAC/
  mismatch и реальными продуктами.

### §X. Слой плана
- [ ] Двухпроходная дизайн-проводка для остальных схем.

### Прочее
- [ ] Windows/R-hub: сборка OligoArrayAux; бейджи/CI (`lifecycle`, R-CMD-check).
- [~] Скорость: устранены горячие точки по `Rprof`, быстрый/полный прогон,
  installer ставит зависимости и пакет (сделано; замерить на чистой машине).

## Z. Внешний научный аудит (окт. 2026) — код, требования и физика

Это не список пожеланий к тексту. Каждый пункт должен закрываться кодом,
тестом и (где есть данные) калибровкой. Приоритеты: **P0** — неверная
биология/вердикт; **P1** — специфичность, популяция, калибровка; **P2** —
расширение методов.

### Z0. Фактические ошибки физики/биологии  [P0]
- [~] **Cas12a collateral — ssDNase, репортёр ssDNA.** Код-комментарии уже
  говорят `ssDNA`; осталось проверить все виньетки и добавить тест на wording в
  `marker_readout()`.
- [ ] **Cas13 — не «PAM-free» без оговорок.** Cas13a: слабое предпочтение PFS
  (избегать G на 3′ фланге); Cas13b: D-PFS. Ввести `cas = "cas13a" | "cas13b"`,
  длину спейсера ~28 nt для LwaCas13a и PFS-поле в кандидатах; тесты на PFS.
- [ ] **dPCR Tm-шкала.** Стандарт: праймеры 58–62 °C, проба +8–10 °C. Текущий
  default 44–50 °C — не менять окном, пока `primer_tm()` не откалиброван против
  вендорского/`TmCalculator`; это блокер Z2/Z8. Временный default помечен как
  калибровочный долг, не как норма.
- [x] **LAMP F1c/B1c > F2/B2.** `lamp_fuse()` теперь фильтрует `t1 > t2`;
  FIP 3′ несёт ALT, а не REF.
- [ ] **KASP: «Taq не продлевает 3′ мисматч» — кинетика, не бинарность.**
  Заменить формулировку и добавить в `marker_design_checks()`/DAG различитель
  типа мисматча (transition/transversion, G:T/A:C, wobble) как `unknown`/риск,
  а не «pass».

### Z1. Вердикт и статусы  [P0]
- [~] `NA` больше не проходит `select_with_ladder()`; DAG получил статус
  `unknown`; но acceptance/verdict всё ещё не различает `unknown` и
  `not_applicable`.
- [ ] Ввести тип `pass | fail | unknown | not_applicable` для каждой проверки.
  `marker_design_checks()` возвращает его; `marker_decision_graph()` и
  `marker_readout()` считают вердикт только из `pass`, а `unknown` не считают
  успехом. Тест: любой `NA` не даёт «all checks pass».
- [ ] Гейты, которые не могут провалиться (`gate_lamp_variant = TRUE`,
  `duplex homologous = NA`), сделать конструкторскими инвариантами или
  `not_applicable`, а не декоративными гейтами.

### Z2. Один источник истины для порогов  [P1]
- [ ] Каждый критерий — data-object: `id`, `scope`, `kind`, `bounds`, `units`,
  `source`, `calibrated = yes/no`. Селектор, checker, plotter и DAG читают один
  объект; расхождение (design GC 0.35–0.70 vs check 0.40–0.60, Tm display vs
  gate, DAG label vs constant) становится невозможным. Тест на отсутствие
  дублирующих констант.
- [ ] `MARKER_CHECK_*`/`MARKER_FAMILY_CHECKS` — начало; перенести в них limits,
  levels и evaluator-id, чтобы семья не собиралась из трёх мест.

### Z3. Условия реакции  [P1]
- [ ] `reaction_conditions` на схему: соли, Mg, dNTP, DMSO/бетаин/PEG,
  температура. `primer_tm()` и `dimer_thermo()` принимают его; ΔG димеров
  считается при температуре реакции (60–65 °C для LAMP, 37–42 °C для RPA),
  а не всегда 37 °C. RPA не получает PCR-Tm как гейт.

### Z4. Движки и capabilities  [P1]
- [ ] Thermo/specificity/scoring — интерфейсы с `capabilities()`. Если движок
  недоступен (OligoArrayAux, BLAST, crisprScore), результат `unknown`, а не
  молчаливый `pass`. `TmCalculator` — независимый cross-check шкалы Tm.
- [ ] **Ecosystem reuse, not reimplementation.** `crisprVerse` (Cas9/Cas12a/
  Cas13d; Bowtie/BWA; `crisprScore`) для PAM/search/off-target — но не как
  trans-cleavage activity predictor и с оговоркой human-only SNP annotation;
  `DECIPHER` `DesignPrimers`/`DesignProbes`/`DesignSignatures` для group/
  homeolog-aware design; `rprimer` для degenerate-панелей; `TmCalculator` (27 NN,
  salt/DMSO) как независимый Tm-cross-check; `fitPoly` для dosage. Внешние
  Primer3/MFEprimer/NUPACK/ViennaRNA — только через capability-интерфейс;
  `openPrimeR` — проверить maintenance перед зависимостью.

### Z5. Mismatch-tolerant specificity  [P1]
- [x] `in_silico_pcr()` ищет IUPAC-праймеры, поддерживает `max_mismatch` и
  склеивает реальные forward/reverse продукты; exact-mode предупреждает о
  near-match продукте.
- [ ] Position-aware off-target: 3′-конец, число и позиция мисматчей; отдельные
  модели для RPA (несколько мисматчей), CRISPR (seed/PAM), паралогов/гомеологов.
  Не считать exact-hit + 15-nt seed достаточной специфичностью.
- [ ] Backends: crisprVerse/Bowtie/BWA для CRISPR, DECIPHER для group-specific
  design, rprimer для degenerate-панелей, Primer3/MFEprimer/NUPACK/ViennaRNA
  как внешние проверки. Каждый — с provenance и capability-флагом.

### Z6. Популяционный слой  [P1]
- [x] `marker_attach_haplotypes()` + `marker_population_checks()` + generic
  haplotype track; CAPS/CRISPR предупреждают без популяции и при вариантах в
  сайте/guide-PAM.
- [ ] Allele-frequency/haplotype object течёт во **все** схемы; отсутствие
  частоты = `NA`, не вес 1. Проверять: CAPS site, CRISPR guide/PAM, LAMP/RPA
  primer sites, SSR null alleles, polyploid dosage. Полиплоидия — отдельная ось,
  не «pop_freq».

### Z7. Намеренные вторые мисматричи  [P1]
- [ ] Параметр deliberate mismatch для KASP/ARMS/AS-LAMP/AS-RPA/CRISPR
  (SHERLOCK synthetic mismatches): позиция, тип, predicted ΔΔG; тесты, что
  слабый transition SNP усиливается, а не просто помечается «риск».

### Z8. Калибровка и benchmark  [P1]
- [ ] Harness/`targets` pipeline: gate-pass → wet-lab outcome; published
  conversion rates; `calibrated`-флаг у каждого порога. Без этого «usable» —
  эвристика, а не вероятность.
- [ ] Analytical performance: LOD, minor-allele fraction, stutter, false-positive
  rate, dPCR partitioning; отдельные поля отчёта, не «not modeled» в тексте.

### Z9. Пробелы по методам  [P1–P2]
- [ ] **LAMP:** paired ALT set или internal control; REF-only ALT = «нет сигнала»
  неотличим от failed reaction; F2–B2 distance, loop length и 3′-stability —
  явные гейты.
- [ ] **HRM:** SNP-class awareness (A/T class 4, Liew 2004), melt-domain/uMELT,
  heteroduplex; A/T warning уже есть, модель плавления — нет.
- [ ] **TaqMan:** MGB/LNA probes (short probes) и probe Tm/MGB-специфичные правила.
- [ ] **Sanger:** ICE/TIDE/Poly Peak deconvolution, minor allele ~15–20% limit,
  primers-over-variants allele dropout.
- [ ] **CRISPR:** screen several spacers; activity scoring; crisprScore/Bowtie;
  reporter/readout kinetics.
- [ ] **dPCR:** Poisson/partition/LOD/minor-fraction model; не только Tm.

### Z10. Структура, плоидность, аналитика  [P2]
- [ ] Amplicon/target folding (особенно LAMP/RPA/Cas13 RNA, GC-rich HRM).
- [ ] Dosage/copy-number (fitPoly/mixture model) для полиплоидных KASP.
- [ ] GRanges/VCF-backed interchange вместо custom tibble там, где это уместно.

### Z11. Документация  [P2]
- [ ] «Аудит решений» переписать из dev-заметок в user-facing ограничения и
  рекомендации; каждое утверждение — с source/calibrated-флагом.
- [ ] Исправить Cas12a/Cas13/dPCR/LAMP/KASP формулировки во всех виньетках и
  комментариях; тест на запрещённые термины (`dsDNase`, «PAM-free» без PFS).

### Уже закрыто в этом проходе
- [x] 29 family-specific checks (`MARKER_CHECK_*`), DAG-узлы и статус `unknown`.
- [x] `NA` не проходит отбор; CAPS/CRISPR population warnings.
- [x] LAMP F1c/B1c > F2/B2, FIP 3′ ALT; RPA THF off variant.
- [x] `in_silico_pcr()` IUPAC/mismatch/real products; near-match warning.
- [x] Generic haplotype track/population checks; examples with known proportions.
- [x] Template/plot fixes (noses, short-feature letters, boxes, Sanger roles).

## ZA. Целевые интерфейсы и швы (codebase-design)

Это план перестройки **модулей и швов**, а не список функций. Принцип:
табличные данные — `tibble`/`vctrs`; поведение — S3-generic; граф —
`tidygraph`; рисование — грамматика `ggplot2`; внешний движок — adapter с
`capabilities()`. Не тащить `R6`/`data.table`; не превращать поведение в
data-frame только ради «tidy».

### ZA1. Целевой путь вызова

```r
context <- marker_context(
  genome, variants = NULL, population = NULL,
  conditions = reaction_conditions("pcr"), engines = marker_engines()
)
scheme <- marker_scheme("caps")
m  <- scheme$design(snp, context)
ev <- scheme$evaluate(m, context)
tr <- scheme$tracks(m, context)
dc <- scheme$decision()
ro <- scheme$readout(m, context)
```

### ZA2. Модули и швы

| Модуль | Tidy/ggplot-форма | Заменяет | Зависимость |
|---|---|---|---|
| `MarkerContext` | classed list, `rlang::list2`, `withr` | ручной `...`-плумбинг | существующие |
| `Scheme` | S3 + `tibble`-реестр | `design_<scheme>_marker` + `switch(scheme)` + naming convention | нет |
| `Criterion`/`CriterionSet` | `tibble` + list-column evaluator, `vctrs::vec_rbind` | дубли порогов в design/check/DAG | `vctrs` (Imports) |
| `Evaluation` | `tibble` + `forcats` status | три evaluator-пути | `vctrs`, `forcats` |
| `Population` | `tibble` + `GenomicRanges`/`VariantAnnotation` adapters | мутация `marker_attach_haplotypes()` | `fitPoly` позже (Suggests) |
| `Engines` | S3-adapter + `capabilities()` | глобальный `PCR_COND`, молчаливые skip | существующие + внешние |
| `DecisionSchema`/`DecisionState` | `tidygraph` + Graphviz | ручная раскладка и статусы | `ggraph` (Suggests) |
| `TrackSpec`/`Figure` | `ggplot2` + `ggenes`/`patchwork` | стиль внутри схем | `scales` (Imports) |
| `Readout`/`Calibration` | model-result + `ggplot2`/`yardstick` | строки «not modeled» | `yardstick`, `rsample`, `targets` (Suggests) |
| YAML-технология | `yaml` | ручная сериализация | `yaml` |
| Отчёты/таблицы | `knitr::kable`/`gt` | ручной текст | опционально Suggests |

### ZA3. Home-rolled → проверенный импорт

- `spec_bind()`/ручная сборка строк → `purrr::list_rbind()`/`vctrs::vec_rbind()`.
- `format_check_value()` → `scales::percent()`, `scales::number()`.
- progress/messages → `cli::cli_progress_*()`, `cli::cli_abort()`.
- cache key → `rlang::hash()`.
- IUPAC/regex-праймеры → `Biostrings::matchPattern(fixed = FALSE)`, `PDict`.
- интервалы → `IRanges`/`GenomicRanges`.
- раскладка графа → `tidygraph` + Graphviz `dot`; `ggraph` только для R-view.
- test doubles → `testthat::local_mocked_bindings()`.
- таблицы в доках → `knitr::kable()`/`gt`.
- метрики/калибровка → `yardstick`; resampling → `rsample`; pipeline → `targets`.

### ZA4. Политика зависимостей

- **Imports (добавить при прямом использовании):** `vctrs`, `scales`, `yaml`.
- **Suggests:** `yardstick`, `rsample`, `targets`, `fitPoly`, `crisprBase`,
  `crisprDesign`, `crisprScore`, `crisprBowtie`, `crisprBwa`, `vdiffr`, `gt`.
- **Уже есть и переиспользовать:** `Biostrings`, `IRanges`, `GenomicRanges`,
  `Rsamtools`, `BiocParallel`, `DECIPHER`, `TmCalculator`, `rBLAST`, `Rsubread`,
  `rprimer`, `dplyr`, `tidyr`, `purrr`, `forcats`, `tibble`, `ggplot2`,
  `ggenes`, `ggforce`, `ggrepel`, `patchwork`, `tidygraph`, `cli`, `withr`,
  `rlang`, `glue`, `processx`.
- **Не добавлять без второго adapter:** `R6`, `data.table`, `DiagrammeR`,
  `openPrimeR` (сначала проверить maintenance).
- Тяжёлые `crispr*`/`fitPoly` — только при реальной реализации
  CRISPR-search/dosage, не «на будущее».

### ZA5. Что не делать tidy/ggplot

- `Scheme` и `Engines` — поведение, S3-generic, не data-frame.
- `DecisionSchema` — графовые данные + state; рисование отдельным adapter.
- `Readout` — typed result (`simulated`, `lod`, `maf`, `stderr`), не ggplot.
- Не вводить seam для одного adapter (deletion test): если удаление модуля
  просто возвращает сложность одному caller, это pass-through.

### ZA6. Порядок миграции

1. `vctrs`-typed `Criterion`/`Evaluation`; статус
   `pass/fail/unknown/not_applicable`.
2. `MarkerContext` + `ReactionConditions` + `Engines` (без тяжёлых deps).
3. `Scheme` S3-adapter вокруг существующих функций; схемы мигрируют по одной.
4. `Population` в context; `marker_attach_haplotypes()` — compatibility adapter.
5. `DecisionSchema`/`DecisionState`; `tidygraph` + Graphviz.
6. `TrackSpec`/`Figure`; `scales`; стиль только в renderer.
7. `Readout` + `yardstick`/`targets`.
8. `yaml` export; затем optional `fitPoly`/`crispr*`.

### ZA7. Acceptance / test surface

- Вызовы в тестах идут тем же швом: `marker_context()` → `scheme$design()` →
  `scheme$evaluate()` → `scheme$tracks()`/`scheme$readout()`.
- Contract-test: каждая схема реализует `design/evaluate/tracks/decision/readout`.
- `NA` никогда не становится `pass`; `unknown` не считается успехом.
- У каждого criterion один bounds-object и `source`/`calibrated`.
- Fake thermo/specificity engines дают deterministic тесты без OligoArrayAux.
- Golden-тесты: criterion-id, evaluation-table и track-rows; не пиксели.
- После каждого шага — `graphify update .` (если CLI установлен) и `dev/check.R`.

### ZA8. Docs/artifacts graph и freshness

Graphify сейчас индексирует код и doc-headings, но **не** файлы документов и **не**
артефакты: `graphify explain vignettes/markers-gel.Rmd` и
`vignettes/figures/caps.svg` не находят узлов. Попытка sidecar-merge также
провалилась: `graph.json` хранит `nodes` как dict, а `graphify merge-graphs`
ожидает node-link list (проверено). Поэтому docs/artifacts слой делаем через
manifest + generated Markdown, а не через sidecar JSON.

- [x] **`docs/INDEX.md` — единый вход.** Все README/AGENTS/CONTRIBUTING/ROADMAP/
  NEWS/INVENTORY/vignettes/dev-docs; каждая ссылка проверяется на существование.
- [x] **`docs/artifacts.yml` — источник истины.** Для каждого артефакта:
  `path`, `kind` (figure/decision/table/benchmark), `generator`, `inputs`,
  `input_hash`, `output_hash`, `docs`, `code`, `status`.
- [x] **`dev/check-docs.R`.** Нет orphan-документа; нет broken relative link;
  каждый `include_graphics()` имеет manifest-entry; обратные ссылки совпадают.
- [x] **`dev/check-artifacts.R`.** Нет orphan-файла в `vignettes/figures`,
  `man/figures`, `results/plot-smoke`; все пути manifest существуют; входной хеш
  не изменился; выходной хеш совпадает — иначе **stale** и gate fail.
- [x] **`docs/generated/artifacts.md`.** Один heading на артефакт со ссылками
  на generator/source/docs; `graphify update . --no-cluster` после генерации.
  Heading label = путь артефакта, чтобы `graphify query "caps.svg"` находил
  docs/code, даже если самого файла нет в AST-графе.
- [ ] **`targets` (Suggests).** `tar_outdated()`/`tar_make()` для gallery,
  decision graphs и vignettes; manifest обновляется из `tar_meta()`. Сейчас
  freshness держат `dev/render-gallery.R` + `dev/artifacts.R --write/--check`;
  `targets` — следующий шаг, не блокер.
- [x] **Gate.** `dev/check.R` вызывает docs/artifacts checks; graphify остаётся
  developer tool (`--graphify` opt-in), не CI-зависимость.
- [ ] **Graphify health.** `diagnose multigraph` сейчас показывает 41
  same-endpoint collapse и 291 external node без `source_file`; учитывать как
  ограничение ориентации, не как доказательство корректности.
- **Состояние:** реализовано в Wave 1 (коммиты `9193e6c`/`e5a581a`).
  `docs/artifacts.yml` хеширует весь `R/`, поэтому **любая** правка `R/` делает
  20 SVG stale: после неё — `dev/render-gallery.R`, затем `dev/artifacts.R --write`.
  PNG галереи git-ignored; SVG коммитятся.

## ZP. Параллельный rollout (coordination)

Правило: **ветки изолируют историю, не рабочее дерево.** Один clone = один
checkout. Никаких N агентов с `git checkout` в одном каталоге.

### ZP0. Base snapshot — блокер  [done]
- [x] Snapshot сделан: ветка `wip/roadmap-base` (коммит `2785097`, 163 файла),
  затем Wave 1 (`e5a581a`) и Wave 0 (`5e31a1a`).
- [x] `main` не трогали (`2077c46`); слияние в `main` — отдельным решением.
- [x] Worktree-only правило соблюдено: агентские ветки слиты и удалены,
  остался один worktree (`main`-clone на `wip/roadmap-base`).
- [x] **Grace-период:** 9 untracked `R/`-файлов, 5 vignettes и 24 man/dev файла
  были реальным исходником и закоммичены; PNG галереи — local-only, SVG — tracked.
- [ ] Новые worktree — только от `wip/roadmap-base`.

### ZP1. Волны и ownership

**Wave 0 — сделано (`5e31a1a`), последовательно (shared files):**
`Evaluation`/status (`marker_check_status()`, `evaluation_verdict()`),
`MarkerContext`/`ReactionConditions`/`Engines` (`marker_run_context()`,
`reaction_conditions()`, `marker_engines()`), KASP-plot нейтральный статус
(`1bc1798`). Файлы: `R/marker_checks.R`, `R/marker_plot_methods.R`,
`R/marker_core.R`, `R/marker_thermo.R`, `R/marker_haplotypes.R`,
`tests/testthat/helper-*`, shared tests. Один владелец.

**Lean test loop (сделано):**
- `Rscript dev/test.R` — 1374 pass, **68 s** (shared contracts; дефолт).
- `Rscript dev/test.R --scheme crispr` — ~66 s.
- `Rscript dev/test.R --full` — 3164 pass, **9m45**, external checks on, 0 skips.
- `Rscript dev/check.R --quick` — **38 s**; `--no-render` — **43 s**.

**Wave 1 — сделано (`e5a581a`), disjoint paths:**

| Agent | Worktree | Branch | Owns | Forbidden |
|---|---|---|---|---|
| docs-artifacts | `../tidymarkers-agent-docs` | `feat/docs-artifacts` | `docs/**`, `dev/check-docs.R`, `dev/check-artifacts.R`, `dev/artifacts.R` | `R/**`, `vignettes/**` |
| crispr-physics | `../tidymarkers-agent-crispr` | `feat/crispr-physics` | `R/marker_crispr.R`, `tests/testthat/test-crispr.R` | `vignettes/**`, other `R/**` |
| sanger-physics | `../tidymarkers-agent-sanger` | `feat/sanger-physics` | `R/marker_sanger.R`, `tests/testthat/test-sanger.R` | `vignettes/**`, other `R/**` |
| lamp-physics | `../tidymarkers-agent-lamp` | `feat/lamp-physics` | `R/marker_lamp.R`, `tests/testthat/test-lamp.R` | `vignettes/**`, other `R/**` |

**Wave 2 — после Wave 0/1:**
- dPCR Tm calibration — только после engine/conditions seam и калибровочной базы.
- Vignettes — один агент на файл, строго после code-агентов; никто не делит
  `markers-kasp.Rmd`/`markers-gel.Rmd`/`markers-probe.Rmd`/
  `markers-isothermal.Rmd`/`markers-crispr-sequencing.Rmd`.
- Decision/plot split — после `Evaluation`, тот же owner, что и Wave 0.

### ZP2. Merge order

1. `wip/roadmap-base` → `main` (или integration worktree).
2. `feat/docs-artifacts` — новые пути, низкий конфликт.
3. `feat/crispr-physics`, `feat/sanger-physics`, `feat/lamp-physics` — по одному.
4. Vignette branches — последними.
5. После каждого merge: `Rscript dev/check.R --no-render`; полный gate — в
   integration worktree. Конфликты решать в integration worktree, не в агентских.

### ZP3. Agent brief

```text
Worktree: /abs/path/tidymarkers-<slug>
Branch: feat/<slug> (already checked out; do NOT checkout elsewhere)
Owns: <paths>
Forbidden: <paths>
Tests: Rscript dev/test.R <filter>; Rscript dev/check.R --no-render
Commit only on feat/<slug>; verify git branch --show-current before commit.
```

### ZP4. Anti-patterns

- N branches in one clone without worktrees.
- Worktrees from stale `HEAD` while 144 paths are dirty.
- Two agents editing the same vignette or shared helper.
- Parallel dPCR calibration before Tm engine/benchmark.
- Merge without the gate, or deleting a worktree before merge.

## Y. Следующий слой: конструктор технологий, цели по области, специфичность  [план]

Основание — два разбора коллег (окт. 2026): (1) API в духе ggplot
(`targets` + `context` + `technologies` + `solve` + `exports`) и (2) требования к
специфичности у гексаплоида (пшеница, гомеологи A/B/D). Это — прямое
развитие уже сделанного: «перебор технологий на **сайт**» расширяется до «перебор
**сайтов** и технологий на **область**». Движок (дизайн, гейты, тиры, DAG) после
перестройки §W готов; основной новый код — слой оценки и ранжирования кандидатов.

### Y0. Имя пакета и API: `tidykasp` → `tidymarkers`  [сделано]
- **Зачем.** Пакет давно не «только KASP»: 20 химий, маркер-агностичный слой (§X) и
  `marker`-обобщение; имя `tidykasp` вводило в заблуждение.
- **Сделано (0.4.0, ломающее):** имя — **`tidymarkers`** (мн. ч.); ренейм пакета
  (файлы/URL/CITATION/бейджи/домен переводов) и обобщение API `kasp_*` → `marker_*`
  (KASP-специфичное оставлено). Детали — в `NEWS.md` (0.4.0); алиасов старых имён
  нет (пакет < 1.0).

### Y1. Конструктор технологий (композиция из блоков)  [фундамент]
- **Сейчас:** схема = `scheme_def_<m>()` + `design_<m>_marker()`, гейты/тиры/проверки
  декларативны, семейства `KASP_CHECK_PRIMER`/`_SITE`/`_AS` (`marker_criteria()`);
  `marker_tech(schemes|preset|exclude, costs=)` уже параметризует набор схем в пуске.
- **Цель:** технология собирается пользователем из **блоков-фильтров и параметров**,
  а не пишется с нуля: готовые блоки (состав/шпильки/димеры, повторы, специфичность,
  3'-якорь, ампликон, readout), и из них — «гига-KASP» = KASP + доп. ограничение
  (напр. «праймеры подлиннее», «шпильки строже»). Набор доступных методов — параметр
  пуска, а не зашитая логика (уже принцип `marker_tech`).
- **Артефакт:** `marker_technology()` (конструктор: семейство + блоки + параметры),
  реестр блоков `marker_filter_blocks()`, `marker_technologies()`; валидация, что собранная
  технология — корректный DAG (`expect_dag_invariants()`).
- **Формат обмена — декларативный, не код.** Технологию для шаринга сериализовать в
  **YAML** (блоки + параметры как данные); R-файл — только для локальных кастомных
  дизайнеров и не является форматом обмена (шаринг R-файла = выполнение чужого кода).
- **Зависимость:** перестройка DAG (§W) — реестр проверок/гейтов/семейств уже единый.
- **Генерализация до `marker`.** Модуль давно шире «просто KASP»: там, где семантика
  не KASP-специфична (проверки, критерии, фрагмент-маркеры, размерные маркеры),
  функции/файлы/описания следует обобщать до `marker` (маркер-агностичный API и доки),
  оставляя `marker_`-префикс только действительно KASP-специфичному. Ломающих
  переименований не делать: вводить общие имена с алиасами/депрекейтом.

### Y2. Цели по области (goal = что угодно из группы сцепления)
- **Сейчас:** `pick_target_snp()` берёт **один** SNP/MNP; `marker_solve()` перебирает
  технологии на сайт; `marker_chunk()`/`haplotype_cover()` дают локальное окно.
- **Цель:** цель — **область** (фазированный блок freebayes / группа сцепления /
  интервал): перебрать сайты × технологии, вернуть лучший по «цена-качество» маркер
  **и его тип**; артефакт пуска показывает, как и из чего выбиралось; автоматическая
  таблица **«гаплотип → какой сигнал маркера»**.
- **Дешёвые префильтры «проклятых мест» — по типу маркера, не хардкод.** «Плохой
  контекст» задаёт сама технология: KASP/HRM вредят GC-экстремумы, гомополимеры,
  повторы; для SSR повторы — *цель*; для FISH тандемные повторы нужны; для CAPS важнее
  наличие/отсутствие сайта. Префильтр — блок технологии со своими порогами. Для
  протяжённых участков — «3 маркера топ на области, между ними 1 kb мрака».
- **Осторожно:** сцепление врёт в обе стороны (потеря и приобретение ложного) при
  ошибках генотипирования и полиплоидии → надёжно лишь **в пределах фазированного
  блока / длины рида**.
- **Зависимости:** Y1, Y4; готовые `marker_plan`/`marker_solve`/`marker_chunk`.

### Y3. Специфичность/off-target — политика технологии  [начато]
- **Сейчас:** `in_silico_pcr()`, `specificity_seed()`/`primer_hits()`,
  `KASP_UNGATED_CHECKS`; off-target гейтится лишь у части схем; `select_with_ladder()`
  ослабляет мягкие гейты.
- **Не «с ноги»:** off-target не всегда порок — для qPCR/мультиплексной детекции
  внеплановая амплификация бывает частью задачи, для FISH неспецифичность ожидаема.
  Поэтому специфичность — **политика технологии**
  (`specificity_mode = "unique" | "homeolog-aware" | "multicopy" | "off"`), а не
  глобальный жёсткий гейт.
- **Цель:**
  1. **Нет чистой пары — нет маркера.** Сейчас при отсутствии чистой пары берётся
     ранговый лидер с предупреждением (§H); для специфичности это недопустимо —
     технология, требующая уникальности, возвращает **NO marker с причиной**, а не
     помеченный off-target'ами маркер. Инвариант «off-target-гейт не в relax-наборе»
     проверяется тестом, как `expect_dag_invariants()`.
  2. **`homeolog-aware` (пшеница).** Гомеологи A/B/D амплифицируются, но AS-праймер
     различает их 3'-концом: хиты на гомеологах допустимы, **только если 3'-мисматч к
     ним гарантирован**, и артефакт показывает, на каком SNP держится субгеномная
     специфичность (`unique | multicopy | off` этот случай не выражают).
  3. **3'-end-aware off-target:** позиция и число несовпадений у 3'-конца, а не только
     совпадения по длине.
  4. **Ленивая оценка поверх ранга:** идти по списку и гонять in silico PCR **пока не
     наберётся k прошедших**; перед ним — дешёвый префильтр уникальности 3'-k-меров.
  5. **Индекс — решить до старта (масштаб).** §O откалиброван под Arabidopsis (120 Mb);
     для пшеницы (~15 Gb) построчный скан по хромосомам через `Biostrings` не годится, а
     `PDict(max.mismatch=)` медленный. Выбрать заранее: хеш 3'-12-меров / BLAST-индекс /
     **FM-индекс** (уже есть `marker_fm_index()`/`specificity_fm()`, подтверждение — BLAST).
     Иначе «ленивая оценка» повиснет.
  6. **Опция доп. сборок/пангенома:** прогонять и по ним (хотя бы финальных кандидатов) —
     одна сборка не видит паралогов/вставок у линий пользователя.
- **Сделано (0.4.0):** `marker_specificity()` + `KASP_SPECIFICITY_MODES` и «нет чистой
  пары → нет маркера» (`assemble_kasp_assay(on_offtarget=)`) для `unique`-схем.
  Остаётся: 3'-end-aware проверка, ленивая оценка + выбор индекса (под 15 Gb),
  опция доп. сборок, логика `homeolog-aware`.
- **Зависимость:** Y1 (технология несёт режим специфичности).

### Y4. Многоцелевое ранжирование и Парето-фронт
- **Сейчас:** жёсткие гейты + лексикографические тиры (узкое место `min`);
  `marker_risk()` low/medium/high; `marker_solve()` идёт по стоимости.
- **Цель:** трёхмерный фронт — **цена**, **качество ассея** (включая информативность в
  панели) и **риск потери сцепления** (нужна генетическая/физическая дистанция или
  карта рекомбинации); выдавать **Парето-фронт** с рекомендованной точкой, а не одно
  «лучшее» решение. Четыре оси не вводить — фронт плохо читается.
- **Артефакт:** `marker_pareto()`; `marker_informativeness()` (по доступным гаплотипам/VCF —
  маркер, различающий двух родителей, может быть бесполезен на популяции; **входит в ось
  качества**); модель дистанции/рекомбинации; интеграция в отчёт `marker_solve`.

### Y5. Инпуты: что уже есть и что доработать
- **Уже есть и работает индексно (переписывать не нужно):** `read_variants_vcf()`
  (vroom / `VariantAnnotation`, регион по tabix), `read_genome_fasta()` + `marker_genome()`
  (индексированный `Rsamtools` / `BSgenome`, `chrom_map`, `seqids`), `read_transcripts()` +
  `haplotypes_from_alignment()`, `marker_chunk()`, аннотация (`marker_annotate.R`).
- **Уточнённый список доработок:**
  1. ридер **GenBank** (`read_genbank()`) — быстрый, с индексацией;
  2. импорт **файлов заказа праймеров** (IDT и родственные) — `read_primer_orders()`;
  3. **глубина из VCF** (DP/AD) как вес уверенности в генотипе: низкопокровные
     гетерозиготы штрафовать/исключать (`marker_context(depth=)`);
  4. обратная «сборка пары»: `find_mate_primer()` — по заданному первому праймеру и
     мишени найти второй;
  5. выравненный фасад разнородных входов к одной паре `targets` + `context`.
- **Фасад для внешних пайплайнов:** вход VCF+DP → готовые кандидаты маркеров.
- **Стабильность фасада.** Для использования пакета как зависимости нужен стабильный
  фасад (`marker_context(depth=)` и т.п.) и честная пометка о нестабильности API
  (`lifecycle`/`experimental` с явным предупреждением в доке).
- **Гексаплоид и глубина.** Низкая глубина и картирование на паралоги портят не только
  гетерозиготы, но и сами вызовы гомозигот: слабую уверенность штрафовать/исключать, а
  не ограничиваться фильтром гетерозигот.

### Y6. Валидация, бенчмарки и обратная задача
- **Минимальный регресс-харнесс — уже в фазе A.** Без набора опубликованных маркеров
  не понять, улучшил ли новый off-target-гейт и ранжирование результат или просто
  изменил его. Первый источник ground truth — уже отобранные/опубликованные маркеры из
  смежного проекта; полный набор — здесь.
- Набор известных случаев, где «правильный» маркер уже **опубликован** (KASP/CAPS/…),
  и проверка, что алгоритм хотя бы не хуже человека; регресс-харнесс; аргумент для
  публикации. `inst/extdata/validated_markers/`, `dev/bench-markers.R`.
- **Обратная задача — «проверить текущий набор маркеров»:** на вход подаётся готовый
  маркер (явная последовательность, фрагмент или один праймер из пары + мишень), тип
  маркера задан пользователем; тул прогоняет его по приёмке (состав, структура,
  специфичность, дискриминация), **не проектируя заново**. Тесты — на опубликованных парах.

### Y7. Артефакт «почему этот маркер» и документация
- Довести начатое в §X: слой объяснения (какие ограничения ослаблены, почему не CAPS,
  почему этот сайт), таблица «гаплотип → сигнал», шарящиеся технологии. Главная боль
  таких тулов — непонятно, почему выдан именно этот маркер; это снимает её.

**Фазы:**
- **A (разблокирует пшеницу):** Y3 (off-target вне relax + «нет чистой пары → нет
  маркера» + `homeolog-aware` + выбор индекса под 15 Gb) + Y1 (конструктор, YAML-обмен)
  + **минимальный регресс-харнесс на опубликованных маркерах** (из Y6).
- **B:** Y5 (инпуты, стабильный фасад) + Y2 (цели по области).
- **C:** Y4 (Парето 3D) + остальная валидация (Y6).

