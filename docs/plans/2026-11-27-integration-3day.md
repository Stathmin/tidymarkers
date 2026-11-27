# Интегративная программа на три дня — 2026-11-27

- Ветка плана: `plan/3day-20261127`, база `main` @ `15a3c41`.
- Статус: **план**; этот коммит меняет только документы, кода не касается.
- Решения: [ADR 0008](../adr/0008-user-scenarios-as-contract.md) (J1–J7),
  [ADR 0009](../adr/0009-panel-orchestrator-contract.md) (панель),
  [ADR 0010](../adr/0010-index-strategy-at-scale.md) (масштаб),
  [ADR 0011](../adr/0011-release-artefact-and-mcp-contract.md) (продукт).
- Очередь: [`ROADMAP.md` §ZV](../../ROADMAP.md), пункты ZV1–ZV10.
- Срез поверхности: [`docs/user-api-surface.md`](../user-api-surface.md) —
  **207 экспортов**, **48 S3-методов**, 255 строк CSV (137 пользовательских),
  **934 типизированных ребра**, 58 точек входа, 116 листьев, 29 изолятов.

## 0. Что уже есть и что не входит

Опора программы — не новая архитектура, а уже принятые решения и замеры:

- ADR 0005 (индексная стратегия), 0006 (швы панели), 0007 (швы модулей и
  MCP-контракт); `ZT T1`–`T13` закрыли индекс, шарды, BiocParallel, фикстуру
  фазовых гаплотипов, оркестратор панели, реестры `plot`/`print`, разрезание
  god-файлов, валидацию `tools/call`, `digest`/`processx` и vctrs-класс.
- Замеры: TAIR10 memory vs per-chromosome bowtie — медиана 18.64 с vs 18.32 с,
  индекс/память 0.98×, максимум 1.18×; `VmHWM` — 1.7 Б/п.н. над baseline,
  консервативный бюджет 5 Б/п.н.; монолит 16 Гб = 21.5 ГБ индекса не влезает
  в 16 ГБ RAM (`docs/decisions/0001-16gb-fixpur.md`).
- Дефекты, которые программа закрывает: ZU1 (`panel` — тупик), ZU2 (`outcomes`
  — тупик), ZU3 (`report` — тупик), ZU4 (`primers` — два контракта), ZU5
  (справочники без моста), ZU6 (нет генотипов и строки `variant`), ZU7
  (дубликаты/дрейф описания), ZU8 (граф и CSV гниют).

**Не входит:** новая наука `Z7`/`Z8`, Rust-адаптер, DuckDB, новые химии,
GUI, CRAN-сабмит, Python-клиент MCP, второй диалект схемы. Всё это — отдельные
решения, а не трёхдневный объём.

### Правила программы

1. Этот коммит — только документы: ADR 0008–0011, §ZV, этот план, ссылки в
   `docs/INDEX.md` и `docs/adr/README.md`.
2. День заканчивается зелёными гейтами. Красный гейт не переносится на
   следующий день; срезается объём.
3. Вердикт на золотом наборе не меняется ни одним рефакторингом (`D1`–`D5`).
4. Новый экспорт допустим только вместе со строкой в
   `docs/user-api-surface.csv`, обновлённым DOT и названным владельцем (ADR или
   ZV). Иначе — не экспорт.
5. Пользовательский шаг не требует `$`, `[[` или ручного вытаскивания колонки.
   Мост — экспортированный аксессор (`as_marker_solved()`,
   `as_primer_set()`), а не знание внутренней формы.
6. Срез при отставании: сначала `ZV3`/`ZV4`, затем `ZV5`; `ZV1`, `ZV2`,
   `ZV6`, `ZV9`, `ZV10` не срезаются.

## 0.1 Контракт метрик — с первого дня, не в конце

Метрики — такой же контракт, как сценарии и артефакт. Базовый прогон
снимается на `main` **до** первой правки кода; каждый следующий прогон пишет
тот же формат. Прогон без строки метрик не считается прогоном.

**Файл:** `data/metrics/<date>-<run>.jsonl` — один JSON-объект на строку, одна
строка на стадию. Поля фиксированы:

| Поле | Тип | Смысл |
| ---- | --- | ----- |
| `run_id` | string | стабильный id прогона; для reproducibility-прогона детерминирован |
| `stage` | enum | `scenario`, `ingest`, `design`, `specificity`, `thermo`, `checks`, `artefact`, `mcp`, `gate` |
| `object_id` | string | id, `artefact_id` или `scenario:JN`; не свободный текст |
| `duration_ms` | integer | длительность стадии |
| `status` | enum | `ok`, `warn`, `fail`, `skip`, `error` |
| `retries` | integer | повторы после первой попытки |
| `tokens_in` | integer/null | входные токены LLM; null, если модель не вызывалась |
| `tokens_out` | integer/null | выходные токены LLM; null, если модель не вызывалась |
| `cost_usd` | number/null | оплаченная стоимость; `0` для локального счёта, null если неизвестна |
| `cache_hit` | boolean | стадия переиспользовала шард, артефакт или кэш |
| `queue_depth` | integer | задач в очереди на старте стадии |
| `language` | string | `en`, `ru` или `und`; язык прогона, не текст |

**Содержимого нет.** Промпты, последовательности, метки целей, пути и тексты
ошибок в метрики не пишутся. `object_id` — id или хеш, не метка. Валидатор
`dev/metrics.R --check` отвергает неизвестные ключи, поэтому `prompt`,
`content` или `sequence` — красный гейт, а не «строчка для отладки».

**Кто пишет.** Стадии пакета пишут `tokens_in = null`, `tokens_out = null`,
`cost_usd = 0`; это локальный R. Обёртка агента/LLM, которая ведёт программу,
пишет те же поля с реальными числами. Смешанный маршрут — один файл и одна
схема.

**Дневной срез.** `Rscript dev/metrics.R --report <date>` пишет
`data/metrics/metrics-<date>.md`:

- p50/p95 `duration_ms` по стадиям (nearest-rank); стадия с `n < 20` —
  `insufficient_n`, не процентиль;
- пропускная способность: завершённых прогонов в час и в день;
- стоимость на маршрут J1–J7: `sum(cost_usd) / число прогонов маршрута`;
- доля провалов `(fail + error) / total` по стадии и маршруту;
- cache-hit rate и p95 `queue_depth`.

Итоговый отчёт программы — `data/metrics/metrics-<date>.md` с дельтой к
baseline по каждой стадии и маршруту.

**Гейты метрик:**

1. `p95(stage) <= p95_baseline(stage)` для каждой стадии с baseline и `n >= 20`;
   иначе `not_run` с числом сэмплов и причиной.
2. `cost_per_route(JN) <= cost_per_route_baseline(JN)` для J1–J7; если baseline
   равен нулю, любая billed-стоимость — красный.
3. Повторный прогон байт-идентичен: выход маршрута (`artefact.json`, CSV)
   сравнивается байтами при фиксированном `SOURCE_DATE_EPOCH`; сырой JSONL
   байтами не сравнивается, потому что `duration_ms` и `queue_depth` разные.
   Гейт сравнивает канонический срез без `duration_ms`, `queue_depth` и
   `cache_hit`; тайминги проверяются p50/p95.
4. Доля провалов выше baseline — красный дневной отчёт даже при зелёном p95;
   отчёт называет стадию и маршрут.

## 1. День 1 — контракт и поверхность

**Цель:** J1 и J2 исполняются как тесты; батч-отчёт и праймеры перестают быть
тупиками; вход варианта и справочники получают типизированный мост.

| Поток | Пункты | Работа | Файлы | Старт |
| ----- | ------ | ------ | ----- | ----- |
| **A. Контракт** | ZV1 | J1/J2 тесты, фикстуры сценариев, статическая проверка шагов в `check-user-api.R`, ADR 0008 | `tests/testthat/test-scenario-j{1,2}.R`, `tests/testthat/fixtures/scenarios/`, `dev/check-user-api.R` | сразу |
| **B. Типизированный выход** | ZV2, ZV3 | `as_marker_solved(report)`; `primer_candidates`/`primer_set` + `as_primer_set()` | `R/marker_batch.R`, `R/marker_strategy.R`, `R/marker_primers.R`, `R/marker_tubes.R`, `R/marker_assemble*.R` | сразу; мерж после A |
| **C. Вход и реестры** | ZV4, ZV5 | мост `marker_schemes()`/`marker_criteria()` в дизайн; `read_variant_genotypes()`, `as_variant_row()` | `R/marker_io.R`, `R/marker_decision.R`, `R/marker_technology.R`, `R/marker_priors.R`, `R/marker_checks.R` | сразу |
| **D. Метрики** | ZV1 (часть) | `dev/metrics.R` (writer, validator, report), baseline на `main`, JSONL с первого прогона | `dev/metrics.R`, `data/metrics/` | **до A/B/C** |

D — общий слой, а не четвёртая фича: A, B и C не начинают прогоны, пока
writer и baseline не готовы. Иначе первый день даёт код без базовой линии, и
p95-гейт становится не с чем сравнивать.

Потоки не пересекаются по файлам. B зависит от A только на мерже: тест J2
пишется против фикстуры, которую вводит A; разработка идёт параллельно.

**Гейт дня:**

```sh
export PATH=/opt/R-4.6.1/bin:/root/development/tidymarkers/.tools/bin:$PATH
export UNAFOLDDAT=/root/development/tidymarkers/.tools/share/oligoarrayaux
export R_LIBS_SITE=/opt/R-4.6.1/lib/R/library

Rscript dev/test.R "scenario|batch|io|plan|dimer|tubes"
Rscript dev/check-user-api.R
Rscript dev/metrics.R --check <date>
Rscript dev/metrics.R --report <date>
Rscript dev/check.R --quick --no-render
```

**Демо:** J1 (один SNP → ассе → проверки → файлы) и J2 (VCF → батч →
`as_marker_solved()` → проверки/экспорт) зелёные на коммиченных фикстурах.

**Провал дня:** J2 требует `$report`; `report -> marker_solved` нет в графе;
`design_as_primers -> dimer_thermo` всё ещё ложное ребро; `check-user-api.R`
красный; новый экспорт без строки CSV; baseline снят после первой правки кода;
в JSONL попало содержимое или промпт.

## 2. День 2 — панель, петля исходов, масштаб

**Цель:** J3, J4 и J5. Панель становится промежуточным объектом с наблюдаемым
порядком швов; исходы влияют на порядок попыток, не на вердикт; пангеном
A/B/D не смешивается.

| Поток | Пункты | Работа | Файлы | Зависит от |
| ----- | ------ | ------ | ----- | ---------- |
| **A. Панель** | ZV6 | `as_marker_solved(panel)` для Pareto/артефакта/экспорта/replay; порядок GWAS → LD → haplotypes → chemistry → artefact; единый `plan$ld`; блок `panel` в `artefact.json` | `R/marker_panel.R`, `R/marker_pareto.R`, `R/marker_artefact.R`, `R/marker_solve.R` | ZV2 |
| **B. Петля** | ZV7 | `marker_priors(outcomes = )` со счётчиками `works`/`fails`; приор меняет порядок попыток, не вердикт | `R/marker_priors.R`, `R/marker_plan.R`, `R/marker_solve.R` | ZV5 (мягко) |
| **C. Масштаб** | ZV8 | `build =` у genome/index; ключ шарда `(genome_build, chrom)`; отказ от смеси сборок; бюджет 5 Б/п.н. и 2× в `dev/perf-gate.R` | `R/marker_genome.R`, `R/marker_index.R`, `dev/perf-gate.R`, `tests/testthat/fixtures/` | независимо; `T11` |

C не зависит от A/B и может мержиться в любой момент после гейта. A мержится
первым: J3 — демо дня и условие `experimental = FALSE` у панели.

**Гейт дня:**

```sh
Rscript dev/test.R "panel|outcomes|priors|genome|index|scenario"
Rscript dev/bench-index-memory.R
Rscript dev/metrics.R --check <date>
Rscript dev/metrics.R --report <date>
Rscript dev/check.R --quick --no-render
```

**Демо:** J3 (targets + haplotypes + LD + costs → `marker_panel` → Pareto и
артефакт), J4 (три сборки, 6 шардов, один вердикт на сборку), J5 (артефакт →
исходы → приоры → следующий порядок попыток).

**Провал дня:** J3 использует `panel$solved`; нефазированная строка на
LD-отброшенной цели всё ещё роняет панель; LD-record не называет dropped
targets; запрос J4 смешивает A/B/D; `VmHWM` выше `baseline + 5 Б/п.н.`;
`ld = NULL` и `ld = 0.8` дают один `artefact_id`; вердикт J5 изменился; p95
стадии вырос против baseline; стоимость маршрута J3/J4/J5 выросла.

## 3. День 3 — продукт и выпуск

**Цель:** J6 и J7. Артефакт становится релизным продуктом с
`contract.version` и `scenario`; MCP-контракт получает один parity-гейт;
поверхность и граф не гниют.

| Поток | Пункты | Работа | Файлы | Зависит от |
| ----- | ------ | ------ | ----- | ---------- |
| **A. Релиз** | ZV9 | `artefact@2`: `contract.version`, `scenario`, `targets` в `inputs`, `files.solved` с причиной; replay сравнивает вердикты | `R/marker_artefact.R`, `R/marker_export_rocrate.R`, `R/marker_patent.R` | ZV6 |
| **B. MCP parity** | ZV10 | `dev/check-mcp-parity.R`: имена, required/формалы, enum/границы, аннотации, классы в CSV, покрытие J1–J7; шаг в `dev/check.R` | `dev/check-mcp-parity.R`, `dev/check.R`, `docs/user-api-surface.csv` | ZV1 |
| **C. Документы и прогон** | — | README/виньетка панели, `docs/INDEX.md`, `NEWS.md`, сухой релиз | `README.md`, `vignettes/`, `docs/`, `NEWS.md` | A и B |

B не зависит от A по коду, но оба трогают релизный гейт; A мержится первым,
B вторым, C — после обоих. `docs/INDEX.md` и `NEWS.md` не правятся в A/B,
чтобы не ловить конфликт.

**Гейт дня (и релизный гейт):**

```sh
Rscript dev/check-user-api.R
Rscript dev/check-mcp-parity.R
Rscript dev/metrics.R --check <date>
Rscript dev/metrics.R --report <date>
Rscript dev/test.R --full
Rscript dev/check.R --full
```

**Демо:** J6 (два прогона с одним `SOURCE_DATE_EPOCH` байт-в-байт, RO-Crate,
ST.26, replay) и J7 (MCP `tools/list` → `tools/call`, read-only и
`marker_technology`), собранный тарбол.

**Провал дня:** текстовые файлы артефакта различаются; `solved.rds` сравнивают
байтами; `scenario`/`contract.version` отсутствуют; parity-гейт не ловит
удалённый инструмент или незадекларированный аргумент; `R CMD check` красный;
`docs/check-docs.R` находит сироту; `check-user-api.R` находит stale DOT;
нет `data/metrics/metrics-<date>.md`; p95 или стоимость на маршрут выросли,
а отчёт этого не показал.

## 4. Что делается параллельно

| День | Параллельно | Последовательно | Критический путь |
| ---- | ----------- | --------------- | ---------------- |
| 1 | D — первым, затем A (ZV1), B (ZV2/ZV3), C (ZV4/ZV5) | baseline D → мерж A → C → B | ZV1 → ZV2 |
| 2 | A (ZV6), B (ZV7), C (ZV8) | мерж A → B → C | ZV2 → ZV6 |
| 3 | A (ZV9), B (ZV10) | C после A и B | ZV6 → ZV9 |

D — общий слой метрик: он стартует до A/B/C и не мержится отдельно от
baseline-прогона.

Общий критический путь программы: **ZV1 → ZV2 → ZV6 → ZV9**. `ZV10` идёт
параллельно `ZV9` и не удлиняет путь, если `ZV1` закрыт в день 1.

Нельзя параллелить:

- `ZV6` и `ZV9` — оба меняют `artefact.json` и `marker_panel`/`marker_solved`;
- `ZV10` и любой шаг, который правит `dev/check.R` (один файл — один
  владелец);
- `docs/INDEX.md` и `NEWS.md` — только поток C дня 3.

## 5. Гейты: что каждый доказывает

| Гейт | Команда | Доказывает | Красный означает |
| ---- | ------- | ---------- | ---------------- |
| Быстрый локальный | `Rscript dev/check.R --quick --no-render` | roxygen, styler, lintr, критерии, docs/artifacts | код или документы разошлись |
| Сценарии | `Rscript dev/test.R scenario` | J1–J7 исполняются на фикстурах | пользовательский путь сломан |
| Поверхность | `Rscript dev/check-user-api.R` | 207/48, CSV, DOT, подтверждённые рёбра, шаги сценариев | экспорт или ребро не описаны |
| MCP parity | `Rscript dev/check-mcp-parity.R` | 38 инструментов, схема, формалы, аннотации, CSV, J1–J7 | агент видит не тот контракт |
| Память/время | `Rscript dev/bench-index-memory.R` + `dev/perf-gate.R` | 5 Б/п.н., 2×, отсутствие монолита выше 21.5 ГБ | масштаб не доказан |
| Метрики | `Rscript dev/metrics.R --check <date>` и `--report <date>` | JSONL-схема, p50/p95, стоимость маршрута, доля провалов, байт-идентичность выхода | регрессия не измерена или в метрики попало содержимое |
| Полный локальный | `Rscript dev/check.R` | всё выше + фигуры + install smoke | релиз не собран |
| Релизный | `Rscript dev/check.R --full` | `R CMD check --as-cran` + полный сьют | пакет не готов |

`not_run` допустим только с названной причиной (нет 16 ГБ фикстуры, нет
движка). Молчаливый пропуск — красный.

## 6. Критерии провала

### Провал дня

1. Любой гейт дня красный на конец дня.
2. Сценарий дня не запускается без сети, без коммиченной фикстуры или без
   интерактивной сессии.
3. Сценарий проходит только через `skip_if_no_*()` и не проверяет ни одного
   контрактного шага.
4. Пользовательский шаг требует `$`, `[[` или ручного вытаскивания колонки.
5. Новый экспорт не имеет строки в `docs/user-api-surface.csv` или владельца
   в ADR/ZV.
6. Вердикт на золотом наборе изменился при том же входе.
7. Baseline снят после первой правки кода, или хотя бы один прогон не оставил
   JSONL.
8. В JSONL есть `prompt`, `content`, `sequence`, метка цели или путь; схема
   отвергла бы это, а гейт пропустил.
9. p95 стадии выше baseline, стоимость маршрута J1–J7 выросла, или выход
   повторного прогона не байт-идентичен (канонический срез метрик — тоже).

### Провал программы

1. J1–J7 не все зелёные к концу дня 3.
2. `artefact@2` не воспроизводится байт-в-байт при фиксированном
   `SOURCE_DATE_EPOCH`.
3. MCP parity-гейт не ловит хотя бы один из шести классов дрейфа (имена,
   аргументы, границы, аннотации, классы CSV, покрытие сценариев).
4. J4 заявляет 16 Гб-замер без записи в
   `docs/decisions/0001-16gb-fixpur.md`; фикстура выдаётся за операторский
   прогон.
5. `panel$experimental` снят до зелёного J3.
6. Закрытый ZV не имеет теста и гейта; закрытие сделано ссылкой на ADR.
7. `dev/check.R --full`, `check-user-api.R` и `check-mcp-parity.R` не все
   зелёные на одном коммите.
8. Нет итогового `data/metrics/metrics-<date>.md` с p50/p95, стоимостью на
   маршрут, долей провалов и дельтой к baseline; или baseline снят не на
   `main` до первой правки.
9. Хотя бы один маршрут J1–J7 не имеет ни одного JSONL-прогона за три дня.

### Срез объёма (не провал, если объявлен)

- День 1: первым `ZV3`, затем `ZV4`, затем `ZV5`.
- День 2: `ZV7` можно отложить, если J5 остаётся красным по данным, а не по
  коду; `ZV6` и `ZV8` — нет.
- День 3: `ZV9` и `ZV10` — нет; документация потока C сокращается до
  `docs/INDEX.md` + `NEWS.md`.
- Срезанный пункт остаётся `[ ]` в §ZV с названной причиной; он не
  превращается в `[x]`.

## 7. Merge, ветки и откат

- Ветка плана `plan/3day-20261127` содержит только этот план, ADR 0008–0011,
  §ZV и ссылки в индексах. Она мержится первой.
- Каждый ZV — отдельный коммит (или PR) в ветку
  `integration/3day-20261127` от `main`; порядок мержа — из таблицы выше.
- `docs/user-api-surface.csv` и `docs/user-api-graph.dot` пересобираются в том
  же коммите, что и экспорт: `Rscript dev/check-user-api.R --write`, затем
  `git diff --exit-code`.
- Откат ZV — revert его коммита; ADR остаётся `proposed`, пункт §ZV остаётся
  `[ ]`. Частичный мерж, оставляющий CSV/DOT stale, запрещён.
- `NEWS.md` и `docs/INDEX.md` правятся последним коммитом дня, чтобы не
  конфликтовать с потоками.

## 8. Определение готовности

- J1–J7 зелёные на коммиченных фикстурах, без сети.
- ZV1–ZV10 закрыты или срезаны с названной причиной; ни один не закрыт
  ссылкой.
- ADR 0008–0011: `accepted`, если их гейт зелёный; иначе `proposed` с
  названным блокером. Статус в `docs/adr/README.md` совпадает.
- `check-user-api.R` и `check-mcp-parity.R` зелёные; DOT не stale; 207/48
  обновлены, если экспорт добавлен.
- `data/metrics/` содержит baseline, JSONL каждого прогона и итоговый
  `metrics-<date>.md`; p95 не хуже baseline, стоимость на маршрут не выросла,
  повторный прогон байт-идентичен.
- `dev/check.R --full` зелёный; тарбол собирается и ставится.
- `NEWS.md` содержит пользовательские изменения (`artefact@2`, порядок LD,
  новые аксессоры); план помечен `done`, §ZV — по факту.

## 9. Риски и как их снимать

| Риск | Симптом | Снятие |
| ---- | ------- | ------ |
| Нет OligoArrayAux/bowtie | J1/J3/J5/J6 только skip | контрактный шаг без движка всё равно проверяется; движок — операторская установка, не заглушка |
| 16 ГБ fixpur недоступен | J4 не доказывает масштаб | три сборки + synthetic 250 Мб + decision 0001; 16 ГБ — операторский прогон |
| `artefact@2` ломает id | старые файлы не читаются | `marker_artefact_read()` с миграционной заметкой; `NEWS.md` |
| `check-user-api.R` статичен | сценарий зелёный в графе, красный в R | пара: статическая проверка + `test-scenario-jN.R` |
| MCP parity слишком широкая | гейт падает на каждом новом поле | шесть классов дрейфа, названных в ADR 0011; additive-поля не ломают |
| Три дня не хватает | ZV3/ZV4/ZV5 в работе | срез по правилу; гейты и J1/J2/J3/J6/J7 не двигаются |
| Метрики шумят | p95 скачет на малых выборках | `n >= 20` и nearest-rank; иначе `not_run` с числом сэмплов |
| JSONL разрастается | файлы на каждый прогон | только числа; один файл на прогон; ротация по дате, содержимое запрещено |

## 10. Ссылки

- [`ROADMAP.md` §ZU](../../ROADMAP.md) — инвентарь и дефекты поверхности.
- [`ROADMAP.md` §ZV](../../ROADMAP.md) — очередь ZV1–ZV10.
- [`docs/user-api-surface.md`](../user-api-surface.md) — 207/48, классы, граф,
  §8.
- [`docs/user-api-surface.csv`](../user-api-surface.csv) — машинный срез.
- [`docs/adr/0008`](../adr/0008-user-scenarios-as-contract.md),
  [`0009`](../adr/0009-panel-orchestrator-contract.md),
  [`0010`](../adr/0010-index-strategy-at-scale.md),
  [`0011`](../adr/0011-release-artefact-and-mcp-contract.md).
- [`dev/check.R`](../../dev/check.R), [`dev/check-user-api.R`](../../dev/check-user-api.R),
  [`dev/check-docs.R`](../../dev/check-docs.R).
- [`docs/index-perf.md`](../index-perf.md) и
  [`docs/decisions/0001-16gb-fixpur.md`](../decisions/0001-16gb-fixpur.md) —
  числа масштаба.
