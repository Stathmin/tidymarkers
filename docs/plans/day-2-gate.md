# День 2 — отчёт гейта

Дата: 2026-11-27 (гейт пройден раньше календарного дня; правило §0.4).
Ветка: `program/3day-20261127`.

## Что закрыто

| Пункт | Доказательство |
| ----- | -------------- |
| ZV6, панель | `test-scenario-j3.R`: LD → haplotypes → chemistry; `plan$ld` одной формы; блок `panel` в `artefact.json`; `experimental = FALSE`; повторный прогон байт-в-байт |
| ZV7, петля исходов | `test-scenario-j5.R`: артефакт → исход → приоры → порядок попыток; вердикт не меняется |
| ZV8, масштаб | `test-scenario-j4.R`: три сборки, шесть шардов, отсутствующий шард называет сборку и хромосому; `specificity_verdict()` отказывает смеси |
| ZV9, релиз | `test-scenario-j6.R`: `artefact@2`, contract/scenario, RO-Crate, ST.26, replay, два прогона байт-в-байт |
| ZV10, поверхность | `dev/check-mcp-parity.R` и `dev/check-user-api.R` — шаги релизного гейта |

## Прогоны

```sh
TIDYMARKERS_METRICS=day2-gate Rscript dev/test.R \
  "panel|outcomes|priors|genome|index|scenario"
# FAIL 0 | WARN 1 | SKIP 1 | PASS 916

Rscript dev/bench-index-memory.R memory tests/testthat/fixtures/scenarios/j1_genome.fa
# memory  peak_rss_mb=470.8
Rscript dev/bench-index-memory.R shard tests/testthat/fixtures/scenarios/j4_a.fa chr1
# shard   peak_rss_mb=563.4

Rscript dev/metrics.R --check 2026-11-27
# ok: 5 file(s), 282 row(s)
Rscript dev/metrics.R --report 2026-11-27
# exit 0; Gate: core cells ok, scenario cells not_run until both runs share run_id

Rscript dev/check.R --quick --no-render
# all checks passed
```

Память шардов на фикстуре ниже порога 1 Мб, поэтому `dev/perf-gate.R`
классифицирует запись как `not_run` — это ожидаемо и названо в
`docs/index-perf.md`.

## Открыто

- ZV12 / T11: шаг `dev/bench.R --gate` в CI стоит; baseline для
  `github-ubuntu-latest` снимается первым прогоном на раннере (ADR 0014).
- ZU7: остаётся свести три батч-входа к одному с параметром `keep`; алиасы
  `design_marker()`/`write_kasp_results()` не удаляются до мажорной версии.
- Блокеров нет; `docs/plans/self-repair-log.md` пуст.
