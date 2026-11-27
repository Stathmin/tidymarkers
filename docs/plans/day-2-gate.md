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

## Закалка после гейта

- `Rscript dev/check.R --full`: `R CMD check` чистый, 2 NOTE (CRAN incoming
  feasibility и время примеров); BiocCheck чистый; 210/210 примеров без
  OligoArrayAux; полный набор — `PASS 40963`, `FAIL 0`, `WARN 3`, `SKIP 28`.
- Повторный полный прогон после поздних правок: `PASS 41079`, `FAIL 0`,
  `WARN 4`, `SKIP 28`; `R CMD check` чистый, те же 2 NOTE; 210/210 примеров.
  Прогон поймал регрессию: `marker_costs()` без аргументов падал на пустом
  векторе (`e53946a`), теперь проверка применяется только к переданным стоимостям.
- `Rscript dev/bench.R --gate`: все кейсы в пределах 2x, peak RSS 689 MB.
- `ZU7` закрыт: `design_marker_batch(keep=)` — один батч, `design_markers()` —
  его сводка; `design_kasp_markers()` остаётся отдельным типизированным
  вердиктом J2.
- Найденные и исправленные дефекты ввода: VCF POS вне `integer` (ошибка, не
  `NA`), неизвестная буква FASTA (ошибка, не сдвиг координат), имя файла из
  `artefact.json` не выходит из каталога.

## Поздние проверки входа

После гейта добавлены отказы вместо тихих `NA`: пустая хромосома в VCF
(`ChrNA` больше не появляется), символический ALT по одному правилу обоих
читателей (ADR 0016), дубли имён FASTA, неразбираемое время артефакта, имя гена
как путь, олиго и число до `system2()`, а также цель, гаплотип, стоимость и
порог с недопустимыми значениями.

## Открыто

- ZV12 / T11: шаг `dev/bench.R --gate` в CI стоит; baseline для
  `github-ubuntu-latest` снимается первым прогоном на раннере (ADR 0014).
- Известный предел: ZV13.3 — `marker_replay()` читает RDS, названный записью;
  безопасного `readRDS()` нет, предел назван в ADR 0015.
- Журнал саморемонта не пуст: 12 строк (10 fixed, 1 bypassed, 1 quarantined);
  новые блокеры пишутся через `dev/repair.R` и несут измеренную длительность.
