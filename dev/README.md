# dev/ — внутренние проверки и измерения

Скрипты разработчика: то, что в обычном пакете делает CI, здесь запускается
локально (CI в репозиторий намеренно не входит — см. `CONTRIBUTING.md`).
Ни один из них не входит в сборку пакета (`.Rbuildignore`).

## `check.R` — приёмка перед коммитом

```sh
Rscript dev/check.R             # быстро: sync доков, styler, lintr
Rscript dev/check.R --full      # + R CMD check (--as-cran)
Rscript dev/check.R --spelling  # + орфография (англ. словарь — шумно на рус. докках)
```

Шаги по порядку: `NAMESPACE`/`man/` синхронны с roxygen → `R/` и `tests/`
форматированы (`styler`) и без линтов (`lintr`) → орфография (`spelling`, только с
`--spelling`) → сборка и проверка (`R CMD check`). Отсутствующий инструмент —
это сообщение «skip», а не падение; провал любого шага даёт ненулевой код возврата.

Орфография вынесена за флаг: проза (README/ROADMAP/NEWS/виньетки) русская, и
англ. словарь помечает каждое слово; английские roxygen/комментарии и так
проверяет `lintr`.

## `bench.R` — время и память по функциям

```sh
Rscript dev/bench.R
```

Прогоняет по представительному вызову горячих функций и печатает медианное время
и объём аллокаций. С `bench` в библиотеке берётся `mem_alloc`; на голой машине —
`system.time()` и базовый `Rprofmem()`. Нужен `OligoArrayAux` (см.
`data-raw/install-oligoarrayaux.R`) — без него термодинамические строки пропускаются.

## Инструменты (в `Suggests`)

`styler`, `lintr`, `spelling`, `bench`, `profmem`, `lobstr`, `pkgdown` (и их
зависимости) запинены в `renv.lock` — ставятся `renv::restore()`. `bench`/`profmem`/
`lobstr` не обязательны: `bench.R` откатывается на `system.time()` + `Rprofmem()`.
