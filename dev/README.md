# dev/ — внутренние проверки и измерения

Скрипты разработчика: то, что в обычном пакете делает CI, здесь запускается
локально (CI в репозиторий намеренно не входит — см. `CONTRIBUTING.md`).
Ни один из них не входит в сборку пакета (`.Rbuildignore`).

## `check.R` — приёмка перед коммитом

```sh
Rscript dev/check.R          # быстро: sync доков, styler, lintr, spelling
Rscript dev/check.R --full   # + R CMD check (--as-cran)
```

Шаги по порядку: `NAMESPACE`/`man/` синхронны с roxygen → `R/` и `tests/`
форматированы (`styler`) и без линтов (`lintr`) → орфография (`spelling`, словарь
`inst/WORDLIST`) → сборка и проверка (`R CMD check`). Отсутствующий инструмент —
это сообщение «skip», а не падение; провал любого шага даёт ненулевой код возврата.

## `bench.R` — время и память по функциям

```sh
Rscript dev/bench.R
```

Прогоняет по представительному вызову горячих функций и печатает медианное время
и объём аллокаций. С `bench` в библиотеке берётся `mem_alloc`; на голой машине —
`system.time()` и базовый `Rprofmem()`. Нужен `OligoArrayAux` (см.
`data-raw/install-oligoarrayaux.R`) — без него термодинамические строки пропускаются.

## Инструменты (в `Suggests`)

`styler`, `lintr`, `spelling` (и опционально `bench`, `profmem`, `lobstr`) ставятся
через `renv`:

```r
renv::install(c("styler", "lintr", "spelling", "bench"))
```
